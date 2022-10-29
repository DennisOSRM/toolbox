mod command_line;
mod deserialize;

use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufWriter,
};

use command_line::Arguments;
use env_logger::{Builder, Env};
use geojson::{feature::Id, Feature, FeatureWriter, Geometry, Value};
use itertools::Itertools;
use log::info;
use toolbox_rs::{
    bounding_box::BoundingBox, convex_hull::monotone_chain, geometry::primitives::FPCoordinate, io,
    partition::PartitionID, space_filling_curve::zorder_cmp, one_iterator::OneIter,
};

use crate::deserialize::binary_partition_file;

// TODO: tool that generate all the runtime data

pub fn main() {
    Builder::from_env(Env::default().default_filter_or("info")).init();

    println!(r#"   ___                       __      __             _        _ "#);
    println!(r#"  / __|    __     __ _      / _|    / _|   ___     | |    __| |"#);
    println!(r#"  \__ \   / _|   / _` |    |  _|   |  _|  / _ \    | |   / _` |"#);
    println!(r#"  |___/   \__|_  \__,_|   _|_|_   _|_|_   \___/   _|_|_  \__,_|"#);
    println!(r#"_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|"#);
    println!(r#""`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"#);
    println!("build from revision: {}", env!("GIT_HASH"));

    // parse and print command line parameters
    let args = <Arguments as clap::Parser>::parse();
    info!("{args}");

    let partition_ids = binary_partition_file(&args.partition_file);
    info!("loaded {} partitions", partition_ids.len());

    let coordinates = io::read_coordinates(&args.coordinates_file);
    info!("loaded {} coordinates", coordinates.len());

    info!("creating and sorting proxy vector");
    let mut known_ids = HashSet::new();
    let mut proxy_vector = Vec::new();
    for (i, partition_id) in partition_ids.iter().enumerate() {
        if !known_ids.contains(partition_id) {
            proxy_vector.push(i);
            known_ids.insert(partition_id);
        }
    }

    info!("number of unique cell ids is {}", proxy_vector.len());

    let mut cells: HashMap<PartitionID, Vec<usize>> = HashMap::new();
    for (i, partition_id) in partition_ids.iter().enumerate() {
        if !cells.contains_key(partition_id) {
            cells.insert(*partition_id, Vec::new());
        }
        cells.get_mut(partition_id).unwrap().push(i);
    }
    let mut hulls = cells
        .iter()
        .map(|(id, indexes)| {
            let cell_coordinates = indexes.iter().map(|i| coordinates[*i]).collect_vec();
            let convex_hull = monotone_chain(&cell_coordinates);
            let bbox = BoundingBox::from_coordinates(&convex_hull);

            (convex_hull, bbox, id)
        })
        .collect_vec();

    if !args.convex_cells_geojson.is_empty() {
        info!("sorting convex cell hulls by Z-order");
        hulls.sort_by(|a, b| zorder_cmp(a.1.center(), b.1.center()));
        info!("writing to {}", &args.convex_cells_geojson);
        serialize_convex_cell_hull_geojson(&hulls, &args.convex_cells_geojson);
    }

    if !args.boundary_nodes_geojson.is_empty() {
        info!("computing geometry of boundary nodes");
        let edges = io::read_graph_into_trivial_edges(&args.graph);
        info!("loaded {} edges", edges.len());
        let boundary_coordinates = edges
            .iter()
            .filter(|edge| partition_ids[edge.source] != partition_ids[edge.target])
            .map(|edge| coordinates[edge.source])
            .collect_vec();
        info!("detection {} boundary nodes", boundary_coordinates.len());

        serialize_boundary_geometry_geojson(&boundary_coordinates, &args.boundary_nodes_geojson);
    }

    // parse level information from integer representation
    info!("decoded level definition: {}", args.level_definition.one_iter().format(" "));

    // TODO: any assertions on the levels possible
    // instantiate graph, extract sub graphs, process cells



    info!("done.");
}

fn serialize_convex_cell_hull_geojson(
    hulls: &[(Vec<FPCoordinate>, BoundingBox, &PartitionID)],
    filename: &str,
) {
    let file = BufWriter::new(File::create(filename).expect("output file cannot be opened"));
    let mut writer = FeatureWriter::from_writer(file);
    for (convex_hull, bbox, id) in hulls {
        // map n + 1 points of the closed polygon into a format that is geojson compliant
        let convex_hull = convex_hull
            .iter()
            .cycle()
            .take(convex_hull.len() + 1)
            .map(|c| {
                // TODO: should this be implemented via the Into<> trait?
                c.to_lon_lat_vec()
            })
            .collect_vec();

        // serialize convex hull polygons as geojson
        let geometry = Geometry::new(Value::Polygon(vec![convex_hull]));

        writer
            .write_feature(&Feature {
                bbox: Some(bbox.into()),
                geometry: Some(geometry),
                id: Some(Id::String(id.to_string())),
                // Features tbd
                properties: None,
                foreign_members: None,
            })
            .unwrap_or_else(|_| panic!("error writing feature: {id}"));
    }
    writer.finish().expect("error writing file");
}

fn serialize_boundary_geometry_geojson(coordinates: &[FPCoordinate], filename: &str) {
    let file = BufWriter::new(File::create(filename).expect("output file cannot be opened"));
    let mut writer = FeatureWriter::from_writer(file);
    for coordinate in coordinates {
        // serialize convex hull polygons as geojson
        let geometry = Geometry::new(Value::Point(coordinate.to_lon_lat_vec()));

        writer
            .write_feature(&Feature {
                bbox: None,
                geometry: Some(geometry),
                id: None,
                // Features tbd
                properties: None,
                foreign_members: None,
            })
            .unwrap_or_else(|_| panic!("error writing feature: {}", coordinate));
    }
    writer.finish().expect("error writing file");
}
