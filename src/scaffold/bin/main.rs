mod command_line;

use std::{collections::HashSet, fs::File, io::BufWriter};

use command_line::Arguments;
use env_logger::{Builder, Env};
use fxhash::FxHashMap;
use geojson::{feature::Id, Feature, FeatureWriter, Geometry, Value};
use itertools::Itertools;
use log::info;
use toolbox_rs::{
    bounding_box::BoundingBox, cell::BaseCell, convex_hull::monotone_chain, edge::InputEdge,
    geometry::primitives::FPCoordinate, io, one_iterator::OneIter, partition::PartitionID,
    space_filling_curve::zorder_cmp,
};

// TODO: tool that generate all the runtime data

pub fn main() {
    Builder::from_env(Env::default().default_filter_or("info")).init();

    println!(r#"   ___                       __      __             _        _ "#);
    println!(r#"  / __|    __     __ _      / _|    / _|   ___     | |    __| |"#);
    println!(r#"  \__ \   / _|   / _` |    |  _|   |  _|  / _ \    | |   / _` |"#);
    println!(r#"  |___/   \__|_  \__,_|   _|_|_   _|_|_   \___/   _|_|_  \__,_|"#);
    println!(r#"_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|"#);
    println!(r#""`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"#);
    println!("build: {}", env!("GIT_HASH"));

    // parse and print command line parameters
    let args = <Arguments as clap::Parser>::parse();
    info!("{args}");

    let partition_ids = io::read_vec_from_file::<PartitionID>(&args.partition_file);
    info!("loaded {} partition ids", partition_ids.len());

    let coordinates = io::read_vec_from_file::<FPCoordinate>(&args.coordinates_file);
    info!("loaded {} coordinates", coordinates.len());

    let edges = io::read_vec_from_file::<InputEdge<usize>>(&args.graph);
    info!("loaded {} edges", edges.len());

    info!("creating and sorting proxy vector");
    let mut known_ids = HashSet::new();
    let mut proxy_vector = Vec::new();
    for (i, partition_id) in partition_ids.iter().enumerate() {
        if !known_ids.contains(partition_id) {
            proxy_vector.push(i);
            known_ids.insert(partition_id);
        }
    }

    proxy_vector.sort();
    info!("number of unique cell ids is {}", proxy_vector.len());

    if !args.convex_cells_geojson.is_empty() {
        info!("generating convex hulls");
        let mut cells: FxHashMap<PartitionID, Vec<usize>> = FxHashMap::default();
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

        info!("sorting convex cell hulls by Z-order");
        hulls.sort_by(|a, b| zorder_cmp(a.1.center(), b.1.center()));
        info!("writing to {}", &args.convex_cells_geojson);
        serialize_convex_cell_hull_geojson(&hulls, &args.convex_cells_geojson);
    }

    if !args.boundary_nodes_geojson.is_empty() {
        info!("computing geometry of boundary nodes");
        let boundary_coordinates = edges
            .iter()
            .filter(|edge| partition_ids[edge.source] != partition_ids[edge.target])
            .map(|edge| coordinates[edge.source])
            .collect_vec();
        info!("detection {} boundary nodes", boundary_coordinates.len());

        serialize_boundary_geometry_geojson(&boundary_coordinates, &args.boundary_nodes_geojson);
    }

    // parse level information from integer representation
    info!(
        "decoded level definition: [{}]",
        args.level_definition.one_iter().format(", ")
    );

    info!("creating all BaseCells");
    // extract subgraphs
    // TODO: can this be done in a faster way without a hash map?
    let mut cells: FxHashMap<PartitionID, BaseCell> = FxHashMap::default();
    for edge in edges {
        let s = edge.source;
        let t = edge.target;
        let source_id = partition_ids[s];
        let target_id = partition_ids[t];
        // edges belong to the cell of their source

        cells.entry(source_id).or_default().edges.push(edge);
        if source_id != target_id {
            // crossing edge
            // sketch of two cells with a crossing directed edge
            //  ..________   ________..
            //           |   |
            //         ~~s-->t~~          i.outgoing_nodes.push(t)
            //    cell i |   | cell k     k.incoming_nodes.push(t)
            //  .._______|   |_______..

            // incoming and outgoing ids in source, target cells
            cells.entry(source_id).or_default().outgoing_nodes.push(t);
            cells.entry(target_id).or_default().incoming_nodes.push(t);
        }
    }
    info!("created {} base cells", cells.len());

    // process cells
    let matrices = cells
        .iter()
        .map(|(key, cell)| {
            // println!("processing cell {}", key);
            (key, cell.process())
        })
        .collect_vec();

    // if let Some(cell) = cells.get(&PartitionID::new(1978734)) {
    //     cell.process();
    // }

    info!("done processing {} base cells.", matrices.len());

    // TODO: process matrix layers

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
