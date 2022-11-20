mod command_line;
mod serialize;

use command_line::Arguments;
use env_logger::{Builder, Env};
use fxhash::{FxHashMap, FxHashSet};
use itertools::Itertools;
use log::{error, info};
use rayon::prelude::*;
use toolbox_rs::{
    bounding_box::BoundingBox, cell::BaseCell, convex_hull::monotone_chain, edge::InputEdge,
    geometry::primitives::FPCoordinate, io, level_directory::LevelDirectory, one_iterator::OneIter,
    partition::PartitionID, space_filling_curve::zorder_cmp,
};

// TODO: tool to generate all the runtime data

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
    let mut known_ids = FxHashSet::default();
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
        let mut hulls: Vec<_> = cells
            .par_iter()
            .map(|(id, indexes)| {
                let cell_coordinates = indexes.iter().map(|i| coordinates[*i]).collect_vec();
                let convex_hull = monotone_chain(&cell_coordinates);
                let bbox = BoundingBox::from_coordinates(&convex_hull);

                (convex_hull, bbox, id)
            })
            .collect();

        info!("sorting convex cell hulls by Z-order");
        hulls.sort_by(|a, b| zorder_cmp(a.1.center(), b.1.center()));
        info!("writing to {}", &args.convex_cells_geojson);
        serialize::convex_cell_hull_geojson(&hulls, &args.convex_cells_geojson);
    }

    if !args.boundary_nodes_geojson.is_empty() {
        info!("computing geometry of boundary nodes");
        let boundary_coordinates = edges
            .iter()
            .filter(|edge| partition_ids[edge.source] != partition_ids[edge.target])
            .map(|edge| coordinates[edge.source])
            .collect_vec();
        info!("detection {} boundary nodes", boundary_coordinates.len());

        serialize::boundary_geometry_geojson(&boundary_coordinates, &args.boundary_nodes_geojson);
    }

    // parse level information from integer representation, level 0 is implicit
    let mut levels = args.level_definition.one_iter().collect_vec();
    levels.push(0);
    levels.sort();
    let levels = levels;
    if levels.len() > 6 {
        error!("number of levels > 6");
    }

    info!(
        "decoded matrix level definition: [{:?}]",
        levels.iter().format(", ")
    );

    let max_level = *levels.iter().max().unwrap();
    for id in &partition_ids {
        debug_assert!(id.level() as u32 <= max_level);
    }

    info!("instantiate level directory");
    let level_directory = LevelDirectory::new(&partition_ids, &levels);

    info!("creating base cells on level 0 and boundary nodes on all levels");
    // extract subgraphs
    let mut cells = Vec::new();
    cells.resize(levels.len(), FxHashMap::<PartitionID, BaseCell>::default());
    // TODO: can this be done in a faster way without a hash map, perhaps in parallel?

    // create first overlay layer and note all boundary nodes at their cells
    for edge in edges {
        let s = edge.source;
        let t = edge.target;
        let source_id = partition_ids[s];
        let target_id = partition_ids[t];

        // edges belong to the cell of their source, always
        cells[0].entry(source_id).or_default().edges.push(edge);

        level_directory
            .get_crossing_levels(s, t)
            .iter()
            .enumerate()
            .for_each(|(index, level)| {
                // TODO: is the level needed?
                // crossing edge
                // sketch of two cells with a crossing directed edge
                //  ..________   ________..
                //           |   |
                //         ~~s-->t~~          i.outgoing_nodes.push(t)
                //    cell i |   | cell k     k.incoming_nodes.push(t)
                //  .._______|   |_______..

                // incoming and outgoing ids in source, target cells
                cells[index]
                    .entry(source_id.parent_at_level(*level))
                    .or_default()
                    .outgoing_nodes
                    .push(t);
                cells[index]
                    .entry(target_id.parent_at_level(*level))
                    .or_default()
                    .incoming_nodes
                    .push(t);
            });
    }
    levels.iter().enumerate().for_each(|(index, level)| {
        info!(
            "[{:02}] created {} base cells on level {level}",
            index,
            cells[index].len()
        );
    });

    let matrices = (0..levels.len())
        .map(|index| {
            // process cells on level i
            let level_matrices: Vec<_> = cells[index]
                .par_iter()
                .map(|(key, cell)| (*key, cell.process()))
                .collect();

            // the number of matrices needs to equal the number of cells on a given level
            debug_assert_eq!(level_matrices.len(), cells[index].len());

            info!(
                "[{:02}-{:02}] created {} matrices on",
                index,
                levels[index],
                level_matrices.len(),
            );

            if index < levels.len() - 1 {
                // no need to extract edges for the very last level
                level_matrices.iter().for_each(|(key, cell)| {
                    let parent = key.parent_at_level(levels[index + 1]);
                    cells[index + 1]
                        .entry(parent)
                        .or_default()
                        .edges
                        .extend(cell.overlay_edges());
                });
            }
            level_matrices
        })
        .collect_vec();

    let sum: usize = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices.iter().fold(0, |acc, (_id, cell)| {
                acc + 4
                    * (1 + cell.incoming_nodes.len()
                        + cell.outgoing_nodes.len()
                        + cell.matrix.len())
            })
        })
        .sum();
    println!("bytes: {sum}");

    let avg_incoming = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices
                .iter()
                .fold(0., |acc, (_, cell)| acc + cell.incoming_nodes.len() as f32)
                / level_matrices.len() as f32
        })
        .collect_vec();
    info!("avg_incoming: [{:?}]", avg_incoming.iter().format(", "));

    let max_incoming = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices
                .iter()
                .max_by_key(|(_, cell)| cell.incoming_nodes.len())
        })
        .collect_vec();
    let max_incoming = max_incoming
        .iter()
        .map(|elem| elem.unwrap().1.incoming_nodes.len())
        .collect_vec();
    info!("max_incoming: [{:?}]", max_incoming.iter().format(", "));
    let min_incoming = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices
                .iter()
                .min_by_key(|(_, cell)| cell.incoming_nodes.len())
        })
        .collect_vec();
    let min_incoming = min_incoming
        .iter()
        .map(|elem| elem.unwrap().1.incoming_nodes.len())
        .collect_vec();
    info!("min_incoming: [{:?}]", min_incoming.iter().format(", "));

    let avg_outgoing = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices
                .iter()
                .fold(0., |acc, (_, cell)| acc + cell.outgoing_nodes.len() as f32)
                / level_matrices.len() as f32
        })
        .collect_vec();
    info!("avg_outgoing: [{:?}]", avg_outgoing.iter().format(", "));
    let max_outgoing = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices
                .iter()
                .max_by_key(|(_, cell)| cell.outgoing_nodes.len())
        })
        .collect_vec();
    let max_outgoing = max_outgoing
        .iter()
        .map(|elem| elem.unwrap().1.outgoing_nodes.len())
        .collect_vec();
    info!("max_outgoing: [{:?}]", max_outgoing.iter().format(", "));
    let min_incoming = matrices
        .iter()
        .map(|level_matrices| {
            level_matrices
                .iter()
                .min_by_key(|(_, cell)| cell.incoming_nodes.len())
        })
        .collect_vec();
    let min_incoming = min_incoming
        .iter()
        .map(|elem| elem.unwrap().1.incoming_nodes.len())
        .collect_vec();
    info!("min_incoming: [{:?}]", min_incoming.iter().format(", "));

    // TODO: sum up the size of the structures
    info!("done.");
}
