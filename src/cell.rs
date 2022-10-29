use std::collections::HashMap;

use itertools::Itertools;
use log::debug;

use crate::{
    edge::InputEdge, graph::NodeID, static_graph::StaticGraph,
    unidirectional_dijkstra::UnidirectionalDijkstra,
};

pub struct BaseCell {
    incoming_nodes: Vec<NodeID>,
    outgoing_nodes: Vec<NodeID>,
    edges: Vec<InputEdge<usize>>,
    // TODO: add renumbering table to support unpacking edges
}

impl BaseCell {
    pub fn process(&self) -> MatrixCell {
        // renumber nodes to be in contiguous range"
        // [sources..targets..other nodes]
        // [0,1,2         ...         n-1]
        let mut seen_nodes = HashMap::new();
        self.incoming_nodes.iter().for_each(|node| {
            debug!("source {}, idx: {}", *node, seen_nodes.len());
            seen_nodes.insert(*node, seen_nodes.len());
        });
        self.outgoing_nodes.iter().for_each(|node| {
            debug!("target {}, idx: {}", *node, seen_nodes.len());
            seen_nodes.insert(*node, seen_nodes.len());
        });

        assert_eq!(
            // assert incoming outgoing nodes ids have empty set intersection
            seen_nodes.len(),
            self.incoming_nodes.len() + self.outgoing_nodes.len()
        );

        let new_edges = self
            .edges
            .iter()
            .map(|edge| {
                // renumber source/target nodes of edges to be in range [0..k-1]
                let mut new_edge = *edge;

                let seen_nodes_len = seen_nodes.len();
                seen_nodes.entry(edge.source).or_insert(seen_nodes_len);
                new_edge.source = *seen_nodes.get(&edge.source).expect("renumbering broken");

                let seen_nodes_len = seen_nodes.len();
                seen_nodes.entry(edge.target).or_insert(seen_nodes_len);
                new_edge.target = *seen_nodes.get(&edge.target).expect("renumbering broken");

                new_edge
            })
            .collect_vec();

        // instantiate subgraph
        let graph = StaticGraph::new(new_edges);
        let mut dijkstra = UnidirectionalDijkstra::new();
        let mut matrix = vec![usize::MAX; self.incoming_nodes.len() * self.outgoing_nodes.len()];

        let source_range = 0..self.incoming_nodes.len();
        for source in source_range {
            // TODO: compute clique information by one-to-many or many-to-many
            let target_range =
                self.incoming_nodes.len()..self.incoming_nodes.len() + self.outgoing_nodes.len();
            for target in target_range {
                let distance = dijkstra.run(&graph, source, target);
                debug!(
                    "matrix[{}] distance({source},{target})={distance}",
                    (source * self.outgoing_nodes.len() + target - self.incoming_nodes.len())
                );

                matrix[source * self.outgoing_nodes.len() + target - self.incoming_nodes.len()] =
                    distance;
            }
        }

        // return MatrixCell
        MatrixCell {
            matrix,
            incoming_nodes: self.incoming_nodes.clone(),
            outgoing_nodes: self.outgoing_nodes.clone(),
        }
    }
}

pub struct MatrixCell {
    incoming_nodes: Vec<NodeID>,
    outgoing_nodes: Vec<NodeID>,
    // matrix of pairwise distances between boundary nodes
    matrix: Vec<usize>,
}

impl MatrixCell {
    // TODO: iterator for row and column access
    pub fn get_distance_row(&self, u: NodeID) -> &[usize] {
        // get row in matrix for node u
        let index = self
            .incoming_nodes
            .binary_search(&u)
            .unwrap_or_else(|_| panic!("node {u} not found in node boundary"));
        // return the row of the matrix
        &self.matrix[index * self.incoming_nodes.len()..(index + 1) * self.outgoing_nodes.len()]
    }
}

#[cfg(test)]
mod tests {
    use super::BaseCell;
    use crate::{edge::InputEdge, graph::UNREACHABLE};

    #[test]
    fn process_base_cell1() {
        let edges: Vec<InputEdge<usize>> = vec![
            InputEdge::new(0, 1, 3),
            InputEdge::new(1, 2, 3),
            InputEdge::new(4, 2, 1),
            InputEdge::new(2, 3, 6),
            InputEdge::new(0, 4, 2),
            InputEdge::new(4, 5, 2),
            InputEdge::new(5, 3, 7),
            InputEdge::new(1, 5, 2),
        ];

        // check first set of source, target nodes
        let incoming_nodes = vec![0, 4];
        let outgoing_nodes = vec![3, 5];

        let base_cell = BaseCell {
            incoming_nodes: incoming_nodes.clone(),
            outgoing_nodes: outgoing_nodes.clone(),
            edges: edges.clone(),
        };
        let matrix_cell = base_cell.process();

        assert_eq!(incoming_nodes, matrix_cell.incoming_nodes);
        assert_eq!(outgoing_nodes, matrix_cell.outgoing_nodes);
        assert_eq!(matrix_cell.matrix, vec![9, 4, 7, 2]);

        assert_eq!(matrix_cell.get_distance_row(0), vec![9, 4]);
        assert_eq!(matrix_cell.get_distance_row(4), vec![7, 2]);

        // check second set of source, target nodes
        let incoming_nodes = vec![0, 1];
        let outgoing_nodes = vec![4, 5];

        let base_cell = BaseCell {
            incoming_nodes: incoming_nodes.clone(),
            outgoing_nodes: outgoing_nodes.clone(),
            edges: edges.clone(),
        };
        let matrix_cell = base_cell.process();

        assert_eq!(incoming_nodes, matrix_cell.incoming_nodes);
        assert_eq!(outgoing_nodes, matrix_cell.outgoing_nodes);
        assert_eq!(matrix_cell.matrix, vec![2, 4, UNREACHABLE, 2]);

        assert_eq!(matrix_cell.get_distance_row(0), vec![2, 4]);
        assert_eq!(matrix_cell.get_distance_row(1), vec![UNREACHABLE, 2]);
    }

    #[test]
    fn process_base_cell2() {
        let edges = vec![
            InputEdge::new(0, 1, 7),
            InputEdge::new(0, 2, 3),
            InputEdge::new(1, 2, 1),
            InputEdge::new(1, 3, 6),
            InputEdge::new(2, 4, 8),
            InputEdge::new(3, 5, 2),
            InputEdge::new(3, 2, 3),
            InputEdge::new(4, 3, 2),
            InputEdge::new(4, 5, 8),
        ];

        // check first set of source, target nodes
        let incoming_nodes = vec![0, 1];
        let outgoing_nodes = vec![4, 5];

        let base_cell = BaseCell {
            incoming_nodes: incoming_nodes.clone(),
            outgoing_nodes: outgoing_nodes.clone(),
            edges: edges.clone(),
        };
        let matrix_cell = base_cell.process();

        assert_eq!(incoming_nodes, matrix_cell.incoming_nodes);
        assert_eq!(outgoing_nodes, matrix_cell.outgoing_nodes);
        assert_eq!(matrix_cell.matrix, vec![11, 15, 9, 8]);

        assert_eq!(matrix_cell.get_distance_row(0), vec![11, 15]);
        assert_eq!(matrix_cell.get_distance_row(1), vec![9, 8]);

        // check second set of source, target nodes
        let incoming_nodes = vec![0, 2];
        let outgoing_nodes = vec![3, 5];

        let base_cell = BaseCell {
            incoming_nodes: incoming_nodes.clone(),
            outgoing_nodes: outgoing_nodes.clone(),
            edges: edges.clone(),
        };
        let matrix_cell = base_cell.process();

        assert_eq!(incoming_nodes, matrix_cell.incoming_nodes);
        assert_eq!(outgoing_nodes, matrix_cell.outgoing_nodes);
        assert_eq!(matrix_cell.matrix, vec![13, 15, 10, 12]);

        assert_eq!(matrix_cell.get_distance_row(0), vec![11, 15]);
        assert_eq!(matrix_cell.get_distance_row(3), vec![9, 8]);
    }

    #[test]
    #[should_panic]
    fn matrix_cell_row_invalid() {
        let edges: Vec<InputEdge<usize>> = vec![
            InputEdge::new(0, 1, 3),
            InputEdge::new(1, 2, 3),
            InputEdge::new(4, 2, 1),
            InputEdge::new(2, 3, 6),
            InputEdge::new(0, 4, 2),
            InputEdge::new(4, 5, 2),
            InputEdge::new(5, 3, 7),
            InputEdge::new(1, 5, 2),
        ];

        // check first set of source, target nodes
        let incoming_nodes = vec![0, 4];
        let outgoing_nodes = vec![3, 5];

        let base_cell = BaseCell {
            incoming_nodes: incoming_nodes.clone(),
            outgoing_nodes: outgoing_nodes.clone(),
            edges: edges.clone(),
        };
        let matrix_cell = base_cell.process();
        matrix_cell.get_distance_row(1);
    }
}
