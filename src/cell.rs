use crate::{edge::SimpleEdge, graph::NodeID};

pub struct UnprocessedCell {
    boundary: Vec<NodeID>,
    edges: Vec<SimpleEdge>,
}

impl UnprocessedCell {
    pub fn process(&self) -> ProcessedCell {
        // sort unique erase boundary nodes
        let mut boundary = self.boundary.clone();
        boundary.sort_unstable();
        boundary.dedup();
        // TODO: renumber graph
        // TODO: compute clique information (first naive point to point, next one-to-many)
        let matrix = vec![usize::MAX; boundary.len() * boundary.len()];
        // return ProcessedCell
        ProcessedCell { boundary, matrix }
    }
}

pub struct ProcessedCell {
    // sorted list of of original boundary node ids
    boundary: Vec<NodeID>,
    // matrix of pairwise distances between boundary nodes
    matrix: Vec<usize>,
}

impl ProcessedCell {
    // TODO: iterator for row and column access
    pub fn get_distance_row(&self, u: NodeID) -> &[usize] {
        // get row in matrix for node u
        let index = self
            .boundary.binary_search(&u)
            .expect(&format!("node {u} not found in node boundary"));
        let node_count = self.boundary.len();
        // return the row of the matrix
        &self.matrix[index * node_count..(index + 1) * node_count]
    }
}
