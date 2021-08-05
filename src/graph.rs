use std::ops::Range;

use crate::static_graph::NodeArrayEntry;

pub type NodeID = u32;
pub type EdgeID = u32;

pub trait Graph<T> {
    fn node_range(&self) -> Range<NodeID>;
    fn edge_range(&self, n: NodeID) -> Range<EdgeID>;
    fn number_of_nodes(&self) -> usize;
    fn number_of_edges(&self) -> usize;
    fn begin_edges(&self, n: NodeID) -> EdgeID;
    fn end_edges(&self, n: NodeID) -> EdgeID;
    fn get_out_degree(&self, n: NodeID) -> usize;
    fn target(&self, e: EdgeID) -> NodeID;
    fn data(&self, e: EdgeID) -> &T;
    fn data_mut(&mut self, e: EdgeID) -> &mut T;
    fn find_edge(&self, s: NodeID, t: NodeID) -> Option<EdgeID>;
}
