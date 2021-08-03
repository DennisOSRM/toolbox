use std::cmp::max;

pub type NodeID = u32;
pub type EdgeID = u32;

#[derive(Eq, PartialOrd, Ord, PartialEq)]
pub struct InputEdge<EdgeDataT: Eq> {
    source: NodeID,
    target: NodeID,
    edge_data: EdgeDataT,
}

impl<EdgeDataT: Eq> InputEdge<EdgeDataT> {
    pub fn new(s: NodeID, t: NodeID, d: EdgeDataT) -> Self {
        Self {
            source: s,
            target: t,
            edge_data: d,
        }
    }
}

pub struct NodeArrayEntry {
    first_edge: EdgeID,
}

impl NodeArrayEntry {
    pub fn new(e: EdgeID) -> NodeArrayEntry {
        NodeArrayEntry { first_edge: e }
    }
}

pub struct EdgeArrayEntry<EdgeDataT> {
    target: NodeID,
    data: EdgeDataT,
}

pub struct StaticGraph<T: Ord> {
    node_array: Vec<NodeArrayEntry>,
    edge_array: Vec<EdgeArrayEntry<T>>,
}

impl<T: Ord + Copy> StaticGraph<T> {
    pub const INVALID_NODE_ID: NodeID = NodeID::MAX;
    pub const INVALID_EDGE_ID: EdgeID = EdgeID::MAX;

    pub fn default() -> Self {
        Self {
            node_array: Vec::new(),
            edge_array: Vec::new(),
        }
    }
    pub fn new(mut input: Vec<InputEdge<T>>) -> Self {
        // TODO: renumber IDs if necessary
        let number_of_edges = input.len();
        let mut number_of_nodes = 0;
        for edge in &input {
            number_of_nodes = max(edge.source, number_of_nodes);
            number_of_nodes = max(edge.target, number_of_nodes);
        }

        let mut graph = Self::default();
        // +1 as we are going to add one sentinel node at the end
        graph.node_array.reserve(number_of_nodes as usize + 1);
        graph.edge_array.reserve(number_of_edges);

        // sort input edges by source/target/data
        // TODO(dl): sorting by source suffices to construct adjacency array
        input.sort();

        // add first entry manually, rest will be computed
        graph.node_array.push(NodeArrayEntry::new(0));
        let mut offset = 0;
        for i in 0..(number_of_nodes) {
            while offset != input.len() && input[offset].source == i {
                offset += 1;
            }
            graph.node_array.push(NodeArrayEntry::new(offset as EdgeID));
        }

        //add sentinel at the end of the node array
        graph
            .node_array
            .push(NodeArrayEntry::new((input.len()) as EdgeID));

        for edge in &input {
            graph.edge_array.push(EdgeArrayEntry {
                target: edge.target,
                data: edge.edge_data,
            });
        }
        graph
    }

    pub fn node_range(&self) -> std::ops::Range<NodeID> {
        std::ops::Range {
            start: 0,
            end: self.number_of_nodes() as NodeID,
        }
    }

    pub fn edge_range(&self, n: NodeID) -> std::ops::Range<EdgeID> {
        std::ops::Range {
            start: self.begin_edges(n),
            end: self.end_edges(n),
        }
    }

    pub fn number_of_nodes(&self) -> usize {
        // minus two because of off-by-one _and_ a sentinel at the end
        // TODO(DL): there must be a more elegant way
        self.node_array.len() - 1
    }

    pub fn number_of_edges(&self) -> usize {
        self.edge_array.len()
    }

    pub fn begin_edges(&self, n: NodeID) -> EdgeID {
        self.node_array[n as usize].first_edge
    }

    pub fn end_edges(&self, n: NodeID) -> EdgeID {
        self.node_array[(n + 1) as usize].first_edge
    }

    pub fn get_out_degree(&self, n: NodeID) -> usize {
        let up = self.end_edges(n);
        let down = self.begin_edges(n);
        (up - down ) as usize
    }

    pub fn target(&self, e: EdgeID) -> NodeID {
        self.edge_array[e as usize].target
    }

    pub fn data(&self, e: EdgeID) -> &T {
        &self.edge_array[e as usize].data
    }

    /// Returns whether the graph contains a cycle by running a node
    /// coloring DFS
    pub fn cycle_check(&self) -> bool {
        #[derive(Clone, PartialEq)]
        enum Colors {
            White,
            Grey,
            Black,
        }

        let mut node_colors = Vec::new();
        node_colors.resize(self.number_of_nodes(), Colors::White);
        let mut stack = Vec::new();

        for root in self.node_range() {
            if node_colors[root as usize] != Colors::White {
                continue;
            }

            stack.push(root);
            while let Some(&node) = stack.last() {
                // pre-order traversal
                if node_colors[node as usize] != Colors::Grey {
                    node_colors[node as usize] = Colors::Grey;

                    for edge in self.edge_range(node) {
                        // push unvisited children to stack
                        let target = self.target(edge);

                        if node_colors[target as usize] == Colors::White {
                            stack.push(target);
                        } else if node_colors[target as usize] == Colors::Grey {
                            // cycle detected
                            return true;
                        }
                    }
                } else if node_colors[node as usize] == Colors::Grey {
                    // post-order traversal
                    stack.pop();
                    node_colors[node as usize] = Colors::Black;
                }
            }
        }
        false
    }
}