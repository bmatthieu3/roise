use std::{collections::{BinaryHeap, HashMap, Vec}, vec};

struct Graph {
    /// Adjacency map
    adj: HashMap<usize, Vec<usize>>,
    nodes: Vec<Node>,
}

struct Node {
    /// cost of the node, can be a difficulty to walk through it
    pub cost: f32,
    /// global heuristic to indicates which node to check 
    pub heuristic: f32,
    /// index
    pub idx: usize,
}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.heuristic.partial_cmp(other.heuristic)
    }
}

impl Graph {
    fn find_path<T>(&mut self, start: usize, end: usize, nodes: &[T]) -> Option<Vec<usize>> {
        let mut came_from: HashMap<usize, usize> = vec![];
        let closed_set = HAhs


        let mut open_list = BinaryHeap::new();
        open_list.push(Node {
            cost: 0.0,
            heuristic: ,// euclidean distance
            idx: start
        });

        while let Some(u) = open_list.pop() {
            if u.idx == end {
                return Some(
                    closed_list.into_iter()
                        .map(|n| n.idx)
                        .collect::<Vec<_>>()
                );
            } else {
                for v in self.adj.get(&u.idx).unwrap() {
                    let cost = u.cost + 1;

                    if !closed_list.contains_key(&v.idx) {
                        open_list.(Node {
                            idx: *v,
                            cost,
                            heuristic: cost + // eucledean distance                            
                        });
                    }
                }
                closed_list.insert(u.idx,);
            }
        }
        // no path to end found
        None
    }
}


