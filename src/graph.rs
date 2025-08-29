use std::collections::{BinaryHeap, HashMap, HashSet};
use std::cmp::Ordering;

use crate::Point2;

struct Graph {
    /// Adjacency map
    adj: HashMap<usize, Vec<usize>>,
}

#[derive(Debug)]
#[derive(PartialEq, Clone)]
struct Node {
    /// cost of the node, can be a difficulty to walk through it
    pub cost: f32,
    /// global heuristic to indicates which node to check 
    pub heuristic: f32,
    /// index
    pub idx: usize,
}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.heuristic.partial_cmp(&self.heuristic)
    }
}


impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        // no NaN expected to be compared
        self.partial_cmp(other).unwrap()
    }
}

impl Eq for Node {}


trait Metric {
    fn barycenter(&self) -> Point2<f32>;

    fn distance2_to<T>(&self, other: &T) -> f32
    where
        T: Metric
    {
        let dp = (self.barycenter() - other.barycenter());

        dp.dot(&dp)
    }
}

impl Metric for Point2<f32> {
    fn barycenter(&self) -> Point2<f32> {
        *self
    }
}

impl Graph {
    fn find_path<T: Metric>(&self, start: usize, end: usize, geometries: &[T]) -> Option<Vec<usize>> {
        let mut came_from: HashMap<usize, usize> = HashMap::new();

        let mut open_list = BinaryHeap::new();
        let mut open_list_set = HashMap::new();
        let node = Node {
            cost: 0.0,
            heuristic: geometries[start].distance2_to(&geometries[end]), // euclidean distance
            idx: start
        };
        open_list_set.insert(start, node.clone());
        open_list.push(node);

        let mut prev_idx = start;

        while let Some(curr) = open_list.pop() {
            open_list_set.remove(&curr.idx).unwrap();

            //came_from.insert(curr.idx, prev_idx);
            //prev_idx = curr.idx;

            if curr.idx == end {
                // We arrived to the end node
                // therefore we can reconstruct the path
                let mut curr_idx = curr.idx;
                let mut path = vec![];
                while let Some(prev_idx) = came_from.get(&curr_idx) {
                    path.push(curr_idx);

                    if curr_idx == start {
                        break;
                    }

                    curr_idx = *prev_idx;
                }
                
                path.reverse();

                return Some(path);
            } else {
                // We are on our way
                for &neigh_idx in self.adj.get(&curr.idx).unwrap() {
                    // Neighbor already included in path are discarded
                    if !came_from.contains_key(&neigh_idx) {
                        let curr_neigh_cost = curr.cost + 1.0;
                        match open_list_set.get(&neigh_idx) {
                            Some(neigh_node) if neigh_node.cost > curr_neigh_cost => {
                                // update the node in the open list if the cost is less
                                let node = Node {
                                    idx: neigh_idx,
                                    cost: curr_neigh_cost,
                                    heuristic: curr_neigh_cost + geometries[neigh_idx].distance2_to(&geometries[end])
                                };
                                open_list_set.insert(neigh_idx, node.clone());
                                open_list.push(node);

                                came_from.insert(neigh_idx, curr.idx);
                            },
                            None => {
                                // update the node in the open list if the cost is less
                                let node = Node {
                                    idx: neigh_idx,
                                    cost: curr_neigh_cost,
                                    heuristic: curr_neigh_cost + geometries[neigh_idx].distance2_to(&geometries[end])
                                };
                                open_list_set.insert(neigh_idx, node.clone());
                                open_list.push(node);

                                came_from.insert(neigh_idx, curr.idx);
                            },
                            _ => ()
                        }
                    }
                }
            }
        }
        // no path to end found
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::Point2;
    use std::collections::HashMap;

    use super::Graph;
    #[test]
    fn test_astar_simple_square() {
        // Build geometry for 4 nodes in a square
        let geometries = vec![
            Point2::new(0.0, 0.0), // node 0
            Point2::new(1.0, 0.0), // node 1
            Point2::new(1.0, 1.0), // node 2
            Point2::new(0.0, 1.0), // node 3
        ];

        // Build adjacency
        let mut adj = HashMap::new();
        adj.insert(0, vec![1, 3]);
        adj.insert(1, vec![0, 2]);
        adj.insert(2, vec![1, 3]);
        adj.insert(3, vec![0, 2]);

        let graph = Graph { adj };

        // Call your A* (assuming signature like):
        // fn astar(graph: &Graph, geometries: &[Point2<f32>], start: usize, goal: usize) -> Option<Vec<usize>>
        let path = graph.find_path(0, 2, &geometries).expect("no path found");

        // Path should be 0 -> 1 -> 2 or 0 -> 3 -> 2
        assert!(
            path == vec![0, 1, 2] || path == vec![0, 3, 2],
            "unexpected path: {:?}",
            path
        );
    }

    #[test]
    fn test_astar_3x3_grid() {
        // 3x3 grid positions (row-major order)
        let geometries = vec![
            Point2::new(0.0, 0.0), // 0
            Point2::new(1.0, 0.0), // 1
            Point2::new(2.0, 0.0), // 2
            Point2::new(0.0, 1.0), // 3
            Point2::new(1.0, 1.0), // 4
            Point2::new(2.0, 1.0), // 5
            Point2::new(0.0, 2.0), // 6
            Point2::new(1.0, 2.0), // 7
            Point2::new(2.0, 2.0), // 8
        ];

        // Build adjacency for 4-connected grid (up, down, left, right)
        let mut adj = HashMap::new();
        let neighbors = [
            vec![(0, 1), (0, 3)],
            vec![(1, 0), (1, 2), (1, 4)],
            vec![(2, 1), (2, 5)],
            vec![(3, 0), (3, 4), (3, 6)],
            vec![(4, 1), (4, 3), (4, 5), (4, 7)],
            vec![(5, 2), (5, 4), (5, 8)],
            vec![(6, 3), (6, 7)],
            vec![(7, 4), (7, 6), (7, 8)],
            vec![(8, 5), (8, 7)],
        ];

        for (node, neighs) in neighbors.iter().enumerate() {
            adj.insert(node, neighs.iter().map(|(_, n)| *n).collect());
        }

        let graph = Graph { adj };

        let path = graph.find_path(0, 8, &geometries).expect("no path found");

        // Path should start at 0 and end at 8
        assert_eq!(*path.first().unwrap(), 0);
        assert_eq!(*path.last().unwrap(), 8);

        // Print path for debugging
        println!("A* path from 0 to 8: {:?}", path);

        // Path length should be minimal (4 steps in Manhattan distance)
        assert_eq!(path.len(), 5);
    }
}


