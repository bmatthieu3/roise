use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

use crate::triangulation::DelaunayTriangulation;
use crate::Point2;
use crate::VertexIdx;

#[derive(Debug)]
pub struct Graph {
    /// Adjacency map
    adj: HashMap<usize, Vec<usize>>,
}

#[derive(Debug, PartialEq, Clone)]
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
        Some(other.cmp(self))
    }
}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        // no NaN expected to be compared
        self.heuristic.partial_cmp(&other.heuristic).unwrap()
    }
}

impl Eq for Node {}

pub trait Metric {
    fn barycenter(&self) -> Point2<f32>;

    fn distance2_to<T>(&self, other: &T) -> f32
    where
        T: Metric,
    {
        let dp = self.barycenter() - other.barycenter();

        dp.dot(&dp)
    }
}

impl Metric for Point2<f32> {
    fn barycenter(&self) -> Point2<f32> {
        *self
    }
}

impl Graph {
    pub fn from_triangulation(triangulation: &DelaunayTriangulation) -> Graph {
        // Build an edge map
        let mut edges: HashMap<(VertexIdx, VertexIdx), usize> = HashMap::new();

        for (i, (u, v, w)) in triangulation.triangles.iter().enumerate() {
            edges.insert((*u, *v), i);
            edges.insert((*v, *w), i);
            edges.insert((*w, *u), i);
        }

        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, (u, v, w)) in triangulation.triangles.iter().enumerate() {
            let mut neigh = vec![];
            if let Some(a) = edges.get(&(*v, *u)) {
                neigh.push(*a);
            }
            if let Some(b) = edges.get(&(*w, *v)) {
                neigh.push(*b);
            }
            if let Some(c) = edges.get(&(*u, *w)) {
                neigh.push(*c);
            }

            adj.insert(i, neigh);
        }

        Graph { adj }
    }

    pub fn find_path<T: Metric>(
        &self,
        start: usize,
        end: usize,
        geometries: &[T],
    ) -> Option<Vec<usize>> {
        let mut came_from: HashMap<usize, usize> = HashMap::new();

        let mut open_list = BinaryHeap::new();
        let mut open_list_set = HashMap::new();
        let node = Node {
            cost: 0.0,
            heuristic: geometries[start].distance2_to(&geometries[end]), // euclidean distance
            idx: start,
        };
        open_list_set.insert(start, node.clone());
        open_list.push(node);

        while let Some(curr) = open_list.pop() {
            let curr_node = open_list_set.remove(&curr.idx);
            if curr_node.is_none() {
                continue;
            }

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
                    if let std::collections::hash_map::Entry::Vacant(e) = came_from.entry(neigh_idx)
                    {
                        let curr_neigh_cost =
                            curr.cost + geometries[curr.idx].distance2_to(&geometries[neigh_idx]);
                        //let curr_neigh_cost = curr.cost + 1.0;

                        let needs_update = match open_list_set.get(&neigh_idx) {
                            Some(neigh_node) => neigh_node.cost > curr_neigh_cost,
                            None => true,
                        };

                        if needs_update {
                            // update the node in the open list if the cost is less
                            let node = Node {
                                idx: neigh_idx,
                                cost: curr_neigh_cost,
                                heuristic: curr_neigh_cost
                                    + geometries[neigh_idx].distance2_to(&geometries[end]),
                            };
                            open_list_set.insert(neigh_idx, node.clone());
                            open_list.push(node);

                            e.insert(curr.idx);
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
    use crate::geometry::closed_polyline::ClosedPolyline;
    use crate::noise::Gradient;
    use image::Rgb;
    use image::RgbImage;
    use imageproc::drawing::draw_line_segment_mut;

    use crate::triangulation::DelaunayTriangulation;
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

    #[test]
    fn test_astar_on_cdt() {
        let gradient = Gradient::new();
        let contours = crate::triangulation::marching_square::extract_isocontours_from_heightmap(
            200,
            |x: Point2<f32>| {
                let noise =
                    gradient.fbm(&(x * 2.1), 0.6, 3.01) * std::f32::consts::FRAC_1_SQRT_2 + 0.5; // in [0, 1]
                noise >= 0.45
            },
        );

        let vertices = contours
            .iter()
            .cloned()
            .flat_map(|ClosedPolyline { mut vertices }| {
                let _ = vertices.pop();
                vertices
            })
            .collect::<Vec<_>>();

        let triangulation = DelaunayTriangulation::from_contours(&contours);
        let barycenters = triangulation
            .triangles
            .iter()
            .map(|(u, v, w)| {
                (u.get_vertex(&vertices) + v.get_vertex(&vertices) + w.get_vertex(&vertices)) / 3.0
            })
            .collect::<Vec<Point2<f32>>>();
        let graph = Graph::from_triangulation(&triangulation);

        let path = graph
            .find_path(100, 600, &barycenters)
            .expect("no path found");

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation {
            for (&idx1, &idx2) in t.iter().zip(t.iter().cycle().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h), // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h), // end point
                    Rgb([69u8, 203u8, 133u8]),                    // RGB colors
                );
            }
        }

        for (&idx1, &idx2) in path.iter().zip(path.iter().skip(1)) {
            draw_line_segment_mut(
                &mut img,
                (barycenters[idx1].x * w, barycenters[idx1].y * h), // start point
                (barycenters[idx2].x * w, barycenters[idx2].y * h), // end point
                Rgb([255u8, 0u8, 0u8]),                             // RGB colors
            );
        }

        img.save("nav_mesh_pathfinding.png").unwrap();
    }

    #[test]
    fn test_graph_from_cdt() {
        let gradient = Gradient::new();
        let contours = crate::triangulation::marching_square::extract_isocontours_from_heightmap(
            200,
            |x: Point2<f32>| {
                let noise =
                    gradient.fbm(&(x * 2.1), 0.6, 3.01) * std::f32::consts::FRAC_1_SQRT_2 + 0.5; // in [0, 1]
                noise >= 0.45
            },
        );

        let vertices = contours
            .iter()
            .cloned()
            .flat_map(|ClosedPolyline { mut vertices }| {
                let _ = vertices.pop();
                vertices
            })
            .collect::<Vec<_>>();

        let triangulation = DelaunayTriangulation::from_contours(&contours);
        let barycenters = triangulation
            .triangles
            .iter()
            .map(|(u, v, w)| {
                (u.get_vertex(&vertices) + v.get_vertex(&vertices) + w.get_vertex(&vertices)) / 3.0
            })
            .collect::<Vec<Point2<f32>>>();

        let Graph { adj } = Graph::from_triangulation(&triangulation);

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation {
            for (&idx1, &idx2) in t.iter().zip(t.iter().cycle().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h), // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h), // end point
                    Rgb([69u8, 203u8, 133u8]),                    // RGB colors
                );
            }
        }

        for (t_idx, neigh_indices) in adj.iter() {
            for n_idx in neigh_indices {
                let neigh_neigh_indices = adj.get(n_idx).unwrap();
                assert!(neigh_neigh_indices.contains(t_idx));
                //if neigh_neigh_indices.contains(t_idx) {
                draw_line_segment_mut(
                    &mut img,
                    (barycenters[*t_idx].x * w, barycenters[*t_idx].y * h), // start point
                    (barycenters[*n_idx].x * w, barycenters[*n_idx].y * h), // end point
                    Rgb([255u8, 0u8, 255u8]),                               // RGB colors
                );
                /*} else {
                    dbg!(t_idx, neigh_neigh_indices);

                    draw_line_segment_mut(
                        &mut img,
                        (barycenters[*t_idx].x * w, barycenters[*t_idx].y * h),              // start point
                        (barycenters[*n_idx].x * w, barycenters[*n_idx].y * h),            // end point
                        Rgb([255u8, 0u8, 0u8]), // RGB colors
                    );
                }*/
                //

                /*for nn_idx in neigh_neigh_indices {
                    draw_line_segment_mut(
                        &mut img,
                        (barycenters[*n_idx].x * w, 1.0 + barycenters[*n_idx].y * h),              // start point
                        (barycenters[*nn_idx].x * w, 1.0 + barycenters[*nn_idx].y * h),            // end point
                        Rgb([255u8, 255u8, 0u8]), // RGB colors
                    );
                }*/
            }
        }

        img.save("nav_mesh_graph.png").unwrap();
    }
}
