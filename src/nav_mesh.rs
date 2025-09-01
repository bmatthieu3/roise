use crate::geometry::coord::{Normed, Vertex};
use crate::{graph, Point2};
use std::collections::HashMap;
use crate::VertexIdx;
use crate::graph::{Graph, Metric};
use crate::triangulation::DelaunayTriangulation;

use crate::Shape;

pub struct NavMesh {
    pub triangulation: DelaunayTriangulation,
    /// portal edges
    pub portals_id: HashMap<[VertexIdx; 2], PortalId>,
    pub portals: Vec<[VertexIdx; 2]>,
    /// the graph containing the adjency map and allowing a* algorithm
    pub graph: Graph,
    // TODO: spatial index
}

pub type PortalId = usize;

impl NavMesh {
    pub fn from_triangulation(triangulation: DelaunayTriangulation) -> NavMesh {
        // Build an edge map
        let mut portals_id: HashMap<[VertexIdx; 2], PortalId> = HashMap::new();
        let mut edge_idx = 0;
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut portals = vec![];

        for ((mut u, mut v), w) in triangulation.vertices.iter() {
            // always treat a edge 1 time
            if u.id() < v.id() && portals_id.get(&[u, v]).is_none() {
                portals_id.insert([u, v], edge_idx);
                portals.push([u, v]);
                edge_idx += 1;
            } else if u.id() > v.id() && portals_id.get(&[v, u]).is_none() {
                portals_id.insert([v, u], edge_idx);
                portals.push([v, u]);
                edge_idx += 1;
            }
        }

        for ([u, v], edge_idx) in portals_id.iter() {
            match (triangulation.adjacent(*u, *v), triangulation.adjacent(*v, *u)) {
                (Some(a), Some(b)) => {
                    let e1 = if a.id() < u.id() {
                        [a, *u]
                    } else {
                        [*u, a]
                    };

                    let e2 = if a.id() < v.id() {
                        [a, *v]
                    } else {
                        [*v, a]
                    };

                     let e3 = if b.id() < u.id() {
                        [b, *u]
                    } else {
                        [*u, b]
                    };

                    let e4 = if b.id() < v.id() {
                        [b, *v]
                    } else {
                        [*v, b]
                    };

                    let neigh_edge_idx1 = *portals_id.get(&e1).unwrap();
                    let neigh_edge_idx2 = *portals_id.get(&e2).unwrap();
                    let neigh_edge_idx3 = *portals_id.get(&e3).unwrap();
                    let neigh_edge_idx4 = *portals_id.get(&e4).unwrap();

                    adj.insert(*edge_idx, vec![
                        neigh_edge_idx1,
                        neigh_edge_idx2,
                        neigh_edge_idx3,
                        neigh_edge_idx4
                    ]);

                    adj.entry(neigh_edge_idx1).and_modify(|e| e.push(*edge_idx)).or_insert(vec![*edge_idx]);
                    adj.entry(neigh_edge_idx2).and_modify(|e| e.push(*edge_idx)).or_insert(vec![*edge_idx]);
                    adj.entry(neigh_edge_idx3).and_modify(|e| e.push(*edge_idx)).or_insert(vec![*edge_idx]);
                    adj.entry(neigh_edge_idx4).and_modify(|e| e.push(*edge_idx)).or_insert(vec![*edge_idx]);

                },
                _ => () // we do not put it the adjacency map
            }
        }

        let graph = Graph::from_adjancy_map(adj);

        Self {
            triangulation,
            portals,
            portals_id,
            graph
        }
    }

    fn find_path_through_portals(&self, start: Point2<f32>, end: Point2<f32>, vertices: &[Point2<f32>]) -> Option<Vec<PortalId>> {
        // find the triangle in which the start and end points are
        use crate::Triangle;

        let mut start_portal_id = None;
        let mut end_portal_id = None;

        for (u, v, w) in &self.triangulation.triangles {
            let u = *u;
            let v = *v;
            let w = *w;
            if Triangle([u, v, w]).contains(&start, vertices) {
                // choose the starting portal
                let e1 = if u.id() < v.id() {
                    [u, v]
                } else {
                    [v, u]
                };
                let e2 = if u.id() < w.id() {
                    [u, w]
                } else {
                    [w, u]
                };
                let e3 = if v.id() < w.id() {
                    [v, w]
                } else {
                    [w, v]
                };

                let a = (e1.barycenter(vertices) - start).magnitude_squared();
                let b = (e2.barycenter(vertices) - start).magnitude_squared();
                let c = (e3.barycenter(vertices) - start).magnitude_squared();

                match (self.portals_id.get(&e1), self.portals_id.get(&e2), self.portals_id.get(&e3)) {
                    (Some(i), Some(j), Some(k)) => {
                        if a < b && a < c {
                            start_portal_id = Some(*i);
                        } else if b < a && b < c {
                            start_portal_id = Some(*j);
                        } else {
                            start_portal_id = Some(*k);
                        }
                    }
                    (Some(i), Some(j), None) => {
                        if a < b {
                            start_portal_id = Some(*i);
                        } else {
                            start_portal_id = Some(*j);
                        }
                    }
                    (None, Some(j), Some(k)) => {
                        if b < c {
                            start_portal_id = Some(*j);
                        } else {
                            start_portal_id = Some(*k);
                        }
                    }
                    (Some(i), None, Some(k)) => {
                        if a < c {
                            start_portal_id = Some(*i);
                        } else {
                            start_portal_id = Some(*k);
                        }
                    }
                    (Some(i), None, None) => {
                        start_portal_id = Some(*i);
                    }
                    (None, Some(j), None) => {
                        start_portal_id = Some(*j);
                    }
                    (None, None, Some(k)) => {
                        start_portal_id = Some(*k);
                    }
                    _ => ()
                }
            }

            if Triangle([u, v, w]).contains(&end, vertices) {
                // choose the starting portal
                let e1 = if u.id() < v.id() {
                    [u, v]
                } else {
                    [v, u]
                };
                let e2 = if u.id() < w.id() {
                    [u, w]
                } else {
                    [w, u]
                };
                let e3 = if v.id() < w.id() {
                    [v, w]
                } else {
                    [w, v]
                };

                let a = (e1.barycenter(vertices) - end).magnitude_squared();
                let b = (e2.barycenter(vertices) - end).magnitude_squared();
                let c = (e3.barycenter(vertices) - end).magnitude_squared();

                match (self.portals_id.get(&e1), self.portals_id.get(&e2), self.portals_id.get(&e3)) {
                    (Some(i), Some(j), Some(k)) => {
                        if a < b && a < c {
                            end_portal_id = Some(*i);
                        } else if b < a && b < c {
                            end_portal_id = Some(*j);
                        } else {
                            end_portal_id = Some(*k);
                        }
                    }
                    (Some(i), Some(j), None) => {
                        if a < b {
                            end_portal_id = Some(*i);
                        } else {
                            end_portal_id = Some(*j);
                        }
                    }
                    (None, Some(j), Some(k)) => {
                        if b < c {
                            end_portal_id = Some(*j);
                        } else {
                            end_portal_id = Some(*k);
                        }
                    }
                    (Some(i), None, Some(k)) => {
                        if a < c {
                            end_portal_id = Some(*i);
                        } else {
                            end_portal_id = Some(*k);
                        }
                    }
                    (Some(i), None, None) => {
                        end_portal_id = Some(*i);
                    }
                    (None, Some(j), None) => {
                        end_portal_id = Some(*j);
                    }
                    (None, None, Some(k)) => {
                        end_portal_id = Some(*k);
                    }
                    _ => ()
                }
            }
        }

        if let (Some(start_portal_id), Some(end_portal_id)) = (start_portal_id, end_portal_id) {
            self.graph
                .find_path(start_portal_id, end_portal_id, &self.portals, &vertices)
        } else {
            None
        }
    }

    pub fn find_path(&self, start: Point2<f32>, end: Point2<f32>, vertices: &[Point2<f32>]) -> Option<Vec<Point2<f32>>> {
        self.find_path_through_portals(start, end, &vertices).and_then(|portal_ids| {
            struct Portal<'a> {
                pub left: &'a Point2<f32>,
                pub right: &'a Point2<f32>
            }

            // funnel algorithm
            // we first need to order all the vertices indexes following the portals from start to the end of path
            let mut vertices_idx_strip = vec![];
            
            let p1 = self.portals[portal_ids[0]];
            let p2 = self.portals[portal_ids[1]];

            let mid_portal = p1.barycenter(&vertices);
            let path_dir = mid_portal - start;

            if p1[0] != p2[0] && p1[0] != p2[1] {
                // 0 is not in common between the 2 portals
                let is_left = path_dir.det(&(*p1[0].get_vertex(&vertices) - start)) > 0.0;
                vertices_idx_strip.push((p1[0], dbg!(is_left)));

                let is_left = path_dir.det(&(*p1[1].get_vertex(&vertices) - start)) > 0.0;
                vertices_idx_strip.push((p1[1], dbg!(is_left)));
            } else {
                let is_left = path_dir.det(&(*p1[1].get_vertex(&vertices) - start)) > 0.0;
                vertices_idx_strip.push((p1[1], dbg!(is_left)));

                let is_left = path_dir.det(&(*p1[0].get_vertex(&vertices) - start)) > 0.0;
                vertices_idx_strip.push((p1[0], dbg!(is_left)));
            }

            for i in 1..(portal_ids.len() - 1) {
                let prev_portal_id = portal_ids[i - 1];
                let curr_portal_id = portal_ids[i];
                let next_portal_id = portal_ids[i + 1];

                let p_prev = self.portals[prev_portal_id];
                let p = self.portals[curr_portal_id];
                let p_next = self.portals[next_portal_id];

                let mid_portal = p.barycenter(&vertices);
                let path_dir = p_next.barycenter(&vertices) - mid_portal;

                if p[0] == p_prev[0] || p[0] == p_prev[1] {
                    // add p[1]
                    let is_left = path_dir.det(&(*p[1].get_vertex(&vertices) - mid_portal)) > 0.0;
                    vertices_idx_strip.push((p[1], dbg!(is_left)));
                } else {
                    // add p[0]
                    let is_left = path_dir.det(&(*p[0].get_vertex(&vertices) - mid_portal)) > 0.0;
                    vertices_idx_strip.push((p[0], dbg!(is_left)));
                }
            }

            // append the last portal vertex
            let prev_portal_id = portal_ids[portal_ids.len() - 2];
            let curr_portal_id = portal_ids[portal_ids.len() - 1];

            let p_prev = self.portals[prev_portal_id];
            let p = self.portals[curr_portal_id];
            let mid_portal = p.barycenter(&vertices);
            let path_dir = end - mid_portal;

            if p[0] == p_prev[0] || p[0] == p_prev[1] {
                // add p[1]
                let is_left = path_dir.det(&(*p[1].get_vertex(&vertices) - mid_portal)) > 0.0;
                vertices_idx_strip.push((p[1], dbg!(is_left)));
            } else {
                // add p[0]
                let is_left = path_dir.det(&(*p[0].get_vertex(&vertices) - mid_portal)) > 0.0;
                vertices_idx_strip.push((p[0], dbg!(is_left)));
            }

            Some(self.apply_funnel(&start, &end, &vertices_idx_strip, vertices))
        })
    }

    fn apply_funnel<'a>(&self, mut apex: &'a Point2<f32>, end: &Point2<f32>, funnel: &[(VertexIdx, bool)], vertices: &'a [Point2<f32>]) -> Vec<Point2<f32>> {
        let (mut funnel_left_idx, mut funnel_right_idx) = if funnel[0].1 { (0, 1) } else { (1, 0) };

        let mut left = funnel[funnel_left_idx].0.get_vertex(vertices);
        let mut right = funnel[funnel_right_idx].0.get_vertex(vertices);

        let mut path = vec![*apex];

        let mut i = 2;
        while i < funnel.len() {
            let v = funnel[i].0.get_vertex(vertices);
            let is_cur_left = funnel[i].1;

            if is_cur_left {
                // v is a 'left' vertex
                if (*v - *apex).det(&(*left - *apex)) > 0.0 {
                    // tighten the funnel by the left
                    left = v;
                    funnel_left_idx = i;

                    // check if the new left is still on the left of the right
                    if (*v - *apex).det(&(*right - *apex)) > 0.0 {
                        // right is left to left!
                        apex = right;
                        path.push(*right);

                        i = funnel_right_idx;
                        //path.extend_from_slice(&self.apply_funnel(right, end, &funnel[(funnel_right_idx+1)..], vertices));
                    }
                }
            } else {
                // v is a 'right' vertex
                if (*right - *apex).det(&(*v - *apex)) > 0.0 {
                    // tighten the funnel by the right
                    right = v;
                    funnel_right_idx = i;

                    // check if the new right is still on the right of the left
                    if (*left - *apex).det(&(*v - *apex)) > 0.0 {
                        // right is left to left!
                        apex = left;
                        path.push(*left);

                        i = funnel_left_idx;
                        //path.extend_from_slice(&self.apply_funnel(left, end, &funnel[(funnel_left_idx+1)..], vertices));
                        //return path;
                    }
                }
            }
            i += 1;
        }

        path.push(*end);

        path
    }
}


#[cfg(test)]
mod tests {
    use image::Rgb;
    use image::RgbImage;
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;

    use crate::nav_mesh::NavMesh;

    use crate::triangulation::DelaunayTriangulation;
    use crate::geometry::closed_polyline::ClosedPolyline;
    use crate::Point2;
    use crate::noise::Gradient;
    use crate::graph::Metric;
    #[test]
    fn test_navmesh() {
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

        let nav_mesh = NavMesh::from_triangulation(triangulation);

        let path = nav_mesh
            .find_path_through_portals(Point2 { x: 0.1, y: 0.1 }, Point2 { x: 0.1, y: 0.6 }, &vertices)
            .expect("no path found");

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in nav_mesh.triangulation {
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
            let b1 = nav_mesh.portals[idx1].barycenter(&vertices);
            let b2 = nav_mesh.portals[idx2].barycenter(&vertices);
            draw_line_segment_mut(
                &mut img,
                (b1.x * w, b1.y * h), // start point
                (b2.x * w, b2.y * h), // end point
                Rgb([255u8, 0u8, 0u8]),                             // RGB colors
            );
        }

        img.save("nav_mesh_portals.png").unwrap();
    }

        #[test]
    fn test_navmesh_funnel() {
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

        let nav_mesh = NavMesh::from_triangulation(triangulation);

        let start = Point2 { x: 0.1, y: 0.1 };
        let end = Point2 { x: 0.8, y: 0.6 };

        let path = nav_mesh
            .find_path(start, end, &vertices)
            .expect("no path found");

        let path_portals = nav_mesh
            .find_path_through_portals(start, end, &vertices)
            .expect("no path found");

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in nav_mesh.triangulation {
            for (&idx1, &idx2) in t.iter().zip(t.iter().cycle().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h), // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h), // end point
                    Rgb([69u8, 203u8, 133u8]),                    // RGB colors
                );
            }
        }

        for &idx1 in path_portals.iter() {
            let [u, v] = nav_mesh.portals[idx1];

            let p1 = u.get_vertex(&vertices);
            let p2 = v.get_vertex(&vertices);

            draw_line_segment_mut(
                &mut img,
                (p1.x * w, p1.y * h), // start point
                (p2.x * w, p2.y * h), // end point
                Rgb([0u8, 255u8, 255u8]),                             // RGB colors
            );
        }

        for (p1, p2) in path.iter().zip(path.iter().skip(1)) {
            draw_line_segment_mut(
                &mut img,
                (p1.x * w, p1.y * h), // start point
                (p2.x * w, p2.y * h), // end point
                Rgb([255u8, 0u8, 0u8]),                             // RGB colors
            );
        }

        img.save("nav_mesh_funnel.png").unwrap();
    }
}
