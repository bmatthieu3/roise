use crate::geometry::coord::{Normed, Vertex};
use crate::{graph, Point2, Triangle};
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

struct EndPointPortal {
    id: usize,
    triangle: Triangle,
}

#[derive(Debug)]
struct Portal<'a> {
    pub left: &'a Point2<f32>,
    pub right: &'a Point2<f32>
}

fn is_left_to(s: &Point2<f32>, a: &Point2<f32>, b: &Point2<f32>) -> bool {
    (*b - *s).det(&(*a - *s)) < 0.0
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

    fn find_triangle_containing_point(&self, p: &Point2<f32>, vertices: &[Point2<f32>]) -> Option<Triangle> {
        // TODO use spatial index here
        for (u, v, w) in &self.triangulation.triangles {
            let u = *u;
            let v = *v;
            let w = *w;
            let tri = Triangle([u, v, w]);
            if tri.contains(p, vertices) {
                return Some(tri);
            }
        }

        None
    }

   

    fn find_nearest_portal(&self, p: &Point2<f32>, vertices: &[Point2<f32>]) -> Option<EndPointPortal> {
        self.find_triangle_containing_point(p, vertices).and_then(|t| {
            let Triangle([u, v, w]) = t.clone();

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

            let a = (e1.barycenter(vertices) - *p).magnitude_squared();
            let b = (e2.barycenter(vertices) - *p).magnitude_squared();
            let c = (e3.barycenter(vertices) - *p).magnitude_squared();

            match (self.portals_id.get(&e1), self.portals_id.get(&e2), self.portals_id.get(&e3)) {
                (Some(i), Some(j), Some(k)) => {
                    if a < b && a < c {
                        Some(EndPointPortal {
                            id: *i,
                            triangle: t
                        })
                    } else if b < a && b < c {
                        Some(EndPointPortal {
                            id: *j,
                            triangle: t
                        })
                    } else {
                        Some(EndPointPortal {
                            id: *k,
                            triangle: t
                        })
                    }
                }
                (Some(i), Some(j), None) => {
                    if a < b {
                        Some(EndPointPortal {
                            id: *i,
                            triangle: t
                        })
                    } else {
                        Some(EndPointPortal {
                            id: *j,
                            triangle: t
                        })
                    }
                }
                (None, Some(j), Some(k)) => {
                    if b < c {
                        Some(EndPointPortal {
                            id: *j,
                            triangle: t
                        })
                    } else {
                        Some(EndPointPortal {
                            id: *k,
                            triangle: t
                        })
                    }
                }
                (Some(i), None, Some(k)) => {
                    if a < c {
                        Some(EndPointPortal {
                            id: *i,
                            triangle: t
                        })
                    } else {
                        Some(EndPointPortal {
                            id: *k,
                            triangle: t
                        })
                    }
                }
                (Some(i), None, None) => {
                    Some(EndPointPortal {
                        id: *i,
                        triangle: t
                    })
                }
                (None, Some(j), None) => {
                    Some(EndPointPortal {
                        id: *j,
                        triangle: t
                    })
                }
                (None, None, Some(k)) => {
                    Some(EndPointPortal {
                        id: *k,
                        triangle: t
                    })
                }
                _ => None
            }
        })
    }

    fn find_path_through_portals(&self, start: Point2<f32>, end: Point2<f32>, vertices: &[Point2<f32>]) -> Option<Vec<PortalId>> {
        // find the triangle in which the start and end points are
        use crate::Triangle;

        let start_portal = self.find_nearest_portal(&start, vertices);
        let end_portal = self.find_nearest_portal(&end, vertices);

        if let (Some(EndPointPortal { id: start_id, triangle: Triangle([su, sv, sw]) }), Some(EndPointPortal { id: end_id, triangle: Triangle([eu, ev, ew]) })) = (start_portal, end_portal) {
            let mut portals = self.graph
                .find_path(start_id, end_id, &self.portals, &vertices);

            if let Some(mut portals) = portals.as_mut() {
                if portals.len() > 1 {
                    let p1 = self.portals[portals[0]];
                    let p2 = self.portals[portals[1]];

                    // check if p1 and p2 belongs to the starting triangle
                    let t1_in_portal = su == p1[0] || su == p1[1] || su == p2[0] || su == p2[1];
                    let t2_in_portal = sv == p1[0] || sv == p1[1] || sv == p2[0] || sv == p2[1];
                    let t3_in_portal = sw == p1[0] || sw == p1[1] || sw == p2[0] || sw == p2[1];

                    if t1_in_portal && t2_in_portal && t3_in_portal {
                        portals.remove(0);
                    }
                }

                if portals.len() > 1 {
                    let p1 = self.portals[portals[portals.len() - 1]];
                    let p2 = self.portals[portals[portals.len() - 2]];

                    // check if p1 and p2 belongs to the ending triangle
                    let t1_in_portal = eu == p1[0] || eu == p1[1] || eu == p2[0] || eu == p2[1];
                    let t2_in_portal = ev == p1[0] || ev == p1[1] || ev == p2[0] || ev == p2[1];
                    let t3_in_portal = ew == p1[0] || ew == p1[1] || ew == p2[0] || ew == p2[1];

                    if t1_in_portal && t2_in_portal && t3_in_portal {
                        portals.pop();
                    }
                }
            }

            portals
        } else {
            None
        }
    }

    pub fn find_path_portals<'a>(&self, start: Point2<f32>, end: Point2<f32>, vertices: &'a [Point2<f32>]) -> Option<Vec<Portal<'a>>> {
        self.find_path_through_portals(start, end, &vertices).and_then(move |portal_ids| {
            // funnel algorithm
            // we first need to order all the vertices indexes following the portals from start to the end of path
            let mut portals = Vec::with_capacity(portal_ids.len());

            let p = self.portals[portal_ids[0]];

            let mut left = p[0].get_vertex(&vertices);
            let mut right = p[1].get_vertex(&vertices);
            if is_left_to(&start, right, left) {
                std::mem::swap(&mut right, &mut left);
            }
            portals.push(Portal {
                left,
                right
            });

            for i in 1..portal_ids.len() {
                let curr_portal_id = portal_ids[i];
                let prev_portal_id = portal_ids[i - 1];

                let p = self.portals[curr_portal_id];
                let p_prev = self.portals[prev_portal_id];

                let prev_mid_portal = p_prev.barycenter(&vertices);

                let mut left = p[0].get_vertex(&vertices);
                let mut right = p[1].get_vertex(&vertices);
                if is_left_to(&prev_mid_portal, right, left) {
                    std::mem::swap(&mut right, &mut left);
                }
                portals.push(Portal {
                    left,
                    right
                });
            }

            Some(portals)

            // append the last portal vertex
            /*let p = self.portals[portal_ids[portal_ids.len() - 1]];
            let mid_portal = p.barycenter(&vertices);
            let path_dir = end - mid_portal;

            let mut left = p[0].get_vertex(&vertices);
            let mut right = p[1].get_vertex(&vertices);
            if path_dir.det(&(*right - mid_portal)) < 0.0 {
                std::mem::swap(&mut right, &mut left);
            }
            portals.push(Portal {
                left,
                right
            });*/
        })
    }


    pub fn find_path(&self, start: Point2<f32>, end: Point2<f32>, vertices: &[Point2<f32>]) -> Option<Vec<Point2<f32>>> {
        self.find_path_through_portals(start, end, &vertices).and_then(|portal_ids| {
            // funnel algorithm
            // we first need to order all the vertices indexes following the portals from start to the end of path
            let mut portals = Vec::with_capacity(portal_ids.len());

            let p = self.portals[portal_ids[0]];

            let mut left = p[0].get_vertex(&vertices);
            let mut right = p[1].get_vertex(&vertices);
            if is_left_to(&start, right, left) {
                std::mem::swap(&mut right, &mut left);
            }
            portals.push(Portal {
                left,
                right
            });

            for i in 1..portal_ids.len() {
                let curr_portal_id = portal_ids[i];
                let prev_portal_id = portal_ids[i - 1];

                let p = self.portals[curr_portal_id];
                let p_prev = self.portals[prev_portal_id];

                let prev_mid_portal = p_prev.barycenter(&vertices);

                let mut left = p[0].get_vertex(&vertices);
                let mut right = p[1].get_vertex(&vertices);
                if is_left_to(&prev_mid_portal, right, left) {
                    std::mem::swap(&mut right, &mut left);
                }
                portals.push(Portal {
                    left,
                    right
                });
            }

            // append the last portal vertex
            /*let p = self.portals[portal_ids[portal_ids.len() - 1]];
            let mid_portal = p.barycenter(&vertices);
            let path_dir = end - mid_portal;

            let mut left = p[0].get_vertex(&vertices);
            let mut right = p[1].get_vertex(&vertices);
            if path_dir.det(&(*right - mid_portal)) < 0.0 {
                std::mem::swap(&mut right, &mut left);
            }
            portals.push(Portal {
                left,
                right
            });*/

            Some(self.apply_funnel(&start, &end, portals, vertices))
        })
    }

    fn apply_funnel<'a>(&self, mut apex: &'a Point2<f32>, end: &'a Point2<f32>, mut portals: Vec<Portal<'a>>, vertices: &[Point2<f32>]) -> Vec<Point2<f32>> {
        portals.push(Portal { left: end, right: end });
        
        let mut apex = *apex;
        let mut path = vec![apex];

        let mut left = *portals[0].left;
        let mut right = *portals[0].right;

        let mut left_portal_id = 0;
        let mut right_portal_id = 0;

        let mut portal_id = 1;
        while portal_id < portals.len() {
            let Portal { left: new_left, right: new_right } = portals[portal_id];

           if ((apex != *new_left && left != *new_left && is_left_to(&apex, &left, new_left)) || left == apex) {
                // tighten the funnel by the left
                left = *new_left;
                left_portal_id = portal_id;

                if is_left_to(&apex, &right, &left) {
                    path.push(right);
                    apex = right;

                    left = apex;
                    left_portal_id = right_portal_id;
                    portal_id = right_portal_id;

                    continue;
                }
            }

            // v is a 'right' vertex
            if ((apex != *new_right && right != *new_right && is_left_to(&apex, new_right, &right)) || right == apex) {
                // tighten the funnel by the right
                right = *new_right;
                right_portal_id = portal_id;

                if is_left_to(&apex, &right, &left) {
                    path.push(left);
                    apex = left;

                    right = apex;
                    right_portal_id = left_portal_id;
                    portal_id = left_portal_id;

                    continue;
                }
            }

            portal_id += 1;
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

        let start = Point2 { x: rand::random::<f32>(), y: rand::random::<f32>() };
        let end = Point2 { x: rand::random::<f32>(), y: rand::random::<f32>() };

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

    use crate::nav_mesh::Portal;
    #[test]
    fn test_navmesh_portals() {
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

        //let start = Point2 { x: 0.2, y: 0.2 };
        //let end = Point2 { x: 0.4, y: 0.8 };

        //let start = Point2 { x: 0.17776436, y: 0.08432096 };
        //let end =  Point2 { x: 0.56069654, y: 0.43880218 };

        let start = Point2 { x: rand::random::<f32>(), y: rand::random::<f32>() };
        let end = Point2 { x: rand::random::<f32>(), y: rand::random::<f32>() };

        let portals = nav_mesh
            .find_path_portals(start, end, &vertices)
            .expect("no path found");

        let path_portals = nav_mesh
            .find_path_through_portals(start, end, &vertices)
            .expect("no path found");

        let path = nav_mesh
            .find_path(start, end, &vertices)
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

        for (i, &idx1) in path_portals.iter().enumerate() {
            let [u, v] = nav_mesh.portals[idx1];

            let p1 = u.get_vertex(&vertices);
            let p2 = v.get_vertex(&vertices);

            if i == 0 {
                draw_line_segment_mut(
                    &mut img,
                    (p1.x * w, p1.y * h), // start point
                    (p2.x * w, p2.y * h), // end point
                    Rgb([0u8, 255u8, 0u8]),                             // RGB colors
                );
            } else {
                draw_line_segment_mut(
                    &mut img,
                    (p1.x * w, p1.y * h), // start point
                    (p2.x * w, p2.y * h), // end point
                    Rgb([0u8, 255u8, 255u8]),                             // RGB colors
                );
            }
        }

        for Portal { left, right } in portals.iter() {
            draw_cross_mut(
                &mut img,
                Rgb([255u8, 20u8, 13u8]),
                (left.x * w) as i32,              // start point
                (left.y * h) as i32,            // end point
            );
            draw_cross_mut(
                &mut img,
                Rgb([25u8, 20u8, 255u8]),
                (right.x * w) as i32,              // start point
                (right.y * h) as i32,            // end point
            );
        }

        for (p1, p2) in path.iter().zip(path.iter().skip(1)) {
            draw_line_segment_mut(
                &mut img,
                (p1.x * w, p1.y * h), // start point
                (p2.x * w, p2.y * h), // end point
                Rgb([255u8, 255u8, 0u8]),                             // RGB colors
            );
        }

        img.save("nav_mesh_portals.png").unwrap();
    }
}
