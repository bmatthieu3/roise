use crate::{graph, Point2};
use std::collections::HashMap;
use crate::VertexIdx;
use crate::graph::Graph;
use crate::triangulation::DelaunayTriangulation;



pub struct NavMesh {
    /// portal edges
    pub portals: Vec<[VertexIdx; 2]>,
    /// the graph containing the adjency map and allowing a* algorithm
    pub graph: Graph,
}

impl NavMesh {
    pub fn from_triangulation(triangulation: &DelaunayTriangulation) -> NavMesh {
        // Build an edge map
        let mut edges: HashMap<(VertexIdx, VertexIdx), usize> = HashMap::new();
        let mut edge_idx = 0;
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut portals = vec![];

        for ((mut u, mut v), w) in triangulation.vertices.iter() {
            // always treat a edge 1 time
            if u.id() < v.id() && edges.get(&(u, v)).is_none() {
                edges.insert((u, v), edge_idx);
                portals.push([u, v]);
                edge_idx += 1;
            } else if u.id() > v.id() && edges.get(&(v, u)).is_none() {
                edges.insert((v, u), edge_idx);
                portals.push([v, u]);
                edge_idx += 1;
            }
        }

        //let edges = dbg!(edges);

        for ((u, v), edge_idx) in edges.iter() {
            match (triangulation.adjacent(*u, *v), triangulation.adjacent(*v, *u)) {
                (Some(a), Some(b)) => {
                    let e1 = if a.id() < u.id() {
                        (a, *u)
                    } else {
                        (*u, a)
                    };

                    let e2 = if a.id() < v.id() {
                        (a, *v)
                    } else {
                        (*v, a)
                    };

                     let e3 = if b.id() < u.id() {
                        (b, *u)
                    } else {
                        (*u, b)
                    };

                    let e4 = if b.id() < v.id() {
                        (b, *v)
                    } else {
                        (*v, b)
                    };

                    //dbg!(a, b, e1, edges.get(&e1), e2, edges.get(&e2), e3, edges.get(&e3), e4, edges.get(&e4));

                    let neigh_edge_idx1 = *edges.get(&e1).unwrap();
                    let neigh_edge_idx2 = *edges.get(&e2).unwrap();
                    let neigh_edge_idx3 = *edges.get(&e3).unwrap();
                    let neigh_edge_idx4 = *edges.get(&e4).unwrap();

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
            portals,
            graph
        }
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

        let NavMesh { portals, graph } = NavMesh::from_triangulation(&triangulation);

        let path = graph
            .find_path(100, 600, &portals, &vertices)
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

        /*for (t_idx, neigh_indices) in adj.iter() {
            for n_idx in neigh_indices {
                let neigh_neigh_indices = adj.get(n_idx).unwrap();
                assert!(neigh_neigh_indices.contains(t_idx));
                draw_line_segment_mut(
                    &mut img,
                    (barycenters[*t_idx].x * w, barycenters[*t_idx].y * h), // start point
                    (barycenters[*n_idx].x * w, barycenters[*n_idx].y * h), // end point
                    Rgb([255u8, 0u8, 255u8]),                               // RGB colors
                );
            }
        }*/

        for (&idx1, &idx2) in path.iter().zip(path.iter().skip(1)) {
            let b1 = portals[idx1].barycenter(&vertices);
            let b2 = portals[idx2].barycenter(&vertices);
            draw_line_segment_mut(
                &mut img,
                (b1.x * w, b1.y * h), // start point
                (b2.x * w, b2.y * h), // end point
                Rgb([255u8, 0u8, 0u8]),                             // RGB colors
            );
        }

        img.save("nav_mesh_portals.png").unwrap();
    }
}
