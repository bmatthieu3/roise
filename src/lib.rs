// externing crate for test-only use
#[cfg(test)]
extern crate image;
#[cfg(test)]
extern crate imageproc;
#[cfg(test)]
extern crate rand;

pub mod graph;

pub mod sampling;
pub mod geometry;
pub mod triangulation;
pub mod noise;
pub mod nav_mesh;

pub use crate::geometry::coord::Point2;

fn lies_on_positive_half_plane(p: &Point2<f32>, a: &Point2<f32>, b: &Point2<f32>) -> bool {
    let ap = p - a;
    let ab = b - a;
    ab.det(&ap) >= 0.0
}

#[derive(Clone, Copy)]
#[derive(PartialEq, Eq, Hash)]
#[derive(Debug)]
pub enum VertexIdx {
    Super(usize),
    Vertices(usize)
}


impl PartialOrd for VertexIdx {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (VertexIdx::Super(i), VertexIdx::Super(j)) => i.partial_cmp(j),
            (VertexIdx::Super(_), VertexIdx::Vertices(_)) => Some(std::cmp::Ordering::Greater),
            (VertexIdx::Vertices(_), VertexIdx::Super(_)) => Some(std::cmp::Ordering::Less),
            (VertexIdx::Vertices(i), VertexIdx::Vertices(j)) => i.partial_cmp(j),
        }
    }
}

impl VertexIdx {
    fn get_vertex<'a>(&self, vertices: &'a [Point2<f32>]) -> &'a Point2<f32> {
        match self {
            VertexIdx::Super(idx) => &crate::triangulation::SUPER_VERTICES[*idx],
            VertexIdx::Vertices(idx) => &vertices[*idx]
        }
    }

    fn id(&self) -> usize {
        match self {
            VertexIdx::Super(idx) => *idx,
            VertexIdx::Vertices(idx) => 3 + *idx
        }
    }
}

use std::convert::TryInto;
use std::collections::HashSet;
pub fn triangulate(vertices: &[Point2<f32>]) -> Box<[[usize; 3]]> {
    let mut triangulation = Vec::new();
    triangulation.push(
        Triangle([VertexIdx::Super(0), VertexIdx::Super(1), VertexIdx::Super(2)])
    );
    let mut bad_tri = vec![];

    //let mut connexity = Connexity::default();

    for (idx_vertex, p) in vertices.iter().enumerate() {
        bad_tri.clear();

        triangulation.retain(|t| {
            if t.in_circumcircle(p, vertices) {
                bad_tri.push(t.clone());
                false
            } else {
                true
            }
        });

        let poly_edges = get_star_shaped_polygon(&bad_tri);

        // Add triangles
        for e in poly_edges.into_iter() {
            let a = e.start().get_vertex(vertices);
            let b = e.end().get_vertex(vertices);

            let ab = b - a;
            let ap = p - a;

            let d = ab.det(&ap);
            let new_tri = if d > 0.0 {
                // ab is well-defined (counter clockwise for the star shaped polygon)
                Triangle([
                    e.start(),
                    e.end(),
                    VertexIdx::Vertices(idx_vertex),
                ])
            } else {
                // ab is not well-defined (clockwise for the star shaped polygon)
                Triangle([
                    e.end(),
                    e.start(),
                    VertexIdx::Vertices(idx_vertex),
                ])
            };

            triangulation.push(new_tri);
        }
    }

    let result = triangulation.iter()
        .filter_map(|t| {
            let mut indices = vec![];
            for idx in &t.0 {
                match idx {
                    VertexIdx::Vertices(i) => {
                        indices.push(*i)
                    },
                    _ => {
                        return None;
                    }
                }
            }
            let indices = indices.into_boxed_slice();
            let boxed_array: Box<[usize; 3]> = indices.try_into().unwrap();
            Some(*boxed_array)
        })
        .collect::<Vec<_>>();

    result.into_boxed_slice()
}

fn get_star_shaped_polygon(triangles: &[Triangle]) -> HashSet<Edge> {
    let num_triangles = triangles.len();

    if num_triangles == 1 {
        triangles[0].edge_iter().collect()
    } else {
        let mut edges = triangles.iter()
            .flat_map(|t| t.edge_iter().collect::<Vec<_>>())
            .collect::<HashSet<_>>();

        for i in 1..num_triangles {
            let t1 = &triangles[i];
            for t2 in triangles.iter().take(i) {
                for (e1, e2) in t1.edge_iter().zip(t2.edge_iter()) {
                    if t1.contains_edge(&e2) {
                        edges.remove(&e2);
                    }

                    if t2.contains_edge(&e1) {
                        edges.remove(&e1);
                    }
                }
            }
        }
        edges
    }
}

pub trait Shape {
    /// Triangle vertices are given in counter-clockwise order
    fn contains(&self, p: &Point2<f32>, vertices: &[Point2<f32>]) -> bool;

    /// Triangle vertices are given in counter-clockwise order
    fn in_circumcircle(&self, p: &Point2<f32>, vertices: &[Point2<f32>]) -> bool;

    /// Check whether a vertex belongs to the shape
    fn contains_vertex(&self, idx: VertexIdx) -> bool;

    /// Check whether an edge belongs to the shape
    fn contains_edge(&self, e: &Edge) -> bool;

    /// An iterator of the edges of the shape
    fn edge_iter(&self) -> EdgeIterator;
}

#[derive(Clone)]
#[derive(PartialEq, Eq, Hash)]
#[derive(Debug)]
struct Triangle([VertexIdx; 3]);

impl Triangle {
    fn get_vertices<'a>(&self, vertices: &'a [Point2<f32>]) -> [&'a Point2<f32>; 3] {
        [
            self.0[0].get_vertex(vertices),
            self.0[1].get_vertex(vertices),
            self.0[2].get_vertex(vertices),
        ]
    }
}

impl Shape for Triangle {
    fn contains(&self, p: &Point2<f32>, vertices: &[Point2<f32>]) -> bool {
        let vertices = self.get_vertices(vertices);

        let pos_e1 = lies_on_positive_half_plane(p, vertices[0], vertices[1]);
        let pos_e2 = lies_on_positive_half_plane(p, vertices[1], vertices[2]);
        if pos_e1 != pos_e2 {
            false
        } else {
            let pos_e3 = lies_on_positive_half_plane(p, vertices[2], vertices[0]);

            pos_e1 == pos_e3
        }
    }

    /// Triangle vertices are given in counter-clockwise order
    fn in_circumcircle(&self, p: &Point2<f32>, vertices: &[Point2<f32>]) -> bool {
        let vertices = self.get_vertices(vertices);
        
        // p is inside the triangle defined by (a, b, c) (given in counter-clockwise order) if:
        //       | ax-px, ay-py, (ax-px)² + (ay-py)² |
        // det = | bx-px, by-py, (bx-px)² + (by-py)² | > 0.0
        //       | cx-px, cy-py, (cx-px)² + (cy-py)² |
        let a1 = vertices[0].x - p.x;
        let b1 = vertices[1].x - p.x;
        let c1 = vertices[2].x - p.x;

        let a2 = vertices[0].y - p.y;
        let b2 = vertices[1].y - p.y;
        let c2 = vertices[2].y - p.y;

        let a3 = a1*a1 + a2*a2;
        let b3 = b1*b1 + b2*b2;
        let c3 = c1*c1 + c2*c2;

        a1*b2*c3 + a2*b3*c1 + b1*c2*a3 - c1*b2*a3 - c2*b3*a1 - b1*a2*c3 > 0.0
    }

    fn contains_vertex(&self, idx: VertexIdx) -> bool {
        self.0[0] == idx || self.0[1] == idx || self.0[2] == idx
    }

    fn contains_edge(&self, e: &Edge) -> bool {
        self.contains_vertex(e.start()) && self.contains_vertex(e.end())
    }

    fn edge_iter(&self) -> EdgeIterator {
        EdgeIterator {
            edges: &self.0,
            cur_idx: 0
        }
    }
}

#[derive(Debug)]
#[derive(Eq, Clone)]
pub struct Edge(VertexIdx, VertexIdx);

impl Edge {
    fn start(&self) -> VertexIdx {
        self.0
    }

    fn end(&self) -> VertexIdx {
        self.1
    }
}

impl std::hash::Hash for Edge {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
        self.1.hash(state);
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.0 == other.0 && self.1 == other.1) ||
        (self.0 == other.1 && self.1 == other.0)
    }
}

pub struct EdgeIterator<'a> {
    edges: &'a [VertexIdx],
    cur_idx: usize,
}
impl Iterator for EdgeIterator<'_> {
    type Item = Edge;

    fn next(&mut self) -> Option<Edge> {
        let num_edges = self.edges.len();
        if self.cur_idx >= num_edges {
            None
        } else {
            let cur = self.cur_idx;
            let next = (self.cur_idx + 1) % num_edges;

            let edge = Edge(self.edges[cur], self.edges[next]);
            self.cur_idx += 1;

            Some(edge)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{lies_on_positive_half_plane, Triangle, Shape, VertexIdx, triangulate};
    use crate::geometry::coord::Point2;
    #[test]
    fn half_plane_position() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.0, 1.0);

        assert!(!lies_on_positive_half_plane(&Point2::new(10.0, -1.0), &a, &b));
        assert!(!lies_on_positive_half_plane(&Point2::new(-10.0, -1.0), &a, &b));
        assert!(lies_on_positive_half_plane(&Point2::new(-10.0, 1.0), &a, &b));

        assert!(!lies_on_positive_half_plane(&Point2::new(10.0, -1.0), &b, &c));
        assert!(lies_on_positive_half_plane(&Point2::new(10.0, -1.0), &c, &a));
    }

    #[test]
    fn test_contains_triangle() {
        let vertices = [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0)
        ];
        let t = Triangle([VertexIdx::Vertices(0), VertexIdx::Vertices(1), VertexIdx::Vertices(2)]);

        assert!(!t.contains(&Point2::new(10.0, -1.0), &vertices));
        assert!(!t.contains(&Point2::new(-1e-5, 0.5), &vertices));
        assert!(t.contains(&Point2::new(0.5, 0.5), &vertices));
        assert!(t.contains(&Point2::new(0.25, 0.25), &vertices));
    }

    #[test]
    fn test_triangulation() {
        let vertices = [
            Point2::new(0.25, 0.25),
            Point2::new(0.75, 0.25),
            Point2::new(0.25, 0.75),
            Point2::new(0.75, 0.75),
        ];

        let triangulation = triangulate(&vertices);
        assert!(triangulation.len() == 2);
    }

    
    #[test]
    fn test_triangulation2() {
        let vertices = [
            Point2::new(0.25, 0.25),
            Point2::new(0.75, 0.25),
            Point2::new(0.25, 0.75),
            Point2::new(0.75, 0.75),
            Point2::new(0.6, 0.3),
        ];

        let triangulation = triangulate(&vertices);
        println!("{:?}", triangulation);
        assert!(triangulation.len() == 4);
    }

    use image::{Rgb, RgbImage};
    use imageproc::drawing::draw_line_segment_mut;
    #[test]
    fn test_triangulation_image() {
        let num_vertices = 1000;
        let vertices = (0..num_vertices)
            .map(|_| Point2::new(rand::random::<f32>(), rand::random::<f32>()))
            .collect::<Vec<_>>();

        let triangulation = triangulate(&vertices);
        let (w, h) = (512.0, 512.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation.iter() {
            for (&idx1, &idx2) in t.iter().zip(t.iter().skip(1).cycle()) {
                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h),              // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h),            // end point
                    Rgb([69u8, 203u8, 133u8]), // RGB colors
                );
            }
        }

        img.save("delaunay1000.png").unwrap();
    }
}
