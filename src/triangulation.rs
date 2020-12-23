use super::{VertexIdx, Edge};
use crate::coord::Point2;

type TriangleIdx = usize;

use std::collections::{HashMap, HashSet};
#[derive(Debug)]
pub struct DelaunayTriangulation {
    vertices: HashMap<(VertexIdx, VertexIdx), VertexIdx>,
    triangles: HashSet<(VertexIdx, VertexIdx)>,
}

/// Compute the intersection given by the two segments
/// [v1; v2] and [v3; v4]
pub fn intersection(v1: &Point2, v2: &Point2,  v3: &Point2, v4: &Point2) -> Option<Point2> {
    let r = v2 - v1;
    let s = v4 - v3;

    let denom = det(&r, &s);
    let v3_minus_v1 = v3 - v1;
    let num = det(&v3_minus_v1, &r);

    // the segments are colinear
    if num == 0.0 && denom == 0.0 {
        let rr = r.dot(&r);
        let mut t0 = v3_minus_v1.dot(&r) / rr;
        let mut t1 = t0 + s.dot(&r) / rr;

        if s.dot(&r) < 0.0 {
            std::mem::swap(&mut t0, &mut t1);
        }

        let is_overlapping = (t1 >= 0.0 && t0 <= 1.0);
        if is_overlapping {
            // Give one point 
            if t0 >= 0.0 {
                Some(v1 + r*t0)
            } else {
                Some(v1 + r*t1)
            }
        } else {
            None
        }
    } else if denom == 0.0 && num != 0.0 {
        // the segments are parallel and not intersecting
        None
    } else if num != 0.0 {
        let u = num / denom;
        let t = det(&v3_minus_v1, &s) / denom;

        // the segments are not parallel and intersecting
        if u >= 0.0 && u <= 1.0 && t >= 0.0 && t <= 1.0 {
            Some(v1 + r*t)
        } else {
            None
        }
    } else {
        // The segments are not parallel and not intersecting
        None
    }
}

fn det(u: &Point2, v: &Point2) -> f32 {
    u.x * v.y - u.y * v.x
}

fn barycenter(u: &Point2, v: &Point2, w: &Point2) -> Point2 {
    (u + v + w)/3.0
}

fn edge_intersect(u: VertexIdx, v: VertexIdx, w: VertexIdx, x: VertexIdx, vertices: &[Point2], super_vertices: &[Point2]) -> (VertexIdx, VertexIdx) {
    // Get the ending vertex
    let up = u.get_vertex(vertices, super_vertices);
    // Get the triangle vertices
    let vp = v.get_vertex(vertices, super_vertices);
    let wp = w.get_vertex(vertices, super_vertices);
    let xp = x.get_vertex(vertices, super_vertices);
    let b = barycenter(vp, wp, xp);

    if let Some(_) = intersection(&b, up, vp, wp) {
        (v, w)
    } else if let Some(_) = intersection(&b, up, wp, xp) {
        (w, x)
    } else {
        assert!(intersection(&b, up, xp, vp).is_some());
        // If it does not intersect (v, w) nor (w, x)
        // it has to intersect (x, v)
        (x, v)
    }
}

impl DelaunayTriangulation {
    pub fn new() -> Self {
        let vertices = HashMap::new();
        let triangles = HashSet::new();

        let mut c = Self {
            vertices,
            triangles,
        };

        // Insert the first super triangle
        c.add_triangle(VertexIdx::Super(0), VertexIdx::Super(1), VertexIdx::Super(2));
        c
    }

    /// u is the vertex to insert in the triangulation
    /// This methods walks in the triangulation to find
    /// one triangle whose circumcircle encloses u
    pub fn get_triangle_whose_circle_encloses_u(&self, u: VertexIdx, vertices: &[Point2], super_vertices: &[Point2]) -> Option<(VertexIdx, VertexIdx, VertexIdx)> {
        if let Some(((mut v, mut w), x)) = self.vertices.iter().next() {
            let mut x = *x;
            // First triangle (vwx) is positively defined
            let mut outside = false;
            while !in_circumcircle(u, v, w, x, vertices, super_vertices) && !outside {
                let (a, b) = edge_intersect(u, v, w, x, vertices, super_vertices);
                // Get the adjacent triangle of swap(a, b) = (b, a)
                if let Some(c) = self.adjacent(b, a) {
                    v = b;
                    w = a;
                    x = c;
                } else {
                    // We go out of the triangulation!
                    outside = true;
                }
            }

            if outside {
                None
            } else {
                Some((v, w, x))
            }
        } else {
            // Empty triangulation
            None
        }
    }

    /// Insert the vertex u in the triangulation
    /// given a positively oriented triangle vwx whose
    /// circumcircle encloses u
    pub fn insert_vertex(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx, x: VertexIdx, vertices: &[Point2], super_vertices: &[Point2]) {
        self.delete_triangle(v, w, x);
        self.dig_cavity(u, v, w, vertices, super_vertices);
        self.dig_cavity(u, w, x, vertices, super_vertices);
        self.dig_cavity(u, x, v, vertices, super_vertices);
    }

    fn dig_cavity(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx, vertices: &[Point2], super_vertices: &[Point2]) {
        if let Some(x) = self.adjacent(w, v) {
            if in_circumcircle(u, w, v, x, vertices, super_vertices) {
                self.delete_triangle(w, v, x);
                self.dig_cavity(u, v, x, vertices, super_vertices);
                self.dig_cavity(u, x, w, vertices, super_vertices);
            } else {
                self.add_triangle(u, v, w);
            }
        } else {
            // TEST that!!!!
            self.add_triangle(u, v, w);
        }
    }

    fn add_triangle(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx) {
        // Reject triangles containing an edge given in the same order
        if self.vertices.contains_key(&(u, v)) || self.vertices.contains_key(&(v, w)) || self.vertices.contains_key(&(w, u)) {
            return;
        }

        // Do not add the triangles at the border of the triangulation
        // i.e. those containing super vertices.
        match (u, v, w) {
            (VertexIdx::Vertices(_), VertexIdx::Vertices(_), VertexIdx::Vertices(_)) => {
                self.triangles.insert((u, v));
            },
            _ => ()
        }
    
        self.vertices.insert((u, v), w);
        self.vertices.insert((v, w), u);
        self.vertices.insert((w, u), v);
    }

    fn delete_triangle(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx) {
        self.triangles.remove(&(u, v)); 
        self.triangles.remove(&(v, w));
        self.triangles.remove(&(w, u));

        self.vertices.remove(&(u, v));
        self.vertices.remove(&(v, w));
        self.vertices.remove(&(w, u));
    }

    fn adjacent(&self, u: VertexIdx, v: VertexIdx) -> Option<VertexIdx> {
        self.vertices.get(&(u, v)).cloned()
    }
}

impl IntoIterator for DelaunayTriangulation {
    type Item = [usize; 3];
    type IntoIter = TriangleIntoIterator;

    fn into_iter(self) -> Self::IntoIter {
        TriangleIntoIterator {
            triangles: self.triangles.into_iter(),
            vertices: self.vertices,
        }
    }
}

pub struct TriangleIntoIterator {
    triangles: std::collections::hash_set::IntoIter<(VertexIdx, VertexIdx)>,
    vertices: HashMap<(VertexIdx, VertexIdx), VertexIdx>,
}
impl Iterator for TriangleIntoIterator {
    type Item = [usize; 3];

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((u, v)) = self.triangles.next() {
            let w = self.vertices.get(&(u, v)).cloned().unwrap();

            let (u, v, w) = match (u, v, w) {
                (VertexIdx::Vertices(u), VertexIdx::Vertices(v), VertexIdx::Vertices(w)) => {
                    (u, v, w)
                },
                _ => unreachable!()
            };
            Some([u, v, w])
        } else {
            // no triangles left
            None
        }
    }
}

/// u is the vertex to test
/// v, w, x defines a positively oriented triangle
fn in_circumcircle(u: VertexIdx, v: VertexIdx, w: VertexIdx, x: VertexIdx, vertices: &[Point2], super_vertices: &[Point2]) -> bool {
    let u = u.get_vertex(vertices, super_vertices);
    let v = v.get_vertex(vertices, super_vertices);
    let w = w.get_vertex(vertices, super_vertices);
    let x = x.get_vertex(vertices, super_vertices);

    // p is inside the triangle defined by (a, b, c) (given in counter-clockwise order) if:
    //       | ax-px, ay-py, (ax-px)² + (ay-py)² |
    // det = | bx-px, by-py, (bx-px)² + (by-py)² | > 0.0
    //       | cx-px, cy-py, (cx-px)² + (cy-py)² |
    let a1 = v.x - u.x;
    let b1 = w.x - u.x;
    let c1 = x.x - u.x;

    let a2 = v.y - u.y;
    let b2 = w.y - u.y;
    let c2 = x.y - u.y;

    let a3 = a1*a1 + a2*a2;
    let b3 = b1*b1 + b2*b2;
    let c3 = c1*c1 + c2*c2;

    a1*b2*c3 + a2*b3*c1 + b1*c2*a3 - c1*b2*a3 - c2*b3*a1 - b1*a2*c3 > 0.0
}

pub fn triangulate2(vertices: &[na::Point2<f32>]) -> DelaunayTriangulation {
    let super_vertices = &[
        Point2::new(-3.0, -1.0),
        Point2::new(3.0, -1.0),
        Point2::new(0.0, 3.0),
    ];

    let vertices = unsafe {
        std::slice::from_raw_parts(vertices.as_ptr() as *const Point2, vertices.len())
    };
    let mut c = DelaunayTriangulation::new();

    for (idx_vertex, vertex) in vertices.iter().enumerate() {
        let u = VertexIdx::Vertices(idx_vertex);
        //c = dbg!(c);

        if let Some((v, w, x)) = c.get_triangle_whose_circle_encloses_u(u, vertices, super_vertices) {
            c.insert_vertex(u, v, w, x, vertices, super_vertices);
        } else {
            panic!("the vertex {:?} is outside the triangulation", vertex);
        }
    }

    c
}
#[cfg(test)]
mod tests {
    use super::{DelaunayTriangulation, triangulate2};
    use crate::coord::Point2;
    use crate::VertexIdx;

    /*#[test]
    fn test_triangulate() {
        let vertices: &[na::Point2<f32>] = &[
            na::Point2::new(0.3, 0.1),
            na::Point2::new(0.5, 0.4),
        ];
        
        for t in triangulate2(vertices).into_iter() {
            println!("{:?}", t);
        }

        //assert_eq!(connexity.triangles.len(), 3);
    }*/

    #[test]
    fn test_triangulate_3_points() {
        println!("aaaaa");
        let vertices: &[na::Point2<f32>] = &[
            na::Point2::new(0.76809347, 0.17880994),
            na::Point2::new(0.14206064, 0.8896956),
            na::Point2::new(0.007408619, 0.17449331),
        ];
        for t in triangulate2(&vertices).into_iter() {
            println!("{:?}", t);
        }
    }
    use image::{Rgb, RgbImage};
    use imageproc::drawing::draw_line_segment_mut;
    #[test]
    fn test_triangulate_complex() {
        let num_vertices = 1000;
        let vertices = (0..num_vertices)
            .map(|_| na::Point2::new(rand::random::<f32>(), rand::random::<f32>()))
            .collect::<Vec<_>>();
        
        /*for t in triangulate2(&vertices).into_iter() {
            println!("{:?}", t);
        }*/

        let triangulation = triangulate2(&vertices);
        let (w, h) = (512.0, 512.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation {
            for (&idx1, &idx2) in t.iter().zip(t.iter().skip(1).cycle()) {
                //let v1 = idx.get_vertex(super_triangle)

                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h),              // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h),            // end point
                    Rgb([69u8, 203u8, 133u8]), // RGB colors
                );
            }
        }

        img.save("delaunay10.png").unwrap();
        //assert_eq!(connexity.triangles.len(), 3);
    }
}