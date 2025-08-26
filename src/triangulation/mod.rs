use super::{VertexIdx, Edge};
use crate::geometry::coord::Point2;

mod marching_square;
mod sweep_line_triangulation;

use crate::geometry::closed_polyline::ClosedPolyline;

type TriangleIdx = usize;

use std::collections::{HashMap, HashSet};
#[derive(Debug)]
pub struct DelaunayTriangulation {
    vertices: HashMap<(VertexIdx, VertexIdx), VertexIdx>,
    triangles: HashSet<(VertexIdx, VertexIdx, VertexIdx)>,
}

/// Compute the intersection given by the two segments
/// [v1; v2] and [v3; v4]
pub fn intersection(v1: &Point2<f32>, v2: &Point2<f32>,  v3: &Point2<f32>, v4: &Point2<f32>) -> Option<Point2<f32>> {
    let r = v2 - v1;
    let s = v4 - v3;

    let denom = (&r).det(&s);
    let v3_minus_v1 = v3 - v1;
    let num = v3_minus_v1.det(&r);
    let eps = 1e-9;

    // the segments are colinear
    if denom.abs() < eps {
        if num.abs() < eps {
            // colinear
            let rr = r.dot(&r);
            if rr.abs() < eps {
                // v1==v2 degenerate
                return None;
            }
            let mut t0 = v3_minus_v1.dot(&r) / rr;
            let mut t1 = t0 + s.dot(&r) / rr;
            if t0 > t1 {
                std::mem::swap(&mut t0, &mut t1);
            }

            if t1 < 0.0 || t0 > 1.0 {
                None
            } else {
                Some(v1 + r * t0.clamp(0.0, 1.0))  // or return (segment)
            }
        } else {
            None
        }
    } else {
        let t = v3_minus_v1.det(&s) / denom;
        let u = num / denom;

        if 0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0 {
            Some(v1 + r * t)
        } else {
            None
        }
    }

    /*if num.abs() < eps && denom.abs() < eps {
        let rr = r.dot(&r);
        if rr.abs() < eps {
            return None;
        }

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
    } else if denom.abs() < eps && num != 0.0 {
        // the segments are parallel and not intersecting
        None
    } else if num != 0.0 {
        let u = num / denom;
        let t = (&v3_minus_v1).det(&s) / denom;

        // the segments are not parallel and intersecting
        if u >= 0.0 && u <= 1.0 && t >= 0.0 && t <= 1.0 {
            Some(v1 + r*t)
        } else {
            None
        }
    } else {
        // The segments are not parallel and not intersecting
        None
    }*/
}

pub fn barycenter(u: &Point2<f32>, v: &Point2<f32>, w: &Point2<f32>) -> Point2<f32> {
    (u + v + w)/3.0
}

fn edge_intersect(u: VertexIdx, v: VertexIdx, w: VertexIdx, x: VertexIdx, vertices: &[Point2<f32>]) -> (VertexIdx, VertexIdx) {
    // Get the ending vertex
    let up = u.get_vertex(vertices);
    // Get the triangle vertices
    let vp = v.get_vertex(vertices);
    let wp = w.get_vertex(vertices);
    let xp = x.get_vertex(vertices);
    let b = barycenter(vp, wp, xp);

    if let Some(_) = intersection(up, &b, vp, wp) {
        (v, w)
    } else if let Some(_) = intersection(up,  &b, wp, xp) {
        (w, x)
    } else {
        //assert!(intersection(&b, up, xp, vp).is_some());
        // If it does not intersect (v, w) nor (w, x)
        // it has to intersect (x, v)
        (x, v)
    }
}

fn is_triangle_convex(a: VertexIdx, b: VertexIdx, c: VertexIdx, vertices: &[Point2<f32>]) -> bool {
    let pa = a.get_vertex(vertices);
    let pb = b.get_vertex(vertices);
    let pc = c.get_vertex(vertices);

    (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x) > 0.0
}

fn point_in_triangle(p: VertexIdx, a: VertexIdx, b: VertexIdx, c: VertexIdx, vertices: &[Point2<f32>]) -> bool {
    let o1 = is_triangle_convex(a, b, p, vertices);
    let o2 = is_triangle_convex(b, c, p, vertices);
    let o3 = is_triangle_convex(c, a, p, vertices);

    let has_pos = o1 || o2 || o3;
    let has_neg = !o1 || !o2 || !o3;

    !(has_pos && has_neg)
}

fn is_ear(prev: VertexIdx, curr: VertexIdx, next: VertexIdx, polygon: &[VertexIdx], vertices: &[Point2<f32>]) -> bool {
    if !is_triangle_convex(prev, curr, next, vertices) {
        return false;
    }

    for p in polygon {
        if *p != prev && *p != curr && *p != next {
            if point_in_triangle(*p, prev, curr, next, vertices) {
                return false;
            }
        }
    }
        
    true
}

pub(crate) const SUPER_VERTICES: &[Point2<f32>] = &[
    Point2::new(-3.0, -1.0),
    Point2::new(3.0, -1.0),
    Point2::new(0.0, 3.0),
];

impl DelaunayTriangulation {
    pub fn from_vertices(vertices: &[Point2<f32>]) -> Self {
        let mut triangulation = Self {
            vertices: HashMap::new(),
            triangles: HashSet::new(),
        };

        // Insert the first super triangle
        triangulation.add_triangle(VertexIdx::Super(0), VertexIdx::Super(1), VertexIdx::Super(2));

        for (idx_vertex, vertex) in vertices.iter().enumerate() {
            let u = VertexIdx::Vertices(idx_vertex);

            if let Some((v, w, x)) = triangulation.get_triangle_whose_circle_encloses_u(u, vertices) {
                triangulation.insert_vertex(u, v, w, x, vertices);
            } else {
                panic!("the vertex {:?} is outside the triangulation", vertex);
            }
        }

        triangulation
    }

    pub fn from_contours(contours: &[ClosedPolyline]) -> Self {
        let mut num_vertices_per_contour = vec![];

        let vertices = contours
            .into_iter()
            .cloned()
            .flat_map(|ClosedPolyline { mut vertices }| {
                let _ = vertices.pop();

                num_vertices_per_contour.push(vertices.len());
                vertices
            })
            .collect::<Vec<_>>();

        let mut triangulation = DelaunayTriangulation::from_vertices(&vertices);

        let mut off = 0;
        for num_vertices_in_contour in num_vertices_per_contour {
            let mut i = num_vertices_in_contour - 1;
            for j in 0..num_vertices_in_contour {
                triangulation.enforce_edge(VertexIdx::Vertices(off + i), VertexIdx::Vertices(off + j), &vertices);

                i = j;
            }

            off += num_vertices_in_contour;
        }

        let mut t_delete = vec![];
        for (u, v, w) in triangulation.triangles.iter() {
            //let w = triangulation.adjacent(u, v).unwrap();
            let b = barycenter(u.get_vertex(&vertices), v.get_vertex(&vertices), w.get_vertex(&vertices));

            let mut is_in_hole = !contours[0].contains(&b);

            for hole in &contours[1..] {
                if hole.contains(&b) {
                    is_in_hole = !is_in_hole;
                }
            }

            if is_in_hole {
                t_delete.push((*u, *v, *w));
            }
        }

        for (u, v, w) in t_delete {
            triangulation.delete_triangle(u, v, w);
        }

        triangulation
    }

    /// u is the vertex to insert in the triangulation
    /// This methods walks in the triangulation to find
    /// one triangle whose circumcircle encloses u
    pub fn get_triangle_whose_circle_encloses_u(&self, u: VertexIdx, vertices: &[Point2<f32>]) -> Option<(VertexIdx, VertexIdx, VertexIdx)> {
        if let Some(((mut v, mut w), x)) = self.vertices.iter().next() {
            let mut x = *x;
            // First triangle (vwx) is positively defined
            let mut outside = false;
            while !in_circumcircle(u, v, w, x, vertices) && !outside {
                let (a, b) = edge_intersect(u, v, w, x, vertices);
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
    pub fn insert_vertex(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx, x: VertexIdx, vertices: &[Point2<f32>]) {
        self.delete_triangle(v, w, x);
        self.dig_cavity(u, v, w, vertices);
        self.dig_cavity(u, w, x, vertices);
        self.dig_cavity(u, x, v, vertices);
    }

    fn dig_cavity(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx, vertices: &[Point2<f32>]) {
        if let Some(x) = self.adjacent(w, v) {
            if in_circumcircle(u, w, v, x, vertices) {
                self.delete_triangle(w, v, x);
                self.dig_cavity(u, v, x, vertices);
                self.dig_cavity(u, x, w, vertices);
            } else {
                self.add_triangle(u, v, w);
            }
        } else {
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
                self.triangles.insert((u, v, w));
                self.triangles.insert((v, w, u));
                self.triangles.insert((w, u, v));
            },
            _ => ()
        }
    
        self.vertices.insert((u, v), w);
        self.vertices.insert((v, w), u);
        self.vertices.insert((w, u), v);
    }

    pub fn delete_triangle(&mut self, u: VertexIdx, v: VertexIdx, w: VertexIdx) {
        self.triangles.remove(&(u, v, w));
        self.triangles.remove(&(v, w, u));
        self.triangles.remove(&(w, u, v));

        self.vertices.remove(&(u, v));
        self.vertices.remove(&(v, w));
        self.vertices.remove(&(w, u));
    }

    fn adjacent(&self, u: VertexIdx, v: VertexIdx) -> Option<VertexIdx> {
        self.vertices.get(&(u, v)).cloned()
    }

    // find a triangle containing u and directed towards v
    fn find_triangle_containing_u(&self, u: VertexIdx, v: VertexIdx, vertices: &[Point2<f32>]) -> Option<(VertexIdx, VertexIdx, VertexIdx)> {
        for ((a, b), c) in &self.vertices {

            if u == *a && intersection(u.get_vertex(vertices), v.get_vertex(vertices), b.get_vertex(vertices), c.get_vertex(vertices)).is_some() {
                return Some((*a, *b, *c));
            } else if u == *b && intersection(u.get_vertex(vertices), v.get_vertex(vertices), c.get_vertex(vertices), a.get_vertex(vertices)).is_some() {
                return Some((*b, *c, *a));
            } else if u == *c && intersection(u.get_vertex(vertices), v.get_vertex(vertices), a.get_vertex(vertices), b.get_vertex(vertices)).is_some() {
                return Some((*c, *a, *b));
            }
        }

        None
    }

    pub fn enforce_edge(&mut self, u: VertexIdx, v: VertexIdx, vertices: &[Point2<f32>]) {
        if let Some(_) = self.adjacent(u, v) {
            return;
        }

        if let Some(_) = self.adjacent(v, u) {
            return;
        }

        // Find a triangle containing u
        if let Some((mut a, mut b, mut c)) = self.find_triangle_containing_u(u, v, vertices) {
            // u == a
            // We know here that uv intersects (b, c)

            let mut triangles_intersecting= vec![(a, b)];

            let mut caveat_polyline = HashMap::new();

            /*let mut insert_triangle = |a: VertexIdx, b: VertexIdx, c: VertexIdx| {
                let e1 = (std::cmp::min(a, b), std::cmp::max(a, b));
                let e2 = (std::cmp::min(b, c), std::cmp::max(b, c));
                let e3 = (std::cmp::min(c, a), std::cmp::max(c, a));

                if let Some(p) = caveat_polyline.get_mut(&e1) {
                    *p += 1;
                } else {
                    caveat_polyline.insert(e1, 1);
                }
                if let Some(p) = caveat_polyline.get_mut(&e2) {
                    *p += 1;
                } else {
                    caveat_polyline.insert(e2, 1);
                }
                if let Some(p) = caveat_polyline.get_mut(&e3) {
                    *p += 1;
                } else {
                    caveat_polyline.insert(e3, 1);
                }
            };*/

            //insert_triangle(a, b, c);

            caveat_polyline.insert(u, b);
            caveat_polyline.insert(c, u);

            while b != v && c != v {
                if let Some(d) = self.adjacent(c, b) {
                    triangles_intersecting.push((c, b));
                    //insert_triangle(c, b, d);

                    if v == d {
                        caveat_polyline.insert(b, d);
                        caveat_polyline.insert(d, c);

                        break;
                    } else {
                        let r = (v.get_vertex(vertices) - u.get_vertex(vertices)).dot(&(d.get_vertex(vertices) - u.get_vertex(vertices)));
                        let is_colinear = r.abs() < 1e-9;

                        dbg!(is_colinear);
                        if is_colinear {
                            self.enforce_edge(d, v, vertices);
                        } else if intersection(u.get_vertex(vertices), v.get_vertex(vertices), b.get_vertex(vertices), d.get_vertex(vertices)).is_some() {
                            caveat_polyline.insert(d, c);
                            c = d;
                        } else if intersection(u.get_vertex(vertices), v.get_vertex(vertices), c.get_vertex(vertices), d.get_vertex(vertices)).is_some() {
                            caveat_polyline.insert(b, d);
                            b = d;
                        } else {
                            panic!("should intersect a edge!: {:?}", r);
                        }
                    }
                } else {
                    panic!("could not find the triangle containing v, border of triangulation");
                }
            }

            //caveat_polyline = caveat_polyline.into_iter().filter(|(k, v)| *v == 1).collect();
            
            dbg!(&caveat_polyline, u, v);

            for t in triangles_intersecting {
                self.delete_triangle(t.0, t.1, self.adjacent(t.0, t.1).unwrap());
            }


            // At this point caveat_polyline form an ordered polygon of the caveat
            // we will create 2 closed polylines in ccw and sharing (u, v)
            let mut pa = vec![u];
            let mut curr = u;
            while let Some(next) = caveat_polyline.get(&curr) {
                pa.push(*next);

                if *next == v {
                    break;
                }

                curr = *next;
            }

            let mut pb = vec![u, v];
            let mut curr = v;
            while let Some(next) = caveat_polyline.get(&curr) {
                if *next == u {
                    break;
                }

                pb.push(*next);

                curr = *next;
            }

            let ta = self.triangulate_polygon(pa, vertices);
            for (a, b, c) in ta {
                self.add_triangle(a, b, c);
            }

            let tb = self.triangulate_polygon(pb, vertices);
            for (a, b, c) in tb {
                self.add_triangle(a, b, c);
            }
        }
    }

    fn triangulate_polygon(&self, mut polyline: Vec<VertexIdx>, vertices: &[Point2<f32>]) -> Vec<(VertexIdx, VertexIdx, VertexIdx)> {
        let mut triangles = vec![];

        while polyline.len() > 3 {
            for i in 0..polyline.len() {
                let prev = polyline[(i-1) % polyline.len()];
                let curr = polyline[i % polyline.len()];
                let next = polyline[(i+1) % polyline.len()];

                // Check ear: no other point inside triangle
                if is_ear(prev, curr, next, &polyline, vertices) {
                    // This is an ear
                    triangles.push( (prev, curr, next) );
                    polyline.remove(i);
                    break;
                }
            }
        }
        // Add the last triangle
        triangles.push( (polyline[0], polyline[1], polyline[2]) );

        triangles
    }

    fn flip_edge(&mut self, u: VertexIdx, v: VertexIdx) {
        match (self.adjacent(u, v), self.adjacent(v, u)) {
            (Some(w), Some(x)) => {
                self.delete_triangle(u, v, w);
                self.delete_triangle(v, u, x);

                self.add_triangle(x, w, u);
                self.add_triangle(w, x, v);
            },
            // no quad found
            _ => ()
        }
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
    triangles: std::collections::hash_set::IntoIter<(VertexIdx, VertexIdx, VertexIdx)>,
    vertices: HashMap<(VertexIdx, VertexIdx), VertexIdx>,
}
impl Iterator for TriangleIntoIterator {
    type Item = [usize; 3];

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((u, v, w)) = self.triangles.next() {
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
fn in_circumcircle(u: VertexIdx, v: VertexIdx, w: VertexIdx, x: VertexIdx, vertices: &[Point2<f32>]) -> bool {
    let uv = u.get_vertex(vertices);
    let vv = v.get_vertex(vertices);
    let wv = w.get_vertex(vertices);
    let xv = x.get_vertex(vertices);

    // p is inside the triangle defined by (a, b, c) (given in counter-clockwise order) if:
    //       | ax-px, ay-py, (ax-px)² + (ay-py)² |
    // det = | bx-px, by-py, (bx-px)² + (by-py)² | > 0.0
    //       | cx-px, cy-py, (cx-px)² + (cy-py)² |
    let a1 = vv.x - uv.x;
    let b1 = wv.x - uv.x;
    let c1 = xv.x - uv.x ;

    let a2 = vv.y - uv.y;
    let b2 = wv.y - uv.y;
    let c2 = xv.y - uv.y;

    let a3 = a1*a1 + a2*a2;
    let b3 = b1*b1 + b2*b2;
    let c3 = c1*c1 + c2*c2;

    let det = a1*b2*c3 + a2*b3*c1 + b1*c2*a3 - c1*b2*a3 - c2*b3*a1 - b1*a2*c3;
    const EPS: f32 = 1e-12;
    if det > EPS {
        true
    } else if det < -EPS {
        false
    } else {
        // 4 cocircular point case
        // tie breaker to deterministically consider it as in circumcicle or not
        let min_idx = u.id().min(v.id().min(w.id().min(x.id())));
        min_idx == w.id() || min_idx == x.id()
    }
}

#[cfg(test)]
mod tests {
    use super::DelaunayTriangulation;
    use crate::geometry::coord::Point2;
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
        let vertices: &[Point2<f32>] = &[
            Point2::new(0.76809347, 0.17880994),
            Point2::new(0.14206064, 0.8896956),
            Point2::new(0.007408619, 0.17449331),
        ];
        for t in DelaunayTriangulation::from_vertices(&vertices).into_iter() {
            println!("{:?}", t);
        }
    }
    use image::{Rgb, RgbImage};
    use imageproc::drawing::draw_line_segment_mut;
    #[test]
    fn test_triangulate_complex() {
        let num_vertices = 1000;
        let vertices = (0..num_vertices)
            .map(|_| Point2::new(rand::random::<f32>(), rand::random::<f32>()))
            .collect::<Vec<_>>();

        let triangulation = DelaunayTriangulation::from_vertices(&vertices);
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