struct Triangulation {
    // storing the vertex indices in counter clockwise order
    vertices: Vec<usize>,
    // storing the triangles indices neighbours
    triangles: Vec<usize>,
}

/// Doubly circular linked chain
struct Hull {
    // triangle indices following the vertices located on the frontier.
    // triangles has the same size as frontier
    vertices: Vec<usize>,
    empty: Vec<usize>,
    // previous vertex index of the ith vertex on the frontier
    prev: Vec<usize>,
    // next vertex index  of the ith vertex on the frontier
    next: Vec<usize>,
    num: usize,
}

const UNASSIGNED: usize = std::usize::MAX;
use crate::coord::Point2;
impl Hull {
    fn new(num_points: usize, first_vertex_idx: usize) -> Self {
        let n = (num_points as f32).sqrt().ceil() as usize;
        let num = 1;
        let mut vertices = vec![UNASSIGNED; n];
        let mut prev = vec![UNASSIGNED; n];
        let mut next = vec![UNASSIGNED; n];
        let empty = (1..n).collect();

        vertices[0] = first_vertex_idx;
        next[0] = first_vertex_idx;
        prev[0] = first_vertex_idx;

        Self {
            vertices,
            empty,
            next,
            prev,
            num,
        }
    }

    fn insert_after(
        &mut self,
        // The vertex index in the hull after which the
        // new vertex will be inserted
        cur_idx: usize,
        // The vertex idx to insert
        vertex_idx: usize,
    ) {
        // Get the location where to insert vertex idx
        let idx = if self.empty.is_empty() {
            self.num
        } else {
            self.empty.pop().unwrap()
        };

        let next_cur = self.next[cur_idx];

        self.prev[idx] = cur_idx;
        self.next[idx] = next_cur;

        self.next[cur_idx] = idx;
        self.prev[next_cur] = idx;

        self.num += 1;
        self.vertices[idx] = vertex_idx;
    }

    fn remove(
        &mut self,
        // The vertex to remove from the hull
        cur_idx: usize,
    ) {
        let next_cur = self.next[cur_idx];
        let prev_cur = self.prev[cur_idx];

        self.prev[next_cur] = prev_cur;
        self.next[prev_cur] = next_cur;

        self.empty.push(cur_idx);
        self.num -= 1;
    }

    fn next(&self, idx: usize) -> usize {
        let next_idx = self.next[idx];
        self.vertices[next_idx]
    }

    fn prev(&self, idx: usize) -> usize {
        let prev_idx = self.prev[idx];
        self.vertices[prev_idx]
    }

    fn edges<'a>(&'a self, starting_idx: usize) -> EdgeIterator<'a> {
        EdgeIterator {
            num: 0,
            hull: self,
        }
    }
}

struct EdgeIterator<'a> {
    hull: &'a Hull,
    // num edges processed
    num: usize,
}

/*impl<'a> Iterator for EdgeIterator<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        // hull.num >= 1 by construction
        if self.num == self.hull.num - 1 {
            None
        } else {
            let u = self.hull.vertices[]


            self.num += 1;
        }
    }
}*/

trait Vertex {
    /// Check if the vertex is in the circumcircle defined
    /// by the vertices a, b and c
    fn in_circumcircle(&self, a: &Self, b: &Self, c: &Self) -> bool;

    /// Check if the vertex facing the edge defined by
    /// vertices a and b
    /// a.y < b.y
    fn is_facing_edge(&self, a: &Self, b: &Self) -> bool;
}
impl Vertex for Point2 {
    /// Check whether the vertex is in the circumcircle
    /// defined by the vertices a, b and c
    fn in_circumcircle(&self, a: &Self, b: &Self, c: &Self) -> bool {    
        // p is inside the triangle defined by (a, b, c) (given in counter-clockwise order) if:
        //       | ax-px, ay-py, (ax-px)² + (ay-py)² |
        // det = | bx-px, by-py, (bx-px)² + (by-py)² | > 0.0
        //       | cx-px, cy-py, (cx-px)² + (cy-py)² |
        let a1 = a.x - self.x;
        let b1 = b.x - self.x;
        let c1 = c.x - self.x;
    
        let a2 = a.y - self.y;
        let b2 = b.y - self.y;
        let c2 = c.y - self.y;
    
        let a3 = a1*a1 + a2*a2;
        let b3 = b1*b1 + b2*b2;
        let c3 = c1*c1 + c2*c2;
    
        a1*b2*c3 + a2*b3*c1 + b1*c2*a3 - c1*b2*a3 - c2*b3*a1 - b1*a2*c3 > 0.0
    }

    fn is_facing_edge(&self, a: &Self, b: &Self) -> bool {
        let ab = b - a;
        let av = self - a;

        av.det(&ab) > 0.0
    }
}

pub fn triangulate(vertices: &[na::Point2<f32>]) {}
#[cfg(test)]
mod tests {
    use crate::coord::Point2;
    use crate::VertexIdx;
    use super::Vertex;

    #[test]
    fn test_edge_facing_hull_vertex() {
        let v = Point2::new(1.0, 0.0);
        let a = Point2::new(0.0, 1.0);
        let b = Point2::new(-0.5, 2.0);

        assert!(v.is_facing_edge(&a, &b));
        assert!(!v.is_facing_edge(&a, &Point2::new(-1.0, 2.0)));
        assert!(!v.is_facing_edge(&a, &Point2::new(-1.0, 1.9)));
        assert!(!v.is_facing_edge(&Point2::new(0.0, 0.0), &Point2::new(-1.0, 0.0)));
        assert!(v.is_facing_edge(&Point2::new(0.0, -15.0), &Point2::new(0.0, -0.5)));
    }

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
        /*for t in triangulate2(&vertices).into_iter() {
            println!("{:?}", t);
        }*/
    }
    // use image::{Rgb, RgbImage};
    // use imageproc::drawing::draw_line_segment_mut;
    // #[test]
    // fn test_triangulate_complex() {
    //     let num_vertices = 1000;
    //     let vertices = (0..num_vertices)
    //         .map(|_| na::Point2::new(rand::random::<f32>(), rand::random::<f32>()))
    //         .collect::<Vec<_>>();

    //     /*for t in triangulate2(&vertices).into_iter() {
    //         println!("{:?}", t);
    //     }*/

    //     let triangulation = triangulate2(&vertices);
    //     let (w, h) = (512.0, 512.0);
    //     let mut img = RgbImage::new(w as u32, h as u32);
    //     for t in triangulation {
    //         for (&idx1, &idx2) in t.iter().zip(t.iter().skip(1).cycle()) {
    //             //let v1 = idx.get_vertex(super_triangle)

    //             draw_line_segment_mut(
    //                 &mut img,
    //                 (vertices[idx1].x * w, vertices[idx1].y * h), // start point
    //                 (vertices[idx2].x * w, vertices[idx2].y * h), // end point
    //                 Rgb([69u8, 203u8, 133u8]),                    // RGB colors
    //             );
    //         }
    //     }

    //     img.save("delaunay10.png").unwrap();
    //     //assert_eq!(connexity.triangles.len(), 3);
    // }
}
