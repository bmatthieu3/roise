Roise
=====

Roise targets game developers and aims at providing core features for developing a game engine

What it contains
----------------

* [x] A collection of noise algorithms and probability density methods
* [x] Contours marching square extraction from a heightmap
* [x] Constrained Delaunay Triangulation (CDT) algorithm
* [x] 2D Polygon API (contains vertex, compute area)
* [ ] Pathfinding (A*) on graph representation
* [ ] Navmesh from transforming a triangulation to a navigable graph
* [ ] Build a 3D mesh of a terrain from a heightmap

API examples
------------

1. Generating a triangulation from a heightmap, first step to build a nav mesh

This uses the image crate for rendering the 2D triangulation

```rust
use image::Rgb;
use imageproc::drawing::draw_cross_mut;
use imageproc::drawing::draw_line_segment_mut;
use image::RgbImage;

use crate::marching_square::extract_isocontours_from_heightmap;
use crate::marching_square::ClosedPolyline;
use crate::triangulation::DelaunayTriangulation;
use crate::Point2;
use crate::noise::Gradient;

let gradient = Gradient::new();
let contours = extract_isocontours_from_heightmap(200, |x: Point2<f32>| {
    let noise = gradient.fbm(&(x * 2.0), 0.6, 3.0)*0.707107 + 0.5; // in [0, 1]
    noise >= 0.45
});

// The contours are all closed polylines with their first and last vertex being the same
// Feeding doublons in a triangulation algorithm is a bad idea
// That is why we have to remove the last vertex before building the triangulation
let vertices = contours
    .iter()
    .cloned()
    .flat_map(|ClosedPolyline { mut vertices }| {
        let _ = vertices.pop();
        vertices
    })
    .collect::<Vec<_>>();

// Build the triangulation
let triangulation = DelaunayTriangulation::from_contours(&contours);

// Plot the triangulation
let (w, h) = (1024.0, 1024.0);
let mut img = RgbImage::new(w as u32, h as u32);
for t in triangulation {
    for (&idx1, &idx2) in t.iter().zip(t.iter().skip(1).cycle()) {
        draw_line_segment_mut(
            &mut img,
            (vertices[idx1].x * w, vertices[idx1].y * h),              // start point
            (vertices[idx2].x * w, vertices[idx2].y * h),            // end point
            Rgb([69u8, 203u8, 133u8]), // RGB colors
        );
    }
}

img.save("coutours_triangulated2.png").unwrap();
```
