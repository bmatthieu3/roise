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
use image::RgbImage;
use imageproc::drawing::draw_cross_mut;
use imageproc::drawing::draw_line_segment_mut;

use roise::{
    nav_mesh::NavMesh,
    triangulation::DelaunayTriangulation,
    geometry::closed_polyline::ClosedPolyline,
    Point2,
    noise::Gradient,
    graph::Metric
};

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

```

<img width="1024" height="1024" alt="nav_mesh_portals" src="https://github.com/user-attachments/assets/55f628f5-c97d-44b7-8bb8-208d68ec5716" />


