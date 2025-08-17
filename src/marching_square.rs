use std::collections::HashMap;

use crate::coord::Point2;

// See this article for the figure giving the 16 possible cases
// https://nils-olovsson.se/articles/marching_squares/
const MARCHING_SQUARE_TABLE: &[&[(Point2<i8>, Point2<i8>)]] = &[
    // TL TR BR BL
    // 0  0  0  0
    &[],
    // 0  0  0  1
    &[(Point2::new(-1, 0), Point2::new(0, 1))],
    // 0  0  1  0
    &[(Point2::new(0, 1), Point2::new(1, 0))],
    // 0  0  1  1
    &[(Point2::new(-1, 0), Point2::new(1, 0))],
    // 0  1  0  0
    &[(Point2::new(1, 0), Point2::new(0, -1))],
    // 0  1  0  1
    &[
        (Point2::new(-1, 0), Point2::new(0, -1)),
        (Point2::new(1, 0), Point2::new(0, 1)),
    ],
    // 0  1  1  0
    &[(Point2::new(0, 1), Point2::new(0, -1))],
    // 0  1  1  1
    &[(Point2::new(-1, 0), Point2::new(0, -1))],
    // 1  0  0  0
    &[(Point2::new(0, -1), Point2::new(-1, 0))],
    // 1  0  0  1
    &[(Point2::new(0, -1), Point2::new(0, 1))],
    // 1  0  1  0
    &[
        (Point2::new(0, -1), Point2::new(1, 0)),
        (Point2::new(0, 1), Point2::new(-1, 0)),
    ],
    // 1  0  1  1
    &[(Point2::new(0, -1), Point2::new(1, 0))],
    // 1  1  0  0
    &[(Point2::new(1, 0), Point2::new(-1, 0))],
    // 1  1  0  1
    &[(Point2::new(1, 0), Point2::new(0, 1))],
    // 1  1  1  0
    &[(Point2::new(0, 1), Point2::new(-1, 0))],
    // 1  1  1  1
    &[],
];

fn extract_isocontours_from_heightmap<F>(num_sampling_vertices: usize, inside_area: F) -> Vec<Polygon>
where
    F: Fn(Point2<f32>) -> bool,
{
    let square_size = 1.0 / (num_sampling_vertices as f32 - 1.0);

    let mut edges: HashMap<(u16, u16), (u16, u16)> = HashMap::new();
    for i in 0..(num_sampling_vertices - 1) {
        let y = (i as f32) * square_size;
        for j in 0..(num_sampling_vertices - 1) {
            let x = (j as f32) * square_size;

            let tl_in = inside_area(Point2::new(x, y));
            let tr_in = inside_area(Point2::new(x + square_size, y));
            let bl_in = inside_area(Point2::new(x, y + square_size));
            let br_in = inside_area(Point2::new(x + square_size, y + square_size));

            let code = ((tl_in as usize) << 3) |((tr_in as usize) << 2) | ((br_in as usize) << 1) | (bl_in as usize);

            let i_half_cell = i << 1;
            let j_half_cell = j << 1;

            for (p1, p2) in MARCHING_SQUARE_TABLE[code] {
                let p1_x = ((j_half_cell as i32) + (p1.x as i32) + 1) as u16;
                let p1_y = ((i_half_cell as i32) + (p1.y as i32) + 1) as u16;

                let p2_x = ((j_half_cell as i32) + (p2.x as i32) + 1) as u16;
                let p2_y = ((i_half_cell as i32) + (p2.y as i32) + 1) as u16;

                edges.insert((p1_x, p1_y), (p2_x, p2_y));
            }
        }
    }

    let mut contours = vec![];
    while !edges.is_empty() {
        let mut c = vec![];
        // Extract one arbitrary edge to start the contour
        if let Some(mut start) = edges.keys().next().cloned() {
            let p1 = Point2::new(
                (start.0 as f32),
                (start.1 as f32),
            ) / ((num_sampling_vertices * 2) as f32);
            c.push(p1);

            while let Some(next) = edges.remove(&start) {
                let p2 = Point2::new(
                    (next.0 as f32),
                    (next.1 as f32),
                ) / ((num_sampling_vertices * 2) as f32);
                c.push(p2);

                start = next;
            }

            contours.push(Polygon {vertices: c});
        }
    }

    contours
}

struct Polygon {
    vertices: Vec<Point2<f32>>
}

impl Polygon {
    fn signed_area(&self) -> f32 {
        let mut i = self.vertices.len() - 1;
        let mut area = 0.0;
        for j in 0..self.vertices.len() {
            area += self.vertices[i].det(&self.vertices[j]);

            i = j;
        }

        area * 0.5
    }

    // Does not work for self intersecting polygons
    fn is_inside(&self, p: &Point2<f32>) -> bool {
        let mut i = self.vertices.len() - 1;
        let mut inside = false;
        for j in 0..self.vertices.len() {
            let Point2 { x: x1, y: y1 } = self.vertices[i];
            let Point2 { x: x2, y: y2 } = self.vertices[j];

            // Check if the edge crossed the line y = p.y
            if y1 <= p.y != y2 <= p.y {
                let xi = x2 + (p.y - y2) * (x2 - x1) / (y2 - y1);

                if xi <= p.x {
                    inside = !inside;
                }
            }

            i = j;
        }

        inside
    }
}

#[cfg(test)]
mod tests {
    use image::Rgb;
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;
    use image::RgbImage;

    use crate::marching_square::extract_isocontours_from_heightmap;
    use crate::marching_square::Polygon;
    use crate::nav_mesh::NavMesh;
    use crate::Point2;
    use crate::noise::Gradient;
    #[test]
    fn test_contour_extract() {
        let gradient = Gradient::new();
        let polygons = extract_isocontours_from_heightmap(200, |x: Point2<f32>| {
            let noise = gradient.fbm(&(x * 2.0), 0.6, 5.1)*0.707107 + 0.5; // in [0, 1]
            noise >= 0.45
        });

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for Polygon { vertices } in polygons {
            let color = Rgb([(rand::random::<f32>() * 255.0) as u8, (rand::random::<f32>() * 255.0) as u8, 255u8]);
            for (p1, p2) in vertices.iter().zip(vertices.iter().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    (((p1.x as f32) * w).round(), ((p1.y as f32) * h).round()),              // start point
                    (((p2.x as f32) * w).round(), ((p2.y as f32) * h).round()),            // end point
                    color, // RGB colors
                );
            }
        }

        img.save("countours.png").unwrap();
    }

    use crate::triangulate2;
    use std::collections::HashSet;
    #[test]
    fn test_triangulate_contours() {
        let gradient = Gradient::new();
        let polygons = extract_isocontours_from_heightmap(128, |x: Point2<f32>| {
            let noise = gradient.fbm(&(x * 2.0), 0.6, 5.1)*0.707107 + 0.5; // in [0, 1]
            noise >= 0.45
        });

        let vertices = polygons
            .into_iter()
            .flat_map(|Polygon { vertices }| vertices)
            .map(|p| ((p.x * 1e6_f32) as u64, (p.y * 1e6_f32) as u64))
            .collect::<HashSet<_>>()
            .into_iter()
            .map(|(x, y)| Point2::new(x as f32 / 1e6, y as f32 / 1e6))
            .collect::<Vec<_>>();

        let triangulation = triangulate2(&vertices);

        //panic!("jjj");


        let (w, h) = (1024.0, 1024.0);
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

        img.save("coutours_triangulated.png").unwrap();
    }

    #[test]
    fn point_inside_polygon() {
        let polygon = Polygon {
            vertices: vec![
                Point2::new()
                    #[test]
            ]
        };
    }
}