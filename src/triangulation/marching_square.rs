use std::collections::HashMap;

use crate::geometry::coord::Point2;
use crate::geometry::closed_polyline::ClosedPolyline;

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

pub fn extract_isocontours_from_heightmap<F>(num_sampling_vertices: i32, inside_area: F) -> Vec<ClosedPolyline>
where
    F: Fn(Point2<f32>) -> bool,
{
    let square_size = 1.0 / (num_sampling_vertices as f32 - 1.0);

    let mut edges: HashMap<(i32, i32), (i32, i32, usize)> = HashMap::new();
    for i in (-1)..(num_sampling_vertices+1) {
        let y = (i as f32) * square_size;
        for j in (-1)..(num_sampling_vertices+1) {
            let x = (j as f32) * square_size;

            let tl_in = if i == -1 || j == -1 {
                false
            } else {
                inside_area(Point2::new(x, y))
            };
            let tr_in = if i == -1 || j == num_sampling_vertices {
                false
            } else {
                inside_area(Point2::new(x + square_size, y))
            };
            let bl_in = if i == num_sampling_vertices || j == -1 {
                false
            } else {
                inside_area(Point2::new(x, y + square_size))
            };
            let br_in = if i == num_sampling_vertices || j == num_sampling_vertices {
                false
            } else {
                inside_area(Point2::new(x + square_size, y + square_size))
            };

            let code = ((tl_in as usize) << 3) | ((tr_in as usize) << 2) | ((br_in as usize) << 1) | (bl_in as usize);

            let i_half_cell = i << 1;
            let j_half_cell = j << 1;

            for (p1, p2) in MARCHING_SQUARE_TABLE[code] {
                let p1_x = j_half_cell + (p1.x as i32) + 1;
                let p1_y = i_half_cell + (p1.y as i32) + 1;

                let p2_x = j_half_cell + (p2.x as i32) + 1;
                let p2_y = i_half_cell + (p2.y as i32) + 1;

                edges.insert((p1_x, p1_y), (p2_x, p2_y, code));
            }
        }
    }

    let mut contours = vec![];
    let mut idx: usize = 0;
    let mut outer_contour = 0;
    let mut max_area = 0.0;

    while !edges.is_empty() {
        let mut vertices = vec![];
        // Extract one arbitrary edge to start the contour
        if let Some(mut start) = edges.keys().next().cloned() {
            let p1 = Point2::new(
                (start.0) as f32,
                (start.1) as f32,
            ) / (((num_sampling_vertices) * 2) as f32);
            vertices.push(p1);

            while let Some((cx, cy, ccode)) = edges.remove(&start) {
                let cur_vertex = Point2::new(
                    cx as f32,
                    cy as f32,
                ) / (((num_sampling_vertices) * 2) as f32);

                // peek the next one to see if its code matches.
                // If so we should we do not need to add the cur vertex to the contour
                match edges.get(&(cx, cy)) {
                    Some((_, _, ncode)) if *ncode == ccode => (),
                    _ => vertices.push(cur_vertex)
                }
                start = (cx, cy);
            }

            assert_eq!(vertices[0], vertices[vertices.len() - 1]);

            let polygon = ClosedPolyline {vertices};

            let area = polygon.signed_area().abs();
            // find the outer contour
            if area > max_area {
                max_area = area;
                outer_contour = idx;
            }

            contours.push(polygon);

            idx += 1;
        }
    }

    // put the outer contour at the first position
    contours.swap(0, outer_contour);

    // simplify the contours using douglas-peucker rec algo
    

    contours.into_iter()
        .filter_map(|mut contour| {
            let last = contour.vertices.pop().unwrap();
            contour.vertices = crate::geometry::closed_polyline::douglas_peucker(&contour.vertices[0..(contour.vertices.len() - 1)], square_size * 0.5);
            contour.vertices.push(last);

            if contour.vertices.len() == 3 {
                None
            } else {
                Some(contour)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use image::Rgb;
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;
    use image::RgbImage;

    use crate::triangulation::marching_square::extract_isocontours_from_heightmap;
    use crate::geometry::closed_polyline::ClosedPolyline;
    
    use crate::triangulation::DelaunayTriangulation;
    use crate::geometry::coord::Point2;
    use crate::noise::Gradient;
    #[test]
    fn test_contour_extract() {
        let gradient = Gradient::new();
        let polygons = extract_isocontours_from_heightmap(200, |x: Point2<f32>| {
            let noise = gradient.fbm(&(x * 2.0), 0.6, 3.0)*std::f32::consts::FRAC_1_SQRT_2 + 0.5; // in [0, 1]
            noise >= 0.45
        });

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for ClosedPolyline { vertices } in polygons {
            let color = Rgb([(rand::random::<f32>() * 255.0) as u8, (rand::random::<f32>() * 255.0) as u8, 255u8]);
            for (p1, p2) in vertices.iter().zip(vertices.iter().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    ((p1.x * w).round(), (p1.y * h).round()),              // start point
                    ((p2.x * w).round(), (p2.y * h).round()),            // end point
                    color, // RGB colors
                );
            }

            for p in vertices.iter() {
                draw_cross_mut(
                    &mut img,
                    Rgb([200u8, 203u8, 133u8]),
                    (p.x * 1024.0) as i32,              // start point
                    (p.y * 1024.0) as i32,            // end point
                );
            }
        }

        img.save("contours.png").unwrap();
    }

    
    #[test]
    fn test_triangulate_contours() {
        let gradient = Gradient::new();
        let polygons = extract_isocontours_from_heightmap(200, |x: Point2<f32>| {
            let noise = gradient.fbm(&(x * 2.0), 0.6, 3.0)*std::f32::consts::FRAC_1_SQRT_2 + 0.5; // in [0, 1]
            noise >= 0.45
        });

        let vertices = polygons
            .into_iter()
            .flat_map(|ClosedPolyline { mut vertices }| {
                let _ = vertices.pop();
                vertices
            })
            .collect::<Vec<_>>();

        let triangulation = DelaunayTriangulation::from_vertices(&vertices);

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation {
            for (&idx1, &idx2) in t.iter().zip(t.iter().cycle().skip(1)) {
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
    fn test_triangulate_cdt_from_contours() {
        let gradient = Gradient::new();
        let contours = extract_isocontours_from_heightmap(200, |x: Point2<f32>| {
            let noise = gradient.fbm(&(x * 2.1), 0.6, 3.01)*std::f32::consts::FRAC_1_SQRT_2 + 0.5; // in [0, 1]
            noise >= 0.45
        });

        let vertices = contours
            .iter()
            .cloned()
            .flat_map(|ClosedPolyline { mut vertices }| {
                let _ = vertices.pop();
                vertices
            })
            .collect::<Vec<_>>();

        let triangulation = DelaunayTriangulation::from_contours(&contours);

        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation {
            for (&idx1, &idx2) in t.iter().zip(t.iter().cycle().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h),              // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h),            // end point
                    Rgb([69u8, 203u8, 133u8]), // RGB colors
                );
            }
        }

        img.save("coutours_triangulated2.png").unwrap();
    }

    #[test]
    fn point_inside_polygon() {
        let polygon = ClosedPolyline {
            vertices: vec![
                Point2::new(0.0, 0.0),
                Point2::new(3.0, 0.0),
                Point2::new(3.0, 3.0),
                Point2::new(1.5, 1.5),
                Point2::new(0.0, 3.0),
                Point2::new(0.0, 0.0),
            ]
        };

        assert!(polygon.contains(&Point2::new(1.5, 1.0)));
        assert!(!polygon.contains(&Point2::new(-1.5, 1.0)));
        assert!(polygon.contains(&Point2::new(1.0, 1.0)));
    }
}