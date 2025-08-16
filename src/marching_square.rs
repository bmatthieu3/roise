use std::collections::HashMap;

use crate::coord::Point2;

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

fn extract_edges<F>(num_sampling_vertices: usize, inside_area: F) -> HashMap<(u16, u16), (u16, u16)>
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

    
    /*while !g.is_empty() {
        // Get one 
        if let Some((&key, &value)) = g.iter().next() {
            println!("Random pick: {} → {}", key, value);
    
            // Remove it
            g.remove(key);
        }

        let a = g.remove(&());
    }*/

    edges
}

#[cfg(test)]
mod tests {
    use image::Rgb;
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;
    use image::RgbImage;

    use crate::marching_square::extract_edges;
    use crate::nav_mesh::NavMesh;
    use crate::Point2;
    use crate::noise::Gradient;
    #[test]
    fn test_contour_extract() {
        let gradient = Gradient::new();
        let mut edges = extract_edges(200, |x: Point2<f32>| {
            let noise = gradient.fbm(&(x * 2.0), 0.6, 5.1)*0.707107 + 0.5; // in [0, 1]
            noise >= 0.45
        });

        let mut contours = vec![];
        while !edges.is_empty() {
            let mut c = vec![];
            // Extract one arbitrary edge to start the contour
            if let Some(mut start) = edges.keys().next().cloned() {
                c.push(start);

                while let Some(next) = edges.remove(&start) {
                    c.push(next);

                    start = next;
                }

                contours.push(c);
            }
        }


        let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for c in contours {
            let color = Rgb([(c[0].0 & 0xff) as u8, (c[0].1 & 0xff) as u8, 255u8]);
            for (p1, p2) in c.iter().zip(c.iter().skip(1)) {
                draw_line_segment_mut(
                    &mut img,
                    ((((p1.0 as f32) / 400.0) * w).round(), (((p1.1 as f32) / 400.0) * h).round()),              // start point
                    ((((p2.0 as f32) / 400.0) * w).round(), (((p2.1 as f32) / 400.0) * h).round()),            // end point
                    color, // RGB colors
                );
            }
        }

        img.save("countours.png").unwrap();
    }
}