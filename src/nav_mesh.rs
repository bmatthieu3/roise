use crate::noise::Gradient;
use crate::sampling::{ConstantStepUniform, Space};
use crate::triangulation::DelaunayTriangulation;
use crate::sampling::Sampler;
use crate::Point2;
use crate::sampling::TwoDim;

pub(crate) struct NavMesh {
    pub triangulation: Vec<[usize; 3]>,
    pub vertices: Vec<Point2>
}
fn build_nav_mesh() -> NavMesh {
    let gradient = Gradient::new();

    let constant_sampler = ConstantStepUniform::new(0.02);
    let s = TwoDim::new(|x: &Point2| {
        let noise = gradient.fbm(&(x * 2.0), 0.6, 2.1)*0.707107 + 0.5; // in [0, 1]
        noise >= 0.45
    });

    let vertices = constant_sampler.sample(&s);
    let mut triangulation = super::triangulation::triangulate2(&vertices);

    let triangulation = triangulation.into_iter()
        .filter(|[a, b, c]| {
            let v1 = &vertices[*a];
            let v2 = &vertices[*b];
            let v3 = &vertices[*c];

            let b = super::triangulation::barycenter(v1, v2, v3);
            //s.inside(&b)
            s.inside(&v1) && s.inside(&v2) && s.inside(&v3)
        }).collect();

    NavMesh { 
        vertices,
        triangulation
    }
}

#[cfg(test)]
mod tests {
    use image::Rgb;
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;
    use image::RgbImage;

    use crate::nav_mesh::NavMesh;
    #[test]
    fn test_build_nav_mesh() {
        
        let NavMesh { triangulation, vertices } = super::build_nav_mesh();
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
        for p in vertices.iter() {
            draw_cross_mut(
                &mut img,
                Rgb([200u8, 203u8, 133u8]),
                (p.x * 1024.0) as i32,              // start point
                (p.y * 1024.0) as i32,            // end point
            );
        }

        img.save("nav_mesh.png").unwrap();
    }
}