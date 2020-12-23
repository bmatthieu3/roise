pub struct RandUniform {
    num_samples: usize
}

use super::{Sampler, Space};
impl<Sp> Sampler<Sp> for RandUniform
where
    Sp: Space
{
    fn sample(&self, space: &Sp) -> Vec<Sp::Sample> {
        (0..self.num_samples)
            .map(|_| space.random())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{Sampler, RandUniform};
    use crate::{
        sampling::space::TwoDim,
        coord::{Point2, Normed}
    };
    use image::{Rgb, RgbImage};
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;
    #[test]
    fn test_poisson_disc() {
        let amplitude = 0.01;
        let offset = 0.005;

        let s = RandUniform {
            num_samples: 1000,
        };

        let vertices = s.sample(
            &TwoDim::new(
                |p| {
                    let p = *p - Point2::new(0.5, 0.5);
                    let r = p.magnitude();

                    if r <= 1.0 && r >= 0.2 {
                        true
                    } else {
                        false
                    }
                }
            )
        );
        let vertices = vertices
            .into_iter()
            .map(|p| <na::Point2<f32> as From<Point2>>::from(p))
            .collect::<Vec<_>>();

        /*let (w, h) = (512.0, 512.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for p in vertices.iter() {
            draw_cross_mut(
                &mut img,
                Rgb([69u8, 203u8, 133u8]),
                (p.x * 512.0) as i32,              // start point
                (p.y * 512.0) as i32,            // end point
            );
        }*/

        let triangulation = crate::triangulate(vertices.as_slice());
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

        img.save("uniform.png").unwrap();
    }
}
