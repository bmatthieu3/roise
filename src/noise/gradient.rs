/// A simple gradient noise
pub struct Gradient;

impl Gradient {
    fn new() -> Self {
        Self {}
    }

    fn fbm(&self, p: &Point2, amplitude_factor: f32, freq_factor: f32) -> f32 {
        let octave = 4;
        let mut amplitude = 1.0;
        let mut freq = 1.0;
        let mut noise = 0.0;
        for i in 0..octave {
            noise += amplitude * self.noise(&(p * freq));
            amplitude *= amplitude_factor;
            freq *= freq_factor;
        }

        noise
    }
}

use crate::coord::{Point2, Vertex};
use crate::sampling::Space;

use super::Noise;
impl Noise<Point2> for Gradient {
    /// p given as coordinates between 0 and 1
    fn noise(&self, p: &Point2) -> f32 {
        let x0 = (p.x as usize) as f32;
        let y0 = (p.y as usize) as f32;

        let u = p.x - x0;
        let v = p.y - y0;

        let a = random_vect(&Point2::new(x0, y0)).dot(&Point2::new(u, v));
        let b = random_vect(&Point2::new(x0 + 1.0, y0)).dot(&Point2::new(u - 1.0, v));
        let c = random_vect(&Point2::new(x0 + 1.0, y0 + 1.0)).dot(&Point2::new(u - 1.0, v - 1.0));
        let d = random_vect(&Point2::new(x0, y0 + 1.0)).dot(&Point2::new(u, v - 1.0));

        let u = u * u * (3.0 - 2.0 * u);
        let v = v * v * (3.0 - 2.0 * v);

        lerp(v, lerp(u, a, b), lerp(u, d, c))
    }
}

fn random_f(p: &Point2) -> f32 {
    let t = p.dot(&Point2::new(12.9898, 78.233));
    (t.sin() * 43758.5453123).fract()
}

fn lerp(x: f32, a0: f32, a1: f32) -> f32 {
    (1.0 - x) * a0 + x * a1
}

fn random_vect(p: &Point2) -> Point2 {
    let theta = 2.0 * std::f32::consts::PI * random_f(p);

    Point2::new(theta.cos(), theta.sin())
}

#[cfg(test)]
mod tests {
    use super::{Gradient};
    use crate::{
        noise::Noise,
        coord::Point2
    };
    use image::{ImageBuffer, Luma, Rgb};
    use crate::sampling::{Sampler, PoissonDisc, TwoDim, CustomDensity};
    use crate::triangulation::triangulate2;
    use image::RgbImage;
    use imageproc::drawing::draw_cross_mut;
    use imageproc::drawing::draw_line_segment_mut;
    #[test]
    fn test_gradient() {
        let gradient = Gradient::new();
        let size = 512;
        let mut pixels = Vec::with_capacity(size * size);
        for i in 0..size {
            for j in 0..size {
                let p = Point2::new(j as f32 / ((size - 1) as f32), i as f32 / ((size - 1) as f32));
                let noise = gradient.fbm(&p, 0.6, 3.0);
                // gradient algorithm has a range between [-sqrt(N/4); sqrt(N/4)]
                // where N is the number of dimension.
                // For N = 2, we need to scale by sqrt(2)/2 and offset by 0.5 to have result between
                // [0; 1]
                let noise = noise*0.707107 + 0.5;
                let color = if noise < 0.5 {
                    0
                } else {
                    (noise * 255.0) as u8
                };
                //let color = (noise * 255.0) as u8;
                pixels.push(color);
            }
        }
        let mut image = ImageBuffer::<Luma<u8>, Vec<u8>>::from_raw(size as u32, size as u32, pixels).unwrap();
        image.save("gradient.png");
    }

    #[test]
    fn test_nav_mesh() {
        let gradient = Gradient::new();

        // Init a poisson disc sampling for the trees
        let s = PoissonDisc::new(
            CustomDensity::new(
                |x: &Point2| {
                    let noise = gradient.fbm(&(x * 2.0), 0.6, 3.0)*0.707107 + 0.5; // in [0, 1]
                    if noise < 0.45 {
                        0.01
                    } else {
                        0.05
                    }
                },
                0.05
            )        
        );

        let vertices = s.sample(
            &TwoDim::new(
                |p| true
            )
        );

        let vertices = vertices
            .into_iter()
            .map(|p| <na::Point2<f32> as From<Point2>>::from(p))
            .collect::<Vec<_>>();

        /*let (w, h) = (1024.0, 1024.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for p in vertices.iter() {
            draw_cross_mut(
                &mut img,
                Rgb([69u8, 203u8, 133u8]),
                (p.x * 1024.0) as i32,              // start point
                (p.y * 1024.0) as i32,            // end point
            );
        }
        img.save("sampling_trees.png").unwrap();*/

        let triangulation = triangulate2(&vertices);
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
