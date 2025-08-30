/// A simple gradient noise
pub struct Gradient;

impl Default for Gradient {
    fn default() -> Self {
        Self::new()
    }
}

impl Gradient {
    pub fn new() -> Self {
        Self {}
    }

    pub fn fbm(&self, p: &Point2<f32>, amplitude_factor: f32, freq_factor: f32) -> f32 {
        let octave = 4;
        let mut amplitude = 1.0;
        let mut freq = 1.0;
        let mut noise = 0.0;
        for _ in 0..octave {
            noise += amplitude * self.noise(&(p * freq));
            amplitude *= amplitude_factor;
            freq *= freq_factor;
        }

        noise
    }
}

use crate::geometry::coord::Point2;

use super::Noise;
impl Noise<Point2<f32>> for Gradient {
    /// p given as coordinates between 0 and 1
    fn noise(&self, p: &Point2<f32>) -> f32 {
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

fn random_f(p: &Point2<f32>) -> f32 {
    let t = p.dot(&Point2::new(12.9898, 78.233));
    (t.sin() * 43_758.547).fract()
}

fn lerp(x: f32, a0: f32, a1: f32) -> f32 {
    (1.0 - x) * a0 + x * a1
}

fn random_vect(p: &Point2<f32>) -> Point2<f32> {
    let theta = 2.0 * std::f32::consts::PI * random_f(p);

    Point2::new(theta.cos(), theta.sin())
}

#[cfg(test)]
mod tests {
    use super::{Gradient};
    use crate::geometry::coord::Point2;
    

    use image::{ImageBuffer, Luma};
    
    
    
    #[test]
    fn test_gradient() {
        let gradient = Gradient::new();
        let size = 512;
        let mut pixels = Vec::with_capacity(size * size);
        for i in 0..size {
            for j in 0..size {
                let p = Point2::new(j as f32 / ((size - 1) as f32), i as f32 / ((size - 1) as f32));
                let noise = gradient.fbm(&p, 0.9, 10.0);
                // gradient algorithm has a range between [-sqrt(N/4); sqrt(N/4)]
                // where N is the number of dimension.
                // For N = 2, we need to scale by sqrt(2)/2 and offset by 0.5 to have result between
                // [0; 1]
                let noise = noise*std::f32::consts::FRAC_1_SQRT_2 + 0.5;
                
                let noise = noise.powf(4.0).abs() * (-1.0) + 1.0;
                let color = (noise * 255.0) as u8;
                pixels.push(color);
            }
        }
        let image = ImageBuffer::<Luma<u8>, Vec<u8>>::from_raw(size as u32, size as u32, pixels).unwrap();
        let _ = image.save("gradient.png");
    }
}
