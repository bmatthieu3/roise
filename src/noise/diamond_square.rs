const SIZE: usize = 513;
/// A pre-computed diamond square noise
/// This type of noise is efficient at designing
/// heightmap terrains
///
/// An amplitude factor between [0, 1] is multiplied each
/// to the amplitude. A value near 1 gives more importance to
/// the higher frequencies.
pub struct DiamondSquare {
    pixels: Vec<Vec<f32>>,
}

impl DiamondSquare {
    fn new(amplitude_factor: f32) -> Self {
        let mut pixels = vec![vec![0.0; SIZE]; SIZE];

        // top left
        pixels[0][0] = rand::random::<f32>();
        // top right
        pixels[0][SIZE - 1] = rand::random::<f32>();
        // bottom left
        pixels[SIZE - 1][0] = rand::random::<f32>();
        // bottom right
        pixels[SIZE - 1][SIZE - 1] = rand::random::<f32>();
        let mut step_size = SIZE - 1; // a power of two
        let mut amplitude = 0.5;
        while step_size > 1 {
            let half_step = (step_size >> 1) as i32;
            let x_start = half_step;
            let y_start = half_step;

            /// 2. Diamond step
            for y in (y_start..(SIZE as i32)).step_by(step_size) {
                for x in (x_start..(SIZE as i32)).step_by(step_size) {
                    let r = (rand::random::<f32>() - 0.5) * amplitude;
                    pixels[y as usize][x as usize] = 0.25 * (
                        pixels[(y - half_step) as usize][(x - half_step) as usize] +
                        pixels[(y - half_step) as usize][(x + half_step) as usize] +
                        pixels[(y + half_step) as usize][(x - half_step) as usize] + 
                        pixels[(y + half_step) as usize][(x + half_step) as usize]
                    ) + r;
                }
            }

            let mut offset = 0;
            /// 3. Square step
            for y in (0_i32..(SIZE as i32)).step_by(half_step as usize) {
                if offset == 0 {
                    offset = half_step;
                } else {
                    offset = 0;
                }
                for x in (offset..(SIZE as i32)).step_by(step_size) {
                    let mut num_acc = 0;
                    let x_off = dbg!(x) as i32;
                    pixels[y as usize][x_off as usize] = 0.0;

                    if x_off - half_step >= 0 {
                        pixels[y as usize][x_off as usize] += pixels[y as usize][(x_off - half_step) as usize];
                        num_acc += 1;
                    }
                    if x_off + half_step < SIZE as i32 {
                        pixels[y as usize][x_off as usize] += pixels[y as usize][(x_off + half_step) as usize];
                        num_acc += 1;
                    }

                    if y - half_step >= 0 {
                        pixels[y as usize][x_off as usize] += pixels[(y - half_step) as usize][x_off as usize];
                        num_acc += 1;
                    }
                    if y + half_step < SIZE as i32 {
                        pixels[y as usize][x_off as usize] += pixels[(y + half_step) as usize][x_off as usize];
                        num_acc += 1;
                    }
                    let r = (rand::random::<f32>() - 0.5) * amplitude;

                    pixels[y as usize][x_off as usize] /= (num_acc as f32);
                    pixels[y as usize][x_off as usize] += r;
                }
            }

            amplitude *= amplitude_factor;
            step_size = step_size >> 1;
        }

        Self {
            pixels,
        }
    }

    const fn get_size() -> usize {
        SIZE
    }
}

use crate::coord::{Point2, Vertex};
use crate::sampling::Space;

use super::Noise;
impl Noise<Point2> for DiamondSquare {
    /// p given as coordinates between 0 and 1
    fn noise(&self, p: &Point2) -> f32 {
        let x = ((SIZE as f32) * p.x) as usize;
        let y = ((SIZE as f32) * p.y) as usize;
        self.pixels[y][x]
    }
}

#[cfg(test)]
mod tests {
    use super::{DiamondSquare, SIZE};
    use image::{ImageBuffer, Luma, Rgb};
    #[test]
    fn test_diamond_square2() {
        let green = [0, 255, 0];
        let red = [255, 0, 0];

        let mut image = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(200, 200);
        image.put_pixel(5, 5, Rgb(green));
        image.save("output1.png");
    }

    #[test]
    fn test_diamond_square() {
        let diamond = DiamondSquare::new(0.9);
        let pixels = diamond.pixels.into_iter()
            .flatten()
            .map(|c| (255.0*c) as u8)
            .collect();
        let mut image = ImageBuffer::<Luma<u8>, Vec<u8>>::from_raw(SIZE as u32, SIZE as u32, pixels).unwrap();

        image.save("output1.jpg");
    }
}

