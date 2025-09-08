use crate::indexing::spatial::{Aabb, SpatialGrid};

pub struct Worley {
    grid: SpatialGrid<Point2<f32>>,
    num_cell_per_side: usize,
}

impl Worley {
    pub fn new(num_cell_per_side: usize) -> Self {
        let mut grid = SpatialGrid::new(1.0 / (num_cell_per_side as f32));

        for i in 0..num_cell_per_side {
            for j in 0..num_cell_per_side {
                let x = (i as f32 + rand::random::<f32>()) / (num_cell_per_side as f32);
                let y = (j as f32 + rand::random::<f32>()) / (num_cell_per_side as f32);

                let p = Point2::new(x, y);

                grid.insert(p);
            }
        }

        Self { grid, num_cell_per_side }
    }

    /*pub fn fbm(&self, p: &Point2<f32>, amplitude_factor: f32, freq_factor: f32) -> f32 {
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
    }*/
}

use crate::geometry::coord::{Normed, Point2};

use super::Noise;
impl Noise<Point2<f32>> for Worley {
    /// p given as coordinates between 0 and 1
    fn noise(&self, p: &Point2<f32>) -> f32 {
        let cell_size = 1.0 / (self.num_cell_per_side as f32);
        let mut nearest_points = self.grid.query_aabb(Aabb {
            xy_min: Point2::new(p.x - cell_size, p.y - cell_size),
            xy_max: Point2::new(p.x + cell_size, p.y + cell_size),
        }).collect::<Vec<_>>();

        nearest_points.sort_unstable_by(|&a, &b| {
            let ma = (*a - *p).magnitude_squared();
            let mb = (*b - *p).magnitude_squared();

            ma.partial_cmp(&mb).unwrap()
        });

        let f1 = (nearest_points[0] - p).magnitude();
        if nearest_points.len() >= 2 {
            let f2 = (nearest_points[1] - p).magnitude();

            f2 - f1
        } else {
            f1
        }
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
    use super::Worley;
    use crate::{geometry::coord::Point2, noise::Noise};

    use image::{ImageBuffer, Luma};

    #[test]
    fn test_worley() {
        let num_cells_per_side = 10;
        let cell_size = 1.0 / (num_cells_per_side as f32);
        let gradient = Worley::new(num_cells_per_side);
        let size = 512;
        let mut pixels = Vec::with_capacity(size * size);
        for i in 0..size {
            for j in 0..size {
                let p = Point2::new(
                    j as f32 / ((size - 1) as f32),
                    i as f32 / ((size - 1) as f32),
                );
                let noise = gradient.noise(&p) / cell_size;
                let color = (noise * 255.0) as u8;
                pixels.push(color);
            }
        }
        let image =
            ImageBuffer::<Luma<u8>, Vec<u8>>::from_raw(size as u32, size as u32, pixels).unwrap();
        let _ = image.save("worley.png");
    }
}
