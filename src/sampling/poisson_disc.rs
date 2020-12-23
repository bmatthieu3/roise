use rand;
//use super::TorusTwoDim;
use crate::coord::{Normed, Point2};
use std::ops::Sub;

pub struct PoissonDisc<Sp, D>
where
    Sp: Space,
    D: Density<Sp>
{
    min_dist: D,
    s: std::marker::PhantomData<Sp>
}

impl<Sp, D> PoissonDisc<Sp, D>
where
    Sp: Space,
    D: Density<Sp>
{
    pub fn new(d: D) -> Self {
        PoissonDisc {
            min_dist: d,
            s: std::marker::PhantomData
        }
    }

    fn is_sampler_valid(grid: &<Sp as Space>::Partition, p: &<Sp as Space>::Sample, samples: &[<Sp as Space>::Sample], min_dist2: f32) -> bool {
        let mut neighbors = grid.neighbors(p, samples);
        while let Some(vertex_idx) = neighbors.pop() {
            for &id in vertex_idx {
                // Get the sample
                let s = &samples[id];
                let pt = *s - *p;
                if pt.magnitude_squared() < min_dist2 {
                    return false;
                }
            }
        }
        
        true
    }
}

use super::partition::EqualSizedGrid;
use crate::coord::Vertex;
const NUM_NEW_POINTS: usize = 10;
use super::{Sampler, Space};
impl<Sp, D> Sampler<Sp> for PoissonDisc<Sp, D>
where
    Sp: Space,
    D: Density<Sp>
{
    fn sample(&self, space: &Sp) -> Vec<Sp::Sample> {
        let max = self.min_dist.max();
        let cell_size = max / 2_f32.sqrt();

        let mut grid = Sp::Partition::new(cell_size);

        let sample = space.random();

        let idx_sample = grid.insert(&sample);

        let mut processed_points_idx = vec![idx_sample];
        let mut samples = vec![sample];
    
        while !processed_points_idx.is_empty() {
            let p_idx = processed_points_idx.pop().unwrap();
    
            for i in 0..NUM_NEW_POINTS {
                let dist = self.min_dist.get(&samples[p_idx]);
                let samp = samples[p_idx].sample_around(dist, 2.0*dist);

                if space.inside(&samp) && Self::is_sampler_valid(&grid, &samp, &samples, dist*dist) {
                    let idx_sample = grid.insert(&samp);
                    samples.push(samp);
                    processed_points_idx.push(idx_sample);

                    processed_points_idx.len();
                }
            }
        }
    
        samples
    }
}

pub trait Density<Sp: Space> {
    fn get(&self, s: &Sp::Sample) -> f32;
    fn max(&self) -> f32;
}

pub struct CustomDensity<Sp, C>
where
    Sp: Space,
    C: Fn(&Sp::Sample) -> f32
{
    constraint: C,
    max_dist: f32,
    space: std::marker::PhantomData<Sp>
}

impl<Sp, C> CustomDensity<Sp, C>
where
    Sp: Space,
    C: Fn(&Sp::Sample) -> f32
{
    pub fn new(constraint: C, max_dist: f32) -> Self {
        Self {
            constraint,
            max_dist,
            space: std::marker::PhantomData
        }
    }
}

impl<Sp, C> Density<Sp> for CustomDensity<Sp, C>
where
    Sp: Space,
    C: Fn(&Sp::Sample) -> f32
{
    fn get(&self, s: &Sp::Sample) -> f32 {
        (self.constraint)(s)
    }

    fn max(&self) -> f32 {
        self.max_dist
    }
}

#[cfg(test)]
mod tests {
    use super::{Sampler, PoissonDisc, CustomDensity};
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
        let offset = 0.01;

        let s = PoissonDisc::new(
            CustomDensity::new(
                |x: &Point2| {
                    //let p = *x - Point2::new(0.5, 0.5);
                    //let r = p.magnitude();

                    //(r.sqrt()*0.03).max(0.0005).min(0.05)


                    let alpha = ((x.x * 40.0).cos() * 0.5 + 0.5)*amplitude + offset;
                    alpha.min(amplitude + offset).max(offset)
                },
                amplitude + offset
                /*|x: &Point2| {
                    0.02
                },
                0.02*/
            )
        );

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
        }
        img.save("poisson.png").unwrap();
        */

        let triangulation = crate::triangulate2(vertices.as_slice());
        let (w, h) = (512.0, 512.0);
        let mut img = RgbImage::new(w as u32, h as u32);
        for t in triangulation.into_iter() {
            for (&idx1, &idx2) in t.iter().zip(t.iter().skip(1).cycle()) {
                draw_line_segment_mut(
                    &mut img,
                    (vertices[idx1].x * w, vertices[idx1].y * h),              // start point
                    (vertices[idx2].x * w, vertices[idx2].y * h),            // end point
                    Rgb([69u8, 203u8, 133u8]), // RGB colors
                );
            }
        }

        img.save("poisson.png").unwrap();
    }
}
