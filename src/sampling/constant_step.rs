use crate::sampling::TwoDim;
pub struct ConstantStepUniform {
    cell_size: f32
}
use crate::coord::{Vertex, Point2};

impl ConstantStepUniform {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size
        }
    }
}
use super::{Sampler, Space};
impl<F> Sampler<TwoDim<F>> for ConstantStepUniform
where
    F: Fn(&Point2) -> bool
{
    fn sample(&self, space: &TwoDim<F>) -> Vec<<TwoDim<F> as Space>::Sample> {
        let num_cell_side = (1.0 / self.cell_size) as usize + 1;

        (0..num_cell_side)
            .map(|i| (0..num_cell_side).map(move |j| Point2::new(i as f32 / ((num_cell_side - 1) as f32), j as f32 / ((num_cell_side - 1) as f32))))
            .flatten()
            .filter(|p| {
                space.inside(p)
            })
            .collect()
    }
}