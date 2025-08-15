mod space;
mod poisson;
mod random_uniform;
mod constant_step;

pub use poisson::{PoissonDisc, CustomDensity};
pub use random_uniform::RandUniform;
pub use constant_step::ConstantStepUniform;

pub use space::{Space, TwoDim};
pub trait Sampler<Sp>
where
    Sp: Space,
{
    fn sample(&self, f: &Sp) -> Vec<Sp::Sample>;
}