mod constant_step;
mod poisson;
mod random_uniform;
mod space;

pub use constant_step::ConstantStepUniform;
pub use poisson::{CustomDensity, PoissonDisc};
pub use random_uniform::RandUniform;

pub use space::{Space, TwoDim};
pub trait Sampler<Sp>
where
    Sp: Space,
{
    fn sample(&self, f: &Sp) -> Vec<Sp::Sample>;
}
