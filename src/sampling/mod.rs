mod space;
mod partition;
mod poisson_disc;
mod random_uniform;

pub use poisson_disc::{PoissonDisc, CustomDensity};
pub use random_uniform::RandUniform;

pub use space::{Space, TwoDim};
pub trait Sampler<Sp>
where
    Sp: Space,
{
    fn sample(&self, f: &Sp) -> Vec<Sp::Sample>;
}