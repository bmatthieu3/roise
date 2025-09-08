use crate::geometry::coord::Vertex;
trait Noise<Sample: Vertex> {
    fn noise(&self, point: &Sample) -> f32;
}

mod diamond_square;
mod gradient;
mod worley;
pub use diamond_square::DiamondSquare;
pub use gradient::Gradient;
