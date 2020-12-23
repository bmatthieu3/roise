use crate::coord::{Point2, Vertex};
trait Noise<Sample: Vertex> {
    fn noise(&self, point: &Sample) -> f32;
}

mod diamond_square;
mod gradient;
pub use diamond_square::DiamondSquare;
pub use gradient::Gradient;