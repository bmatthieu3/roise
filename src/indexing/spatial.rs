use std::hash::Hash;
use crate::Point2;

trait Spatial {
    type Id: Copy + Eq + Hash;

    /// Bounding box (for fast insertion/search)
    //fn aabb(&self) -> Aabb;

    /// Optional: exact geometry check if needed
    fn contains(&self, point: Point2<f32>) -> bool;

    fn id(&self) -> Self::Id;
}