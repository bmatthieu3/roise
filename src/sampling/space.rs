use std::marker::Sized;
use super::partition::EqualSizedGrid;
use crate::coord::{Vertex, Point2};
pub trait Space: Sized {
    type Sample: Vertex;
    type Partition: EqualSizedGrid<Sp=Self>;

    fn inside(&self, p: &Self::Sample) -> bool;

    fn random(&self) -> Self::Sample {
        let mut s = Self::Sample::default();
        let mut inside = false;
        while !inside {
            s = Self::random_unconstrained();
            inside = self.inside(&s);
        }

        s
    }

    fn random_unconstrained() -> Self::Sample;
}

pub struct TwoDim<F>
where
    F: Fn(&Point2) -> bool
{
    constraint: F
}

impl<F> TwoDim<F>
where
    F: Fn(&Point2) -> bool
{
    pub fn new(constraint: F) -> Self {
        Self {
            constraint
        }
    }
}

use std::ops::Fn;
use super::partition::Equal2DSizedGrid;
impl<F> Space for TwoDim<F>
where F: Fn(&Point2) -> bool {
    type Sample = Point2;
    type Partition = Equal2DSizedGrid<F>;

    fn random_unconstrained() -> Self::Sample {
        Point2::new(
            rand::random::<f32>(),
            rand::random::<f32>()
        )
    }

    fn inside(&self, p: &Self::Sample) -> bool {
        // 1. chech whether it is in [0, 1] x [0, 1]
        let inside_full_space = p.x >= 0.0 && p.x < 1.0 &&
            p.y >= 0.0 && p.y < 1.0;
        if inside_full_space {
            // 2. check the constraint
            (self.constraint)(p)
        } else {
            false
        }
    }
}
