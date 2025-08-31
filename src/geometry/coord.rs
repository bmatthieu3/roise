use core::fmt::Debug;
pub trait Vertex: Sub<Output = Self> + Normed + Add<Output = Self> + Clone + Copy + Debug {
    fn sample_around(&self, min_rad: f32, max_rad: f32) -> Self;
    fn default() -> Self;
}

pub trait Normed {
    fn magnitude(&self) -> f32;
    fn magnitude_squared(&self) -> f32;
}

impl Normed for Point2<f32> {
    fn magnitude(&self) -> f32 {
        self.magnitude_squared().sqrt()
    }

    fn magnitude_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Point2<T> {
    pub const fn new(x: T, y: T) -> Self {
        Point2 { x, y }
    }
}

impl Point2<f32> {
    pub fn dot(&self, other: &Self) -> f32 {
        self.x * other.x + other.y * self.y
    }

    pub fn det(&self, other: &Self) -> f32 {
        self.x * other.y - other.x * self.y
    }
}

use core::ops::Add;
impl Add for Point2<f32> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y)
    }
}

use std::iter::Sum;
impl Sum<Point2<f32>> for Point2<f32> {
    fn sum<I: Iterator<Item = Point2<f32>>>(iter: I) -> Self {
        iter.fold(Point2::new(0.0, 0.0), |acc, p| Point2::new(acc.x + p.x, acc.y + p.y))
    }
}

impl Add for &Point2<f32> {
    type Output = Point2<f32>;

    fn add(self, other: Self) -> Self::Output {
        Point2::new(self.x + other.x, self.y + other.y)
    }
}
impl Add<&Point2<f32>> for Point2<f32> {
    type Output = Self;

    fn add(self, other: &Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y)
    }
}
impl Add<Point2<f32>> for &Point2<f32> {
    type Output = Point2<f32>;

    fn add(self, other: Point2<f32>) -> Self::Output {
        Point2::new(self.x + other.x, self.y + other.y)
    }
}

use core::ops::Sub;
impl Sub for Point2<f32> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y)
    }
}
impl Sub for &Point2<f32> {
    type Output = Point2<f32>;

    fn sub(self, other: Self) -> Self::Output {
        Point2::new(self.x - other.x, self.y - other.y)
    }
}

use core::ops::Mul;
impl Mul<f32> for Point2<f32> {
    type Output = Self;

    fn mul(self, other: f32) -> Self::Output {
        Point2::new(self.x * other, self.y * other)
    }
}
impl Mul<f32> for &Point2<f32> {
    type Output = Point2<f32>;

    fn mul(self, other: f32) -> Self::Output {
        Point2::new(self.x * other, self.y * other)
    }
}

use core::ops::Div;
impl Div<f32> for Point2<f32> {
    type Output = Self;

    fn div(self, other: f32) -> Self::Output {
        Self::new(self.x / other, self.y / other)
    }
}
impl Div<f32> for &Point2<f32> {
    type Output = Point2<f32>;

    fn div(self, other: f32) -> Self::Output {
        Point2::new(self.x / other, self.y / other)
    }
}

impl Vertex for Point2<f32> {
    fn sample_around(&self, min_rad: f32, max_rad: f32) -> Self {
        let r = min_rad + rand::random::<f32>() * (max_rad - min_rad);
        let theta = 2.0 * std::f32::consts::PI * rand::random::<f32>();

        Point2::new(self.x + r * theta.cos(), self.y + r * theta.sin())
    }

    fn default() -> Self {
        Point2::new(0.0, 0.0)
    }
}
