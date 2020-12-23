use core::fmt::Debug;
pub trait Vertex: Sub<Output=Self> + Normed + Add<Output=Self> + Clone + Copy + Debug {
    fn sample_around(&self, min_rad: f32, max_rad: f32) -> Self;
    fn default() -> Self;
}

pub trait Normed {
    fn magnitude(&self) -> f32;
    fn magnitude_squared(&self) -> f32;
}

impl Normed for Point2 {
    fn magnitude(&self) -> f32 {
        (self.x*self.x + self.y*self.y).sqrt()
    }

    fn magnitude_squared(&self) -> f32 {
        self.x*self.x + self.y*self.y
    }
}

#[derive(Copy, Clone, Debug)]
#[derive(PartialEq)]
pub struct Point2(na::Point2<f32>);

impl Point2 {
    pub fn new(x: f32, y: f32) -> Self {
        Point2(na::Point2::new(x, y))
    }

    pub fn dot(&self, other: &Point2) -> f32 {
        self.x*other.x + other.y*self.y
    }

    pub fn det(&self, other: &Point2) -> f32 {
        self.x*other.y - other.x*self.y
    }
}
use core::ops::Deref;
impl Deref for Point2 {
    type Target = na::Point2<f32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

use core::ops::Add;
impl Add for Point2 {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Point2::new(
            self.x + other.x,
            self.y + other.y
        )
    }
}
impl Add for &Point2 {
    type Output = Point2;

    fn add(self, other: Self) -> Self::Output {
        Point2::new(
            self.x + other.x,
            self.y + other.y
        )
    }
}
impl Add<&Point2> for Point2 {
    type Output = Point2;

    fn add(self, other: &Self) -> Self::Output {
        Point2::new(
            self.x + other.x,
            self.y + other.y
        )
    }
}
impl Add<Point2> for &Point2 {
    type Output = Point2;

    fn add(self, other: Point2) -> Self::Output {
        Point2::new(
            self.x + other.x,
            self.y + other.y
        )
    }
}

use core::ops::Sub;
impl Sub for Point2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Point2::new(
            self.x - other.x,
            self.y - other.y
        )
    }
}
impl Sub for &Point2 {
    type Output = Point2;

    fn sub(self, other: Self) -> Self::Output {
        Point2::new(
            self.x - other.x,
            self.y - other.y
        )
    }
}

use core::ops::Mul;
impl Mul<f32> for Point2 {
    type Output = Self;

    fn mul(self, other: f32) -> Self::Output {
        Point2::new(
            self.x * other,
            self.y * other
        )
    }
}
impl Mul<f32> for &Point2 {
    type Output = Point2;

    fn mul(self, other: f32) -> Self::Output {
        Point2::new(
            self.x * other,
            self.y * other
        )
    }
}

use core::ops::Div;
impl Div<f32> for Point2 {
    type Output = Self;

    fn div(self, other: f32) -> Self::Output {
        Point2::new(
            self.x / other,
            self.y / other
        )
    }
}
impl Div<f32> for &Point2 {
    type Output = Point2;

    fn div(self, other: f32) -> Self::Output {
        Point2::new(
            self.x / other,
            self.y / other
        )
    }
}

impl From<Point2> for na::Point2<f32> {
    fn from(p: Point2) -> Self {
        *p
    }
}

impl Vertex for Point2 {
    fn sample_around(&self, min_rad: f32, max_rad: f32) -> Self {
        let r = min_rad + rand::random::<f32>() * (max_rad - min_rad);
        let theta = 2.0 * std::f32::consts::PI * rand::random::<f32>();

        Point2::new(
            self.x + r * theta.cos(),
            self.y + r * theta.sin()
        )
    }

    fn default() -> Self {
        Point2::new(0.0, 0.0)
    }
}