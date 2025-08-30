use crate::geometry::coord::Point2;

#[derive(Clone)]
pub struct ClosedPolyline {
    pub vertices: Vec<Point2<f32>>,
}

impl ClosedPolyline {
    pub(crate) fn signed_area(&self) -> f32 {
        let mut i = self.vertices.len() - 1;
        let mut area = 0.0;
        for j in 0..self.vertices.len() {
            area += self.vertices[i].det(&self.vertices[j]);

            i = j;
        }

        area * 0.5
    }

    // Does not work for self intersecting polygons
    pub fn contains(&self, p: &Point2<f32>) -> bool {
        let mut i = self.vertices.len() - 1;
        let mut inside = false;
        for j in 0..self.vertices.len() {
            let Point2 { x: x1, y: y1 } = self.vertices[i];
            let Point2 { x: x2, y: y2 } = self.vertices[j];

            // Check if the edge crossed the line y = p.y
            if (y1 <= p.y) != (y2 <= p.y) {
                let xi = x2 + (p.y - y2) * (x2 - x1) / (y2 - y1);

                if xi < p.x {
                    inside = !inside;
                }
            }

            i = j;
        }

        inside
    }
}

/// Compute perpendicular distance from point P to segment [A,B]
fn perpendicular_distance(p: Point2<f32>, a: Point2<f32>, b: Point2<f32>) -> f32 {
    let ab = b - a;
    let ap = p - a;
    let len = (ab.x * ab.x + ab.y * ab.y).sqrt();
    if len < 1e-9 {
        (ap.x * ap.x + ap.y * ap.y).sqrt()
    } else {
        (ab.y * p.x - ab.x * p.y + b.x * a.y - b.y * a.x).abs() / len
    }
}

pub(crate) fn douglas_peucker(points: &[Point2<f32>], epsilon: f32) -> Vec<Point2<f32>> {
    if points.len() < 3 {
        return points.to_vec(); // nothing to simplify
    }

    // Find the point with the max distance from the line [first, last]
    let (mut index, mut max_dist) = (0, 0.0);
    for i in 1..(points.len() - 1) {
        let d = perpendicular_distance(points[i], points[0], points[points.len() - 1]);
        if d > max_dist {
            index = i;
            max_dist = d;
        }
    }

    // If the max distance is greater than epsilon, recursively simplify
    if max_dist > epsilon {
        let mut left = douglas_peucker(&points[0..=index], epsilon);
        let right = douglas_peucker(&points[index..], epsilon);

        left.pop(); // avoid duplicating the split point
        left.extend(right);
        left
    } else {
        // Otherwise just keep the endpoints
        vec![points[0], points[points.len() - 1]]
    }
}
