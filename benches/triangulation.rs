use criterion::{black_box, criterion_group, criterion_main, Criterion};

use roise::{triangulate, triangulate2};
use nalgebra as na;

fn criterion_benchmark(c: &mut Criterion) {
    let num_vertices = 10000;
    let vertices = (0..num_vertices)
        .map(|_| na::Point2::new(rand::random::<f32>(), rand::random::<f32>()))
        .collect::<Vec<_>>();
    c.bench_function("triangulation", |b| b.iter(|| triangulate(black_box(&vertices))));
    c.bench_function("triangulation2", |b| b.iter(|| triangulate2(black_box(&vertices))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);