use criterion::{black_box, criterion_group, criterion_main, Criterion};

use roise::{triangulate, Point2};

fn criterion_benchmark(c: &mut Criterion) {
    let num_vertices = 10000;
    let vertices = (0..num_vertices)
        .map(|_| Point2::new(rand::random::<f32>(), rand::random::<f32>()))
        .collect::<Vec<_>>();
    c.bench_function("triangulation", |b| {
        b.iter(|| triangulate(black_box(&vertices)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
