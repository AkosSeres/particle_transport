use criterion::{black_box, criterion_group, criterion_main, Criterion};
use particle_transport::{
    photon::{Photon, F},
    vec3::Vector,
};

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("simulate_default_settings", |b| {
        b.iter(|| {
            let mut random_photon = Photon {
                energy: black_box(661.7),
                pos: Vector::<F>::new(black_box(3.0), black_box(-3.0), black_box(2.0)),
                dir: Vector::<F>::random_isotropic_normed(),
            };
            let mut energy_hit_size = random_photon.simulate();
            if energy_hit_size <= 0.0 {
                return;
            }
            for _ in 0..12 {
                energy_hit_size += (fastrand::f64() - 0.5) * black_box(8.0);
            }
            let channel_width = black_box(700.0) / (black_box(1024) as f64);
            let idx = (energy_hit_size / channel_width).floor() as usize;
            if idx >= black_box(1024) {
                return;
            }
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
