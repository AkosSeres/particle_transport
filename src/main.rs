use crate::photon::set_detector;
use clap::Parser;

mod photon;
pub mod rand_gen;
mod vec3;

/// Simulates the transport of the photons of a monoenergetic gamma-source inside a scintillation detector.
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Radius of the detector in cm
    #[clap(long, default_value_t = 2.5)]
    radius: f64,

    /// Height of the detector in cm
    #[clap(short, long, default_value_t = 3.0)]
    height: f64,

    /// Energy of the emitted photons in keV
    #[clap(short, long, default_value_t = 661.7)]
    energy: f64,

    /// The density of the detector material, in g/cm³
    #[clap(short, long, default_value_t = 3.67)]
    density: f64,

    /// X coordinate of the emitter
    #[clap(long, default_value_t = 3.0)]
    rx: f64,
    /// Y coordinate of the emitter
    #[clap(long, default_value_t = -3.0)]
    ry: f64,
    /// Z coordinate of the emitter
    #[clap(long, default_value_t = 2.0)]
    rz: f64,

    /// Number of particles to evaluate
    #[clap(short, long, default_value_t = 10_000)]
    particlenum: u64,
}

fn main() {
    let args = Args::parse();
    set_detector(args.radius, args.height, args.density);

    println!("Hello, world! {:?}", args);
}
