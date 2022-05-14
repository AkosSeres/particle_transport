use crate::vec3::Vector;

type F = f64;

/// Radius of the detector
static mut DET_R: F = 2.5;
/// Height of the detector
static mut DET_HEIGHT: F = 3.0;
/// Density of the detector material, in g/cm³
static mut DET_DENSITY: F = 3.67;
/// Sets the properties of the detector
pub fn set_detector(r: F, height: F, density: F) {
    unsafe {
        DET_R = r;
        DET_HEIGHT = height;
        DET_DENSITY = density;
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Photon {
    energy: F,
    pos: Vector<F>,
    dir: Vector<F>,
}

impl Photon {
    pub fn new(energy: F, pos: Vector<F>, dir: Vector<F>) -> Self {
        Self { energy, pos, dir }
    }
}
