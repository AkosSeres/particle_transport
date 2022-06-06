#[cfg(target_arch = "wasm32")]
use nanorand::{Rng, WyRand};

#[cfg(target_arch = "wasm32")]
static mut rng: WyRand = WyRand::new_seed(6454353);

#[cfg(target_arch = "wasm32")]
pub fn set_rng_seed(seed: u64) {
    unsafe {
        rng = WyRand::new_seed(seed);
    }
}

pub trait RandGen {
    fn rand() -> Self;
}

impl RandGen for f32 {
    #[inline(always)]
    fn rand() -> f32 {
        #[cfg(not(target_arch = "wasm32"))]
        return fastrand::f32();
        #[cfg(target_arch = "wasm32")]
        unsafe {
            return rng.generate::<f32>();
        }
    }
}

impl RandGen for f64 {
    #[inline(always)]
    fn rand() -> f64 {
        #[cfg(not(target_arch = "wasm32"))]
        return fastrand::f64();
        #[cfg(target_arch = "wasm32")]
        unsafe {
            return rng.generate::<f64>();
        }
    }
}
