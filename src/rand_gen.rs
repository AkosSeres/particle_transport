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

    fn rand_normal_sum12() -> Self
    where
        Self: num::Float,
        Self: std::iter::Sum<Self>,
        f64: Into<Self>,
        f32: Into<Self>,
    {
        (0..12).map(|_| Self::rand()).sum::<Self>() - 6.0.into()
    }

    fn rand_normal_box_muller() -> (Self, Self)
    where
        Self: num::Float,
        Self: std::ops::MulAssign<Self>,
        f64: Into<Self>,
        f32: Into<Self>,
    {
        let usq = (Self::rand().ln() * (-2.0).into()).sqrt();
        let v2p = Self::rand() * (2.0 * 3.141592653589793).into();
        let mut ret = v2p.sin_cos();
        ret.0 *= usq;
        ret.1 *= usq;
        ret
    }
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
