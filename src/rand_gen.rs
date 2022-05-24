pub trait RandGen {
    fn rand() -> Self;
}

impl RandGen for f32 {
    #[inline(always)]
    fn rand() -> f32 {
        #[cfg(not(target_arch = "wasm32"))]
        return fastrand::f32();
        #[cfg(target_arch = "wasm32")]
        return js_sys::Math::random() as f32;
    }
}

impl RandGen for f64 {
    #[inline(always)]
    fn rand() -> f64 {
        #[cfg(not(target_arch = "wasm32"))]
        return fastrand::f64();
        #[cfg(target_arch = "wasm32")]
        return js_sys::Math::random();
    }
}
