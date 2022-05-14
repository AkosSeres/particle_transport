pub trait RandGen {
    fn rand() -> Self;
}

impl RandGen for f32 {
    #[inline(always)]
    fn rand() -> f32 {
        fastrand::f32()
    }
}

impl RandGen for f64 {
    #[inline(always)]
    fn rand() -> f64 {
        fastrand::f64()
    }
}
