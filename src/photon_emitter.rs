use crate::{rand_gen::RandGen, vec3::Vector};

/// This is used to emit photons towards the detector,
/// in a controlled manner, instead of isotropically.
/// This saves a lot of computing time, since we won't shoot
/// unnecessary photons in the opposide direction.
pub struct PhotonEmitter {
    cos_alpha: f64,
    mid_dir: Vector<f64>,
    solid_angle: f64,
}

impl PhotonEmitter {
    /// Creates a new [PhotonEmitter] from the size of the detector
    /// and the desired coordinates of the point-like emitter.
    /// The detector is assumed to be centered at the origin.
    pub fn from_params(
        detector_r: f64,
        detector_h: f64,
        rx: f64,
        ry: f64,
        rz: f64,
    ) -> PhotonEmitter {
        // The bounding radius of the detector
        let bounding_radius = (detector_h * detector_h / 4.0 + detector_r * detector_r).sqrt();
        // Distance of the emmitter from the origin
        let dist_from_origo = (rx * rx + ry * ry + rz * rz).sqrt();
        // The sine of the maximum angle of the emitter
        let sin_alpha = bounding_radius / dist_from_origo;
        // The cossine of that angle
        // If the emitter is inside the bounding sphere of the detector,
        // then the maximum angle shall be 180 degrees, thus `cos_alpha` is -1.0
        let cos_alpha = if sin_alpha >= 1.0 {
            -1.0
        } else {
            (1.0 - sin_alpha * sin_alpha).sqrt()
        };
        Self {
            cos_alpha,
            // We save the middle direction of emission, since it will be a
            // good starting point for generating random directions
            mid_dir: (Vector::new(0., 0., 0.) - Vector::new(rx, ry, rz)).normalized(),
            solid_angle: 2.0 * std::f64::consts::PI * (1.0 - cos_alpha),
        }
    }

    /// Generates a random photon's direction towards the detector.
    pub fn gen_photon_dir(&self) -> Vector<f64> {
        // If we are inside the bounding sphere,
        // generate an isotropic direction
        if self.cos_alpha == -1.0 {
            return Vector::random_isotropic_normed();
        }
        // Else, generate a direction towards the detector
        let cos_theta = self.cos_alpha + (1.0 - self.cos_alpha) * f64::rand();
        let mut ret_photon = self.mid_dir.clone();
        ret_photon.rotate_random_by_angle_cosine(cos_theta);
        ret_photon
    }

    /// Returns the solid angle, in which the emission can happen.
    #[inline(always)]
    pub fn get_solid_angle(&self) -> f64 {
        self.solid_angle
    }
}
