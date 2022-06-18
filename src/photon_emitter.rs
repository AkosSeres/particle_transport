use crate::{rand_gen::RandGen, vec3::Vector};

pub struct PhotonEmitter {
    cos_alpha: f64,
    mid_dir: Vector<f64>,
}

impl PhotonEmitter {
    pub fn from_params(
        detector_r: f64,
        detector_h: f64,
        rx: f64,
        ry: f64,
        rz: f64,
    ) -> PhotonEmitter {
        let bounding_radius = (detector_h * detector_h / 4.0 + detector_r * detector_r).sqrt();
        let dist_from_origo = (rx * rx + ry * ry + rz * rz).sqrt();
        let sin_alpha = bounding_radius / dist_from_origo;
        Self {
            cos_alpha: (1.0 - sin_alpha * sin_alpha).sqrt(),
            mid_dir: (Vector::new(0., 0., 0.) - Vector::new(rx, ry, rz)).normalized(),
        }
    }

    pub fn gen_photon_dir(&self) -> Vector<f64> {
        let cos_theta = self.cos_alpha + (1.0 - self.cos_alpha) * f64::rand();
        let mut ret_photon = self.mid_dir.clone();
        ret_photon.rotate_random_by_angle_cosine(cos_theta);
        ret_photon
    }
}
