use crate::vec3::Vector;

/// The float type to use in calculations
type F = f64;

/// Radius of the detector
static mut DET_R: F = 2.5;
static mut DET_RSQ: F = 6.25;
/// Height of the detector
static mut DET_HEIGHT: F = 3.0;
static mut DET_ZBOTTOM: F = -1.5;
static mut DET_ZTOP: F = 1.5;
/// Density of the detector material, in g/cm³
static mut DET_DENSITY: F = 3.67;
/// Sets the properties of the detector
pub fn set_detector(r: F, height: F, density: F) {
    unsafe {
        DET_R = r;
        DET_RSQ = r * r;
        DET_HEIGHT = height;
        DET_ZBOTTOM = -0.5 * height;
        DET_ZTOP = 0.5 * height;
        DET_DENSITY = density;
    }
}
#[inline(always)]
fn get_detector_rsq() -> F {
    unsafe { DET_RSQ }
}
#[inline(always)]
fn get_detector_zbottom() -> F {
    unsafe { DET_ZBOTTOM }
}
#[inline(always)]
fn get_detector_ztop() -> F {
    unsafe { DET_ZTOP }
}

/// Describes what happens to the photon in the simulation step
pub enum PhotonStep {
    /// Completely misses the detector
    Miss,
    /// Enters the detector
    Enter,
    /// The photon transfers energy inside the detector
    TransferEnergy(F),
    /// Exits the detector
    Exit,
}

#[derive(PartialEq)]
enum PhotonLocation {
    OutsideMisses,
    OutsideInto(F),
    Inside(F),
}

#[derive(Debug, Copy, Clone)]
pub struct Photon {
    /// The current energy of the photon in keV
    energy: F,
    /// The current position of the photon
    pos: Vector<F>,
    /// The movement direction of the photon, expressed as a unit vector
    dir: Vector<F>,
}

impl Photon {
    pub fn step(&mut self) -> PhotonStep {
        PhotonStep::Miss
    }

    /// Intersects the path of the photon with an x-y plane and if they intersect, it returns the distance.
    /// # Arguments
    /// * `plane_z` - The z-level of the plane
    fn intersect_plane(&self, plane_z: F) -> Option<F> {
        if self.pos.z == plane_z {
            return Some(0.0);
        }
        if self.dir.z == 0.0 {
            return None;
        }
        Some((plane_z - self.pos.z) / self.dir.z)
    }

    /// Intersects the path of the photon with the bottom base of the detector and if they intersect, it returns the distance.
    fn intersect_bottom_base(&self) -> Option<F> {
        if self.pos.z == get_detector_zbottom() {
            return Some(0.0);
        }
        if self.dir.z == 0.0 {
            return None;
        }
        let dist = (get_detector_zbottom() - self.pos.z) / self.dir.z;
        let sect_point = self.pos + self.dir * dist;
        if sect_point.x * sect_point.x + sect_point.y * sect_point.y >= get_detector_rsq() {
            return None;
        }
        Some(dist)
    }
    /// Intersects the path of the photon with the top base of the detector and if they intersect, it returns the distance.
    fn intersect_top_base(&self) -> Option<F> {
        if self.pos.z == get_detector_ztop() {
            return Some(0.0);
        }
        if self.dir.z == 0.0 {
            return None;
        }
        let dist = (get_detector_ztop() - self.pos.z) / self.dir.z;
        let sect_point = self.pos + self.dir * dist;
        if sect_point.x * sect_point.x + sect_point.y * sect_point.y >= get_detector_rsq() {
            return None;
        }
        Some(dist)
    }

    /// Intersects the path of the photon with the infinitely extended cylinder of the detector and if they intersect, it returns the distances (both the entering distance and the exit distance as well).
    fn intersect_infinite_cylinder(&self) -> Option<(F, F)> {
        let a = self.dir.x * self.dir.x + self.dir.y * self.dir.y;
        let b = 2.0 * (self.pos.x * self.dir.x + self.pos.y * self.dir.y);
        let c = self.pos.x * self.pos.x + self.pos.y * self.pos.y - get_detector_rsq();
        let d = b * b - 4.0 * a * c;
        if d <= 0.0 {
            return None;
        }
        let diff = d.sqrt();
        Some(((-b - diff) / 2.0 / a, (-b + diff) / 2.0 / a))
    }

    /// Intersects the path of the photon with the cylinder of the detector and if they intersect, it returns the distances (both the entering distance and the exit distance as well).
    fn intersect_cylinder(&self) -> (Option<F>, Option<F>) {
        let a = self.dir.x * self.dir.x + self.dir.y * self.dir.y;
        let b = 2.0 * (self.pos.x * self.dir.x + self.pos.y * self.dir.y);
        let c = self.pos.x * self.pos.x + self.pos.y * self.pos.y - get_detector_rsq();
        let d = b * b - 4.0 * a * c;
        if d <= 0.0 {
            return (None, None);
        }
        let diff = d.sqrt();
        let sectz_0 = self.pos.z + self.dir.z * ((-b - diff) / 2.0 / a);
        let sectz_1 = self.pos.z + self.dir.z * ((-b + diff) / 2.0 / a);
        (
            if sectz_0 <= get_detector_zbottom() || sectz_0 >= get_detector_ztop() {
                None
            } else {
                Some(sectz_0)
            },
            if sectz_1 <= get_detector_zbottom() || sectz_1 >= get_detector_ztop() {
                None
            } else {
                Some(sectz_1)
            },
        )
    }

    fn intersect_detector(&self) -> PhotonLocation {
        let mut behind_dist: Option<F> = None;
        let mut forward_dist: Option<F> = None;
        let mut process_any_dist = |d_opt: Option<F>| match d_opt {
            None => {}
            Some(new_d) => {
                if new_d <= 0.0 {
                    if behind_dist.is_some() && behind_dist.unwrap() < new_d {
                        behind_dist = Some(new_d);
                    }
                } else {
                    if forward_dist.is_some() && forward_dist.unwrap() > new_d {
                        forward_dist = Some(new_d);
                    }
                }
            }
        };

        // Make intersections and process them
        process_any_dist(self.intersect_bottom_base());
        process_any_dist(self.intersect_top_base());
        let cylinder_dists = self.intersect_cylinder();
        process_any_dist(cylinder_dists.0);
        process_any_dist(cylinder_dists.1);

        if forward_dist.is_none() {
            return PhotonLocation::OutsideMisses;
        }
        if behind_dist.is_none() {
            return PhotonLocation::OutsideInto(forward_dist.unwrap());
        }
        PhotonLocation::Inside(forward_dist.unwrap())
    }

    fn sample_free_path(&self) -> F {
        // TODO
        0.5
    }

    fn sample_new_direction(&mut self) {
        // TODO
        self.dir = Vector::random_isotropic_normed();
    }

    fn move_by(&mut self, dist: F) {
        self.pos += self.dir * dist;
    }

    fn move_inside(&mut self, inside_wall_dist: F) -> (bool, F) {
        let free_path = self.sample_free_path();
        self.move_by(free_path);
        if free_path < inside_wall_dist {
            self.sample_new_direction();
            return (true, 30.0);
        } else {
            return (false, 0.0);
        }
    }

    /// Simulates the path of the photon, and we return how much energy it gives off to the detector
    fn simulate(&mut self) -> F {
        let mut location = self.intersect_detector();
        let mut energy_transfered = 0.0;
        while location != PhotonLocation::OutsideMisses {
            match location {
                PhotonLocation::OutsideInto(dist) => self.move_by(dist),
                PhotonLocation::Inside(dist) => {
                    let move_result = self.move_inside(dist);
                    energy_transfered += move_result.1;
                    if !move_result.0 {
                        return energy_transfered;
                    }
                }
                PhotonLocation::OutsideMisses => return energy_transfered,
            }
            location = self.intersect_detector();
        }
        0.0
    }
}

// Energy dependent cross sections
const ENERGIES: [F; 58] = [
    1.0, 1.035, 1.072, 1.072, 1.5, 2.0, 3.0, 4.0, 4.557, 4.557, 4.702, 4.852, 4.852, 5.0, 5.188,
    5.188, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 33.17, 33.17, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0,
    200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1022.0, 1250.0, 1500.0, 2000.0, 2044.0,
    3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 11000.0, 12000.0, 13000.0,
    14000.0, 15000.0, 16000.0, 18000.0, 20000.0,
];
const COMPTON: [F; 58] = [
    0.005916, 0.006244, 0.006586, 0.006587, 0.01065, 0.01541, 0.02459, 0.03292, 0.03712, 0.03712,
    0.03816, 0.03923, 0.03923, 0.04026, 0.04155, 0.04155, 0.04679, 0.05799, 0.06731, 0.0841,
    0.09463, 0.1067, 0.109, 0.109, 0.1127, 0.1157, 0.1169, 0.1165, 0.1144, 0.1071, 0.1, 0.0886,
    0.08009, 0.07351, 0.06822, 0.06012, 0.05413, 0.05355, 0.04846, 0.04407, 0.03764, 0.03717,
    0.02963, 0.02472, 0.02135, 0.01887, 0.01697, 0.01544, 0.0142, 0.01315, 0.01227, 0.0115,
    0.01084, 0.01025, 0.009737, 0.009271, 0.008471, 0.007811,
];
const FOTOEFFECT: [F; 58] = [
    7794.0, 7245.0, 6736.0, 7924.0, 3801.0, 1917.0, 700.3, 335.1, 238.7, 658.5, 616.6, 577.4,
    771.9, 727.7, 661.1, 760.4, 529.7, 248.9, 137.5, 45.7, 20.62, 6.607, 4.972, 29.76, 18.24,
    10.06, 6.111, 2.746, 1.462, 0.4592, 0.2019, 0.06481, 0.02987, 0.01684, 0.01079, 0.005588,
    0.003491, 0.003323, 0.00224, 0.001604, 0.0009761, 0.0009415, 0.0005178, 0.0003428, 0.0002535,
    0.0002, 0.0001647, 0.0001398, 0.0001212, 0.0001069, 0.00009565, 0.00008648, 0.00007888,
    0.00007249, 0.00006706, 0.0000624, 0.00005472, 0.00004873,
];
const PAIR_PRODUCTION: [F; 58] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0001634, 0.0007594, 0.002575, 0.002751, 0.00652836, 0.01005228, 0.01311417, 0.0157591,
    0.0181239, 0.0202573, 0.0221889, 0.0239583, 0.0255956, 0.0271008, 0.028514, 0.0298055,
    0.0310054, 0.0321239, 0.0341662, 0.0359942,
];
const SUM_OF_CROSS_SECTIONS: [F; 58] = [
    7794.005916,
    7245.006244,
    6736.006586,
    7924.006587,
    3801.01065,
    1917.01541,
    700.32459,
    335.13292,
    238.73712,
    658.53712,
    616.63816,
    577.43923,
    771.93923,
    727.74026,
    661.14155,
    760.44155,
    529.74679,
    248.95799,
    137.56731,
    45.7841,
    20.71463,
    6.7137,
    5.081,
    29.869,
    18.3527,
    10.1757,
    6.2279,
    2.8625,
    1.5764,
    0.5663,
    0.3019,
    0.15341,
    0.10996,
    0.09035,
    0.07901,
    0.065708,
    0.057621,
    0.056873,
    0.0508634,
    0.0464334,
    0.0411911,
    0.0408625,
    0.03667616,
    0.03511508,
    0.03471767,
    0.0348291,
    0.0352586,
    0.0358371,
    0.0365101,
    0.0372152,
    0.03796125,
    0.03868728,
    0.03943288,
    0.04012799,
    0.04080946,
    0.0414573,
    0.04269192,
    0.04385393,
];
