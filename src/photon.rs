use crate::{rand_gen::RandGen, vec3::Vector};

/// Radius of the detector
static mut DET_R: f64 = 2.5;
static mut DET_RSQ: f64 = 6.25;
/// Height of the detector
static mut DET_HEIGHT: f64 = 3.0;
static mut DET_ZBOTTOM: f64 = -1.5;
static mut DET_ZTOP: f64 = 1.5;
/// Density of the detector material, in g/cm³
static mut DET_DENSITY: f64 = 3.67;
/// Default photon energy
static mut DEFAULT_ENERGY: f64 = ENERGIES[0];
/// Default cross sections corresponding to the default [DEFAULT_ENERGY]
static mut DEFAULT_CROSS_SECTIONS: (f64, f64, f64, f64) = (
    COMPTON[0],
    FOTOEFFECT[0],
    PAIR_PRODUCTION[0],
    SUM_OF_CROSS_SECTIONS[0],
);
/// Sets the properties of the detector
pub fn set_detector(r: f64, height: f64, density: f64) {
    unsafe {
        DET_R = r;
        DET_RSQ = r * r;
        DET_HEIGHT = height;
        DET_ZBOTTOM = -0.5 * height;
        DET_ZTOP = 0.5 * height;
        DET_DENSITY = density;
    }
}
pub fn set_default_energy(energy: f64) {
    unsafe {
        DEFAULT_ENERGY = energy;
        DEFAULT_CROSS_SECTIONS = Photon::get_cross_sections(energy);
    }
}
#[inline(always)]
fn get_detector_rsq() -> f64 {
    unsafe { DET_RSQ }
}
#[inline(always)]
fn get_detector_zbottom() -> f64 {
    unsafe { DET_ZBOTTOM }
}
#[inline(always)]
fn get_detector_ztop() -> f64 {
    unsafe { DET_ZTOP }
}
#[inline(always)]
fn get_detector_density() -> f64 {
    unsafe { DET_DENSITY }
}

/// Describes what happens to the photon in the simulation step
pub enum PhotonStep {
    /// Completely misses the detector
    Miss,
    /// Enters the detector
    Enter,
    /// The photon transfers energy inside the detector
    TransferEnergy(f64),
    /// Exits the detector
    Exit,
}

enum InteractionType {
    ComptonScatter,
    Photoeffect,
    PairProduction,
}

#[derive(PartialEq, Debug)]
pub enum PhotonLocation {
    OutsideMisses,
    OutsideInto(f64),
    Inside(f64),
}

#[derive(Debug, Copy, Clone)]
pub struct Photon {
    /// The current energy of the photon in keV
    pub energy: f64,
    /// The current position of the photon
    pub pos: Vector<f64>,
    /// The movement direction of the photon, expressed as a unit vector
    pub dir: Vector<f64>,
}

impl Photon {
    /// Intersects the path of the photon with an x-y plane and if they intersect, it returns the distance.
    /// # Arguments
    /// * `plane_z` - The z-level of the plane
    fn intersect_plane(&self, plane_z: f64) -> Option<f64> {
        if self.pos.z == plane_z {
            return Some(0.0);
        }
        if self.dir.z == 0.0 {
            return None;
        }
        Some((plane_z - self.pos.z) / self.dir.z)
    }

    /// Intersects the path of the photon with the bottom base of the detector and if they intersect, it returns the distance.
    pub fn intersect_bottom_base(&self) -> Option<f64> {
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
    pub fn intersect_top_base(&self) -> Option<f64> {
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
    pub fn intersect_infinite_cylinder(&self) -> Option<(f64, f64)> {
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

    /// Intersects the path of the photon with the cylinder of the detector and if they intersect, it returns the frontal distance
    pub fn intersect_cylinder(&self) -> Option<f64> {
        let a = self.dir.x * self.dir.x + self.dir.y * self.dir.y;
        let b = 2.0 * (self.pos.x * self.dir.x + self.pos.y * self.dir.y);
        let c = self.pos.x * self.pos.x + self.pos.y * self.pos.y - get_detector_rsq();
        let d = b * b - 4.0 * a * c;
        if d <= 0.0 {
            return None;
        }
        let diff = d.sqrt();
        let dist0 = (-b - diff) / 2.0 / a;
        let dist1 = (-b + diff) / 2.0 / a;
        let sectz_0 = self.pos.z + self.dir.z * dist0;
        let sectz_1 = self.pos.z + self.dir.z * dist1;
        let mut ret_dist = None;
        if dist0 > 0.0 && sectz_0 > get_detector_zbottom() && sectz_0 < get_detector_ztop() {
            ret_dist = Some(dist0);
        }
        if dist1 > 0.0 && sectz_1 > get_detector_zbottom() && sectz_1 < get_detector_ztop() {
            if ret_dist.is_none() {
                ret_dist = Some(dist1);
            } else {
                unsafe {
                    if ret_dist.unwrap_unchecked() > dist1 {
                        ret_dist = Some(dist1);
                    }
                }
            }
        }
        ret_dist
    }

    pub fn intersect_detector(&self) -> PhotonLocation {
        let is_inside = self.pos.x * self.pos.x + self.pos.y * self.pos.y < get_detector_rsq()
            && self.pos.z > get_detector_zbottom()
            && self.pos.z < get_detector_ztop();
        let mut forward_dist: Option<f64> = None;
        let mut process_any_dist = |d_opt: Option<f64>| match d_opt {
            None => {}
            Some(new_d) => {
                if new_d > 0.0 {
                    if forward_dist.is_none() {
                        forward_dist = Some(new_d);
                    } else {
                        unsafe {
                            if forward_dist.unwrap_unchecked() > new_d {
                                forward_dist = Some(new_d);
                            }
                        }
                    }
                }
            }
        };

        // Make intersections and process them
        process_any_dist(self.intersect_bottom_base());
        process_any_dist(self.intersect_top_base());
        process_any_dist(self.intersect_cylinder());

        if forward_dist.is_some() {
            if is_inside {
                unsafe {
                    return PhotonLocation::Inside(forward_dist.unwrap_unchecked());
                }
            } else {
                unsafe {
                    return PhotonLocation::OutsideInto(forward_dist.unwrap_unchecked());
                }
            }
        }
        return PhotonLocation::OutsideMisses;
    }

    fn sample_free_path_and_interaction_type(&self) -> (f64, InteractionType)
    where
        f64: RandGen,
    {
        let cross_sections = Photon::get_cross_sections_cached(self.energy);
        let cross_section_sum_macroscopic = cross_sections.3 * get_detector_density();
        let free_path = -(1.0 / cross_section_sum_macroscopic) * f64::rand().ln();

        let interaction_coeff = cross_sections.3 * f64::rand();
        let interaction_type = if interaction_coeff <= cross_sections.0 {
            InteractionType::ComptonScatter
        } else if interaction_coeff <= (cross_sections.1 + cross_sections.0) {
            InteractionType::Photoeffect
        } else {
            InteractionType::PairProduction
        };
        (free_path, interaction_type)
    }

    /// Returns all of the cross sections corresponding to the given energy.
    /// The energy must be given in keV, and the tuple conatins the following values in this order:
    /// (cross section of Compton scattering, cross section for fotoeffect, cross section for pair production, sum of all cross sections)
    ///
    /// The cross sections are given in cm²/g
    fn get_cross_sections(energy: f64) -> (f64, f64, f64, f64) {
        let mut idx_low: usize = 0;
        let mut idx_high: usize = CROSS_SECTION_DATA_POINTS - 1;
        if energy <= ENERGIES[0] {
            return (
                COMPTON[idx_low],
                FOTOEFFECT[idx_low],
                PAIR_PRODUCTION[idx_low],
                SUM_OF_CROSS_SECTIONS[idx_low],
            );
        }
        if energy >= ENERGIES[CROSS_SECTION_DATA_POINTS - 1] {
            return (
                COMPTON[idx_high],
                FOTOEFFECT[idx_high],
                PAIR_PRODUCTION[idx_high],
                SUM_OF_CROSS_SECTIONS[idx_high],
            );
        }
        while idx_high - idx_low > 1 {
            let idx_mid = (idx_high + idx_low) / 2;
            let energy_mid = ENERGIES[idx_mid];
            if energy >= energy_mid {
                idx_low = idx_mid;
            } else {
                idx_high = idx_mid;
            }
        }
        if idx_high == idx_low {
            return (
                COMPTON[idx_high],
                FOTOEFFECT[idx_high],
                PAIR_PRODUCTION[idx_high],
                SUM_OF_CROSS_SECTIONS[idx_high],
            );
        }
        let energy_low = ENERGIES[idx_low];
        let energy_high = ENERGIES[idx_high];
        let hcoeff = (energy - energy_low) / (energy_high - energy_low);
        let lcoeff = 1.0 - hcoeff;
        (
            COMPTON[idx_low] * lcoeff + COMPTON[idx_high] * hcoeff,
            FOTOEFFECT[idx_low] * lcoeff + FOTOEFFECT[idx_high] * hcoeff,
            PAIR_PRODUCTION[idx_low] * lcoeff + PAIR_PRODUCTION[idx_high] * hcoeff,
            SUM_OF_CROSS_SECTIONS[idx_low] * lcoeff + SUM_OF_CROSS_SECTIONS[idx_high] * hcoeff,
        )
    }

    /// Returns all of the cross sections corresponding to the given energy, but this function caches these values
    /// for the default energy value, unlike _get_cross_sections_
    /// The energy must be given in keV, and the tuple conatins the following values in this order:
    /// (cross section of Compton scattering, cross section for fotoeffect, cross section for pair production, sum of all cross sections)
    ///
    /// The cross sections are given in cm²/g
    fn get_cross_sections_cached(energy: f64) -> (f64, f64, f64, f64) {
        unsafe {
            if energy == DEFAULT_ENERGY {
                return DEFAULT_CROSS_SECTIONS;
            }
            // Pair producted photon energy
            if energy == 511.0 {
                return (0.07292810000000001, 0.0161745, 0.0, 0.0891026);
            }
        }
        Self::get_cross_sections(energy)
    }

    fn move_by(&mut self, dist: f64) {
        self.pos = self.pos + self.dir * dist;
    }

    fn move_inside(&mut self, inside_wall_dist: f64) -> (bool, f64) {
        let next_move = self.sample_free_path_and_interaction_type();
        let interaction_type = next_move.1;
        let free_path = next_move.0;
        self.move_by(free_path);
        if free_path < inside_wall_dist {
            match interaction_type {
                InteractionType::ComptonScatter => {
                    // Kahn method
                    let lambda = 511.0 / self.energy;
                    let R;
                    loop {
                        let (r1, r2, r3) = (f64::rand(), f64::rand(), f64::rand());
                        if r1 <= (1.0 + 2.0 / lambda) / (9.0 + 2.0 / lambda) {
                            let R_ = 1.0 + 2.0 * r2 / lambda;
                            if r3 <= 4.0 * (1.0 / R_ + 1.0 / (R_ * R_)) {
                                R = R_;
                                break;
                            } else {
                                continue;
                            }
                        } else {
                            let R_ = (1.0 + 2.0 / lambda) / (1.0 + 2.0 * r2 / lambda);
                            if r3 <= 0.5 * ((lambda - R_ * lambda + 1.0).powi(2) + (1.0 / R_)) {
                                R = R_;
                                break;
                            } else {
                                continue;
                            }
                        }
                    }
                    let costheta = 1.0 + lambda - lambda * R;
                    self.dir.rotate_random_by_angle(costheta.acos());
                    let post_energy = self.energy / R;
                    let energy_d = self.energy - post_energy;
                    self.energy = post_energy;
                    //(self.energy > 30.0, 30.0)
                    (true, energy_d)
                }
                InteractionType::Photoeffect => {
                    // The photon gets entirely absorbed
                    (false, self.energy)
                }
                InteractionType::PairProduction => {
                    let mut photon1 = Photon {
                        energy: 511.0,
                        pos: self.pos,
                        dir: Vector::<f64>::random_isotropic_normed(),
                    };
                    let mut photon2 = Photon {
                        energy: 511.0,
                        pos: self.pos,
                        dir: -photon1.dir,
                    };

                    (
                        false,
                        self.energy - 1022.0
                            + photon1.simulate().energy_transfered
                            + photon2.simulate().energy_transfered,
                    )
                }
            }
        } else {
            return (false, 0.0);
        }
    }

    /// Simulates the path of the photon, and we return how much energy it gives off to the detector
    pub fn simulate(&mut self) -> SimulationResults {
        let mut location = self.intersect_detector();
        let start_energy = self.energy;
        let mut energy_transfered = 0.0;

        let mut results = SimulationResults {
            energy_transfered: 0.0,
            hit_detector: false,
        };

        while location != PhotonLocation::OutsideMisses {
            results.hit_detector = true;
            match location {
                PhotonLocation::OutsideInto(dist) => self.move_by(dist),
                PhotonLocation::Inside(dist) => {
                    let move_result = self.move_inside(dist);
                    results.energy_transfered += move_result.1;
                    if !move_result.0 {
                        return results;
                    }
                }
                _ => return results,
            }
            location = self.intersect_detector();
        }
        results
    }
}

pub struct SimulationResults {
    pub energy_transfered: f64,
    pub hit_detector: bool,
}

const CROSS_SECTION_DATA_POINTS: usize = 58;
// Energy dependent cross sections
const ENERGIES: [f64; CROSS_SECTION_DATA_POINTS] = [
    1.0, 1.035, 1.072, 1.072, 1.5, 2.0, 3.0, 4.0, 4.557, 4.557, 4.702, 4.852, 4.852, 5.0, 5.188,
    5.188, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 33.17, 33.17, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0,
    200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1022.0, 1250.0, 1500.0, 2000.0, 2044.0,
    3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 11000.0, 12000.0, 13000.0,
    14000.0, 15000.0, 16000.0, 18000.0, 20000.0,
];
const COMPTON: [f64; CROSS_SECTION_DATA_POINTS] = [
    0.005916, 0.006244, 0.006586, 0.006587, 0.01065, 0.01541, 0.02459, 0.03292, 0.03712, 0.03712,
    0.03816, 0.03923, 0.03923, 0.04026, 0.04155, 0.04155, 0.04679, 0.05799, 0.06731, 0.0841,
    0.09463, 0.1067, 0.109, 0.109, 0.1127, 0.1157, 0.1169, 0.1165, 0.1144, 0.1071, 0.1, 0.0886,
    0.08009, 0.07351, 0.06822, 0.06012, 0.05413, 0.05355, 0.04846, 0.04407, 0.03764, 0.03717,
    0.02963, 0.02472, 0.02135, 0.01887, 0.01697, 0.01544, 0.0142, 0.01315, 0.01227, 0.0115,
    0.01084, 0.01025, 0.009737, 0.009271, 0.008471, 0.007811,
];
const FOTOEFFECT: [f64; CROSS_SECTION_DATA_POINTS] = [
    7794.0, 7245.0, 6736.0, 7924.0, 3801.0, 1917.0, 700.3, 335.1, 238.7, 658.5, 616.6, 577.4,
    771.9, 727.7, 661.1, 760.4, 529.7, 248.9, 137.5, 45.7, 20.62, 6.607, 4.972, 29.76, 18.24,
    10.06, 6.111, 2.746, 1.462, 0.4592, 0.2019, 0.06481, 0.02987, 0.01684, 0.01079, 0.005588,
    0.003491, 0.003323, 0.00224, 0.001604, 0.0009761, 0.0009415, 0.0005178, 0.0003428, 0.0002535,
    0.0002, 0.0001647, 0.0001398, 0.0001212, 0.0001069, 0.00009565, 0.00008648, 0.00007888,
    0.00007249, 0.00006706, 0.0000624, 0.00005472, 0.00004873,
];
const PAIR_PRODUCTION: [f64; CROSS_SECTION_DATA_POINTS] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0001634, 0.0007594, 0.002575, 0.002751, 0.00652836, 0.01005228, 0.01311417, 0.0157591,
    0.0181239, 0.0202573, 0.0221889, 0.0239583, 0.0255956, 0.0271008, 0.028514, 0.0298055,
    0.0310054, 0.0321239, 0.0341662, 0.0359942,
];
const SUM_OF_CROSS_SECTIONS: [f64; CROSS_SECTION_DATA_POINTS] = [
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
