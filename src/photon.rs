use crate::vec3::Vector;

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
