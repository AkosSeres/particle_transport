/// Struct for representing a 3D vector in Cartesian coordinates.
#[derive(Debug, Copy, Clone)]
pub struct Vector<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

/// Struct for representing a 3D vector in spherical coordinates.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Spherical<T> {
    /// The magnitude (radius) of the vector.
    pub radius: T,
    /// Angle of inclination
    pub theta: T,
    /// Azimuthal angle
    pub phi: T,
}

impl<T> Vector<T> {
    /// Creates a new Vector from the given Cartesian coordinates.
    pub fn new(x_: T, y_: T, z_: T) -> Vector<T> {
        Vector::<T> {
            x: x_,
            y: y_,
            z: z_,
        }
    }

    /// Creates a 3D vector from the given spherical coordinates.
    /// # Arguments
    /// * `radius` - The radius (magnitude) of the vector, in range \[0, ∞\]
    /// * `theta` - Angle of inclination, in the range \[0, 𝜋\]
    /// * `phi` - Azimuth angle, in the range \[0, 2𝜋\]
    /// # Comments
    /// * `radius` is allowed to be negative, but if it is, then the vector will face in the opposite direction
    pub fn from_spherical(radius: T, theta: T, phi: T) -> Vector<T>
    where
        T: num::traits::Float,
    {
        let sintheta = theta.sin();
        Vector {
            x: phi.cos() * sintheta,
            y: phi.sin() * sintheta,
            z: theta.cos(),
        } * radius
    }

    /// Returns a [Spherical] with the spherical coordinates of the vector.
    pub fn to_spherical(&self) -> Spherical<T>
    where
        T: num::traits::Float,
    {
        let r = self.mag();
        Spherical::<T> {
            radius: r,
            theta: (self.z / r).acos(),
            phi: self.y.atan2(self.x),
        }
    }

    /// Dot product of two vectors.
    pub fn dot(&self, other: &Self) -> T
    where
        T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
    {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// The angle between two vectors.
    pub fn angle(&self, other: &Self) -> T
    where
        T: num::traits::Float,
    {
        (self.dot(other) / (self.mag() * other.mag())).acos()
    }

    /// The square of the magnitude of the vector.
    /// This is useful for comparing the magnitude of two vectors,
    /// because we can avoid the square root function.
    pub fn mag_sq(&self) -> T
    where
        T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns the magnitude of the vector.
    pub fn mag(&self) -> T
    where
        T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + num::traits::Float,
    {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Sets the magnitude of the vector to the given value.
    /// # Arguments
    /// * `new_mag` - The new magnitude of the vector
    /// # Comments
    /// * `new_mag` can be negative, in which case the vector will face in the opposite direction.
    pub fn set_mag(&mut self, new_mag: T)
    where
        T: num::traits::Float,
    {
        *self = *self * (new_mag / self.mag());
    }

    /// Sets the magnitude of the vector to 1.
    pub fn normalize(&mut self)
    where
        T: num::traits::Float,
    {
        *self = *self / self.mag();
    }

    /// Returns a [normalize][Self::normalize]d copy of the vector.
    pub fn normalized(&self) -> Self
    where
        T: num::traits::Float,
    {
        *self / self.mag()
    }

    /// Returns a [normalize][Self::normalize]d random vector,
    /// with _isotropic_ distribution. In other words, we get a
    /// uniform random point on the surface of the unit sphere.
    ///
    /// The vector is generated using the Marsaglia method.
    pub fn random_isotropic_normed() -> Vector<T>
    where
        T: crate::rand_gen::RandGen + num::traits::Float,
        f64: Into<T>,
        f32: Into<T>,
    {
        let mut rhosq: T;
        let mut u: T;
        let mut v: T;
        while {
            u = T::rand() * 2.0.into() - 1.0.into();
            v = T::rand() * 2.0.into() - 1.0.into();
            rhosq = u * u + v * v;
            rhosq > 1.0.into()
        } {}
        let sqrt_part = (-rhosq + 1.0.into()).sqrt();
        Vector::<T> {
            x: 1.0.into() - 2.0.into() * rhosq,
            y: 2.0.into() * u * sqrt_part,
            z: 2.0.into() * v * sqrt_part,
        }
    }

    /// Rotates the vector by the given angle, but this rotation happens
    /// around a random axis that is perpendicular to the vector.
    pub fn rotate_random_by_angle(&mut self, angle: T)
    where
        T: crate::rand_gen::RandGen + num::traits::Float,
        f64: Into<T>,
        f32: Into<T>,
    {
        let theta = angle;
        let phi = T::rand() * 3.141592653589793.into() * 2.0.into();
        let (u, v, w) = (self.x, self.y, self.z);
        let thetasin = theta.sin();
        let sincos = thetasin * phi.cos();
        let sinsin = thetasin * phi.sin();
        let coss = theta.cos();
        let sqpart = (u * u + v * v).sqrt();
        let new_v = Self {
            x: sincos * v / sqpart + u * w * sinsin / sqpart + u * coss,
            y: -u * sincos / sqpart + v * w * sinsin / sqpart + v * coss,
            z: -sqpart * sinsin + w * coss,
        };
        *self = new_v;
    }

    /// The same as the [rotate_random_by_angle][Self::rotate_random_by_angle] function,
    /// but here we give the function the _cosine_ of the angle, so the execution time is
    /// slightly faster.
    pub fn rotate_random_by_angle_cosine(&mut self, angle_cosine: T)
    where
        T: crate::rand_gen::RandGen + num::traits::Float,
        f64: Into<T>,
        f32: Into<T>,
    {
        let theta_cos = angle_cosine;
        let phi = T::rand() * 3.141592653589793.into() * 2.0.into();
        let (u, v, w) = (self.x, self.y, self.z);
        let thetasin = (-theta_cos * theta_cos + 1.0.into()).sqrt();
        let sincos = thetasin * phi.cos();
        let sinsin = thetasin * phi.sin();
        let coss = theta_cos;
        let sqpart = (u * u + v * v).sqrt();
        let new_v = Self {
            x: sincos * v / sqpart + u * w * sinsin / sqpart + u * coss,
            y: -u * sincos / sqpart + v * w * sinsin / sqpart + v * coss,
            z: -sqpart * sinsin + w * coss,
        };
        *self = new_v;
    }
}

impl<T> std::ops::Add for Vector<T>
where
    T: std::ops::Add<Output = T>,
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T> std::ops::AddAssign for Vector<T>
where
    T: std::ops::AddAssign,
{
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T> std::ops::Sub for Vector<T>
where
    T: std::ops::Sub<Output = T>,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T> std::ops::SubAssign for Vector<T>
where
    T: std::ops::SubAssign,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl<T> std::cmp::PartialEq<Vector<T>> for Vector<T>
where
    T: std::cmp::PartialEq,
{
    fn eq(&self, other: &Vector<T>) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}

impl<T> std::ops::Div<T> for Vector<T>
where
    T: std::ops::Div<Output = T> + core::marker::Copy,
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T> std::ops::Mul<T> for Vector<T>
where
    T: std::ops::Mul<Output = T> + core::marker::Copy,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> std::ops::Neg for Vector<T>
where
    T: std::ops::Neg<Output = T>,
{
    type Output = Vector<T>;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

/// Vector of [f32] coordinates.
pub type Vectorf = Vector<f32>;
/// Vector of [f64] coordinates.
pub type Vectord = Vector<f64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equals() {
        assert_eq!(
            Vector::<f32> {
                x: 6.75,
                y: 15.825,
                z: 32.5
            },
            Vector::<f32> {
                x: 6.75,
                y: 15.825,
                z: 32.5
            }
        );
    }

    #[test]
    fn create() {
        assert_eq!(
            Vector::<f32> {
                x: 0.34,
                y: 2.5,
                z: 9.25
            },
            Vector::<f32>::new(0.34, 2.5, 9.25)
        );
    }
}
