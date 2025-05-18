/*
Copyright 2023 Norbert Fröhlich


Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


use crate::proportion::TwoProportion;
use num::{Zero,One};
use std::ops::{Mul, Add, Sub, Div};
use crate::matrix::{Matrix, RowVector, ColumnVector};


/// Represents a projective 1D point
#[derive(Debug, Clone)]
pub struct ProjOnePoint<T> {
    pub x : TwoProportion<T>,
}


impl<T> ProjOnePoint<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(x: T, y: T) -> ProjOnePoint<T> {
        ProjOnePoint { x: TwoProportion::new(x,y) }
    }
}


impl<T> ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero + One,
    T: Clone,
{
    fn metric_blue() -> Matrix<T> {
        Matrix::new(vec![vec![T::one(), T::zero()],
                    vec![T::zero(), T::one()]])
    }

    fn metric_red() -> Matrix<T> {
        Matrix::new(vec![vec![T::one(), T::zero()],
                         vec![T::zero(), T::zero() - T::one()]])
    }

    fn metric_green() -> Matrix<T> {
        Matrix::new(vec![vec![T::zero(), T::one()],
                    vec![T::one(), T::zero()]])
    }

    pub fn dot_blue(&self, other: &Self) -> T {
        self.dot_metric(&other, &ProjOnePoint::metric_blue())
    }

    pub fn dot_red(&self, other: &Self) -> T {
        self.dot_metric(&other, &ProjOnePoint::metric_red())
    }

    pub fn dot_green(&self, other: &Self) -> T {
        self.dot_metric(&other, &ProjOnePoint::metric_green())
    }

    pub fn dot_metric(&self, other: &Self, metric: &Matrix<T>) -> T {
        let a = self.x.a.clone();
        let b = self.x.b.clone();
        let c = other.x.a.clone();
        let d = other.x.b.clone();

        let v1 = RowVector::new(vec![a, b]);
        let v2 = ColumnVector::new(vec![c, d]);

        v1*(*metric).clone()*v2
    }

    pub fn cross(&self, other: &Self) -> T {
        let a = self.x.a.clone();
        let b = self.x.b.clone();
        let c = other.x.a.clone();
        let d = other.x.b.clone();

        let v1 = RowVector::new(vec![a, b]);
        let m1 = Matrix::new(vec![vec![T::zero(), T::one()],
                                  vec![T::zero() - T::one(), T::zero()]]);
        let v2 = ColumnVector::new(vec![c, d]);

        v1*m1*v2
    }

    pub fn quadrance(&self, other: &Self) -> T {
        let n = self.cross(&other);
        let numer = n.clone() * n;
        let d1 = self.dot_blue(&self);
        let d2 = other.dot_blue(&other);
        return numer / (d1*d2);
    }

    pub fn quadrance_red(&self, other: &Self) -> T {
        let n = self.cross(&other);
        let numer = T::zero() - n.clone() * n;
        let d1 = self.dot_red(&self);
        let d2 = other.dot_red(&other);
        return numer / (d1*d2);
    }
    pub fn quadrance_green(&self, other: &Self) -> T {
        let n = self.cross(&other);
        let numer = T::zero() - n.clone() * n;
        let d1 = self.dot_green(&self);
        let d2 = other.dot_green(&other);
        return numer / (d1*d2);
    }

}

impl<T> ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    pub fn is_perpendicular(&self, other: &Self) -> bool {
        let x1 = self.x.a.clone();
        let y1 = self.x.b.clone();
        let x2 = other.x.a.clone();
        let y2 = other.x.b.clone();
        let r = x1*x2+y1*y2;
        r.is_zero()
    }
    pub fn perpendicular(&self) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        ProjOnePoint::new(T::zero() - y, x)
    }
    pub fn is_perpendicular_red(&self, other: &Self) -> bool {
        let x1 = self.x.a.clone();
        let y1 = self.x.b.clone();
        let x2 = other.x.a.clone();
        let y2 = other.x.b.clone();
        let r = x1*x2-y1*y2;
        r.is_zero()
    }
    pub fn perpendicular_red(&self) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        ProjOnePoint::new(y, x)
    }
    pub fn is_perpendicular_green(&self, other: &Self) -> bool {
        let x1 = self.x.a.clone();
        let y1 = self.x.b.clone();
        let x2 = other.x.a.clone();
        let y2 = other.x.b.clone();
        let r = x1*y2+x2*y1;
        r.is_zero()
    }
    pub fn perpendicular_green(&self) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        ProjOnePoint::new(x, T::zero() - y)
    }
}

impl<T> PartialEq for ProjOnePoint<T>
where
    TwoProportion<T>: PartialEq,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x
    }
}


/// Rotation
#[derive(Debug,Clone)]
pub struct Rotation<T> {
    vector: TwoProportion<T>,
}

impl<T> Rotation<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> Rotation<T> {
        Rotation { vector: TwoProportion::new(a, b) }
    }
}

impl<T> Mul<Rotation<T>> for ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: Rotation<T>) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        let a = r.vector.a.clone();
        let b = r.vector.b.clone();
        let rx = a.clone()*x.clone() - b.clone()*y.clone();
        let ry = b*x + a*y;
        ProjOnePoint::new(rx,ry)
    }
}

impl<T> Mul<Rotation<T>> for Rotation<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = Rotation<T>;

    fn mul(self, other: Self) -> Rotation<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let ra = a.clone()*c.clone() - b.clone()*d.clone();
        let rb = a*d + b*c;
        Rotation::new(ra,rb)
    }
}

impl<T> One for Rotation<T>
where
    T: Sub<T, Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn one() -> Rotation<T> {
        Rotation::new(T::one(),T::zero())
    }
}

/// Reflection
#[derive(Debug,Clone)]
pub struct Reflection<T> {
    vector: TwoProportion<T>,
}

impl<T> Reflection<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> Reflection<T> {
        Reflection { vector: TwoProportion::new(a, b) }
    }
}

impl<T> Mul<Reflection<T>> for ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: Reflection<T>) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        let a = r.vector.a.clone();
        let b = r.vector.b.clone();
        let rx = a.clone()*x.clone() + b.clone()*y.clone();
        let ry = b*x - a*y;
        ProjOnePoint::new(rx,ry)
    }
}

impl<T> Mul<Reflection<T>> for Reflection<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = Rotation<T>;

    fn mul(self, other: Self) -> Rotation<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let ra = a.clone()*c.clone() + b.clone()*d.clone();
        let rb = a*d - b*c;
        Rotation::new(ra,rb)
    }
}

impl<T> Mul<Reflection<T>> for Rotation<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = Reflection<T>;

    fn mul(self, other: Reflection<T>) -> Reflection<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let sa = a.clone()*c.clone() + b.clone()*d.clone();
        let sb = a*d - b*c;
        Reflection::new(sa,sb)
    }
}

impl<T> Mul<Rotation<T>> for Reflection<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = Reflection<T>;

    fn mul(self, other: Rotation<T>) -> Reflection<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let sa = a.clone()*c.clone() - b.clone()*d.clone();
        let sb = a*d + b*c;
        Reflection::new(sa,sb)
    }
}

/// red Rotation
#[derive(Debug,Clone)]
pub struct RotationRed<T> {
    vector: TwoProportion<T>,
}

impl<T> RotationRed<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> RotationRed<T> {
        RotationRed { vector: TwoProportion::new(a, b) }
    }
}

impl<T> Mul<RotationRed<T>> for ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: RotationRed<T>) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        let a = r.vector.a.clone();
        let b = r.vector.b.clone();
        let rx = a.clone()*x.clone() + b.clone()*y.clone();
        let ry = b*x + a*y;
        ProjOnePoint::new(rx,ry)
    }
}

impl<T> Mul<RotationRed<T>> for RotationRed<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = RotationRed<T>;

    fn mul(self, other: Self) -> RotationRed<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let ra = a.clone()*c.clone() + b.clone()*d.clone();
        let rb = a*d + b*c;
        RotationRed::new(ra,rb)
    }
}

impl<T> One for RotationRed<T>
where
    T: Sub<T, Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn one() -> RotationRed<T> {
        RotationRed::new(T::one(),T::zero())
    }
}


/// red Reflection
#[derive(Debug,Clone)]
pub struct ReflectionRed<T> {
    vector: TwoProportion<T>,
}

impl<T> ReflectionRed<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> ReflectionRed<T> {
        ReflectionRed { vector: TwoProportion::new(a, b) }
    }
}

impl<T> Mul<ReflectionRed<T>> for ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: ReflectionRed<T>) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        let a = r.vector.a.clone();
        let b = r.vector.b.clone();
        let rx = a.clone()*x.clone() - b.clone()*y.clone();
        let ry = b*x - a*y;
        ProjOnePoint::new(rx,ry)
    }
}

impl<T> Mul<ReflectionRed<T>> for ReflectionRed<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = RotationRed<T>;

    fn mul(self, other: Self) -> RotationRed<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let ra = a.clone()*c.clone() - b.clone()*d.clone();
        let rb = a*d - b*c;
        RotationRed::new(ra,rb)
    }
}

impl<T> Mul<ReflectionRed<T>> for RotationRed<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ReflectionRed<T>;

    fn mul(self, other: ReflectionRed<T>) -> ReflectionRed<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let sa = a.clone()*c.clone() - b.clone()*d.clone();
        let sb = a*d - b*c;
        ReflectionRed::new(sa,sb)
    }
}

impl<T> Mul<RotationRed<T>> for ReflectionRed<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ReflectionRed<T>;

    fn mul(self, other: RotationRed<T>) -> ReflectionRed<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let sa = a.clone()*c.clone() + b.clone()*d.clone();
        let sb = a*d + b*c;
        ReflectionRed::new(sa,sb)
    }
}


/// green Rotation
#[derive(Debug,Clone)]
pub struct RotationGreen<T> {
    vector: TwoProportion<T>,
}

impl<T> RotationGreen<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> RotationGreen<T> {
        RotationGreen { vector: TwoProportion::new(a, b) }
    }
}

impl<T> Mul<RotationGreen<T>> for ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: RotationGreen<T>) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        let a = r.vector.a.clone();
        let b = r.vector.b.clone();
        let rx = a.clone()*x.clone() + b.clone()*y.clone();
        let ry = b*x + a*y;
        ProjOnePoint::new(rx,ry)
    }
}

impl<T> Mul<RotationGreen<T>> for RotationGreen<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = RotationGreen<T>;

    fn mul(self, other: Self) -> RotationGreen<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let ra = a.clone()*c.clone() + b.clone()*d.clone();
        let rb = a*d + b*c;
        RotationGreen::new(ra,rb)
    }
}

impl<T> One for RotationGreen<T>
where
    T: Sub<T, Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn one() -> RotationGreen<T> {
        RotationGreen::new(T::one(),T::zero())
    }
}

/// green Reflection
#[derive(Debug,Clone)]
pub struct ReflectionGreen<T> {
    vector: TwoProportion<T>,
}

impl<T> ReflectionGreen<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> ReflectionGreen<T> {
        ReflectionGreen { vector: TwoProportion::new(a, b) }
    }
}

impl<T> Mul<ReflectionGreen<T>> for ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: ReflectionGreen<T>) -> ProjOnePoint<T> {
        let x = self.x.a.clone();
        let y = self.x.b.clone();
        let a = r.vector.a.clone();
        let b = r.vector.b.clone();
        let rx = a.clone()*x.clone() - b.clone()*y.clone();
        let ry = b*x - a*y;
        ProjOnePoint::new(rx,ry)
    }
}

impl<T> Mul<ReflectionGreen<T>> for ReflectionGreen<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = RotationGreen<T>;

    fn mul(self, other: Self) -> RotationGreen<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let ra = a.clone()*c.clone() - b.clone()*d.clone();
        let rb = a*d - b*c;
        RotationGreen::new(ra,rb)
    }
}

impl<T> Mul<ReflectionGreen<T>> for RotationGreen<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ReflectionGreen<T>;

    fn mul(self, other: ReflectionGreen<T>) -> ReflectionGreen<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let sa = a.clone()*c.clone() - b.clone()*d.clone();
        let sb = a*d - b*c;
        ReflectionGreen::new(sa,sb)
    }
}

impl<T> Mul<RotationGreen<T>> for ReflectionGreen<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = ReflectionGreen<T>;

    fn mul(self, other: RotationGreen<T>) -> ReflectionGreen<T> {
        let a = self.vector.a.clone();
        let b = self.vector.b.clone();
        let c = other.vector.a.clone();
        let d = other.vector.b.clone();
        let sa = a.clone()*c.clone() + b.clone()*d.clone();
        let sb = a*d + b*c;
        ReflectionGreen::new(sa,sb)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::Ratio;
    #[test]
    fn one_dimensional_euclidean_projective_quadrance() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(1,2));
    }
    #[test]
    fn one_dimensional_euclidean_projective_quadrance_2() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(4,1));
        let q = Ratio::new(1,1)-Ratio::new(11*11,1)/Ratio::new(13*17,1);
        assert_eq!(a1.quadrance(&a2),q);
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_quadrance() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        assert_eq!(a1.quadrance_red(&a2),Ratio::new(25,24));
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_quadrance_2() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance_red(&a2),Ratio::new(-49,15));
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_quadrance_3() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance_red(&a2),Ratio::new(121,40));
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_quadrance() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        assert_eq!(a1.quadrance_green(&a2),Ratio::new(-25,24));
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_quadrance_2() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance_green(&a2),Ratio::new(49,48));
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_quadrance_3() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance_green(&a2),Ratio::new(121,72));
    }
    #[test]
    fn one_dimensional_euclidean_projective_perpendicularity() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(4,1));
        let a2 = ProjOnePoint::new(Ratio::new(-8,1),Ratio::new(6,1));
        assert!(a1.is_perpendicular(&a2));
    }
    #[test]
    fn one_dimensional_euclidean_projective_perpendicularity_2() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(0,1));
        let a2 = ProjOnePoint::new(Ratio::new(0,1),Ratio::new(1,1));
        assert!(a1.is_perpendicular(&a2));
    }
    #[test]
    fn one_dimensional_euclidean_projective_perpendicularity_3() {
        let a1 = ProjOnePoint::new(Ratio::new(2,3),Ratio::new(1,5));
        let a2 = ProjOnePoint::new(Ratio::new(-3,1),Ratio::new(10,1));
        assert!(a1.is_perpendicular(&a2));
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_perpendicularity() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(-1,1));
        assert!(a1.is_perpendicular_red(&a1));
        assert!(a2.is_perpendicular_red(&a2));
    }
    #[test]
    fn perpendicular_of_one_dimensional_projective_point_red() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(2,1),Ratio::new(1,1));
        assert_eq!(a1.perpendicular_red(),a2);
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_perpendicularity() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(-1,1));
        assert!(a1.is_perpendicular_green(&a2));
    }
    #[test]
    fn perpendicular_of_one_dimensional_projective_point_green() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(-2,1));
        assert_eq!(a1.perpendicular_green(),a2);
    }
    #[test]
    fn rotation_of_one_dimensional_projective_points() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(0,1));
        let a2 = ProjOnePoint::new(Ratio::new( 3,1),Ratio::new( 4,1));
        let a3 = ProjOnePoint::new(Ratio::new( 0,1),Ratio::new( 1,1));
        let a4 = ProjOnePoint::new(Ratio::new(-4,1),Ratio::new( 3,1));
        let r = Rotation::new(Ratio::new( 3,1),Ratio::new( 4,1));
        assert_eq!(a1.clone()*r.clone(),a2.clone());
        assert_eq!(a3.clone()*r,a4.clone());
        assert_eq!(a1.quadrance(&a3),a2.quadrance(&a4));
    }
    #[test]
    fn rotation_of_one_dimensional_projective_points_2() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a3 = ProjOnePoint::new(Ratio::new(-2,1),Ratio::new(1,1));
        let r = Rotation::new(Ratio::new(-1,1),Ratio::new(1,1));
        assert_eq!(a1.clone()*r.clone(),a2.clone());
        assert_eq!(a2*r,a3);
    }
    #[test]
    fn multiplication_of_rotations_of_one_dimensional_projective_points() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(5,1));
        let r1 = Rotation::new(Ratio::new( 7,1),Ratio::new( 4,1));
        let r2 = Rotation::new(Ratio::new( 11,1),Ratio::new( 4,1));
        assert_eq!((a1.clone()*r1.clone())*r2.clone(),a1*(r1*r2));
    }
    #[test]
    fn reflection_of_one_dimensional_projective_points() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(0,1));
        let a2 = ProjOnePoint::new(Ratio::new(-1,1),Ratio::new(1,1));
        let a3 = ProjOnePoint::new(Ratio::new(0,1),Ratio::new(1,1));
        let a4 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(1,1));
        let a5 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a6 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(3,1));
        let a7 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a8 = ProjOnePoint::new(Ratio::new(-1,1),Ratio::new(2,1));

        let s = Reflection::new(Ratio::new(-1,1),Ratio::new(1,1));

        assert_eq!(a1.clone()*s.clone(),a2.clone());
        assert_eq!(a3.clone()*s.clone(),a4.clone());
        assert_eq!(a2.clone()*s.clone(),a1.clone());
        assert_eq!(a4.clone()*s.clone(),a3.clone());
        assert_eq!(a5.clone()*s.clone(),a6.clone());
        assert_eq!(a7.clone()*s.clone(),a8.clone());
    }
    #[test]
    fn perpendicular_of_one_dimensional_projective_point() {
        let a1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(-2,1),Ratio::new(1,1));
        assert_eq!(a1.perpendicular(),a2);
    }
    #[test]
    fn projective_quadruple_quad_formula() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let a3 = ProjOnePoint::new(Ratio::new(-1,1),Ratio::new(1,1));
        let a4 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(2,1));

        let q1 = a1.quadrance(&a2);
        let q2 = a2.quadrance(&a3);
        let q3 = a3.quadrance(&a4);
        let q4 = a4.quadrance(&a1);

        assert_eq!(q1,Ratio::new(1,2));
        assert_eq!(q2,Ratio::new(9,10));
        assert_eq!(q3,Ratio::new(25,26));
        assert_eq!(q4,Ratio::new(9,130));

        let r1 = q1+q2+q3+q4;
        let r2 = r1*r1 - (q1*q1+q2*q2+q3*q3+q4*q4)*2
                 - (q1*q2*q3 + q1*q2*q4+q1*q3*q4+q2*q3*q4)*4
                 + q1*q2*q3*q4*8;
        let r = r2*r2 - q1*q2*q3*q4*64*(q1-1)*(q2-1)*(q3-1)*(q4-1);
        assert!(r.is_zero());
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_rotation() {
        let v1 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(3,4));
        let v2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(1,5));
        let v3 = RotationRed::new(Ratio::new(1,1),Ratio::new(1,2));
        assert_eq!(v1.quadrance_red(&v2),Ratio::new(-121,168));
        assert_eq!((v1*v3.clone()).quadrance_red(&(v2*v3)),Ratio::new(-121,168));
    }
}

