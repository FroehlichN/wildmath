/*
Copyright 2023 - 2025 Norbert Fröhlich


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


use crate::proportion::Proportion;
use num::{Zero,One};
use std::ops::{Mul, Add, Sub, Div, Neg};
use crate::complex::Complex;
use crate::matrix::Matrix;



/// Represents a projective 1D point using a complex number
#[derive(Debug, Clone)]
pub struct ProjOnePoint<T> {
    c : Complex<T>,
}


impl<T> ProjOnePoint<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new(x: T, y: T) -> ProjOnePoint<T> {
        let _p = Proportion::new(vec![x.clone(),y.clone()]);
        let c = Complex::new(x,y);
        ProjOnePoint { c: c }
    }

    pub fn new_blue(x: T, y: T) -> ProjOnePoint<T> {
        let _p = Proportion::new(vec![x.clone(),y.clone()]);
        let c = Complex::new_blue(x,y);
        ProjOnePoint { c: c }
    }

    pub fn new_red(x: T, y: T) -> ProjOnePoint<T> {
        let _p = Proportion::new(vec![x.clone(),y.clone()]);
        let c = Complex::new_red(x,y);
        ProjOnePoint { c: c }
    }

    pub fn new_green(x: T, y: T) -> ProjOnePoint<T> {
        let two = T::one() + T::one();
        let a = (x.clone()+y.clone())/two.clone();
        let b = (x.clone()-y.clone())/two;
        let _p = Proportion::new(vec![x,y]);
        let c = Complex::new_green(a,b);
        ProjOnePoint { c: c }
    }
}


impl<T> ProjOnePoint<T>
where
    T: PartialEq,
    T: Neg<Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero + One,
    T: Clone,
{
    pub fn quadrance(&self, other: &Self) -> T {
        let n = self.c.clone().quadrance(other.c.clone());
        let numer = n.clone() * n;
        let d1 = self.c.clone().quadrance(self.c.clone());
        let d2 = other.c.clone().quadrance(other.c.clone());
        return T::one() - numer / (d1*d2)
    }
}

impl<T> ProjOnePoint<T>
where
    T: PartialEq,
    T: Neg<Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero + One,
    T: Clone,
{
    pub fn is_perpendicular(&self, other: &Self) -> bool {
        let q = self.c.clone().quadrance(other.c.clone());
        q.is_zero()
    }
    pub fn perpendicular(&self) -> ProjOnePoint<T> {
        let s = self.c.clone();
        let cc = s.clone().complex_conjugate();
        let p = s.clone()*(s.clone() - cc);
        ProjOnePoint { c: p }
    }
}

impl<T> PartialEq for ProjOnePoint<T>
where
    T: Zero + One,
    T: PartialEq,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let sp = Proportion::new(
            vec![self.c.matrix.get(0,0),self.c.matrix.get(0,1),
                 self.c.matrix.get(1,0),self.c.matrix.get(1,1)]);

        let op = Proportion::new(
            vec![other.c.matrix.get(0,0),other.c.matrix.get(0,1),
                 other.c.matrix.get(1,0),other.c.matrix.get(1,1)]);

        sp == op
    }
}


/// Rotation
#[derive(Debug,Clone)]
pub struct Rotation<T> {
    c: Complex<T>,
}

impl<T> Rotation<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new(a: T, b: T) -> Rotation<T> {
        let _p = Proportion::new(vec![a.clone(), b.clone()]);
        let c = Complex::new(a,b);
        Rotation { c: c }
    }

    pub fn new_red(a: T, b: T) -> Rotation<T> {
        let _p = Proportion::new(vec![a.clone(), b.clone()]);
        let c = Complex::new_red(a,b);
        Rotation { c: c }
    }

    pub fn new_blue(a: T, b: T) -> Rotation<T> {
        let _p = Proportion::new(vec![a.clone(), b.clone()]);
        let c = Complex::new_blue(a,b);
        Rotation { c: c }
    }

    pub fn new_green(a: T, b: T) -> Rotation<T> {
        let two = T::one() + T::one();
        let x = (a.clone()+b.clone())/two.clone();
        let y = (a.clone()-b.clone())/two;
        let _p = Proportion::new(vec![a.clone(), b.clone()]);
        let c = Complex::new_green(x,y);
        Rotation { c: c }
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
        let pr = self.c * r.c;
        ProjOnePoint { c: pr }
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
        let r = self.c * other.c;
        Rotation { c: r }
    }
}

impl<T> One for Rotation<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    fn one() -> Rotation<T> {
        Rotation::new(T::one(),T::zero())
    }
}

/// Reflection
#[derive(Debug,Clone)]
pub struct Reflection<T> {
    c: Complex<T>,
}

impl<T> Reflection<T>
where
    T: PartialEq,
    T: Zero + One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new(a: T, b: T) -> Reflection<T> {
        let _p = Proportion::new(vec![a.clone(), b.clone()]);
        let c = Complex::new(a,-b);
        Reflection { c: c }
    }
}

impl<T> Mul<Reflection<T>> for ProjOnePoint<T>
where
    T: PartialEq,
    T: Zero + One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    type Output = ProjOnePoint<T>;

    fn mul(self, r: Reflection<T>) -> ProjOnePoint<T> {
        let r = self.c * r.c;
        ProjOnePoint { c: r.complex_conjugate() }
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
        let r = self.c * other.c;
        Rotation { c: r }
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
        let r = self.c * other.c;
        Reflection { c: r }
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
        let r = self.c * other.c;
        Reflection { c: r }
    }
}

/// Represents a projective vector consisting of a
/// scalar u and a
/// vector v
#[derive(Debug,Clone)]
pub struct ProjVector<T> {
    u: T,
    v: Vec<T>,
}

impl<T> ProjVector<T>
where
    T: Zero,
{
    pub fn new(u: T, v: Vec<T>) -> ProjVector<T> {
        if u.is_zero() {
            let mut v_is_zero = true;
            for (_, e) in v.iter().enumerate() {
                if !e.is_zero() {
                    v_is_zero = false;
                    break;
                }
            }
            if v_is_zero {
                panic!("Projective vector must not be all zero");
            }
        }
        ProjVector{ u: u, v: v }
    }
}

impl<T> PartialEq for ProjVector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let mut sv = Vec::new();
        sv.push(self.u.clone());
        for (_, e) in self.v.iter().enumerate() {
            sv.push(e.clone());
        }

        let mut ov = Vec::new();
        ov.push(other.u.clone());
        for (_, e) in other.v.iter().enumerate() {
            ov.push(e.clone());
        }

        let m = Matrix::new(vec![sv,ov]);
        m.rank() == 1
    }
}

/// Represents reflections and rotations
#[derive(Debug, Clone)]
pub struct Isometry<T> {
    matrix: Matrix<T>,
}

impl<T> Isometry<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Div<Output = T>,
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn anchor(&self) -> Option<Matrix<T>> {
        let ao = self.matrix.cayley();
        match ao {
            None => return None,
            Some(a) => return Some(-a),
        }
    }

    pub fn from_anchor(u: T, anchor: Matrix<T>) -> Option<Isometry<T>> {
        if anchor.rows != anchor.cols {
            return None;
        }

        let q = anchor.dot(&anchor);
        let denom = u.clone()*u.clone() - q;

        if denom.is_zero() {
            return None;
        }

        let i = Matrix::identity(anchor.rows,anchor.cols);
        let two = T::one() + T::one();
        let r = i + (anchor.clone()*u + anchor.clone()*anchor)*two/denom;
        Some(Isometry{ matrix: r })
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
        let a1 = ProjOnePoint::new_red(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(2,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(25,24));
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_quadrance_2() {
        let a1 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new_red(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(-49,15));
    }
    #[test]
    fn one_dimensional_red_relativistic_projective_quadrance_3() {
        let a1 = ProjOnePoint::new_red(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new_red(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(121,40));
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_quadrance() {
        let a1 = ProjOnePoint::new_green(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new_green(Ratio::new(1,1),Ratio::new(2,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(-25,24));
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_quadrance_2() {
        let a1 = ProjOnePoint::new_green(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new_green(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(49,48));
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_quadrance_3() {
        let a1 = ProjOnePoint::new_green(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new_green(Ratio::new(2,1),Ratio::new(-3,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(121,72));
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
        let a1 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(-1,1));
        assert!(a1.is_perpendicular(&a1));
        assert!(a2.is_perpendicular(&a2));
    }
    #[test]
    fn perpendicular_of_one_dimensional_projective_point_red() {
        let a1 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new_red(Ratio::new(2,1),Ratio::new(1,1));
        assert_eq!(a1.perpendicular(),a2);
    }
    #[test]
    fn one_dimensional_green_relativistic_projective_perpendicularity() {
        let a1 = ProjOnePoint::new_green(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new_green(Ratio::new(3,1),Ratio::new(-1,1));
        assert!(a1.is_perpendicular(&a2));
    }
    #[test]
    fn perpendicular_of_one_dimensional_projective_point_green() {
        let a1 = ProjOnePoint::new_green(Ratio::new(1,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new_green(Ratio::new(1,1),Ratio::new(-2,1));
        assert_eq!(a1.perpendicular(),a2);
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
        let v1 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(3,4));
        let v2 = ProjOnePoint::new_red(Ratio::new(1,1),Ratio::new(1,5));
        let v3 = Rotation::new_red(Ratio::new(1,1),Ratio::new(1,2));
        assert_eq!(v1.quadrance(&v2),Ratio::new(-121,168));
        assert_eq!((v1*v3.clone()).quadrance(&(v2*v3)),Ratio::new(-121,168));
    }
    #[should_panic]
    fn invalid_projective_vector() {
        ProjVector::new(0,vec![0,0]);
    }
    #[test]
    fn equality_of_projective_vectors() {
        let v1 = ProjVector::new(1,vec![2,3,4]);
        let v2 = ProjVector::new(2,vec![4,6,8]);
        let v3 = ProjVector::new(3,vec![3,1,2]);
        assert_eq!(v1,v2);
        assert!(!(v1==v3));
    }
}

