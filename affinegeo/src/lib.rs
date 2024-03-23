/*
Copyright 2024 Norbert Fröhlich


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

use num::{Zero, One};
use std::ops::{Mul, Add, Sub, Div};
use linalg::{ColumnVector, RowVector, Matrix};



/// Represents an n-D point
#[derive(Debug, Clone)]
pub struct Point<T> {
    coords: Vec<T>,
}

impl<T> Point<T>
{
    pub fn new(coords: Vec<T>) -> Point<T> {
        Point { coords: coords }
    }
}

impl<T> Point<T>
where
    T: Zero,
    T: PartialEq,
    T: Add,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn lies_on(&self, pi: &Plane<T>) -> bool {
        let pointv = RowVector::new(self.coords.clone());
        let planev = ColumnVector::new(pi.coords.clone());
        let r = pointv * planev;
        r == pi.d.clone()
    }
}

impl<T> Point<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn quadrance(&self, other: &Self) -> T {
        let v = Vector::new((*self).clone(), (*other).clone());
        v.quadrance_blue()
    }
}

impl<T> Add<Vector<T>> for Point<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    type Output = Point<T>;

    fn add(self, other: Vector<T>) -> Point<T> {
        let selfvec = ColumnVector::new(self.coords);
        let sumvec = selfvec + other.columnvector();
        Point::new(sumvec.elem)
    }
}

impl<T> PartialEq for Point<T>
where
    T: Zero,
    T: PartialEq,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let svec = ColumnVector::new(self.coords.clone());
        let ovec = ColumnVector::new(other.coords.clone());
        svec == ovec
    }
}


/// Represents an n-D vector
#[derive(Debug, Clone)]
pub struct Vector<T> {
    start: Point<T>,
    end: Point<T>,
}

impl<T> Vector<T>
{
    pub fn new(start: Point<T>, end: Point<T>) -> Vector<T> {
        Vector {start: start, end: end}
    }
}

impl<T> Vector<T>
where
    T: Zero,
{
    pub fn from(coords: Vec<T>) -> Vector<T> {
        let mut zero = Vec::new();
        for (_, _) in coords.iter().enumerate() {
            zero.push(T::zero());
        }
        let start = Point::new(zero);
        let end = Point::new(coords);
        Vector::new(start, end)
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Clone,
{
    pub fn columnvector(&self) -> ColumnVector<T> {
        let selfs = ColumnVector::new(self.start.coords.clone());
        let selfe = ColumnVector::new(self.end.coords.clone());
        selfe - selfs
    }
    pub fn rowvector(&self) -> RowVector<T> {
        let selfs = RowVector::new(self.start.coords.clone());
        let selfe = RowVector::new(self.end.coords.clone());
        selfe - selfs
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn dot_blue(&self, other: &Self) -> T {
        let v1 = self.rowvector();
        let v2 = other.columnvector();

        v1*v2
    }

    pub fn quadrance_blue(&self) -> T {
        self.dot_blue(&self)
    }
}

impl<T> PartialEq for Vector<T>
where
    T: PartialEq,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let s = self.columnvector();
        let o = other.columnvector();
        s == o
    }
}

impl<T> Add for Vector<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn add(self, other: Self) -> Vector<T> {
        let selfendvec = ColumnVector::new(self.end.coords);
        let sumendvec = selfendvec + other.columnvector();
        let sumend = Point::new(sumendvec.elem);
        Vector::new(self.start, sumend)
    }
}

impl<T> Mul<T> for Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn mul(self, scalar: T) -> Vector<T> {
        let vec = self.columnvector();
        let startvec = ColumnVector::new(self.start.coords.clone());
        let endvec = startvec + vec * scalar;
        let end = Point::new(endvec.elem);
        Vector::new(self.start, end)
    }
}


/// Represents an n-D plane
/// a1*x1 + a2*x2 + a3*x3 + ... = d
#[derive(Debug, Clone)]
pub struct Plane<T> {
    coords: Vec<T>,
    d: T,
}

impl<T> Plane<T>
{
    pub fn new(coords: Vec<T>, d: T) -> Plane<T> {
        Plane { coords: coords, d: d }
    }
}

impl<T> Plane<T>
where
    T: Zero + One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn projection_matrix(&self, direction: &Vector<T>) -> Matrix<T> {
        let pv = self.coords.clone();
        let dv = direction.columnvector().elem;

        let s = T::one() / (RowVector::<T>::new(pv.clone()) * ColumnVector::<T>::new(dv.clone()));

        let ident = Matrix::<T>::identitiy(dv.len(),pv.len());

        let mut dmtv = Vec::new();
        dmtv.push(dv);
        let dmt = Matrix::<T>::new(dmtv);
        let dm = dmt.transpose();

        let mut pmv = Vec::new();
        pmv.push(pv);
        let pm = Matrix::<T>::new(pmv);

        ident - dm*pm*s
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};

    #[test]
    fn points_lie_on_a_plane() {
        let p1 = Point::new(vec![3,0,0]);
        let p2 = Point::new(vec![0,1,0]);
        let p3 = Point::new(vec![0,0,2]);
        let pi = Plane::new(vec![2,6,3],6);
        assert!(p1.lies_on(&pi));
        assert!(p2.lies_on(&pi));
        assert!(p3.lies_on(&pi));
    }

    #[test]
    fn equal_vectors() {
        let a1 = Point::new(vec![Ratio::new(1,1), Ratio::new(1,1)]);
        let a2 = Point::new(vec![Ratio::new(3,1), Ratio::new(4,1)]);
        let a3 = Point::new(vec![Ratio::new(4,1), Ratio::new(2,1)]);
        let a4 = Point::new(vec![Ratio::new(6,1), Ratio::new(5,1)]);
        let v1 = Vector::new(a1, a2);
        let v2 = Vector::new(a3, a4);
        assert_eq!(v1,v2);
    }

    #[test]
    fn add_vectors() {
        let v1 = Vector::from(vec![Ratio::new(2,1), Ratio::new(3,1)]);
        let v2 = Vector::from(vec![Ratio::new(2,1), Ratio::new(-4,1)]);
        let v3 = Vector::from(vec![Ratio::new(4,1), Ratio::new(-1,1)]);
        assert_eq!(v1+v2,v3);
    }

    #[test]
    fn add_point_and_vector() {
        let p1 = Point::new(vec![1,2,3]);
        let v1 = Vector::from(vec![4,5,6]);
        let p2 = Point::new(vec![5,7,9]);
        assert_eq!(p1+v1,p2);
    }

    #[test]
    fn scalar_vector_mul() {
        let s = Ratio::new(3,1);
        let v1 = Vector::from(vec![Ratio::new(2,1), Ratio::new(-4,1)]);
        let v2 = Vector::from(vec![Ratio::new(6,1), Ratio::new(-12,1)]);
        assert_eq!(v1*s,v2);
    }

    #[test]
    fn projections_onto_a_plane() {
        let pi1 = Plane::new(vec![Ratio::from(1),Ratio::from(-5),Ratio::from(2)],Ratio::from(0));
        let v1 = Vector::from(vec![Ratio::from(1),Ratio::from(1),Ratio::from(3)]);
        let v2 = Vector::from(vec![Ratio::from(1),Ratio::from(-5),Ratio::from(2)]);
        let pm1 = pi1.projection_matrix(&v1);
        let pm2 = pi1.projection_matrix(&v2);
        assert_eq!(pm1.clone()*pm1.clone(),pm1.clone());
        assert_eq!(pm2.clone()*pm2.clone(),pm2.clone());

        let u1 = ColumnVector::new(vec![Ratio::from(12),Ratio::from(-341),Ratio::from(35)]);
        let u2 = ColumnVector::new(vec![Ratio::from(-123),Ratio::from(-55),Ratio::from(32)]);

        let pu11 = pm1.clone()*u1.clone();
        let pu12 = pm1*u2.clone();
        let pu21 = pm2.clone()*u1;
        let pu22 = pm2*u2;

        let vpu11 = Vector::from(pu11.elem);
        let vpu12 = Vector::from(pu12.elem);
        let vpu21 = Vector::from(pu21.elem);
        let vpu22 = Vector::from(pu22.elem);

        let o = Point::new(vec![Ratio::from(0),Ratio::from(0),Ratio::from(0)]);
        assert!((o.clone() + vpu11).lies_on(&pi1));
        assert!((o.clone() + vpu12).lies_on(&pi1));
        assert!((o.clone() + vpu21).lies_on(&pi1));
        assert!((o + vpu22).lies_on(&pi1));
    }

    #[test]
    fn quadrance() {
        let p1 = Point::new(vec![Ratio::new(2,1), Ratio::new(1,1)]);
        let p2 = Point::new(vec![Ratio::new(6,1), Ratio::new(2,1)]);
        let q  = Ratio::new(17,1);
        assert_eq!(p1.quadrance(&p2),q);
    }

    #[test]
    fn quadrance_between_one_points() {
        let a1 = Point::new(vec![2]);
        let a2 = Point::new(vec![5]);
        let a3 = Point::new(vec![7]);
        assert_eq!(a2.quadrance(&a3),4);
        assert_eq!(a1.quadrance(&a3),25);
        assert_eq!(a1.quadrance(&a2),9);
    }

}

