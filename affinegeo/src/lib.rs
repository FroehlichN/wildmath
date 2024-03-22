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

use num::{Zero};
use std::ops::{Mul, Add, Sub};
use linalg::{ColumnVector, RowVector};



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

}

