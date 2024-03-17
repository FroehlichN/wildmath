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
use std::ops::{Mul, Div, Add, Sub, Neg};
use linalg::{Matrix, ColumnVector};


/// Represents some geometric object or nothing
#[derive(Debug, Clone, PartialEq)]
pub enum GeoObj<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    Nothing,
    Point(Point<T>),
    Line(Line<T>),
}


/// Represents a 3D point
#[derive(Debug, Clone)]
pub struct Point<T> {
    pub x : T,
    pub y : T,
    pub z : T,
}


impl<T> Point<T>
{
    pub fn new(x: T, y: T, z: T) -> Point<T> {
        Point { x: x, y: y, z: z }
    }
}


impl<T> Point<T>
where
    T: Zero,
    T: Add,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn lies_on(&self, pi: &Plane<T>) -> bool {
        let r = pi.a.clone() * self.x.clone()
              + pi.b.clone() * self.y.clone()
              + pi.c.clone() * self.z.clone()
              + pi.d.clone();
        r.is_zero()
    }
}

impl<T> PartialEq for Point<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}


/// Represents a line in 3D
#[derive(Debug, Clone)]
pub struct Line<T> {
    a: Point<T>,
    b: Point<T>,
}

impl<T> Line<T>
{
    pub fn new(a: Point<T>, b: Point<T>) -> Line<T> {
        Line { a:a, b:b }
    }
}

impl<T> PartialEq for Line<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        other.passes_through(&self.a) && other.passes_through(&self.b)
    }
}

impl<T> Line<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn passes_through(&self, point: &Point<T>) -> bool {
        let mut mv = Vec::new();
        mv.push(vec![self.a.x.clone(), self.a.y.clone(), self.a.z.clone()]);
        mv.push(vec![self.b.x.clone(), self.b.y.clone(), self.b.z.clone()]);
        mv.push(vec![point.x.clone(), point.y.clone(), point.z.clone()]);
        let m = Matrix::new(mv);
        m.determinant().is_zero()
    }
}


/// Represents a plane in 3D
/// a*x + b*y + c*z + d = 0
#[derive(Debug, Clone)]
pub struct Plane<T> {
    a: T,
    b: T,
    c: T,
    d: T,
}


impl<T> Plane<T>
{
    pub fn new(a: T, b: T, c: T, d: T) -> Plane<T> {
        Plane { a: a, b: b, c: c, d: d }
    }
}

impl<T> Plane<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Mul<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn meet(self, others: &[Self]) -> GeoObj<T> {
        if others.len() == 0 {
            return GeoObj::<T>::Nothing;
        }

        if others.len() == 2 {
            let mut mv = Vec::new();
            let mut dv = Vec::new();
            mv.push(vec![self.a.clone(), self.b.clone(), self.c.clone()]);
            dv.push( - self.d.clone());
            for (_, pi) in others.iter().enumerate() {
                mv.push(vec![pi.a.clone(), pi.b.clone(), pi.c.clone()]);
                dv.push(T::zero() - pi.d.clone());
            }
            let m = Matrix::new(mv);
            let oinv = m.inverse();
            let inv = match oinv {
                None => return GeoObj::<T>::Nothing,
                Some(inv) => inv,
            };
            let cv = ColumnVector::new(dv);
            let rv = inv*cv;
            return GeoObj::<T>::Point(Point::new(rv.get(0),rv.get(1),rv.get(2)));
        }

        GeoObj::<T>::Nothing
    }
}


/// Represents a 3D tetrahedon
/// T is type of the coordinates of the points
pub struct Tetrahedron<T> {
    points: [Point<T>; 3],
}


impl<T> Tetrahedron<T>
where
    T: One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new(a1: Point<T>, a2: Point<T>, a3: Point<T>)
        -> Tetrahedron<T> {
/*        if a1.is_collinear(&a2, &a3) {
            panic!("Tetrahedron cannot have collinear points");
        }*/
        Tetrahedron { points: [a1, a2, a3] }
    }

    pub fn volume(&self) -> T {
        let six = T::one() + T::one() + T::one()
                + T::one() + T::one() + T::one();
        let x1 = self.points[0].x.clone();
        let x2 = self.points[1].x.clone();
        let x3 = self.points[2].x.clone();
        let y1 = self.points[0].y.clone();
        let y2 = self.points[1].y.clone();
        let y3 = self.points[2].y.clone();
        let z1 = self.points[0].z.clone();
        let z2 = self.points[1].z.clone();
        let z3 = self.points[2].z.clone();

        let n = x1.clone()*y2.clone()*z3.clone()
              - x1*y3.clone()*z2.clone()
              + x2.clone()*y3*z1.clone()
              - x3.clone()*y2*z1
              + x3*y1.clone()*z2
              - x2*y1*z3;
        return n / six;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::Ratio;

    #[test]
    fn points_lie_on_a_plane() {
        let p1 = Point::new(3,0,0);
        let p2 = Point::new(0,1,0);
        let p3 = Point::new(0,0,2);
        let pi = Plane::new(2,6,3,-6);
        assert!(p1.lies_on(&pi));
        assert!(p2.lies_on(&pi));
        assert!(p3.lies_on(&pi));
    }

    #[test]
    fn three_plains_meet_in_one_point() {
        let pi1 = Plane::new(Ratio::from(1),Ratio::from(2),Ratio::from(-1),Ratio::from(3));
        let pi2 = Plane::new(Ratio::from(3),Ratio::from(7),Ratio::from(2),Ratio::from(-1));
        let pi3 = Plane::new(Ratio::from(4),Ratio::from(-2),Ratio::from(1),Ratio::from(2));

        let p1 = pi1.clone().meet(&[pi2.clone(),pi3.clone()]);
        let p2 = GeoObj::Point(Point::new(Ratio::from(-1),Ratio::from(0),Ratio::from(2)));
        assert_eq!(p1,p2);
        _ = match p1 {
            GeoObj::Point(p) => {assert!(p.lies_on(&pi1));
                                 assert!(p.lies_on(&pi2));
                                 assert!(p.lies_on(&pi3));},
            _ => assert!(false),
        }
    }

    #[test]
    fn line_passes_through_point() {
        let p1 = Point::new(3,1,1);
        let p2 = Point::new(2,2,-1);
        let p3 = Point::new(1,3,-3);
        let p4 = Point::new(0,4,3);
        let l = Line::new(p1,p2);
        assert!(l.passes_through(&p3));
        assert!(!l.passes_through(&p4));
    }
}

