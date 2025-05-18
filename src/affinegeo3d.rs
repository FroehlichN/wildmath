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
use crate::matrix::{Matrix, ColumnVector, RowVector};
use crate::create_polynumber_one;
use crate::create_polynumber_var;


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

impl<T> Point<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn is_collinear(&self, p2: &Self, p3: &Self) -> bool {
        let l = Line::new(p2.clone(),p3.clone());
        l.passes_through(&self)
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

impl<T> Add<Vector<T>> for Point<T>
where
    T: Add<Output = T>,
{
    type Output = Point<T>;

    fn add(self, v: Vector<T>) -> Point<T> {
        let x = self.x + v.x;
        let y = self.y + v.y;
        let z = self.z + v.z;
        Point { x: x, y: y, z: z }
    }
}


/// Represents a 3D vector
#[derive(Debug, Clone)]
pub struct Vector<T> {
    pub x : T,
    pub y : T,
    pub z : T,
}


impl<T> Vector<T>
{
    pub fn new(x: T, y: T, z: T) -> Vector<T> {
        Vector { x: x, y: y, z: z }
    }
}


impl<T> Vector<T>
where
    T: Sub<Output = T>,
{
    pub fn newp(start: Point<T>, end: Point<T>) -> Vector<T> {
        let x = end.x-start.x;
        let y = end.y-start.y;
        let z = end.z-start.z;

        Vector { x: x, y: y, z: z }
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: PartialEq,
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    pub fn is_collinear(&self, other: &Self) -> bool {
        self.clone().cross(other.clone()).is_zero()
    }
}

impl<T> Vector<T>
where
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    pub fn cross(self, other: Self) -> Vector<T> {
        let x = self.y.clone() * other.z.clone() - self.z.clone() * other.y.clone();
        let y = self.z.clone() * other.x.clone() - self.x.clone() * other.z.clone();
        let z = self.x.clone() * other.y.clone() - self.y.clone() * other.x.clone();
        Vector { x:x, y:y, z:z }
    }
}

impl<T> Zero for Vector<T>
where
    T: Zero,
    T: PartialEq,
{
    fn zero() -> Vector<T> {
        Vector { x: T::zero(), y: T::zero(), z: T::zero() }
    }

    fn is_zero(&self) -> bool {
        self.x == T::zero() && self.y == T::zero() && self.z == T::zero()
    }
}

impl<T> Add for Vector<T>
where
    T: Add<Output = T>,
{
    type Output = Vector<T>;

    fn add(self, other: Self) -> Vector<T> {
        let x = self.x + other.x;
        let y = self.y + other.y;
        let z = self.z + other.z;
        Vector { x:x, y:y, z:z }
    }
}

impl<T> Mul<T> for Vector<T>
where
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn mul(self, scalar: T) -> Vector<T> {
        let x = scalar.clone() * self.x;
        let y = scalar.clone() * self.y;
        let z = scalar*self.z;

        Vector::new(x,y,z)
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

impl<T> Line<T>
where
    T: Add<Output = T>,
    T: Clone,
{
    pub fn newv(p: Point<T>, v: Vector<T>) -> Line<T> {
        Line { a:p.clone(), b:p+v }
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

impl<T> Line<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn meet(&self, plane: &Plane<T>) -> Option<Point<T>> {
        let a = plane.a.clone();
        let b = plane.b.clone();
        let c = plane.c.clone();
        let d = plane.d.clone();
        let x1 = self.a.x.clone();
        let x2 = self.b.x.clone();
        let y1 = self.a.y.clone();
        let y2 = self.b.y.clone();
        let z1 = self.a.z.clone();
        let z2 = self.b.z.clone();
        let denom = a.clone()*(x1.clone()-x2.clone())
                   +b.clone()*(y1.clone()-y2.clone())
                   +c.clone()*(z1.clone()-z2.clone());
        if denom.is_zero() {
            return None;
        }
        let numer = a*x1 + b*y1 + c*z1 + d;
        let lambda = numer / denom;
        Some(self.a.clone() + Vector::newp(self.a.clone(),self.b.clone())*lambda)
    }
}

impl<T> Line<T>
where
    T: Zero,
    T: Add,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn lies_on(&self, plane: &Plane<T>) -> bool {
        self.a.lies_on(&plane) && self.b.lies_on(&plane)
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
    T: PartialEq,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn newp(a: Point<T>, b: Point<T>, c: Point<T>) -> Plane<T> {
        if a.is_collinear(&b, &c) {
            panic!("Plane cannot be defined by collinear points");
        }
        let x = create_polynumber_var!(x; x,y,z; T);
        let y = create_polynumber_var!(y; x,y,z; T);
        let z = create_polynumber_var!(z; x,y,z; T);
        let one = create_polynumber_one!(x,y,z; T);

        let mut mv = Vec::new();
        mv.push(vec![x - one.clone()*a.x.clone(),
                     y - one.clone()*a.y.clone(),
                     z - one.clone()*a.z.clone()]);
        mv.push(vec![one.clone()*(b.x - a.x.clone()),
                     one.clone()*(b.y - a.y.clone()),
                     one.clone()*(b.z - a.z.clone())]);
        mv.push(vec![one.clone()*(c.x - a.x),
                     one.clone()*(c.y - a.y),
                     one*(c.z - a.z)]);
        let m = Matrix::new(mv);
        let det = m.determinant();
        let a = det.get(1).get(0).get(0);
        let b = det.get(0).get(1).get(0);
        let c = det.get(0).get(0).get(1);
        let d = det.get(0).get(0).get(0);
        Plane { a: a, b: b, c: c, d: d }
    }

    pub fn newpv(a: Point<T>, u: Vector<T>, v: Vector<T>) -> Plane<T> {
        if u.is_collinear(&v) {
            panic!("Plane cannot be defined by collinear vectors");
        }
        Plane::newp(a.clone(), a.clone()+u, a+v)
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
    pub fn meet(&self, others: &[Self]) -> GeoObj<T> {
        if others.len() == 0 {
            return GeoObj::<T>::Nothing;
        }

        if others.len() == 1 {
            let pi1 = Plane::new(T::one(),T::zero(),T::zero(),T::zero());
            let pi2 = Plane::new(T::zero(),T::one(),T::zero(),T::zero());
            let pi3 = Plane::new(T::zero(),T::zero(),T::one(),T::zero());

            let go1 = self.meet(&[others[0].clone(), pi1]);
            let go2 = self.meet(&[others[0].clone(), pi2]);
            let go3 = self.meet(&[others[0].clone(), pi3]);

            match (go1, go2, go3) {
                (GeoObj::<T>::Nothing, GeoObj::<T>::Nothing, GeoObj::<T>::Nothing) =>
                    return GeoObj::<T>::Nothing,
                (GeoObj::<T>::Nothing, GeoObj::<T>::Nothing, GeoObj::<T>::Point(_p3)) =>
                    return GeoObj::<T>::Nothing,
                (GeoObj::<T>::Nothing, GeoObj::<T>::Point(_p2), GeoObj::<T>::Nothing) =>
                    return GeoObj::<T>::Nothing,
                (GeoObj::<T>::Nothing, GeoObj::<T>::Point(p2), GeoObj::<T>::Point(p3)) =>
                    return GeoObj::<T>::Line(Line::new(p2,p3)),
                (GeoObj::<T>::Point(_p1), GeoObj::<T>::Nothing, GeoObj::<T>::Nothing) =>
                    return GeoObj::<T>::Nothing,
                (GeoObj::<T>::Point(p1), GeoObj::<T>::Nothing, GeoObj::<T>::Point(p3)) =>
                    return GeoObj::<T>::Line(Line::new(p1,p3)),
                (GeoObj::<T>::Point(p1), GeoObj::<T>::Point(p2), GeoObj::<T>::Nothing) =>
                    return GeoObj::<T>::Line(Line::new(p1,p2)),
                (GeoObj::<T>::Point(p1), GeoObj::<T>::Point(p2), GeoObj::<T>::Point(_p3)) =>
                    return GeoObj::<T>::Line(Line::new(p1,p2)),
                _ => return GeoObj::<T>::Nothing,
            };
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
        let pv = vec![self.a.clone(), self.b.clone(), self.c.clone()];
        let dv = vec![direction.x.clone(), direction.y.clone(), direction.z.clone()];
        let zerov = vec![T::zero(), T::zero(), T::zero()];

        let s = T::one() / (RowVector::<T>::new(pv.clone()) * ColumnVector::<T>::new(dv.clone()));

        let ident = Matrix::<T>::identity(3,3);

        let mut dmtv = Vec::new();
        dmtv.push(dv);
        dmtv.push(zerov.clone());
        dmtv.push(zerov.clone());
        let dmt = Matrix::<T>::new(dmtv);
        let dm = dmt.transpose();

        let mut pmv = Vec::new();
        pmv.push(pv);
        pmv.push(zerov.clone());
        pmv.push(zerov.clone());
        let pm = Matrix::<T>::new(pmv);

        ident - dm*pm*s
    }
}


/// Represents a 3D tetrahedon
/// T is type of the coordinates of the points
pub struct Tetrahedron<T> {
    points: [Point<T>; 3],
}


impl<T> Tetrahedron<T>
where
    T: Zero + One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn new(a1: Point<T>, a2: Point<T>, a3: Point<T>)
        -> Tetrahedron<T> {
        if a1.is_collinear(&a2, &a3) {
            panic!("Tetrahedron cannot have collinear points");
        }
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

    #[test]
    fn meet_of_line_and_plane() {
        let p = Point::new(Ratio::from(1),Ratio::from(2),Ratio::from(-3));
        let v = Vector::new(Ratio::from(-1),Ratio::from(4),Ratio::from(1));
        let l = Line::newv(p,v);
        let pi = Plane::new(Ratio::from(1),Ratio::from(3),Ratio::from(-1),Ratio::from(-5));
        let pm = l.meet(&pi);
        assert!(l.passes_through(&pm.as_ref().unwrap()));
        assert!(pm.as_ref().unwrap().lies_on(&pi));
        let pr = Point::new(Ratio::new(3,2),Ratio::from(0),Ratio::new(-7,2));
        assert_eq!(*pm.as_ref().unwrap(),pr);
    }

    #[test]
    fn two_planes_meet_in_one_line() {
        let pi1 = Plane::new(Ratio::from(1),Ratio::from(2),Ratio::from(-1),Ratio::from(3));
        let pi2 = Plane::new(Ratio::from(3),Ratio::from(7),Ratio::from(2),Ratio::from(-1));
        let go = pi1.meet(&[pi2.clone()]);
        match go {
            GeoObj::Line(l) => {
                    assert!(l.lies_on(&pi1));
                    assert!(l.lies_on(&pi2));
                    let p = Point::new(Ratio::from(-23),Ratio::from(10),Ratio::from(0));
                    let v = Vector::new(Ratio::from(11),Ratio::from(-5),Ratio::from(1));
                    let l2 = Line::newv(p,v);
                    assert_eq!(l,l2);
                },
            _ => assert!(false),
        }
    }

    #[test]
    fn two_planes_meet_in_one_line_2() {
        let pi1 = Plane::new(Ratio::from(1),Ratio::from(3),Ratio::from(-2),Ratio::from(-2));
        let pi2 = Plane::new(Ratio::from(2),Ratio::from(6),Ratio::from(-5),Ratio::from(-3));
        let go = pi1.meet(&[pi2.clone()]);
        match go {
            GeoObj::Line(l) => {
                    assert!(l.lies_on(&pi1));
                    assert!(l.lies_on(&pi2));
                    let p = Point::new(Ratio::from(4),Ratio::from(0),Ratio::from(1));
                    let v = Vector::new(Ratio::from(-3),Ratio::from(1),Ratio::from(0));
                    let l2 = Line::newv(p,v);
                    assert_eq!(l,l2);
                },
            _ => assert!(false),
        }
    }

    #[test]
    fn three_points_define_a_plane() {
        let p1 = Point::new(1,0,-1);
        let p2 = Point::new(2,1,3);
        let p3 = Point::new(-4,2,5);
        let pi = Plane::newp(p1.clone(),p2.clone(),p3.clone());
        assert!(p1.lies_on(&pi));
        assert!(p2.lies_on(&pi));
        assert!(p3.lies_on(&pi));
    }

    #[test]
    fn one_point_and_two_vectors_define_a_plane() {
        let a = Point::new(1,0,-1);
        let u = Vector::new(1,1,4);
        let v = Vector::new(-5,2,6);
        let pi = Plane::newpv(a.clone(),u.clone(),v.clone());
        assert!(a.lies_on(&pi));
        assert!((a.clone()+u).lies_on(&pi));
        assert!((a+v).lies_on(&pi));
    }

    #[test]
    fn projections_onto_a_plane() {
        let pi1 = Plane::new(Ratio::from(1),Ratio::from(-5),Ratio::from(2),Ratio::from(0));
        let v1 = Vector::new(Ratio::from(1),Ratio::from(1),Ratio::from(3));
        let v2 = Vector::new(Ratio::from(1),Ratio::from(-5),Ratio::from(2));
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

        let vpu11 = Vector::new(pu11.get(0), pu11.get(1), pu11.get(2));
        let vpu12 = Vector::new(pu12.get(0), pu12.get(1), pu12.get(2));
        let vpu21 = Vector::new(pu21.get(0), pu21.get(1), pu21.get(2));
        let vpu22 = Vector::new(pu22.get(0), pu22.get(1), pu22.get(2));

        let o = Point::new(Ratio::from(0),Ratio::from(0),Ratio::from(0));
        assert!((o.clone() + vpu11).lies_on(&pi1));
        assert!((o.clone() + vpu12).lies_on(&pi1));
        assert!((o.clone() + vpu21).lies_on(&pi1));
        assert!((o + vpu22).lies_on(&pi1));
    }
}

