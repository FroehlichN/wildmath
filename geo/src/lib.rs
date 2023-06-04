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


use num::{Num, Zero};
use std::ops::{Mul, Add, Sub};

pub trait Point {
    fn is_collinear(&self, a2: &Self, a3: &Self) -> bool;
}

/// Represents a 2D point
#[derive(Debug, Clone)]
pub struct TwoPoint<T> {
    x : T,
    y : T,
}

impl<T> TwoPoint<T>
where
    T: Num,
    T: Clone,
{
    pub fn lies_on(&self, l: &TwoLine<T>) -> bool {
        let r = l.a.clone() * self.x.clone()
              + l.b.clone() * self.y.clone() + l.c.clone();
        r.is_zero()
    }

    pub fn join(&self, other: &TwoPoint<T>) -> TwoLine<T> {
        let x1 = self.x.clone();
        let x2 = other.x.clone();
        let y1 = self.y.clone();
        let y2 = other.y.clone();
        let a = y1.clone() - y2.clone();
        let b = x2.clone() - x1.clone();
        let c = x1 * y2 - x2 * y1;
        TwoLine::new(a, b, c)
    }

    pub fn quadrance(&self, other: &Self) -> T {
        let dx = other.x.clone() - self.x.clone();
        let dy = other.y.clone() - self.y.clone();
        dx.clone()*dx + dy.clone()*dy
    }
}

impl<T> Point for TwoPoint<T>
where
    T: Num,
    T: Clone,
{
    fn is_collinear(&self, a2: &Self, a3: &Self) -> bool {
        let l = self.join(&a2);
        a3.lies_on(&l)
    }
}

impl<T> PartialEq for TwoPoint<T>
where
    T: Num,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.quadrance(&other).is_zero()
    }
}

impl<T> Add<TwoVector<T>> for TwoPoint<T>
where
    T: Num,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn add(self, v: TwoVector<T>) -> TwoPoint<T> {
        TwoPoint {x: self.x.clone() + v.dx(), y: self.y.clone() + v.dy() }
    }
}

impl<T> Mul<Translation<T>> for TwoPoint<T>
where
    T: Num,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn mul(self, t: Translation<T>) -> TwoPoint<T> {
        self + t.vector
    }
}

/// Represents a 2D line a*x+b*y+c=0
#[derive(Debug)]
pub struct TwoLine<T> {
    a : T,
    b : T,
    c : T,
}

impl<T> TwoLine<T> {
    /// Creates a 'TwoLine' without checking for validity of the proportion
    pub const fn new_raw(a: T, b: T, c: T) -> TwoLine<T> {
        TwoLine { a, b, c }
    }
}

impl<T> TwoLine<T>
where
    T: Num,
    T: Clone,
{
    pub fn new(a: T, b: T, c: T) -> TwoLine<T> {
        if a.is_zero() && b.is_zero() {
            panic!("Line has proportion of <0:0:c>");
        }
        TwoLine::new_raw(a, b, c)
    }

    pub fn newpv(p: &TwoPoint<T>, v: &TwoVector<T>) -> TwoLine<T> {
        let vx = v.dx();
        let vy = v.dy();
        let c1 = p.x.clone() * vy.clone() - p.y.clone() * vx.clone();
        let c2 = p.y.clone() * vx.clone() - p.x.clone() * vy.clone();
        TwoLine::new(c1.clone() * vy, c2.clone() * vx, c1*c2)
    }

    pub fn is_parallel(self, other: &Self) -> bool {
        let r = self.a * other.b.clone() - other.a.clone() * self.b;
        r.is_zero()
    }

    pub fn is_perpendicular(self, other: &Self) -> bool {
        let r = self.a * other.a.clone() + self.b * other.b.clone();
        r.is_zero()
    }

    pub fn meet(&self, other: &Self) -> TwoPoint<T> {
        let xn = self.b.clone() * other.c.clone() - other.b.clone() * self.c.clone();
        let yn = self.c.clone() * other.a.clone() - other.c.clone() * self.a.clone();
        let d  = self.a.clone() * other.b.clone() - other.a.clone() * self.b.clone();
        let x = xn / d.clone();
        let y = yn / d;
        TwoPoint { x: x, y: y }
    }

    pub fn is_concurrent(&self, l2: &Self, l3: &Self) -> bool {
        let m = self.meet(&l2);
        m.lies_on(&l3)
    }

    pub fn spread(&self, other: &Self) -> T {
        let cross = self.a.clone() * other.b.clone()
                  - other.a.clone() * self.b.clone();
        let q1 = self.a.clone() * self.a.clone()
               + self.b.clone() * self.b.clone();
        let q2 = other.a.clone() * other.a.clone()
               + other.b.clone() * other.b.clone();
        (cross.clone() * cross) / (q1 * q2)
    }
}

impl<T> PartialEq for TwoLine<T>
where
    T: Num,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let d1 = self.a.clone()*other.b.clone() - other.a.clone()*self.b.clone();
        let d2 = self.b.clone()*other.c.clone() - other.b.clone()*self.c.clone();
        d1.is_zero() & d2.is_zero()
    }
}

impl<T> Mul<Translation<T>> for TwoLine<T>
where
    T: Num,
    T: Clone,
{
    type Output = TwoLine<T>;

    fn mul(self, t: Translation<T>) -> TwoLine<T> {
        TwoLine::new(self.a.clone(), self.b.clone(),
            self.c.clone()-self.a.clone()*t.vector.dx()-self.b.clone()*t.vector.dy())
    }
}

/// Represents a 2D circle
/// given in terms of the proportion (1:0:1:-2a:-2b:c)
/// for the eq. dx²+exy+fy²+ax+by+c=0
/// namly x²+y²-2ax-2by+c=0
/// which can be rewritten as (x-a)²+(y-b)²=K
/// with K = a²+b²-c
/// T is type of the parameters a,b,c
pub struct TwoCircle<T> {
    a: T,
    b: T,
    c: T,
}

impl<T> TwoCircle<T>
where
    T: Num,
    T: Clone,
{
    pub fn center(&self) -> TwoPoint<T> {
        TwoPoint{x: self.a.clone(), y: self.b.clone()}
    }

    pub fn quadrance(&self) -> T {
        let k = self.a.clone()*self.a.clone()
              + self.b.clone()*self.b.clone() - self.c.clone();
        return k;
    }

    pub fn new(center: TwoPoint<T>, quadrance: T) -> TwoCircle<T> {
        let cc = center.x.clone()*center.x.clone()
               + center.y.clone()*center.y.clone() - quadrance;
        TwoCircle {a: center.x, b: center.y, c: cc}
    }

    pub fn lies_on(&self, point: &TwoPoint<T>) -> bool {
        let dx = point.x.clone() - self.a.clone();
        let dy = point.y.clone() - self.b.clone();
        let k = self.quadrance();
        let r = dx.clone() * dx + dy.clone() * dy - k;
        r.is_zero()
    }

}

/// Represents a triangle
/// T is type of the points
pub struct Triangle<T> {
    points: [T; 3],
}

impl<T: Point> Triangle<T> {
    pub fn new(a1: T, a2: T, a3: T) -> Triangle<T> {
        if a1.is_collinear(&a2, &a3) {
            panic!("Triangle cannot have collinear points");
        }
        Triangle { points: [a1, a2, a3] }
    }
}

/// Represents a 2D triangle
/// T is type of the coordinates of the points
pub struct TwoTriangle<T> {
    points: [TwoPoint<T>; 3],
}

impl<T> TwoTriangle<T>
where
    T: Num,
    T: Clone,
    TwoPoint<T>: Point,
    TwoPoint<T>: Clone,
{
    pub fn new(a1: TwoPoint<T>, a2: TwoPoint<T>, a3: TwoPoint<T>)
        -> TwoTriangle<T> {
        if a1.is_collinear(&a2, &a3) {
            panic!("Triangle cannot have collinear points");
        }
        TwoTriangle { points: [a1, a2, a3] }
    }

    pub fn area(&self) -> T {
        let v1 = TwoVector {start: self.points[0].clone(),
                            end: self.points[1].clone()};
        let v2 = TwoVector {start: self.points[0].clone(),
                            end: self.points[2].clone()};
        (v1.dx() * v2.dy() - v2.dx() * v1.dy()) / (T::one() + T::one())
    }
}

/// Represents a 2D vector
#[derive(Debug, Clone)]
pub struct TwoVector<T> {
    start: TwoPoint<T>,
    end: TwoPoint<T>,
}

impl<T> PartialEq for TwoVector<T>
where
    T: Num,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let r = self.end.x.clone() - self.start.x.clone()
              - other.end.x.clone() + other.start.x.clone();
        r.is_zero()
    }
}


impl<T> Mul<T> for TwoVector<T>
where
    T: Mul<T, Output = T>,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, rhs: T) -> TwoVector<T> {
        let s = TwoPoint {x: rhs.clone() * self.start.x.clone(),
                          y: rhs.clone() * self.start.y.clone()};
        let e = TwoPoint {x: rhs.clone() * self.end.x.clone(),
                          y: rhs.clone() * self.end.y.clone()};
        TwoVector {start: s, end: e}
    }
}

impl<T> Add<TwoVector<T>> for TwoVector<T>
where
    T: Add<T, Output = T>,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn add(self, rhs: TwoVector<T>) -> TwoVector<T> {
        let s = TwoPoint {x: rhs.start.x.clone() + self.start.x.clone(),
                          y: rhs.start.y.clone() + self.start.y.clone()};
        let e = TwoPoint {x: rhs.end.x.clone() + self.end.x.clone(),
                          y: rhs.end.y.clone() + self.end.y.clone()};
        TwoVector {start: s, end: e}
    }
}

impl<T> TwoVector<T>
where
    T: Clone,
    T: Zero,
    T: Sub<T, Output = T>,
{
    pub fn new(a: T, b: T) -> TwoVector<T> {
        TwoVector {start: TwoPoint {x: T::zero(), y: T::zero()},
                   end: TwoPoint {x: a, y: b}}
    }

    pub fn dx(&self) -> T {
        let dx = self.end.x.clone() - self.start.x.clone();
        return dx;
    }
    pub fn dy(&self) -> T {
        let dy = self.end.y.clone() - self.start.y.clone();
        return dy;
    }
}

impl<T> TwoVector<T>
where
    T: Num,
    T: Clone,
{
    pub fn area(&self) -> T {
        (self.start.x.clone() * self.end.y.clone()
         - self.end.x.clone() * self.start.y.clone()) / (T::one() + T::one())
    }

}

/// Represents a polygon
pub struct Polygon<T> {
    points: Vec<TwoPoint<T>>,
}

impl<T> Polygon<T>
where
    T: Num,
    T: Clone,
{
    fn area(&self) -> T {
        let mut a = T::zero();

        if self.points.len() < 2 {
            return a;
        }

        for i in 0..(self.points.len()-1) {
            let v = TwoVector {start: self.points[i].clone(),
                               end: self.points[i+1].clone()};
            a = a + v.area();
        }

        let lp = self.points.last();
        match lp {
            Some(p) => {
                    let v = TwoVector {start: p.clone(),
                                       end: self.points[0].clone()};
                    return a + v.area();
                },
            None => a,
        }
    }
}

/// Translation
#[derive(Clone, Debug)]
pub struct Translation<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<Translation<T>> for Translation<T>
where
    T: Num,
    T: Clone,
{
    type Output = Translation<T>;

    fn mul(self, other: Translation<T>) -> Translation<T> {
        let v = self.vector.clone() + other.vector.clone();
        Translation {vector: v}
    }
}

impl<T> PartialEq for Translation<T>
where
    T: Num,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.vector == other.vector
    }
}

/// Rotation
#[derive(Clone)]
pub struct Rotation<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<Rotation<T>> for TwoVector<T>
where
    T: Num,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: Rotation<T>) -> TwoVector<T> {
        let a = r.vector.dx()*self.dx() - r.vector.dy()*self.dy();
        let b = r.vector.dy()*self.dx() + r.vector.dx()*self.dy();
        TwoVector::new(a, b)
    }
}

/// Reflection
#[derive(Clone)]
pub struct Reflection<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<Reflection<T>> for TwoVector<T>
where
    T: Num,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: Reflection<T>) -> TwoVector<T> {
        let a = r.vector.dx()*self.dx() + r.vector.dy()*self.dy();
        let b = r.vector.dy()*self.dx() - r.vector.dx()*self.dy();
        TwoVector::new(a, b)
    }
}

/// Represents a quidrilateral
/// T is type of the points
pub struct Quadrilateral<T> {
    points: [T; 4],
}

impl<T: Point> Quadrilateral<T> {
    pub fn new(a1: T, a2: T, a3: T, a4: T) -> Quadrilateral<T> {
        if a1.is_collinear(&a2, &a3) {
            panic!("Quadrilateral cannot have collinear points");
        }
        if a2.is_collinear(&a3, &a4) {
            panic!("Quadrilateral cannot have collinear points");
        }
        if a3.is_collinear(&a4, &a1) {
            panic!("Quadrilateral cannot have collinear points");
        }
        if a4.is_collinear(&a1, &a1) {
            panic!("Quadrilateral cannot have collinear points");
        }

        Quadrilateral { points: [a1, a2, a3, a4] }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio, Rational};

    #[test]
    #[should_panic]
    fn invalid_two_line() {
        TwoLine::new(0,0,1);
    }
    #[test]
    fn point_lies_on_line() {
        let a = TwoPoint {x: Ratio::new( 6, 1),
            y: Ratio::new( -3, 2)};
        let l = TwoLine::new( Ratio::new(1, 1),
            Ratio::new(2, 1), Ratio::new(-3, 1) );
        assert!(a.lies_on(&l));
    }
    #[test]
    fn points_lie_on_join() {
        let a = TwoPoint {x: Ratio::new( 6, 1),
            y: Ratio::new( -3, 2) };
        let b = TwoPoint {x: Ratio::new( 4, 3),
            y: Ratio::new( 3, 5) };
        let l = a.join(&b);
        assert!(a.lies_on(&l));
        assert!(b.lies_on(&l));
    }
    #[test]
    fn parallel_lines() {
        let l1 = TwoLine::new(3, 4, -1);
        let l2 = TwoLine::new(6, 8, 5);
        assert!(l1.is_parallel(&l2));
    }
    #[test]
    fn perpendicular_lines() {
        let l1 = TwoLine::new(3, 4, -1);
        let l2 = TwoLine::new(-4, 3, 2);
        assert!(l1.is_perpendicular(&l2));
    }
    #[test]
    fn meet_lies_on_lines() {
        let l1 = TwoLine::new(Rational::new(3, 1),
            Rational::new(4, 1), Rational::new(-1, 1));
        let l2 = TwoLine::new(Rational::new(-4, 1),
            Rational::new(3, 1), Rational::new(2, 1));
        let a = l1.meet(&l2);
        assert!(a.lies_on(&l1));
        assert!(a.lies_on(&l2));
    }
    #[test]
    fn create_triangle() {
        let p1 = TwoPoint {x: Ratio::new(-3,1), y: Ratio::new(4,1)};
        let p2 = TwoPoint {x: Ratio::new(4,1), y: Ratio::new(5,1)};
        let p3 = TwoPoint {x: Ratio::new(1,1), y: Ratio::new(1,1)};
        let _t = Triangle::new(p1, p2, p3);
    }
    #[test]
    fn equal_points() {
        let p1 = TwoPoint {x: Ratio::new(-3,1), y: Ratio::new(4,1)};
        let p2 = TwoPoint {x: Ratio::new(-3,1), y: Ratio::new(4,1)};
        assert_eq!(p1,p2);
    }
    #[test]
    fn quadrance() {
        let p1 = TwoPoint {x: Ratio::new(2,1), y: Ratio::new(1,1)};
        let p2 = TwoPoint {x: Ratio::new(6,1), y: Ratio::new(2,1)};
        let q  = Ratio::new(17,1);
        assert_eq!(p1.quadrance(&p2),q);
    }
    #[test]
    fn spread() {
        let l1 = TwoLine::new(Rational::new(3,1),
            Rational::new(2,1), Rational::new(1,1));
        let l2 = TwoLine::new(Rational::new(1,1),
            Rational::new(-4,1), Rational::new(2,1));
        let s  = Rational::new(196, 221);
        assert_eq!(l1.spread(&l2),s);
    }
    #[test]
    fn point_lies_on_circle() {
        let c = TwoCircle {a: Rational::new(0,1), b: Rational::new(0,1),
            c: Rational::new(-2,1)};
        let p = TwoPoint { x: Rational::new(1,1), y: Rational::new(1,1) };
        assert!(c.lies_on(&p));
    }
    #[test]
    fn equal_vectors() {
        let a1 = TwoPoint {x: Ratio::new(1,1), y: Ratio::new(1,1)};
        let a2 = TwoPoint {x: Ratio::new(3,1), y: Ratio::new(4,1)};
        let a3 = TwoPoint {x: Ratio::new(4,1), y: Ratio::new(2,1)};
        let a4 = TwoPoint {x: Ratio::new(6,1), y: Ratio::new(5,1)};
        let v1 = TwoVector {start: a1, end: a2};
        let v2 = TwoVector {start: a3, end: a4};
        assert_eq!(v1,v2);
    }
    #[test]
    fn scalar_vector_mul() {
        let s = Ratio::new(3,1);
        let v1 = TwoVector::new(Ratio::new(2,1), Ratio::new(-4,1));
        let v2 = TwoVector::new(Ratio::new(6,1), Ratio::new(-12,1));
        assert_eq!(v1*s,v2);
    }
    #[test]
    fn add_vectors() {
        let v1 = TwoVector::new(Ratio::new(2,1), Ratio::new(3,1));
        let v2 = TwoVector::new(Ratio::new(2,1), Ratio::new(-4,1));
        let v3 = TwoVector::new(Ratio::new(4,1), Ratio::new(-1,1));
        assert_eq!(v1+v2,v3);
    }
    #[test]
    fn point_vector_line() {
        let a1 = TwoPoint {x: Ratio::new(3,1), y: Ratio::new(2,1)};
        let v1 = TwoVector::new(Ratio::new(1,1), Ratio::new(2,1));
        let a2 = a1.clone() + v1.clone();
        let a3 = a1.clone() + v1.clone() * Ratio::new(-12,5);
        let l1 = TwoLine::newpv(&a1, &v1);
        assert!(a1.lies_on(&l1));
        assert!(a2.lies_on(&l1));
        assert!(a3.lies_on(&l1));
    }
    #[test]
    fn triangle_area() {
        let a1 = TwoPoint {x: Ratio::new(1,1), y: Ratio::new(1,1)};
        let a2 = TwoPoint {x: Ratio::new(7,1), y: Ratio::new(3,1)};
        let a3 = TwoPoint {x: Ratio::new(2,1), y: Ratio::new(5,1)};
        let t = TwoTriangle::new(a1, a2, a3);
        assert_eq!(t.area(),Ratio::new(11,1));
    }
    #[test]
    fn vector_area() {
        let a1 = TwoPoint {x: Ratio::new(3,1), y: Ratio::new(6,1)};
        let a2 = TwoPoint {x: Ratio::new(7,1), y: Ratio::new(5,1)};
        let v1 = TwoVector {start: a1, end: a2};
        assert_eq!(v1.area(),Ratio::new(-27,2));
    }
    #[test]
    fn polygon_area() {
        let a1 = TwoPoint {x: Ratio::new(3,1), y: Ratio::new(6,1)};
        let a2 = TwoPoint {x: Ratio::new(7,1), y: Ratio::new(5,1)};
        let a3 = TwoPoint {x: Ratio::new(7,1), y: Ratio::new(1,1)};
        let a4 = TwoPoint {x: Ratio::new(2,1), y: Ratio::new(2,1)};
        let a5 = TwoPoint {x: Ratio::new(1,1), y: Ratio::new(4,1)};
        let p1 = Polygon {points: vec![a1, a2, a3, a4, a5]};
        assert_eq!(p1.area(), Ratio::new(-43,2));
    }
    #[test]
    fn translation_of_2d_point() {
        let a1 = TwoPoint {x: Ratio::new(3,1), y: Ratio::new(6,1)};
        let v1 = TwoVector::new(Ratio::new(1,1), Ratio::new(1,1));
        let t1 = Translation {vector: v1.clone()};
        assert_eq!(a1.clone()*t1,a1+v1);
    }
    #[test]
    fn translation_of_2d_line() {
        let p1 = TwoPoint {x: Ratio::new(8,1), y: Ratio::new(3,1)};
        let p2 = TwoPoint {x: Ratio::new(9,1), y: Ratio::new(1,1)};
        let l1 = p1.join(&p2);
        let v1 = TwoVector::new(Ratio::new(3,1), Ratio::new(1,1));
        let t1 = Translation {vector: v1.clone()};
        let p3 = p1*t1.clone();
        let p4 = p2*t1.clone();
        let m1 = p3.join(&p4);
        assert_eq!(l1*t1,m1);
    }
    #[test]
    fn addition_of_translations() {
        let v1 = TwoVector::new(Ratio::new(3,1), Ratio::new(1,1));
        let v2 = TwoVector::new(Ratio::new(2,1), Ratio::new(6,1));
        let v3 = v1.clone() + v2.clone();
        let t1 = Translation {vector: v1};
        let t2 = Translation {vector: v2};
        let t3 = Translation {vector: v3};
        assert_eq!(t1*t2,t3);
    }
    #[test]
    fn rotation() {
        let v1 = TwoVector::new(Ratio::new(1,1), Ratio::new(0,1));
        let v2 = TwoVector::new(Ratio::new(0,1), Ratio::new(1,1));
        let v3 = TwoVector::new(Ratio::new(2,1), Ratio::new(6,1));
        let r1 = Rotation {vector: v3.clone()};
        assert_eq!(v1*r1.clone(),v3);
        let v4 = TwoVector::new(Ratio::new(-6,1), Ratio::new(2,1));
        assert_eq!(v2*r1,v4);
    }
    #[test]
    fn reflection() {
        let v1 = TwoVector::new(Ratio::new(1,1), Ratio::new(0,1));
        let v2 = TwoVector::new(Ratio::new(0,1), Ratio::new(1,1));
        let v3 = TwoVector::new(Ratio::new(2,1), Ratio::new(6,1));
        let r1 = Reflection {vector: v3.clone()};
        assert_eq!(v1*r1.clone(),v3);
        let v4 = TwoVector::new(Ratio::new(6,1), Ratio::new(-2,1));
        assert_eq!(v2*r1,v4);
    }
}
