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


use num::{Num, Zero, One, Signed};
use std::ops::{Mul, Add, Sub, Div};
use linalg::{Matrix, RowVector, ColumnVector};


pub trait Point {
    fn is_collinear(&self, a2: &Self, a3: &Self) -> bool;
}


pub fn archimedes<T>(x: T, y: T, z: T) -> T
where
    T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> + One,
    T: Clone,
{
    let two = T::one() + T::one();
    let xyz = x.clone() + y.clone() + z.clone();
    let xyz2 = xyz.clone() * xyz;
    xyz2 - two * (x.clone()*x + y.clone()*y + z.clone()*z)
}


pub fn half_slope<T>(h: T) -> TwoPoint<T>
where
    T: One + Zero,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    let two = T::one() + T::one();
    let h2 = h.clone() * h.clone();
    let x = (T::one() - h2.clone()) / (T::one() + h2.clone());
    let y = two*h.clone() / (T::one() + h2);
    TwoPoint::new(x,y)
}

pub fn circle_sum<T>(h1: T, h2: T) -> T
where
    T: One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    (h1.clone() + h2.clone()) / (T::one() - h1*h2)
}

/// Represents a 1D point
#[derive(Debug, Clone)]
pub struct OnePoint<T> {
    pub x : T,
}

impl<T> OnePoint<T>
{
    pub fn new(x : T) -> OnePoint<T> {
        OnePoint{ x: x }
    }
}

impl<T> OnePoint<T>
where
    T: Sub<T, Output = T>,
    T: Signed,
    T: Clone,
{
    pub fn length(&self, other: &Self) -> T {
        let d = self.x.clone() - other.x.clone();
        return d.abs();
    }
}

impl<T> OnePoint<T>
where
    T: Sub<T, Output = T>,
    T: Clone,
{
    /// displacement
    pub fn disp(&self, other: &Self) -> T {
        other.x.clone() - self.x.clone()
    }
}

impl<T> OnePoint<T>
where
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Clone,
{
    pub fn quadrance(&self, other: &Self) -> T {
        let dx = other.x.clone() - self.x.clone();
        dx.clone()*dx
    }
}

/// Represents a 2D point
#[derive(Debug, Clone)]
pub struct TwoPoint<T> {
    pub x : T,
    pub y : T,
}

impl<T> TwoPoint<T>
{
    pub fn new(x: T, y: T) -> TwoPoint<T> {
        TwoPoint {x: x, y: y}
    }
}

impl<T> TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
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
        let v = TwoVector { start: self.clone(), end: other.clone() };
        v.quadrance_blue()
    }

    pub fn quadrance_red(&self, other: &Self) -> T {
        let v = TwoVector { start: self.clone(), end: other.clone() };
        v.quadrance_red()
    }

    pub fn quadrance_green(&self, other: &Self) -> T {
        let v = TwoVector { start: self.clone(), end: other.clone() };
        v.quadrance_green()
    }

    pub fn quadrance_metric(&self, other: &Self, metric: &Matrix<T>) -> T {
        let v = TwoVector { start: self.clone(), end: other.clone() };
        v.quadrance_metric(&metric)
    }
}

impl<T> Point for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn is_collinear(&self, a2: &Self, a3: &Self) -> bool {
        let l = self.join(&a2);
        a3.lies_on(&l)
    }
}

impl<T> PartialEq for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.quadrance(&other).is_zero()
    }
}

impl<T> Add<TwoVector<T>> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn add(self, v: TwoVector<T>) -> TwoPoint<T> {
        TwoPoint {x: self.x.clone() + v.dx(), y: self.y.clone() + v.dy() }
    }
}

impl<T> Add<TwoPoint<T>> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn add(self, other: TwoPoint<T>) -> TwoPoint<T> {
        let nx =   self.x.clone() * other.y.clone()
                 + self.y.clone() * other.x.clone();
        let ny = self.y.clone() * other.y.clone();
        TwoPoint {x: nx, y: ny }
    }
}

impl<T> Mul<Translation<T>> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn mul(self, t: Translation<T>) -> TwoPoint<T> {
        self + t.vector
    }
}

impl<T> Mul<T> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn mul(self, scalar: T) -> TwoPoint<T> {
        let nx = self.x.clone() * scalar.clone();
        let ny = self.y.clone() * scalar;
        TwoPoint {x: nx, y: ny }
    }
}

impl<T> Mul<TwoPoint<T>> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn mul(self, other: TwoPoint<T>) -> TwoPoint<T> {
        let u = TwoVector { start: TwoPoint::zero(), end: other };
        let v = TwoVector { start: TwoPoint::zero(), end: self };
        let r = Rotation { vector: u };
        let w = v * r;
        w.end
    }
}

impl<T> Sub<TwoPoint<T>> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn sub(self, other: TwoPoint<T>) -> TwoPoint<T> {
        let nx =   self.x.clone() * other.y.clone()
                 - self.y.clone() * other.x.clone();
        let ny = self.y.clone() * other.y.clone();
        TwoPoint {x: nx, y: ny }
    }
}

impl<T> Div<TwoPoint<T>> for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    type Output = TwoPoint<T>;

    fn div(self, other: TwoPoint<T>) -> TwoPoint<T> {
        let nx = self.x.clone() * other.y.clone();
        let ny = self.y.clone() * other.x.clone();
        TwoPoint {x: nx, y: ny }
    }
}

impl<T> Zero for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn zero() -> TwoPoint<T> {
        TwoPoint { x: T::zero(), y: T::zero() }
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero()
    }
}

impl<T> One for TwoPoint<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    fn one() -> TwoPoint<T> {
        TwoPoint { x: T::one(), y: T::one() }
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
    T: Zero,
{
    pub fn new(a: T, b: T, c: T) -> TwoLine<T> {
        if a.is_zero() && b.is_zero() {
            panic!("Line has proportion of <0:0:c>");
        }
        TwoLine::new_raw(a, b, c)
    }
}

impl<T> TwoLine<T>
where
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    pub fn newpv(p: &TwoPoint<T>, v: &TwoVector<T>) -> TwoLine<T> {
        let vx = v.dx();
        let vy = v.dy();
        let c1 = p.x.clone() * vy.clone() - p.y.clone() * vx.clone();
        let c2 = p.y.clone() * vx.clone() - p.x.clone() * vy.clone();
        TwoLine::new(c1.clone() * vy, c2.clone() * vx, c1*c2)
    }
}

impl<T> TwoLine<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
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
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
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

    pub fn curvature(&self) -> T {
        T::one() / self.quadrance()
    }

    pub fn quadrance_red(&self) -> T {
        let k = self.a.clone()*self.a.clone()
              - self.b.clone()*self.b.clone() - self.c.clone();
        return k;
    }

    pub fn new(center: TwoPoint<T>, quadrance: T) -> TwoCircle<T> {
        let cc = center.x.clone()*center.x.clone()
               + center.y.clone()*center.y.clone() - quadrance;
        TwoCircle {a: center.x, b: center.y, c: cc}
    }

    pub fn new_red(center: TwoPoint<T>, quadrance: T) -> TwoCircle<T> {
        let cc = center.x.clone()*center.x.clone()
               - center.y.clone()*center.y.clone() - quadrance;
        TwoCircle {a: center.x, b: center.y, c: cc}
    }

    pub fn lies_on(&self, point: &TwoPoint<T>) -> bool {
        let dx = point.x.clone() - self.a.clone();
        let dy = point.y.clone() - self.b.clone();
        let k = self.quadrance();
        let r = dx.clone() * dx + dy.clone() * dy - k;
        r.is_zero()
    }

    pub fn lies_on_red(&self, point: &TwoPoint<T>) -> bool {
        let dx = point.x.clone() - self.a.clone();
        let dy = point.y.clone() - self.b.clone();
        let k = self.quadrance_red();
        let r = dx.clone() * dx - dy.clone() * dy - k;
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

    pub fn quadrea(&self) -> T {
        let two = T::one() + T::one();
        let four = two.clone() * two;
        let sixteen = four.clone() * four;
        let a = self.area();
        sixteen * a.clone() * a
    }

    pub fn circumquadrance(&self) -> T {
        let q0 = self.points[1].quadrance(&self.points[2]);
        let q1 = self.points[0].quadrance(&self.points[2]);
        let q2 = self.points[0].quadrance(&self.points[1]);

        q0*q1*q2/self.quadrea()
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
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
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

impl<T> Sub<TwoVector<T>> for TwoVector<T>
where
    T: Sub<T, Output = T>,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn sub(self, rhs: TwoVector<T>) -> TwoVector<T> {
        let s = TwoPoint {x: self.start.x.clone() - rhs.start.x.clone(),
                          y: self.start.y.clone() - rhs.start.y.clone()};
        let e = TwoPoint {x: self.end.x.clone() - rhs.end.x.clone(),
                          y: self.end.y.clone() - rhs.end.y.clone()};
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

    pub fn new0e(end: TwoPoint<T>) -> TwoVector<T> {
        TwoVector {start: TwoPoint {x: T::zero(), y: T::zero()}, end: end }
    }

    pub fn newse(start: TwoPoint<T>, end: TwoPoint<T>) -> TwoVector<T> {
        TwoVector {start: start, end: end }
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
    T: Clone,
    T: Zero,
    T: One,
    T: Sub<T, Output = T>,
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
}

impl<T> TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    pub fn dot_blue(&self, other: &Self) -> T {
        self.dot_metric(&other, &TwoVector::metric_blue())
    }
    pub fn dot_red(&self, other: &Self) -> T {
        self.dot_metric(&other, &TwoVector::metric_red())
    }
    pub fn dot_green(&self, other: &Self) -> T {
        self.dot_metric(&other, &TwoVector::metric_green())
    }
    pub fn dot_metric(&self, other: &Self, metric: &Matrix<T>) -> T {
        let a = self.dx();
        let b = self.dy();
        let c = other.dx();
        let d = other.dy();

        let v1 = RowVector::new(vec![a, b]);
        let v2 = ColumnVector::new(vec![c, d]);

        v1*(*metric).clone()*v2
    }
    pub fn cross(&self, other: &Self) -> T {
        let a = self.dx();
        let b = self.dy();
        let c = other.dx();
        let d = other.dy();

        let v1 = RowVector::new(vec![a, b]);
        let m1 = Matrix::new(vec![vec![T::zero(), T::one()],
                                  vec![T::zero() - T::one(), T::zero()]]);
        let v2 = ColumnVector::new(vec![c, d]);

        v1*m1*v2
    }
    pub fn area(&self) -> T {
        (self.start.x.clone() * self.end.y.clone()
         - self.end.x.clone() * self.start.y.clone()) / (T::one() + T::one())
    }
    pub fn spread(&self, other: &Self) -> T {
        self.spread_blue(&other)
    }
    pub fn quadrance_blue(&self) -> T {
        self.dot_blue(&self)
    }
    pub fn quadrance_red(&self) -> T {
        self.dot_red(&self)
    }
    pub fn quadrance_green(&self) -> T {
        self.dot_green(&self)
    }
    pub fn quadrance_metric(&self, metric: &Matrix<T>) -> T {
        self.dot_metric(&self, &metric)
    }
    pub fn spread_blue(&self, other: &Self) -> T {
        self.spread_metric(&other, &TwoVector::metric_blue())
    }
    pub fn spread_red(&self, other: &Self) -> T {
        self.spread_metric(&other, &TwoVector::metric_red())
    }
    pub fn spread_green(&self, other: &Self) -> T {
        self.spread_metric(&other, &TwoVector::metric_green())
    }
    pub fn spread_metric(&self, other: &Self, metric: &Matrix<T>) -> T {
        let dot = self.dot_metric(&other, &metric);
        let dot2 = dot.clone() * dot;
        T::one() - dot2 / (self.quadrance_metric(&metric) * other.quadrance_metric(&metric))
    }
    pub fn proj(&self, other: &Self, metric: &Matrix<T>) -> TwoVector<T> {
        let vu = self.dot_metric(&other, &metric);
        let u = self.quadrance_metric(&metric);
        let s = vu/(u.clone() * u);
        self.clone() * s
    }
    pub fn refl(&self, other: &Self, metric: &Matrix<T>) -> TwoVector<T> {
        let two = T::one() + T::one();
        let vu = self.dot_metric(&other, &metric);
        let u = self.quadrance_metric(&metric);
        let s = vu/(u.clone() * u);
        self.clone() * two * s - other.clone()
    }
}


impl<T> TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: PartialEq,
    T: Clone,
{
    pub fn is_perpendicular_metric(&self, other: &Self, metric: &Matrix<T>) -> bool {
        self.dot_metric(&other, &metric) == T::zero()
    }
    pub fn is_parallel(&self, other: &Self) -> bool {
        self.is_parallel_metric(&other, &TwoVector::metric_blue())
    }
    pub fn is_parallel_metric(&self, other: &Self, metric: &Matrix<T>) -> bool {
        let qu = self.quadrance_metric(&metric);
        let qv = other.quadrance_metric(&metric);
        let uv = (*self).clone() - (*other).clone();
        let quv = uv.quadrance_metric(&metric);
        archimedes(qu,qv,quv) == T::zero()
    }
}

/// Represents a polygon
pub struct Polygon<T> {
    points: Vec<TwoPoint<T>>,
}

impl<T> Polygon<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: PartialEq,
    T: Clone,
{
    pub fn new(p: Vec<TwoPoint<T>>) -> Polygon<T> {
        Polygon { points: p }
    }

    pub fn area(&self) -> T {
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
#[derive(Clone,Debug)]
pub struct Rotation<T> {
    vector: TwoVector<T>,
}

impl<T> Rotation<T>
{
    pub fn new(vector: TwoVector<T>) -> Rotation<T> {
        Rotation { vector: vector }
    }
}

impl<T> PartialEq for Rotation<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.vector == other.vector
    }
}

impl<T> Mul<Rotation<T>> for TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: Rotation<T>) -> TwoVector<T> {
        let a = r.vector.dx();
        let b = r.vector.dy();
        let x = self.dx();
        let y = self.dy();
        let v = RowVector::new(vec![x,y]);
        let m = Matrix::new(vec![vec![a.clone(),b.clone()],
                                 vec![T::zero()-b,a]]);
        let u = v*m;
        TwoVector::new(u.get(0), u.get(1))
    }
}

impl<T> Mul<Rotation<T>> for Rotation<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = Rotation<T>;

    fn mul(self, other: Rotation<T>) -> Rotation<T> {
        let a = self.vector.dx();
        let b = self.vector.dy();
        let c = other.vector.dx();
        let d = other.vector.dy();
        let m1 = Matrix::new(vec![vec![a.clone(),b.clone()],
                                  vec![T::zero()-b,a]]);
        let m2 = Matrix::new(vec![vec![c.clone(),d.clone()],
                                  vec![T::zero()-d,c]]);
        let m3 = m1*m2;
        Rotation { vector: TwoVector::new(m3.get(0,0), m3.get(0,1)) }
    }
}

impl<T> Mul<Reflection<T>> for Rotation<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = Reflection<T>;

    fn mul(self, other: Reflection<T>) -> Reflection<T> {
        let a = self.vector.dx();
        let b = self.vector.dy();
        let c = other.vector.dx();
        let d = other.vector.dy();
        let m1 = Matrix::new(vec![vec![a.clone(),b.clone()],
                                  vec![T::zero()-b,a]]);
        let m2 = Matrix::new(vec![vec![c.clone(),d.clone()],
                                  vec![d,T::zero()-c]]);
        let m3 = m1*m2;
        Reflection { vector: TwoVector::new(m3.get(0,0), m3.get(0,1)) }
    }
}

/// Red Rotation
#[derive(Clone,Debug)]
pub struct RotationRed<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<RotationRed<T>> for TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: RotationRed<T>) -> TwoVector<T> {
        let a = r.vector.dx();
        let b = r.vector.dy();
        let x = self.dx();
        let y = self.dy();
        let v = RowVector::new(vec![x,y]);
        let m = Matrix::new(vec![vec![a.clone(),b.clone()],
                                 vec![b,a]]);
        let u = v*m;
        TwoVector::new(u.get(0), u.get(1))
    }
}

/// Green Rotation
#[derive(Clone,Debug)]
pub struct RotationGreen<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<RotationGreen<T>> for TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: RotationGreen<T>) -> TwoVector<T> {
        let a = r.vector.dx();
        let b = r.vector.dy();
        let x = self.dx();
        let y = self.dy();
        let v = RowVector::new(vec![x,y]);
        let m = Matrix::new(vec![vec![a,T::zero()],
                                 vec![T::zero(),b]]);
        let u = v*m;
        TwoVector::new(u.get(0), u.get(1))
    }
}

/// Reflection
#[derive(Clone,Debug)]
pub struct Reflection<T> {
    vector: TwoVector<T>,
}

impl<T> Reflection<T>
{
    pub fn new(vector: TwoVector<T>) -> Reflection<T> {
        Reflection { vector: vector }
    }
}

impl<T> PartialEq for Reflection<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.vector == other.vector
    }
}

impl<T> Mul<Reflection<T>> for TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: Reflection<T>) -> TwoVector<T> {
        let a = r.vector.dx();
        let b = r.vector.dy();
        let x = self.dx();
        let y = self.dy();
        let v = RowVector::new(vec![x,y]);
        let m = Matrix::new(vec![vec![a.clone(),b.clone()],
                                 vec![b,T::zero()-a]]);
        let u = v*m;
        TwoVector::new(u.get(0), u.get(1))
    }
}

impl<T> Mul<Reflection<T>> for Reflection<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = Rotation<T>;

    fn mul(self, other: Reflection<T>) -> Rotation<T> {
        let a = self.vector.dx();
        let b = self.vector.dy();
        let c = other.vector.dx();
        let d = other.vector.dy();
        let m1 = Matrix::new(vec![vec![a.clone(),b.clone()],
                                  vec![b,T::zero()-a]]);
        let m2 = Matrix::new(vec![vec![c.clone(),d.clone()],
                                  vec![d,T::zero()-c]]);
        let m3 = m1*m2;
        Rotation { vector: TwoVector::new(m3.get(0,0), m3.get(0,1)) }
    }
}

impl<T> Mul<Rotation<T>> for Reflection<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = Reflection<T>;

    fn mul(self, other: Rotation<T>) -> Reflection<T> {
        let a = self.vector.dx();
        let b = self.vector.dy();
        let c = other.vector.dx();
        let d = other.vector.dy();
        let m1 = Matrix::new(vec![vec![a.clone(),b.clone()],
                                  vec![b,T::zero()-a]]);
        let m2 = Matrix::new(vec![vec![c.clone(),d.clone()],
                                  vec![T::zero()-d,c]]);
        let m3 = m1*m2;
        Reflection { vector: TwoVector::new(m3.get(0,0), m3.get(0,1)) }
    }
}

/// Red Reflection
#[derive(Clone,Debug)]
pub struct ReflectionRed<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<ReflectionRed<T>> for TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: ReflectionRed<T>) -> TwoVector<T> {
        let a = r.vector.dx();
        let b = r.vector.dy();
        let x = self.dx();
        let y = self.dy();
        let v = RowVector::new(vec![x,y]);
        let m = Matrix::new(vec![vec![a.clone(),b.clone()],
                                 vec![T::zero()-b,T::zero()-a]]);
        let u = v*m;
        TwoVector::new(u.get(0), u.get(1))
    }
}

/// Green Reflection
#[derive(Clone,Debug)]
pub struct ReflectionGreen<T> {
    vector: TwoVector<T>,
}

impl<T> Mul<ReflectionGreen<T>> for TwoVector<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    type Output = TwoVector<T>;

    fn mul(self, r: ReflectionGreen<T>) -> TwoVector<T> {
        let a = r.vector.dx();
        let b = r.vector.dy();
        let x = self.dx();
        let y = self.dy();
        let v = RowVector::new(vec![x,y]);
        let m = Matrix::new(vec![vec![T::zero(),b],
                                 vec![a,T::zero()]]);
        let u = v*m;
        TwoVector::new(u.get(0), u.get(1))
    }
}

/// Represents a quadrilateral
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
        if a4.is_collinear(&a1, &a2) {
            panic!("Quadrilateral cannot have collinear points");
        }

        Quadrilateral { points: [a1, a2, a3, a4] }
    }
}

impl<T> Quadrilateral<TwoPoint<T>>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Zero,
    T: One,
    T: Clone,
{
    pub fn quadrea(&self) -> T {
        let four = T::one() + T::one() + T::one() + T::one();
        let q01 = self.points[0].quadrance(&self.points[1]);
        let q02 = self.points[0].quadrance(&self.points[2]);
        let q03 = self.points[0].quadrance(&self.points[3]);
        let q12 = self.points[1].quadrance(&self.points[2]);
        let q13 = self.points[1].quadrance(&self.points[3]);
        let q23 = self.points[2].quadrance(&self.points[3]);
        let s = q01+q23-q03-q12;
        return four * q02*q13 - s.clone()*s;
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
    fn addition_of_points() {
        let a1 = TwoPoint {x: 2, y: 3};
        let a2 = TwoPoint {x: 4, y: 5};
        let a3 = TwoPoint {x: 22, y: 15};
        assert_eq!(a1+a2,a3);

        let b1 = TwoPoint {x: 3, y: -1};
        let b2 = TwoPoint {x: 0, y: 2};
        let b3 = TwoPoint {x: 6, y: -2};
        assert_eq!(b1+b2,b3);

        let c1 = TwoPoint {x: 6, y: 0};
        let c2 = TwoPoint {x: 3, y: 1};
        let c3 = TwoPoint {x: 6, y: 0};
        assert_eq!(c1+c2,c3);
    }
    #[test]
    fn multiplication_of_points() {
        let zero = Ratio::new(0,1);
        let one = Ratio::new(1,1);
        let pi = TwoPoint {x: one, y: zero};
        let a = Ratio::new(15,17);
        let b = Ratio::new(8,17);
        let c = Ratio::new(-7,25);
        let d = Ratio::new(24,25);
        let pa = TwoPoint {x: a, y: b};
        let pb = TwoPoint {x: c, y: d};
        let pc = pa.clone()*pb.clone();
        let u = TwoVector { start: pa, end: pb };
        let v = TwoVector { start: pi, end: pc };
        let ic = TwoVector::new(a*c-b*d-one, a*b+b*c);
        assert_eq!(v,ic);
        assert!(u.is_parallel(&v));
    }
    #[test]
    fn subtraction_of_points() {
        let a1 = TwoPoint {x: 1, y: -1};
        let a2 = TwoPoint {x: 2, y: 0};
        let a3 = TwoPoint {x: -1, y: 3};
        let a4 = TwoPoint {x: 6, y: 0};
        assert_eq!(a1.clone()-(a2.clone()-a3.clone()),a4.clone());
        assert_eq!(a1.clone()-a4.clone(),a4.clone());
        assert_eq!((a1.clone()-a2.clone())+a3.clone(),a4.clone());
        assert_eq!(a2.clone()+a3.clone(),a4.clone());

        let b1 = TwoPoint {x: 2, y: 3};
        let b2 = TwoPoint {x: 1, y: 4};
        let b3 = TwoPoint {x: -1, y: 1};
        let b4 = TwoPoint {x: 17, y: 12};
        assert_eq!(b1.clone()-(b2.clone()+b3.clone()),b4.clone());
        assert_eq!((b1.clone()-b2.clone())-b3.clone(),b4.clone());
    }
    #[test]
    fn distributive_laws_for_points2() {
        let a1 = TwoPoint {x: 2, y: 3};
        let a2 = TwoPoint {x: 1, y: -1};
        let a3 = TwoPoint {x: 5, y: 4};
        let a4 = TwoPoint {x: -100, y: -75};
        assert_eq!(((a1.clone()-a2.clone())/a3.clone())*5,a4.clone());
        assert_eq!((a1.clone()/a3.clone())-(a2.clone()/a3.clone()),a4.clone());
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
    fn spread_between_vectors() {
        let o = TwoPoint {x: Ratio::new(0,1), y: Ratio::new(0,1)};
        let a = TwoPoint {x: Ratio::new(4,1), y: Ratio::new(2,1) };
        let b = TwoPoint {x: Ratio::new(-1,1), y: Ratio::new(3,1)};
        let oa = TwoVector {start: o.clone(), end: a.clone()};
        let ob = TwoVector {start: o, end: b.clone()};
        let ab = TwoVector {start: a, end: b};
        let s = oa.spread(&ob);
        let t = ob.spread(&ab);
        let r = oa.spread(&ab);
        assert_eq!(s,Ratio::new(49,50));
        assert_eq!(t,Ratio::new(49,65));
        assert_eq!(r,Ratio::new(49,130));
    }
    #[test]
    fn blue_red_green_quadrance_of_vectors() {
        let o = TwoPoint {x: Ratio::new(0,1), y: Ratio::new(0,1)};
        let a = TwoPoint {x: Ratio::new(4,1), y: Ratio::new(2,1) };
        let oa = TwoVector {start: o.clone(), end: a.clone()};
        let qb2 = oa.quadrance_blue() * oa.quadrance_blue();
        let qr2 = oa.quadrance_red() * oa.quadrance_red();
        let qg2 = oa.quadrance_green() * oa.quadrance_green();
        assert_eq!(qb2,qr2+qg2);
    }
    #[test]
    fn blue_red_green_spread_of_vectors() {
        let v = TwoVector::new(Ratio::new(3,1),Ratio::new(1,1));
        let u = TwoVector::new(Ratio::new(1,1),Ratio::new(4,1));
        let sb = v.spread_blue(&u);
        let sr = v.spread_red(&u);
        let sg = v.spread_green(&u);
        let one = Ratio::new(1,1);
        let two = Ratio::new(2,1);
        let si = one/sb + one/sr + one/sg;
        assert_eq!(sb,Ratio::new(121,170));
        assert_eq!(sr,Ratio::new(121,120));
        assert_eq!(sg,Ratio::new(-121,48));
        assert_eq!(si,two);
    }
    #[test]
    fn quadrance_with_general_metric() {
        let m = Matrix::new(vec![vec![Ratio::new(3,8),Ratio::new(-1,4)],
                                 vec![Ratio::new(-1,4),Ratio::new(1,2)]]);
        let u = TwoVector::new(Ratio::new(2,1),Ratio::new(1,1));
        let v = TwoVector::new(Ratio::new(-2,1),Ratio::new(-1,1));
        let u2 = TwoVector::new(Ratio::new(4,1),Ratio::new(2,1));
        let v2 = TwoVector::new(Ratio::new(-1,1),Ratio::new(-1,2));
        let qu = u.quadrance_metric(&m);
        let qv = v.quadrance_metric(&m);
        let qu2 = u2.quadrance_metric(&m);
        let qv2 = v2.quadrance_metric(&m);
        assert_eq!(qu,Ratio::new(1,1));
        assert_eq!(qv,Ratio::new(1,1));
        assert_eq!(qu2,Ratio::new(4,1));
        assert_eq!(qv2,Ratio::new(1,4));
    }
    #[test]
    fn perpendicular_vectors_with_general_metric() {
        let m = Matrix::new(vec![vec![Ratio::new(3,8),Ratio::new(-1,4)],
                                 vec![Ratio::new(-1,4),Ratio::new(1,2)]]);
        let u = TwoVector::new(Ratio::new(2,1),Ratio::new(1,1));
        let v = TwoVector::new(Ratio::new(0,1),Ratio::new(1,1));
        assert!(u.is_perpendicular_metric(&v, &m));
    }
    #[test]
    fn spread_between_vectors_with_general_metric() {
        let m = Matrix::new(vec![vec![Ratio::new(3,8),Ratio::new(-1,4)],
                                 vec![Ratio::new(-1,4),Ratio::new(1,2)]]);

        let a1 = TwoPoint::new(Ratio::new(0,1),Ratio::new(0,1));
        let a2 = TwoPoint::new(Ratio::new(2,1),Ratio::new(2,1));
        let a3 = TwoPoint::new(Ratio::new(-3,1),Ratio::new(1,1));

        let v12 = TwoVector {start: a1.clone(), end: a2.clone()};
        let v23 = TwoVector {start: a2, end: a3.clone()};
        let v31 = TwoVector {start: a3, end: a1};

        let q1 = v23.quadrance_metric(&m);
        let q2 = v31.quadrance_metric(&m);
        let q3 = v12.quadrance_metric(&m);

        assert_eq!(q1,Ratio::new(59,8));
        assert_eq!(q2,Ratio::new(43,8));
        assert_eq!(q3,Ratio::new(3,2));

        let s1 = v12.spread_metric(&v31,&m);
        let s2 = v12.spread_metric(&v23,&m);
        let s3 = v23.spread_metric(&v31,&m);

        assert_eq!(s1,Ratio::new(128,129));
        assert_eq!(s2,Ratio::new(128,177));
        assert_eq!(s3,Ratio::new(512,2537));

        assert_eq!(s1/q1,s2/q2);
        assert_eq!(s2/q2,s3/q3);
        assert_eq!(s3/q3,Ratio::new(1024,7611));
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
    fn composition_of_rotations() {
        let ra = Rotation { vector: TwoVector::new(Ratio::new(15,17),Ratio::new(8,17)) };
        let rb = Rotation { vector: TwoVector::new(Ratio::new(-7,25),Ratio::new(24,25)) };
        let rc = Rotation { vector: TwoVector::new(Ratio::new(-297,425),Ratio::new(304,425)) };
        assert_eq!(ra*rb,rc);
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
    #[test]
    fn length_between_one_points() {
        let a = OnePoint::new(Ratio::new(3,2));
        let b = OnePoint::new(Ratio::new(1,5));
        assert_eq!(a.length(&b),Ratio::new(13,10));
    }
    #[test]
    fn quadrance_between_one_points() {
        let a1 = OnePoint::new(2);
        let a2 = OnePoint::new(5);
        let a3 = OnePoint::new(7);
        assert_eq!(a2.quadrance(&a3),4);
        assert_eq!(a1.quadrance(&a3),25);
        assert_eq!(a1.quadrance(&a2),9);
    }
    #[test]
    fn bretschneider_von_staudt() {
        let p1 = TwoPoint::new(Ratio::new(1,1), Ratio::new(1,2));
        let p2 = TwoPoint::new(Ratio::new(2,1), Ratio::new(3,2));
        let p3 = TwoPoint::new(Ratio::new(-1,1), Ratio::new(4,3));
        let p4 = TwoPoint::new(Ratio::new(-3,1), Ratio::new(-4,3));

        let q = Quadrilateral::new(p1.clone(),p2.clone(),p3.clone(),p4.clone());
        let p = Polygon {points: vec![p1, p2, p3, p4]};
        let a = p.area();
        assert_eq!(Ratio::new(16,1)*a*a,q.quadrea());
    }
    #[test]
    fn unit_circle_in_red_geometry() {
        let center = TwoPoint::new(Ratio::new(0,1), Ratio::new(0,1));
        let quadrance = Ratio::new(1,1);
        let a1 = TwoPoint::new(Ratio::new(1,1), Ratio::new(0,1));
        let a2 = TwoPoint::new(Ratio::new(13,12), Ratio::new(5,12));
        let a3 = TwoPoint::new(Ratio::new(-5,3), Ratio::new(4,3));
        let a4 = TwoPoint::new(Ratio::new(-17,15), Ratio::new(-8,15));
        let circle = TwoCircle::new_red(center, quadrance);
        assert!(circle.lies_on_red(&a1));
        assert!(circle.lies_on_red(&a2));
        assert!(circle.lies_on_red(&a3));
        assert!(circle.lies_on_red(&a4));
        assert_eq!(a1.quadrance_red(&a2), Ratio::new(-1,6));
        assert_eq!(a2.quadrance_red(&a3), Ratio::new(121,18));
        assert_eq!(a3.quadrance_red(&a4), Ratio::new(-16,5));
        assert_eq!(a1.quadrance_red(&a4), Ratio::new(64,15));
        assert_eq!(a1.quadrance_red(&a3), Ratio::new(16,3));
        assert_eq!(a2.quadrance_red(&a4), Ratio::new(361,90));
    }
    #[test]
    fn spreads_in_a_triangle() {
        let a1 = TwoPoint::new(Ratio::new(-1,1),Ratio::new(-2,1));
        let a2 = TwoPoint::new(Ratio::new(4,1),Ratio::new(1,1));
        let a3 = TwoPoint::new(Ratio::new(1,1),Ratio::new(3,1));

        let l1 = a2.join(&a3);
        let l2 = a1.join(&a3);
        let l3 = a1.join(&a2);

        let s1 = l2.spread(&l3);
        let s2 = l1.spread(&l3);
        let s3 = l1.spread(&l2);

        assert_eq!(s1,Ratio::new(361,986));
        assert_eq!(s2,Ratio::new(361,442));
        assert_eq!(s3,Ratio::new(361,377));
    }
    #[test]
    fn circumquadrance_of_a_triangle() {
        let a1 = TwoPoint::new(Ratio::new(1,1),Ratio::new(3,1));
        let a2 = TwoPoint::new(Ratio::new(2,1),Ratio::new(2,1));
        let a3 = TwoPoint::new(Ratio::new(-1,1),Ratio::new(-2,1));

        let t = TwoTriangle::new(a1,a2,a3);
        assert_eq!(t.quadrea(),Ratio::new(196,1));
        assert_eq!(t.circumquadrance(),Ratio::new(725,98));
    }
    #[test]
    fn curvature_of_a_circle() {
        let a = TwoPoint::new(Ratio::new(1,1),Ratio::new(2,1));
        let r = Ratio::new(5,1);
        let c = TwoCircle::new(a,r);
        assert_eq!(c.curvature(),Ratio::new(1,5));
    }
}
