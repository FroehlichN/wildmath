use num::Zero;
use std::ops::{Mul, Add, Sub, Div};

pub trait Point {
    fn is_collinear(&self, a2: &Self, a3: &Self) -> bool;
}

/// Represents a 2D point
#[derive(Debug)]
pub struct TwoPoint<T> {
    x : T,
    y : T,
}

impl<T> TwoPoint<T>
where
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero,
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
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero,
    T: Clone,
{
    fn is_collinear(&self, a2: &Self, a3: &Self) -> bool {
        let l = self.join(&a2);
        a3.lies_on(&l)
    }
}

impl<T> PartialEq for TwoPoint<T>
where
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.quadrance(&other).is_zero()
    }
}

/// Represents a 2D line a*x+b*y+c=0
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
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T, c: T) -> TwoLine<T> {
        if a.is_zero() && b.is_zero() {
            panic!("Line has proportion of <0:0:c>");
        }
        TwoLine::new_raw(a, b, c)
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
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Div<T, Output = T>,
    T: Zero,
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
}
