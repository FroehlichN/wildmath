use num::Zero;
use std::ops::{Mul, Add, Sub};

/// Represents a 2D point
pub struct TwoPoint<T> {
    x : T,
    y : T,
}

impl<T> TwoPoint<T>
where
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Zero,
    T: Clone,
{
    pub fn lies_on(self, l: &TwoLine<T>) -> bool {
        let r = l.a.clone() * self.x + l.b.clone() * self.y + l.c.clone();
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

impl<T: Clone + Zero> TwoLine<T> {
    pub fn new(a: T, b: T, c: T) -> TwoLine<T> {
        if a.is_zero() && b.is_zero() {
            panic!("Line has proportion of <0:0:c>");
        }
        TwoLine::new_raw(a, b, c)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::Ratio;

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
}
