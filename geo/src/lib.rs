use num::Zero;
use std::ops::{Mul, Add};

/// Represents a 2D point
pub struct TwoPoint<T> {
    x : T,
    y : T,
}

impl<T> TwoPoint<T>
where
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Zero,
{
    pub fn lies_on(self, l: TwoLine<T>) -> bool {
        let r = l.a * self.x + l.b * self.y + l.c;
        r.is_zero()
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
        assert!(a.lies_on(l));
    }
}
