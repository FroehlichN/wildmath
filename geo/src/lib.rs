use num::Zero;

/// Represents a 2D point
pub struct TwoPoint<T> {
    x : T,
    y : T,
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

    #[test]
    #[should_panic]
    fn invalid_two_line() {
        TwoLine::new(0,0,1);
    }
}
