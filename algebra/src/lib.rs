use num::{Num};
use std::ops::{Div, Mul, Add, Sub, Neg};

#[derive(Debug)]
enum LinSummand<T> {
    Number(T),
    Unknown(T),
}

/// Represents a linear equation a*x + b + c*x + ... = ... + w + y * x + z
pub struct LinEq<T> {
    lhs: Vec<LinSummand<T>>,
    rhs: Vec<LinSummand<T>>,
}

impl<T> LinEq<T>
where
    T: Num,
    T: Copy,
{
    fn x(&self) -> T {
       let mut m = T::zero();
       let mut n = T::zero();

       for (_, s) in self.lhs.iter().enumerate() {
            match s {
                LinSummand::Number(a) => n = n - *a,
                LinSummand::Unknown(b) => m = m + *b,
            }
        }
        for (_, s) in self.rhs.iter().enumerate() {
            match s {
                LinSummand::Number(a) => n = n + *a,
                LinSummand::Unknown(b) => m = m - *b,
            }
        }
        n/m
    }
}

/// Represents a linear term in x
#[derive(Debug)]
pub struct LinTermInX<T> {
    summands: Vec<LinSummand<T>>,
}

impl<T> Div<T> for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn div(self, other: T) -> LinTermInX<T> {
        let mut q: Vec<LinSummand<T>> = Vec::new();

        for (_, s) in self.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => q.push( LinSummand::Number(*a / other) ),
                LinSummand::Unknown(b) => q.push( LinSummand::Unknown(*b / other) ),
            }
        }

        LinTermInX {summands: q}
    }
}

impl<T> Mul<T> for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn mul(self, other: T) -> LinTermInX<T> {
        let mut p: Vec<LinSummand<T>> = Vec::new();

        for (_, s) in self.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => p.push( LinSummand::Number(*a * other) ),
                LinSummand::Unknown(b) => p.push( LinSummand::Unknown(*b * other) ),
            }
        }

        LinTermInX {summands: p}
    }
}

impl<T> Add<T> for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn add(self, other: T) -> LinTermInX<T> {
        let mut s: Vec<LinSummand<T>> = self.summands;
        s.push( LinSummand::Number(other) );
        LinTermInX {summands: s}
    }
}

impl<T> Sub<T> for LinTermInX<T>
where
    T: Num,
    T: Neg<Output = T>,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn sub(self, other: T) -> LinTermInX<T> {
        let mut d: Vec<LinSummand<T>> = self.summands;
        d.push( LinSummand::Number(-other) );
        LinTermInX {summands: d}
    }
}

impl<T> Sub<LinSummand<T>> for LinTermInX<T>
where
    T: Num,
    T: Neg<Output = T>,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn sub(self, other: LinSummand<T>) -> LinTermInX<T> {
        let mut d: Vec<LinSummand<T>> = self.summands;
        match other {
            LinSummand::Number(a) => d.push( LinSummand::Number(-a) ),
            LinSummand::Unknown(b) => d.push( LinSummand::Unknown(-b) ),
        }
        LinTermInX {summands: d}
    }
}

impl<T> PartialEq for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    fn eq(&self, other: &Self) -> bool {
        let mut m = T::zero();
        let mut n = T::zero();
        let mut o = T::zero();
        let mut p = T::zero();

        for (_, s) in self.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => m = m + *a,
                LinSummand::Unknown(b) => n = n + *b,
            }
        }

        for (_, s) in other.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => o = o + *a,
                LinSummand::Unknown(b) => p = p + *b,
            }
        }

        return (m == o) && (n == p);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};

    #[test]
    fn xplus2eq7() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(1), LinSummand::Number(2)],
            rhs: vec![LinSummand::Number(7)]};
        assert_eq!(eq1.x(), 5);
    }
    #[test]
    fn xplus10eq2() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(1), LinSummand::Number(10)],
            rhs: vec![LinSummand::Number(2)]};
        assert_eq!(eq1.x(), -8);
    }
    #[test]
    fn threexeq12() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(3)],
            rhs: vec![LinSummand::Number(12)]};
        assert_eq!(eq1.x(), 4);
    }
    #[test]
    fn fivexeq11() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(Ratio::new(5, 1))],
            rhs: vec![LinSummand::Number(Ratio::new(11, 1))]};
        assert_eq!(eq1.x(), Ratio::new(11, 5));
    }
    #[test]
    fn solve_lin_eq() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(3), LinSummand::Number(4)],
            rhs: vec![LinSummand::Number(11), LinSummand::Unknown(-2), LinSummand::Number(8)]};
        assert_eq!(eq1.x(), 3);
    }
    #[test]
    fn terms_in_x() {
        let lhs = LinTermInX {summands: vec![LinSummand::Unknown(Ratio::new(6, 1)),
            LinSummand::Number(Ratio::new(-1, 1))]} / Ratio::new(2, 1);
        let rhs = LinTermInX {summands: vec![LinSummand::Unknown(Ratio::new(3, 1)),
            LinSummand::Number(Ratio::new(5, 1))]} / Ratio::new(7, 1);
        let eq1 = LinEq {lhs: lhs.summands, rhs: rhs.summands};
        assert_eq!(eq1.x(), Ratio::new(17, 36));
    }
    #[test]
    fn arithmatic_with_terms_in_x() {
        let t1 = LinTermInX {summands: vec![LinSummand::Unknown(1),
            LinSummand::Number(3)]};
        let t2 = t1 * 2;
        let t3 = t2 + 10;
        let t4 = t3 - 4;
        let t5 = t4 / 2;
        let t6 = t5 - LinSummand::Unknown(1);
        assert_eq!(t6, LinTermInX {summands: vec![LinSummand::Number(6)]});
    }
}
