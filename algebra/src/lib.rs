use num::{Num};

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

       for (index, s) in self.lhs.iter().enumerate() {
            match s {
                LinSummand::Number(a) => n = n - *a,
                LinSummand::Unknown(b) => m = m + *b,
            }
        }
        for (index, s) in self.rhs.iter().enumerate() {
            match s {
                LinSummand::Number(a) => n = n + *a,
                LinSummand::Unknown(b) => m = m - *b,
            }
        }
        n/m
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
}
