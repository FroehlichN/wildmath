use num::{Num};

/// Represents the equation x + m = n, x=?
pub struct BabyEq<T> {
    m: T,
    n: T,
}

impl<T> BabyEq<T>
where
    T: Num,
    T: Copy,
{
    fn x(&self) -> T {
        self.n - self.m
    }
}

/// Represents the equation m*x = n, x=?
pub struct RatEq<T> {
    m: T,
    n: T,
}

impl<T> RatEq<T>
where
    T: Num,
    T: Copy,
{
    fn x(&self) -> T {
        self.n/self.m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};

    #[test]
    fn xplus2eq7() {
        let eq1 = BabyEq {m: 2, n: 7};
        assert_eq!(eq1.x(), 5);
    }
    #[test]
    fn xplus10eq2() {
        let eq1 = BabyEq {m: 10, n: 2};
        assert_eq!(eq1.x(), -8);
    }
    #[test]
    fn threexeq12() {
        let eq1 = RatEq {m: 3, n: 12};
        assert_eq!(eq1.x(), 4);
    }
    #[test]
    fn fivexeq11() {
        let eq1 = RatEq {m: Ratio::new(5, 1), n: Ratio::new(11, 1) };
        assert_eq!(eq1.x(), Ratio::new(11, 5));
    }
}
