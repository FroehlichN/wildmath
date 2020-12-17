use num::{Num};

/// Represents the equation x + m = n
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xplus2eq7() {
        let eq1 = BabyEq {m: 2, n: 7};
        assert_eq!(eq1.x(), 5);
    }
    fn xplus10eq2() {
        let eq1 = BabyEq {m: 10, n: 2};
        assert_eq!(eq1.x(), -8);
    }
}
