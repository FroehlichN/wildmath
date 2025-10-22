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


use num::{Num, Integer, Zero, One};
use std::ops::{Mul, Add, Sub, Div, Neg};
use std::fmt;



/// RatInf
/// extended rational numbers
/// rational numbers with infinity
#[derive(Clone,Copy)]
pub struct RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    a: T,
    b: T,
}

impl<T> RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    pub fn new(a: T, b: T) -> RatInf<T> {
        RatInf {a: a, b: b}
    }

    pub fn is_finite(&self) -> bool {
        !self.b.is_zero()
    }

    pub fn is_standard(&self) -> bool {
        !self.a.is_zero() || !self.b.is_zero()
    }

    pub fn is_infinite(&self) -> bool {
        !self.a.is_zero() && self.b.is_zero()
    }

    pub fn is_nil(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

}

impl<T> PartialEq for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        if self.is_nil() && other.is_nil() {
            return true;
        }
        let lhs = self.a.clone() * other.b.clone();
        let rhs = self.b.clone() * other.a.clone();
        let s_non_zero = !self.a.is_zero() || !self.b.is_zero();
        let o_non_zero = !other.a.is_zero() || !other.b.is_zero();

        return (lhs == rhs) && s_non_zero && o_non_zero ;
    }
}

impl<T> Add for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn add(self, other: Self) -> RatInf<T> {
        let na = self.a.clone() * other.b.clone()
                + self.b.clone() * other.a.clone();
        let nb = self.b.clone() * other.b.clone();
        return RatInf {a: na, b: nb };
    }
}

impl<T> Sub for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn sub(self, other: Self) -> RatInf<T> {
        let na = self.a.clone() * other.b.clone()
                - self.b.clone() * other.a.clone();
        let nb = self.b.clone() * other.b.clone();
        return RatInf {a: na, b: nb };
    }
}

impl<T> Mul for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn mul(self, other: Self) -> RatInf<T> {
        let na = self.a.clone() * other.a.clone();
        let nb = self.b.clone() * other.b.clone();
        return RatInf {a: na, b: nb };
    }
}

impl<T> Div for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn div(self, other: Self) -> RatInf<T> {
        let na = self.a.clone() * other.b.clone();
        let nb = self.b.clone() * other.a.clone();
        return RatInf {a: na, b: nb };
    }
}

impl<T> Zero for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    fn zero() -> RatInf<T> {
        RatInf {a: T::zero(), b: T::one()}
    }

    fn is_zero(&self) -> bool {
        self.a.is_zero()
    }
}

impl<T> One for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    fn one() -> RatInf<T> {
        RatInf {a: T::one(), b: T::one()}
    }
}

impl<T> Neg for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn neg(self) -> RatInf<T> {
        RatInf::<T>::zero() - self
    }
}

impl<T> Mul<T> for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn mul(self, other: T) -> RatInf<T> {
        self * RatInf::new(other,T::one())
    }
}

impl<T> Div<T> for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    type Output = RatInf<T>;

    fn div(self, other: T) -> RatInf<T> {
        self / RatInf::new(other,T::one())
    }
}

impl<T> fmt::Display for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
    T: std::fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({})/({})", self.a, self.b)
    }
}

impl<T> fmt::Debug for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
    T: std::fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RatInf")
         .field("numerator", &self.a)
         .field("denominator", &self.b)
         .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn equality_of_extended_rational_numbers() {
        let a1 = RatInf::new(2,3);
        let a2 = RatInf::new(4,6);
        assert_eq!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers2() {
        let a1 = RatInf::new(3,0);
        let a2 = RatInf::new(-2,0);
        assert_eq!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers3() {
        let a1 = RatInf::new(0,7);
        let a2 = RatInf::new(0,-24);
        assert_eq!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers4() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(0,0);
        assert_eq!(a1, a2);
    }
    #[test]
    fn inequality_of_extended_rational_numbers1() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(2,3);
        assert!(!(a1==a2));
    }
    #[test]
    fn inequality_of_extended_rational_numbers2() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(5,0);
        assert!(!(a1==a2));
    }
    #[test]
    fn inequality_of_extended_rational_numbers3() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(0,4);
        assert!(!(a1==a2));
    }
    #[test]
    fn finite_extednded_ration_numbers() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(0,4);
        assert!(!a1.is_finite());
        assert!(a2.is_finite());
    }
    #[test]
    fn standard_extednded_ration_numbers() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(0,4);
        let a3 = RatInf::new(1,4);
        let a4 = RatInf::new(1,0);
        assert!(!a1.is_standard());
        assert!(a2.is_standard());
        assert!(a3.is_standard());
        assert!(a4.is_standard());
    }
}

