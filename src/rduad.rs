/*
Copyright 2025 Norbert Fröhlich


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


use num::{Num, Integer, Zero, One, Signed, signum};
use std::ops::{Mul, Add, Sub, Div, Neg};
use std::fmt;
use crate::algebra::prime_factors;



/// RDuad
/// extended rational numbers
/// ray-rational numbers with +/-infinity
#[derive(Clone,Copy)]
pub struct RDuad<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    a: T,
    b: T,
}

impl<T> RDuad<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    pub fn new(a: T, b: T) -> RDuad<T> {
        RDuad {a: a, b: b}
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

impl<T> RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    fn reduce(&self) -> RDuad<T> {
        if self.a.clone().is_zero() {
            return self.clone();
        }

        if self.b.clone().is_zero() {
            return self.clone();
        }


        let mut af = prime_factors(self.a.clone());
        let mut bf = prime_factors(self.b.clone());

        let mut common_factor_found = false;
        for (ai, av) in af.clone().iter().enumerate() {
            if *av == -T::one() {
                continue;
            } else {
                for (bi, bv) in bf.clone().iter().enumerate() {
                    if *av == *bv {
                        af.remove(ai);
                        bf.remove(bi);
                        common_factor_found = true;
                        break;
                    }
                }
            }
            if common_factor_found {
                break;
            }
        }

        let mut ar = T::one();
        for (_,v) in af.iter().enumerate() {
            ar = ar*v.clone();
        }

        let mut br = T::one();
        for (_,v) in bf.iter().enumerate() {
            br = br*v.clone();
        }

        let result = RDuad::new(ar,br);

        if common_factor_found {
            return result.reduce();
        } else {
            return result;
        }
    }
}

impl<T> PartialEq for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Signed,
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

        return (lhs == rhs) && s_non_zero && o_non_zero
               && (signum(self.a.clone()) == signum(other.a.clone()))
               && (signum(self.b.clone()) == signum(other.b.clone()));
    }
}

impl<T> Add for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn add(self, other: Self) -> RDuad<T> {
        let na = self.a.clone() * other.b.clone()
                + self.b.clone() * other.a.clone();
        let nb = self.b.clone() * other.b.clone();
        let n = RDuad {a: na, b: nb };
        return n.reduce();
    }
}

impl<T> Sub for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn sub(self, other: Self) -> RDuad<T> {
        let na = self.a.clone() * other.b.clone()
                - self.b.clone() * other.a.clone();
        let nb = self.b.clone() * other.b.clone();
        let n = RDuad {a: na, b: nb };
        return n.reduce();
    }
}

impl<T> Mul for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn mul(self, other: Self) -> RDuad<T> {
        let na = self.a.clone() * other.a.clone();
        let nb = self.b.clone() * other.b.clone();
        let n = RDuad {a: na, b: nb };
        return n.reduce();
    }
}

impl<T> Div for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn div(self, other: Self) -> RDuad<T> {
        let na = self.a.clone() * other.b.clone();
        let nb = self.b.clone() * other.a.clone();
        let n = RDuad {a: na, b: nb };
        return n.reduce();
    }
}

impl<T> Zero for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    fn zero() -> RDuad<T> {
        RDuad {a: T::zero(), b: T::one()}
    }

    fn is_zero(&self) -> bool {
        self.a.is_zero()
    }
}

impl<T> One for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    fn one() -> RDuad<T> {
        RDuad {a: T::one(), b: T::one()}
    }
}

impl<T> Neg for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn neg(self) -> RDuad<T> {
        RDuad::<T>::zero() - self
    }
}

impl<T> RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn subgation(self) -> RDuad<T> {
        RDuad {a: self.a, b: -(self.b)}
    }

    pub fn antigation(self) -> RDuad<T> {
        RDuad {a: -(self.a), b: -(self.b)}
    }
}

impl<T> Mul<T> for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn mul(self, other: T) -> RDuad<T> {
        self * RDuad::new(other,T::one())
    }
}

impl<T> Div<T> for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    type Output = RDuad<T>;

    fn div(self, other: T) -> RDuad<T> {
        self / RDuad::new(other,T::one())
    }
}

impl<T> fmt::Display for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Clone,
    T: std::fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({})//({})", self.a, self.b)
    }
}

impl<T> fmt::Debug for RDuad<T>
where
    T: Num,
    T: Integer,
    T: Clone,
    T: std::fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RDuad")
         .field("numerator", &self.a)
         .field("denominator", &self.b)
         .finish()
    }
}

impl<T> RDuad<T>
where
    T: Num,
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn turn(&self, other: &Self) -> RDuad<T> {
        let a = self.a.clone();
        let b = self.b.clone();
        let c = other.a.clone();
        let d = other.b.clone();
        RDuad::new(a.clone()*d.clone()-b.clone()*c.clone(), a*c+b*d)
    }
    pub fn coturn(&self, other: &Self) -> RDuad<T> {
        let a = self.a.clone();
        let b = self.b.clone();
        let c = other.a.clone();
        let d = other.b.clone();
        RDuad::new(a.clone()*c.clone()+b.clone()*d.clone(), a*d-b*c)
    }
    pub fn turn_red(&self, other: &Self) -> RDuad<T> {
        let a = self.a.clone();
        let b = self.b.clone();
        let c = other.a.clone();
        let d = other.b.clone();
        RDuad::new(a.clone()*d.clone()-b.clone()*c.clone(), a*c-b*d)
    }
    pub fn coturn_red(&self, other: &Self) -> RDuad<T> {
        let a = self.a.clone();
        let b = self.b.clone();
        let c = other.a.clone();
        let d = other.b.clone();
        RDuad::new(a.clone()*c.clone()-b.clone()*d.clone(), a*d-b*c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn equality_of_extended_rational_numbers() {
        let a1 = RDuad::new(2,3);
        let a2 = RDuad::new(4,6);
        assert_eq!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers2() {
        let a1 = RDuad::new(3,0);
        let a2 = RDuad::new(-2,0);
        assert_ne!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers3() {
        let a1 = RDuad::new(0,7);
        let a2 = RDuad::new(0,-24);
        assert_ne!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers4() {
        let a1 = RDuad::new(0,0);
        let a2 = RDuad::new(0,0);
        assert_eq!(a1, a2);
    }
    #[test]
    fn inequality_of_extended_rational_numbers1() {
        let a1 = RDuad::new(0,0);
        let a2 = RDuad::new(2,3);
        assert_ne!(a1, a2);
    }
    #[test]
    fn inequality_of_extended_rational_numbers2() {
        let a1 = RDuad::new(0,0);
        let a2 = RDuad::new(5,0);
        assert_ne!(a1, a2);
    }
    #[test]
    fn inequality_of_extended_rational_numbers3() {
        let a1 = RDuad::new(0,0);
        let a2 = RDuad::new(0,4);
        assert_ne!(a1, a2);
    }
    #[test]
    fn finite_extednded_ration_numbers() {
        let a1 = RDuad::new(0,0);
        let a2 = RDuad::new(0,4);
        assert!(!a1.is_finite());
        assert!(a2.is_finite());
    }
    #[test]
    fn standard_extednded_ration_numbers() {
        let a1 = RDuad::new(0,0);
        let a2 = RDuad::new(0,4);
        let a3 = RDuad::new(1,4);
        let a4 = RDuad::new(1,0);
        assert!(!a1.is_standard());
        assert!(a2.is_standard());
        assert!(a3.is_standard());
        assert!(a4.is_standard());
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_1() {
        let a1 = RDuad::new(2,3);
        let a2 = RDuad::new(4,5);
        let a3 = RDuad::new(22,15);
        assert_eq!(a1+a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_2() {
        let a1 = RDuad::new(-5,2);
        let a2 = RDuad::new(-3,-7);
        let a3 = RDuad::new(29,-14);
        assert_eq!(a1+a2,a3);
        let a4 = RDuad::new(58,-28);
        assert_eq!(a3,a4);
        let a5 = RDuad::new(-29,14);
        assert_ne!(a4,a5);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_3() {
        let a1 = RDuad::new(1,-3);
        let a2 = RDuad::new(-6,5);
        let a3 = RDuad::new(-6,-15);
        assert_eq!(a1*a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_4() {
        let a1 = RDuad::new(1,-3);
        let a2 = RDuad::new(-6,5);
        let a3 = RDuad::new(5,18);
        assert_eq!(a1/a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_5() {
        let a1 = RDuad::new(0,-1);
        let a2 = RDuad::new(2,8);
        let a3 = RDuad::new(-2,-8);
        assert_eq!(a1+a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_6() {
        let a1 = RDuad::new(0,-1);
        let a2 = RDuad::new(2,8);
        let a3 = RDuad::new(0,-8);
        assert_eq!(a1*a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_7() {
        let a1 = RDuad::new(1,-0);
        let a2 = RDuad::new(-2,3);
        let a3 = RDuad::new(3,0);
        assert_eq!(a1+a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_8() {
        let a1 = RDuad::new(1,-0);
        let a2 = RDuad::new(-2,3);
        let a3 = RDuad::new(-2,0);
        assert_eq!(a1*a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_9() {
        let a1 = RDuad::new(0,2);
        let a2 = RDuad::new(-3,0);
        let a3 = RDuad::new(-6,0);
        assert_eq!(a1+a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_10() {
        let a1 = RDuad::new(-5,-2);
        let a2 = RDuad::new(4,-3);
        let a3 = RDuad::new(23,6);
        assert_eq!(a1-a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_11() {
        let a1 = RDuad::new(-5,-2);
        let a2 = RDuad::new(-4,3);
        let a3 = RDuad::new(-23,-6);
        assert_eq!(a1-a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_12() {
        let a1 = RDuad::new(5,2);
        let a2 = RDuad::new(4,-3);
        let a3 = RDuad::new(-23,-6);
        assert_eq!(a1-a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_13() {
        let a1 = RDuad::new(5,2);
        let a2 = RDuad::new(-4,3);
        let a3 = RDuad::new(23,6);
        assert_eq!(a1-a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_14() {
        let a1 = RDuad::new(1,0);
        let a2 = RDuad::new(-1,0);
        let a3 = RDuad::new(0,0);
        assert_eq!(a1+a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_15() {
        let a1 = RDuad::new(0,-1);
        let a2 = RDuad::new(-2,-0);
        let a3 = RDuad::new(2,0);
        assert_eq!(a1+a2,a3);
    }
    #[test]
    fn arithmetic_with_ray_rational_numbers_16() {
        let a1 = RDuad::new(2,0);
        let a2 = RDuad::new(-3,0);
        let a3 = RDuad::new(-6,0);
        assert_eq!(a1*a2,a3);
    }
    #[test]
    fn negation_of_ray_rational_numbers() {
        let a1 = RDuad::new(2,1);
        let a2 = RDuad::new(-2,1);
        assert_eq!(-a1,a2);
    }
    #[test]
    fn subgation_of_ray_rational_numbers() {
        let a1 = RDuad::new(4,3);
        let a2 = RDuad::new(4,-3);
        let a3 = RDuad::new(-4,3);
        assert_eq!(a1.subgation(),a2);
        assert_ne!(a2,a3);
    }
    #[test]
    fn antigation_of_ray_rational_numbers() {
        let a1 = RDuad::new(2,1);
        let a2 = RDuad::new(4,2);
        let a3 = RDuad::new(-2,-1);
        assert_eq!(a1,a2);
        assert_ne!(a1,a3);
        assert_eq!(a1.antigation(),a3);
    }
    #[test]
    fn triple_turn_formula() {
        let v1 = RDuad::new(3,4);
        let v2 = RDuad::new(7,2);
        let v3 = RDuad::new(2,-3);

        let ru1 = v2.turn(&v3);
        let ru2 = v3.turn(&v1);
        let ru3 = v1.turn(&v2);

        assert_eq!(ru1+ru2+ru3,ru1*ru2*ru3);
    }
    #[test]
    fn triple_coturn_formula() {
        let v1 = RDuad::new(3,4);
        let v2 = RDuad::new(7,2);
        let v3 = RDuad::new(2,-3);

        let ro1 = v2.coturn(&v3);
        let ro2 = v3.coturn(&v1);
        let ro3 = v1.coturn(&v2);

        assert_eq!(ro1*ro2+ro2*ro3+ro3*ro1,RDuad::new(1,1));
    }
    #[test]
    fn red_triple_turn_formula() {
        let v1 = RDuad::new(3,4);
        let v2 = RDuad::new(7,2);
        let v3 = RDuad::new(2,-3);

        let ru1 = v2.turn_red(&v3);
        let ru2 = v3.turn_red(&v1);
        let ru3 = v1.turn_red(&v2);

        assert_eq!(ru1+ru2+ru3,-ru1*ru2*ru3);
    }
    #[test]
    fn red_triple_coturn_formula() {
        let v1 = RDuad::new(3,4);
        let v2 = RDuad::new(7,2);
        let v3 = RDuad::new(2,-3);

        let ro1 = v2.coturn_red(&v3);
        let ro2 = v3.coturn_red(&v1);
        let ro3 = v1.coturn_red(&v2);

        assert_eq!(ro1*ro2+ro2*ro3+ro3*ro1,RDuad::new(-1,1));
    }
}

