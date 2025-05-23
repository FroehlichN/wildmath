/*
Copyright 2024 Norbert Fröhlich


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

use num::{Integer,Zero,One};
use num::rational::{Ratio};
use std::ops::{Mul, Add, Sub, Neg, Div};
use crate::algebra::prime_factors;

/// Represents natural numbers extended with their nth roots
#[derive(Debug, Clone)]
pub struct Root<T,P>
where
    T: Clone,
    P: Integer,
    P: Copy,
{
    sum: Vec<RootProduct<T,P>>,
}

/// Represents a product of natural numbers extended with their nth roots
#[derive(Debug, Clone)]
struct RootProduct<T,P>
where
    P: Integer,
    P: Copy,
{
    factor: T,
    product: Vec<RootPrime<P>>,
}

/// Represents the nth root of a prime number
#[derive(Debug, Clone,PartialEq)]
struct RootPrime<P>
where
    P: Integer,
    P: Copy,
{
    exponent: Ratio::<u32>,
    radicand: P,
}

impl<T,P> Root<T,P>
where
    T: Clone,
    P: Integer,
    P: Copy,
{
    pub fn from(num: T) -> Root<T,P> {
        let p = Vec::new();
        let rp = RootProduct { factor: num, product: p };
        Root { sum: vec![rp] }
    }
}

impl<T,P> Root<T,P>
where
    T: Zero + One,
    T: Mul<P, Output = T>,
    T: Clone,
    P: Integer,
    P: Neg<Output = P>,
    P: Copy,
    Root<T,P>: Mul<Output = Root<T,P>>,
{
    pub fn root(degree: u32, radicand: P) -> Root<T,P> {
        if radicand.is_zero() {
            return Root::from(T::zero());
        }
        let pf = prime_factors(radicand);
        let mut product = Root::from(T::one());
        for (_, f) in pf.iter().enumerate() {
            let rprimef = RootPrime {exponent: Ratio::new(1, degree), radicand: *f};
            let rprodf = RootProduct {factor: T::one(), product: vec![rprimef]};
            let rootf = Root{sum: vec![rprodf]};
            product = product * rootf;
        }
        product
    }
}

impl<P> Root<Ratio<P>,P>
where
    P: Integer,
    P: Neg<Output = P>,
    P: Copy,
{
    pub fn root_of_ratio(degree: u32, radicand: Ratio<P>) -> Root<Ratio<P>,P> {
        let pf_numer = prime_factors(radicand.numer().clone());
        let pf_denom = prime_factors(radicand.denom().clone());
        let mut product = Root::from(Ratio::<P>::one()/radicand.denom());
        for (_, f) in pf_numer.iter().enumerate() {
            let rprimef = RootPrime {exponent: Ratio::new(1, degree), radicand: *f};
            let rprodf = RootProduct {factor: Ratio::<P>::one(), product: vec![rprimef]};
            let rootf = Root{sum: vec![rprodf]};
            product = product * rootf;
        }
        for (_, f) in pf_denom.iter().enumerate() {
            let rprimef = RootPrime {exponent: Ratio::new(degree - 1, degree), radicand: *f};
            let rprodf = RootProduct {factor: Ratio::<P>::one(), product: vec![rprimef]};
            let rootf = Root{sum: vec![rprodf]};
            product = product * rootf;
        }
        product
    }
}

impl<T,P> PartialEq for Root<T,P>
where
    T: PartialEq,
    T: Clone,
    P: Integer,
    P: Copy,
{
    fn eq(&self, other: &Self) -> bool {
        for (_, sval) in self.sum.iter().enumerate() {
            let mut value_found = false;
            for (_, oval) in other.sum.iter().enumerate() {
                if sval == oval {
                    value_found = true;
                    break;
                }
            }
            if !value_found {
                return false;
            }
        }

        for (_, oval) in other.sum.iter().enumerate() {
            let mut value_found = false;
            for (_, sval) in self.sum.iter().enumerate() {
                if sval == oval {
                    value_found = true;
                    break;
                }
            }
            if !value_found {
                return false;
            }
        }

        true
    }
}

impl<T,P> Mul<Root<T,P>> for Root<T,P>
where
    T: Zero,
    T: Clone,
    P: Integer,
    P: Copy,
    Root<T,P>: Mul<RootProduct<T,P>, Output = Root<T,P>>,
{
    type Output = Root<T,P>;

    fn mul(self, other: Self) -> Root<T,P> {
        let mut sum = Root::from(T::zero());
        for (_, s) in self.sum.iter().enumerate() {
            sum = sum.clone() + other.clone() * s.clone();
        }
        sum
    }
}

impl<T,P> Mul<RootProduct<T,P>> for Root<T,P>
where
    T: Zero,
    T: Clone,
    P: Integer,
    P: Copy,
    RootProduct<T,P>: Mul<Output = RootProduct<T,P>>,
{
    type Output = Root<T,P>;

    fn mul(self, other: RootProduct<T,P>) -> Root<T,P> {
        let mut sum = Root::from(T::zero());
        for (_, s) in self.sum.iter().enumerate() {
            let p = s.clone() * other.clone();
            let rp = Root { sum: vec![p] };
            sum = sum + rp;
        }
        sum
    }
}

impl<T,P> Mul<T> for Root<T,P>
where
    T: Zero,
    T: Clone,
    P: Integer,
    P: Copy,
    RootProduct<T,P>: Mul<Output = RootProduct<T,P>>,
{
    type Output = Root<T,P>;

    fn mul(self, other: T) -> Root<T,P> {
        let rp = RootProduct{ factor: other, product: vec![] };
        self * rp
    }
}

impl<T,P> PartialEq for RootProduct<T,P>
where
    T: PartialEq,
    T: Clone,
    P: Integer,
    P: Copy,
{
    fn eq(&self, other: &Self) -> bool {
        if self.factor != other.factor {
            return false;
        }

        for (_, sval) in self.product.iter().enumerate() {
            let mut value_found = false;
            for (_, oval) in other.product.iter().enumerate() {
                if sval == oval {
                    value_found = true;
                    break;
                }
            }
            if !value_found {
                return false;
            }
        }

        for (_, oval) in other.product.iter().enumerate() {
            let mut value_found = false;
            for (_, sval) in self.product.iter().enumerate() {
                if sval == oval {
                    value_found = true;
                    break;
                }
            }
            if !value_found {
                return false;
            }
        }

        true
    }
}


impl<T,P> Mul<RootProduct<T,P>> for RootProduct<T,P>
where
    T: One,
    T: Mul<P, Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
{
    type Output = RootProduct<T,P>;

    fn mul(self, other: RootProduct<T,P>) -> RootProduct<T,P> {
        let mut factor = self.factor.clone() * other.factor.clone();

        // search for common root primes and reduce, if possible
        let mut sproduct = self.product.clone();
        let mut oproduct = other.product.clone();
        let mut common_root_prime_found = false;
        for (srpi, srp) in self.product.iter().enumerate() {
            if common_root_prime_found {break;}

            for (orpi, orp) in other.product.iter().enumerate() {
                if common_root_prime_found {break;}

                if srp.radicand == orp.radicand {
                    common_root_prime_found = true;
                    sproduct.remove(srpi);
                    oproduct.remove(orpi);

                    let mut exponent = srp.exponent + orp.exponent;
                    if exponent >= Ratio::from(1) {
                        factor = factor * srp.radicand;
                        exponent = exponent - Ratio::from(1);
                    }

                    if exponent > Ratio::from(0) {
                        let rp = RootPrime { exponent: exponent, radicand: srp.radicand };
                        sproduct.push(rp);
                    }
                }
            }
        }

        if common_root_prime_found {
            let srpr = RootProduct {factor: factor, product: sproduct}; // self root product reduced
            let orpr = RootProduct {factor: T::one(), product: oproduct}; // other root product reduced
            // recursively reduce until no common factors are left
            return srpr * orpr;
        }

        // after reducing there are no common factors in sproduct and oproduct
        // combine both vectors of factors to calculate the resulting product
        for (_, orp) in oproduct.iter().enumerate() {
            sproduct.push(orp.clone());
        }

        RootProduct {factor: factor, product: sproduct}
    }
}

impl<T,P> Add<Root<T,P>> for Root<T,P>
where
    T: Zero,
    T: Add<Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
{
    type Output = Root<T,P>;

    fn add(self, other: Self) -> Root<T,P> {

        let mut ssum = self.sum.clone();
        let mut osum = other.sum.clone();
        let mut common_factor_found = false;

        for (sidx, sval) in self.sum.iter().enumerate() {
            if common_factor_found {break;}

            for (oidx, oval) in other.sum.iter().enumerate() {
                if common_factor_found {break;}

                if sval.eq_product(oval) {
                    common_factor_found = true;
                    ssum.remove(sidx);
                    osum.remove(oidx);
                    let factor = sval.factor.clone() + oval.factor.clone();
                    if !factor.is_zero() {
                        let rp = RootProduct {factor: factor, product: sval.product.clone()};
                        ssum.push(rp);
                    }
                }
            }
        }

        if common_factor_found {
            let ssumr = Root {sum: ssum};
            let osumr = Root {sum: osum};
            return ssumr + osumr;
        }

        for (_, oval) in osum.iter().enumerate() {
            ssum.push(oval.clone());
        }
        Root {sum: ssum}
    }
}

impl<T,P> Sub<Root<T,P>> for Root<T,P>
where
    T: Neg<Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
    Root<T,P>: Add<Output = Root<T,P>>,
{
    type Output = Root<T,P>;

    fn sub(self, other: Self) -> Root<T,P> {
        self + (-other.clone())
    }
}

impl<T,P> Neg for Root<T,P>
where
    T: Neg<Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
{
    type Output = Root<T,P>;

    fn neg(self) -> Root<T,P> {
        let mut sum = Vec::new();
        for (_, val) in self.sum.iter().enumerate() {
            sum.push(-val.clone());
        }
        Root {sum: sum}
    }
}

impl<T,P> Div for Root<T,P>
where
    T: One + Zero,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Div<P, Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
    Root<T,P>: Mul<Output = Root<T,P>>,
    Root<T,P>: Sub<Output = Root<T,P>>,
{
    type Output = Root<T,P>;

    fn div(self, other: Self) -> Root<T,P> {
        self * other.inverse()
    }
}

impl<T,P> Root<T,P>
where
    T: One + Zero,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Div<P, Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
    Root<T,P>: Mul<Output = Root<T,P>>,
    Root<T,P>: Sub<Output = Root<T,P>>,
{
    pub fn inverse(&self) -> Root<T,P> {
        let s0 = self.sum[0].clone();
        let rpone = RootProduct {factor: T::one(), product: vec![]};
        let mut r0 = Root {sum: vec![rpone / s0.clone()]};
        loop {
            let rone = Root::from(T::one());
            let t0 = self.clone() * r0.clone();
            let d0 = rone - t0;
            if d0.sum.len().is_zero() {
                return r0;
            }
            if !d0.has_roots() {
                let m = Root::from(T::one() / (T::one() - d0.sum[0].factor.clone()));
                return r0 * m;
            }
            r0 = r0.clone() + Root {sum: vec![d0.sum[0].clone() / s0.clone()]};
        }
    }
}

impl<T,P> Root<T,P>
where
    T: Zero,
    T: Clone,
    P: Integer,
    P: Copy,
{
    pub fn has_roots(&self) -> bool {
        for (_, s) in self.sum.iter().enumerate() {
            if !s.factor.is_zero() && s.product.len() > 0 {
                return true;
            }
        }
        false
    }
}

impl<T,P> RootProduct<T,P>
where
    T: Clone,
    P: Integer,
    P: Copy,
{
    pub fn eq_product(&self, other: &Self) -> bool {
        for (_, sval) in self.product.iter().enumerate() {
            let mut value_found = false;
            for (_, oval) in other.product.iter().enumerate() {
                if sval == oval {
                    value_found = true;
                    break;
                }
            }
            if !value_found {
                return false;
            }
        }

        for (_, oval) in other.product.iter().enumerate() {
            let mut value_found = false;
            for (_, sval) in self.product.iter().enumerate() {
                if sval == oval {
                    value_found = true;
                    break;
                }
            }
            if !value_found {
                return false;
            }
        }

        true
    }
}

impl<T,P> Neg for RootProduct<T,P>
where
    T: Neg<Output = T>,
    T: Clone,
    P: Integer,
    P: Copy,
{
    type Output = RootProduct<T,P>;

    fn neg(self) -> RootProduct<T,P> {
        RootProduct {factor: -self.factor, product: self.product}
    }
}

impl<T,P> Div for RootProduct<T,P>
where
    T: Div<Output = T>,
    T: Div<P, Output = T>,
    P: Integer,
    P: Copy,
{
    type Output = RootProduct<T,P>;

    fn div(self, other: Self) -> RootProduct<T,P> {
        let mut factor = self.factor / other.factor;
        let mut product = self.product.clone();
        for (_, p) in other.product.iter().enumerate() {
            let exp = Ratio::from(1) - p.exponent;
            let rprime = RootPrime {exponent: exp, radicand: p.radicand};
            factor = factor / p.radicand;
            product.push(rprime);
        }
        RootProduct {factor: factor, product: product}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_root_squared() {
        let sq_root_two = Root::root(2,2);
        assert_eq!(sq_root_two.clone() * sq_root_two, Root::from(2));

        let sq_root_three = Root::root(2,3);
        assert_eq!(sq_root_three.clone() * sq_root_three, Root::from(3));
    }

    #[test]
    fn cube_root_cubed() {
        let cb_root_two = Root::root(3,2);
        assert_eq!(cb_root_two.clone() * cb_root_two.clone() * cb_root_two, Root::from(2));

        let sq_root_three = Root::root(2,3);
        let cq_root_three = Root::root(3,3);
        assert_eq!(sq_root_three.clone() * cq_root_three.clone() * cq_root_three.clone() * cq_root_three,
                   Root::from(3) * sq_root_three);
    }

    #[test]
    fn root_of_rational_number() {
        let r1 = Ratio::new(4,3);
        let root1 = Root::root_of_ratio(2,r1);
        assert_eq!(root1.clone() * root1, Root::from(r1));

        let r2 = Ratio::new(4,3);
        let root2 = Root::root_of_ratio(3,r2);
        assert_eq!(root2.clone() * root2.clone() * root2, Root::from(r2));
    }

    #[test]
    fn invers_of_rational_number() {
        let r1 = Root::from(Ratio::from(1)) + Root::root(3,2);
        let r2 = Root::from(Ratio::new(1,3)) *
            (Root::from(Ratio::from(1)) - Root::root(3,2) + Root::root(3,2) * Root::root(3,2));
        assert_eq!(r1.inverse(),r2);
        assert_eq!(r1.inverse()*r1.clone(),Root::from(Ratio::from(1)));
    }

    #[test]
    fn multiplication_with_one() {
        let one = Root::from(1);
        let val = Root::from(5) + Root::root(2, 2) * 10;
        assert_eq!(one.clone() * val.clone(), val);
        assert_eq!(val.clone() * one, val);
    }

    #[test]
    fn commute_add_mul() {
        let val1 = Root::from(1234) + Root::root(2, 2) * 78;
        let val2 = Root::from(329587) + Root::root(2, 2) * 10294;
        assert_eq!(val1.clone() + val2.clone(), val2.clone() + val1.clone());
        assert_eq!(val1.clone() * val2.clone(), val2.clone() * val1.clone());
    }

    #[test]
    fn inverse() {
        let a = 198;
        let b = 897;
        let one = Ratio::new(1,1);
        let val1 = Root::from(one*a) + Root::root_of_ratio(2, one*2) * (one * b);
        let val2 = val1.inverse();
        assert_eq!(val1 * val2, Root::from(one));
    }

    #[test]
    fn square_root_two() {
        let rt = Root::from(0) + Root::root(2, 2);
        let two = Root::from(2) + Root::root(2, 0);
        assert_eq!(rt.clone() * rt, two);
    }

    #[test]
    fn addition_of_exponents() {
        let two = Root::from(Ratio::from(2));
        let sqrt2 = Root::root_of_ratio(2, Ratio::from(2));
        let root6th2 = Root::root_of_ratio(6, Ratio::from(2));
        let root3rd2 = Root::root_of_ratio(3, Ratio::from(2));
        let r1 = sqrt2.clone()*root3rd2.clone()*sqrt2.clone();
        let r2 = two.clone()*root6th2.clone()*root6th2.clone();
        assert_eq!(r1,r2);
    }
}

