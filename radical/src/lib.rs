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

use num::{Integer};
use std::ops::{Mul, Add};
use algebra::prime_factors;

/// Represents natural numbers extended with their nth roots
#[derive(Debug, Clone)]
pub struct Root<T> {
    sum: Vec<RootProduct<T>>,
}

/// Represents a product of natural numbers extended with their nth roots
#[derive(Debug, Clone)]
struct RootProduct<T> {
    factor: T,
    product: Vec<RootPrime<T>>,
}

/// Represents the nth root of a prime number
#[derive(Debug, Clone,PartialEq)]
struct RootPrime<T> {
    degree: u32,
    exponent: u32,
    radicand: T,
}

impl<T> Root<T>
{
    pub fn from_int(num: T) -> Root<T> {
        let p = Vec::new();
        let rp = RootProduct { factor: num, product: p };
        Root { sum: vec![rp] }
    }
}

impl<T> Root<T>
where
    T: Integer,
    T: Copy,
{
    pub fn root(degree: u32, radicand: T) -> Root<T> {
        let pf = prime_factors(radicand);
        let mut product = Root::from_int(T::one());
        for (_, f) in pf.iter().enumerate() {
            let rprimef = RootPrime {degree: degree, exponent: 1, radicand: *f};
            let rprodf = RootProduct {factor: T::one(), product: vec![rprimef]};
            let rootf = Root{sum: vec![rprodf]};
            product = product * rootf;
        }
        product
    }
}

impl<T> PartialEq for Root<T>
where
    T: Integer,
    T: Copy,
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

impl<T> Mul<Root<T>> for Root<T>
where
    T: Integer,
    T: Copy,
{
    type Output = Root<T>;

    fn mul(self, other: Self) -> Root<T> {
        let mut sum = Root::from_int(T::zero());
        for (_, s) in self.sum.iter().enumerate() {
            sum = sum.clone() + other.clone() * s.clone();
        }
        sum
    }
}

impl<T> Mul<RootProduct<T>> for Root<T>
where
    T: Integer,
    T: Copy,
{
    type Output = Root<T>;

    fn mul(self, other: RootProduct<T>) -> Root<T> {
        let mut sum = Root::from_int(T::zero());
        for (_, s) in self.sum.iter().enumerate() {
            let p = s.clone() * other.clone();
            let rp = Root { sum: vec![p] };
            sum = sum + rp;
        }
        sum
    }
}

impl<T> PartialEq for RootProduct<T>
where
    T: Integer,
    T: Copy,
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


impl<T> Mul<RootProduct<T>> for RootProduct<T>
where
    T: Integer,
    T: Copy,
{
    type Output = RootProduct<T>;

    fn mul(self, other: RootProduct<T>) -> RootProduct<T> {
        let mut factor = self.factor.clone() * other.factor.clone();

        // search for common root primes and reduce, if possible
        let mut sproduct = self.product.clone();
        let mut oproduct = other.product.clone();
        let mut common_root_prime_found = false;
        for (srpi, srp) in self.product.iter().enumerate() {
            if common_root_prime_found {break;}

            for (orpi, orp) in other.product.iter().enumerate() {
                if common_root_prime_found {break;}

                if srp.degree == orp.degree && srp.radicand == orp.radicand {
                    common_root_prime_found = true;
                    sproduct.remove(srpi);
                    oproduct.remove(orpi);

                    let mut exponent = srp.exponent + orp.exponent;
                    if exponent >= srp.degree {
                        factor = factor * srp.radicand;
                        exponent = exponent - srp.degree;
                    }

                    if exponent > 0 {
                        let rp = RootPrime { degree: srp.degree, exponent: exponent, radicand: srp.radicand };
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

impl<T> Add<Root<T>> for Root<T>
where
    T: Integer,
    T: Copy,
{
    type Output = Root<T>;

    fn add(self, other: Self) -> Root<T> {

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
                    let factor = sval.factor + oval.factor;
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

impl<T> RootProduct<T>
where
    T: Integer,
    T: Copy,
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_root_squared() {
        let sq_root_two = Root::root(2,2);
        assert_eq!(sq_root_two.clone() * sq_root_two, Root::from_int(2));

        let sq_root_three = Root::root(2,3);
        assert_eq!(sq_root_three.clone() * sq_root_three, Root::from_int(3));
    }

    #[test]
    fn cube_root_cubed() {
        let cb_root_two = Root::root(3,2);
        assert_eq!(cb_root_two.clone() * cb_root_two.clone() * cb_root_two, Root::from_int(2));

        let sq_root_three = Root::root(2,3);
        let cq_root_three = Root::root(3,3);
        assert_eq!(sq_root_three.clone() * cq_root_three.clone() * cq_root_three.clone() * cq_root_three,
                   Root::from_int(3) * sq_root_three);
    }

}

