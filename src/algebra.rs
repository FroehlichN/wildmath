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

use num::{Num, Integer,One};
use std::ops::{Div, Mul, Add, Sub, Neg};
use std::fmt::Debug;


#[derive(Debug)]
pub enum LinSummand<T> {
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
    pub fn x(&self) -> Option<T> {
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
        if m.is_zero() {
            return None;
        } else {
            return Some(n/m);
        }
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

pub fn pentagonal_nr<T>(n: T) -> Option<T>
where
    T: Num,
    T: Copy,
{
    let two = T::one() + T::one();
    let three = T::one() + T::one() + T::one();

    if two.is_zero() {
        return None;
    } else {
        return Some(n*(three*n-T::one())/two);
    }
}

pub fn sum_of_divisors<T>(i: T) -> T
where
    T: Integer,
    T: Copy,
{
    let mut s = T::zero();
    let mut n = T::zero() - i;

    while n <= i {
        let p_n = pentagonal_nr(n);

        match p_n {
            None => return T::one(),
            Some(p) => {
                let d = i - p;
                let sd: T;

                if d.is_zero() {
                    sd = i;
                } else if d > T::zero() && d < i {
                    sd = sum_of_divisors(d);
                } else {
                    sd = T::zero();
                }

                if n.is_odd() {
                    s = s + sd;
                } else {
                    s = s - sd;
                }
            },
        }

        n = n + T::one();
    }

    return s;
}

pub fn nr_of_partitions<T>(i: T) -> T
where
    T: Integer,
    T: Copy,
{
    let mut s = T::zero();
    let mut n = T::zero() - i;

    while n <= i {
        let p_n = pentagonal_nr(n);

        match p_n {
            None => return T::one(),
            Some(p) => {
                let d = i - p;
                let sd: T;

                if d.is_zero() {
                    sd = T::one();
                } else if d > T::zero() && d < i {
                    sd = nr_of_partitions(d);
                } else {
                    sd = T::zero();
                }

                if n.is_odd() {
                    s = s + sd;
                } else {
                    s = s - sd;
                }
            },
        }

        n = n + T::one();
    }

    return s;
}

pub fn factorial<T>(n: T) -> T
where
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: One,
    T: PartialOrd,
    T: Clone,
{
    let mut f = T::one();
    let mut p = T::one();
    while f <= n {
        p = p.clone() * f.clone();
        f = f.clone() + T::one();
    }
    return p;
}

pub fn pascal_array<T>(m: T, k: T) -> T
where
    T: Integer,
    T: Copy,
{
    factorial(m+k)/(factorial(m)*factorial(k))
}

pub fn choose<T>(n: T, k: T) -> T
where
    T: Integer,
    T: Copy,
{
    factorial(n)/(factorial(k)*factorial(n-k))
}

pub fn catalan<T>(n: T) -> T
where
    T: Integer,
    T: Copy,
{
    let two = T::one() + T::one();
    choose(two*n, n)/(n + T::one())
}

pub fn hyper_catalan<T>(n: Vec<T>) -> T
where
    T: Integer,
    T: Copy,
{
    let two = T::one() + T::one();
    let mut f1 = two;
    let mut f2 = T::one();
    let mut n1 = T::zero();
    let mut n2 = T::one();
    for m in &n {
        n1 = n1 + f1 * *m;
        n2 = n2 + f2 * *m;
        f1 = f1 + T::one();
        f2 = f2 + T::one();
    }
    let mut hc = factorial(n1)/factorial(n2);
    for m in n {
        hc = hc/factorial(m);
    }
    hc
}

pub fn prime_factors<T>(number: T) -> Vec<T>
where
    T: Integer,
    T: Neg<Output = T>,
    T: Clone,
{
    let mut pf = Vec::new();

    if number.is_zero() {
        return pf;
    }

    let mut rest;
    let mut factor = T::one() + T::one();

    if number < T::zero() {
        pf.push(-T::one());
        rest = - number;
    } else {
        rest = number;
    }

    while rest.clone() != T::one() {
        loop {
            let n = rest.clone() / factor.clone();
            if n.clone() * factor.clone() == rest.clone() {
                pf.push(factor.clone());
                rest = n.clone();
            } else {
                break;
            }
        }
        factor = factor.clone() + T::one();
    }
    pf
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};

    #[test]
    fn xplus2eq7() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(1), LinSummand::Number(2)],
            rhs: vec![LinSummand::Number(7)]};
        assert_eq!(eq1.x(), Some(5));
    }
    #[test]
    fn xplus10eq2() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(1), LinSummand::Number(10)],
            rhs: vec![LinSummand::Number(2)]};
        assert_eq!(eq1.x(), Some(-8));
    }
    #[test]
    fn threexeq12() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(3)],
            rhs: vec![LinSummand::Number(12)]};
        assert_eq!(eq1.x(), Some(4));
    }
    #[test]
    fn fivexeq11() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(Ratio::new(5, 1))],
            rhs: vec![LinSummand::Number(Ratio::new(11, 1))]};
        assert_eq!( eq1.x(), Some(Ratio::new(11, 5)) );
    }
    #[test]
    fn solve_lin_eq() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(3), LinSummand::Number(4)],
            rhs: vec![LinSummand::Number(11), LinSummand::Unknown(-2), LinSummand::Number(8)]};
        assert_eq!(eq1.x(), Some(3));
    }
    #[test]
    fn terms_in_x() {
        let lhs = LinTermInX {summands: vec![LinSummand::Unknown(Ratio::new(6, 1)),
            LinSummand::Number(Ratio::new(-1, 1))]} / Ratio::new(2, 1);
        let rhs = LinTermInX {summands: vec![LinSummand::Unknown(Ratio::new(3, 1)),
            LinSummand::Number(Ratio::new(5, 1))]} / Ratio::new(7, 1);
        let eq1 = LinEq {lhs: lhs.summands, rhs: rhs.summands};
        assert_eq!( eq1.x(), Some(Ratio::new(17, 36)) );
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
    #[test]
    fn pentagonal_numbers() {
        assert_eq!(pentagonal_nr(5),Some(35));
        assert_eq!(pentagonal_nr(6),Some(51));
    }
    #[test]
    fn sums_of_divisors() {
        assert_eq!(sum_of_divisors(1),1);
        assert_eq!(sum_of_divisors(2),3);
        assert_eq!(sum_of_divisors(3),4);
        assert_eq!(sum_of_divisors(4),7);
        assert_eq!(sum_of_divisors(5),6);
        assert_eq!(sum_of_divisors(6),12);
        assert_eq!(sum_of_divisors(12),28);
        assert_eq!(sum_of_divisors(14),24);
    }
    #[test]
    fn test_nr_of_partitions() {
        assert_eq!(nr_of_partitions(1),1);
        assert_eq!(nr_of_partitions(2),2);
        assert_eq!(nr_of_partitions(12),77);
        assert_eq!(nr_of_partitions(13),101);
        assert_eq!(nr_of_partitions(14),135);
    }
    #[test]
    fn factorials() {
        assert_eq!(factorial(1),1);
        assert_eq!(factorial(2),2);
        assert_eq!(factorial(4),24);
        assert_eq!(factorial(5),120);
        assert_eq!(factorial(6),720);
    }
    #[test]
    fn pascal_number() {
        assert_eq!(pascal_array(0,0),1);
        assert_eq!(pascal_array(1,1),2);
        assert_eq!(pascal_array(3,4),35);
        assert_eq!(pascal_array(3,5),56);
        assert_eq!(pascal_array(3,6),84);
        assert_eq!(pascal_array(4,5),126);
        assert_eq!(pascal_array(4,6),210);
        assert_eq!(choose(3,1),3);
        assert_eq!(choose(4,2),6);
    }
    #[test]
    fn prime_factorization() {
        assert_eq!(prime_factors(6),vec![2,3]);
        assert_eq!(prime_factors(10+23),vec![3,11]);
        assert_eq!(prime_factors(10*10+23),vec![3,41]);
        assert_eq!(prime_factors(10*10*10+23),vec![3,11,31]);
        assert_eq!(prime_factors(10*10*10*10+23),vec![3,13,257]);
        assert_eq!(prime_factors(10*10*10*10*10+23),vec![3,7,11,433]);
        assert_eq!(prime_factors(10*10*10*10*10*10+23),vec![3,333341]);
        assert_eq!(prime_factors(10*10*10*10*10*10*10+23),vec![3,11,19,41,389]);
    }
    #[test]
    fn prime_factorization_of_negative_numbers() {
        assert_eq!(prime_factors(-6),vec![-1,2,3]);
    }
    #[test]
    fn catalan_numbers() {
        assert_eq!(catalan(0),1);
        assert_eq!(catalan(1),1);
        assert_eq!(catalan(2),2);
        assert_eq!(catalan(3),5);
        assert_eq!(catalan(4),14);
        assert_eq!(catalan(5),42);
        assert_eq!(catalan(6i64),132);
        assert_eq!(catalan(7i64),429);
        assert_eq!(catalan(8i64),1430);
        assert_eq!(catalan(9i64),4862);
        assert_eq!(catalan(10i64),16796);
        assert_eq!(catalan(11i128),58786);
    }
    #[test]
    fn hyper_catalan_numbers() {
        assert_eq!(hyper_catalan(vec![0]),1);
        assert_eq!(hyper_catalan(vec![1]),1);
        assert_eq!(hyper_catalan(vec![2]),2);
        assert_eq!(hyper_catalan(vec![3]),5);
        assert_eq!(hyper_catalan(vec![4]),14);
        assert_eq!(hyper_catalan(vec![5]),42);
        assert_eq!(hyper_catalan(vec![6i64]),132);
        assert_eq!(hyper_catalan(vec![7i64]),429);
        assert_eq!(hyper_catalan(vec![8i64]),1430);
        assert_eq!(hyper_catalan(vec![9i64]),4862);
        assert_eq!(hyper_catalan(vec![10i128]),16796);
        assert_eq!(hyper_catalan(vec![11i128]),58786);
        assert_eq!(hyper_catalan(vec![2,0,1]),28);
        assert_eq!(hyper_catalan(vec![7i128,5i128]),2283421140);
    }
}

