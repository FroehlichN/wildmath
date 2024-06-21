/*
Copyright 2023 - 2024 Norbert Fröhlich


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

use num::{Zero,One,Integer};
use std::ops::{Mul, Add, Sub, Neg, Div};
use num::rational::Ratio;
use polynumber::*;
use algebra::factorial;


type Sequence<T> = Vec<T>;


pub fn farey<T>(n : T) -> Sequence<Ratio<T>>
where
    T: Zero,
    T: One,
    T: Integer,
    T: PartialEq,
    T: Copy,
{
    if n == T::zero() {
        return Sequence::new();
    }

    let mut f = Sequence::from([Ratio::new(T::zero(),T::one()),
                                Ratio::new(T::one() ,T::one())]);

    let mut denom = T::one();

    while denom <= n {
        let mut numer = T::one();
        while numer < denom {
            let new = Ratio::new(numer,denom);
            for (i, n) in f.iter().enumerate() {
                if new == *n {
                    break;
                }
                if new < *n {
                    f.insert(i, new);
                    break;
                }
            }
            numer = numer + T::one();
        }
        denom = denom + T::one();
    }
    return f;
}

pub fn forward_diff<T>(s: Sequence<T>) -> Sequence<T>
where
    T: Sub<Output = T>,
    T: Clone,
{
    let mut ds = Sequence::new();
    if s.len() < 2 {
        return ds;
    }

    for i in 1..s.len() {
        ds.push(s[i].clone() - s[i-1].clone());
    }
    ds
}


pub fn polynumber<T>(s0: Sequence<T>) -> PolyNumber<T>
where
    T: PartialEq + PartialOrd + Zero + One + Mul + Add + Clone,
    T: Sub<Output = T>,
    T: Neg<Output = T>,
    T: Div<Output = T>,
    T: std::fmt::Display,
{
    let one = create_polynumber_one!(alpha ; T);
    if s0.len() == 0 {
        return one.clone() * T::zero();
    }

    let mut p = one.clone() * s0[0].clone();
    let mut s1 = s0.clone();
    let alpha = create_polynumber_var!(alpha; alpha ; T);
    let mut index : usize = 1;
    let mut k = T::one();

    loop {
        s1 = forward_diff(s1);
        if s1.len() == 0 {
            return p;
        }

        let denom = factorial(k.clone());
        p = p.clone() + alpha.lowering_factorial_power(index) * (s1[0].clone()/denom);

        index = index + 1;
        k = k.clone() + T::one();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn farey_sequence() {
        let f1 = farey(1);
        let s1 = Sequence::from([Ratio::new(0,1), Ratio::new(1,1)]);
        assert_eq!(f1,s1);
        let f2 = farey(2);
        let s2 = Sequence::from([Ratio::new(0,1),Ratio::new(1,2),Ratio::new(1,1)]);
        assert_eq!(f2,s2);
        let f3 = farey(3);
        let s3 = Sequence::from([Ratio::new(0,1),Ratio::new(1,3),Ratio::new(1,2),
                                 Ratio::new(2,3),Ratio::new(1,1)]);
        assert_eq!(f3,s3);
        let f4 = farey(4);
        let s4 = Sequence::from([Ratio::new(0,1),Ratio::new(1,4),Ratio::new(1,3),
                    Ratio::new(1,2),Ratio::new(2,3),Ratio::new(3,4),Ratio::new(1,1)]);
        assert_eq!(f4,s4);
        let f5 = farey(5);
        let s5 = Sequence::from([Ratio::new(0,1),Ratio::new(1,5),Ratio::new(1,4),
                    Ratio::new(1,3),Ratio::new(2,5),Ratio::new(1,2),Ratio::new(3,5),
                    Ratio::new(2,3),Ratio::new(3,4),Ratio::new(4,5),Ratio::new(1,1)]);
        assert_eq!(f5,s5);
        let f6 = farey(6);
        let s6 = Sequence::from([Ratio::new(0,1),Ratio::new(1,6),Ratio::new(1,5),
                    Ratio::new(1,4),Ratio::new(1,3),Ratio::new(2,5),Ratio::new(1,2),
                    Ratio::new(3,5),Ratio::new(2,3),Ratio::new(3,4),Ratio::new(4,5),
                    Ratio::new(5,6),Ratio::new(1,1)]);
        assert_eq!(f6,s6);
    }

    #[test]
    fn polynumber_from_sequence_of_square_pyramidal_numbers() {
        let s1 = Sequence::from([Ratio::from(0),Ratio::from(1),Ratio::from(5),
            Ratio::from(14),Ratio::from(30),Ratio::from(55),Ratio::from(91),
            Ratio::from(140)]);
        let p1 = polynumber(s1);
        let alpha = create_polynumber_var!(alpha; alpha ; Ratio::<i64>);
        let one = create_polynumber_one!(alpha ; Ratio::<i64>);
        let p2 = alpha.clone()*(alpha.clone()+one.clone())
            *(alpha.clone()*Ratio::from(2)+one)*Ratio::new(1,6);
        assert_eq!(p1,p2);
    }

    #[test]
    fn polynumber_from_sequence_of_cubes() {
        let s1 = Sequence::from([Ratio::from(0),Ratio::from(1),Ratio::from(9),
            Ratio::from(36),Ratio::from(100),Ratio::from(225),Ratio::from(441)]);
        let p1 = polynumber(s1);
        let alpha = create_polynumber_var!(alpha; alpha ; Ratio::<i64>);
        let one = create_polynumber_one!(alpha ; Ratio::<i64>);
        let p2 = alpha.clone()*(alpha.clone()+one.clone())*Ratio::new(1,2);
        let p3 = p2.clone() * p2;
        assert_eq!(p1,p3);
    }
}

