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

use num::{Zero,One};
use std::ops::{Add,Sub,Mul};
use std::cmp::PartialOrd;
use polynumber::*;


/// Represents a formal power series in several variables of form
/// sum(f(n1,n2,n3,...)*x1^n1*x2^n2*x3^n3... ; n1,n2,3=0..infinity).
/// Coefficients are calculated by a function f taking the vector
/// of indices n1,n2,n3,... from the coefficients of a polynumber as an input.
/// The resulting sum is an infinite polynumber in as many variables
/// as the index vector, implemented as polynumber, has elements.
type FormalPowerFnRaw<T> = dyn Fn(PolyNumber<T>) -> T;
type FormalPowerFn<T> = Box<FormalPowerFnRaw<T>>;

pub struct FormalPowerSeries<T>
{
    f: FormalPowerFn<T>,
}


impl<T> FormalPowerSeries<T>
{
    pub fn new(f: FormalPowerFn<T>) -> Self {
        FormalPowerSeries{ f: f }
    }

    pub fn coef(&self, n: PolyNumber<T>) -> T {
        (self.f)(n)
    }
}


impl<T> Add for FormalPowerSeries<T>
where
    T: Add<Output = T>,
    T: 'static,
    T: Clone,
{
    type Output = FormalPowerSeries<T>;

    fn add(self, other: Self) -> FormalPowerSeries<T> {
        let sf = std::boxed::Box::leak(self.f);
        let of = std::boxed::Box::leak(other.f);
        let f = Box::new(|n: PolyNumber<T>| (sf)(n.clone()) + (of)(n));
        FormalPowerSeries::new(f)
    }
}

impl<T> Mul for FormalPowerSeries<T>
where
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: PartialOrd,
    T: Zero,
    T: One,
    T: 'static,
    T: Clone,
{
    type Output = FormalPowerSeries<T>;

    fn mul(self, other: Self) -> FormalPowerSeries<T> {
        let sf = std::boxed::Box::leak(self.f);
        let of = std::boxed::Box::leak(other.f);

        let sum = |f: &dyn Fn(PolyNumber<T>)->T, k: PolyNumber<T>| {
                let mut a = T::zero();
                let mut l = k.clone() * T::zero();
                let alpha = create_polynumber_var!(alpha; alpha ; T);
                let mut alpha_to_the_i = create_polynumber_one!(alpha; T);
                for i in 0 .. k.order() {
                    while l.get(i) <= k.get(i) {
                        a = a + f(l.clone());
                        l = l + alpha_to_the_i.clone();
                    }
                    alpha_to_the_i = alpha_to_the_i.clone() * alpha.clone();
                }
                return a;
            };
        let lsum = std::boxed::Box::leak(Box::new(sum));
        let cauchy = Box::new(|k: PolyNumber<T>| lsum(&|l: PolyNumber<T>|(sf)(l.clone()) * (of)(k.clone()-l),k.clone()));
        FormalPowerSeries::new(cauchy)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coefficient_extraction_from_formal_power_series() {
        let leibniz = |n: PolyNumber<i32>| (2*n.get(0)+1) * if n.get(0)%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));
        let one = create_polynumber_one!(alpha; i32);

        assert_eq!(a.coef(one.clone()*0),1);
        assert_eq!(a.coef(one.clone()*1),-3);
        assert_eq!(a.coef(one.clone()*2),5);
        assert_eq!(a.coef(one.clone()*3),-7);
        assert_eq!(a.coef(one.clone()*4),9);
        assert_eq!(a.coef(one.clone()*5),-11);
    }
    #[test]
    fn adding_formal_power_series() {
        let leibniz = |n: PolyNumber<i32>| (2*n.get(0)+1) * if n.get(0)%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));
        let ifodd = |n: PolyNumber<i32>| if n.get(0)%2 > 0 {n.get(0)+1} else {0};
        let b = FormalPowerSeries::new(Box::new(ifodd));
        let c = a + b;
        let one = create_polynumber_one!(alpha; i32);

        assert_eq!(c.coef(one.clone()*0),1);
        assert_eq!(c.coef(one.clone()*1),-1);
        assert_eq!(c.coef(one.clone()*2),5);
        assert_eq!(c.coef(one.clone()*3),-3);
        assert_eq!(c.coef(one.clone()*4),9);
        assert_eq!(c.coef(one.clone()*5),-5);
    }
    #[test]
    fn multiplying_formal_power_series() {
        let leibniz = |n: PolyNumber<i32>| (2*n.get(0)+1) * if n.get(0)%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));
        let ifodd = |n: PolyNumber<i32>| if n.get(0)%2 > 0 {n.get(0)+1} else {0};
        let b = FormalPowerSeries::new(Box::new(ifodd));
        let c = a * b;
        let one = create_polynumber_one!(alpha; i32);

        assert_eq!(c.coef(one.clone()*0),0);
        assert_eq!(c.coef(one.clone()*1),2);
        assert_eq!(c.coef(one.clone()*2),-6);
        assert_eq!(c.coef(one.clone()*3),14);
        assert_eq!(c.coef(one.clone()*4),-26);
        assert_eq!(c.coef(one.clone()*5),44);
    }
}

