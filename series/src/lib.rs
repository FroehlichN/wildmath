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


/// Represents a formal power series of form sum(f(n)*x^n,n=0..infinity).
/// Coefficients are calculated by a function f taking the index n as an input.
/// The resulting sum is an infinite polynumber.
type FormalPowerFn<T> = Box<dyn Fn(T) -> T>;

pub struct FormalPowerSeries<T>
{
    f: FormalPowerFn<T>,
}


impl<T> FormalPowerSeries<T>
{
    pub fn new(f: FormalPowerFn<T>) -> Self {
        FormalPowerSeries{ f: f }
    }

    pub fn coef(&self, n: T) -> T {
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
        let f = Box::new(|n: T| (sf)(n.clone()) + (of)(n));
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

        let sum = |f: &dyn Fn(T)->T, k: T| {
                let mut a = T::zero();
                let mut l = T::zero();
                while l <= k {
                    a = a + f(l.clone());
                    l = l + T::one();
                }
                return a;
            };
        let lsum = std::boxed::Box::leak(Box::new(sum));
        let cauchy = Box::new(|k: T| lsum(&|l: T|(sf)(l.clone()) * (of)(k.clone()-l),k.clone()));
        FormalPowerSeries::new(cauchy)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn adding_formal_power_series() {
        let a = FormalPowerSeries::new(Box::new(|n: i32| n*n));
        let b = FormalPowerSeries::new(Box::new(|n: i32| n*n*n));
        let _c = a + b;
    }
    #[test]
    fn multiplying_formal_power_series() {
        let a = FormalPowerSeries::new(Box::new(|n: i32| n));
        let b = FormalPowerSeries::new(Box::new(|n: i32| n*n));
        let _c = a * b;
    }
    #[test]
    fn coefficient_extraction_from_formal_power_series() {
        let leibniz = |n: i32| (2*n+1) * if n%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));

        assert_eq!(a.coef(0),1);
        assert_eq!(a.coef(1),-3);
        assert_eq!(a.coef(2),5);
        assert_eq!(a.coef(3),-7);
        assert_eq!(a.coef(4),9);
        assert_eq!(a.coef(5),-11);
    }
}

