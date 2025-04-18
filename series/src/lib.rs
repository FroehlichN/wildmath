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

use std::ops::{Add};


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


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn adding_formal_power_series() {
        let a = FormalPowerSeries::new(Box::new(|n: i32| n*n));
        let b = FormalPowerSeries::new(Box::new(|n: i32| n*n*n));
        let _c = a + b;
    }
}

