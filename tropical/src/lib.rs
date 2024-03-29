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


use std::ops::{Add, Mul, Div, Sub};
use core::cmp::max;


/// Represents a number in max-plus algebra
#[derive(PartialEq, Debug, Clone)]
pub struct MaxPlus<T> {
    pub n : T,
}

impl<T> MaxPlus<T>
{
    pub fn new(n : T) -> MaxPlus<T> {
        MaxPlus{ n: n }
    }
}


impl<T> Add<MaxPlus<T>> for MaxPlus<T>
where
    T: Ord,
{
    type Output = MaxPlus<T>;

    fn add(self, other: MaxPlus<T>) -> MaxPlus<T> {
        MaxPlus {n: max(self.n, other.n) }
    }
}

impl<T> Mul<MaxPlus<T>> for MaxPlus<T>
where
    T: Add<T, Output=T>,
{
    type Output = MaxPlus<T>;

    fn mul(self, other: MaxPlus<T>) -> MaxPlus<T> {
        MaxPlus {n: (self.n + other.n) }
    }
}

impl<T> MaxPlus<T>
where
    T: Mul<T, Output=T>,
{
    pub fn pow(self, exp: Self) -> MaxPlus<T> {
        MaxPlus {n: (self.n * exp.n) }
    }
}

impl<T> Div<MaxPlus<T>> for MaxPlus<T>
where
    T: Sub<T, Output=T>,
{
    type Output = MaxPlus<T>;

    fn div(self, other: MaxPlus<T>) -> MaxPlus<T> {
        MaxPlus {n: (self.n - other.n) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn max_plus_addition() {
        let t1 = MaxPlus::new(5);
        let t2 = MaxPlus::new(7);
        assert_eq!(t1+t2.clone(), t2);
    }
    #[test]
    fn max_plus_multiplication() {
        let t1 = MaxPlus::new(5);
        let t2 = MaxPlus::new(7);
        let t3 = MaxPlus::new(12);
        assert_eq!(t1*t2, t3);
    }
    #[test]
    fn max_plus_exponentiation() {
        let t1 = MaxPlus::new(3);
        let t2 = MaxPlus::new(4);
        let t3 = MaxPlus::new(12);
        assert_eq!(t1.pow(t2), t3);
    }
    #[test]
    fn max_plus_division() {
        let t1 = MaxPlus::new(12);
        let t2 = MaxPlus::new(7);
        let t3 = MaxPlus::new(5);
        assert_eq!(t1/t2, t3);
    }
}

