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

use num::{Zero,One};
use std::ops::{Mul, Add, Sub, Div};


/// Represents a number whose square can be two
#[derive(Debug,Clone)]
pub struct RootTwo<T> {
    a: T,
    b: T,
}

impl<T> PartialEq for RootTwo<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        return self.a == other.a && self.b == other.b;
    }
}

impl<T> Zero for RootTwo<T>
where
    T: Zero,
{
    fn zero() -> RootTwo<T> {
        return RootTwo{ a: T::zero(), b: T::zero() };
    }

    fn is_zero(&self) -> bool {
        return self.a.is_zero() && self.b.is_zero();
    }
}

impl<T> One for RootTwo<T>
where
    T: Zero,
    T: One,
    T: Clone,
{
    fn one() -> RootTwo<T> {
        return RootTwo{ a: T::one(), b: T::zero() };
    }
}

impl<T> Add for RootTwo<T>
where
    T: Add<Output = T>,
{
    type Output = RootTwo<T>;

    fn add(self, other: Self) -> RootTwo<T> {
        let sa = self.a + other.a;
        let sb = self.b + other.b;
        return RootTwo{ a: sa, b: sb };
    }
}

impl<T> Sub for RootTwo<T>
where
    T: Sub<Output = T>,
{
    type Output = RootTwo<T>;

    fn sub(self, other: Self) -> RootTwo<T> {
        let da = self.a - other.a;
        let db = self.b - other.b;
        return RootTwo{ a: da, b: db };
    }
}

impl<T> Mul for RootTwo<T>
where
    T: One,
    T: Add<T, Output = T>,
    T: Mul,
    T: Clone,
{
    type Output = RootTwo<T>;

    fn mul(self, other: Self) -> RootTwo<T> {
        let two = T::one() + T::one();
        let pa = self.a.clone() * other.a.clone() + two * self.b.clone() * other.b.clone();
        let pb = self.a * other.b + other.a * self.b;
        return RootTwo{ a: pa, b: pb };
    }
}

impl<T> RootTwo<T>
where
    T: Zero,
    T: One,
    T: Mul,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    fn inverse(&self) -> RootTwo<T> {
        let a2 = self.a.clone() * self.a.clone();
        let b2 = self.b.clone() * self.b.clone();
        let two = T::one() + T::one();
        let denom = a2 - two * b2;
        let ia = self.a.clone() / denom.clone();
        let ib = (T::zero() - self.b.clone()) / denom;
        return RootTwo{ a: ia, b: ib };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
    #[test]
    fn multiplication_with_one() {
        let one = RootTwo::<i32>::one();
        let val = RootTwo{ a: 5, b: 10 };
        assert_eq!(one.clone() * val.clone(), val);
        assert_eq!(val.clone() * one, val);
    }
    #[test]
    fn commute_add_mul() {
        let val1 = RootTwo{ a: 1234, b: 78 };
        let val2 = RootTwo{ a: 329587, b: 10294 };
        assert_eq!(val1.clone() + val2.clone(), val2.clone() + val1.clone());
        assert_eq!(val1.clone() * val2.clone(), val2.clone() * val1.clone());
    }
    #[test]
    fn inverse() {
        let a = 198;
        let b = 897;
        let one = Ratio::new(1,1);
        let val1 = RootTwo{ a: one*a, b: one*b };
        let val2 = val1.inverse();
        assert_eq!(val1 * val2, RootTwo::one());
    }
    #[test]
    fn square_root_two() {
        let rt = RootTwo{ a: 0, b: 1};
        let two = RootTwo{ a: 2, b: 0};
        assert_eq!(rt.clone() * rt, two);
    }
}

