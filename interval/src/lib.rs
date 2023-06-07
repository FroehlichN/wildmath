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


use std::ops::{Mul, Add};


/// Represents the interval between two numbers
#[derive(Debug,Clone)]
pub struct Interval<T> {
    l: T,
    u: T,
}

impl<T> Interval<T>
where
    T: PartialOrd,
{
    fn new(m : T, n : T) -> Interval<T> {
        if m < n {
            return Interval{ l: m, u: n };
        } else {
            return Interval{ l: n, u: m };
        }
    }
}

impl<T> Add for Interval<T>
where
    T: PartialOrd,
    T: Add<Output = T>,
{
    type Output = Interval<T>;

    fn add(self, other : Self) -> Interval<T> {
        let sl = self.l + other.l;
        let su = self.u + other.u;
        return Interval::new(sl, su);
    }
}

impl<T> Mul for Interval<T>
where
    T: PartialOrd,
    T: Mul<Output = T>,
{
    type Output = Interval<T>;

    fn mul(self, other : Self) -> Interval<T> {
        let pl = self.l * other.l;
        let pu = self.u * other.u;
        return Interval::new(pl, pu);
    }
}

impl<T> PartialEq for Interval<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        return self.l == other.l && self.u == other.u;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn new_interval_is_ordered() {
        let m = 10;
        let n = 9;
        let i = Interval::new(m, n); 
        assert!(i.l < i.u);
    }
    #[test]
    fn addition_of_intervals() {
        let i1 = Interval::new(1,2);
        let i2 = Interval::new(4,7);
        let i3 = Interval::new(5,9);
        assert_eq!(i1+i2,i3);
    }
    #[test]
    fn addition_of_intervals2() {
        let i1 = Interval::new(7,7);
        let i2 = Interval::new(8,8);
        let i3 = Interval::new(15,15);
        assert_eq!(i1+i2,i3);
    }
    #[test]
    fn multiplication_of_intervals() {
        let i1 = Interval::new(1,2);
        let i2 = Interval::new(4,7);
        let i3 = Interval::new(4,14);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn multiplication_of_intervals2() {
        let i1 = Interval::new(7,7);
        let i2 = Interval::new(8,8);
        let i3 = Interval::new(56,56);
        assert_eq!(i1*i2,i3);
    }
}

