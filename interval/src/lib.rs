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

use num::{Zero};
use std::ops::{Mul, Add, Sub, Div};
use std::cmp::{min, max, Ordering};


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
    T: Ord,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Interval<T>;

    fn mul(self, other : Self) -> Interval<T> {
        let p1 = self.l.clone() * other.l.clone();
        let p2 = self.l * other.u.clone();
        let p3 = self.u.clone() * other.l;
        let p4 = self.u * other.u;
        let pl = min(min(p1.clone(),p2.clone()),min(p3.clone(),p4.clone()));
        let pu = max(max(p1,p2),max(p3,p4));
        return Interval{ l: pl, u: pu};
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

impl<T> Interval<T>
where
    T: PartialOrd,
{
    fn is_between(&self, other : &T) -> bool {
        return self.l <= *other && *other <= self.u;
    }
}

impl<T> Sub for Interval<T>
where
    T: PartialOrd,
    T: Sub<Output = T>,
{
    type Output = Interval<T>;

    fn sub(self, other : Self) -> Interval<T> {
        let dl = self.l - other.u;
        let du = self.u - other.l;
        return Interval::new(dl, du);
    }
}

impl<T> PartialOrd for Interval<T>
where
    T: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.l > other.l && self.u < other.u {
            return Some(Ordering::Less);
        }
        if self.l < other.l && self.u > other.u {
            return Some(Ordering::Greater);
        }
        if self.l == other.l && self.u == other.u {
            return Some(Ordering::Equal);
        }
        return None;
    }

    // Provided methods
    fn le(&self, other: &Self) -> bool {
        return self.l >= other.l && self.u <= other.u;
    }

    fn ge(&self, other: &Self) -> bool {
        return self.l <= other.l && self.u >= other.u;
    }
}

impl<T> Div for Interval<T>
where
    T: Div<Output = T>,
    T: Zero,
    T: Clone,
    T: Ord,
{
    type Output = Interval<T>;

    fn div(self, other: Self) -> Interval<T> {
        if other.is_between(&T::zero()) {
            panic!("Division by Interval containing zero not allowed!");
        }
        let i1 = self.u.clone() / other.u.clone();
        let i2 = self.u / other.l.clone();
        let i3 = self.l.clone() / other.u;
        let i4 = self.l / other.l;
        let qu = max(max(i1.clone(),i2.clone()),max(i3.clone(),i4.clone()));
        let ql = min(min(i1,i2),min(i3,i4));
        return Interval::new(ql,qu);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
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
    #[test]
    fn is_between() {
        let i = Interval::new(-3,4);
        assert!(i.is_between(&-3));
        assert!(!i.is_between(&-4));
        assert!(i.is_between(&4));
        assert!(!i.is_between(&5));
    }
    #[test]
    fn multiplication_of_intervals3() {
        let i1 = Interval::new(-7,7);
        let i2 = Interval::new(8,10);
        let i3 = Interval::new(-70,70);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn multiplication_of_intervals4() {
        let i1 = Interval::new(-7,-5);
        let i2 = Interval::new(-10,8);
        let i3 = Interval::new(-56,70);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn subtraction_of_intervals() {
        let i1 = Interval::new(-1,3);
        let i2 = Interval::new(-2,4);
        let i3 = Interval::new(-5,5);
        assert_eq!(i1-i2,i3);
    }
    #[test]
    fn multiplication_of_intervals5() {
        let i1 = Interval::new(-1,2);
        let i2 = Interval::new(3,5);
        let i3 = Interval::new(-5,10);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn multiplication_of_intervals6() {
        let i1 = Interval::new(-4,-3);
        let i2 = Interval::new(1,2);
        let i3 = Interval::new(-8,-3);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn multiplication_of_intervals7() {
        let i1 = Interval::new(-5,1);
        let i2 = Interval::new(-2,3);
        let i3 = Interval::new(-15,10);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn multiplication_of_intervals8() {
        let i1 = Interval::new(-1,2);
        let i2 = Interval::new(-6,-4);
        let i3 = Interval::new(-12,6);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn multiplication_of_intervals9() {
        let i1 = Interval::new(-7,-5);
        let i2 = Interval::new(-1,0);
        let i3 = Interval::new(0,7);
        assert_eq!(i1*i2,i3);
    }
    #[test]
    fn interval_contains_interval() {
        let i1 = Interval::new(-6,12);
        let i2 = Interval::new(-11,13);
        assert!(i1<=i2);
    }
    #[test]
    fn addition_of_rational_intervals() {
        let one = Ratio::new(1,1);
        let i1 = Interval::new(one/2,one*3/4);
        let i2 = Interval::new(-one/3,one/4);
        let i3 = Interval::new(one/6,one);
        assert!(i1+i2==i3);
    }
    #[test]
    fn multiplication_of_rational_intervals() {
        let one = Ratio::new(1,1);
        let i1 = Interval::new(one/2,one*3/4);
        let i2 = Interval::new(-one/3,one/4);
        let i3 = Interval::new(-one/4,one*3/16);
        assert!(i1*i2==i3);
    }
    #[test]
    #[should_panic]
    fn division_of_rational_intervals_by_zero() {
        let one = Ratio::new(1,1);
        let i1 = Interval::new(one/2,one*3/4);
        let i2 = Interval::new(-one/3,one/4);
        let i3 = i1/i2;
    }
    #[test]
    fn division_of_rational_intervals() {
        let one = Ratio::new(1,1);
        let i1 = Interval::new(-one/3,one/4);
        let i2 = Interval::new(one/2,one*3/4);
        let i3 = Interval::new(-one*2/3,one/2);
        assert!(i1/i2==i3);
    }
}

