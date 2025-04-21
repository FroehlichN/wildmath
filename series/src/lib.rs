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



/// Index set for formal power series
#[derive(Clone)]
struct IndexSet<T> {
    m: Vec<T>,
}

impl<T> IndexSet<T>
{
    pub fn new(v : Vec<T>) -> IndexSet<T> {
        return IndexSet { m: v };
    }
}

impl<T> Mul<T> for IndexSet<T>
where
    T: Mul<T, Output = T>,
    T: Clone,
{
    type Output = IndexSet<T>;

    fn mul(self, other: T) -> IndexSet<T> {
        let mut p: Vec<T> = Vec::new();

        for a in self.m {
            p.push(a * other.clone());
        }

        return IndexSet::new(p);
    }
}

impl<T> Sub for IndexSet<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = IndexSet<T>;

    fn sub(self, other: Self) -> IndexSet<T> {
        let mut index: usize = 0;
        let mut s: Vec<T> = Vec::new();

        loop {
            let a = self.m.get(index);
            let b = other.m.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push((*aa).clone() - (*bb).clone()),
                (Some(aa), None)     => s.push((*aa).clone()),
                (None, Some(bb))     => s.push(T::zero() - (*bb).clone()),
                (None, None)         => return IndexSet::new(s),
            }

            index += 1;
        }
    }
}



/// Represents a formal power series in several variables of form
/// sum(f(n1,n2,n3,...)*x1^n1*x2^n2*x3^n3... ; n1,n2,3=0..infinity).
/// Coefficients are calculated by a function f taking the vector
/// of indices n1,n2,n3,... as an input.
/// The resulting sum is an infinite polynumber in as many variables
/// as the index vector has elements.
type FormalPowerFnRaw<T> = dyn Fn(Vec<T>) -> T;
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

    pub fn coef(&self, n: Vec<T>) -> T {
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
        let f = Box::new(|n: Vec<T>| (sf)(n.clone()) + (of)(n));
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

        let sum = |f: &dyn Fn(Vec<T>)->T, k: Vec<T>| {
                let mut a = T::zero();
                let mut l = IndexSet::<T>::new(k.clone()) * T::zero();
                for i in 0 .. k.len() {
                    while l.m[i] <= k[i] {
                        a = a + f(l.m.clone());
                        l.m[i] = (l.m[i]).clone() + T::one();
                    }
                }
                return a;
            };
        let lsum = std::boxed::Box::leak(Box::new(sum));
        let cauchy = Box::new(|k: Vec<T>|
                lsum(&|l: Vec<T>| {
                    let isl = IndexSet::<T>::new(l.clone());
                    let isk = IndexSet::<T>::new(k.clone());
                    let isd = isk - isl;
                    (sf)(l.clone()) * (of)(isd.m)
                },k.clone())
            );
        FormalPowerSeries::new(cauchy)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coefficient_extraction_from_formal_power_series() {
        let leibniz = |n: Vec<i32>| (2*n[0]+1) * if n[0]%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));

        assert_eq!(a.coef(vec![0]),1);
        assert_eq!(a.coef(vec![1]),-3);
        assert_eq!(a.coef(vec![2]),5);
        assert_eq!(a.coef(vec![3]),-7);
        assert_eq!(a.coef(vec![4]),9);
        assert_eq!(a.coef(vec![5]),-11);
    }
    #[test]
    fn adding_formal_power_series() {
        let leibniz = |n: Vec<i32>| (2*n[0]+1) * if n[0]%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));
        let ifodd = |n: Vec<i32>| if n[0]%2 > 0 {n[0]+1} else {0};
        let b = FormalPowerSeries::new(Box::new(ifodd));
        let c = a + b;

        assert_eq!(c.coef(vec![0]),1);
        assert_eq!(c.coef(vec![1]),-1);
        assert_eq!(c.coef(vec![2]),5);
        assert_eq!(c.coef(vec![3]),-3);
        assert_eq!(c.coef(vec![4]),9);
        assert_eq!(c.coef(vec![5]),-5);
    }
    #[test]
    fn multiplying_formal_power_series() {
        let leibniz = |n: Vec<i32>| (2*n[0]+1) * if n[0]%2 > 0 {-1} else {1};
        let a = FormalPowerSeries::new(Box::new(leibniz));
        let ifodd = |n: Vec<i32>| if n[0]%2 > 0 {n[0]+1} else {0};
        let b = FormalPowerSeries::new(Box::new(ifodd));
        let c = a * b;

        assert_eq!(c.coef(vec![0]),0);
        assert_eq!(c.coef(vec![1]),2);
        assert_eq!(c.coef(vec![2]),-6);
        assert_eq!(c.coef(vec![3]),14);
        assert_eq!(c.coef(vec![4]),-26);
        assert_eq!(c.coef(vec![5]),44);
    }
}

