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
use std::ops::{Add,Sub,Mul,Div,Neg};


/// Represents a red complex number
#[derive(Debug,Clone)]
pub struct ComplexRed<T> {
    t: T,
    x: T,
}

impl<T> PartialEq for ComplexRed<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t && self.x == other.x
    }
}

impl<T> Add for ComplexRed<T>
where
    T: Add<Output = T>,
{
    type Output = ComplexRed<T>;

    fn add(self, other: Self) -> Self {
        let sumt = self.t + other.t;
        let sumx = self.x + other.x;
        ComplexRed{ t: sumt, x: sumx }
    }
}

impl<T> Mul for ComplexRed<T>
where
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = ComplexRed<T>;

    fn mul(self, other: Self) -> Self {
        let mult = self.clone().t * other.clone().t + self.clone().x * other.clone().x;
        let mulx = self.clone().t * other.clone().x + self.x * other.t;
        ComplexRed{ t: mult, x: mulx }
    }
}

impl<T> Mul<T> for ComplexRed<T>
where
    T: Clone,
    T: Mul<Output = T>,
{
    type Output = ComplexRed<T>;

    fn mul(self, other: T) -> Self {
        let mult = self.clone().t * other.clone();
        let mulx = self.x * other;
        ComplexRed{ t: mult, x: mulx }
    }
}

impl<T> ComplexRed<T>
where
    T: Neg<Output = T>,
{
    pub fn complex_conjugate(self) -> Self {
        ComplexRed{ t: self.t, x: -self.x }
    }
}

impl<T> ComplexRed<T>
where
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    pub fn quadrance(self) -> T {
        self.clone().t*self.clone().t - self.clone().x*self.x
    }
}

impl<T> ComplexRed<T>
where
    T: One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    ComplexRed<T>: Mul<T, Output = ComplexRed<T>>,
    T: Clone,
{
    pub fn inverse(self) -> Self {
        self.clone().complex_conjugate()*(T::one()/self.quadrance())
    }
}

impl<T> ComplexRed<T>
where
    T: One,
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new_from_half_slope(h: T) -> Self {
        let h2 = h.clone()*h.clone();
        let two = T::one() + T::one();
        let t = (T::one()+h2.clone())/(T::one()-h2.clone());
        let x = two*h/(T::one()-h2);
        ComplexRed{ t: t, x: x }
    }
}

impl<T> ComplexRed<T>
where
    T: Div<Output = T>,
{
    pub fn velocity(self) -> T {
        self.x / self.t
    }
}

impl<T> ComplexRed<T>
{
    pub fn new(t: T, x: T) -> Self {
        ComplexRed{ t: t, x: x }
    }
}

pub fn lemmermeyer_product<T>(h1: T, h2: T) -> T
where
    T: One,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    (h1.clone() + h2.clone())/(T::one() + h1*h2)
}

/// Represents a blue complex number
#[derive(Debug,Clone)]
pub struct ComplexBlue<T> {
    t: T,
    x: T,
}

impl<T> PartialEq for ComplexBlue<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.t == other.t && self.x == other.x
    }
}

impl<T> Zero for ComplexBlue<T>
where
    T: Zero,
{
    fn zero() -> ComplexBlue<T> {
        ComplexBlue{ t: T::zero(), x: T::zero() }
    }

    fn is_zero(&self) -> bool {
        self.t.is_zero() && self.x.is_zero()
    }
}

impl<T> One for ComplexBlue<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Clone,
{
    fn one() -> ComplexBlue<T> {
        ComplexBlue{ t: T::one(), x: T::zero() }
    }
}

impl<T> Add for ComplexBlue<T>
where
    T: Add<Output = T>,
{
    type Output = ComplexBlue<T>;

    fn add(self, other: Self) -> Self {
        let sumt = self.t + other.t;
        let sumx = self.x + other.x;
        ComplexBlue{ t: sumt, x: sumx }
    }
}

impl<T> Sub for ComplexBlue<T>
where
    T: Sub<Output = T>,
{
    type Output = ComplexBlue<T>;

    fn sub(self, other: Self) -> Self {
        let dift = self.t - other.t;
        let difx = self.x - other.x;
        ComplexBlue{ t: dift, x: difx }
    }
}

impl<T> Mul for ComplexBlue<T>
where
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = ComplexBlue<T>;

    fn mul(self, other: Self) -> Self {
        let mult = self.clone().t * other.clone().t - self.clone().x * other.clone().x;
        let mulx = self.clone().t * other.clone().x + self.x * other.t;
        ComplexBlue{ t: mult, x: mulx }
    }
}

impl<T> Div for ComplexBlue<T>
where
    T: Neg<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    type Output = ComplexBlue<T>;

    fn div(self, other: Self) -> Self {
        let div = self * other.clone().complex_conjugate() / other.quadrance();
        return div;
    }
}

impl<T> Mul<T> for ComplexBlue<T>
where
    T: Clone,
    T: Mul<Output = T>,
{
    type Output = ComplexBlue<T>;

    fn mul(self, other: T) -> Self {
        let mult = self.clone().t * other.clone();
        let mulx = self.x * other;
        ComplexBlue{ t: mult, x: mulx }
    }
}

impl<T> Div<T> for ComplexBlue<T>
where
    T: Clone,
    T: Div<Output = T>,
{
    type Output = ComplexBlue<T>;

    fn div(self, other: T) -> Self {
        let divt = self.clone().t / other.clone();
        let divx = self.x / other;
        ComplexBlue{ t: divt, x: divx }
    }
}

impl<T> ComplexBlue<T>
where
    T: Neg<Output = T>,
{
    pub fn complex_conjugate(self) -> Self {
        ComplexBlue{ t: self.t, x: -self.x }
    }
}

impl<T> ComplexBlue<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Clone,
{
    pub fn quadrance(self) -> T {
        self.clone().t*self.clone().t + self.clone().x*self.x
    }
}

impl<T> ComplexBlue<T>
where
    T: One,
    T: Neg<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    ComplexBlue<T>: Mul<T, Output = ComplexBlue<T>>,
    T: Clone,
{
    pub fn inverse(self) -> Self {
        self.clone().complex_conjugate()*(T::one()/self.quadrance())
    }
}

impl<T> ComplexBlue<T>
where
    T: Div<Output = T>,
{
    pub fn velocity(self) -> T {
        self.x / self.t
    }
}

impl<T> ComplexBlue<T>
{
    pub fn new(t: T, x: T) -> Self {
        ComplexBlue{ t: t, x: x }
    }

    pub fn real(self) -> T {
        self.t
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
    #[test]
    fn rational_parametrization_of_unit_circle() {
        let h = Ratio::new(1,2);
        let z = ComplexRed::new_from_half_slope(h);
        let z_ = ComplexRed::new(Ratio::new(5,3),Ratio::new(4,3));
        let w = ComplexRed::new(Ratio::new(1,1),Ratio::new(1,2));
        let w2 = w.clone()*w.clone();
        let w2_ = ComplexRed::new(Ratio::new(5,4),Ratio::new(1,1));
        let q = w.clone().quadrance();
        let v = z.clone().velocity();

        assert_eq!(q,Ratio::new(3,4));
        assert_eq!(z,z_);
        assert_eq!(v,Ratio::new(4,5));
        assert_eq!(w2,w2_);
    }
    #[test]
    fn stereographic_projection() {
        let h1 = Ratio::new(1,2);
        let u1 = ComplexRed::new_from_half_slope(h1);
        let h2 = Ratio::new(-3,2);
        let u2 = ComplexRed::new_from_half_slope(h2);
        let u3 = u1*u2;
        let h3 = lemmermeyer_product(h1,h2);
        assert_eq!(u3,ComplexRed::new_from_half_slope(h3));
    }
    #[test]
    fn collision() {
        let b_a = ComplexRed::new(Ratio::new(1,1),Ratio::new(1,2));
        let a_b = b_a.clone().inverse();
        assert_eq!(a_b.clone(),ComplexRed::new(Ratio::new(4,3),Ratio::new(-2,3)));

        let p_1a = ComplexRed::new(Ratio::new(-1,1),Ratio::new(-1,3));
        let p_1b = p_1a.clone() * a_b.clone();
        assert_eq!(p_1b.clone(),ComplexRed::new(Ratio::new(-10,9),Ratio::new(2,9)));

        let p_2a = ComplexRed::new(Ratio::new(-1,1),Ratio::new(1,4));
        let p_2b = p_2a.clone() * a_b.clone();
        assert_eq!(p_2b.clone(),ComplexRed::new(Ratio::new(-3,2),Ratio::new(1,1)));

        let r_1a = ComplexRed::new(Ratio::new(1,1),Ratio::new(-7,12));
        let r_1b = r_1a.clone() * a_b.clone();
        assert_eq!(r_1b.clone(),ComplexRed::new(Ratio::new(31,18),Ratio::new(-13,9)));

        let r_2a = ComplexRed::new(Ratio::new(1,1),Ratio::new(0,1));
        let r_2b = r_2a.clone() * a_b.clone();
        assert_eq!(r_2b.clone(),ComplexRed::new(Ratio::new(4,3),Ratio::new(-2,3)));

        let v_1b = p_1b.clone().velocity();
        assert_eq!(v_1b,Ratio::new(-1,5));

        let v_2b = p_2b.clone().velocity();
        assert_eq!(v_2b,Ratio::new(-2,3));

        let u_1b = r_1b.clone().velocity();
        assert_eq!(u_1b,Ratio::new(-26,31));

        let u_2b = r_2b.clone().velocity();
        assert_eq!(u_2b,Ratio::new(-1,2));
    }
}

