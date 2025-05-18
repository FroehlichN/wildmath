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
use crate::matrix::*;


/// Represents a complex number as 2x2 matrix
#[derive(Debug,Clone)]
pub struct Complex<T> {
    matrix: Matrix<T>,
}

impl<T> PartialEq for Complex<T>
where
    T: PartialEq,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        self.matrix == other.matrix
    }
}

impl<T> Zero for Complex<T>
where
    T: Zero,
    T: Clone,
{
    fn zero() -> Complex<T> {
        let m = Matrix::new(vec![vec![T::zero(),T::zero()],
                                 vec![T::zero(),T::zero()]]);
        Complex{ matrix: m }
    }

    fn is_zero(&self) -> bool {
        self.matrix.is_zero()
    }
}

impl<T> One for Complex<T>
where
    T: Zero + One,
    T: Clone,
{
    fn one() -> Complex<T> {
        let m = Matrix::new(vec![vec![T::one(),T::zero()],
                                 vec![T::zero(),T::one()]]);
        Complex{ matrix: m }
    }
}

impl<T> Add for Complex<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Clone,
{
    type Output = Complex<T>;

    fn add(self, other: Self) -> Self {
        let sum = self.matrix + other.matrix;
        Complex{ matrix: sum }
    }
}

impl<T> Sub for Complex<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Clone,
{
    type Output = Complex<T>;

    fn sub(self, other: Self) -> Self {
        let sum = self.matrix - other.matrix;
        Complex{ matrix: sum }
    }
}

impl<T> Mul for Complex<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Complex<T>;

    fn mul(self, other: Self) -> Self {
        let product = self.matrix * other.matrix;
        Complex{ matrix: product }
    }
}

impl<T> Mul<T> for Complex<T>
where
    T: Zero,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Complex<T>;

    fn mul(self, other: T) -> Self {
        let product = self.matrix * other.clone();
        Complex{ matrix: product }
    }
}

impl<T> Complex<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: PartialEq,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn complex_conjugate(self) -> Self {
        Complex{ matrix: self.matrix.adjugate() }
    }

    pub fn quadrance(self) -> T {
        let two = T::one() + T::one();
        let q = self.clone() * self.complex_conjugate();
        return q.matrix.trace()/two;
    }

    pub fn inverse(self) -> Self {
        self.clone().complex_conjugate()*(T::one()/self.quadrance())
    }
}

impl<T> Div for Complex<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: PartialEq,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    type Output = Complex<T>;

    fn div(self, other: Self) -> Self {
        let div = self * other.inverse();
        return div;
    }
}

impl<T> Complex<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{

    pub fn new_blue_param(h: T) -> Self {
        let h2 = h.clone()*h.clone();
        let two = T::one() + T::one();
        let re = (T::one()-h2.clone())/(T::one()+h2.clone());
        let im = two*h/(T::one()+h2);
        Complex::new_blue(re,im)
    }

    pub fn new_red_param(h: T) -> Self {
        let h2 = h.clone()*h.clone();
        let two = T::one() + T::one();
        let re = (T::one()+h2.clone())/(T::one()-h2.clone());
        let im = two*h/(T::one()-h2);
        Complex::new_red(re,im)
    }

    pub fn new_green_param(h: T) -> Self {
        let m = Matrix::new(vec![vec![h.clone(),T::zero()],
                                 vec![T::zero(),T::one()/h]]);
        Complex{ matrix: m }
    }

}

impl<T> Complex<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    fn matrix_one() -> Matrix<T> {
        Matrix::new(vec![vec![T::zero(), T::one()],
                         vec![T::one(), T::zero()]])
    }
   fn matrix_blue() -> Matrix<T> {
        Matrix::new(vec![vec![T::zero(), T::one()],
                    vec![-T::one(), T::zero()]])
    }
    fn matrix_red() -> Matrix<T> {
        Matrix::new(vec![vec![T::zero(), T::one()],
                         vec![T::one(), T::zero()]])
    }
    fn matrix_green() -> Matrix<T> {
        Matrix::new(vec![vec![T::one(), T::zero()],
                    vec![T::zero(), -T::zero()]])
    }

    pub fn new(re: T, im: T) -> Self {
        let m = Matrix::new(vec![vec![re.clone(),im.clone()],
                                 vec![-im,re]]);
        Complex{ matrix: m }
    }

    pub fn new_blue(re: T, im: T) -> Self {
        let m = Matrix::new(vec![vec![re.clone(),im.clone()],
                                 vec![-im,re]]);
        Complex{ matrix: m }
    }

    pub fn new_red(re: T, im: T) -> Self {
        let m = Matrix::new(vec![vec![re.clone(),im.clone()],
                                 vec![im,re]]);
        Complex{ matrix: m }
    }

    pub fn new_green(re: T, im: T) -> Self {
        let m = Matrix::new(vec![vec![re.clone() - im.clone(), T::zero()],
                                 vec![T::zero(), im+re]]);
        Complex{ matrix: m }
    }

    pub fn real(self) -> T {
        self.matrix.get(0,0)
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


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
    #[test]
    fn rational_parametrization_of_unit_circle() {
        let h = Ratio::new(1,2);
        let z = Complex::new_red_param(h);
        let z_ = Complex::new_red(Ratio::new(5,3),Ratio::new(4,3));
        let w = Complex::new_red(Ratio::new(1,1),Ratio::new(1,2));
        let w2 = w.clone()*w.clone();
        let w2_ = Complex::new_red(Ratio::new(5,4),Ratio::new(1,1));
        let q = w.clone().quadrance();

        assert_eq!(q,Ratio::new(3,4));
        assert_eq!(z,z_);
        assert_eq!(w2,w2_);
    }
    #[test]
    fn stereographic_projection() {
        let h1 = Ratio::new(1,2);
        let u1 = Complex::new_red_param(h1);
        let h2 = Ratio::new(-3,2);
        let u2 = Complex::new_red_param(h2);
        let u3 = u1*u2;
        let h3 = lemmermeyer_product(h1,h2);
        assert_eq!(u3,Complex::new_red_param(h3));
    }
    #[test]
    fn collision() {
        let b_a = Complex::new_red(Ratio::new(1,1),Ratio::new(1,2));
        let a_b = b_a.clone().inverse();
        assert_eq!(a_b.clone(),Complex::new_red(Ratio::new(4,3),Ratio::new(-2,3)));

        let p_1a = Complex::new_red(Ratio::new(-1,1),Ratio::new(-1,3));
        let p_1b = p_1a.clone() * a_b.clone();
        assert_eq!(p_1b.clone(),Complex::new_red(Ratio::new(-10,9),Ratio::new(2,9)));

        let p_2a = Complex::new_red(Ratio::new(-1,1),Ratio::new(1,4));
        let p_2b = p_2a.clone() * a_b.clone();
        assert_eq!(p_2b.clone(),Complex::new_red(Ratio::new(-3,2),Ratio::new(1,1)));

        let r_1a = Complex::new_red(Ratio::new(1,1),Ratio::new(-7,12));
        let r_1b = r_1a.clone() * a_b.clone();
        assert_eq!(r_1b.clone(),Complex::new_red(Ratio::new(31,18),Ratio::new(-13,9)));

        let r_2a = Complex::new_red(Ratio::new(1,1),Ratio::new(0,1));
        let r_2b = r_2a.clone() * a_b.clone();
        assert_eq!(r_2b.clone(),Complex::new_red(Ratio::new(4,3),Ratio::new(-2,3)));
    }
}

