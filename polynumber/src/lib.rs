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
use std::fmt;

/// Construct variables of a polynumber, commonly known as polynomial.
/// first argument is the wanted variable.
/// second argument is the list of all variables sepearated by commas ,.
/// ctype is the type of the coefficients. The type has to implement Zero and One.
/// E.g. for a polynomial in two variables like 1 + 3x - 5.7xy^2, you can
/// create the x with
///   let x = create_polynumber_var(x; x,y ; Ratio::<i64>)
/// and the y with
///   let y = create_polynumber_var(y; x,y ; Ratio::<i64>)
#[macro_export]
macro_rules! create_polynumber_var {
    ( $var: ident ; $first: ident ; $ctype: ty) => {
        $crate::PolyNumber::new_var(vec![ <$ctype>::zero(), <$ctype>::one() ], stringify!($first))
    };
    ( $var: ident ; $first: ident, $( $tail: ident ),+ ; $ctype: ty) => {
        {
            let mut n = Vec::new();
            if stringify!($var) == stringify!($first) {
                n.push( $crate::create_polynumber_zero!( $($tail),* ; $ctype ) );
                n.push( $crate::create_polynumber_one!( $($tail),* ; $ctype ) );
            } else {
                n.push( $crate::create_polynumber_var!( $var; $($tail),* ; $ctype ) );
            }
            $crate::PolyNumber::new_var(n, stringify!($first))
        }
    };
}

/// Construct a polynumber, commonly known as polynomial, equal one.
/// first argument is the list of all variables.
/// ctype is the type of the coefficients. The type has to implement Zero and One.
#[macro_export]
macro_rules! create_polynumber_one {
    ( $first: ident ; $ctype: ty) => {
        $crate::PolyNumber::new_var(vec![ <$ctype>::one() ], stringify!($first))
    };

    ( $first: ident, $( $tail: ident ),+ ; $ctype: ty) => {
        $crate::PolyNumber::new_var(vec![ $crate::create_polynumber_one!( $($tail),* ; $ctype ) ], stringify!($first))
    };

}

/// Construct a polynumber, commonly known as polynomial, equal zero.
/// first argument is the list of all variables.
/// ctype is the type of the coefficients. The type has to implement Zero and One.
#[macro_export]
macro_rules! create_polynumber_zero {
    ( $first: ident ; $ctype: ty) => {
        $crate::PolyNumber::new(vec![ <$ctype>::zero() ])
    };

    ( $first: ident , $( $tail: ident ),+ ; $ctype: ty) => {
        $crate::PolyNumber::new(vec![ $crate::create_polynumber_zero!( $($tail),* ; $ctype) ])
    };

}

#[derive(Clone)]
pub struct PolyNumber<T> {
    n: Vec<T>,
    var: String,
}

impl<T> PolyNumber<T>
{
    pub fn order(&self) -> usize {
        return self.n.len();
    }
}

impl<T> PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    pub fn lowest_order(&self) -> Option<usize> {

        for (di, dv) in self.n.iter().enumerate() {
            if *dv != T::zero() {
                return Some(di);
            }
        }
        return None;
    }
}

impl<T> PartialEq for PolyNumber<T>
where
    T: PartialEq,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let mut index: usize = 0;

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => if *aa != *bb       { return false; },
                (Some(aa), None)     => if *aa != T::zero() { return false; },
                (None, Some(bb))     => if *bb != T::zero() { return false; },
                (None, None)         => return true,
            }

            index += 1;
        }
    }
}

impl<T> Add for PolyNumber<T>
where
    T: Add<Output = T>,
    T: Clone,
{
    type Output = PolyNumber<T>;

    fn add(self, other: Self) -> PolyNumber<T> {
        let mut index: usize = 0;
        let mut s: Vec<T> = Vec::new();
        let var = if self.var == "" {other.var} else {self.var};

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push((*aa).clone() + (*bb).clone()),
                (Some(aa), None)     => s.push((*aa).clone()),
                (None, Some(bb))     => s.push((*bb).clone()),
                (None, None)         => return PolyNumber::new_var(s, &var),
            }

            index += 1;
        }
    }
}

impl<T> Sub for PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = PolyNumber<T>;

    fn sub(self, other: Self) -> PolyNumber<T> {
        let mut index: usize = 0;
        let mut s: Vec<T> = Vec::new();
        let var = if self.var == "" {other.var} else {self.var};

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push((*aa).clone() - (*bb).clone()),
                (Some(aa), None)     => s.push((*aa).clone()),
                (None, Some(bb))     => s.push(T::zero() - (*bb).clone()),
                (None, None)         => return PolyNumber::new_var(s , &var),
            }

            index += 1;
        }
    }
}

impl<T> Mul for PolyNumber<T>
where
    T: Mul<Output = T>,
    T: Add<Output = T>,
    T: Clone,
{
    type Output = PolyNumber<T>;

    fn mul(self, other: Self) -> PolyNumber<T> {
        let mut p: Vec<T> = Vec::new();
        let var = if self.var == "" {other.var} else {self.var};

        for (ai, av) in self.n.iter().enumerate() {
            for (bi, bv) in other.n.iter().enumerate() {
                let ci = ai+bi;
                let c = p.get(ci);

                match c {
                    Some(cv) => p[ci] = (*cv).clone() + (*av).clone() * (*bv).clone(),
                    None     => p.push((*av).clone() * (*bv).clone()),
                }
            }
        }

        return PolyNumber::new_var(p, &var);
    }
}

impl<T> Mul<T> for PolyNumber<T>
where
    T: Mul<T, Output = T>,
    T: Clone,
{
    type Output = PolyNumber<T>;

    fn mul(self, other: T) -> PolyNumber<T> {
        let mut p: Vec<T> = Vec::new();

        for a in self.n {
            p.push(a * other.clone());
        }

        return PolyNumber::new_var(p, &self.var);
    }
}

impl<T> Mul<T> for PolyNumber<PolyNumber<T>>
where
    T: Mul,
    PolyNumber<T>: Mul<PolyNumber<T>, Output = PolyNumber<T>>,
    T: Clone,
{
    type Output = PolyNumber<PolyNumber<T>>;

    fn mul(self, other: T) -> PolyNumber<PolyNumber<T>> {
        let mut p: Vec<PolyNumber<T>> = Vec::new();

        for a in self.n {
            p.push(a.clone() * PolyNumber::new_var(vec![other.clone()], &a.var));
        }

        return PolyNumber::new_var(p, &self.var);
    }
}

impl<T> PolyNumber<T>
where
    T: Clone,
{
    pub fn eval<O>(&self, c: O) -> O
    where
        O: One,
        O: Zero,
        O: Mul,
        O: Mul<T, Output = O>,
        O: Add,
        O: Clone,
    {
        let mut ck = O::one();
        let mut pc = O::zero();

        for a in &self.n {
            pc = pc + ck.clone() * (*a).clone();
            ck = ck.clone() * c.clone();
        }

        return pc;
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: Clone,
{
    pub fn eval2<O>(&self, c: O) -> PolyNumber<O>
    where
        O: One,
        O: Zero,
        O: Mul,
        O: Mul<T, Output = O>,
        O: Add,
        O: Clone,
    {
        let mut v : Vec<O> = Vec::new();

        for a in &self.n {
            v.push((*a).eval(c.clone()));
        }

        return PolyNumber::new_var(v, &self.var);
    }
}

impl<T> PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn ltrans(&self, c: T) -> PolyNumber<T> {
        let q = PolyNumber::new_var(vec![c,T::one()], &self.var);
        return self.eval(q);
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn ltrans2(&self, r: T, s: T) -> PolyNumber<PolyNumber<T>> {
        let a = PolyNumber::new_var(vec![r,T::one()], &self.var); // r + a
        let dvar = "d".to_string() + &self.var.to_string();
        let b = PolyNumber::new_var(vec![PolyNumber::new_var(vec![s], &dvar), // s
                                     PolyNumber::new_var(vec![T::one()], &dvar) ], &self.var); // b
        let s2 = self.eval2(a); // set a = a + r
        return s2.eval(b); // set b = b + s
    }
}

impl<T> Zero for PolyNumber<T>
where
    T: Zero,
    T: PartialEq,
    T: Clone,
{
    fn zero() -> PolyNumber<T> {
        return PolyNumber::new_var(vec![T::zero()], "");
    }

    fn is_zero(&self) -> bool {
        return *self == Self::zero();
    }
}

impl<T> One for PolyNumber<T>
where
    T: One,
    T: Add<Output = T>,
    T: Clone,
{
    fn one() -> PolyNumber<T> {
        return PolyNumber::new_var(vec![T::one()], "");
    }
}

impl<T> PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn truncate(self, k: usize) -> PolyNumber<T> {
        let mut v : Vec<T> = Vec::new();
        for n in 0..k+1 {
            let a = self.n.get(n);
            match a {
                Some(aa) => v.push((*aa).clone()),
                None     => (),
            }
        }
        return PolyNumber::new_var(v, &self.var);
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn truncate2(self, k: usize) -> PolyNumber<PolyNumber<T>> {
        let mut v : Vec<PolyNumber<T>> = Vec::new();
        for i in 0..k+1 {
            let a = self.n.get(i);
            match a {
                Some(aa) => v.push((*aa).clone().truncate(k-i)),
                None     => (),
            }
        }
        return PolyNumber::new_var(v, &self.var);
    }
}

impl<T> PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    fn tangent(self, k: usize, c: T) -> PolyNumber<T> {
        let p_alpha_plus_c = self.ltrans(c.clone());
        let trunc = p_alpha_plus_c.truncate(k);
        return trunc.ltrans(T::zero() - c);
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    pub fn tangent2(self, k: usize, r: T, s: T) -> PolyNumber<PolyNumber<T>> {
        let p_at_rs = self.ltrans2(r.clone(),s.clone());
        let trunc = p_at_rs.truncate2(k);
        return trunc.ltrans2(T::zero() - r, T::zero() - s);
    }
}

impl<T> PolyNumber<T>
where
    T: Zero,
    T: One,
    T: PartialEq,
    T: Clone,
{
    fn taylor(self) -> PolyNumber<PolyNumber<T>> {
        let dvar = "d".to_string() + &self.var.to_string();
        let p = PolyNumber::new_var(vec! [ PolyNumber::new_var(vec![T::zero(),T::one()], &dvar),
                                       PolyNumber::new_var(vec![T::one(),T::zero()], &dvar) ], &self.var );
        return self.eval(p);
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: Zero,
    T: One,
    T: PartialEq,
    T: Clone,
{
    fn taylor2(self) -> PolyNumber<PolyNumber<PolyNumber<PolyNumber<T>>>> {
        let apg = PolyNumber::new(vec![PolyNumber::new(vec![T::zero(),T::one()]), // a
                                      PolyNumber::new(vec![T::one(),T::zero()]) ]); // g
        let bpd = PolyNumber::new(vec![
                      PolyNumber::new(vec![ // d^0
                          PolyNumber::new(vec![ // b^0
                              PolyNumber::new(vec![T::zero()] ) ] ),
                          PolyNumber::new(vec![ // b^1
                              PolyNumber::new(vec![T::one()] ) ] ) ] ), // b
                      PolyNumber::new(vec![ // d^1
                          PolyNumber::new(vec![ // b^0
                              PolyNumber::new(vec![T::one()] ) ] ) ] ) ] ); // d

        let papg = self.eval2(apg);
        return papg.eval(bpd);
    }
}

impl<T> PolyNumber<T>
where
    T: Zero,
    T: One,
    T: PartialEq,
    T: Clone,
{
    fn derivative(self, grade: usize) -> PolyNumber<T> {
        let tp = self.taylor();
        let  s = tp.n.get(grade);
        match s {
            Some(v) => return (*v).clone(),
            None    => return Self::zero(),
        }
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: Zero,
    T: One,
    T: PartialEq,
    T: Clone,
{
    fn derivative2(self, i : usize, j : usize) -> PolyNumber<PolyNumber<T>> {
        let tp = self.clone().taylor2();

        let bpol : PolyNumber<PolyNumber<PolyNumber<T>>>;
        let s = tp.n.get(j);
        match s {
            Some(p) => bpol = (*p).clone(),
            None    => return Self::zero(),
        }

        let mut v : Vec<PolyNumber<T>> = Vec::new();
        for b in bpol.n {
            let go = b.n.get(i);
            match go {
                Some(g) => v.push((*g).clone()),
                None    => v.push(PolyNumber::<T>::zero()),
            }
        }
        return PolyNumber::new_var(v, &self.var);
    }
}

impl<T> PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
    T: Div<Output = T>,
{
    fn newton_approx(self, k : T, r1 : T) -> T {
        let t = self.tangent(1,r1);
        let o = &t.n[0];
        let s = &t.n[1];
        return (k - (*o).clone()) / (*s).clone();
    }
}

impl<T> PolyNumber<T>
where
    T: Zero,
    T: Clone,
{
    pub fn get(&self, i : usize) -> T {
        let ni = self.n.get(i);
        match ni {
            Some(nni) => return (*nni).clone(),
            None      => return T::zero(),
        }
    }
}

impl<T> PolyNumber<T>
{
    pub fn new(v : Vec<T>) -> PolyNumber<T> {
        return PolyNumber::new_var(v, "x");
    }
    pub fn new_var(v : Vec<T>, var: &str) -> PolyNumber<T> {
        return PolyNumber { n: v, var: var.to_string() };
    }
}

impl<T> fmt::Display for PolyNumber<T>
where
    T: std::fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "(")?;

        for (i, v) in self.n.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{}", v)?;
            if i > 0 {
                write!(f, "*{}", self.var)?;
            }
            if i > 1 {
                write!(f, "^{}", i)?;
            }
        }
        write!(f, ")")
    }
}

impl<T> fmt::Debug for PolyNumber<T>
where
    T: std::fmt::Debug,
    T: std::fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
    use num::{BigInt};
    use assert_approx_eq::assert_approx_eq;
    #[test]
    fn eq_poly_numbers() {
        assert_eq!(PolyNumber::new_var(vec![1, 2, 0], "a" ), PolyNumber::new_var(vec![1, 2, 0, 0], "a" ));
        assert_ne!(PolyNumber::new(vec![5, 2, 0, 3] ), PolyNumber::new(vec![5, 2, 3] ));
    }
    #[test]
    fn addition_of_poly_numbers() {
        let p1 = PolyNumber::new_var(vec![1, 0, 2], "a" );
        let p2 = PolyNumber::new_var(vec![3, 4, 7, 1], "a" );
        let p3 = PolyNumber::new_var(vec![4, 4, 9, 1], "a" );
        assert_eq!(p1+p2,p3);
    }
    #[test]
    fn multiplication_of_poly_number() {
        let p1 = PolyNumber::new_var(vec![1, 0, 2], "a" );
        let p2 = PolyNumber::new_var(vec![3, 4, 7, 1], "a" );
        let p3 = PolyNumber::new_var(vec![3, 4, 13, 9, 14, 2], "a" );
        assert_eq!(p1*p2,p3);
        let p4 = PolyNumber::new_var(vec![0, 1, 5, 2], "a" );
        let p5 = PolyNumber::new_var(vec![4, 3], "a" );
        let p6 = PolyNumber::new_var(vec![0, 4, 23, 23, 6], "a" );
        assert_eq!(p4*p5,p6);
        let p7 = PolyNumber::new_var(vec![2, 1, 3], "a" );
        let p8 = PolyNumber::new_var(vec![5, 7], "a" );
        let p9 = PolyNumber::new_var(vec![10, 19, 22, 21], "a" );
        assert_eq!(p7*p8,p9);
    }
    #[test]
    fn scalars_for_poly_number() {
        let p1 = PolyNumber::new_var(vec![2, 3], "a" );
        let s1 = 5;
        let p2 = PolyNumber::new_var(vec![10, 15], "a" );
        assert_eq!(p1*s1,p2);
    }
    #[test]
    fn evaluate_poly_number() {
        let p1 = PolyNumber::new_var(vec![2, -3, 1], "a" );
        let c1 = 4;
        assert_eq!(p1.eval(c1),6);
        let p2 = PolyNumber::new_var(vec![1, -5], "a" );
        let c2 = 4;
        assert_eq!(p2.eval(c2),-19);
    }
    #[test]
    fn composing_poly_numbers() {
        let p = PolyNumber::new_var(vec![1,3,4], "a" );
        let q = PolyNumber::new_var(vec![2,1], "a" );
        let pq = PolyNumber::new_var(vec![23,19,4], "a" );
        assert_eq!(p.eval(q),pq);
    }
    #[test]
    fn composing_poly_numbers2() {
        let p = PolyNumber::new_var(vec![1,3,4], "a" );
        let q = PolyNumber::new_var(vec![37,1], "a" );
        let pq = PolyNumber::new_var(vec![5588,299,4], "a" );
        assert_eq!(p.eval(q),pq);
    }
    #[test]
    fn left_translate_poly_numbers() {
        let p = PolyNumber::new_var(vec![1,3,4], "a" );
        let L3p = PolyNumber::new_var(vec![46,27,4], "a" );
        assert_eq!(p.ltrans(3),L3p);
    }
    #[test]
    fn adding_bi_poly_numbers() {
        let bp1 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![1,  0, 5, 4], "a" ),
                                        PolyNumber::new_var(vec![3, -1, 7], "a" ),
                                        PolyNumber::new_var(vec![4,  2], "a" ) ], "b");
        let bp2 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![2, -1, 3], "a" ),
                                        PolyNumber::new_var(vec![1,  0, 2, -3], "a" ) ], "b");
        let bp3 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![3, -1, 8,  4], "a" ),
                                        PolyNumber::new_var(vec![4, -1, 9, -3], "a" ) ,
                                        PolyNumber::new_var(vec![4,  2], "a" ) ], "b");
        assert_eq!(bp1+bp2,bp3);
    }
    #[test]
    fn multiplication_of_bi_poly_numbers() {
        let bp1 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![1, 2], "a" ),
                                        PolyNumber::new_var(vec![5, 3], "a" ) ], "b");
        let bp2 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![4, 5, 1], "a" ),
                                        PolyNumber::new_var(vec![2, 1, 3], "a" ) ], "b");
        let bp3 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![ 4, 13, 11, 2], "a" ),
                                        PolyNumber::new_var(vec![22, 42, 25, 9], "a" ),
                                        PolyNumber::new_var(vec![10, 11, 18, 9], "a" ) ], "b");
        assert_eq!(bp1*bp2,bp3);
    }
    #[test]
    fn scalar_multiplication_of_bi_poly_numbers() {
        let bp1 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![1,  2], "a" ),
                                        PolyNumber::new_var(vec![5,  3], "a" ) ], "b");
        let bp2 = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![3,  6], "a" ),
                                        PolyNumber::new_var(vec![15, 9], "a" ) ], "b");
        assert_eq!(bp1*3,bp2);
    }
    #[test]
    fn evaluate_poly_number_at_bi_poly_number() {
        let alpha = create_polynumber_var!(alpha; alpha,beta ; i32);
        let beta  = create_polynumber_var!(beta;  alpha,beta ; i32);
        let p = alpha + beta;
        let q = PolyNumber::new_var(vec![-4,7,10,-6,2], "a" );
        let t = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![-4,  7, 10,-6,2], "a" ),
                                      PolyNumber::new_var(vec![ 7, 20,-18, 8], "a" ),
                                      PolyNumber::new_var(vec![10,-18,12], "a" ),
                                      PolyNumber::new_var(vec![-6,  8], "a" ),
                                      PolyNumber::new_var(vec![2], "a" ) ], "b");
        assert_eq!(q.eval(p),t);
    }
    #[test]
    fn taylor_bi_poly_number() {
        let q = PolyNumber::new(vec![-4,7,10,-6,2]);
        let t = PolyNumber::new_var(vec![ PolyNumber::new_var(vec![-4,  7, 10,-6,2], "a" ),
                                      PolyNumber::new_var(vec![ 7, 20,-18, 8], "a" ),
                                      PolyNumber::new_var(vec![10,-18,12], "a" ),
                                      PolyNumber::new_var(vec![-6,  8], "a" ),
                                      PolyNumber::new_var(vec![2], "a" ) ], "b");
        assert_eq!(q.taylor(),t);
    }
    #[test]
    fn subderivatives_of_poly_numbers() {
        let q  = PolyNumber::new_var(vec![-4,7,10,-6,2], "a" );
        let d1 = PolyNumber::new_var(vec![7,20,-18,8], "a" );
        let d2 = PolyNumber::new_var(vec![10,-18,12], "a" );
        let d3 = PolyNumber::new_var(vec![-6,8], "a" );
        let d4 = PolyNumber::new_var(vec![2], "a" );
        assert_eq!(q.clone().derivative(1),d1);
        assert_eq!(q.clone().derivative(2),d2);
        assert_eq!(q.clone().derivative(3),d3);
        assert_eq!(q.clone().derivative(4),d4);
    }
    #[test]
    fn truncation_of_polynumbers() {
        let p = PolyNumber::new_var(vec![8,-5,0,4,-1], "a" );
        let t1p = PolyNumber::new_var(vec![8,-5], "a" );
        let t2p = PolyNumber::new_var(vec![8,-5], "a" );
        let t3p = PolyNumber::new_var(vec![8,-5,0,4], "a" );
        assert_eq!(p.clone().truncate(1),t1p);
        assert_eq!(p.clone().truncate(2),t2p);
        assert_eq!(p.clone().truncate(3),t3p);
    }
    #[test]
    fn tangents_of_polynumbers() {
        let p = PolyNumber::new_var(vec![8,-5,0,4,-1], "a" );
        let t0p1 = PolyNumber::new_var(vec![6], "a" );
        let t1p1 = PolyNumber::new_var(vec![3,3], "a" );
        let t2p1 = PolyNumber::new_var(vec![9,-9,6], "a" );
        assert_eq!(p.clone().tangent(0,1),t0p1);
        assert_eq!(p.clone().tangent(1,1),t1p1);
        assert_eq!(p.clone().tangent(2,1),t2p1.clone());
        assert_eq!(p.clone().tangent(3,1),t2p1);
    }
    #[test]
    fn taylor_expansion_for_bi_polynumbers() {
        let p = PolyNumber::new_var(vec![PolyNumber::new_var(vec![-1,0,1], "a" ), // -1 + a^2
                                    PolyNumber::new_var(vec![0], "a" ),
                                    PolyNumber::new_var(vec![1], "a" ) ], "b"); // b^2
        let t = PolyNumber::new_var(vec![
                    PolyNumber::new_var(vec![ // d^0
                        PolyNumber::new_var(vec![ // b^0
                            PolyNumber::new_var(vec![-1,0,1], "a" ), // -1 + a^2
                            PolyNumber::new_var(vec![0,2], "a" ), // 2ag
                            PolyNumber::new_var(vec![1], "a" ) ], "g"),  // g^2
                        PolyNumber::new_var(vec![ // b^1
                            PolyNumber::new_var(vec![0], "a" ) ], "g"),
                        PolyNumber::new_var(vec![ // b^2
                            PolyNumber::new_var(vec![1], "a" ) ], "g") ], "b"), // b^2
                    PolyNumber::new_var(vec![ // d^1
                        PolyNumber::new_var(vec![ // b^0
                            PolyNumber::new_var(vec![0], "a" ) ], "g"),
                        PolyNumber::new_var(vec![ // b^1
                            PolyNumber::new_var(vec![2], "a" ) ], "g") ], "b" ), // 2bd
                    PolyNumber::new_var(vec![ // d^2
                        PolyNumber::new_var(vec![ // b^0
                            PolyNumber::new_var(vec![1], "a" ) ], "g") ], "b") ], "d"); // d^2
        assert_eq!(p.taylor2(),t);
    }
    #[test]
    fn taylor_expansion_for_bi_polynumbers_with_macros() {
        // create bi-polynumbers
        let ba = create_polynumber_var!(a; a,b ; i32);
        let bb  = create_polynumber_var!(b;  a,b ; i32);
        let bone = create_polynumber_one!(a,b ; i32);

        let p = ba.clone()*ba.clone() - bone.clone() + bb.clone()*bb.clone();

        // create quad-polynumbers
        let qa = create_polynumber_var!(a; a,g,b,d ; i32);
        let qb  = create_polynumber_var!(b;  a,g,b,d ; i32);
        let qd = create_polynumber_var!(d; a,g,b,d ; i32);
        let qg  = create_polynumber_var!(g;  a,g,b,d ; i32);
        let qone = create_polynumber_one!(a,g,b,d ; i32);
        let qtwo = qone.clone() + qone.clone();

        let t = qa.clone()*qa.clone() - qone.clone()
                + qtwo.clone()*qa.clone()*qg.clone() + qg.clone()*qg.clone()
                + qb.clone()*qb.clone() + qtwo.clone()*qb.clone()*qd.clone()
                + qd.clone()*qd.clone();

        assert_eq!(p.taylor2(),t);
    }
    #[test]
    fn derivatives_of_by_polynumbers() {
        let p = PolyNumber::new_var(vec![PolyNumber::new_var(vec![-1,0,1], "a" ), // -1 + a^2
                                    PolyNumber::new_var(vec![0], "a" ),
                                    PolyNumber::new_var(vec![1], "a" ) ], "b"); // b^2
        let d01 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![0], "a" ),
                                      PolyNumber::new_var(vec![2], "a" ) ], "b"); // 2b
        let d02 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![1], "a" ) ], "b"); // 1
        let d10 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![0,2], "a" ) ], "b"); // 2a
        let d11 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![0], "a" ) ], "b"); // 0
        let d20 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![1], "a" ) ], "b"); // 1
        assert_eq!(p.clone().derivative2(0,0),p.clone());
        assert_eq!(p.clone().derivative2(0,1),d01);
        assert_eq!(p.clone().derivative2(0,2),d02);
        assert_eq!(p.clone().derivative2(1,0),d10);
        assert_eq!(p.clone().derivative2(1,1),d11);
        assert_eq!(p.clone().derivative2(2,0),d20);
    }
    #[test]
    fn tangent_plane() {
        let p = PolyNumber::new_var(vec![PolyNumber::new_var(vec![-1,0,1], "a" ), // -1 + a^2
                                    PolyNumber::new_var(vec![0], "a" ),
                                    PolyNumber::new_var(vec![1], "a" ) ], "b"); // b^2
        let p35 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![33, 6, 1], "a" ), // 23 + 6a + a^2
                                      PolyNumber::new_var(vec![10], "a" ), // 10b
                                      PolyNumber::new_var(vec![1], "a" ) ], "b"); // b^2
        // first tangent,i.e. tangent plane, of p at [3,5]
        let t35 = PolyNumber::new_var(vec![PolyNumber::new_var(vec![-35, 6], "a" ), // -35 + 6a
                                      PolyNumber::new_var(vec![10], "a" ) ], "b"); // 10b
        assert_eq!(p.clone().ltrans2(3,5),p35);
        assert_eq!(p.tangent2(1,3,5),t35);
    }
    #[test]
    fn newton_approximation() {
        let one = Ratio::new(1,1);
        let p = PolyNumber::new(vec![one*0, one*0, one]);
        let a1 = p.clone().newton_approx(one*2,one*2);
        assert_eq!(a1,one*3/2);
        let a2 = p.clone().newton_approx(one*2,a1);
        assert_eq!(a2,one*17/12);
        let a3 = p.clone().newton_approx(one*2,a2);
        assert_eq!(a3,one*577/408);
        let a4 = p.clone().newton_approx(one*2,a3);
        assert_eq!(a4,one*665857/470832);
    }
    #[test]
    fn cube_root_of_five() {
        let bone = BigInt::from(1);
        let bzero = BigInt::from(0);
        let one : Ratio<BigInt> = Ratio::new(bone.clone(),bone.clone());
        let zero : Ratio<BigInt> = Ratio::new(bzero.clone(),bone.clone());
        let five : Ratio<BigInt> = one.clone()*(bone.clone()*5);
        let p = PolyNumber::new(vec![zero.clone(), zero.clone(), zero.clone(), one.clone()]);
        let a1 = p.clone().newton_approx(five.clone(),one.clone());
        assert_eq!(a1,one.clone()*(bone.clone()*7)/(bone.clone()*3));
        let a2 = p.clone().newton_approx(five.clone(),a1);
        assert_eq!(a2,one.clone()*(bone.clone()*821)/(bone.clone()*441));
        let a3 = p.clone().newton_approx(five.clone(),a2);
        assert_eq!(a3,one.clone()*(bone.clone()*1535605927)/(bone.clone()*891756243));
        let a4 = p.clone().newton_approx(five.clone(),a3.clone());
    }
    #[test]
    fn float_cube_root_of_five() {
        let p = PolyNumber::new_var(vec![0.0f64, 0.0f64, 0.0f64, 1.0f64], "a" );
        let a1 = p.clone().newton_approx(5.0f64,1.0f64);
        assert_approx_eq!(a1, 7.0f64/3.0f64, 0.0000001f64);
        let a2 = p.clone().newton_approx(5.0f64,a1);
        assert_approx_eq!(a2, 821.0f64/441.0f64, 0.0000001f64);
        let a3 = p.clone().newton_approx(5.0f64,a2);
        assert_approx_eq!(a3, 1535605927.0f64/891756243.0f64, 0.0000001f64);
        let a4 = p.clone().newton_approx(5.0f64,a3);
    }
}

