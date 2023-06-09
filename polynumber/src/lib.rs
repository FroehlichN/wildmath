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
use std::fmt::Debug;


#[derive(Debug, Clone)]
pub struct PolyNumber<T> {
    n: Vec<T>,
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

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push((*aa).clone() + (*bb).clone()),
                (Some(aa), None)     => s.push((*aa).clone()),
                (None, Some(bb))     => s.push((*bb).clone()),
                (None, None)         => return PolyNumber { n: s },
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

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push((*aa).clone() - (*bb).clone()),
                (Some(aa), None)     => s.push((*aa).clone()),
                (None, Some(bb))     => s.push(T::zero() - (*bb).clone()),
                (None, None)         => return PolyNumber { n: s },
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

        return PolyNumber { n: p };
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

        return PolyNumber { n: p };
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
            p.push(a * PolyNumber{ n: vec![other.clone()] });
        }

        return PolyNumber { n: p };
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

        return PolyNumber{ n: v };
    }
}

impl<T> PolyNumber<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn ltrans(&self, c: T) -> PolyNumber<T> {
        let q = PolyNumber { n: vec![c,T::one()] };
        return self.eval(q);
    }
}

impl<T> PolyNumber<PolyNumber<T>>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn ltrans2(&self, r: T, s: T) -> PolyNumber<PolyNumber<T>> {
        let a = PolyNumber { n: vec![r,T::one()] }; // r + a
        let b = PolyNumber { n: vec![PolyNumber{ n: vec![s] }, // s
                                     PolyNumber{ n: vec![T::one()] } ] }; // b
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
        return PolyNumber { n: vec![T::zero()] };
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
        return PolyNumber { n: vec![T::one()] };
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
        return PolyNumber{ n: v };
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
        return PolyNumber{ n: v };
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
        let p = PolyNumber { n: vec! [ PolyNumber { n: vec![T::zero(),T::one()] },
                                       PolyNumber { n: vec![T::one(),T::zero()] } ] };
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
        let apg = PolyNumber{ n: vec![PolyNumber{ n: vec![T::zero(),T::one()] }, // a
                                      PolyNumber{ n: vec![T::one(),T::zero()] } ] }; // g
        let bpd = PolyNumber{ n: vec![
                      PolyNumber{ n: vec![ // d^0
                          PolyNumber{ n: vec![ // b^0
                              PolyNumber{ n: vec![T::zero()] } ] },
                          PolyNumber{ n: vec![ // b^1
                              PolyNumber{ n: vec![T::one()] } ] } ] }, // b
                      PolyNumber{ n: vec![ // d^1
                          PolyNumber{ n: vec![ // b^0
                              PolyNumber{ n: vec![T::one()] } ] } ] } ] }; // d

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
        let tp = self.taylor2();

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
        return PolyNumber{ n: v };
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
        return PolyNumber{ n: v };
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
        assert_eq!(PolyNumber { n: vec![1, 2, 0] }, PolyNumber { n: vec![1, 2, 0, 0] });
        assert_ne!(PolyNumber { n: vec![5, 2, 0, 3] }, PolyNumber { n: vec![5, 2, 3] });
    }
    #[test]
    fn addition_of_poly_numbers() {
        let p1 = PolyNumber { n: vec![1, 0, 2] };
        let p2 = PolyNumber { n: vec![3, 4, 7, 1] };
        let p3 = PolyNumber { n: vec![4, 4, 9, 1] };
        assert_eq!(p1+p2,p3);
    }
    #[test]
    fn multiplication_of_poly_number() {
        let p1 = PolyNumber { n: vec![1, 0, 2] };
        let p2 = PolyNumber { n: vec![3, 4, 7, 1] };
        let p3 = PolyNumber { n: vec![3, 4, 13, 9, 14, 2] };
        assert_eq!(p1*p2,p3);
        let p4 = PolyNumber { n: vec![0, 1, 5, 2] };
        let p5 = PolyNumber { n: vec![4, 3] };
        let p6 = PolyNumber { n: vec![0, 4, 23, 23, 6] };
        assert_eq!(p4*p5,p6);
        let p7 = PolyNumber { n: vec![2, 1, 3] };
        let p8 = PolyNumber { n: vec![5, 7] };
        let p9 = PolyNumber { n: vec![10, 19, 22, 21] };
        assert_eq!(p7*p8,p9);
    }
    #[test]
    fn scalars_for_poly_number() {
        let p1 = PolyNumber { n: vec![2, 3] };
        let s1 = 5;
        let p2 = PolyNumber { n: vec![10, 15] };
        assert_eq!(p1*s1,p2);
    }
    #[test]
    fn evaluate_poly_number() {
        let p1 = PolyNumber { n: vec![2, -3, 1] };
        let c1 = 4;
        assert_eq!(p1.eval(c1),6);
        let p2 = PolyNumber { n: vec![1, -5] };
        let c2 = 4;
        assert_eq!(p2.eval(c2),-19);
    }
    #[test]
    fn composing_poly_numbers() {
        let p = PolyNumber { n: vec![1,3,4] };
        let q = PolyNumber { n: vec![2,1] };
        let pq = PolyNumber { n: vec![23,19,4] };
        assert_eq!(p.eval(q),pq);
    }
    #[test]
    fn composing_poly_numbers2() {
        let p = PolyNumber { n: vec![1,3,4] };
        let q = PolyNumber { n: vec![37,1] };
        let pq = PolyNumber { n: vec![5588,299,4] };
        assert_eq!(p.eval(q),pq);
    }
    #[test]
    fn left_translate_poly_numbers() {
        let p = PolyNumber { n: vec![1,3,4] };
        let L3p = PolyNumber { n: vec![46,27,4] };
        assert_eq!(p.ltrans(3),L3p);
    }
    #[test]
    fn adding_bi_poly_numbers() {
        let bp1 = PolyNumber { n: vec![ PolyNumber { n: vec![1,  0, 5, 4] },
                                        PolyNumber { n: vec![3, -1, 7] },
                                        PolyNumber { n: vec![4,  2] } ] };
        let bp2 = PolyNumber { n: vec![ PolyNumber { n: vec![2, -1, 3] },
                                        PolyNumber { n: vec![1,  0, 2, -3] } ] };
        let bp3 = PolyNumber { n: vec![ PolyNumber { n: vec![3, -1, 8,  4] },
                                        PolyNumber { n: vec![4, -1, 9, -3] } ,
                                        PolyNumber { n: vec![4,  2] } ] };
        assert_eq!(bp1+bp2,bp3);
    }
    #[test]
    fn multiplication_of_bi_poly_numbers() {
        let bp1 = PolyNumber { n: vec![ PolyNumber { n: vec![1, 2] },
                                        PolyNumber { n: vec![5, 3] } ] };
        let bp2 = PolyNumber { n: vec![ PolyNumber { n: vec![4, 5, 1] },
                                        PolyNumber { n: vec![2, 1, 3] } ] };
        let bp3 = PolyNumber { n: vec![ PolyNumber { n: vec![ 4, 13, 11, 2] },
                                        PolyNumber { n: vec![22, 42, 25, 9] },
                                        PolyNumber { n: vec![10, 11, 18, 9] } ] };
        assert_eq!(bp1*bp2,bp3);
    }
    #[test]
    fn scalar_multiplication_of_bi_poly_numbers() {
        let bp1 = PolyNumber { n: vec![ PolyNumber { n: vec![1,  2] },
                                        PolyNumber { n: vec![5,  3] } ] };
        let bp2 = PolyNumber { n: vec![ PolyNumber { n: vec![3,  6] },
                                        PolyNumber { n: vec![15, 9] } ] };
        assert_eq!(bp1*3,bp2);
    }
    #[test]
    fn evaluate_poly_number_at_bi_poly_number() {
        let q = PolyNumber { n: vec![-4,7,10,-6,2] };
        let p = PolyNumber { n: vec![ PolyNumber { n: vec![0,1] },
                                      PolyNumber { n: vec![1,0] } ] };
        let t = PolyNumber { n: vec![ PolyNumber { n: vec![-4,  7, 10,-6,2] },
                                      PolyNumber { n: vec![ 7, 20,-18, 8] },
                                      PolyNumber { n: vec![10,-18,12] },
                                      PolyNumber { n: vec![-6,  8] },
                                      PolyNumber { n: vec![2] } ] };
        assert_eq!(q.eval(p),t);
    }
    #[test]
    fn taylor_bi_poly_number() {
        let q = PolyNumber { n: vec![-4,7,10,-6,2] };
        let t = PolyNumber { n: vec![ PolyNumber { n: vec![-4,  7, 10,-6,2] },
                                      PolyNumber { n: vec![ 7, 20,-18, 8] },
                                      PolyNumber { n: vec![10,-18,12] },
                                      PolyNumber { n: vec![-6,  8] },
                                      PolyNumber { n: vec![2] } ] };
        assert_eq!(q.taylor(),t);
    }
    #[test]
    fn subderivatives_of_poly_numbers() {
        let q  = PolyNumber { n: vec![-4,7,10,-6,2] };
        let d1 = PolyNumber { n: vec![7,20,-18,8] };
        let d2 = PolyNumber { n: vec![10,-18,12] };
        let d3 = PolyNumber { n: vec![-6,8] };
        let d4 = PolyNumber { n: vec![2] };
        assert_eq!(q.clone().derivative(1),d1);
        assert_eq!(q.clone().derivative(2),d2);
        assert_eq!(q.clone().derivative(3),d3);
        assert_eq!(q.clone().derivative(4),d4);
    }
    #[test]
    fn truncation_of_polynumbers() {
        let p = PolyNumber { n: vec![8,-5,0,4,-1] };
        let t1p = PolyNumber { n: vec![8,-5] };
        let t2p = PolyNumber { n: vec![8,-5] };
        let t3p = PolyNumber { n: vec![8,-5,0,4] };
        assert_eq!(p.clone().truncate(1),t1p);
        assert_eq!(p.clone().truncate(2),t2p);
        assert_eq!(p.clone().truncate(3),t3p);
    }
    #[test]
    fn tangents_of_polynumbers() {
        let p = PolyNumber{ n: vec![8,-5,0,4,-1] };
        let t0p1 = PolyNumber{ n: vec![6] };
        let t1p1 = PolyNumber{ n: vec![3,3] };
        let t2p1 = PolyNumber{ n: vec![9,-9,6] };
        assert_eq!(p.clone().tangent(0,1),t0p1);
        assert_eq!(p.clone().tangent(1,1),t1p1);
        assert_eq!(p.clone().tangent(2,1),t2p1.clone());
        assert_eq!(p.clone().tangent(3,1),t2p1);
    }
    #[test]
    fn taylor_expansion_for_bi_polynumbers() {
        let p = PolyNumber{ n: vec![PolyNumber{ n: vec![-1,0,1] }, // -1 + a^2
                                    PolyNumber{ n: vec![0] },
                                    PolyNumber{ n: vec![1] } ] }; // b^2
        let t = PolyNumber{ n: vec![
                    PolyNumber{ n: vec![ // d^0
                        PolyNumber{ n: vec![ // b^0
                            PolyNumber{ n: vec![-1,0,1] }, // -1 + a^2
                            PolyNumber{ n: vec![0,2] }, // 2ag
                            PolyNumber{ n: vec![1] } ] },  // g^2
                        PolyNumber{ n: vec![ // b^1
                            PolyNumber{ n: vec![0] } ] },
                        PolyNumber{ n: vec![ // b^2
                            PolyNumber{ n: vec![1] } ] } ] }, // b^2
                    PolyNumber{ n: vec![ // d^1
                        PolyNumber{ n: vec![ // b^0
                            PolyNumber{ n: vec![0] } ] },
                        PolyNumber{ n: vec![ // b^1
                            PolyNumber{ n: vec![2] } ] } ] }, // 2bd
                    PolyNumber{ n: vec![ // d^2
                        PolyNumber{ n: vec![ // b^0
                            PolyNumber{ n: vec![1] } ] } ] } ] }; // d^2
        assert_eq!(p.taylor2(),t);
    }
    #[test]
    fn derivatives_of_by_polynumbers() {
        let p = PolyNumber{ n: vec![PolyNumber{ n: vec![-1,0,1] }, // -1 + a^2
                                    PolyNumber{ n: vec![0] },
                                    PolyNumber{ n: vec![1] } ] }; // b^2
        let d01 = PolyNumber{ n: vec![PolyNumber{ n: vec![0] },
                                      PolyNumber{ n: vec![2] } ] }; // 2b
        let d02 = PolyNumber{ n: vec![PolyNumber{ n: vec![1] } ] }; // 1
        let d10 = PolyNumber{ n: vec![PolyNumber{ n: vec![0,2] } ] }; // 2a
        let d11 = PolyNumber{ n: vec![PolyNumber{ n: vec![0] } ] }; // 0
        let d20 = PolyNumber{ n: vec![PolyNumber{ n: vec![1] } ] }; // 1
        assert_eq!(p.clone().derivative2(0,0),p.clone());
        assert_eq!(p.clone().derivative2(0,1),d01);
        assert_eq!(p.clone().derivative2(0,2),d02);
        assert_eq!(p.clone().derivative2(1,0),d10);
        assert_eq!(p.clone().derivative2(1,1),d11);
        assert_eq!(p.clone().derivative2(2,0),d20);
    }
    #[test]
    fn tangent_plane() {
        let p = PolyNumber{ n: vec![PolyNumber{ n: vec![-1,0,1] }, // -1 + a^2
                                    PolyNumber{ n: vec![0] },
                                    PolyNumber{ n: vec![1] } ] }; // b^2
        let p35 = PolyNumber{ n: vec![PolyNumber{ n: vec![33, 6, 1] }, // 23 + 6a + a^2
                                      PolyNumber{ n: vec![10] }, // 10b
                                      PolyNumber{ n: vec![1] } ] }; // b^2
        // first tangent,i.e. tangent plane, of p at [3,5]
        let t35 = PolyNumber{ n: vec![PolyNumber{ n: vec![-35, 6] }, // -35 + 6a
                                      PolyNumber{ n: vec![10] } ] }; // 10b
        assert_eq!(p.clone().ltrans2(3,5),p35);
        assert_eq!(p.tangent2(1,3,5),t35);
    }
    #[test]
    fn newton_approximation() {
        let one = Ratio::new(1,1);
        let p = PolyNumber{ n: vec![one*0, one*0, one] };
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
        let p = PolyNumber{ n: vec![zero.clone(), zero.clone(), zero.clone(), one.clone()] };
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
        let p = PolyNumber{ n: vec![0.0f64, 0.0f64, 0.0f64, 1.0f64] };
        let a1 = p.clone().newton_approx(5.0f64,1.0f64);
        assert_approx_eq!(a1, 7.0f64/3.0f64, 0.0000001f64);
        let a2 = p.clone().newton_approx(5.0f64,a1);
        assert_approx_eq!(a2, 821.0f64/441.0f64, 0.0000001f64);
        let a3 = p.clone().newton_approx(5.0f64,a2);
        assert_approx_eq!(a3, 1535605927.0f64/891756243.0f64, 0.0000001f64);
        let a4 = p.clone().newton_approx(5.0f64,a3);
    }
}

