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


use num::{Zero,One,Integer};
use std::ops::{Mul, Add, Sub, Div};
use std::fmt::Debug;
use polynumber::PolyNumber;
use extrational::RatInf;



/// Represents the ratio between two poly numbers
#[derive(Debug,Clone)]
pub struct PolyRatio<T> {
    numer: PolyNumber<T>,
    denom: PolyNumber<T>,
}

impl<T> PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    PolyNumber<T>: Sub<Output = PolyNumber<T>>,
{
    pub fn new(numer: PolyNumber<T>, denom: PolyNumber<T>) -> PolyRatio<T> {
        let dlo = denom.lowest_order();
        let dlow = match dlo {
            Some(i) => i,
            None => panic!("Rational poly number has a denominator equal zero."),
        };
        let d1 = denom.get(dlow);

        let mut res = numer.clone();
        let mut q: Vec<T> = Vec::new();


        if numer.order() < denom.order() {
            return PolyRatio{numer: numer, denom: denom};
        }

        let order = numer.order()-denom.order()+1;

        for _i in 0..order {
            let nlo = res.lowest_order();

            let nlow = match nlo {
                Some(i) => i,
                None => 0,
            };


            if nlow < dlow {
                return PolyRatio{numer: numer.clone(), denom: denom};
            }

            if nlow > 0 {
                for _j in q.clone().len()..nlow {
                    q.push(T::zero());
                }
            }

            let n1 = res.get(nlow);

            let q1 = n1.clone()/d1.clone();

            if q1.clone()*d1.clone() == n1 {
                q.push(q1);
                let pq1 = PolyNumber::new( q.clone() );
                let s1 = denom.clone() * PolyNumber::new( q.clone() );
                res = numer.clone() - s1;
                let zero : PolyNumber<T> =  PolyNumber::zero();

                if res == zero {
                    return PolyRatio{numer: pq1, denom: PolyNumber::new( vec![T::one()] ) };
                }

            } else {
                return PolyRatio{numer: numer, denom: denom};
            }

        }
        PolyRatio{numer: numer, denom: denom}
    }
}

impl<T> PartialEq for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let p = self.numer.clone();
        let q = self.denom.clone();
        let r = other.numer.clone();
        let s = other.denom.clone();
        p * s == r * q
    }
}

impl<T> Add for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn add(self, other: Self) -> PolyRatio<T> {
        let p = self.numer.clone();
        let q = self.denom.clone();
        let r = other.numer.clone();
        let s = other.denom.clone();
        PolyRatio::new( p*s.clone() + r*q.clone(), q*s )
    }
}

impl<T> Mul for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn mul(self, other: Self) -> PolyRatio<T> {
        let p = self.numer.clone();
        let q = self.denom.clone();
        let r = other.numer.clone();
        let s = other.denom.clone();
        PolyRatio::new( p*r, q*s )
    }
}

impl<T> Mul<T> for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn mul(self, other: T) -> PolyRatio<T> {
        let p = self.numer.clone();
        let q = self.denom.clone();
        PolyRatio::new( p*other, q )
    }
}

impl<T> Sub for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn sub(self, other: Self) -> PolyRatio<T> {
        let p = self.numer.clone();
        let q = self.denom.clone();
        let r = other.numer.clone();
        let s = other.denom.clone();
        PolyRatio::new( p*s.clone() - r*q.clone(), q*s )
    }
}

impl<T> Div for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn div(self, other: Self) -> PolyRatio<T> {
        let p = self.numer.clone();
        let q = self.denom.clone();
        let r = other.numer.clone();
        let s = other.denom.clone();
        PolyRatio::new( p*s, r*q )
    }
}

impl<T> Zero for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    fn zero() -> PolyRatio<T> {
        return PolyRatio{ numer: PolyNumber::<T>::zero(),
                          denom: PolyNumber::<T>::one() };
    }

    fn is_zero(&self) -> bool {
        return *self == Self::zero();
    }
}

impl<T> One for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    T: Sub<Output = T>,
{
    fn one() -> PolyRatio<T> {
        return PolyRatio{ numer: PolyNumber::<T>::one(),
                          denom: PolyNumber::<T>::one() };
    }
}

impl<T> PolyRatio<RatInf<T>>
where
    T: Integer,
    T: Clone,
{
    pub fn eval(&self, c: RatInf<T>) -> RatInf<T> {
        let mut numer = self.numer.clone();
        let mut denom = self.denom.clone();

        loop {
            let num = numer.eval(c.clone());
            let den = denom.eval(c.clone());

            if den.is_zero() {
                if !num.is_zero() {
                    return num/den;
                }

                let numer_reduced = PolyRatio::new(numer.clone(),
                                    PolyNumber::new(vec![-c.clone(), RatInf::<T>::one()]));
                let denom_reduced = PolyRatio::new(denom.clone(),
                                    PolyNumber::new(vec![-c.clone(), RatInf::<T>::one()]));

                numer = numer_reduced.numer;
                denom = denom_reduced.numer;
            } else {
                return num/den;
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn new_rational_poly_numbers() {
        let p1 = PolyNumber::new( vec![2, 7, 2, -3] );
        let p2 = PolyNumber::new( vec![2, 1, -1] );
        let p3 = PolyNumber::new( vec![1, 3] );
        let p4 = PolyNumber::one();
        let rp1 = PolyRatio::new(p1,p2);
        match rp1 {
            PolyRatio { numer: n, denom: d} => {
                assert_eq!(n,p3);
                assert_eq!(d,p4);
            }
        }
        let p21 = PolyNumber::new( vec![12, 8, -7, -2, 1] );
        let p22 = PolyNumber::new( vec![4, 0, -1] );
        let p23 = PolyNumber::new( vec![3, 2, -1] );
        let p24 = PolyNumber::one();
        let rp21 = PolyRatio::new(p21,p22);
        match rp21 {
            PolyRatio { numer: n, denom: d} => {
                assert_eq!(n,p23);
                assert_eq!(d,p24);
            }
        }
    }
    #[test]
    fn equality_of_rational_poly_numbers() {
        let p1 = PolyNumber::new( vec![1, 0, -2] );
        let p2 = PolyNumber::new( vec![2, -5] );
        let rp1 = PolyRatio::new(p1,p2);
        let p3 = PolyNumber::new( vec![2, 0, -4] );
        let p4 = PolyNumber::new( vec![4, -10] );
        let rp2 = PolyRatio::new(p3,p4);
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn arithmetic_with_rat_poly_numbers() {
        let rp1 = PolyRatio{ numer: PolyNumber::new( vec![2,1] ), denom: PolyNumber::new( vec![3,-1] ) };
        let rp2 = PolyRatio{ numer: PolyNumber::new( vec![4,0,-1]), denom: PolyNumber::new( vec![1, 1] ) };
        let rp3 = PolyRatio{ numer: PolyNumber::new( vec![14,-1,-2,1]), denom: PolyNumber::new( vec![3,2,-1] ) };
        assert_eq!(rp1+rp2,rp3);
    }
    #[test]
    fn arithmetic_with_rat_poly_numbers2() {
        let rp1 = PolyRatio{ numer: PolyNumber::new( vec![5,-1,0,1] ), denom: PolyNumber::new( vec![1,0,0,0,-1] ) };
        let rp2 = PolyRatio{ numer: PolyNumber::new( vec![6,0,-1] ), denom: PolyNumber::new( vec![0,0,0,0,0,1] ) };
        let rp3 = PolyRatio{ numer: PolyNumber::new( vec![30,-6,-5,7,0,-1] ), denom: PolyNumber::new( vec![0,0,0,0,0,1,0,0,0,-1] ) };
        assert_eq!(rp1*rp2,rp3);
    }
    #[test]
    fn equality_of_rational_polynumbers() {
        let rp1 = PolyRatio{ numer: PolyNumber::new( vec![1,0,-1] ), denom: PolyNumber::new( vec![1,-1] ) };
        let rp2 = PolyRatio{ numer: PolyNumber::new( vec![1,1] ), denom: PolyNumber::one() };
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn equality_of_rational_polynumbers2() {
        let rp1 = PolyRatio{ numer: PolyNumber::new( vec![1,0,0,-1] ), denom: PolyNumber::new( vec![1,-1] ) };
        let rp2 = PolyRatio{ numer: PolyNumber::new( vec![1,1,1] ), denom: PolyNumber::one() };
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn equality_of_rational_polynumbers3() {
        let rp1 = PolyRatio{ numer: PolyNumber::new( vec![1,0,1,0,1] ), denom: PolyNumber::new( vec![1,1,1] ) };
        let rp2 = PolyRatio{ numer: PolyNumber::new( vec![1,-2,1,0,-1] ), denom: PolyNumber::new( vec![1,-1,-1] ) };
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn folium_of_descartes() {
        let p = PolyNumber::new( vec![PolyNumber::<i128>::new( vec![0,0,0,1] ), // x^3
                                    PolyNumber::new( vec![0,3] ), // 3xy
                                    PolyNumber::new( vec![0] ),
                                    PolyNumber::new( vec![1] ) ] ); // y^3
        let tx = PolyRatio{ numer: PolyNumber::new( vec![0,-3] ),
                            denom: PolyNumber::new( vec![1,0,0,1] ) }; // -3t/(1+t^3)
        let ty = PolyRatio{ numer: PolyNumber::new( vec![0,0,-3] ),
                            denom: PolyNumber::new( vec![1,0,0,1] ) }; // -3t^2/(1+t^3)
        let zero = PolyRatio{ numer: PolyNumber::new( vec![0] ),
                              denom: PolyNumber::new( vec![1] ) };
        let ptx = p.eval2(tx);
        let pt = ptx.eval(ty);
        assert_eq!(pt,zero);
        assert_eq!(p.eval2(2).eval(1),15);
        let t1 = PolyNumber::new( vec![PolyNumber::new( vec![-24,15] ),
                                     PolyNumber::new( vec![9] ) ] );
        assert_eq!(p.clone().tangent2(1,2,1),t1);
        let t2 = PolyNumber::new( vec![PolyNumber::new( vec![3,-4,2] ),
                                     PolyNumber::new( vec![-1,1] ),
                                     PolyNumber::new( vec![1] ) ] ) * 3;
        assert_eq!(p.clone().tangent2(2,2,1),t2);

        let z = RatInf::new( 0, 1);
        let o = RatInf::new( 1, 1);
        assert_eq!(p.eval2(-o*3/2).eval(-o*3/2),z);
        let pr = PolyNumber::new( vec![PolyNumber::new( vec![z,z,z,o] ), // x^3
                                     PolyNumber::new( vec![z,o*3] ), // 3xy
                                     PolyNumber::new( vec![z] ),
                                     PolyNumber::new( vec![o] ) ] ); // y^3
        let t3 = PolyNumber::new( vec![PolyNumber::new( vec![o*27/4,o*9/4] ),
                                     PolyNumber::new( vec![o*9/4] ) ] );
        assert_eq!(pr.clone().tangent2(1,-o*3/2,-o*3/2),t3);
        let t4 = PolyNumber::new( vec![PolyNumber::new( vec![-o*27/4,-o*27/4,-o*9/2] ),
                                     PolyNumber::new( vec![-o*27/4, o*3] ),
                                     PolyNumber::new( vec![-o*9/2] ) ] );
        assert_eq!(pr.clone().tangent2(2,-o*3/2,-o*3/2),t4);
    }
    #[test]
    fn evaluation_of_rational_poly_numbers() {
        let numer = PolyNumber::new(vec![RatInf::new(0,1),RatInf::new(1,1),
                                         RatInf::new(-1,1)]);
        let denom = PolyNumber::new(vec![RatInf::new(1,1),RatInf::new(0,1),
                                         RatInf::new(-1,1)]);
        let p = PolyRatio::new(numer.clone(),denom.clone());

        assert_eq!(p.eval(RatInf::new(0,1)),RatInf::new(0,1));
        assert_eq!(p.eval(RatInf::new(1,1)),RatInf::new(1,2));
        assert!(p.eval(RatInf::new(-1,1)).is_infinite());
    }
    #[test]
    fn evaluation_of_rational_poly_numbers2() {
        let numer = PolyNumber::new(vec![RatInf::new(6,1),RatInf::new(-5,1),RatInf::new(1,1)]);
        let denom = PolyNumber::new(vec![RatInf::new(-16,1),RatInf::new(20,1),
                                         RatInf::new(-8,1),RatInf::new(1,1)]);
        let p = PolyRatio::new(numer.clone(),denom.clone());

        assert_eq!(p.eval(RatInf::new(1,1)),RatInf::new(-2,3));
        assert!(p.eval(RatInf::new(2,1)).is_infinite());
        assert_eq!(p.eval(RatInf::new(3,1)),RatInf::new(0,1));
        assert!(p.eval(RatInf::new(4,1)).is_infinite());
        assert_eq!(p.eval(RatInf::new(5,1)),RatInf::new(2,3));
    }
}

