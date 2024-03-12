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
use std::fmt;
use polynumber::PolyNumber;
use extrational::RatInf;



/// Represents the ratio between two poly numbers
#[derive(Clone)]
pub struct PolyRatio<T> {
    numer_factors: Vec<PolyNumber<T>>,
    denom_factors: Vec<PolyNumber<T>>,
}

impl<T> PolyRatio<T>
where
    PolyNumber<T>: One,
    PolyNumber<T>: Mul<Output = PolyNumber<T>>,
    T: Clone,
{
    pub fn new(numer: PolyNumber<T>, denom: PolyNumber<T>) -> PolyRatio<T> {
        let numer_factors = vec![numer.clone()];
        let denom_factors = vec![denom.clone()];
        PolyRatio{numer_factors: numer_factors, denom_factors: denom_factors}
    }
    pub fn newf(numer_factors: Vec<PolyNumber<T>>, denom_factors: Vec<PolyNumber<T>>) -> PolyRatio<T> {
        PolyRatio{numer_factors: numer_factors, denom_factors: denom_factors}
    }
    pub fn numer(&self) -> PolyNumber<T> {
        let mut numer = PolyNumber::<T>::one();
        for (_nfi, nfv) in self.numer_factors.iter().enumerate() {
            numer = numer * nfv.clone();
        }
        numer
    }
    pub fn denom(&self) -> PolyNumber<T> {
        let mut denom = PolyNumber::<T>::one();
        for (_dfi, dfv) in self.denom_factors.iter().enumerate() {
            denom = denom * dfv.clone();
        }
        denom
    }
}


impl<T> PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Div + Clone,
    T: Div<Output = T>,
    PolyNumber<T>: Sub<Output = PolyNumber<T>>,
{
    pub fn reduce(numer: PolyNumber<T>, denom: PolyNumber<T>) -> PolyRatio<T> {
        let dlo = denom.lowest_order();
        let dlow = match dlo {
            Some(i) => i,
            None => panic!("Rational poly number has a denominator equal zero."),
        };
        let d1 = denom.get(dlow);

        let mut res = numer.clone();
        let mut q: Vec<T> = Vec::new();


        if numer.order() < denom.order() {
            return PolyRatio::new(numer, denom);
        }

        let order = numer.order()-denom.order()+1;

        for _i in 0..order {
            let nlo = res.lowest_order();

            let nlow = match nlo {
                Some(i) => i,
                None => 0,
            };


            if nlow < dlow {
                return PolyRatio::new(numer.clone(), denom);
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
                    return PolyRatio::new(pq1, PolyNumber::new( vec![T::one()] ) );
                }

            } else {
                return PolyRatio::new(numer, denom);
            }

        }
        PolyRatio::new(numer, denom)
    }
}

impl<T> PartialEq for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let p = self.numer();
        let q = self.denom();
        let r = other.numer();
        let s = other.denom();
        p * s == r * q
    }
}

impl<T> Add for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn add(self, other: Self) -> PolyRatio<T> {
        let mut sdf = self.denom_factors.clone(); // reduced self denominator factors
        let mut odf = other.denom_factors.clone(); // reduced other denominator factors
        let mut cdf = PolyNumber::<T>::one(); // common denominator factor
        let mut denoms_have_common_factor = false;

        for (sfi, sfv) in self.denom_factors.iter().enumerate() {
            for (ofi, ofv) in other.denom_factors.iter().enumerate() {
                if sfv == ofv {
                    cdf = sfv.clone();
                    sdf.remove(sfi);
                    odf.remove(ofi);
                    denoms_have_common_factor = true;
                    break;
                }
            }
            if denoms_have_common_factor {
                break;
            }
        }

        if denoms_have_common_factor {
            let s1 = PolyRatio::newf(self.numer_factors.clone(), sdf);
            let s2 = PolyRatio::newf(other.numer_factors.clone(), odf);
            let s3 = PolyRatio::new(cdf, PolyNumber::<T>::one());
            return (s1+s2)/s3;
        }

        let mut snf = self.numer_factors.clone(); // reduced self numerator factors
        let mut onf = other.numer_factors.clone(); // reduced other numerator factors
        let mut cnf = PolyNumber::<T>::one(); // common numerator factor
        let mut numerators_have_common_factor = false;

        for (sfi, sfv) in self.numer_factors.iter().enumerate() {
            for (ofi, ofv) in other.numer_factors.iter().enumerate() {
                if sfv == ofv {
                    cnf = sfv.clone();
                    snf.remove(sfi);
                    onf.remove(ofi);
                    numerators_have_common_factor = true;
                    break;
                }
            }
            if numerators_have_common_factor {
                break;
            }
        }

        if numerators_have_common_factor {
            let s1 = PolyRatio::newf(snf, self.denom_factors.clone());
            let s2 = PolyRatio::newf(onf, other.denom_factors.clone());
            let s3 = PolyRatio::new(cnf, PolyNumber::<T>::one());
            return (s1+s2)*s3;
        }

        let p = self.numer();
        let q = self.denom();
        let r = other.numer();
        let s = other.denom();

        let mut qf = self.denom_factors.clone();
        let mut sf = other.denom_factors.clone();
        qf.append(&mut sf);

        PolyRatio::newf( vec![p*s.clone() + r*q.clone()], qf )
    }
}

impl<T> Mul for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn mul(self, other: Self) -> PolyRatio<T> {
        let mut p = self.numer_factors.clone();
        let mut q = self.denom_factors.clone();
        let mut r = other.numer_factors.clone();
        let mut s = other.denom_factors.clone();
        p.append(&mut r);
        q.append(&mut s);
        PolyRatio::newf(p, q)
    }
}

impl<T> Mul<T> for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn mul(self, other: T) -> PolyRatio<T> {
        let mut p = self.numer_factors.clone();
        let q = self.denom_factors.clone();
        p.push(PolyNumber::<T>::one() * other);
        PolyRatio::newf( p, q )
    }
}

impl<T> Sub for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn sub(self, other: Self) -> PolyRatio<T> {
        let mone = T::zero() - T::one();
        self.clone() + (other.clone() * mone)
    }
}

impl<T> Div for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    type Output = PolyRatio<T>;

    fn div(self, other: Self) -> PolyRatio<T> {
        let mut p = self.numer_factors.clone();
        let mut q = self.denom_factors.clone();
        let mut r = other.numer_factors.clone();
        let mut s = other.denom_factors.clone();
        p.append(&mut s);
        r.append(&mut q);
        PolyRatio::newf( p, r )
    }
}

impl<T> Zero for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    fn zero() -> PolyRatio<T> {
        return PolyRatio::new(PolyNumber::<T>::zero(),
                              PolyNumber::<T>::one() );
    }

    fn is_zero(&self) -> bool {
        return *self == Self::zero();
    }
}

impl<T> One for PolyRatio<T>
where
    T: PartialEq + Zero + One + Mul + Add + Sub + Clone,
    T: Sub<Output = T>,
{
    fn one() -> PolyRatio<T> {
        return PolyRatio::new(PolyNumber::<T>::one(),
                              PolyNumber::<T>::one() );
    }
}

impl<T> PolyRatio<RatInf<T>>
where
    T: Integer,
    T: Clone,
{
    pub fn eval(&self, c: RatInf<T>) -> RatInf<T> {
        let mut numer = self.numer();
        let mut denom = self.denom();

        loop {
            let num = numer.eval(c.clone());
            let den = denom.eval(c.clone());

            if den.is_zero() {
                if !num.is_zero() {
                    return num/den;
                }

                let numer_reduced = PolyRatio::reduce(numer.clone(),
                                    PolyNumber::new(vec![-c.clone(), RatInf::<T>::one()]));
                let denom_reduced = PolyRatio::reduce(denom.clone(),
                                    PolyNumber::new(vec![-c.clone(), RatInf::<T>::one()]));

                numer = numer_reduced.numer();
                denom = denom_reduced.numer();
            } else {
                return num/den;
            }
        }
    }
}

impl<T> fmt::Debug for PolyRatio<T>
where
    T: std::fmt::Debug,
    T: std::fmt::Display,
    T: Clone,
    PolyNumber<T>: Mul,
    PolyNumber<T>: One,

{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PolyRatio")
         .field("numerator", &self.numer())
         .field("denominator", &self.denom())
         .finish()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reduce_rational_poly_numbers() {
        let p1 = PolyNumber::new( vec![2, 7, 2, -3] );
        let p2 = PolyNumber::new( vec![2, 1, -1] );
        let p3 = PolyNumber::new( vec![1, 3] );
        let p4 = PolyNumber::one();
        let rp1 = PolyRatio::reduce(p1,p2);
        assert_eq!(rp1.numer(),p3);
        assert_eq!(rp1.denom(),p4);

        let p21 = PolyNumber::new( vec![12, 8, -7, -2, 1] );
        let p22 = PolyNumber::new( vec![4, 0, -1] );
        let p23 = PolyNumber::new( vec![3, 2, -1] );
        let p24 = PolyNumber::one();
        let rp21 = PolyRatio::reduce(p21,p22);
        assert_eq!(rp21.numer(),p23);
        assert_eq!(rp21.denom(),p24);
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
        let rp1 = PolyRatio::new(PolyNumber::new( vec![2,1] ), PolyNumber::new( vec![3,-1] ) );
        let rp2 = PolyRatio::new(PolyNumber::new( vec![4,0,-1]), PolyNumber::new( vec![1, 1] ) );
        let rp3 = PolyRatio::new(PolyNumber::new( vec![14,-1,-2,1]), PolyNumber::new( vec![3,2,-1] ) );
        assert_eq!(rp1+rp2,rp3);
    }
    #[test]
    fn arithmetic_with_rat_poly_numbers2() {
        let rp1 = PolyRatio::new(PolyNumber::new( vec![5,-1,0,1] ), PolyNumber::new( vec![1,0,0,0,-1] ) );
        let rp2 = PolyRatio::new(PolyNumber::new( vec![6,0,-1] ), PolyNumber::new( vec![0,0,0,0,0,1] ) );
        let rp3 = PolyRatio::new(PolyNumber::new( vec![30,-6,-5,7,0,-1] ),
                                 PolyNumber::new( vec![0,0,0,0,0,1,0,0,0,-1] ) );
        assert_eq!(rp1*rp2,rp3);
    }
    #[test]
    fn equality_of_rational_polynumbers() {
        let rp1 = PolyRatio::new(PolyNumber::new( vec![1,0,-1] ), PolyNumber::new( vec![1,-1] ) );
        let rp2 = PolyRatio::new(PolyNumber::new( vec![1,1] ), PolyNumber::one() );
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn equality_of_rational_polynumbers2() {
        let rp1 = PolyRatio::new(PolyNumber::new( vec![1,0,0,-1] ), PolyNumber::new( vec![1,-1] ) );
        let rp2 = PolyRatio::new(PolyNumber::new( vec![1,1,1] ), PolyNumber::one() );
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn equality_of_rational_polynumbers3() {
        let rp1 = PolyRatio::new(PolyNumber::new( vec![1,0,1,0,1] ), PolyNumber::new( vec![1,1,1] ) );
        let rp2 = PolyRatio::new(PolyNumber::new( vec![1,-2,1,0,-1] ), PolyNumber::new( vec![1,-1,-1] ) );
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn folium_of_descartes() {
        let p = PolyNumber::new( vec![PolyNumber::<i64>::new( vec![0,0,0,1] ), // x^3
                                    PolyNumber::new( vec![0,3] ), // 3xy
                                    PolyNumber::new( vec![0] ),
                                    PolyNumber::new( vec![1] ) ] ); // y^3
        let tx = PolyRatio::new(PolyNumber::new( vec![0,-3] ),
                                PolyNumber::new( vec![1,0,0,1] ) ); // -3t/(1+t^3)
        let ty = PolyRatio::new(PolyNumber::new( vec![0,0,-3] ),
                                PolyNumber::new( vec![1,0,0,1] ) ); // -3t^2/(1+t^3)
        let zero = PolyRatio::new(PolyNumber::new( vec![0] ),
                                  PolyNumber::new( vec![1] ) );
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

