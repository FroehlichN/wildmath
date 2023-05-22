use num::{Num, Integer};
use std::ops::{Div, Mul, Add, Sub, Neg};
use std::fmt::Debug;


#[derive(Debug)]
enum LinSummand<T> {
    Number(T),
    Unknown(T),
}

/// Represents a linear equation a*x + b + c*x + ... = ... + w + y * x + z
pub struct LinEq<T> {
    lhs: Vec<LinSummand<T>>,
    rhs: Vec<LinSummand<T>>,
}

impl<T> LinEq<T>
where
    T: Num,
    T: Copy,
{
    fn x(&self) -> Option<T> {
       let mut m = T::zero();
       let mut n = T::zero();

       for (_, s) in self.lhs.iter().enumerate() {
            match s {
                LinSummand::Number(a) => n = n - *a,
                LinSummand::Unknown(b) => m = m + *b,
            }
        }
        for (_, s) in self.rhs.iter().enumerate() {
            match s {
                LinSummand::Number(a) => n = n + *a,
                LinSummand::Unknown(b) => m = m - *b,
            }
        }
        if m.is_zero() {
            return None;
        } else {
            return Some(n/m);
        }
    }
}

/// Represents a linear term in x
#[derive(Debug)]
pub struct LinTermInX<T> {
    summands: Vec<LinSummand<T>>,
}

impl<T> Div<T> for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn div(self, other: T) -> LinTermInX<T> {
        let mut q: Vec<LinSummand<T>> = Vec::new();

        for (_, s) in self.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => q.push( LinSummand::Number(*a / other) ),
                LinSummand::Unknown(b) => q.push( LinSummand::Unknown(*b / other) ),
            }
        }

        LinTermInX {summands: q}
    }
}

impl<T> Mul<T> for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn mul(self, other: T) -> LinTermInX<T> {
        let mut p: Vec<LinSummand<T>> = Vec::new();

        for (_, s) in self.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => p.push( LinSummand::Number(*a * other) ),
                LinSummand::Unknown(b) => p.push( LinSummand::Unknown(*b * other) ),
            }
        }

        LinTermInX {summands: p}
    }
}

impl<T> Add<T> for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn add(self, other: T) -> LinTermInX<T> {
        let mut s: Vec<LinSummand<T>> = self.summands;
        s.push( LinSummand::Number(other) );
        LinTermInX {summands: s}
    }
}

impl<T> Sub<T> for LinTermInX<T>
where
    T: Num,
    T: Neg<Output = T>,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn sub(self, other: T) -> LinTermInX<T> {
        let mut d: Vec<LinSummand<T>> = self.summands;
        d.push( LinSummand::Number(-other) );
        LinTermInX {summands: d}
    }
}

impl<T> Sub<LinSummand<T>> for LinTermInX<T>
where
    T: Num,
    T: Neg<Output = T>,
    T: Copy,
{
    type Output = LinTermInX<T>;

    fn sub(self, other: LinSummand<T>) -> LinTermInX<T> {
        let mut d: Vec<LinSummand<T>> = self.summands;
        match other {
            LinSummand::Number(a) => d.push( LinSummand::Number(-a) ),
            LinSummand::Unknown(b) => d.push( LinSummand::Unknown(-b) ),
        }
        LinTermInX {summands: d}
    }
}

impl<T> PartialEq for LinTermInX<T>
where
    T: Num,
    T: Copy,
{
    fn eq(&self, other: &Self) -> bool {
        let mut m = T::zero();
        let mut n = T::zero();
        let mut o = T::zero();
        let mut p = T::zero();

        for (_, s) in self.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => m = m + *a,
                LinSummand::Unknown(b) => n = n + *b,
            }
        }

        for (_, s) in other.summands.iter().enumerate() {
            match s {
                LinSummand::Number(a) => o = o + *a,
                LinSummand::Unknown(b) => p = p + *b,
            }
        }

        return (m == o) && (n == p);
    }
}

fn pentagonal_nr<T>(n: T) -> Option<T>
where
    T: Num,
    T: Copy,
{
    let two = T::one() + T::one();
    let three = T::one() + T::one() + T::one();

    if two.is_zero() {
        return None;
    } else {
        return Some(n*(three*n-T::one())/two);
    }
}

fn sum_of_divisors<T>(i: T) -> T
where
    T: Integer,
    T: Copy,
{
    let mut s = T::zero();
    let mut n = T::zero() - i;

    while n <= i {
        let p_n = pentagonal_nr(n);

        match p_n {
            None => return T::one(),
            Some(p) => {
                let d = i - p;
                let sd: T;

                if d.is_zero() {
                    sd = i;
                } else if d > T::zero() && d < i {
                    sd = sum_of_divisors(d);
                } else {
                    sd = T::zero();
                }

                if n.is_odd() {
                    s = s + sd;
                } else {
                    s = s - sd;
                }
            },
        }

        n = n + T::one();
    }

    return s;
}

fn nr_of_partitions<T>(i: T) -> T
where
    T: Integer,
    T: Copy,
{
    let mut s = T::zero();
    let mut n = T::zero() - i;

    while n <= i {
        let p_n = pentagonal_nr(n);

        match p_n {
            None => return T::one(),
            Some(p) => {
                let d = i - p;
                let sd: T;

                if d.is_zero() {
                    sd = T::one();
                } else if d > T::zero() && d < i {
                    sd = nr_of_partitions(d);
                } else {
                    sd = T::zero();
                }

                if n.is_odd() {
                    s = s + sd;
                } else {
                    s = s - sd;
                }
            },
        }

        n = n + T::one();
    }

    return s;
}

fn factorial<T>(n: T) -> T
where
    T: Integer,
    T: Copy,
{
    let mut f = T::one();
    let mut p = T::one();
    while f <= n {
        p = p * f;
        f = f + T::one();
    }
    return p;
}

fn pascal_array<T>(m: T, k: T) -> T
where
    T: Integer,
    T: Copy,
{
    factorial(m+k)/(factorial(m)*factorial(k))
}

fn choose<T>(n: T, k: T) -> T
where
    T: Integer,
    T: Copy,
{
    factorial(n)/(factorial(k)*factorial(n-k))
}

#[derive(Debug, Clone)]
pub struct PolyNumber<T> {
    n: Vec<T>,
}

impl<T> PolyNumber<T>
where
    T: Num,
{
    fn lowest_order(&self) -> Option<usize> {

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
    T: Num,
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
    T: Num,
    T: Copy,
{
    type Output = PolyNumber<T>;

    fn add(self, other: Self) -> PolyNumber<T> {
        let mut index: usize = 0;
        let mut s: Vec<T> = Vec::new();

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push(*aa + *bb),
                (Some(aa), None)     => s.push(*aa),
                (None, Some(bb))     => s.push(*bb),
                (None, None)         => return PolyNumber { n: s },
            }

            index += 1;
        }
    }
}

impl<T> Sub for PolyNumber<T>
where
    T: Num,
    T: Copy,
{
    type Output = PolyNumber<T>;

    fn sub(self, other: Self) -> PolyNumber<T> {
        let mut index: usize = 0;
        let mut s: Vec<T> = Vec::new();

        loop {
            let a = self.n.get(index);
            let b = other.n.get(index);

            match (a, b) {
                (Some(aa), Some(bb)) => s.push(*aa - *bb),
                (Some(aa), None)     => s.push(*aa),
                (None, Some(bb))     => s.push(T::zero() - *bb),
                (None, None)         => return PolyNumber { n: s },
            }

            index += 1;
        }
    }
}

impl<T> Mul for PolyNumber<T>
where
    T: Num,
    T: Copy,
{
    type Output = PolyNumber<T>;

    fn mul(self, other: Self) -> PolyNumber<T> {
        let mut p: Vec<T> = Vec::new();

        for (ai, av) in self.n.iter().enumerate() {
            for (bi, bv) in other.n.iter().enumerate() {
                let ci = ai+bi;
                let c = p.get(ci);

                match c {
                    Some(cv) => p[ci] = *cv + *av * *bv,
                    None     => p.push(*av * *bv),
                }
            }
        }

        return PolyNumber { n: p };
    }
}

impl<T> Mul<T> for PolyNumber<T>
where
    T: Num,
    T: Copy,
{
    type Output = PolyNumber<T>;

    fn mul(self, other: T) -> PolyNumber<T> {
        let mut p: Vec<T> = Vec::new();

        for a in self.n {
            p.push(a * other);
        }

        return PolyNumber { n: p };
    }
}

impl<T> PolyNumber<T>
where
    T: Num,
    T: Copy,
{
    fn eval(&self, c: T) -> T {
        let mut ck = T::one();
        let mut pc = T::zero();

        for a in &self.n {
            pc = pc + *a * ck;
            ck = ck * c;
        }

        return pc;
    }
}

impl<T> PolyNumber<T>
where
    T: Num,
    T: Copy,
{
    fn zero() -> PolyNumber<T> {
        return PolyNumber { n: vec![T::zero()] };
    }
}

impl<T> Div for PolyNumber<T>
where
    T: Num,
    T: Copy,
    T: Debug,
{
    type Output = Option<PolyRatio<T>>;

    fn div(self, other: PolyNumber<T>) -> Option<PolyRatio<T>> {
        let zero : PolyNumber<T> = PolyNumber::zero();
        if other == zero {
            return None;
        } else {
            return Some( PolyRatio::new(self, other) );
        }
    }
}

/// Represents the ratio between two poly numbers
#[derive(Debug)]
pub struct PolyRatio<T> {
    numer: PolyNumber<T>,
    denom: PolyNumber<T>,
}

impl<T> PolyRatio<T>
where
    T: Num,
    T: Copy,
    T: Debug,
{
    fn new(numer: PolyNumber<T>, denom: PolyNumber<T>) -> PolyRatio<T> {
        let dlo = denom.lowest_order();
        let dlow = match dlo {
            Some(i) => i,
            None => panic!("Rational poly number has a denominator equal zero."),
        };
        let d1 = denom.n[dlow];

        let mut res = numer.clone();
        let mut q: Vec<T> = Vec::new();

        for _i in 0..(numer.n.len()-denom.n.len()+1) {
            let nlo = res.lowest_order();

            let nlow = match nlo {
                Some(i) => i,
                None => 0,
            };


            if nlow < dlow {
                return PolyRatio{numer: numer.clone(), denom: denom};
            }

            if nlow > 0 {
                for _j in q.clone().len()..(nlow -1) {
                    q.push(T::zero());
                }
            }

            let n1 = res.n[nlow];

            let q1 = n1/d1;

            if q1*d1 == n1 {
                q.push(q1);
                let pq1 = PolyNumber { n: q.clone() };
                let s1 = denom.clone() * PolyNumber { n: q.clone() };
                res = numer.clone() - s1;
                let zero : PolyNumber<T> =  PolyNumber::zero();

                if res == zero {
                    return PolyRatio{numer: pq1, denom: PolyNumber{n: vec![T::one()]}};
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
    T: Num,
    T: Copy,
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
    T: Num,
    T: Copy,
    T: Debug,
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
    T: Num,
    T: Copy,
    T: Debug,
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

impl<T> Sub for PolyRatio<T>
where
    T: Num,
    T: Copy,
    T: Debug,
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
    T: Num,
    T: Copy,
    T: Debug,
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

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};

    #[test]
    fn xplus2eq7() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(1), LinSummand::Number(2)],
            rhs: vec![LinSummand::Number(7)]};
        assert_eq!(eq1.x(), Some(5));
    }
    #[test]
    fn xplus10eq2() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(1), LinSummand::Number(10)],
            rhs: vec![LinSummand::Number(2)]};
        assert_eq!(eq1.x(), Some(-8));
    }
    #[test]
    fn threexeq12() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(3)],
            rhs: vec![LinSummand::Number(12)]};
        assert_eq!(eq1.x(), Some(4));
    }
    #[test]
    fn fivexeq11() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(Ratio::new(5, 1))],
            rhs: vec![LinSummand::Number(Ratio::new(11, 1))]};
        assert_eq!( eq1.x(), Some(Ratio::new(11, 5)) );
    }
    #[test]
    fn solve_lin_eq() {
        let eq1 = LinEq {lhs: vec![LinSummand::Unknown(3), LinSummand::Number(4)],
            rhs: vec![LinSummand::Number(11), LinSummand::Unknown(-2), LinSummand::Number(8)]};
        assert_eq!(eq1.x(), Some(3));
    }
    #[test]
    fn terms_in_x() {
        let lhs = LinTermInX {summands: vec![LinSummand::Unknown(Ratio::new(6, 1)),
            LinSummand::Number(Ratio::new(-1, 1))]} / Ratio::new(2, 1);
        let rhs = LinTermInX {summands: vec![LinSummand::Unknown(Ratio::new(3, 1)),
            LinSummand::Number(Ratio::new(5, 1))]} / Ratio::new(7, 1);
        let eq1 = LinEq {lhs: lhs.summands, rhs: rhs.summands};
        assert_eq!( eq1.x(), Some(Ratio::new(17, 36)) );
    }
    #[test]
    fn arithmatic_with_terms_in_x() {
        let t1 = LinTermInX {summands: vec![LinSummand::Unknown(1),
            LinSummand::Number(3)]};
        let t2 = t1 * 2;
        let t3 = t2 + 10;
        let t4 = t3 - 4;
        let t5 = t4 / 2;
        let t6 = t5 - LinSummand::Unknown(1);
        assert_eq!(t6, LinTermInX {summands: vec![LinSummand::Number(6)]});
    }
    #[test]
    fn pentagonal_numbers() {
        assert_eq!(pentagonal_nr(5),Some(35));
        assert_eq!(pentagonal_nr(6),Some(51));
    }
    #[test]
    fn sums_of_divisors() {
        assert_eq!(sum_of_divisors(1),1);
        assert_eq!(sum_of_divisors(2),3);
        assert_eq!(sum_of_divisors(3),4);
        assert_eq!(sum_of_divisors(4),7);
        assert_eq!(sum_of_divisors(5),6);
        assert_eq!(sum_of_divisors(6),12);
        assert_eq!(sum_of_divisors(12),28);
        assert_eq!(sum_of_divisors(14),24);
    }
    #[test]
    fn test_nr_of_partitions() {
        assert_eq!(nr_of_partitions(1),1);
        assert_eq!(nr_of_partitions(2),2);
        assert_eq!(nr_of_partitions(12),77);
        assert_eq!(nr_of_partitions(13),101);
        assert_eq!(nr_of_partitions(14),135);
    }
    #[test]
    fn factorials() {
        assert_eq!(factorial(1),1);
        assert_eq!(factorial(2),2);
        assert_eq!(factorial(4),24);
        assert_eq!(factorial(5),120);
        assert_eq!(factorial(6),720);
    }
    #[test]
    fn pascal_number() {
        assert_eq!(pascal_array(0,0),1);
        assert_eq!(pascal_array(1,1),2);
        assert_eq!(pascal_array(3,4),35);
        assert_eq!(pascal_array(3,5),56);
        assert_eq!(pascal_array(3,6),84);
        assert_eq!(pascal_array(4,5),126);
        assert_eq!(pascal_array(4,6),210);
        assert_eq!(choose(3,1),3);
        assert_eq!(choose(4,2),6);
    }
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
    fn rational_poly_numbers() {
        let p1 = PolyNumber { n: vec![1, 0, -2] };
        let p2 = PolyNumber { n: vec![2, -5] };
        let rp1 = p1/p2;
        let p3 = PolyNumber { n: vec![2, 0, -4] };
        let p4 = PolyNumber { n: vec![4, -10] };
        let rp2 = p3/p4;
        assert_eq!(rp1,rp2);
    }
    #[test]
    fn dividing_poly_numbers() {
        let p1 = PolyNumber { n: vec![2, 7, 2, -3] };
        let p2 = PolyNumber { n: vec![2, 1, -1] };
        let p3 = PolyNumber { n: vec![1, 3] };
        let p4 = PolyNumber { n: vec![1] };
        let rp1 = p1/p2;
        match rp1 {
            Some(PolyRatio { numer: n, denom: d }) => {
                assert_eq!(n,p3);
                assert_eq!(d,p4); },
            None => panic!("Legal division of poly numbers returns None."),
        }
        let p21 = PolyNumber { n: vec![12, 8, -7, -2, 1] };
        let p22 = PolyNumber { n: vec![4, 0, -1] };
        let p23 = PolyNumber { n: vec![3, 2, -1] };
        let p24 = PolyNumber { n: vec![1] };
        let rp21 = p21/p22;
        match rp21 {
            Some(PolyRatio { numer: n, denom: d }) => {
                assert_eq!(n,p23);
                assert_eq!(d,p24); },
            None => panic!("Legal division of poly numbers returns None."),
        }
    }
}
