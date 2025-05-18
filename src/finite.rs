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


pub use num::{Zero,One};
pub use std::ops::{Add,Sub,Mul,Div};
pub use paste::paste;


/// Creates a finite number field
/// The argument N should be some prime number.
/// The resulting arithmetic sets N = 0.
#[macro_export]
macro_rules! create_finite_field {
    ( $N:literal ) => {

        paste! {

            #[derive(Debug, Clone, Copy)]
            pub struct [<Finite $N>] {
                n: i64,
            }

            impl [<Finite $N>] {
                pub fn new(n: i64) -> [<Finite $N>] {
                    let mut m = n % $N;
                    if m < 0 {
                        m += $N;
                    }
                    [<Finite $N>] { n: m }
                }
            }

            impl Add<[<Finite $N>]> for [<Finite $N>]
            {
                type Output = [<Finite $N>];

                fn add(self, other: [<Finite $N>]) -> [<Finite $N>] {
                    Self::new(self.n.clone() + other.n.clone())
                }
            }

            impl Sub<[<Finite $N>]> for [<Finite $N>]
            {
                type Output = [<Finite $N>];

                fn sub(self, other: [<Finite $N>]) -> [<Finite $N>] {
                    Self::new(self.n.clone() - other.n.clone())
                }
            }

            impl Mul<[<Finite $N>]> for [<Finite $N>]
            {
                type Output = [<Finite $N>];

                fn mul(self, other: [<Finite $N>]) -> [<Finite $N>] {
                    Self::new(self.n.clone() * other.n.clone())
                }
            }

            impl Div<[<Finite $N>]> for [<Finite $N>]
            {
                type Output = [<Finite $N>];

                fn div(self, other: [<Finite $N>]) -> [<Finite $N>] {

                    if other.n.is_zero() {
                        panic!("Division by zero not allowed.");
                    }

                    let mut m : i64 = 1;
                    while Self::new(other.n * m) != Self::one() {
                        m += 1;
                    }
                    Self::new(self.n * m)
                }
            }

            impl Zero for [<Finite $N>] {
                fn zero() -> [<Finite $N>] {
                    [<Finite $N>] { n: 0 }
                }

                fn is_zero(&self) -> bool {
                    self.n.is_zero()
                }
            }

            impl One for [<Finite $N>] {
                fn one() -> [<Finite $N>] {
                    [<Finite $N>] { n: 1 }
                }
            }

            impl PartialEq for [<Finite $N>]
            {
                fn eq(&self, other: &Self) -> bool {
                    self.n == other.n
                }
            }

        }
    }
}




#[cfg(test)]
mod tests {
    use super::*;

    create_finite_field!(3);
    create_finite_field!(7);

    #[test]
    fn three_equals_zero() {
        let f = Finite3::new(3);
        assert!(f.is_zero());
    }
    #[test]
    fn division_in_f7() {
        let n1 = Finite7::new(12);
        let d1 = Finite7::new(17);
        let n2 = Finite7::new(5);
        let d2 = Finite7::new(3);
        assert_eq!(n1/d1,n2/d2);
    }
    #[test]
    fn division_in_f7_2() {
        let n1 = Finite7::new(356);
        let d1 = Finite7::new(95);
        let n2 = Finite7::new(3);
        let d2 = Finite7::new(2);
        assert_eq!(n1/d1,n2/d2);
    }
    #[test]
    fn division_in_f7_3() {
        let n1 = Finite7::new(1);
        let d1 = Finite7::new(2);
        let n2 = Finite7::new(4);
        assert_eq!(n1/d1,n2);
    }
    #[test]
    fn division_in_f7_4() {
        let n1 = Finite7::new(1);
        let d1 = Finite7::new(3);
        let n2 = Finite7::new(5);
        assert_eq!(n1/d1,n2);
    }
    #[test]
    fn division_in_f7_5() {
        let n1 = Finite7::new(1);
        let d1 = Finite7::new(6);
        let n2 = Finite7::new(6);
        assert_eq!(n1/d1,n2);
    }
    #[test]
    fn arithmetic_in_f7() {
        let n1 = Finite7::new(32);
        let d1 = Finite7::new(19);
        let n2 = Finite7::new(43);
        let d2 = Finite7::new(101);
        let n3 = Finite7::new(22);
        let d3 = Finite7::new(5);
        let n4 = Finite7::new(47);
        let d4 = Finite7::new(11);
        let n6 = Finite7::new(6);
        assert_eq!(n1/d1*n2/d2+n3/d3*n4/d4,n6);
    }
}

