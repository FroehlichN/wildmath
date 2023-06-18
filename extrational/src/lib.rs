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


use num::{Num, Integer};



/// RatInf
/// extended rational numbers
/// rational numbers with infinity
#[derive(Debug, Clone)]
pub struct RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    a: T,
    b: T,
}

impl<T> RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    pub fn new(a: T, b: T) -> RatInf<T> {
        RatInf {a: a, b: b}
    }
}

impl<T> PartialEq for RatInf<T>
where
    T: Num,
    T: Integer,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let lhs = self.a.clone() * other.b.clone();
        let rhs = self.b.clone() * other.a.clone();
        let s_non_zero = !self.a.is_zero() || !self.b.is_zero();
        let o_non_zero = !other.a.is_zero() || !other.b.is_zero();

        return (lhs == rhs) && s_non_zero && o_non_zero ;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn equality_of_extended_rational_numbers() {
        let a1 = RatInf::new(2,3);
        let a2 = RatInf::new(4,6);
        assert_eq!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers2() {
        let a1 = RatInf::new(3,0);
        let a2 = RatInf::new(-2,0);
        assert_eq!(a1, a2);
    }
    #[test]
    fn equality_of_extended_rational_numbers3() {
        let a1 = RatInf::new(0,7);
        let a2 = RatInf::new(0,-24);
        assert_eq!(a1, a2);
    }
    #[test]
    fn inequality_of_extended_rational_numbers1() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(2,3);
        assert!(!(a1==a2));
    }
    #[test]
    fn inequality_of_extended_rational_numbers2() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(5,0);
        assert!(!(a1==a2));
    }
    #[test]
    fn inequality_of_extended_rational_numbers3() {
        let a1 = RatInf::new(0,0);
        let a2 = RatInf::new(0,4);
        assert!(!(a1==a2));
    }
}

