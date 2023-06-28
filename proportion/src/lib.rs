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

use num::{Zero};
use std::ops::{Mul, Add};



/// Represents a 2-proportion
#[derive(Debug, Clone)]
pub struct TwoProportion<T> {
    pub a : T,
    pub b : T,
}

impl<T> TwoProportion<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(a: T, b: T) -> TwoProportion<T> {
        if a.is_zero() && b.is_zero() {
            panic!("Proportions cannot have all zero entries");
        }
        return TwoProportion{ a: a, b: b };
    }
}

impl<T> PartialEq for TwoProportion<T>
where
    T: PartialEq,
    T: Zero,
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
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
    use num::rational::{Ratio};
    #[test]
    #[should_panic]
    fn invalid_2proportion() {
        TwoProportion::new(0,0);
    }
    #[test]
    fn equality_of_2proportions() {
        let a1 = TwoProportion::new(3,4);
        let b1 = TwoProportion::new(6,8);
        assert_eq!(a1,b1);
        let a2 = TwoProportion::new(Ratio::new(1,2),Ratio::new(3,1));
        let b2 = TwoProportion::new(Ratio::new(3,2),Ratio::new(9,1));
        assert_eq!(a2,b2);
        let a3 = TwoProportion::new(Ratio::new(1,1),Ratio::new(0,1));
        let b3 = TwoProportion::new(Ratio::new(7,1),Ratio::new(0,1));
        let c3 = TwoProportion::new(Ratio::new(-7,4),Ratio::new(0,1));
        assert_eq!(a3,b3);
        assert_eq!(a3,c3);
        let a4 = TwoProportion::new(Ratio::new(0,1),Ratio::new(-5,1));
        let b4 = TwoProportion::new(Ratio::new(0,1),Ratio::new(17,19));
        assert_eq!(a4,b4);
    }
}

