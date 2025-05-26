/*
Copyright 2023 - 2025 Norbert Fröhlich


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



/// Represents a proportion
#[derive(Debug, Clone)]
pub struct Proportion<T> {
    p : Vec<T>,
}

impl<T> Proportion<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(elem: Vec<T>) -> Proportion<T> {
        let mut all_zero = true;
        for (_, e) in elem.iter().enumerate() {
            if !e.is_zero() {
                all_zero = false;
            }
        }
        if all_zero {
            panic!("Proportions cannot have all zero entries");
        }
        return Proportion{ p: elem };
    }
}

impl<T> PartialEq for Proportion<T>
where
    T: PartialEq,
    T: Zero,
    T: Mul<T, Output = T>,
    T: Add<T, Output = T>,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let sl = self.p.len();
        let ol = other.p.len();
        if sl != ol {
            return false;
        }

        for i in 0..(sl-1) {
            let lhs = self.p[i].clone() * other.p[i+1].clone();
            let rhs = self.p[i+1].clone() * other.p[i].clone();

            if lhs != rhs {
                return false;
            }
        }
        return true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
    #[test]
    #[should_panic]
    fn invalid_2proportion() {
        Proportion::new(vec![0,0]);
    }
    #[test]
    fn equality_of_2proportions() {
        let a1 = Proportion::new(vec![3,4]);
        let b1 = Proportion::new(vec![6,8]);
        assert_eq!(a1,b1);
        let a2 = Proportion::new(vec![Ratio::new(1,2),Ratio::new(3,1)]);
        let b2 = Proportion::new(vec![Ratio::new(3,2),Ratio::new(9,1)]);
        assert_eq!(a2,b2);
        let a3 = Proportion::new(vec![Ratio::new(1,1),Ratio::new(0,1)]);
        let b3 = Proportion::new(vec![Ratio::new(7,1),Ratio::new(0,1)]);
        let c3 = Proportion::new(vec![Ratio::new(-7,4),Ratio::new(0,1)]);
        assert_eq!(a3,b3);
        assert_eq!(a3,c3);
        let a4 = Proportion::new(vec![Ratio::new(0,1),Ratio::new(-5,1)]);
        let b4 = Proportion::new(vec![Ratio::new(0,1),Ratio::new(17,19)]);
        assert_eq!(a4,b4);
    }
}

