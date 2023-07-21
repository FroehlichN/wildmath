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


use proportion::TwoProportion;
use num::{Zero};
use std::ops::{Mul, Add, Sub, Div};


/// Represents a projective 1D point
#[derive(Debug, Clone)]
pub struct ProjOnePoint<T> {
    pub x : TwoProportion<T>,
}


impl<T> ProjOnePoint<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(x: T, y: T) -> ProjOnePoint<T> {
        ProjOnePoint { x: TwoProportion::new(x,y) }
    }
}


impl<T> ProjOnePoint<T>
where
    T: Add<T, Output = T>,
    T: Sub<T, Output = T>,
    T: Mul<T, Output = T>,
    T: Div<T, Output = T>,
    T: Clone,
{
    pub fn quadrance(&self, other: &Self) -> T {
        let x1 = self.x.a.clone();
        let y1 = self.x.b.clone();
        let x2 = other.x.a.clone();
        let y2 = other.x.b.clone();

        let n = x1.clone() * y2.clone() - x2.clone() * y1.clone();
        let numer = n.clone() * n;
        let d1 = x1.clone() * x1 + y1.clone() * y1;
        let d2 = x2.clone() * x2 + y2.clone() * y2;
        return numer / (d1*d2);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::Ratio;
    #[test]
    fn one_dimensional_euclidean_projective_quadrance() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(1,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(2,1));
        assert_eq!(a1.quadrance(&a2),Ratio::new(1,2));
    }
    #[test]
    fn one_dimensional_euclidean_projective_quadrance_2() {
        let a1 = ProjOnePoint::new(Ratio::new(3,1),Ratio::new(2,1));
        let a2 = ProjOnePoint::new(Ratio::new(1,1),Ratio::new(4,1));
        let q = Ratio::new(1,1)-Ratio::new(11*11,1)/Ratio::new(13*17,1);
        assert_eq!(a1.quadrance(&a2),q);
    }
}

