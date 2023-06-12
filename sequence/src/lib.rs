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

use std::collections::BTreeSet;
use num::{Zero,One,Integer};
use num::rational::Ratio;


type Sequence<T> = BTreeSet<T>;


pub fn farey<T>(n : T) -> Sequence<Ratio<T>>
where
    T: Zero,
    T: One,
    T: Integer,
    T: PartialEq,
    T: Copy,
{
    if n == T::zero() {
        return Sequence::new();
    }

    let mut f = Sequence::from([Ratio::new(T::zero(),T::one()),
                                Ratio::new(T::one() ,T::one())]);

    let mut denom = T::one();

    while denom <= n {
        let mut numer = T::one();
        while numer < denom {
            f.insert(Ratio::new(numer,denom));
            numer = numer + T::one();
        }
        denom = denom + T::one();
    }
    return f;
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn farey_sequence() {
        let f1 = farey(1);
        let s1 = Sequence::from([Ratio::new(0,1), Ratio::new(1,1)]);
        assert_eq!(f1,s1);
        let f2 = farey(2);
        let s2 = Sequence::from([Ratio::new(0,1),Ratio::new(1,2),Ratio::new(1,1)]);
        assert_eq!(f2,s2);
        let f3 = farey(3);
        let s3 = Sequence::from([Ratio::new(0,1),Ratio::new(1,3),Ratio::new(1,2),
                                 Ratio::new(2,3),Ratio::new(1,1)]);
        assert_eq!(f3,s3);
        let f4 = farey(4);
        let s4 = Sequence::from([Ratio::new(0,1),Ratio::new(1,4),Ratio::new(1,3),
                    Ratio::new(1,2),Ratio::new(2,3),Ratio::new(3,4),Ratio::new(1,1)]);
        assert_eq!(f4,s4);
        let f5 = farey(5);
        let s5 = Sequence::from([Ratio::new(0,1),Ratio::new(1,5),Ratio::new(1,4),
                    Ratio::new(1,3),Ratio::new(2,5),Ratio::new(1,2),Ratio::new(3,5),
                    Ratio::new(2,3),Ratio::new(3,4),Ratio::new(4,5),Ratio::new(1,1)]);
        assert_eq!(f5,s5);
        let f6 = farey(6);
        let s6 = Sequence::from([Ratio::new(0,1),Ratio::new(1,6),Ratio::new(1,5),
                    Ratio::new(1,4),Ratio::new(1,3),Ratio::new(2,5),Ratio::new(1,2),
                    Ratio::new(3,5),Ratio::new(2,3),Ratio::new(3,4),Ratio::new(4,5),
                    Ratio::new(5,6),Ratio::new(1,1)]);
        assert_eq!(f6,s6);
    }
}

