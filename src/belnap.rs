/*
Copyright 2025 Norbert Fröhlich


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


use std::ops::{Not,BitOr,BitAnd};


/// Implementation of Belnap's logic with four possible values
#[derive(Debug,Clone,Copy,PartialEq)]
enum Belnap {
    True,
    False,
    Both,
    Neither,
}

impl Not for Belnap {
    type Output = Belnap;

    fn not(self) -> Belnap {
        match self {
            Belnap::True => Belnap::False,
            Belnap::False => Belnap::True,
            Belnap::Both => Belnap::Both,
            Belnap::Neither => Belnap::Neither,
        }
    }
}

impl BitOr for Belnap {
    type Output = Belnap;

    fn bitor(self,other: Self) -> Belnap {
        match (self,other) {
            (Belnap::True,_) => Belnap::True,
            (_,Belnap::True) => Belnap::True,
            (Belnap::False,_) => other,
            (_,Belnap::False) => self,
            (Belnap::Both,Belnap::Both) => Belnap::Both,
            (Belnap::Both,Belnap::Neither) => Belnap::True,
            (Belnap::Neither,Belnap::Neither) => Belnap::Neither,
            (Belnap::Neither,Belnap::Both) => Belnap::True,
        }
    }
}

impl BitAnd for Belnap {
    type Output = Belnap;

    fn bitand(self,other: Self) -> Belnap {
        match (self,other) {
            (Belnap::True,_) => other,
            (_,Belnap::True) => self,
            (Belnap::False,_) => Belnap::False,
            (_,Belnap::False) => Belnap::False,
            (Belnap::Both,Belnap::Both) => Belnap::Both,
            (Belnap::Both,Belnap::Neither) => Belnap::False,
            (Belnap::Neither,Belnap::Neither) => Belnap::Neither,
            (Belnap::Neither,Belnap::Both) => Belnap::False,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn belnap_logic_not() {
        assert_eq!(!Belnap::True, Belnap::False);
        assert_eq!(!Belnap::False, Belnap::True);
        assert_eq!(!Belnap::Both, Belnap::Both);
        assert_eq!(!Belnap::Neither, Belnap::Neither);
    }
    #[test]
    fn belnap_logic_and_or() {
        let v = [Belnap::True, Belnap::False, Belnap::Both, Belnap::Neither];

        for a in v {
            for b in v {
                if a & b == a {
                    assert_eq!(a | b, b);
                }
                if a | b == b {
                    assert_eq!(a & b, a);
                }
            }
        }
    }
    #[test]
    fn belnap_de_morgan() {
        let v = [Belnap::True, Belnap::False, Belnap::Both, Belnap::Neither];

        for a in v {
            for b in v {
                assert_eq!(!(a & b), !a | !b);
                assert_eq!(!(a | b), !a & !b);
            }
        }
    }
}

