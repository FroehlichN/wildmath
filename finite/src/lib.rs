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
use std::ops::{Add};
use paste::paste;

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

            impl Zero for [<Finite $N>] {
                fn zero() -> [<Finite $N>] {
                    [<Finite $N>] { n: 0 }
                }

                fn is_zero(&self) -> bool {
                    self.n.is_zero()
                }
            }

        }
    }
}




#[cfg(test)]
mod tests {
    use super::*;

    create_finite_field!(3);

    #[test]
    fn three_equals_zero() {
        let f = Finite3::new(3);
        assert!(f.is_zero());
    }
}

