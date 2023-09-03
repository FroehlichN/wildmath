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


use std::ops::Add;
use core::cmp::max;


/// Represents a number in max-plus algebra
#[derive(PartialEq, Debug, Clone)]
pub struct MaxPlus<T> {
    pub n : T,
}

impl<T> MaxPlus<T>
{
    pub fn new(n : T) -> MaxPlus<T> {
        MaxPlus{ n: n }
    }
}


impl<T> Add<MaxPlus<T>> for MaxPlus<T>
where
    T: Ord,
{
    type Output = MaxPlus<T>;

    fn add(self, other: MaxPlus<T>) -> MaxPlus<T> {
        MaxPlus {n: max(self.n, other.n) }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn addition() {
        let t1 = MaxPlus::new(5);
        let t2 = MaxPlus::new(7);
        assert_eq!(t1+t2.clone(), t2);
    }
}

