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


/// Represents the interval between two numbers
#[derive(Debug,Clone)]
pub struct Interval<T> {
    l: T,
    u: T,
}

impl<T> Interval<T>
where
    T: PartialOrd,
{
    fn new(m : T, n : T) -> Interval<T> {
        if m < n {
            return Interval{ l: m, u: n };
        } else {
            return Interval{ l: n, u: m };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn new_interval_is_ordered() {
        let m = 10;
        let n = 9;
        let i = Interval::new(m, n); 
        assert!(i.l < i.u);
    }
}

