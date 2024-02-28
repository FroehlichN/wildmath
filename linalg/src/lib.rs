/*
Copyright 2024 Norbert Fröhlich


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


#[derive(Debug, Clone)]
pub struct Matrix<T> {
    rows: usize,
    cols: usize,
    elem: Vec<Vec<T>>,
}

impl<T> Matrix<T>
where
    T: Zero,
    T: Clone,
{
    fn get(&self, ri: usize, ci: usize) -> T {
        let r = self.elem.get(ri);

        let empty_row = Vec::new();
        let row = match r {
            Some(rv) => rv,
            None => &empty_row,
        };

        let v = row.get(ci);

        let value = match v {
            Some(vv) => vv.clone(),
            None => T::zero(),
        };
        value
    }
}

impl<T> PartialEq for Matrix<T>
where
    T: PartialEq,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let rows = if self.rows > other.rows {self.rows} else {other.rows};
        let cols = if self.cols > other.cols {self.cols} else {other.cols};

        for ri in 0..rows {
            for ci in 0..cols {
                if !(self.get(ri,ci) == other.get(ri,ci)) {
                    return false;
                }
            }
        }
        return true;
    }
}


impl<T> Add for Matrix<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Clone,
{
    type Output = Matrix<T>;

    fn add(self, other: Self) -> Matrix<T> {
        let rows = if self.rows > other.rows {self.rows} else {other.rows};
        let cols = if self.cols > other.cols {self.cols} else {other.cols};

        let mut elem : Vec<Vec<T>> = Vec::new();
        for ri in 0..rows {
            let mut row : Vec<T> = Vec::new();
            for ci in 0..cols {
                row.push(self.get(ri,ci) + other.get(ri,ci));
            }
            elem.push(row);
        }
        Matrix { rows: rows, cols: cols, elem: elem }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn matrix_addition() {
        let m1 = Matrix { rows: 2, cols: 2, elem: vec![vec![-1, 0], vec![0, 1]] };
        let m2 = Matrix { rows: 2, cols: 2, elem: vec![vec![1,  1], vec![0, 1]] };
        let m3 = Matrix { rows: 2, cols: 2, elem: vec![vec![0,  1], vec![0, 2]] };
        assert_eq!(m1+m2, m3);
    }
}

