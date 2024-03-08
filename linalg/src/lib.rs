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

use num::{Zero,One};
use std::ops::{Add,Sub,Mul,Div,Neg};


#[derive(Debug, Clone)]
pub struct RowVector<T> {
    elem: Vec<T>,
}

impl<T> RowVector<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(elem: Vec<T>) -> RowVector<T> {
        RowVector { elem: elem }
    }

    pub fn get(&self, ci: usize) -> T {
        let v = self.elem.get(ci);

        let value = match v {
            Some(vv) => vv.clone(),
            None => T::zero(),
        };
        value
    }
}

impl<T> PartialEq for RowVector<T>
where
    T: PartialEq,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let scols = self.elem.len();
        let ocols = other.elem.len();
        let cols = if scols > ocols {scols} else {ocols};

        for ci in 0..cols {
            if !(self.get(ci) == other.get(ci)) {
                return false;
            }
        }
        return true;
    }
}


impl<T> Mul<Matrix<T>> for RowVector<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = RowVector<T>;

    fn mul(self, other: Matrix<T>) -> RowVector<T> {
        let vcols = self.elem.len();
        let cols = if vcols < other.rows {vcols} else {other.rows};

        let mut elem : Vec<T> = Vec::new();
        for mci in 0..other.cols {
            let mut s = T::zero();
            for vci in 0..cols {
                s = s + self.get(vci) * other.get(vci,mci);
            }
            elem.push(s);
        }
        RowVector { elem: elem }
    }
}

impl<T> Mul<ColumnVector<T>> for RowVector<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = T;

    fn mul(self, other: ColumnVector<T>) -> T {
        let scols = self.elem.len();
        let orows = other.elem.len();
        let len = if scols < orows {scols} else {orows};

        let mut s = T::zero();
        for vi in 0..len {
            s = s + self.get(vi) * other.get(vi);
        }
        s
    }
}


#[derive(Debug, Clone)]
pub struct ColumnVector<T> {
    elem: Vec<T>,
}

impl<T> ColumnVector<T>
where
    T: Zero,
    T: Clone,
{
    pub fn new(elem: Vec<T>) -> ColumnVector<T> {
        ColumnVector { elem: elem }
    }

    pub fn get(&self, ri: usize) -> T {
        let v = self.elem.get(ri);

        let value = match v {
            Some(vv) => vv.clone(),
            None => T::zero(),
        };
        value
    }
}

impl<T> PartialEq for ColumnVector<T>
where
    T: PartialEq,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let srows = self.elem.len();
        let orows = other.elem.len();
        let rows = if srows > orows {srows} else {orows};

        for ri in 0..rows {
            if !(self.get(ri) == other.get(ri)) {
                return false;
            }
        }
        return true;
    }
}


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
    pub fn new(elem: Vec<Vec<T>>) -> Matrix<T> {
        let rows = elem.len();
        let mut cols : usize = 0;
        for (_, rv) in elem.iter().enumerate() {
            if rv.len() >  cols {
                cols = rv.len();
            }
        }
        Matrix { rows: rows, cols: cols, elem: elem }
    }

    pub fn get(&self, ri: usize, ci: usize) -> T {
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

impl<T> Matrix<T>
where
    T: Zero + One,
    T: Div<Output = T>,
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn inverse(&self) -> Option<Matrix<T>> {
        if self.rows == 1 && self.cols == 1 {
            let a = self.get(0,0);
            if a.is_zero() {
                return None;
            } else {
                let b = T::one()/a;
                return Some(Matrix::new(vec![vec![b]]));
            }
        } else if self.rows == 2 && self.cols == 2 {
            let a = self.get(0,0);
            let b = self.get(0,1);
            let c = self.get(0,1);
            let d = self.get(1,1);

            let det = a.clone()*d.clone() - b.clone()*c.clone();
            if det.is_zero() {
                return None;
            } else {
                let ai = d/det.clone();
                let bi = -b/det.clone();
                let ci = -c/det.clone();
                let di = a/det.clone();
                let mi = Matrix::new(vec![vec![ai,bi],vec![ci,di]]);
                return Some(mi);
            }
        }
        return None;
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

impl<T> Mul<Matrix<T>> for Matrix<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Matrix<T>;

    fn mul(self, other: Matrix<T>) -> Matrix<T> {
        let rows = self.rows;
        let cols = other.cols;
        let scount = if self.cols < other.rows {self.cols} else {other.rows};

        let mut elem : Vec<Vec<T>> = Vec::new();
        for mri in 0..rows {
            let mut row : Vec<T> = Vec::new();
            for mci in 0..cols {
                let mut s = T::zero();
                for si in 0..scount {
                    s = s + self.get(mri,si) * other.get(si,mci);
                }
                row.push(s)
            }
            elem.push(row);
        }
        Matrix { rows: rows, cols: cols, elem: elem }
    }
}

impl<T> Mul<ColumnVector<T>> for Matrix<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = ColumnVector<T>;

    fn mul(self, other: ColumnVector<T>) -> ColumnVector<T> {
        let vrows = self.elem.len();
        let rows = if vrows < self.cols {vrows} else {self.cols};

        let mut elem : Vec<T> = Vec::new();
        for mri in 0..self.rows {
            let mut s = T::zero();
            for vri in 0..rows {
                s = s + self.get(mri,vri) * other.get(vri);
            }
            elem.push(s);
        }
        ColumnVector { elem: elem }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};
    #[test]
    fn matrix_addition() {
        let m1 = Matrix::new(vec![vec![-1, 0], vec![0, 1]]);
        let m2 = Matrix::new(vec![vec![1,  1], vec![0, 1]]);
        let m3 = Matrix::new(vec![vec![0,  1], vec![0, 2]]);
        assert_eq!(m1+m2, m3);
    }
    #[test]
    fn matrix_multiplication() {
        let m1 = Matrix::new(vec![vec![1, 2], vec![3, 4]]);
        let m2 = Matrix::new(vec![vec![5,  6], vec![7, 8]]);
        let m3 = Matrix::new(vec![vec![19, 22], vec![43, 50]]);
        assert_eq!(m1*m2, m3);
    }
    #[test]
    fn column_vector_matrix_multiplication() {
        let v1 = RowVector::new(vec![1, 2]);
        let m1 = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        let v2 = RowVector::new(vec![4, 7]);
        assert_eq!(v1*m1,v2);
    }
    #[test]
    fn matrix_row_vector_multiplication() {
        let m1 = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        let v1 = ColumnVector::new(vec![1, 2]);
        let v2 = ColumnVector::new(vec![2, 8]);
        assert_eq!(m1*v1,v2);
    }
    #[test]
    fn inverse_of_1d_matrix() {
        let m1 = Matrix::new(vec![vec![0]]);
        assert_eq!(m1.inverse(),None);
        let m2 = Matrix::new(vec![vec![Ratio::new(2,3)]]);
        let m3 = Matrix::new(vec![vec![Ratio::new(3,2)]]);
        assert_eq!(m2.inverse(),Some(m3));
    }
    #[test]
    fn inverse_of_2d_matrix() {
        let m1 = Matrix::new(vec![vec![Ratio::new(3,1),Ratio::new(1,1)],
                                  vec![Ratio::new(1,1),Ratio::new(2,1)]]);
        let m2 = Matrix::new(vec![vec![Ratio::new(2,5),Ratio::new(-1,5)],
                                  vec![Ratio::new(-1,5),Ratio::new(3,5)]]);
        assert_eq!(m1.inverse(),Some(m2));
    }
}

