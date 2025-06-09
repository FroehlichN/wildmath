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

use num::{Zero, One};
use std::ops::{Mul, Add, Sub, Div, Neg};
use crate::matrix::{ColumnVector, RowVector, Matrix};


/// Represents some geometric object consisting of a single point and
/// none or multiple vectors.
/// A point without a vector is just a point.
/// A point with a vector is a line.
/// A point with two vectors is a plane.
/// A point with three vectors is a volume.
#[derive(Debug, Clone)]
pub struct GeoObj<T>
{
    point: Point<T>,
    vectors: Vec<Vec<T>>,
}

impl<T> GeoObj<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn new(point: Point<T>, vectors: Vec<Vec<T>>) -> GeoObj<T> {
        if vectors.len() < 2 {
            return GeoObj{ point: point, vectors: vectors };
        }

        // If more than one vector is given, check all of them are linear
        // independent. Skip the ones, that are.
        let mut m = Vec::new();
        for (_, v) in vectors.iter().enumerate() {
            m.push(v.clone());
            let matrix = Matrix::new(m.clone());
            if !matrix.has_full_rank() {
                m.pop();
            }
        }
        GeoObj{ point: point, vectors: m }
    }

    pub fn from(points: Vec<Point<T>>) -> GeoObj<T> {
        let point = points[0].clone();

        let mut vectors = Vec::new();
        for i in 1..points.len() {
            let point1 = points[i].clone();
            let v = Vector::new(point.clone(),point1);
            let cv = v.columnvector();
            vectors.push(cv.elem.clone());
        }
        GeoObj::new(point, vectors)
    }
}

impl<T> GeoObj<T>
where
    T: Zero,
    T: PartialEq,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn contains(&self, other: &Self) -> bool {
        if self.vectors.len().is_zero() {
            if other.vectors.len().is_zero() {
                return self.point == other.point;
            } else {
                return false;
            }
        }

        // Lines cannot be contained in points and planes cannot be contained
        // in lines and volumes cannot be contained in planes etc. so ...
        if self.vectors.len() < other.vectors.len() {
            return false;
        }

        // If the vector from self.point to other.point can be described as a
        // linear combination of the self.vectors, then other.point lies in
        // self. This is the case exactly if the matrix formed by the
        // coordinates of these vectors has not full rank.
        let mut m = Vec::new();
        for (_, sv) in self.vectors.iter().enumerate() {
            m.push(sv.clone());
        }

        let spopv = Vector::new(self.point.clone(),other.point.clone());
        let mut pm = m.clone();
        pm.push(spopv.rowvector().elem);
        let pmatrix = Matrix::new(pm);
        if pmatrix.has_full_rank() {
            return false;
        }

        // Additionally, if every other.vector can be described as a linear
        // combination of self.vectors, then other is contained in self.
        for (_, ov) in other.vectors.iter().enumerate() {
            let mut vm = m.clone();
            vm.push(ov.clone());
            let vmatrix = Matrix::new(vm);
            if vmatrix.has_full_rank() {
                return false;
            }
        }
        return true;
    }
}

impl<T> PartialEq for GeoObj<T>
where
    T: Zero,
    T: PartialEq,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        if self.vectors.len() == other.vectors.len() {
            return self.contains(&other);
        } else {
            return false;
        }
    }
}

impl<T> GeoObj<T>
where
    T: Zero + One,
    T: PartialEq,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn meet(&self, other: &Self) -> Option<GeoObj<T>> {

        if self.contains(other) {
            return Some(other.clone());
        }

        if other.contains(self) {
            return Some(self.clone());
        }

        // From the vectors of self and other form the columns of a matrix.
        // If the matrix is not solvable, remomve rows, until the matrix
        // becomes solvable and put those rows in an auxilliary matrix.
        // The found solution also hast to be consistent with the auxilliary
        // matrix. If it is not, return None.
        // If no solvable matrix can be formed, return None.
        // self:  p0 + pv0*mu0 + pv1*mu1 + pv2*mu2 + ...
        // other: q0 + qv0*nu0 + qv1*nu1 + qv2*nu2 + ...
        // matrix: [pv0 pv1 pv2 ... -qv0 -qv1 -qv2 ...]
        // vector of unknowns: (mu0 mu1 mu2 ... nu0 nu1 nu2 ...)^T
        // result vector: (q0-p0)
        let mut m1 = Vec::new();
        for (_, v) in self.vectors.iter().enumerate() {
            m1.push(v.clone());
        }
        for (_, v) in other.vectors.iter().enumerate() {
            let cv = -ColumnVector::new(v.clone());
            m1.push(cv.elem.clone());
        }
        let matrix1 = Matrix::new(m1).transpose();

        let cvo = ColumnVector::new(other.point.coords.clone());
        let cvs = ColumnVector::new(self.point.coords.clone());
        let res1 = cvo - cvs;

        // Check if parallel
        if !matrix1.has_full_rank() {
            return None;
        }

        // If we have more rows than needed, remove them.
        let mut m2 = matrix1.elem.clone();
        let mut r2 = res1.elem.clone();

        if matrix1.rows > matrix1.cols {
            for _n in 0..(matrix1.rows - matrix1.cols) {
                m2.pop();
                r2.pop();
            }
        }
        let matrix2 = Matrix::new(m2);
        let res2 = ColumnVector::new(r2.clone());

        // Fill up missing rows. Every missing row corresponds with a vector in
        // the resulting meet.
        let nr_of_vectors = matrix1.cols - matrix1.rows;

        let mut m3 = matrix2.elem.clone();
        let mut matrix3 = Matrix::new(m3.clone());
        let mut r3 = res2.elem.clone();
        for _n in 0..nr_of_vectors {
            let mut a = T::one();
            let mut row = Vec::new();
            r3.push(T::zero());
            for _c in 0..matrix2.cols {
                row.push(a.clone());
                a = a.clone() + T::one();
            }
            loop {
                m3.push(row.clone());
                matrix3 = Matrix::new(m3.clone());
                if matrix3.has_full_rank() {
                    break;
                }
                row.pop();
                a = a.clone() + T::one();
                row.push(a.clone());
            }
        }
        let res3 = ColumnVector::new(r3);

        let inv = matrix3.inverse().unwrap();
        let sol3 = inv.clone()*res3.clone();

        let residual = matrix1*sol3.clone() - res1;
        if !residual.is_zero() {
            // No consistent solution could be found. There is not meet.
            return None;
        }

        let mut l0 = ColumnVector::new(self.point.coords.clone());
        for (i, v) in self.vectors.iter().enumerate() {
            let mu = sol3.get(i);
            let cv = ColumnVector::new(v.clone());
            l0 = l0 + cv*mu;
        }
        let point = Point::new(l0.elem);

        let mut vectors = Vec::new();
        for n in 0..nr_of_vectors {
            // create more points to create vectors in the line / plane /
            //volume ... or, in general, meet.
            let mut d1 = r2.clone();
            for r in 0..nr_of_vectors {
                if r == n {
                    d1.push(T::one());
                } else {
                    d1.push(T::zero());
                }
            }
            let d1vector = ColumnVector::new(d1);
            let sol4 = inv.clone()*d1vector;

            let mut l1 = ColumnVector::new(self.point.coords.clone());
            for (i, v) in self.vectors.iter().enumerate() {
                let mu = sol4.get(i);
                let cv = ColumnVector::new(v.clone());
                l1 = l1 + cv*mu;
            }
            let point1 = Point::new(l1.elem);

            let v1 = Vector::new(point.clone(),point1.clone());
            vectors.push(v1.columnvector().elem);
        }

        return Some(GeoObj{ point: point, vectors: vectors });
    }
}


/// Represents an n-D point
#[derive(Debug, Clone)]
pub struct Point<T> {
    coords: Vec<T>,
}

impl<T> Point<T>
{
    pub fn new(coords: Vec<T>) -> Point<T> {
        Point { coords: coords }
    }
}

impl<T> Point<T>
where
    T: Zero,
    T: PartialEq,
    T: Add,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn lies_on(&self, pi: &Slice<T>) -> bool {
        let pointv = RowVector::new(self.coords.clone());
        let planev = ColumnVector::new(pi.coords.clone());
        let r = pointv * planev;
        r == pi.d.clone()
    }
}

impl<T> Point<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn quadrance(&self, other: &Self) -> T {
        let v = Vector::new((*self).clone(), (*other).clone());
        v.quadrance()
    }
}

impl<T> Add<Vector<T>> for Point<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Clone,
{
    type Output = Point<T>;

    fn add(self, other: Vector<T>) -> Point<T> {
        let selfvec = ColumnVector::new(self.coords);
        let sumvec = selfvec + other.columnvector();
        Point::new(sumvec.elem)
    }
}

impl<T> PartialEq for Point<T>
where
    T: Zero,
    T: PartialEq,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let svec = ColumnVector::new(self.coords.clone());
        let ovec = ColumnVector::new(other.coords.clone());
        svec == ovec
    }
}


/// Represents an n-D vector with a metric defined by its dot product
#[derive(Debug, Clone)]
pub struct Vector<T> {
    start: Point<T>,
    end: Point<T>,
    dot: fn(Vector<T>, other: Vector<T>) -> T,
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn new(start: Point<T>, end: Point<T>) -> Vector<T> {
        Vector {start: start, end: end, dot: dot_blue}
    }

    pub fn new_red(start: Point<T>, end: Point<T>) -> Vector<T> {
        Vector {start: start, end: end, dot: dot_red}
    }
}


fn dot_blue<T>(lhs: Vector<T>, rhs: Vector<T>) -> T
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    let v1 = lhs.rowvector();
    let v2 = rhs.columnvector();

    v1*v2
}

fn dot_red<T>(lhs: Vector<T>, rhs: Vector<T>) -> T
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    let v1 = lhs.rowvector();
    let v2 = rhs.columnvector();

    let t = v1.get(0)*v2.get(0);
    let t2 = t.clone() + t;
    t2 - v1*v2
}


impl<T> Mul for Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = T;

    fn mul(self, rhs: Self) -> T {
        (self.dot)(self, rhs)
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn from(coords: Vec<T>) -> Vector<T> {
        let mut zero = Vec::new();
        for (_, _) in coords.iter().enumerate() {
            zero.push(T::zero());
        }
        let start = Point::new(zero);
        let end = Point::new(coords);
        Vector::new(start, end)
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Clone,
{
    pub fn columnvector(&self) -> ColumnVector<T> {
        let selfs = ColumnVector::new(self.start.coords.clone());
        let selfe = ColumnVector::new(self.end.coords.clone());
        selfe - selfs
    }
    pub fn rowvector(&self) -> RowVector<T> {
        let selfs = RowVector::new(self.start.coords.clone());
        let selfe = RowVector::new(self.end.coords.clone());
        selfe - selfs
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn is_perpendicular(&self, other: &Self) -> bool {
        (self.clone() * other.clone()).is_zero()
    }

    pub fn quadrance(&self) -> T {
        (self.dot)(self.clone(), self.clone())
    }
}

impl<T> Vector<T>
where
    T: Zero + One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn spread(&self, other: &Self) -> T {
        let v1 = self.rowvector();
        let v2 = other.columnvector();
        let n = v1*v2;
        let n2 = n.clone()*n;
        T::one() - n2 / (self.quadrance() * other.quadrance())
    }
}

impl<T> Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn projection(&self, other: &Self) -> Vector<T> {
        let vw = self.clone() * other.clone();
        let vv = self.clone() * self.clone();
        self.clone() * (vw/vv)
    }
}

impl<T> Vector<T>
where
    T: Zero + One,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn reflection_in_normal(&self, other: &Self) -> Vector<T> {
        let vw = self.clone() * other.clone();
        let vv = self.clone() * self.clone();
        let two = T::one()+T::one();
        other.clone() - self.clone() * (two*vw/vv)
    }

    pub fn reflection_in_vector(&self, other: &Self) -> Vector<T> {
        let vw = self.clone() * other.clone();
        let vv = self.clone() * self.clone();
        let two = T::one()+T::one();
        self.clone() * (two*vw/vv) - other.clone()
    }
}

impl<T> PartialEq for Vector<T>
where
    T: PartialEq,
    T: Sub<Output = T>,
    T: Zero,
    T: Clone,
{
    fn eq(&self, other: &Self) -> bool {
        let s = self.columnvector();
        let o = other.columnvector();
        s == o
    }
}

impl<T> Add for Vector<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn add(self, other: Self) -> Vector<T> {
        let selfendvec = ColumnVector::new(self.end.coords);
        let sumendvec = selfendvec + other.columnvector();
        let sumend = Point::new(sumendvec.elem);
        Vector::new(self.start, sumend)
    }
}

impl<T> Sub for Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn sub(self, other: Self) -> Vector<T> {
        let selfendvec = ColumnVector::new(self.end.coords);
        let subendvec = selfendvec - other.columnvector();
        let subend = Point::new(subendvec.elem);
        Vector::new(self.start, subend)
    }
}

impl<T> Mul<T> for Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn mul(self, scalar: T) -> Vector<T> {
        let vec = self.columnvector();
        let startvec = ColumnVector::new(self.start.coords.clone());
        let endvec = startvec + vec * scalar;
        let end = Point::new(endvec.elem);
        Vector::new(self.start, end)
    }
}

impl<T> Neg for Vector<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn neg(self) -> Vector<T> {
        let mone = -T::one();
        self*mone
    }
}


/// Represents the slicing of n-D
/// a1*x1 + a2*x2 + a3*x3 + ... = d
/// A slice in 1D is just a point, basically a single number.
/// A slice in 2D is a line.
/// A slice in 3D is a plane.
/// A slice in 4D is a volume.
#[derive(Debug, Clone)]
pub struct Slice<T> {
    coords: Vec<T>,
    d: T,
}

impl<T> Slice<T>
where
    T: Zero,
{
    pub fn new(coords: Vec<T>, d: T) -> Slice<T> {
        for (_, v) in coords.iter().enumerate() {
            if !v.is_zero() {
                return Slice { coords: coords, d: d };
            }
        }
        panic!("Slice has proportion of all zeros.");
    }
}

impl<T> Slice<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn from(point: Point<T>, normal: Vector<T>) -> Slice<T> {
        let pv = Vector::from(point.coords.clone());
        let d = normal.clone() * pv;
        let ncv = normal.columnvector();
        Slice::new(ncv.elem, d)
    }
}

impl<T> Slice<T>
where
    T: Zero + One,
    T: PartialEq,
    T: Neg<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn join(points: Vec<Point<T>>) -> Option<Slice<T>> {
        let mut m = Vec::new();
        for (_, v) in points.iter().enumerate() {
            m.push(v.coords.clone());
        }
        let matrix = Matrix::new(m);
        let det = matrix.determinant();

        if det.is_zero() {
            return None;
        }

        let adj = matrix.adjugate();

        let mut d = Vec::new();
        for _n in 0..(points.len()) {
            d.push(T::one());
        }
        let dcv = ColumnVector::new(d);
        let slice = adj*dcv;
        Some(Slice::new(slice.elem,det))
    }
}

impl<T> Slice<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn is_parallel(&self, other: &Self) -> bool {
        let mut m : Vec<Vec<T>> = Vec::new();
        m.push(self.coords.clone());
        m.push(other.coords.clone());
        let matrix = Matrix::new(m);
        if matrix.has_full_rank() {
            return false;
        } else {
            return true;
        }
    }
}

impl<T> Slice<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn is_perpendicular(self, other: &Self) -> bool {
        let sv = RowVector::new(self.coords.clone());
        let ov = ColumnVector::new(other.coords.clone());
        let r = sv*ov;
        r.is_zero()
    }
}

impl<T> Slice<T>
where
    T: Zero + One,
    T: Add<Output = T>,
    T: Div<Output = T>,
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn meet(&self, others: Vec<Self>) -> Option<GeoObj<T>> {
        let mut m : Vec<Vec<T>> = Vec::new();
        let mut d : Vec<T> = Vec::new();
        let mut sumv = RowVector::new(self.coords.clone());

        m.push(self.coords.clone());
        d.push(self.d.clone());

        for (_, o) in others.iter().enumerate() {
            m.push(o.coords.clone());
            d.push(o.d.clone());
            let rv = RowVector::new(o.coords.clone());
            sumv = sumv.clone() + rv;
        }

        let omatrix = Matrix::new(m.clone());

        if omatrix.rows > omatrix.cols {
            // Problem is overdetermined.
            return None;
        }

        if !omatrix.has_full_rank() {
            // If the matrix has not full rank, this means some of the slices
            // are parallel. Either those slices meet in all or no points.
            // In any case ...
            return None;
        }

        // Number of vectors must be zero or positive, otherwise the matrix
        // would have been rejected as overdetermined.
        let nr_of_vectors = omatrix.cols - omatrix.rows;

        for _n in 0..nr_of_vectors {
            // Create new slices to fill up the original matrix based on the
            // sum of all the previous rows and adding the column index to
            // each entry, to ensure linear independence.
            let mut rv = Vec::new();
            let mut a = T::zero();
            for (_, c) in sumv.elem.iter().enumerate() {
                rv.push(c.clone() + a.clone());
                a = a.clone() + T::one();
            }

            loop {
                let mut mt = m.clone();
                mt.push(rv.clone());
                let matrixt = Matrix::new(mt);

                if matrixt.has_full_rank() {
                    // matrix is invertible
                    m.push(rv.clone());
                    break;
                } else {
                    // Increase last entry in new row, to make it linear
                    // independent.
                    let mut c = rv.pop().unwrap();
                    c = c.clone() + T::one();
                    rv.push(c);
                }
            }

            sumv = sumv.clone() + RowVector::new(rv.clone());
        }

        let matrix = Matrix::new(m.clone());

        let mut d0 = d.clone();
        for _n in 0..nr_of_vectors {
            d0.push(T::zero());
        }

        let d0vector = ColumnVector::new(d0);

        let inv = matrix.inverse().unwrap();
        let p0 = inv.clone()*d0vector;
        let point0 = Point::new(p0.elem);

        let mut vectors = Vec::new();
        for n in 0..nr_of_vectors {
            // create more points to create vectors in the line / plane /
            //volume ... or, in general, meet.
            let mut d1 = d.clone();
            for r in 0..nr_of_vectors {
                if r == n {
                    d1.push(T::one());
                } else {
                    d1.push(T::zero());
                }
            }
            let d1vector = ColumnVector::new(d1);
            let p1 = inv.clone()*d1vector;
            let point1 = Point::new(p1.elem);
            let v1 = Vector::new(point0.clone(),point1.clone());
            vectors.push(v1.columnvector().elem);
        }

        return Some(GeoObj{ point: point0, vectors: vectors });
    }
}

impl<T> Slice<T>
where
    T: Zero + One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn projection_on_line_matrix(&self, direction: &Vector<T>) -> Matrix<T> {
        let pv = self.coords.clone();
        let dv = direction.columnvector().elem;

        let s = T::one() / (RowVector::<T>::new(pv.clone()) * ColumnVector::<T>::new(dv.clone()));

        let mut dmtv = Vec::new();
        dmtv.push(dv);
        let dmt = Matrix::<T>::new(dmtv);
        let dm = dmt.transpose();

        let mut pmv = Vec::new();
        pmv.push(pv);
        let pm = Matrix::<T>::new(pmv);
        dm*pm*s
    }

    pub fn projection_on_plane_matrix(&self, direction: &Vector<T>) -> Matrix<T> {
        let pm = self.projection_on_line_matrix(direction);
        let dvlen = direction.columnvector().elem.len();
        let pvlen = self.coords.len();
        let ident = Matrix::<T>::identity(dvlen,pvlen);

        ident - pm
    }

    pub fn reflection_across_plane_matrix(&self, direction: &Vector<T>) -> Matrix<T> {
        let pm = self.projection_on_line_matrix(direction);
        let dvlen = direction.columnvector().elem.len();
        let pvlen = self.coords.len();
        let ident = Matrix::<T>::identity(dvlen,pvlen);
        let two = T::one()+T::one();

        ident - pm*two
    }

    pub fn reflection_across_line_matrix(&self, direction: &Vector<T>) -> Matrix<T> {
        let mone = T::zero() - T::one();
        self.reflection_across_plane_matrix(direction) * mone
    }

    pub fn some_point(&self) -> Point<T> {
        let dim = self.coords.len();
        let mut coords = Vec::new();
        if self.d.is_zero() {
            for _i in 0..dim {
                coords.push(T::zero());
            }
            return Point::new(coords);
        }

        for _i in 0..dim {
            coords.push(T::one());
        }
        let v1 = ColumnVector::new(coords);
        let p1 = RowVector::new(self.coords.clone());
        let d1 = p1*v1.clone();
        let v2 = v1 * (self.d.clone() / d1);

        Point::new(v2.elem)
    }

    pub fn quadrance(&self, point: &Point<T>) -> T {
        let normal = Vector::from(self.coords.clone());
        let cva = ColumnVector::new(point.coords.clone());
        let pv = self.some_point();
        let cvv = ColumnVector::new(pv.coords);
        let cvg = cva - cvv;
        let vg = Vector::from(cvg.elem);
        let vu = normal.projection(&vg);
        vu.quadrance()
    }
}


/// Represents reflections and rotations
#[derive(Debug, Clone)]
pub struct Isometry<T> {
    matrix: Matrix<T>,
}

impl<T> Isometry<T>
where
    T: Zero + One,
    T: PartialEq,
    T: Neg<Output = T>,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new_reflection(coords: Vec<T>) -> Isometry<T> {
        let n = coords.clone().len();
        let i = Matrix::identity(n,n);
        let two = T::one() + T::one();
        let rv = RowVector::new(coords.clone());
        let cv = ColumnVector::new(coords);
        let m = i - cv.clone() * rv.clone() / (rv * cv) * two;
        Isometry{ matrix: m }
    }

    pub fn new_rotation(coords: Vec<T>) -> Isometry<T> {
        let refl= Self::new_reflection(coords);
        Isometry{ matrix: -refl.matrix }
    }

    pub fn is_reflection(&self) -> bool {
        self.matrix.determinant() == -T::one()
    }

    pub fn is_rotation(&self) -> bool {
        self.matrix.determinant() == T::one()
    }
}

impl<T> Mul<Isometry<T>> for Vector<T>
where
    T: Zero,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Vector<T>;

    fn mul(self, other: Isometry<T>) -> Vector<T> {
        let rv = self.rowvector();
        let nv = rv * other.matrix;
        Vector::from(nv.elem)
    }
}

impl<T> Mul for Isometry<T>
where
    T: Zero,
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: Clone,
{
    type Output = Isometry<T>;

    fn mul(self, other: Isometry<T>) -> Isometry<T> {
        let m = self.matrix * other.matrix;
        Isometry{ matrix: m }
    }
}

impl<T> Isometry<T>
where
    T: Zero + One,
    T: Neg<Output = T>,
    T: Div<Output = T>,
    T: Mul<Output = T>,
    T: Sub<Output = T>,
    T: Neg<Output = T>,
    T: Clone,
{
    pub fn anchor(&self) -> Option<Matrix<T>> {
        let ao = self.matrix.cayley();
        match ao {
            None => return None,
            Some(a) => return Some(-a),
        }
    }

    pub fn from_anchor(anchor: Matrix<T>) -> Option<Isometry<T>> {
        let ro = anchor.little_cayley();
        match ro {
            None => return None,
            Some(r) => return Some(Isometry{ matrix: r }),
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num::rational::{Ratio};

    #[test]
    #[should_panic]
    fn invalid_two_line() {
        Slice::new(vec![0,0],-1);
    }

    #[test]
    fn point_lies_on_line() {
        let a = Point::new(vec![Ratio::new( 6, 1), Ratio::new( -3, 2)]);
        let l = Slice::new(vec![Ratio::new(1, 1),  Ratio::new(2, 1)],
                           Ratio::new(3, 1));
        assert!(a.lies_on(&l));
    }

    #[test]
    fn points_lie_on_join() {
        let a = Point::new(vec![Ratio::new( 6, 1), Ratio::new( -3, 2)]);
        let b = Point::new(vec![Ratio::new( 4, 3), Ratio::new( 3, 5)]);
        let l = Slice::join(vec![a.clone(),b.clone()]).unwrap();
        assert!(a.lies_on(&l));
        assert!(b.lies_on(&l));
    }

    #[test]
    fn parallel_lines() {
        let l1 = Slice::new(vec![3, 4], 1);
        let l2 = Slice::new(vec![6, 8], -5);
        assert!(l1.is_parallel(&l2));
    }

    #[test]
    fn perpendicular_lines() {
        let l1 = Slice::new(vec![3, 4], 1);
        let l2 = Slice::new(vec![-4, 3], -2);
        assert!(l1.is_perpendicular(&l2));
    }

    #[test]
    fn points_lie_on_a_plane() {
        let p1 = Point::new(vec![3,0,0]);
        let p2 = Point::new(vec![0,1,0]);
        let p3 = Point::new(vec![0,0,2]);
        let pi = Slice::new(vec![2,6,3],6);
        assert!(p1.lies_on(&pi));
        assert!(p2.lies_on(&pi));
        assert!(p3.lies_on(&pi));
    }

    #[test]
    fn three_plains_meet_in_one_point() {
        let pi1 = Slice::new(vec![Ratio::from(1),Ratio::from(2),Ratio::from(-1)],Ratio::from(-3));
        let pi2 = Slice::new(vec![Ratio::from(3),Ratio::from(7),Ratio::from(2)],Ratio::from(1));
        let pi3 = Slice::new(vec![Ratio::from(4),Ratio::from(-2),Ratio::from(1)],Ratio::from(-2));

        let p1 = pi1.clone().meet(vec![pi2.clone(),pi3.clone()]).unwrap();
        let p2 = Point::new(vec![Ratio::from(-1),Ratio::from(0),Ratio::from(2)]);
        assert_eq!(p1.point,p2);
        let emptyvec : Vec<Vec<Ratio<i64>>> = Vec::new();
        assert_eq!(p1.vectors,emptyvec);
    }

    #[test]
    fn three_plains_contain_the_point_in_which_they_meet() {
        let pi1 = Slice::new(vec![Ratio::from(1),Ratio::from(2),Ratio::from(-1)],Ratio::from(-3));
        let pi2 = Slice::new(vec![Ratio::from(3),Ratio::from(7),Ratio::from(2)],Ratio::from(1));
        let pi3 = Slice::new(vec![Ratio::from(4),Ratio::from(-2),Ratio::from(1)],Ratio::from(-2));

        let point = pi1.clone().meet(vec![pi2.clone(),pi3.clone()]).unwrap();
        let plane1 = pi1.clone().meet(vec![]).unwrap();
        let plane2 = pi2.clone().meet(vec![]).unwrap();
        let plane3 = pi3.clone().meet(vec![]).unwrap();

        assert!(plane1.contains(&point));
        assert!(plane2.contains(&point));
        assert!(plane3.contains(&point));

        let line1 = pi2.clone().meet(vec![pi3.clone()]).unwrap();
        let line2 = pi3.clone().meet(vec![pi1.clone()]).unwrap();
        let line3 = pi1.clone().meet(vec![pi2.clone()]).unwrap();

        assert!(!plane1.contains(&line1));
        assert!(plane1.contains(&line2));
        assert!(plane1.contains(&line3));

        assert!(plane2.contains(&line1));
        assert!(!plane2.contains(&line2));
        assert!(plane2.contains(&line3));

        assert!(plane3.contains(&line1));
        assert!(plane3.contains(&line2));
        assert!(!plane3.contains(&line3));
    }

    #[test]
    fn line_passes_through_point() {
        let p1 = Point::new(vec![3,1,1]);
        let p2 = Point::new(vec![2,2,-1]);
        let p3 = GeoObj::new(Point::new(vec![1,3,-3]),vec![]);
        let p4 = GeoObj::new(Point::new(vec![0,4,3]),vec![]);
        let v = Vector::new(p1.clone(),p2).columnvector().elem;
        let l = GeoObj::new(p1,vec![v]);
        assert!(l.contains(&p3));
        assert!(!l.contains(&p4));
    }

    #[test]
    fn meet_of_line_and_plane() {
        let p = Point::new(vec![Ratio::from(1),Ratio::from(2),Ratio::from(-3)]);
        let v = vec![Ratio::from(-1),Ratio::from(4),Ratio::from(1)];
        let l = GeoObj::new(p,vec![v]);
        let pi = Slice::new(vec![Ratio::from(1),Ratio::from(3),Ratio::from(-1)],Ratio::from(5));
        let plane = pi.meet(vec![]).unwrap();
        let pm = l.meet(&plane).unwrap();
        assert!(l.contains(&pm));
        assert!(plane.contains(&pm));
        let pr = Point::new(vec![Ratio::new(3,2),Ratio::from(0),Ratio::new(-7,2)]);
        assert_eq!(pm.point,pr);
    }

    #[test]
    fn two_planes_meet_in_one_line() {
        let pi1 = Slice::new(vec![Ratio::from(1),Ratio::from(2),Ratio::from(-1)],Ratio::from(-3));
        let pi2 = Slice::new(vec![Ratio::from(3),Ratio::from(7),Ratio::from(2)],Ratio::from(1));
        let plane1 = pi1.meet(vec![]).unwrap();
        let plane2 = pi2.meet(vec![]).unwrap();

        let go1 = pi1.meet(vec![pi2.clone()]).unwrap();
        assert!(plane1.contains(&go1));
        assert!(plane2.contains(&go1));

        let go2 = plane1.meet(&plane2).unwrap();
        assert!(plane1.contains(&go2));
        assert!(plane2.contains(&go2));
        assert_eq!(go1,go2);

        let p = Point::new(vec![Ratio::from(-23),Ratio::from(10),Ratio::from(0)]);
        let v = vec![Ratio::from(11),Ratio::from(-5),Ratio::from(1)];
        let l2 = GeoObj::new(p,vec![v]);
        assert_eq!(go2,l2);
    }

    #[test]
    fn two_planes_meet_in_one_line_2() {
        let pi1 = Slice::new(vec![Ratio::from(1),Ratio::from(3),Ratio::from(-2)],Ratio::from(2));
        let pi2 = Slice::new(vec![Ratio::from(2),Ratio::from(6),Ratio::from(-5)],Ratio::from(3));
        let go = pi1.meet(vec![pi2.clone()]).unwrap();
        let plane1 = pi1.meet(vec![]).unwrap();
        let plane2 = pi2.meet(vec![]).unwrap();
        assert!(plane1.contains(&go));
        assert!(plane2.contains(&go));

        let p = Point::new(vec![Ratio::from(4),Ratio::from(0),Ratio::from(1)]);
        let v = vec![Ratio::from(-3),Ratio::from(1),Ratio::from(0)];
        let l2 = GeoObj::new(p,vec![v]);
        assert_eq!(go,l2);
    }

    #[test]
    fn three_points_define_a_plane() {
        let p1 = Point::new(vec![1,0,-1]);
        let p2 = Point::new(vec![2,1,3]);
        let p3 = Point::new(vec![-4,2,5]);
        let pi = GeoObj::from(vec![p1.clone(),p2.clone(),p3.clone()]);
        assert!(pi.contains(&GeoObj{point: p1, vectors: vec![]}));
        assert!(pi.contains(&GeoObj{point: p2, vectors: vec![]}));
        assert!(pi.contains(&GeoObj{point: p3, vectors: vec![]}));
    }

    #[test]
    fn three_points_define_a_plane_2() {
        let p1 = Point::new(vec![1,0,-1]);
        let p2 = Point::new(vec![2,1,3]);
        let p3 = Point::new(vec![-4,2,5]);
        let pi = Slice::join(vec![p1.clone(),p2.clone(),p3.clone()]).unwrap();
        assert!(p1.lies_on(&pi));
        assert!(p2.lies_on(&pi));
        assert!(p3.lies_on(&pi));
    }

    #[test]
    fn one_point_and_two_vectors_define_a_plane() {
        let a = Point::new(vec![1,0,-1]);
        let u = vec![1,1,4];
        let uv = Vector::from(u.clone());
        let v = vec![-5,2,6];
        let vv = Vector::from(v.clone());
        let pi = GeoObj::new(a.clone(),vec![u.clone(),v.clone()]);
        assert!(pi.contains(&GeoObj{point: a.clone(), vectors: vec![]}));
        assert!(pi.contains(&GeoObj{point: (a.clone()+uv), vectors: vec![]}));
        assert!(pi.contains(&GeoObj{point: (a+vv), vectors: vec![]}));
    }

    #[test]
    fn point_normal_form_of_a_line() {
        // In 2D, slices are lines
        let nv = Vector::from(vec![Ratio::from(2),Ratio::from(-1)]);
        let p = Point::new(vec![Ratio::from(-1),Ratio::from(3)]);
        let l = Slice::from(p.clone(),nv);
        assert!(p.lies_on(&l));
        // Adding a direction vector
        let dv = Vector::from(vec![Ratio::from(1),Ratio::from(2)]);
        assert!((p+dv).lies_on(&l));
    }

    #[test]
    fn equal_vectors() {
        let a1 = Point::new(vec![Ratio::new(1,1), Ratio::new(1,1)]);
        let a2 = Point::new(vec![Ratio::new(3,1), Ratio::new(4,1)]);
        let a3 = Point::new(vec![Ratio::new(4,1), Ratio::new(2,1)]);
        let a4 = Point::new(vec![Ratio::new(6,1), Ratio::new(5,1)]);
        let v1 = Vector::new(a1, a2);
        let v2 = Vector::new(a3, a4);
        assert_eq!(v1,v2);
    }

    #[test]
    fn add_vectors() {
        let v1 = Vector::from(vec![Ratio::new(2,1), Ratio::new(3,1)]);
        let v2 = Vector::from(vec![Ratio::new(2,1), Ratio::new(-4,1)]);
        let v3 = Vector::from(vec![Ratio::new(4,1), Ratio::new(-1,1)]);
        assert_eq!(v1+v2,v3);
    }

    #[test]
    fn add_point_and_vector() {
        let p1 = Point::new(vec![1,2,3]);
        let v1 = Vector::from(vec![4,5,6]);
        let p2 = Point::new(vec![5,7,9]);
        assert_eq!(p1+v1,p2);
    }

    #[test]
    fn scalar_vector_mul() {
        let s = Ratio::new(3,1);
        let v1 = Vector::from(vec![Ratio::new(2,1), Ratio::new(-4,1)]);
        let v2 = Vector::from(vec![Ratio::new(6,1), Ratio::new(-12,1)]);
        assert_eq!(v1*s,v2);
    }

    #[test]
    fn projections_onto_a_plane() {
        let pi1 = Slice::new(vec![Ratio::from(1),Ratio::from(-5),Ratio::from(2)],Ratio::from(0));
        let v1 = Vector::from(vec![Ratio::from(1),Ratio::from(1),Ratio::from(3)]);
        let v2 = Vector::from(vec![Ratio::from(1),Ratio::from(-5),Ratio::from(2)]);
        let pm1 = pi1.projection_on_plane_matrix(&v1);
        let pm2 = pi1.projection_on_plane_matrix(&v2);
        assert_eq!(pm1.clone()*pm1.clone(),pm1.clone());
        assert_eq!(pm2.clone()*pm2.clone(),pm2.clone());

        let u1 = ColumnVector::new(vec![Ratio::from(12),Ratio::from(-341),Ratio::from(35)]);
        let u2 = ColumnVector::new(vec![Ratio::from(-123),Ratio::from(-55),Ratio::from(32)]);

        let pu11 = pm1.clone()*u1.clone();
        let pu12 = pm1*u2.clone();
        let pu21 = pm2.clone()*u1;
        let pu22 = pm2*u2;

        let vpu11 = Vector::from(pu11.elem);
        let vpu12 = Vector::from(pu12.elem);
        let vpu21 = Vector::from(pu21.elem);
        let vpu22 = Vector::from(pu22.elem);

        let o = Point::new(vec![Ratio::from(0),Ratio::from(0),Ratio::from(0)]);
        assert!((o.clone() + vpu11).lies_on(&pi1));
        assert!((o.clone() + vpu12).lies_on(&pi1));
        assert!((o.clone() + vpu21).lies_on(&pi1));
        assert!((o + vpu22).lies_on(&pi1));
    }

    #[test]
    fn projection_onto_a_line_1() {
        // In 2D slices are lines
        let coords = vec![Ratio::from(4),Ratio::from(1)];
        let l = Slice::new(coords.clone(),Ratio::from(0));
        // Projection along the normal to the plane (line)
        let nv = Vector::from(coords);
        let pm = l.projection_on_plane_matrix(&nv);
        assert_eq!(pm.clone()*pm.clone(),pm.clone());

        let v1 = ColumnVector::new(vec![Ratio::from(3),Ratio::from(2)]);

        // Projection vector
        let pu = pm*v1;
        let vpu = Vector::from(pu.elem);

        // Origin
        let o = Point::new(vec![Ratio::from(0),Ratio::from(0)]);
        let pp = o + vpu;
        assert!(pp.lies_on(&l));
    }

    #[test]
    fn projection_onto_a_line_2() {
        // In 2D slices are lines
        let coords = vec![Ratio::from(3),Ratio::from(2)];
        let l = Slice::new(coords.clone(),Ratio::from(0));
        // Projection along the normal to the plane (line)
        let nv = Vector::from(coords);
        let pm = l.projection_on_plane_matrix(&nv);
        assert_eq!(pm.clone()*pm.clone(),pm.clone());

        let v1 = ColumnVector::new(vec![Ratio::from(4),Ratio::from(1)]);

        // Projection vector
        let pu = pm*v1;
        let vpu = Vector::from(pu.elem);

        // Origin
        let o = Point::new(vec![Ratio::from(0),Ratio::from(0)]);
        let pp = o + vpu;
        assert!(pp.lies_on(&l));
    }

    #[test]
    fn quadrance_from_line_to_point_1() {
        // In 2D, slices are lines
        let nv = Vector::from(vec![Ratio::from(3),Ratio::from(4)]);
        let p = Point::new(vec![Ratio::from(0),Ratio::from(1)]);
        let l = Slice::from(p.clone(),nv);

        // Point A
        let a = Point::new(vec![Ratio::from(1),Ratio::from(2)]);

        let q = l.quadrance(&a);
        assert_eq!(q,Ratio::new(49,25));
    }

    #[test]
    fn quadrance_from_line_to_point_2() {
        // In 2D, slices are lines
        let l = Slice::new(vec![Ratio::from(1),Ratio::from(-2)],Ratio::from(5));

        // Point A
        let a = Point::new(vec![Ratio::from(5),Ratio::from(-2)]);

        let q = l.quadrance(&a);
        assert_eq!(q,Ratio::new(16,5));
    }

    #[test]
    fn projection_of_vector_onto_vector() {
        let v1 = Vector::from(vec![Ratio::from(3),Ratio::from(2)]);
        let v2 = Vector::from(vec![Ratio::from(4),Ratio::from(1)]);
        let u1 = v2.projection(&v1);
        let u2 = v1.projection(&v2);
        let wmu1 = v1.clone() - u1.clone();
        let wmu2 = v2.clone() - u2.clone();
        assert!(wmu1.is_perpendicular(&u1));
        assert!(wmu2.is_perpendicular(&u2));
        let lambda1 = v2.columnvector().get(0) / u1.columnvector().get(0);
        let lambda2 = v1.columnvector().get(0) / u2.columnvector().get(0);
        assert_eq!(u1 * lambda1, v2);
        assert_eq!(u2 * lambda2, v1);
    }

    #[test]
    fn quadrance() {
        let p1 = Point::new(vec![Ratio::new(2,1), Ratio::new(1,1)]);
        let p2 = Point::new(vec![Ratio::new(6,1), Ratio::new(2,1)]);
        let q  = Ratio::new(17,1);
        assert_eq!(p1.quadrance(&p2),q);
    }

    #[test]
    fn quadrance_between_one_points() {
        let a1 = Point::new(vec![2]);
        let a2 = Point::new(vec![5]);
        let a3 = Point::new(vec![7]);
        assert_eq!(a2.quadrance(&a3),4);
        assert_eq!(a1.quadrance(&a3),25);
        assert_eq!(a1.quadrance(&a2),9);
    }

    #[test]
    fn quadrances_in_right_triangle() {
        let a = Point::new(vec![1,2,4]);
        let b = Point::new(vec![4,3,2]);
        let c = Point::new(vec![5,4,4]);
        let qab = a.quadrance(&b);
        let qbc = b.quadrance(&c);
        let qca = c.quadrance(&a);
        assert_eq!(qab+qbc,qca);
    }

    #[test]
    fn perpendicular_vectors() {
        let v = Vector::from(vec![1,-2,3]);
        let u = Vector::from(vec![4,-1,-2]);
        assert!(v.is_perpendicular(&u));
    }

    #[test]
    fn reflection_in_normal() {
        let v = Vector::from(vec![Ratio::from(1),Ratio::from(-2),Ratio::from(3)]);
        let v1 = Vector::from(vec![Ratio::from(-4),Ratio::from(5),Ratio::from(6)]);
        let v2 = Vector::from(vec![Ratio::from(7),Ratio::from(-8),Ratio::from(-9)]);

        let mv1 = v.reflection_in_normal(&v1);
        let mv2 = v.reflection_in_normal(&v2);

        assert_eq!(mv1.clone() * mv2.clone(),v1.clone() * v2.clone());

        let mmv1 = v.reflection_in_normal(&mv1);
        let mmv2 = v.reflection_in_normal(&mv2);

        assert_eq!(mmv1,v1);
        assert_eq!(mmv2,v2);
    }

    #[test]
    fn reflection_in_vector() {
        let v = Vector::from(vec![Ratio::from(1),Ratio::from(-2),Ratio::from(3)]);
        let v1 = Vector::from(vec![Ratio::from(-4),Ratio::from(5),Ratio::from(6)]);
        let v2 = Vector::from(vec![Ratio::from(7),Ratio::from(-8),Ratio::from(-9)]);

        let mv1 = v.reflection_in_vector(&v1);
        let mv2 = v.reflection_in_vector(&v2);

        assert_eq!(mv1.clone() * mv2.clone(),v1.clone() * v2.clone());

        let mmv1 = v.reflection_in_vector(&mv1);
        let mmv2 = v.reflection_in_vector(&mv2);

        assert_eq!(mmv1,v1);
        assert_eq!(mmv2,v2);
    }

    #[test]
    fn reflection_matrix() {
        let elem1 = vec![Ratio::from(1),Ratio::from(2),Ratio::from(3)];
        let r = Isometry::new_reflection(elem1.clone());
        let v11 = Vector::from(elem1);
        let v12 = v11.clone() * r.clone();
        assert_eq!(v11,-v12);

        let elem2 = vec![Ratio::from(-2),Ratio::from(1),Ratio::from(0)];
        let v21 = Vector::from(elem2);
        let v22 = v21.clone() * r.clone();
        assert_eq!(v21,v22);

        let elem3 = vec![Ratio::from(-3),Ratio::from(0),Ratio::from(1)];
        let v31 = Vector::from(elem3);
        let v32 = v31.clone() * r;
        assert_eq!(v31,v32);
    }

    #[test]
    fn multiplication_of_reflections() {
        let elemv = vec![Ratio::from(1),Ratio::from(2),Ratio::from(3)];
        let elemw = vec![Ratio::from(2),Ratio::from(1),Ratio::from(-1)];

        let rv = Isometry::new_reflection(elemv);
        let rw = Isometry::new_reflection(elemw);

        let r= rv*rw;

        let elemf = vec![Ratio::from(5),Ratio::from(-7),Ratio::from(3)];
        let vf = Vector::from(elemf);

        assert_eq!(vf.clone()*r,vf);
    }

    #[test]
    fn rotations_in_3d() {
        let elemv1 = vec![Ratio::from(1),Ratio::from(2),Ratio::from(3)];
        let elemv2 = vec![Ratio::from(2),Ratio::from(-1),Ratio::from(5)];

        let sigmav1 = Isometry::new_reflection(elemv1);
        let sigmav2 = Isometry::new_reflection(elemv2);

        let rho = sigmav1*sigmav2;
        assert!(rho.is_rotation());
    }

    #[test]
    fn spread_between_2d_vectors() {
        let o = Point::new(vec![Ratio::new(0,1), Ratio::new(0,1)]);
        let a = Point::new(vec![Ratio::new(4,1), Ratio::new(2,1)]);
        let b = Point::new(vec![Ratio::new(-1,1),Ratio::new(3,1)]);
        let oa = Vector::new(o.clone(), a.clone());
        let ob = Vector::new(o, b.clone());
        let ab = Vector::new(a, b);
        let s = oa.spread(&ob);
        let t = ob.spread(&ab);
        let r = oa.spread(&ab);
        assert_eq!(s,Ratio::new(49,50));
        assert_eq!(t,Ratio::new(49,65));
        assert_eq!(r,Ratio::new(49,130));
    }

    #[test]
    fn spread_between_3d_vectors() {
        let elemv1 = vec![Ratio::from(1),Ratio::from(2),Ratio::from(3)];
        let elemv2 = vec![Ratio::from(2),Ratio::from(-1),Ratio::from(5)];
        let v1 = Vector::from(elemv1);
        let v2 = Vector::from(elemv2);
        assert_eq!(v1.spread(&v2),Ratio::new(13,28));
    }

    #[test]
    fn anchor_of_rotations() {
        let zero = Ratio::from(0);
        let one = Ratio::from(1);
        let two = Ratio::from(2);
        let three = Ratio::from(3);
        let four = Ratio::from(4);
        let five = Ratio::from(5);
        let r12 = Ratio::from(12);
        let r19 = Ratio::from(19);
        let r28 = Ratio::from(28);

        let ma = Matrix::new(vec![vec![zero,five,-one],
                                  vec![-five,zero,four],
                                  vec![one,-four,zero]]);

        let mb = Matrix::new(vec![vec![zero,two,-four],
                                  vec![-two,zero,-three],
                                  vec![four,three,zero]]);

        let mc = Matrix::new(vec![vec![zero,r12,r28],
                                  vec![-r12,zero,-r19],
                                  vec![-r28,r19,zero]]);

        let rho_a = Isometry::from_anchor(ma).unwrap();
        let rho_b = Isometry::from_anchor(mb).unwrap();
        let rho_c = rho_a*rho_b;
        let mc2 = rho_c.anchor().unwrap();

        assert_eq!(mc2,mc);
    }

    #[test]
    fn quadrance_in_red_geometry() {
        let a1 = Point::new(vec![Ratio::new(1,1), Ratio::new(0,1)]);
        let a2 = Point::new(vec![Ratio::new(13,12), Ratio::new(5,12)]);
        let a3 = Point::new(vec![Ratio::new(-5,3), Ratio::new(4,3)]);
        let a4 = Point::new(vec![Ratio::new(-17,15), Ratio::new(-8,15)]);

        let v12 = Vector::new_red(a1.clone(),a2.clone());
        let v23 = Vector::new_red(a2.clone(),a3.clone());
        let v34 = Vector::new_red(a3.clone(),a4.clone());
        let v14 = Vector::new_red(a1.clone(),a4.clone());
        let v13 = Vector::new_red(a1,a3);
        let v24 = Vector::new_red(a2,a4);

        assert_eq!(v12.quadrance(), Ratio::new(-1,6));
        assert_eq!(v23.quadrance(), Ratio::new(121,18));
        assert_eq!(v34.quadrance(), Ratio::new(-16,5));
        assert_eq!(v14.quadrance(), Ratio::new(64,15));
        assert_eq!(v13.quadrance(), Ratio::new(16,3));
        assert_eq!(v24.quadrance(), Ratio::new(361,90));
    }
}

