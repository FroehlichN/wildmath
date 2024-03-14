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
use std::ops::{Mul, Div, Add, Sub};


/// Represents a 3D point
#[derive(Debug, Clone)]
pub struct Point<T> {
    pub x : T,
    pub y : T,
    pub z : T,
}


impl<T> Point<T>
{
    pub fn new(x: T, y: T, z: T) -> Point<T> {
        Point { x: x, y: y, z: z }
    }
}


impl<T> Point<T>
where
    T: Zero,
    T: Add,
    T: Mul<Output = T>,
    T: Clone,
{
    pub fn lies_on(self, pi: &Plane<T>) -> bool {
        let r = pi.a.clone() * self.x.clone()
              + pi.b.clone() * self.y.clone()
              + pi.c.clone() * self.z.clone()
              + pi.d.clone();
        r.is_zero()
    }
}


/// Represents a plane in 3D
/// a*x + b*y + c*z + d = 0
#[derive(Debug, Clone)]
pub struct Plane<T> {
    a: T,
    b: T,
    c: T,
    d: T,
}


impl<T> Plane<T>
{
    pub fn new(a: T, b: T, c: T, d: T) -> Plane<T> {
        Plane { a: a, b: b, c: c, d: d }
    }
}


/// Represents a 3D tetrahedon
/// T is type of the coordinates of the points
pub struct Tetrahedron<T> {
    points: [Point<T>; 3],
}


impl<T> Tetrahedron<T>
where
    T: One,
    T: Add<Output = T>,
    T: Sub<Output = T>,
    T: Mul<Output = T>,
    T: Div<Output = T>,
    T: Clone,
{
    pub fn new(a1: Point<T>, a2: Point<T>, a3: Point<T>)
        -> Tetrahedron<T> {
/*        if a1.is_collinear(&a2, &a3) {
            panic!("Tetrahedron cannot have collinear points");
        }*/
        Tetrahedron { points: [a1, a2, a3] }
    }

    pub fn volume(&self) -> T {
        let six = T::one() + T::one() + T::one()
                + T::one() + T::one() + T::one();
        let x1 = self.points[0].x.clone();
        let x2 = self.points[1].x.clone();
        let x3 = self.points[2].x.clone();
        let y1 = self.points[0].y.clone();
        let y2 = self.points[1].y.clone();
        let y3 = self.points[2].y.clone();
        let z1 = self.points[0].z.clone();
        let z2 = self.points[1].z.clone();
        let z3 = self.points[2].z.clone();

        let n = x1.clone()*y2.clone()*z3.clone()
              - x1*y3.clone()*z2.clone()
              + x2.clone()*y3*z1.clone()
              - x3.clone()*y2*z1
              + x3*y1.clone()*z2
              - x2*y1*z3;
        return n / six;
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn points_lie_on_a_plane() {
        let p1 = Point::new(3,0,0);
        let p2 = Point::new(0,1,0);
        let p3 = Point::new(0,0,2);
        let pi = Plane::new(2,6,3,-6);
        assert!(p1.lies_on(&pi));
        assert!(p2.lies_on(&pi));
        assert!(p3.lies_on(&pi));
    }
}

