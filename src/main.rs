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

use affinegeo::*;
use polynumber::*;

fn main() {
    meet_of_fermat_and_bernoulli();
}

fn meet_of_fermat_and_bernoulli() {
    println!("Intersect of Fermat curve and Lemniscate of Bernoulli:");
    // Fermat curve p = -1 + x^3 + y^3
    let fp = PolyNumber::new(vec![
                PolyNumber::new(vec![-1f64, 0f64, 0f64, 1f64] ),
                PolyNumber::new(vec![ 0f64] ),
                PolyNumber::new(vec![ 0f64] ),
                PolyNumber::new(vec![ 1f64] ) ] );

    // Lemniscate of Bernoulli (x^2 + y^2)^2 = 2*(x^2 - y^2)
    let lp = PolyNumber::new(vec![
                PolyNumber::new(vec![ 0f64, 0f64, -2f64, 0f64, 1f64] ),
                PolyNumber::new(vec![ 0f64] ),
                PolyNumber::new(vec![ 2f64, 0f64, 2f64] ),
                PolyNumber::new(vec![ 0f64] ),
                PolyNumber::new(vec![ 1f64] ) ] );

    let mut x = 2f64;
    let mut y = 1f64;

    for n in 1..7 {
        // Tangent plane to Fermat curve at [2,1]
        let ft = fp.clone().tangent2(1,x,y);
        // Tangent plane to Lemniscate at [2,1]
        let lt = lp.clone().tangent2(1,x,y);

        // convert to TwoLines
        let ftl = TwoLine::new(ft.get(0).get(1),ft.get(1).get(0),ft.get(0).get(0));
        let ltl = TwoLine::new(lt.get(0).get(1),lt.get(1).get(0),lt.get(0).get(0));

        // Calculate meet of two TwoLines
        let m = ftl.meet(&ltl);

        x = m.x;
        y = m.y;

        println!("Iteration {:?} : {:?}",n,m);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use polyratio::*;
    use polynumber::*;
    use extrational::*;
    use affinegeo::*;
    use finite::*;
    use projectivegeo::*;

    create_finite_field!(3);
    create_finite_field!(7);
    create_finite_field!(11);
    create_finite_field!(13);

    #[test]
    fn rational_poly_number_with_infinity() {
        let numer = PolyNumber::new(vec![RatInf::new(1,1),RatInf::new(2,1)]);
        let denom = PolyNumber::new(vec![RatInf::new(1,1),RatInf::new(-1,1)]);
        let r = PolyRatio::new(numer,denom);
        assert_eq!(r.eval(RatInf::new(1,1)),RatInf::new(1,0));
        assert!(r.eval(RatInf::new(1,0)).is_nil());
    }
    #[test]
    fn bretschneider_von_staudt_over_f7() {
        let p1 = TwoPoint::new(Finite7::new(1), Finite7::new(1));
        let p2 = TwoPoint::new(Finite7::new(2), Finite7::new(2));
        let p3 = TwoPoint::new(Finite7::new(3), Finite7::new(4));
        let p4 = TwoPoint::new(Finite7::new(5), Finite7::new(6));

        let q = Quadrilateral::new(p1.clone(),p2.clone(),p3.clone(),p4.clone());
        let p = Polygon::new(vec![p1, p2, p3, p4]);
        let a = p.area();
        assert_eq!(Finite7::new(16)*a*a,q.quadrea());
    }
    #[test]
    fn bretschneider_von_staudt_over_f11() {
        let p1 = TwoPoint::new(Finite11::new(1), Finite11::new(1));
        let p2 = TwoPoint::new(Finite11::new(2), Finite11::new(2));
        let p3 = TwoPoint::new(Finite11::new(3), Finite11::new(4));
        let p4 = TwoPoint::new(Finite11::new(5), Finite11::new(6));

        let q = Quadrilateral::new(p1.clone(),p2.clone(),p3.clone(),p4.clone());
        let p = Polygon::new(vec![p1, p2, p3, p4]);
        let a = p.area();
        assert_eq!(Finite11::new(16)*a*a,q.quadrea());
    }
    #[test]
    fn bretschneider_von_staudt_over_f13() {
        let p1 = TwoPoint::new(Finite13::new(1), Finite13::new(1));
        let p2 = TwoPoint::new(Finite13::new(2), Finite13::new(2));
        let p3 = TwoPoint::new(Finite13::new(3), Finite13::new(4));
        let p4 = TwoPoint::new(Finite13::new(5), Finite13::new(6));

        let q = Quadrilateral::new(p1.clone(),p2.clone(),p3.clone(),p4.clone());
        let p = Polygon::new(vec![p1, p2, p3, p4]);
        let a = p.area();
        assert_eq!(Finite13::new(16)*a*a,q.quadrea());
    }
    #[test]
    fn unit_circle_in_f7() {
        let center = TwoPoint::new(Finite7::new(0), Finite7::new(0));
        let quadrance = Finite7::new(1);
        let a1 = TwoPoint::new(Finite7::new(0), Finite7::new(1));
        let a2 = TwoPoint::new(Finite7::new(2), Finite7::new(2));
        let a3 = TwoPoint::new(Finite7::new(5), Finite7::new(5));
        let a4 = TwoPoint::new(Finite7::new(6), Finite7::new(0));
        let circle = TwoCircle::new(center, quadrance);
        assert!(circle.lies_on(&a1));
        assert!(circle.lies_on(&a2));
        assert!(circle.lies_on(&a3));
        assert!(circle.lies_on(&a4));
        assert_eq!(a1.quadrance(&a2), Finite7::new(5));
        assert_eq!(a2.quadrance(&a3), Finite7::new(4));
        assert_eq!(a3.quadrance(&a4), Finite7::new(5));
        assert_eq!(a1.quadrance(&a4), Finite7::new(2));
        assert_eq!(a1.quadrance(&a3), Finite7::new(6));
        assert_eq!(a2.quadrance(&a4), Finite7::new(6));
    }
    #[test]
    fn quadrance_of_projective_one_points_in_f3() {
        let am1 = ProjOnePoint::new(Finite3::new(1),Finite3::new(-1));
        let a0 = ProjOnePoint::new(Finite3::new(1),Finite3::new(0));
        let a1 = ProjOnePoint::new(Finite3::new(1),Finite3::new(1));
        let ainf = ProjOnePoint::new(Finite3::new(0),Finite3::new(1));

        let half = Finite3::new(1)/Finite3::new(2);
        let one = Finite3::new(1);
        let zero = Finite3::new(0);

        assert_eq!(am1.quadrance(&am1),zero);
        assert_eq!(am1.quadrance(&a0),half);
        assert_eq!(am1.quadrance(&a1),one);
        assert_eq!(am1.quadrance(&ainf),half);

        assert_eq!(a0.quadrance(&am1),half);
        assert_eq!(a0.quadrance(&a0),zero);
        assert_eq!(a0.quadrance(&a1),half);
        assert_eq!(a0.quadrance(&ainf),one);

        assert_eq!(a1.quadrance(&am1),one);
        assert_eq!(a1.quadrance(&a0),half);
        assert_eq!(a1.quadrance(&a1),zero);
        assert_eq!(a1.quadrance(&ainf),half);

        assert_eq!(ainf.quadrance(&am1),half);
        assert_eq!(ainf.quadrance(&a0),one);
        assert_eq!(ainf.quadrance(&a1),half);
        assert_eq!(ainf.quadrance(&ainf),zero);
    }
}

