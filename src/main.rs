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

use wildmath::affinegeo2d::*;
use wildmath::polynumber::*;
use num::rational::{Ratio};
use wildmath::matrix::{Matrix,ColumnVector};
use wildmath::belnap::Belnap;

fn main() {
    meet_of_fermat_and_bernoulli();
    eigenvalues_of_2d_matracies();
    line_trough_two_point();
    parabola_through_three_points();
    belnap_logic_of_penguins_birds_flying();
    belnap_logic_identity();
    belnap_logic_implication();
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

fn eigenvalues_of_2d_matracies() {
    let m1 = Matrix::new(vec![vec![Ratio::from(3),Ratio::from(2)],
                              vec![Ratio::from(0),Ratio::from(5)]]);
    let ev1 = m1.eigenvalues();
    println!("Eigenvalues of {:?} = {:?}",m1,ev1);

    let m2 = Matrix::new(vec![vec![Ratio::from(4),Ratio::from(6)],
                              vec![Ratio::from(6),Ratio::from(4)]]);
    let ev2 = m2.eigenvalues();
    println!("Eigenvalues of {:?} = {:?}",m2,ev2);

    let m3 = Matrix::new(vec![vec![Ratio::from(13),Ratio::from(-9)],
                              vec![Ratio::from(24),Ratio::from(-17)]]);
    let ev3 = m3.eigenvalues();
    println!("Eigenvalues of {:?} = {:?}",m3,ev3);
}

fn line_trough_two_point() {
    let m1 = Matrix::new(vec![vec![Ratio::from(1), Ratio::from(-1)],
                              vec![Ratio::from(1), Ratio::from(5)]]);
    let v1 = ColumnVector::new(vec![Ratio::from(2),Ratio::from(4)]);
    let v2 = m1.inverse().unwrap() * v1;
    println!("Line going through [-1,2] and [5,4]: y = {:?} + {:?} * x",v2.get(0),v2.get(1));
}

fn parabola_through_three_points() {
    let m1 = Matrix::new(vec![vec![Ratio::from(1), Ratio::from(-3), Ratio::from(9)],
                              vec![Ratio::from(1), Ratio::from(0), Ratio::from(0)],
                              vec![Ratio::from(1), Ratio::from(1), Ratio::from(1)]]);
    let v1 = ColumnVector::new(vec![Ratio::from(4), Ratio::from(-2), Ratio::from(5)]);
    let v2 = m1.inverse().unwrap() * v1;
    println!("Parabola going through [-3,4], [0,-2], and [1,5]: y = {:?} + {:?} * x + {:?} * x^2",v2.get(0),v2.get(1),v2.get(2));
}

fn belnap_logic_of_penguins_birds_flying() {
    println!("");
    println!("phi = is it a penguin AND penguins are birds AND birds can fly AND penguins cannot fly.");
    println!("In Belnap logic, there are non-false cases for this normally inconsistent truth statement phi:");
    let v = [Belnap::True, Belnap::False, Belnap::Both, Belnap::Neither];

    // is it a penguin?
    for p in v {
        // it is a bird?
        for b in v {
            // can it fly?
            for f in v {
                // is always false in classical two-valued logic
                let phi = p & p.mimplies(b) & b.mimplies(f) & p.mimplies(!f);
                // list the cases, when it is not
                if !(phi == Belnap::False) {
                    println!("phi = {:?}, P = {:?}, B = {:?}, F = {:?}",phi,p,b,f);
                }
            }
        }
    }
}

fn belnap_logic_identity() {
    println!("");
    println!("Belnap logic of a => b AND b => a");
    let v = [Belnap::True, Belnap::False, Belnap::Both, Belnap::Neither];

    let mut header = format!("        |");
    for a in v {
        header = format!("{header}{:<8}",format!("{:?}",a));
    }

    println!("{header}");
    println!("-----------------------------------------");

    for a in v {
        let mut line = format!("{:<8}|",format!("{:?}",a));
        for b in v {
            let c = a.mimplies(b) & b.mimplies(a);
            line = format!("{line}{:<8}",format!("{:?}",c));
        }
        println!("{line}");
    }
}

fn belnap_logic_implication() {
    println!("");
    println!("Belnap logic of a => b AND -a => -b");
    let v = [Belnap::True, Belnap::False, Belnap::Both, Belnap::Neither];

    let mut header = format!("        |");
    for a in v {
        header = format!("{header}{:<8}",format!("{:?}",a));
    }

    println!("{header}");
    println!("-----------------------------------------");

    for a in v {
        let mut line = format!("{:<8}|",format!("{:?}",a));
        for b in v {
            let c = a.mimplies(b) & (!a).mimplies(!b);
            line = format!("{line}{:<8}",format!("{:?}",c));
        }
        println!("{line}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wildmath::polyratio::*;
    use wildmath::extrational::*;
    use wildmath::finite::*;
    use wildmath::create_finite_field;
    use wildmath::projectivegeo::*;
    use num::rational::Ratio;
    use wildmath::create_polynumber_var;
    use wildmath::create_polynumber_one;

    create_finite_field!(3);
    create_finite_field!(5);
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
    #[test]
    fn rotations_of_projective_one_points_in_f3() {
        let am1 = ProjOnePoint::new(Finite3::new(1),Finite3::new(-1));
        let a0 = ProjOnePoint::new(Finite3::new(1),Finite3::new(0));
        let a1 = ProjOnePoint::new(Finite3::new(1),Finite3::new(1));
        let ainf = ProjOnePoint::new(Finite3::new(0),Finite3::new(1));

        let r1 = wildmath::projectivegeo::Rotation::new(Finite3::new(1),Finite3::new(1));
        let rm1 = wildmath::projectivegeo::Rotation::new(Finite3::new(1),Finite3::new(-1));
        let rinf = wildmath::projectivegeo::Rotation::new(Finite3::new(0),Finite3::new(1));

        assert_eq!(a0.clone()*r1.clone(),a1);
        assert_eq!(a1.clone()*r1.clone(),ainf);
        assert_eq!(ainf.clone()*r1.clone(),am1);
        assert_eq!(am1.clone()*r1.clone(),a0);

        assert_eq!(a0.clone()*rm1.clone(),am1);
        assert_eq!(a1.clone()*rm1.clone(),a0);
        assert_eq!(ainf.clone()*rm1.clone(),a1);
        assert_eq!(am1.clone()*rm1.clone(),ainf);

        assert_eq!(a0.clone()*rinf.clone(),ainf);
        assert_eq!(a1.clone()*rinf.clone(),am1);
        assert_eq!(ainf.clone()*rinf.clone(),a0);
        assert_eq!(am1.clone()*rinf.clone(),a1);
    }
    #[test]
    fn reflections_of_projective_one_points_in_f3() {
        let am1 = ProjOnePoint::new(Finite3::new(1),Finite3::new(-1));
        let a0 = ProjOnePoint::new(Finite3::new(1),Finite3::new(0));
        let a1 = ProjOnePoint::new(Finite3::new(1),Finite3::new(1));
        let ainf = ProjOnePoint::new(Finite3::new(0),Finite3::new(1));

        let s0 = wildmath::projectivegeo::Reflection::new(Finite3::new(1),Finite3::new(0));
        let s1 = wildmath::projectivegeo::Reflection::new(Finite3::new(1),Finite3::new(1));
        let sm1 = wildmath::projectivegeo::Reflection::new(Finite3::new(1),Finite3::new(-1));
        let sinf = wildmath::projectivegeo::Reflection::new(Finite3::new(0),Finite3::new(1));

        assert_eq!(a0.clone()*s0.clone(),a0);
        assert_eq!(a1.clone()*s0.clone(),am1);
        assert_eq!(ainf.clone()*s0.clone(),ainf);
        assert_eq!(am1.clone()*s0.clone(),a1);

        assert_eq!(a0.clone()*s1.clone(),a1);
        assert_eq!(a1.clone()*s1.clone(),a0);
        assert_eq!(ainf.clone()*s1.clone(),am1);
        assert_eq!(am1.clone()*s1.clone(),ainf);

        assert_eq!(a0.clone()*sm1.clone(),am1);
        assert_eq!(a1.clone()*sm1.clone(),ainf);
        assert_eq!(ainf.clone()*sm1.clone(),a1);
        assert_eq!(am1.clone()*sm1.clone(),a0);

        assert_eq!(a0.clone()*sinf.clone(),ainf);
        assert_eq!(a1.clone()*sinf.clone(),a1);
        assert_eq!(ainf.clone()*sinf.clone(),a0);
        assert_eq!(am1.clone()*sinf.clone(),am1);
    }
    #[test]
    #[should_panic]
    fn quadrance_of_projective_null_point_in_f5() {
        let am2 = ProjOnePoint::new(Finite5::new(1),Finite5::new(-2));
        let _q = am2.quadrance(&am2);
    }
    #[test]
    fn quadrance_of_projective_one_points_in_f5() {
        let am1 = ProjOnePoint::new(Finite5::new(1),Finite5::new(-1));
        let a0 = ProjOnePoint::new(Finite5::new(1),Finite5::new(0));
        let a1 = ProjOnePoint::new(Finite5::new(1),Finite5::new(1));
        let ainf = ProjOnePoint::new(Finite5::new(0),Finite5::new(1));

        let half = Finite5::new(1)/Finite5::new(2);
        let one = Finite5::new(1);
        let zero = Finite5::new(0);

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
    #[test]
    fn rotations_of_projective_one_points_in_f5() {
        let am1 = ProjOnePoint::new(Finite5::new(1),Finite5::new(-1));
        let a0 = ProjOnePoint::new(Finite5::new(1),Finite5::new(0));
        let a1 = ProjOnePoint::new(Finite5::new(1),Finite5::new(1));
        let a2 = ProjOnePoint::new(Finite5::new(1),Finite5::new(2));
        let ainf = ProjOnePoint::new(Finite5::new(0),Finite5::new(1));

        let r2 = wildmath::projectivegeo::Rotation::new(Finite5::new(1),Finite5::new(2));

        assert_eq!(am1.clone()*r2.clone(),a2);
        assert_eq!(a0.clone()*r2.clone(),a2);
        assert_eq!(a1.clone()*r2.clone(),a2);
        assert_eq!(a2.clone()*r2.clone(),a2);
        assert_eq!(ainf.clone()*r2.clone(),a2);
    }
    #[test]
    fn proof_multiplication_theorem_for_half_slopes() {
        let pone = PolyNumber::new(vec![PolyNumber::new(vec![Ratio::<i128>::new(1,1)])]);
        let ph = PolyNumber::new(vec![PolyNumber::new(vec![Ratio::new(0,1),Ratio::new(1,1)])]);
        let pk = PolyNumber::new(vec![PolyNumber::new(vec![Ratio::new(0,1)]),
                                     PolyNumber::new(vec![Ratio::new(1,1)])]);
        let one = PolyRatio::new(pone.clone(),pone.clone());
        let h = PolyRatio::new(ph.clone(),pone.clone());
        let k = PolyRatio::new(pk.clone(),pone.clone());
        let hk = (h.clone()+k.clone())/(one - h.clone()*k.clone());
        let eh = half_slope(h.clone());
        let ek = half_slope(k.clone());
        let ehk = half_slope(hk.clone());
        assert_eq!(eh*ek,ehk);
    }
    #[test]
    fn proof_multiplications_of_rotations_and_reflections_with_half_slopes() {
        let pone = PolyNumber::new(vec![PolyNumber::new(vec![Ratio::<i128>::new(1,1)])]);
        let ph = PolyNumber::new(vec![PolyNumber::new(vec![Ratio::new(0,1),Ratio::new(1,1)])]);
        let pk = PolyNumber::new(vec![PolyNumber::new(vec![Ratio::new(0,1)]),
                                     PolyNumber::new(vec![Ratio::new(1,1)])]);
        let one = PolyRatio::new(pone.clone(),pone.clone());

        let h = PolyRatio::new(ph.clone(),pone.clone());
        let k = PolyRatio::new(pk.clone(),pone.clone());
        let hk1 = (k.clone()-h.clone())/(one.clone() + h.clone()*k.clone());
        let hk2 = (k.clone()+h.clone())/(one.clone() - h.clone()*k.clone());

        let eh = half_slope(h.clone());
        let ek = half_slope(k.clone());
        let ehk1 = half_slope(hk1.clone());
        let ehk2 = half_slope(hk2.clone());

        let vh = TwoVector::new0e(eh);
        let vk = TwoVector::new0e(ek);
        let vhk1 = TwoVector::new0e(ehk1);
        let vhk2 = TwoVector::new0e(ehk2);

        let sh = wildmath::affinegeo2d::Reflection::new(vh.clone());
        let sk = wildmath::affinegeo2d::Reflection::new(vk.clone());
        let rh = wildmath::affinegeo2d::Rotation::new(vh.clone());
        let rk = wildmath::affinegeo2d::Rotation::new(vk.clone());

        assert_eq!(sh.clone()*sk.clone(),wildmath::affinegeo2d::Rotation::new(vhk1.clone()));
        assert_eq!(rh.clone()*rk.clone(),wildmath::affinegeo2d::Rotation::new(vhk2.clone()));
        assert_eq!(sh.clone()*rk.clone(),wildmath::affinegeo2d::Reflection::new(vhk2.clone()));
        assert_eq!(rh.clone()*sk.clone(),wildmath::affinegeo2d::Reflection::new(vhk1.clone()));
    }
    #[test]
    fn proof_addition_theorem_for_circle_sum() {
        let pone = create_polynumber_one!(h1,h2,h3; Ratio::<i32>);
        let ph1 = create_polynumber_var!(h1; h1,h2,h3; Ratio::<i32>);
        let ph2 = create_polynumber_var!(h2; h1,h2,h3; Ratio::<i32>);
        let ph3 = create_polynumber_var!(h3; h1,h2,h3; Ratio::<i32>);

        let one = PolyRatio::new(pone.clone(),pone.clone());
        let h1 = PolyRatio::new(ph1,pone.clone());
        let h2 = PolyRatio::new(ph2,pone.clone());
        let h3 = PolyRatio::new(ph3,pone.clone());

        let h123 = (h1.clone() + h2.clone() + h3.clone() - h1.clone()*h2.clone()*h3.clone())
                    /(one - (h1.clone()*h2.clone() + h2.clone()*h3.clone() + h1.clone()*h3.clone()));
        assert_eq!(circle_sum(circle_sum(h1.clone(),h2.clone()),h3.clone()),h123);
    }
    #[test]
    fn proof_diagnoals_of_parallelogram_bisect_each_other() {
        let polone = create_polynumber_one!(ax,ay,bx,by,dx,dy; Ratio::<i32>);
        let polax = create_polynumber_var!(ax; ax,ay,bx,by,dx,dy; Ratio::<i32>);
        let polay = create_polynumber_var!(ay; ax,ay,bx,by,dx,dy; Ratio::<i32>);
        let polbx = create_polynumber_var!(bx; ax,ay,bx,by,dx,dy; Ratio::<i32>);
        let polby = create_polynumber_var!(by; ax,ay,bx,by,dx,dy; Ratio::<i32>);
        let poldx = create_polynumber_var!(dx; ax,ay,bx,by,dx,dy; Ratio::<i32>);
        let poldy = create_polynumber_var!(dy; ax,ay,bx,by,dx,dy; Ratio::<i32>);

        let ax = PolyRatio::new(polax,polone.clone());
        let ay = PolyRatio::new(polay,polone.clone());
        let bx = PolyRatio::new(polbx,polone.clone());
        let by = PolyRatio::new(polby,polone.clone());
        let dx = PolyRatio::new(poldx,polone.clone());
        let dy = PolyRatio::new(poldy,polone.clone());

        let pa = TwoPoint::new(ax,ay);
        let pb = TwoPoint::new(bx,by);
        let pd = TwoPoint::new(dx,dy);

        let vv = TwoVector::newse(pa.clone(),pb.clone());
        let vu = TwoVector::newse(pa.clone(),pd.clone());

        let pc = pa.clone() + vu.clone() + vv.clone();

        let lac = TwoLine::newpv(&pa,&(vu.clone()+vv.clone()));
        let lbd = TwoLine::newpv(&pb,&(vu.clone()-vv.clone()));

        let pe = lac.meet(&lbd);

        let vae = TwoVector::newse(pa.clone(),pe.clone());
        let vec = TwoVector::newse(pe.clone(),pc.clone());
        let vbe = TwoVector::newse(pb.clone(),pe.clone());
        let ved = TwoVector::newse(pe.clone(),pd.clone());

        assert_eq!(vae,vec);
        assert_eq!(vbe,ved);
    }
    #[test]
    fn proof_medians_of_a_triangle() {
        let polone = create_polynumber_one!(ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polax = create_polynumber_var!(ax; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polay = create_polynumber_var!(ay; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polbx = create_polynumber_var!(bx; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polby = create_polynumber_var!(by; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polcx = create_polynumber_var!(dx; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polcy = create_polynumber_var!(dy; ax,ay,bx,by,dx,dy; Ratio::<i64>);

        let half = PolyRatio::new(polone.clone(),polone.clone()+polone.clone());
        let ax = PolyRatio::new(polax,polone.clone());
        let ay = PolyRatio::new(polay,polone.clone());
        let bx = PolyRatio::new(polbx,polone.clone());
        let by = PolyRatio::new(polby,polone.clone());
        let cx = PolyRatio::new(polcx,polone.clone());
        let cy = PolyRatio::new(polcy,polone.clone());

        let pa = TwoPoint::new(ax,ay);
        let pb = TwoPoint::new(bx,by);
        let pc = TwoPoint::new(cx,cy);

        let vu = TwoVector::newse(pa.clone(),pb.clone());
        let vv = TwoVector::newse(pa.clone(),pc.clone());

        let pd = pa.clone() + vu.clone()*half.clone();
        let pe = pa.clone() + vv.clone()*half.clone();

        let lbe = pb.join(&pe);
        let lcd = pc.join(&pd);

        let pg = lbe.meet(&lcd);

        let vbg = TwoVector::newse(pb.clone(),pg.clone());
        let vge = TwoVector::newse(pg.clone(),pe.clone());

        let vcg = TwoVector::newse(pc.clone(),pg.clone());
        let vgd = TwoVector::newse(pg.clone(),pd.clone());

        assert_eq!(vbg*half.clone(),vge);
        assert_eq!(vcg*half.clone(),vgd);
    }
    #[test]
    fn proof_vector_through_midpoints_parallel_third_side() {
        let polone = create_polynumber_one!(ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polax = create_polynumber_var!(ax; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polay = create_polynumber_var!(ay; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polbx = create_polynumber_var!(bx; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polby = create_polynumber_var!(by; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polcx = create_polynumber_var!(dx; ax,ay,bx,by,dx,dy; Ratio::<i64>);
        let polcy = create_polynumber_var!(dy; ax,ay,bx,by,dx,dy; Ratio::<i64>);

        let half = PolyRatio::new(polone.clone(),polone.clone()+polone.clone());
        let ax = PolyRatio::new(polax,polone.clone());
        let ay = PolyRatio::new(polay,polone.clone());
        let bx = PolyRatio::new(polbx,polone.clone());
        let by = PolyRatio::new(polby,polone.clone());
        let cx = PolyRatio::new(polcx,polone.clone());
        let cy = PolyRatio::new(polcy,polone.clone());

        let pa = TwoPoint::new(ax,ay);
        let pb = TwoPoint::new(bx,by);
        let pc = TwoPoint::new(cx,cy);

        let vu = TwoVector::newse(pa.clone(),pb.clone());
        let vv = TwoVector::newse(pa.clone(),pc.clone());
        let vw = TwoVector::newse(pb.clone(),pc.clone());

        let pd = pa.clone() + vu.clone()*half.clone();
        let pe = pa.clone() + vv.clone()*half.clone();

        let vde = TwoVector::newse(pd.clone(),pe.clone());

        assert_eq!(vw*half,vde);
    }
    #[test]
    fn proof_varignons_theorem() {
        let polone = create_polynumber_one!(ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let polax = create_polynumber_var!(ax; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let polay = create_polynumber_var!(ay; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let polbx = create_polynumber_var!(bx; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let polby = create_polynumber_var!(by; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let polcx = create_polynumber_var!(cx; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let polcy = create_polynumber_var!(cy; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let poldx = create_polynumber_var!(dx; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);
        let poldy = create_polynumber_var!(dy; ax,ay,bx,by,cx,cy,dx,dy; Ratio::<i32>);

        let half = PolyRatio::new(polone.clone(),polone.clone()+polone.clone());
        let ax = PolyRatio::new(polax,polone.clone());
        let ay = PolyRatio::new(polay,polone.clone());
        let bx = PolyRatio::new(polbx,polone.clone());
        let by = PolyRatio::new(polby,polone.clone());
        let cx = PolyRatio::new(polcx,polone.clone());
        let cy = PolyRatio::new(polcy,polone.clone());
        let dx = PolyRatio::new(poldx,polone.clone());
        let dy = PolyRatio::new(poldy,polone.clone());

        let pa = TwoPoint::new(ax.clone(),ay);
        let pb = TwoPoint::new(bx,by);
        let pc = TwoPoint::new(cx,cy);
        let pd = TwoPoint::new(dx,dy);

        let vab = TwoVector::newse(pa.clone(),pb.clone());
        let vbc = TwoVector::newse(pb.clone(),pc.clone());
        let vcd = TwoVector::newse(pc.clone(),pd.clone());
        let vda = TwoVector::newse(pd.clone(),pa.clone());

        let pmab = pa.clone() + vab.clone()*half.clone();
        let pmbc = pb.clone() + vbc.clone()*half.clone();
        let pmcd = pc.clone() + vcd.clone()*half.clone();
        let pmda = pd.clone() + vda.clone()*half.clone();

        let vmabmbc = TwoVector::newse(pmab.clone(),pmbc.clone());
        let vmbcmcd = TwoVector::newse(pmbc.clone(),pmcd.clone());
        let vmcdmda = TwoVector::newse(pmcd.clone(),pmda.clone());
        let vmdamab = TwoVector::newse(pmda.clone(),pmab.clone());

        assert!(vmabmbc.is_parallel(&vmcdmda));
        assert!(vmbcmcd.is_parallel(&vmdamab));
    }
}

