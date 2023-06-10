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

use geo::*;
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

