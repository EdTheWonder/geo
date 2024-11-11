use super::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pure_geometric_invariants() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        // Pure geometric verification using only compass operations
        let verify_circle_points = |p: &Point, c1: &Circle, c2: &Circle| -> bool {
            let d1 = distance(&p, &c1.center);
            let d2 = distance(&p, &c2.center);
            (d1 - c1.radius).abs() < 1e-10 && (d2 - c2.radius).abs() < 1e-10
        };
        
        let verify_square_distances = |a: &Point, b: &Point, c: &Point| -> bool {
            let side = distance(a, b);
            let next_side = distance(b, c);
            let ratio = if next_side != 0.0 { side / next_side } else { 0.0 };
            (ratio - 1.0).abs() < 1e-10
        };
        
        // Vesica Piscis verification
        if let Some(p1) = construction.points.p1 {
            assert!(verify_circle_points(
                &p1,
                &construction.circle_a,
                construction.circle_b.as_ref().unwrap()
            ));
        }
        
        // Square verification through pure distance relationships
        if let (Some(a), Some(b), Some(c)) = (
            construction.points.a,
            construction.points.b,
            construction.points.c
        ) {
            // Verify through direct distance comparisons only
            assert!(verify_square_distances(&a, &b, &c));
            
            // Document fundamental ratio
            let side = distance(&a, &b);
            let diagonal = distance(&a, &c);
            println!("\nFundamental Ratio (side:diagonal): {}", side/diagonal);
        }
    }

    #[test]
    fn test_circle_intersection_angles() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(p1), Some(p2)) = (construction.points.p1, construction.points.p2) {
            // Verify P1 and P2 are equidistant from both circle centers
            let d1_center_a = distance(&p1, &construction.circle_a.center);
            let d1_center_b = distance(&p1, &construction.circle_b.as_ref().unwrap().center);
            let d2_center_a = distance(&p2, &construction.circle_a.center);
            let d2_center_b = distance(&p2, &construction.circle_b.as_ref().unwrap().center);
            
            // Both points should be on both circles (distance = radius)
            assert!((d1_center_a - construction.circle_a.radius).abs() < 1e-10);
            assert!((d1_center_b - construction.circle_b.as_ref().unwrap().radius).abs() < 1e-10);
            assert!((d2_center_a - construction.circle_a.radius).abs() < 1e-10);
            assert!((d2_center_b - construction.circle_b.as_ref().unwrap().radius).abs() < 1e-10);
        }
    }

    #[test]
    fn test_square_properties() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(a), Some(b), Some(c), Some(d)) = (
            construction.points.a,
            construction.points.b,
            construction.points.c,
            construction.points.d
        ) {
            // Verify square properties through pure distance relationships
            let side_ab = distance(&a, &b);
            let side_bc = distance(&b, &c);
            let side_cd = distance(&c, &d);
            let side_da = distance(&d, &a);
            
            // All sides should be equal
            assert!((side_ab - side_bc).abs() < 1e-10);
            assert!((side_bc - side_cd).abs() < 1e-10);
            assert!((side_cd - side_da).abs() < 1e-10);
            
            // Diagonals should be equal
            let diag_ac = distance(&a, &c);
            let diag_bd = distance(&b, &d);
            assert!((diag_ac - diag_bd).abs() < 1e-10);
            
            // Pythagorean relationship
            let sides = side_ab.powi(2) + side_bc.powi(2);
            let diagonal = diag_ac.powi(2);
            assert!((sides - diagonal).abs() < 1e-10);
        }
    }

    #[test]
    fn test_construction_invariants() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        // Test if all points lie on their respective circles
        if let Some(circle_b) = &construction.circle_b {
            for point in [
                construction.points.p1,
                construction.points.p3,
                construction.points.p4,
                construction.points.p5,
                construction.points.p6
            ].iter().flatten() {
                assert!(construction.point_on_circle(point, &construction.circle_a) ||
                       construction.point_on_circle(point, circle_b));
            }
        }
    }

    #[test]
    fn test_circle_square_equality() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(a), Some(b), Some(p1)) = (
            construction.points.a,
            construction.points.b,
            construction.points.p1
        ) {
            let vesica_height = construction.distance_to_line(
                &p1,
                &construction.circle_a.center,
                &construction.circle_b.as_ref().unwrap().center
            );
            
            let square_side = distance(&a, &b);
            
            // Compare squared values to avoid floating point errors
            let side_squared = square_side * square_side;
            let height_squared = vesica_height * vesica_height * 4.0;
            assert!((side_squared - height_squared).abs() < 1e-10);
        }
    }
}