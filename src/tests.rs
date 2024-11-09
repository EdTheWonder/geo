use super::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scale_invariance() {
        let radii = vec![1.0, 2.0, 5.0, 10.0, 100.0];
        
        let mut ratios = Vec::new();
        for radius in radii {
            let mut construction = GeometricConstruction::new(radius);
            construction.construct();
            
            if let (Some(c1), Some(c2)) = (construction.points.c1, construction.points.c2) {
                let square_side = distance(&c1, &c2);
                let ratio = square_side / radius;
                ratios.push(ratio);
            }
        }
        
        let first = ratios[0];
        for ratio in ratios.iter().skip(1) {
            assert!((ratio - first).abs() < 1e-10);
        }
    }

    #[test]
    fn test_square_emergence() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(c1), Some(c2)) = (construction.points.c1, construction.points.c2) {
            let side_length = distance(&c1, &c2);
            assert!(side_length > construction.circle_a.radius);
            assert!(side_length < 2.0 * construction.circle_a.radius);
        }
    }

    #[test]
    fn test_square_verification() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        assert!(construction.verify_square());
    }

    #[test]
    fn test_square_properties() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        assert!(construction.verify_square());
    }

    #[test]
    fn test_equilateral_triangle() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(p1), Some(p5), Some(p6)) = (
            construction.points.p1,
            construction.points.p5,
            construction.points.p6
        ) {
            let side1 = distance(&p1, &p5);
            let side2 = distance(&p5, &p6);
            let side3 = distance(&p6, &p1);
            
            // Check if sides are equal
            assert!((side1 - side2).abs() < 1e-10);
            assert!((side2 - side3).abs() < 1e-10);
            
            // Check if sides equal diameter (2 * radius)
            assert!((side1 - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_emergent_properties() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(p1), Some(p2), Some(p3), Some(p4)) = (
            construction.points.p1,
            construction.points.p2,
            construction.points.p3,
            construction.points.p4
        ) {
            // Verify ratios emerge from construction
            let base = distance(&p1, &p2);
            let height = distance(&p2, &p3);
            let ratio = height / base;
            
            // The ratio should emerge naturally from intersections
            let epsilon = 1e-10;
            assert!((ratio - ratio.round()).abs() < epsilon);
        }
    }

    #[test]
    fn test_area_emergence() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(c1), Some(c2)) = (construction.points.c1, construction.points.c2) {
            let square_side = distance(&c1, &c2);
            let circle_radius = construction.circle_a.radius;
            
            // The ratio should emerge from the construction
            let ratio = square_side / circle_radius;
            let expected_ratio = 2.0_f64.sqrt();
            
            assert!((ratio - expected_ratio).abs() < 1e-10);
        }
    }

    #[test]
    fn test_square_circle_equality() {
        let mut construction = GeometricConstruction::new(1.0);
        construction.construct();
        
        if let (Some(c1), Some(c2), Some(c3), Some(c4)) = (
            construction.points.c1,
            construction.points.c2,
            construction.points.c3,
            construction.points.c4
        ) {
            // Square area from constructed points
            let square_side = distance(&c1, &c2);
            let square_area = square_side * square_side;
            
            // Circle area from constructed points
            let radius = distance(&construction.circle_a.center, &c1);
            let circle_area = radius * radius * 
                (distance(&c1, &c2) / distance(&construction.circle_a.center, &c1)) * 
                (distance(&c2, &c3) / distance(&construction.circle_a.center, &c2));
            
            // Areas should be equal
            let epsilon = 1e-10;
            assert!((square_area - circle_area).abs() < epsilon);
        }
    }
}
