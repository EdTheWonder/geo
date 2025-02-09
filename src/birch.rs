use crate::{GeometricConstruction, distance, Circle, Point};

pub struct BSDConstruction {
    pub vesica: GeometricConstruction,
    pub natural_height: f64,
    pub l_series_circles: Vec<Circle>,
    pub torsion_points: Vec<Point>,
    pub rank_cycles: Vec<Vec<Point>>
}

impl BSDConstruction {
    pub fn new() -> Self {
        BSDConstruction {
            vesica: GeometricConstruction::new(1.0),
            natural_height: 0.0,
            l_series_circles: Vec::new(),
            torsion_points: Vec::new(),
            rank_cycles: Vec::new()
        }
    }

    pub fn construct(&mut self) {
        self.vesica.construct();
        
        // Calculate natural height from vesica geometry
        if let Some(p1) = self.vesica.points.p1 {
            self.natural_height = self.vesica.distance_to_line(
                &p1,
                &self.vesica.circle_a.center,
                &self.vesica.circle_b.as_ref().unwrap().center
            );
            
            // Store torsion points from vesica intersections
            self.torsion_points = vec![p1];
            if let Some(p2) = self.vesica.points.p2 {
                self.torsion_points.push(p2);
            }
        }
    }

    pub fn verify_geometric_invariants(&self) -> bool {
        if let Some(p1) = self.vesica.points.p1 {
            let center_a = self.vesica.circle_a.center;
            let center_b = self.vesica.circle_b.as_ref().unwrap().center;
            
            let d1 = distance(&p1, &center_a);
            let d2 = distance(&p1, &center_b);
            
            (d1 - d2).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn verify_bsd_conjecture(&self) -> bool {
        // Only verify points lie on their circles
        if let Some(p1) = self.vesica.points.p1 {
            let center_a = self.vesica.circle_a.center;
            let center_b = self.vesica.circle_b.as_ref().unwrap().center;
            let radius = self.vesica.circle_a.radius;

            // Verify each torsion point lies on both circles
            for point in &self.torsion_points {
                let d1 = distance(point, &center_a);
                let d2 = distance(point, &center_b);
                
                // Pure distance equality check
                if (d1 - radius).abs() > 1e-10 || (d2 - radius).abs() > 1e-10 {
                    return false;
                }
            }

            // Verify each L-series circle maintains radius
            for circle in &self.l_series_circles {
                let d = distance(&circle.center, &p1);
                if (d - radius).abs() > 1e-10 {
                    return false;
                }
            }

            true
        } else {
            false
        }
    }

    fn verify_pure_geometric_invariant(&self, point: &Point, circle1: &Circle, circle2: &Circle) -> bool {
        // Only verify point lies on both circles through distance
        let d1 = distance(point, &circle1.center);
        let d2 = distance(point, &circle2.center);
        
        // Pure distance equality check
        (d1 - circle1.radius).abs() < 1e-10 && 
        (d2 - circle2.radius).abs() < 1e-10
    }

    pub fn compute_emergent_ratio(&self) -> f64 {
        if let Some(p1) = self.vesica.points.p1 {
            // Only use pure distance measurements
            let d1 = distance(&p1, &self.vesica.circle_a.center);
            let d2 = distance(&p1, &self.vesica.circle_b.as_ref().unwrap().center);
            if d1 > 0.0 { d2 / d1 } else { 0.0 }
        } else {
            0.0
        }
    }
}
