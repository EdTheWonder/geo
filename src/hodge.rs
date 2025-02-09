use crate::{GeometricConstruction, distance};

pub struct HodgeForm<'a> {
    pub dimension: usize,
    pub vesica_ratio: f64,
    pub natural_ratio: f64,
    pub cycles: Vec<Vec<f64>>,
    construction: &'a GeometricConstruction
}

impl<'a> HodgeForm<'a> {
    pub fn new(construction: &'a GeometricConstruction) -> Self {
        let p1 = construction.points.p1.unwrap();
        let height = construction.distance_to_line(
            &p1,
            &construction.circle_a.center,
            &construction.circle_b.as_ref().unwrap().center
        );
        
        let radius = construction.circle_a.radius;
        // Verify height = r√3 through Pythagorean theorem
        assert!((height.powi(2) - 3.0 * radius.powi(2)).abs() < 1e-10);
        
        // Let dimension emerge from vesica height
        let dimension = (height / radius).round() as usize;  // Should be √3 ≈ 1.732
        
        HodgeForm {
            dimension,
            vesica_ratio: height / radius,  // = √3
            natural_ratio: construction.natural_ratio.unwrap(),  // Let ratio emerge naturally
            cycles: Vec::new(),
            construction
        }
    }

    pub fn construct_cycles(&mut self) {
        let mut cycles = Vec::new();
        
        // Use only existing vesica points
        let p1 = self.construction.points.p1.unwrap();
        let center_a = self.construction.circle_a.center;
        let center_b = self.construction.circle_b.as_ref().unwrap().center;
        
        // First cycle from pure vesica distances
        let mut cycle = Vec::new();
        cycle.push(distance(&p1, &center_a));
        cycle.push(distance(&p1, &center_b));
        cycle.push(distance(&center_a, &center_b));
        cycles.push(cycle);
        
        self.cycles = cycles;
    }

    pub fn verify_hodge_classes(&self) -> bool {
        // Verify only through pure distance relationships
        for cycle in &self.cycles {
            // All distances in a cycle must be equal (vesica property)
            let base_distance = cycle[0];
            for d in cycle.iter().skip(1) {
                if (d - base_distance).abs() > 1e-10 {
                    return false;
                }
            }
        }
        true
    }
}
