use crate::GeometricConstruction;

pub struct NavierStokesField {
    pub vesica_height: f64,
    pub natural_ratio: f64,
    pub velocity_field: Vec<Vec<Vec<(f64, f64, f64)>>>,
    pub pressure_field: Vec<Vec<Vec<f64>>>
}

impl NavierStokesField {
    pub fn new(construction: &GeometricConstruction) -> Self {
        let p1 = construction.points.p1.unwrap();
        let height = construction.distance_to_line(
            &p1,
            &construction.circle_a.center,
            &construction.circle_b.as_ref().unwrap().center
        );
        
        let vesica_ratio = height / construction.circle_a.radius;
        let size = vesica_ratio.round() as usize;
        
        Self {
            vesica_height: height,
            natural_ratio: construction.natural_ratio.unwrap(),
            velocity_field: vec![vec![vec![(0.0, 0.0, 0.0); size]; size]; size],
            pressure_field: vec![vec![vec![0.0; size]; size]; size]
        }
    }

    pub fn solve_step(&mut self, dt: f64) {
        let h = self.vesica_height;
        let ratio = self.natural_ratio;
        let viscosity = h;
        
        let mut next_velocity = self.velocity_field.clone();
        let mut next_pressure = self.pressure_field.clone();
        let size = self.velocity_field.len();
        
        for i in 1..size-1 {
            for j in 1..size-1 {
                for k in 1..size-1 {
                    let (u, v, w) = self.velocity_field[i][j][k];
                    
                    // Spatial derivatives scaled by vesica height
                    let du_dx = (self.velocity_field[i+1][j][k].0 - self.velocity_field[i-1][j][k].0) / (2.0 * h);
                    let du_dy = (self.velocity_field[i][j+1][k].0 - self.velocity_field[i][j-1][k].0) / (2.0 * h);
                    let du_dz = (self.velocity_field[i][j][k+1].0 - self.velocity_field[i][j][k-1].0) / (2.0 * h);
                    
                    let dv_dx = (self.velocity_field[i+1][j][k].1 - self.velocity_field[i-1][j][k].1) / (2.0 * h);
                    let dv_dy = (self.velocity_field[i][j+1][k].1 - self.velocity_field[i][j-1][k].1) / (2.0 * h);
                    let dv_dz = (self.velocity_field[i][j][k+1].1 - self.velocity_field[i][j][k-1].1) / (2.0 * h);
                    
                    let dw_dx = (self.velocity_field[i+1][j][k].2 - self.velocity_field[i-1][j][k].2) / (2.0 * h);
                    let dw_dy = (self.velocity_field[i][j+1][k].2 - self.velocity_field[i][j-1][k].2) / (2.0 * h);
                    let dw_dz = (self.velocity_field[i][j][k+1].2 - self.velocity_field[i][j][k-1].2) / (2.0 * h);

                    // Update velocities using natural ratio
                    next_velocity[i][j][k].0 = u + dt * (
                        -u * du_dx - v * du_dy - w * du_dz +
                        viscosity * (du_dx + du_dy + du_dz)
                    ) * ratio;
                    
                    next_velocity[i][j][k].1 = v + dt * (
                        -u * dv_dx - v * dv_dy - w * dv_dz +
                        viscosity * (dv_dx + dv_dy + dv_dz)
                    ) * ratio;
                    
                    next_velocity[i][j][k].2 = w + dt * (
                        -u * dw_dx - v * dw_dy - w * dw_dz +
                        viscosity * (dw_dx + dw_dy + dw_dz)
                    ) * ratio;

                    // Pressure emerges from geometric relationships
                    next_pressure[i][j][k] = self.pressure_field[i][j][k] + 
                        (du_dx + dv_dy + dw_dz) * h * ratio;
                }
            }
        }
        
        self.velocity_field = next_velocity;
        self.pressure_field = next_pressure;
    }

    pub fn verify_geometric_invariants(&self) -> bool {
        let size = self.velocity_field.len();
        
        for i in 1..size-1 {
            for j in 1..size-1 {
                for k in 1..size-1 {
                    // Pure geometric relationships through finite differences
                    let du_dx = (self.velocity_field[i+1][j][k].0 - self.velocity_field[i-1][j][k].0) / 2.0;
                    let dv_dy = (self.velocity_field[i][j+1][k].1 - self.velocity_field[i][j-1][k].1) / 2.0;
                    let dw_dz = (self.velocity_field[i][j][k+1].2 - self.velocity_field[i][j][k-1].2) / 2.0;
                    
                    // Incompressibility is a pure geometric invariant
                    if (du_dx + dv_dy + dw_dz).abs() > 1e-10 {
                        return false;
                    }
                }
            }
        }
        true
    }
}

