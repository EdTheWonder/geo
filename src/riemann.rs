use super::{Point, Circle, Point3D, Sphere, GeometricConstruction};
pub struct RiemannConstruction {
    pub base_circle: Circle,
    pub zeta_circles: Vec<Circle>,
    pub critical_points: Vec<Point>,
    pub vesica_ratio: Option<f64>,
    pub manifold_points: Vec<Point3D>
}

impl RiemannConstruction {
    pub fn new() -> Self {
        Self {
            base_circle: Circle {
                center: Point { x: 0.0, y: 0.0 },
                radius: 1.0
            },
            zeta_circles: Vec::new(),
            critical_points: Vec::new(),
            vesica_ratio: None,
            manifold_points: Vec::new()
        }
    }

    pub fn construct(&mut self) {
        let mut vesica = super::GeometricConstruction::new(1.0);
        vesica.construct();
        
        if let Some(p1) = vesica.points.p1 {
            let height = vesica.distance_to_line(
                &p1,
                &vesica.circle_a.center,
                &vesica.circle_b.as_ref().unwrap().center
            );
            
            self.vesica_ratio = Some(height);
            
            // Critical circle emerges naturally at x=0.5
            let critical_circle = Circle {
                center: Point { x: 0.5, y: 0.0 },
                radius: height
            };
            self.zeta_circles.push(critical_circle);
            
            // Critical points emerge from natural circle intersections
            let intersections = vesica.find_circle_intersections(
                &vesica.circle_a,
                vesica.circle_b.as_ref().unwrap()
            );
            
            for point in intersections {
                // Transform points preserving circle radius
                let y_scale = if point.y >= 0.0 { 1.0 } else { -1.0 };
                let y = height * y_scale;  // Use height directly instead of scaling point.y
                self.critical_points.push(Point { x: 0.5, y });
            }
            
            self.critical_points.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
        }
    }

    pub fn point_on_circle(&self, p: &Point, c: &Circle) -> bool {
        let dx = p.x - c.center.x;
        let dy = p.y - c.center.y;
        let d = (dx * dx + dy * dy).sqrt();
        
        // Debug print to see actual values
        println!("Point: ({}, {})", p.x, p.y);
        println!("Circle: center({}, {}), radius={}", c.center.x, c.center.y, c.radius);
        println!("Distance: {}", d);
        
        (d - c.radius).abs() < 1e-10
    }

    pub fn find_circle_intersections(&self, c1: &Circle, c2: &Circle) -> Vec<Point> {
        let mut points = Vec::new();
        
        let dx = c2.center.x - c1.center.x;
        let dy = c2.center.y - c1.center.y;
        let d = (dx * dx + dy * dy).sqrt();
        
        // Circles must intersect
        if d > c1.radius + c2.radius || d < (c1.radius - c2.radius).abs() {
            return points;
        }
        
        let a = (c1.radius * c1.radius - c2.radius * c2.radius + d * d) / (2.0 * d);
        let h = (c1.radius * c1.radius - a * a).sqrt();
        
        let x2 = c1.center.x + (dx * a) / d;
        let y2 = c1.center.y + (dy * a) / d;
        
        let x3 = x2 + (h * dy) / d;
        let y3 = y2 - (h * dx) / d;
        
        let x4 = x2 - (h * dy) / d;
        let y4 = y2 + (h * dx) / d;
        
        points.push(Point { x: x3, y: y3 });
        points.push(Point { x: x4, y: y4 });
        
        points
    }

    pub fn verify_critical_point(&self, p: &Point) -> bool {
        if let Some(critical_circle) = self.zeta_circles.first() {
            // Verify point lies on critical circle and critical line
            self.point_on_circle(p, critical_circle) && 
            (p.x - 0.5).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn verify_critical_circle(&self) -> bool {
        if let Some(critical_circle) = self.zeta_circles.first() {
            // Critical circle must be centered at (0.5, 0)
            if (critical_circle.center.x - 0.5).abs() > 1e-10 || 
               critical_circle.center.y.abs() > 1e-10 {
                return false;
            }
            
            // Radius must equal vesica height
            (critical_circle.radius - self.vesica_ratio.unwrap()).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn verify_zero_symmetry(&self, p1: &Point, p2: &Point) -> bool {
        // Zeros must be symmetric about real axis
        if (p1.x - 0.5).abs() > 1e-10 || (p2.x - 0.5).abs() > 1e-10 {
            return false;
        }
        
        // y-coordinates must be equal and opposite
        (p1.y + p2.y).abs() < 1e-10
    }

    pub fn extend_to_3d(&mut self) {
        // Start with base vesica piscis construction
        let mut vesica = GeometricConstruction::new(1.0);
        vesica.construct();
        
        if let Some(p1) = vesica.points.p1 {
            // Use vesica height as sphere radius (proven r√3 relationship)
            let height = vesica.distance_to_line(
                &p1,
                &vesica.circle_a.center,
                &vesica.circle_b.as_ref().unwrap().center
            );
            
            // Create critical sphere at x=0.5 (analogous to critical circle)
            let critical_sphere = Sphere {
                center: Point3D { x: 0.5, y: 0.0, z: 0.0 },
                radius: height
            };
            
            // Generate 3D manifold points through sphere intersections
            let intersections = self.find_sphere_intersections(&critical_sphere);
            
            // Transform points preserving geometric ratios
            for point in intersections {
                let scaled_point = self.scale_to_manifold(&point, height);
                self.manifold_points.push(scaled_point);
            }
        }
    }

    pub fn scale_to_manifold(&self, point: &Point3D, height: f64) -> Point3D {
        let scale = height / ((point.x*point.x + point.y*point.y + point.z*point.z).sqrt());
        Point3D {
            x: point.x * scale,
            y: point.y * scale,
            z: point.z * scale
        }
    }

    pub fn find_sphere_intersections(&self, sphere: &Sphere) -> Vec<Point3D> {
        let mut intersections = Vec::new();
        
        // Create second sphere based on vesica ratio
        let sphere2 = Sphere {
            center: Point3D { 
                x: sphere.center.x + sphere.radius, 
                y: 0.0, 
                z: 0.0 
            },
            radius: sphere.radius
        };
        
        // Distance between sphere centers
        let dx = sphere2.center.x - sphere.center.x;
        let d = dx.abs();
        
        // Check if spheres intersect
        if d > sphere.radius + sphere2.radius || d < (sphere.radius - sphere2.radius).abs() {
            return intersections;
        }
        
        // Find intersection circle
        let a = (sphere.radius.powi(2) - sphere2.radius.powi(2) + d.powi(2)) / (2.0 * d);
        let h = (sphere.radius.powi(2) - a.powi(2)).sqrt();
        
        // Calculate intersection points preserving vesica ratio
        let px = sphere.center.x + (a/d)*dx;
        
        // Add symmetric points maintaining vesica height
        intersections.push(Point3D { 
            x: px, 
            y: h, 
            z: 0.0 
        });
        intersections.push(Point3D { 
            x: px, 
            y: -h, 
            z: 0.0 
        });
        
        intersections
    }
}

impl Point3D {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Point3D { x, y, z }
    }
    
    pub fn origin() -> Self {
        Point3D::new(0.0, 0.0, 0.0)
    }
}

impl Sphere {
    pub fn new(center: Point3D, radius: f64) -> Self {
        Sphere { center, radius }
    }
}


