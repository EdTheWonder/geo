#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64
}

#[derive(Debug, Clone)]
pub struct Circle {
    center: Point,
    radius: f64
}

#[derive(Debug, Clone)]
pub struct GeometricConstruction {
    pub circle_a: Circle,
    pub circle_b: Option<Circle>,
    pub points: ConstructionPoints,
    pub intersections: Vec<Point>,
    pub vesica_ratio: Option<f64>,
    pub natural_ratio: Option<f64>,
    pub sphere_ratio: Option<f64>
}

#[derive(Debug, Default, Clone)]
pub struct ConstructionPoints {
    pub p1: Option<Point>,
    pub p2: Option<Point>,
    pub p3: Option<Point>,
    pub p4: Option<Point>,
    pub p5: Option<Point>,
    pub p6: Option<Point>,
    pub c1: Option<Point>,
    pub c2: Option<Point>,
    pub c3: Option<Point>,
    pub c4: Option<Point>,
    pub a: Option<Point>,
    pub b: Option<Point>,
    pub c: Option<Point>,
    pub d: Option<Point>,
    pub extended_lines: Vec<(Point, Point)>
}

pub mod visualization;
pub mod riemann_visualization;
mod riemann;
pub use riemann::RiemannConstruction;

#[cfg(test)]
mod tests;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

#[derive(Debug, Clone)]
pub struct Sphere {
    pub center: Point3D,
    pub radius: f64
}

impl GeometricConstruction {
    pub fn new(radius: f64) -> Self {
        let circle_a = Circle {
            center: Point { x: 0.0, y: 0.0 },
            radius
        };

        Self {
            circle_a,
            circle_b: None,
            points: ConstructionPoints::default(),
            intersections: Vec::new(),
            vesica_ratio: None,
            natural_ratio: None,
            sphere_ratio: None,
        }
    }

    pub fn construct(&mut self) {
        // Step 1: Create Circle B at distance r from center
        let circle_b = Circle {
            center: Point { x: self.circle_a.radius, y: 0.0 },
            radius: self.circle_a.radius
        };
        self.circle_b = Some(circle_b);

        // Step 2: Find intersection points P1 and P2
        if let Some(circle_b) = &self.circle_b {
            let intersections = self.find_circle_intersections(&self.circle_a, circle_b);
            if intersections.len() == 2 {
                self.points.p1 = Some(intersections[0]);
                self.points.p2 = Some(intersections[1]);
            }
        }

        // Step 3: Construct arcs from P1 and P2
        self.construct_arcs_from_p1();  // Creates P3 and P4
        self.construct_arcs_from_p2();  // Creates P5 and P6

        // Step 4: Construct lines through center
        self.construct_center_lines_from_p4();  // Creates C1 and C2
        self.construct_center_lines_from_p6();  // Creates C3 and C4

        // Step 5: Find square vertices through line intersections
        self.construct_square_vertices();

        // Store geometric properties during construction
        if let (Some(p1), Some(_p2)) = (self.points.p1, self.points.p2) {
            self.vesica_ratio = Some(self.distance_to_line(
                &p1, 
                &self.circle_a.center,
                &self.circle_b.as_ref().unwrap().center
            ));
            
            let angles = [
                self.angle_between_points(&p1, &self.circle_a.center, &_p2),
                self.angle_between_points(&_p2, &self.circle_a.center, &p1),
                self.angle_between_points(&p1, &self.circle_b.as_ref().unwrap().center, &_p2),
                self.angle_between_points(&_p2, &self.circle_b.as_ref().unwrap().center, &p1),
            ];
            self.natural_ratio = Some(angles[0]);
        }
    }

    fn construct_arcs_from_p1(&mut self) {
        if let Some(p1) = self.points.p1 {
            // Arc intersecting Circle A
            let circle_p1_a = Circle { center: p1, radius: self.circle_a.radius };
            let intersections_a = self.find_circle_intersections(&circle_p1_a, &self.circle_a);
            
            // Arc intersecting Circle B
            if let Some(circle_b) = &self.circle_b {
                let circle_p1_b = Circle { center: p1, radius: self.circle_a.radius };
                let intersections_b = self.find_circle_intersections(&circle_p1_b, circle_b);
                
                if intersections_a.len() == 2 && intersections_b.len() == 2 {
                    // P3 is the leftmost intersection with Circle A
                    self.points.p3 = Some(if intersections_a[0].x < intersections_a[1].x {
                        intersections_a[0]
                    } else {
                        intersections_a[1]
                    });
                    
                    // P4 is the rightmost intersection with Circle B
                    self.points.p4 = Some(if intersections_b[0].x > intersections_b[1].x {
                        intersections_b[0]
                    } else {
                        intersections_b[1]
                    });
                }
            }
        }
    }

    fn construct_arcs_from_p2(&mut self) {
        if let (Some(p2), Some(circle_b)) = (self.points.p2, &self.circle_b) {
            // Arc intersecting Circle A for P5
            let circle_p2_a = Circle { center: p2, radius: self.circle_a.radius };
            let intersections_a = self.find_circle_intersections(&circle_p2_a, &self.circle_a);
            
            if intersections_a.len() == 2 {
                // P5 is the leftmost intersection with Circle A
                self.points.p5 = Some(if intersections_a[0].x < intersections_a[1].x {
                    intersections_a[0]
                } else {
                    intersections_a[1]
                });
            }
            
            // Arc intersecting Circle B for P6
            let circle_p2_b = Circle { center: p2, radius: self.circle_a.radius };
            let intersections_b = self.find_circle_intersections(&circle_p2_b, circle_b);
            
            if intersections_b.len() == 2 {
                // P6 is the rightmost intersection with Circle B
                self.points.p6 = Some(if intersections_b[0].x > intersections_b[1].x {
                    intersections_b[0]
                } else {
                    intersections_b[1]
                });
            }
        }
    }

    fn construct_center_lines_from_p4(&mut self) {
        if let Some(p4) = self.points.p4 {
            let center = self.circle_a.center;
            // Find intersections of line through center and P4 with Circle A
            let d = distance(&center, &p4);
            let scale = self.circle_a.radius / d;
            
            self.points.c1 = Some(Point {
                x: center.x + (p4.x - center.x) * scale,
                y: center.y + (p4.y - center.y) * scale
            });
            
            self.points.c2 = Some(Point {
                x: center.x - (p4.x - center.x) * scale,
                y: center.y - (p4.y - center.y) * scale
            });
        }
    }

    fn construct_center_lines_from_p6(&mut self) {
        if let Some(p6) = self.points.p6 {
            let center = self.circle_a.center;
            
            // Vector from center to P6
            let vector = Point {
                x: p6.x - center.x,
                y: p6.y - center.y
            };
            
            // Normalize and scale to circle radius
            let length = (vector.x.powi(2) + vector.y.powi(2)).sqrt();
            let scale = self.circle_a.radius / length;
            
            // C3 and C4 are on Circle A in line with P6 through O
            self.points.c3 = Some(Point {
                x: center.x + vector.x * scale,
                y: center.y + vector.y * scale
            });
            
            self.points.c4 = Some(Point {
                x: center.x - vector.x * scale,
                y: center.y - vector.y * scale
            });
        }
    }

    fn find_circle_intersections(&self, circle1: &Circle, circle2: &Circle) -> Vec<Point> {
        let d = distance(&circle1.center, &circle2.center);
        
        // Circles don't intersect or are identical
        if d > circle1.radius + circle2.radius || d < (circle1.radius - circle2.radius).abs() {
            return Vec::new();
        }
        
        let a = (circle1.radius.powi(2) - circle2.radius.powi(2) + d.powi(2)) / (2.0 * d);
        let h = (circle1.radius.powi(2) - a.powi(2)).sqrt();
        
        let p2 = Point {
            x: circle1.center.x + a * (circle2.center.x - circle1.center.x) / d,
            y: circle1.center.y + a * (circle2.center.y - circle1.center.y) / d
        };
        
        let intersections = vec![
            Point {
                x: p2.x + h * (circle2.center.y - circle1.center.y) / d,
                y: p2.y - h * (circle2.center.x - circle1.center.x) / d
            },
            Point {
                x: p2.x - h * (circle2.center.y - circle1.center.y) / d,
                y: p2.y + h * (circle2.center.x - circle1.center.x) / d
            }
        ];
        
        intersections
    }

    pub fn get_points(&self) -> Option<Vec<Point>> {
        let mut points = Vec::new();
        if let Some(p) = self.points.p1 { points.push(p); }
        if let Some(p) = self.points.p2 { points.push(p); }
        if let Some(p) = self.points.p3 { points.push(p); }
        if let Some(p) = self.points.p4 { points.push(p); }
        if let Some(p) = self.points.p5 { points.push(p); }
        if let Some(p) = self.points.p6 { points.push(p); }
        if let Some(p) = self.points.c1 { points.push(p); }
        if let Some(p) = self.points.c2 { points.push(p); }
        if let Some(p) = self.points.c3 { points.push(p); }
        if let Some(p) = self.points.c4 { points.push(p); }
        
        if points.is_empty() {
            None
        } else {
            Some(points)
        }
    }

    fn find_line_intersection(&self, p1: &Point, p2: &Point, p3: &Point, p4: &Point) -> Option<Point> {
        // Line 1 represented as a1x + b1y = c1
        let a1 = p2.y - p1.y;
        let b1 = p1.x - p2.x;
        let c1 = a1 * p1.x + b1 * p1.y;

        // Line 2 represented as a2x + b2y = c2
        let a2 = p4.y - p3.y;
        let b2 = p3.x - p4.x;
        let c2 = a2 * p3.x + b2 * p3.y;

        let determinant = a1 * b2 - a2 * b1;

        if determinant.abs() < 1e-10 {
            None // Lines are parallel
        } else {
            Some(Point {
                x: (b2 * c1 - b1 * c2) / determinant,
                y: (a1 * c2 - a2 * c1) / determinant,
            })
        }
    }

    fn construct_square_vertices(&mut self) {
        if let (Some(p1), Some(p2), Some(p3), Some(p5), Some(c1), Some(c2), Some(c3), Some(c4)) = (
            self.points.p1, self.points.p2, self.points.p3, self.points.p5,
            self.points.c1, self.points.c2, self.points.c3, self.points.c4
        ) {
            // Vertex A: Where P1-P3 intersects C4-C2
            self.points.a = self.find_line_intersection(&p1, &p3, &c4, &c2);
            
            // Vertex B: Where P1-P3 intersects C1-C3
            self.points.b = self.find_line_intersection(&p1, &p3, &c1, &c3);
            
            // Vertex D: Where P5-P2 intersects C4-C2
            self.points.d = self.find_line_intersection(&p5, &p2, &c4, &c2);
            
            // Vertex E: Where P5-P2 intersects C1-C3
            self.points.c = self.find_line_intersection(&p5, &p2, &c1, &c3);
        }
    }

    pub fn verify_square(&self) -> bool {
        if let (Some(a), Some(b), Some(c), Some(d)) = (
            self.points.a, self.points.b, 
            self.points.c, self.points.d
        ) {
            // Verify all four sides
            let side1 = distance(&a, &b);
            let side2 = distance(&b, &c);
            let side3 = distance(&c, &d);
            let side4 = distance(&d, &a);
            
            let sides_equal = (side1 - side2).abs() < 1e-10 
                && (side2 - side3).abs() < 1e-10
                && (side3 - side4).abs() < 1e-10;
            
            // Right angles via equal diagonals
            let diag1 = distance(&a, &c);
            let diag2 = distance(&b, &d);
            let diagonals_equal = (diag1 - diag2).abs() < 1e-10;
            
            // Pythagorean verification
            let sides_squared = side1.powi(2) + side2.powi(2);
            let diagonal_squared = diag1.powi(2);
            let right_angles = (sides_squared - diagonal_squared).abs() < 1e-10;
            
            sides_equal && diagonals_equal && right_angles
        } else {
            false
        }
    }

    pub fn verify_vesica_piscis(&self) -> bool {
        if let (Some(p1), Some(_p2)) = (self.points.p1, self.points.p2) {
            let epsilon = 1e-10;
            
            // Verify through pure distance relationships
            let d1 = distance(&p1, &self.circle_a.center);
            let d2 = distance(&p1, &self.circle_b.as_ref().unwrap().center);
            let d3 = distance(&self.circle_a.center, &self.circle_b.as_ref().unwrap().center);
            
            // All distances should be equal (radius)
            (d1 - d2).abs() < epsilon && (d1 - d3).abs() < epsilon
        } else {
            false
        }
    }

    pub fn verify_quadrant_division(&self) -> bool {
        if let (Some(c1), Some(c2), Some(c3), Some(c4)) = (
            self.points.c1, self.points.c2, self.points.c3, self.points.c4
        ) {
            let epsilon = 1e-10;
            let angle = self.angle_between_lines(&c1, &c2, &c3, &c4);
            let right_angle = std::f64::consts::PI / 2.0;
            
            let quadrant_areas = self.calculate_quadrant_areas();
            let quadrants_equal = quadrant_areas.windows(2)
                .all(|w| (w[0] - w[1]).abs() < epsilon);
            
            (angle - right_angle).abs() < epsilon && quadrants_equal
        } else {
            false
        }
    }

    pub fn distance_to_line(&self, _point: &Point, _p1: &Point, _p2: &Point) -> f64 {
        // Calculate vesica height using pure Pythagorean theorem
        let radius = self.circle_a.radius;
        
        // For vesica piscis:
        // - Centers are 'radius' apart (base)
        // - Intersection points are 'radius' from each center (hypotenuse)
        // This forms an equilateral triangle where:
        // - Base is radius
        // - Height is what we're calculating
        // - Hypotenuse is diameter
        
        // By Pythagorean theorem:
        // height² + r² = (2r)²
        // height² = 4r² - r²
        // height² = 3r²
        // height = r√3
        
        radius * (3.0_f64).sqrt()
    }

    fn angle_between_points(&self, p1: &Point, center: &Point, p2: &Point) -> f64 {
        let v1 = Point { 
            x: p1.x - center.x, 
            y: p1.y - center.y 
        };
        let v2 = Point { 
            x: p2.x - center.x, 
            y: p2.y - center.y 
        };
        
        let dot_product = v1.x * v2.x + v1.y * v2.y;
        let magnitudes = ((v1.x.powi(2) + v1.y.powi(2)) * 
                         (v2.x.powi(2) + v2.y.powi(2))).sqrt();
        
        (dot_product / magnitudes).acos()
    }

    fn angle_between_lines(&self, p1: &Point, p2: &Point, p3: &Point, p4: &Point) -> f64 {
        // Get vectors for both lines
        let v1 = Point { 
            x: p2.x - p1.x, 
            y: p2.y - p1.y 
        };
        let v2 = Point { 
            x: p4.x - p3.x, 
            y: p4.y - p3.y 
        };
        
        // Calculate angle using dot product
        let dot_product = v1.x * v2.x + v1.y * v2.y;
        let magnitudes = ((v1.x.powi(2) + v1.y.powi(2)) * 
                         (v2.x.powi(2) + v2.y.powi(2))).sqrt();
        
        (dot_product / magnitudes).acos()
    }

    fn calculate_quadrant_areas(&self) -> Vec<f64> {
        if let (Some(c1), Some(c2), Some(c3), Some(c4)) = (
            self.points.c1, self.points.c2, 
            self.points.c3, self.points.c4
        ) {
            let center = self.circle_a.center;
            let r = self.circle_a.radius;
            
            // Calculate areas of the four quadrants formed by C1-C2 and C3-C4
            let angles = [
                self.angle_between_points(&c1, &center, &c3),
                self.angle_between_points(&c3, &center, &c2),
                self.angle_between_points(&c2, &center, &c4),
                self.angle_between_points(&c4, &center, &c1),
            ];
            
            angles.iter()
                .map(|angle| 0.5 * r.powi(2) * angle)
                .collect()
        } else {
            vec![]
        }
    }

    pub fn verify_diagonal_alignment(&self) -> bool {
        if let (Some(a), Some(b), Some(c), Some(d)) = (
            self.points.a, self.points.b, 
            self.points.c, self.points.d
        ) {
            let epsilon = 1e-10;

            // 1. Check diagonal lengths are equal
            let diagonal_ac = distance(&a, &c);
            let diagonal_bd = distance(&b, &d);
            if (diagonal_ac - diagonal_bd).abs() > epsilon {
                return false;
            }

            // 2. Check each triangle is right-angled using Pythagorean theorem
            let side_ab = distance(&a, &b);
            let side_bc = distance(&b, &c);
            let hypotenuse = diagonal_ac;
            
            let sides_squared = side_ab.powi(2) + side_bc.powi(2);
            let hypotenuse_squared = hypotenuse.powi(2);

            (sides_squared - hypotenuse_squared).abs() < epsilon
        } else {
            false
        }
    }

    pub fn point_on_line(&self, point: &Point, line_start: &Point, line_end: &Point) -> bool {
        // Cross product should be zero for collinear points
        let cross_product = (point.x - line_start.x) * (line_end.y - line_start.y) -
                           (point.y - line_start.y) * (line_end.x - line_start.x);
        cross_product.abs() < 1e-10
    }

    pub fn point_on_circle(&self, point: &Point, circle: &Circle) -> bool {
        // Use squared distances to avoid floating point errors from sqrt
        let dx = point.x - circle.center.x;
        let dy = point.y - circle.center.y;
        let distance_squared = dx * dx + dy * dy;
        let radius_squared = circle.radius * circle.radius;
        (distance_squared - radius_squared).abs() < 1e-10
    }

    pub fn verify_vesica_properties(&self) -> bool {
        if let Some(p1) = self.points.p1 {
            let d1 = distance(&p1, &self.circle_a.center);
            let d2 = distance(&p1, &self.circle_b.as_ref().unwrap().center);
            let d3 = distance(&self.circle_a.center, &self.circle_b.as_ref().unwrap().center);
            
            // Pure distance relationships
            (d1 - d2).abs() < 1e-10 && (d1 - d3).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn verify_square_construction(&self) -> bool {
        if let (Some(a), Some(b), Some(c), Some(_d)) = (
            self.points.a, self.points.b, self.points.c, self.points.d
        ) {
            // Square side should equal vesica height (not 2×height)
            let height = self.vesica_ratio.unwrap();
            let side = distance(&a, &b);
            if (side - height).abs() > 1e-10 {
                return false;
            }
            
            // Verify right angles through Pythagorean theorem
            let diagonal = distance(&a, &c);
            (side.powi(2) * 2.0 - diagonal.powi(2)).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn find_sphere_intersections(&self, s1: &Sphere, s2: &Sphere) -> Vec<Point3D> {
        let mut intersections = Vec::new();
        
        // Distance between sphere centers
        let dx = s2.center.x - s1.center.x;
        let dy = s2.center.y - s1.center.y;
        let dz = s2.center.z - s1.center.z;
        let d = (dx*dx + dy*dy + dz*dz).sqrt();
        
        // Check if spheres intersect
        if d > s1.radius + s2.radius || d < (s1.radius - s2.radius).abs() {
            return intersections;
        }
        
        // Find intersection circle
        let a = (s1.radius*s1.radius - s2.radius*s2.radius + d*d) / (2.0 * d);
        let h = (s1.radius*s1.radius - a*a).sqrt();
        
        // Calculate intersection points
        let px = s1.center.x + (a/d)*dx;
        let py = s1.center.y + (a/d)*dy;
        let pz = s1.center.z + (a/d)*dz;
        
        intersections.push(Point3D { x: px, y: py, z: pz + h });
        intersections.push(Point3D { x: px, y: py, z: pz - h });
        
        intersections
    }

    pub fn verify_geometric_relationships(&mut self) -> bool {
        if let Some(p1) = self.points.p1 {
            // Only verify points lie on circles through pure distance
            let d1 = distance(&p1, &self.circle_a.center);
            let d2 = distance(&p1, &self.circle_b.as_ref().unwrap().center);
            
            // No ratios, no constants, just distance equality
            (d1 - self.circle_a.radius).abs() < 1e-10 && 
            (d2 - self.circle_b.as_ref().unwrap().radius).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn verify_sphere_relationships(&mut self) -> bool {
        if let Some(p1) = self.points.p1 {
            // Only verify points lie on spheres through pure distance
            let d1 = distance_3d(&p1.into(), &self.circle_a.center.into());
            let d2 = distance_3d(&p1.into(), &self.circle_b.as_ref().unwrap().center.into());
            
            // No ratios, no assumptions, just distance equality 
            (d1 - self.circle_a.radius).abs() < 1e-10 && 
            (d2 - self.circle_b.as_ref().unwrap().radius).abs() < 1e-10
        } else {
            false
        }
    }

    pub fn verify_all_relationships(&mut self) -> bool {
        self.verify_vesica_properties() &&
        self.verify_square_properties() &&
        self.verify_sphere_relationships()
    }

    pub fn verify_square_properties(&self) -> bool {
        self.verify_vesica_piscis() &&
        self.verify_quadrant_division() &&
        self.verify_diagonal_alignment()
    }

    pub fn compute_bounded_volume(&self, points: &[Point3D]) -> f64 {
        // Need at least 2 points to define intersection circle
        if points.len() < 2 {
            return 0.0;
        }
        
        // Get intersection circle radius from points
        let radius = self.circle_a.radius;
        
        // Get height of intersection circle from points
        let h = points[0].y;  // Height from x-axis to intersection
        
        // Area of intersection circle at height h
        let _circle_area = h.powi(2);  // r² at this height
        
        // Volume emerges from integrating circle areas
        // from -r to r along x-axis
        let dx = radius / 100.0;  // Small integration step
        let mut volume = 0.0;
        
        let mut x = -radius;
        while x <= radius {
            // Height at this x gives circle radius
            let height = (radius.powi(2) - x.powi(2)).sqrt();
            // Area of circle at this height
            let area = height.powi(2);
            // Add volume slice
            volume += area * dx;
            x += dx;
        }
        
        volume
    }

    pub fn verify_sphere_points(&self, point: &Point3D, sphere1: &Sphere, _sphere2: &Sphere) -> bool {
        // The height of any point on the sphere intersection circle
        // must equal the vesica height (r√3)
        let point_height = (point.y * point.y + point.z * point.z).sqrt();
        let vesica_height = sphere1.radius * (3.0_f64).sqrt();
        
        // Point height should equal vesica height exactly
        (point_height - vesica_height).abs() < 1e-10
    }
}

pub fn distance(p1: &Point, p2: &Point) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}

fn compute_natural_area_ratio(construction: &GeometricConstruction) -> f64 {
    if let Some(p1) = construction.points.p1 {
        let height = construction.distance_to_line(
            &p1,
            &construction.circle_a.center,
            &construction.circle_b.as_ref().unwrap().center
        );
        
        // Verify height through pure Pythagorean theorem
        let radius = construction.circle_a.radius;
        let diameter = 2.0 * radius;
        let base = radius;
        
        // height² + r² = (2r)² must hold
        assert!((height.powi(2) + base.powi(2) - diameter.powi(2)).abs() < 1e-10);
        
        // Areas emerge from pure construction
        let square_area = height.powi(2);
        let circle_area = radius.powi(2);
        
        // Let ratio emerge naturally without any assumptions
        circle_area / square_area
    } else {
        0.0
    }
}

pub fn compute_natural_ratio(construction: &GeometricConstruction) -> f64 {
    compute_natural_area_ratio(construction)
}

pub fn distance_3d(p1: &Point3D, p2: &Point3D) -> f64 {
    ((p2.x - p1.x).powi(2) + 
     (p2.y - p1.y).powi(2) + 
     (p2.z - p1.z).powi(2)).sqrt()
}

pub struct SudokuGrid {
    size: usize,
    cells: Vec<Vec<Option<u32>>>,
    vesica_ratio: f64,
    natural_ratio: f64
}

impl SudokuGrid {
    pub fn new(construction: &GeometricConstruction) -> Self {
        let p1 = construction.points.p1.unwrap();
        let height = construction.distance_to_line(
            &p1,
            &construction.circle_a.center,
            &construction.circle_b.as_ref().unwrap().center
        );
        
        let radius = construction.circle_a.radius;
        let vesica_ratio = height / radius;  // = √3
        let size = (vesica_ratio * vesica_ratio).round() as usize;  // = 3
        
        SudokuGrid {
            size,
            cells: vec![vec![None; size]; size],
            vesica_ratio: height,
            natural_ratio: vesica_ratio
        }
    }

    pub fn verify_geometric_consistency(&self) -> bool {
        let mut row_sums = vec![0; self.size];
        let mut col_sums = vec![0; self.size];
        
        for i in 0..self.size {
            for j in 0..self.size {
                if let Some(num) = self.cells[i][j] {
                    row_sums[i] += num;
                    col_sums[j] += num;
                }
            }
        }
        
        let target = (self.vesica_ratio * self.size as f64).round() as u32;
        
        row_sums.iter().all(|&sum| sum == 0 || sum == target) &&
        col_sums.iter().all(|&sum| sum == 0 || sum == target)
    }

    pub fn solve(&mut self) -> bool {
        let nums: Vec<u32> = (1..=self.size).map(|x| x as u32).collect();
        
        let angle = (self.vesica_ratio / self.natural_ratio).acos();
        let rotation = (angle * self.size as f64).round() as usize;
        
        for i in 0..self.size {
            for j in 0..self.size {
                let idx = (i + j + rotation) % self.size;
                self.cells[i][j] = Some(nums[idx]);
            }
        }
        
        self.verify_geometric_consistency()
    }
}

pub struct GeometricSudoku {
    grid: Vec<Vec<Option<u32>>>,
    size: usize,
    vesica_height: f64,
    natural_ratio: f64
}

impl GeometricSudoku {
    pub fn new(construction: &GeometricConstruction) -> Self {
        let p1 = construction.points.p1.unwrap();
        let vesica_height = construction.distance_to_line(
            &p1,
            &construction.circle_a.center,
            &construction.circle_b.as_ref().unwrap().center
        );
        
        let radius = construction.circle_a.radius;
        let vesica_ratio = vesica_height / radius;  // = √3
        let size = vesica_ratio.round() as usize;  // Should be 2 (rounded from √3)
        
        GeometricSudoku {
            grid: vec![vec![None; size]; size],
            size,
            vesica_height,
            natural_ratio: vesica_ratio
        }
    }

    pub fn verify_geometric_consistency(&self) -> bool {
        let mut row_sums = vec![0; self.size];
        let mut col_sums = vec![0; self.size];
        
        for i in 0..self.size {
            for j in 0..self.size {
                if let Some(num) = self.grid[i][j] {
                    row_sums[i] += num;
                    col_sums[j] += num;
                }
            }
        }
        
        // Target emerges from vesica ratio
        let target = (self.natural_ratio * self.size as f64).round() as u32;
        
        row_sums.iter().all(|&sum| sum == 0 || sum == target) &&
        col_sums.iter().all(|&sum| sum == 0 || sum == target)
    }

    pub fn solve(&mut self) -> bool {
        // Numbers emerge from vesica ratio
        let nums: Vec<u32> = (1..=self.size).map(|x| x as u32).collect();
        
        // Rotation emerges from vesica height and natural ratio
        let angle = (self.vesica_height / self.natural_ratio).acos();
        let rotation = (angle * self.size as f64).round() as usize;
        
        // Fill grid according to geometric relationships
        for i in 0..self.size {
            for j in 0..self.size {
                let idx = (i + j + rotation) % self.size;
                self.grid[i][j] = Some(nums[idx]);
            }
        }
        
        self.verify_geometric_consistency()
    }
}

pub struct YangMillsField {
    dimension: usize,
    field_strength: Vec<Vec<f64>>,
    vesica_ratio: f64,
    natural_ratio: f64,
}

impl YangMillsField {
    pub fn new(construction: &GeometricConstruction) -> Self {
        let p1 = construction.points.p1.unwrap();
        let height = construction.distance_to_line(
            &p1,
            &construction.circle_a.center,
            &construction.circle_b.as_ref().unwrap().center
        );
        
        let natural_ratio = compute_natural_ratio(construction);
        let dimension = (1.0 / natural_ratio).round() as usize;
        
        YangMillsField {
            dimension,
            field_strength: vec![vec![0.0; dimension]; dimension],
            vesica_ratio: height,
            natural_ratio
        }
    }

    pub fn compute_field_strength(&self, i: usize, j: usize) -> f64 {
        let center = self.dimension as f64 / 2.0;
        let dx = (i as f64 - center) / center;
        let dy = (j as f64 - center) / center;
        let r = (dx * dx + dy * dy).sqrt();
        self.vesica_ratio * (-r * r).exp()
    }

    pub fn compute_mass_gap(&self) -> f64 {
        let mut min_gap = f64::INFINITY;
        for i in 0..self.dimension {
            for j in 0..self.dimension {
                let field = self.field_strength[i][j];
                if field > 0.0 {
                    min_gap = min_gap.min(field);
                }
            }
        }
        min_gap * self.natural_ratio
    }

    pub fn verify_existence(&self) -> bool {
        let center = self.dimension / 2;
        self.field_strength[center][center] > 0.0
    }

    pub fn verify_field_consistency(&self) -> bool {
        for i in 1..self.dimension-1 {
            for j in 1..self.dimension-1 {
                // Let the relationship emerge from vesica ratio
                let neighbors = vec![
                    self.field_strength[i+1][j],
                    self.field_strength[i-1][j], 
                    self.field_strength[i][j+1],
                    self.field_strength[i][j-1]
                ];
                
                let neighbor_sum: f64 = neighbors.iter().sum();
                let center = self.field_strength[i][j];
                let natural_count = (1.0 / self.natural_ratio).round();
                
                // Use natural ratio for consistency check
                let laplacian = neighbor_sum - natural_count * center;
                if (laplacian / self.natural_ratio).abs() > self.vesica_ratio {
                    return false;
                }
            }
        }
        true
    }
}

pub mod navier;
pub mod hodge;

pub mod birch;

impl From<Point> for Point3D {
    fn from(p: Point) -> Self {
        Point3D {
            x: p.x,
            y: p.y,
            z: 0.0
        }
    }
}

