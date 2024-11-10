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
    pub intersections: Vec<Point>
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
    pub d: Option<Point>,
    pub e: Option<Point>,
    pub extended_lines: Vec<(Point, Point)>
}

pub mod visualization;

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
            intersections: Vec::new()
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
            // Construct lines for square vertices
            // Line P1-P3 intersects with C4-C2 to form vertex A
            self.points.a = self.find_line_intersection(&p1, &p3, &c4, &c2);
            
            // Line P1-P3 intersects with C1-C3 to form vertex B
            self.points.b = self.find_line_intersection(&p1, &p3, &c1, &c3);
            
            // Line P5-P2 intersects with C4-C2 to form vertex D
            self.points.d = self.find_line_intersection(&p5, &p2, &c4, &c2);
            
            // Line P5-P2 intersects with C1-C3 to form vertex E
            self.points.e = self.find_line_intersection(&p5, &p2, &c1, &c3);

            // Store the extended lines for visualization
            self.points.extended_lines = vec![
                (p1, p3),
                (p5, p2),
                (c1, c3),
                (c4, c2)
            ];
        }
    }

    pub fn verify_square(&self) -> bool {
        if let (Some(a), Some(b), Some(d), Some(e)) = (
            self.points.a, self.points.b, self.points.d, self.points.e
        ) {
            let epsilon = 1e-10;
            
            // Check all sides are equal
            let side1 = distance(&a, &b);
            let side2 = distance(&b, &e);
            let side3 = distance(&e, &d);
            let side4 = distance(&d, &a);
            
            let sides_equal = (side1 - side2).abs() < epsilon
                && (side2 - side3).abs() < epsilon
                && (side3 - side4).abs() < epsilon;
            
            // Check angles are 90 degrees using dot product
            let check_right_angle = |p1: &Point, p2: &Point, p3: &Point| -> bool {
                let v1 = Point { x: p2.x - p1.x, y: p2.y - p1.y };
                let v2 = Point { x: p3.x - p2.x, y: p3.y - p2.y };
                let dot_product = v1.x * v2.x + v1.y * v2.y;
                dot_product.abs() < epsilon
            };
            
            let angles_right = check_right_angle(&a, &b, &e)
                && check_right_angle(&b, &e, &d)
                && check_right_angle(&e, &d, &a)
                && check_right_angle(&d, &a, &b);
            
            sides_equal && angles_right
        } else {
            false
        }
    }
}

pub fn distance(p1: &Point, p2: &Point) -> f64 {
    ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
}

#[cfg(test)]
mod tests;

