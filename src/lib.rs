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
pub struct Line {
    start: Point,
    end: Point
}

#[derive(Debug)]
pub struct GeometricConstruction {
    circle_a: Circle,
    circle_b: Option<Circle>,
    points: ConstructionPoints,
    intersections: Vec<Point>
}

#[derive(Debug, Default)]
struct ConstructionPoints {
    p1: Option<Point>,
    p2: Option<Point>,
    p3: Option<Point>,
    p4: Option<Point>,
    p5: Option<Point>,
    p6: Option<Point>,
    c1: Option<Point>,
    c2: Option<Point>,
    c3: Option<Point>,
    c4: Option<Point>
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
            intersections: Vec::new()
        }
    }

    pub fn construct(&mut self) {
        // Initial circle intersection defines all other points
        let initial_circle = Circle {
            center: Point { x: self.circle_a.radius, y: 0.0 },
            radius: self.circle_a.radius
        };
        self.circle_b = Some(initial_circle);

        if let Some(circle_b) = &self.circle_b {
            let intersections = self.find_circle_intersections(&self.circle_a, circle_b);
            if intersections.len() == 2 {
                self.points.p1 = Some(intersections[0]);
                self.points.p2 = Some(intersections[1]);
            }
        }

        self.construct_arcs();
        self.construct_center_lines();
        self.construct_square();
    }



    fn find_circle_intersections(&self, circle1: &Circle, circle2: &Circle) -> Vec<Point> {
        let dx = circle2.center.x - circle1.center.x;
        let dy = circle2.center.y - circle1.center.y;
        let d = (dx * dx + dy * dy).sqrt();
        
        if d > circle1.radius + circle2.radius || d < (circle1.radius - circle2.radius).abs() {
            return Vec::new();
        }
        
        let a = (circle1.radius * circle1.radius - circle2.radius * circle2.radius + d * d) / (2.0 * d);
        let h = (circle1.radius * circle1.radius - a * a).sqrt();
        
        let x2 = circle1.center.x + (dx * a) / d;
        let y2 = circle1.center.y + (dy * a) / d;
        
        vec![
            Point {
                x: x2 + (h * dy) / d,
                y: y2 - (h * dx) / d
            },
            Point {
                x: x2 - (h * dy) / d,
                y: y2 + (h * dx) / d
            }
        ]
    }

    fn find_line_circle_intersection(&self, line: &Line, circle: &Circle) -> Vec<Point> {
        let dx = line.end.x - line.start.x;
        let dy = line.end.y - line.start.y;
        let a = dx * dx + dy * dy;
        let b = 2.0 * (dx * (line.start.x - circle.center.x) + 
                      dy * (line.start.y - circle.center.y));
        let c = circle.center.x * circle.center.x + 
                circle.center.y * circle.center.y +
                line.start.x * line.start.x +
                line.start.y * line.start.y -
                2.0 * (circle.center.x * line.start.x + 
                       circle.center.y * line.start.y) -
                circle.radius * circle.radius;
        
        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return Vec::new();
        }
        
        let t1 = (-b + discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b - discriminant.sqrt()) / (2.0 * a);
        
        vec![
            Point {
                x: line.start.x + t1 * dx,
                y: line.start.y + t1 * dy
            },
            Point {
                x: line.start.x + t2 * dx,
                y: line.start.y + t2 * dy
            }
        ]
    }

    fn construct_arcs(&mut self) {
        if let (Some(p1), Some(p2)) = (self.points.p1, self.points.p2) {
            // Create circles centered at P1 and P2 with radius equal to original circle
            let circle_p1 = Circle { center: p1, radius: self.circle_a.radius };
            let circle_p2 = Circle { center: p2, radius: self.circle_a.radius };
            
            // Find intersections with original circle
            let p1_intersections = self.find_circle_intersections(&self.circle_a, &circle_p1);
            let p2_intersections = self.find_circle_intersections(&self.circle_a, &circle_p2);
            
            if p1_intersections.len() == 2 {
                self.points.p3 = Some(p1_intersections[0]);
                self.points.p4 = Some(p1_intersections[1]);
            }
            
            if p2_intersections.len() == 2 {
                self.points.p5 = Some(p2_intersections[0]);
                self.points.p6 = Some(p2_intersections[1]);
            }
        }
    }

    fn construct_center_lines(&mut self) {
        if let (Some(p4), Some(p6)) = (self.points.p4, self.points.p6) {
            // Line through P4 and center
            let line1 = Line { start: p4, end: self.circle_a.center };
            let intersections1 = self.find_line_circle_intersection(&line1, &self.circle_a);
            if intersections1.len() == 2 {
                self.points.c1 = Some(intersections1[0]);
                self.points.c2 = Some(intersections1[1]);
            }

            // Line through P6 and center
            let line2 = Line { start: p6, end: self.circle_a.center };
            let intersections2 = self.find_line_circle_intersection(&line2, &self.circle_a);
            if intersections2.len() == 2 {
                self.points.c3 = Some(intersections2[0]);
                self.points.c4 = Some(intersections2[1]);
            }
        }
    }

    fn construct_square(&mut self) {
        if let (Some(p1), Some(p2), Some(p3), Some(_p4), Some(p5), Some(_p6)) = (
            self.points.p1, self.points.p2, 
            self.points.p3, self.points.p4,
            self.points.p5, self.points.p6
        ) {
            // Calculate intersection points between lines
            let _line_p1_p3 = Line { start: p1, end: p3 };
            let _line_p2_p5 = Line { start: p2, end: p5 };
            
            // Store intersections for verification
            if let (Some(c1), Some(c2), Some(c3), Some(c4)) = (
                self.points.c1, self.points.c2,
                self.points.c3, self.points.c4
            ) {
                self.intersections = vec![c1, c2, c3, c4];
            }
        }
    }

    pub fn verify_square(&self) -> bool {
        if let (Some(p1), Some(p2), Some(p3), Some(p4)) = (
            &self.points.c1, &self.points.c2, 
            &self.points.c3, &self.points.c4) {
            
            let d1 = distance(p1, p2);
            let d2 = distance(p2, p3);
            let d3 = distance(p3, p4);
            let d4 = distance(p4, p1);
            
            let diag1 = distance(p1, p3);
            let diag2 = distance(p2, p4);
            
            let epsilon = 1e-10;
            
            (d1 - d2).abs() < epsilon &&
            (d2 - d3).abs() < epsilon &&
            (d3 - d4).abs() < epsilon &&
            (diag1 - diag2).abs() < epsilon
        } else {
            false
        }
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
}

pub fn distance(p1: &Point, p2: &Point) -> f64 {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    (dx * dx + dy * dy).sqrt()
}

#[cfg(test)]
mod tests;