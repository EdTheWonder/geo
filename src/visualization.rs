use plotters::prelude::*;
use plotters::element::Circle as PlottersCircle;
use plotters::coord::types::RangedCoordf32;
use crate::{GeometricConstruction, Point};

fn draw_label<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf32, RangedCoordf32>>,
    point: &Point,
    label: &str
) -> Result<(), Box<dyn std::error::Error>> 
where
    DB::ErrorType: 'static,
{
    chart.draw_series(std::iter::once(Text::new(
        label.to_string(),
        (point.x as f32, point.y as f32),
        ("sans-serif", 15.0).into_font(),
    )))?;
    Ok(())
}

fn draw_circle<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf32, RangedCoordf32>>,
    center: &Point,
    radius: f64,
    color: &RGBColor
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let points: Vec<(f32, f32)> = (0..=100).map(|i| {
        let angle = 2.0 * std::f64::consts::PI * (i as f64 / 100.0);
        (
            (center.x + radius * angle.cos()) as f32,
            (center.y + radius * angle.sin()) as f32
        )
    }).collect();
    
    chart.draw_series(LineSeries::new(points, color))?;
    Ok(())
}

fn draw_arc<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf32, RangedCoordf32>>,
    center: &Point,
    start_point: &Point,
    end_point: &Point,
    radius: f64,
    color: &RGBColor
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let start_angle = (start_point.y - center.y).atan2(start_point.x - center.x);
    let end_angle = (end_point.y - center.y).atan2(end_point.x - center.x);
    
    // Ensure we draw the shorter arc
    let mut angle_diff = end_angle - start_angle;
    if angle_diff > std::f64::consts::PI {
        angle_diff -= 2.0 * std::f64::consts::PI;
    } else if angle_diff < -std::f64::consts::PI {
        angle_diff += 2.0 * std::f64::consts::PI;
    }
    
    let points: Vec<(f32, f32)> = (0..=50).map(|i| {
        let t = i as f64 / 50.0;
        let angle = start_angle + t * angle_diff;
        (
            (center.x + radius * angle.cos()) as f32,
            (center.y + radius * angle.sin()) as f32
        )
    }).collect();
    
    chart.draw_series(LineSeries::new(points, color))?;
    Ok(())
}

pub fn draw_construction(construction: &GeometricConstruction) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("geometric_construction.png", (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let min = -2.0f32;
    let max = 2.0f32;
    
    let mut chart = ChartBuilder::on(&root)
        .margin(5)
        .build_cartesian_2d(min..max, min..max)?;

    // Step 1: Draw Circle A (Unity Circle)
    draw_circle(&mut chart, &construction.circle_a.center, construction.circle_a.radius, &BLUE)?;
    draw_label(&mut chart, &construction.circle_a.center, "O")?;

    // Step 2: Draw Circle B (Iteration Circle)
    if let Some(circle_b) = &construction.circle_b {
        draw_circle(&mut chart, &circle_b.center, circle_b.radius, &RGBColor(100, 100, 255))?;
        draw_label(&mut chart, &circle_b.center, "P")?;
    }

    // Step 3: Draw Arcs of Convergence
    if let (Some(p1), Some(p3), Some(p4)) = (construction.points.p1, construction.points.p3, construction.points.p4) {
        draw_arc(&mut chart, &p1, &p3, &p4, construction.circle_a.radius, &RGBColor(255, 100, 100))?;
    }

    if let (Some(p2), Some(p5), Some(p6)) = (construction.points.p2, construction.points.p5, construction.points.p6) {
        draw_arc(&mut chart, &p2, &p5, &p6, construction.circle_a.radius, &RGBColor(255, 100, 100))?;
    }

    // Step 4: Draw Lines through Center
    if let (Some(p4), Some(c1), Some(c2)) = (
        construction.points.p4,
        construction.points.c1,
        construction.points.c2
    ) {
        // Draw line from P4 to center to verify construction
        chart.draw_series(LineSeries::new(
            vec![(p4.x as f32, p4.y as f32), (0.0, 0.0)],
            &RGBColor(100, 200, 100)
        ))?;
        // Draw C1-C2 line
        chart.draw_series(LineSeries::new(
            vec![(c1.x as f32, c1.y as f32), (c2.x as f32, c2.y as f32)],
            &RGBColor(100, 200, 100)
        ))?;
    }

    if let (Some(p6), Some(c3), Some(c4)) = (
        construction.points.p6,
        construction.points.c3,
        construction.points.c4
    ) {
        // Draw line from P6 to center to verify construction
        chart.draw_series(LineSeries::new(
            vec![(p6.x as f32, p6.y as f32), (0.0, 0.0)],
            &RGBColor(100, 200, 100)
        ))?;
        // Draw C3-C4 line
        chart.draw_series(LineSeries::new(
            vec![(c3.x as f32, c3.y as f32), (c4.x as f32, c4.y as f32)],
            &RGBColor(100, 200, 100)
        ))?;
    }

    // Draw extended construction lines
    for (p1, p2) in &construction.points.extended_lines {
        chart.draw_series(LineSeries::new(
            vec![(p1.x as f32, p1.y as f32), (p2.x as f32, p2.y as f32)],
            &RGBColor(150, 150, 150)
        ))?;
    }

    // Draw the square with thicker lines
    if let (Some(a), Some(b), Some(d), Some(c)) = (
        construction.points.a,
        construction.points.b,
        construction.points.d,
        construction.points.c
    ) {
        let square_points = vec![
            (a.x as f32, a.y as f32),
            (b.x as f32, b.y as f32),
            (c.x as f32, c.y as f32),
            (d.x as f32, d.y as f32),
            (a.x as f32, a.y as f32),
        ];
        
        chart.draw_series(LineSeries::new(
            square_points,
            &RGBColor(0, 255, 0)
        ))?;
    }

    // Draw all construction points with labels
    let construction_points = [
        (construction.points.p1, "P₁", &RED),
        (construction.points.p2, "P₂", &RED),
        (construction.points.p3, "P₃", &RED),
        (construction.points.p4, "P₄", &RED),
        (construction.points.p5, "P₅", &RED),
        (construction.points.p6, "P₆", &RED),
        (construction.points.c1, "C₁", &GREEN),
        (construction.points.c2, "C₂", &GREEN),
        (construction.points.c3, "C₃", &GREEN),
        (construction.points.c4, "C₄", &GREEN),
        (construction.points.a, "A", &RGBColor(0, 255, 0)),
        (construction.points.b, "B", &RGBColor(0, 255, 0)),
        (construction.points.d, "D", &RGBColor(0, 255, 0)),
        (construction.points.c, "C", &RGBColor(0, 255, 0)),
    ];

    for (point_opt, label, color) in construction_points.iter() {
        if let Some(point) = point_opt {
            chart.draw_series(PointSeries::of_element(
                vec![(point.x as f32, point.y as f32)],
                3,
                *color,
                &|coord, size, style| PlottersCircle::new(coord, size, style.filled()),
            ))?;
            draw_label(&mut chart, point, label)?;
        }
    }

    root.present()?;
    Ok(())
}
