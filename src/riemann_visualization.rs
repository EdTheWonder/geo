use plotters::prelude::*;
use plotters::element::Circle as PlottersCircle;
use plotters::coord::types::RangedCoordf32;
use super::{Circle, RiemannConstruction};

pub fn draw_riemann_construction(construction: &RiemannConstruction) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("riemann_construction.png", (800, 800))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .build_cartesian_2d(-2.0f32..2.0f32, -2.0f32..2.0f32)?;

    // Draw critical line
    chart.draw_series(LineSeries::new(
        vec![(0.5f32, -2.0f32), (0.5f32, 2.0f32)],
        &BLACK.mix(0.5),
    ))?;

    // Draw base circle
    if let Some(circle) = construction.zeta_circles.first() {
        draw_circle(&mut chart, circle, &BLUE)?;
    }

    // Draw critical points (zeta zeros)
    for point in &construction.critical_points {
        chart.draw_series(PointSeries::of_element(
            vec![(point.x as f32, point.y as f32)],
            5,
            &RED,
            &|c, s, st| PlottersCircle::new(c, s, st.filled()),
        ))?;
    }

    // Add labels
    chart.draw_series(vec![Text::new(
        "Critical Line Re(s) = 1/2",
        (0.6f32, 1.8f32),
        ("sans-serif", 20.0).into_font(),
    )])?;

    Ok(())
}

fn draw_circle<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf32, RangedCoordf32>>,
    circle: &Circle,
    color: &RGBColor,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static
{
    let points: Vec<(f32, f32)> = (0..=360).map(|i| {
        let angle = i as f64 * std::f64::consts::PI / 180.0;
        (
            (circle.center.x + circle.radius * angle.cos()) as f32,
            (circle.center.y + circle.radius * angle.sin()) as f32
        )
    }).collect();

    chart.draw_series(LineSeries::new(points, color))?;
    Ok(())
}
