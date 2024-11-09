use plotters::prelude::*;
use crate::GeometricConstruction;

pub fn draw_construction(construction: &GeometricConstruction) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("geometric_construction.png", (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let min = -2.0f32;
    let max = 2.0f32;
    
    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(min..max, min..max)?;

    // Draw circles
    let points: Vec<(f32, f32)> = (0..=100).map(|i| {
        let angle = 2.0 * std::f32::consts::PI * (i as f32 / 100.0);
        (angle.cos(), angle.sin())
    }).collect();

    // Draw unit circle
    chart.draw_series(LineSeries::new(
        points.iter().map(|(x, y)| (*x, *y)),
        &BLUE,
    ))?;

    // Draw points
    if let Some(points) = construction.get_points() {
        for point in points {
            chart.draw_series(PointSeries::of_element(
                vec![(point.x as f32, point.y as f32)],
                3,
                &RED,
                &|c, s, st| Circle::new(c, s, st.filled()),
            ))?;
        }
    }

    root.present()?;
    Ok(())
}
