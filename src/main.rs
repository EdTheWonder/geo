use geometric_construction::GeometricConstruction;
mod visualization;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut construction = GeometricConstruction::new(1.0);
    construction.construct();
    visualization::draw_construction(&construction)?;
    println!("Construction completed and saved to geometric_construction.png");
    Ok(())
}