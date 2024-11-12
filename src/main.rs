use geometric_construction::{
    GeometricConstruction, 
    RiemannConstruction, 
    visualization, 
    riemann_visualization
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut construction = GeometricConstruction::new(1.0);
    construction.construct();
    visualization::draw_construction(&construction)?;
    println!("Circle-square construction saved to geometric_construction.png");

    let mut riemann = RiemannConstruction::new();
    riemann.construct();
    riemann_visualization::draw_riemann_construction(&riemann)?;
    println!("Riemann construction saved to riemann_construction.png");
    
    Ok(())
}