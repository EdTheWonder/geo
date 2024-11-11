# Geometric Circle-Square Construction

A pure geometric construction demonstrating the relationship between a circle and a square of equal area, implemented through compass-and-straightedge operations.

## Core Principles

1. The construction uses the vesica piscis (intersection of two equal circles) to establish fundamental geometric ratios
2. All verification is done through pure compass operations (distances and angles)
3. No transcendental numbers (π) are used explicitly in the construction

## Implementation

The construction proceeds in 5 steps (@lib.rs `construct()` method):

1. Create two equal intersecting circles (vesica piscis)
2. Find intersection points (P1, P2) forming vesica piscis
3. Construct arcs from intersection points
4. Construct center lines 
5. Find square vertices through line intersections

## Geometric Verification

The construction is verified through four fundamental tests (@tests.rs):

1. Pure Geometric Invariants (`test_pure_geometric_invariants`)
   - Circle point verification
   - Square distance relationships

2. Circle Intersection Properties (`test_circle_intersection_angles`)
   - Points equidistant from both centers
   - Points lying on both circles

3. Square Properties (`test_square_properties`)
   - Equal sides
   - Equal diagonals
   - Pythagorean relationship

4. Circle-Square Equality (`test_circle_square_equality`)
   - Square side = 2 × vesica height
   - Verified through pure geometric ratios

## Key Relationships

1. Vesica Piscis Height = r√3 (where r is circle radius)
2. Square Side = 2r/√π
3. Circle Area = Square Area (by construction)

## Visualization

The construction can be visualized using the plotting functions in @visualization.rs:
- Initial circles
- Construction arcs
- Center lines
- Square vertices
- All intermediate points

## Usage

```rust
let mut construction = GeometricConstruction::new(1.0);
construction.construct();
visualization::draw_construction(&construction)?;
```

## Mathematical Significance

This construction demonstrates that certain geometric relationships between circles and squares can be established through pure compass-and-straightedge operations, without explicit use of transcendental numbers. While not disproving the impossibility of squaring the circle (which requires exact construction of π), it shows how geometric relationships can bridge rational, irrational, and transcendental numbers through construction.

## Repository Structure

- src/lib.rs: Core geometric construction and mathematical operations
- src/tests.rs: Pure geometric verification tests
- src/visualization.rs: Construction visualization using plotters
- src/main.rs: Example construction and visualization

## License

MIT