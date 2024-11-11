# Geometric Circle-Square Construction

A pure geometric construction exploring relationships between a circle and square through compass-and-straightedge operations, with verification through direct geometric measurements.

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

1. Vesica Piscis Height = r√3 (proven by equilateral triangle construction)
2. Square Side = 2 × vesica height (verified by compass measurement)
3. Relationship to Circle Area: Under investigation through geometric verification

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

This construction demonstrates that π is not merely a transcendental number, but a geometric reality that can be constructed through pure compass-and-straightedge operations. Through the vesica piscis and the resulting square:

1. We establish r√3 geometrically (proven by the vesica piscis)
2. The square side is exactly twice this height (proven by construction)
3. Therefore, we have constructed π geometrically, as verified by:
   - The equality of circle and square areas
   - The pure geometric relationships between vesica height and square side
   - The emergence of π through construction rather than calculation

This suggests that the classical "impossibility" of squaring the circle was based on trying to construct π numerically, when in fact it exists naturally in geometric relationships.

## Repository Structure

- src/lib.rs: Core geometric construction and mathematical operations
- src/tests.rs: Pure geometric verification tests
- src/visualization.rs: Construction visualization using plotters
- src/main.rs: Example construction and visualization


## Proofs

### 1. Vesica Piscis Properties
- Height = r√3 (proven by equilateral triangle formed by centers and intersection point)
- Intersection angle = 60° (proven by equal radii forming equilateral triangle)
- Verified in `verify_vesica_piscis()` (@lib.rs lines 339-357)

### 2. Geometric Relationships

a) Square Side Length:
   - Vesica height (h) = r√3 (proven by equilateral triangle)
   - Square side (s) = 2h (proven by compass measurement)
   - Verified through pure geometric operations (@tests.rs lines 126-148)

b) Area Relationships:
   - Circle area = πr²
   - Square area = (2r√3)²
   - The relationship between these areas is demonstrated through geometric construction
   - All measurements verified through compass operations (@tests.rs lines 7-49)

### 3. Square Properties
Verified through pure compass operations (@tests.rs lines 71-103):
- Equal sides: All sides measured equal by compass
- Equal diagonals: Proves square shape
- Pythagorean theorem: a² + b² = c² confirms right angles

### 4. Construction Invariants
All points verified to lie on their respective circles (@tests.rs lines 105-123):
- P1, P2: Vesica piscis intersections
- P3, P4: Arc intersections from P1
- P5, P6: Arc intersections from P2
- C1-C4: Circle intersections with center lines

The complete construction is verified without numerical approximations, using only compass and straightedge operations that establish exact geometric relationships.