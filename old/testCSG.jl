using ConstructiveSolidGeometry

top =   Plane(Coord(0.0, 0.0, 150.0),  unitize(Coord(0.0, 0.0, 1.0)),  "reflective")
bot =   Plane(Coord(0.0, 0.0, -150.0), unitize(Coord(0.0, 0.0, -1.0)), "reflective")
left =  Plane(Coord(-.63, 0.0, 0.0),   unitize(Coord(-1.0, 0.0, 0.0)), "reflective")
right = Plane(Coord(0.63, 0.0, 0.0),   unitize(Coord(1.0, 0.0, 0.0)),  "reflective")
up =    Plane(Coord(0.0, 0.63, 0.0),   unitize(Coord(0.0, 1.0, 0.0)),  "reflective")
down =  Plane(Coord(0.0, -0.63, 0.0),  unitize(Coord(0.0, -1.0, 0.0)), "reflective");


clad_outer = InfCylinder(Coord(0.0, 0.0, 0.0), unitize(Coord(0.0, 0.0, 1.0)), 0.4750)
clad_inner = InfCylinder(Coord(0.0, 0.0, 0.0), unitize(Coord(0.0, 0.0, 1.0)), 0.4180)
fuel =       InfCylinder(Coord(0.0, 0.0, 0.0), unitize(Coord(0.0, 0.0, 1.0)), 0.4096);

cells = Cell[]
regions = Region[]


push!(regions, Region(top, -1))
push!(regions, Region(bot, -1))
push!(regions, Region(left, -1))
push!(regions, Region(right, -1))
push!(regions, Region(up, -1))
push!(regions, Region(down, -1))
push!(regions, Region(clad_outer, 1));

ex = :(1 ^ 2 ^ 3 ^ 4 ^ 5 ^ 6 ^ 7)

push!(cells, Cell(regions, ex));


# Make the Cladding Cell
regions = Region[]

push!(regions, Region(top, -1))
push!(regions, Region(bot, -1))
push!(regions, Region(clad_outer, -1))
push!(regions, Region(clad_inner, 1))

ex = :(1 ^ 2 ^ 3 ^ 4)

push!(cells, Cell(regions, ex))

# Make the Gap Cell
regions = Region[]

push!(regions, Region(top, -1))
push!(regions, Region(bot, -1))
push!(regions, Region(clad_inner, -1))
push!(regions, Region(fuel, 1))

ex = :(1 ^ 2 ^ 3 ^ 4)

push!(cells, Cell(regions, ex))

# Make the Fuel Cell
regions = Region[]

push!(regions, Region(top, -1))
push!(regions, Region(bot, -1))
push!(regions, Region(fuel, -1))

ex = :(1 ^ 2 ^ 3)

push!(cells, Cell(regions, ex));

bounding_box = Box(Coord(-.63, -.63, -150), Coord(.63, .63, 150))

geometry = Geometry(cells, bounding_box);

plot_geometry_2D(geometry, Box(Coord(-0.63, -0.63, 0), Coord(.63, 0.63, 0)), 1000)

plot_cell_2D(geometry, Box(Coord(-0.63, -0.63, 0), Coord(.63, 0.63, 0)), 1000, 2)

ray = generate_random_ray(geometry.bounding_box);
