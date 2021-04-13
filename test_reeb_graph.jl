using FileIO
using MeshIO
using WGLMakie
using Rotations

include("mesh_utils.jl")
include("crit_pts.jl")
include("mesh_io.jl")


X, T = MeshUtils.readoff("moomoo.off")
X = X * RotX(pi/2)

model = load("moomoo.off")
vertices = decompose(Point3f0, model)
faces = decompose(TriangleFace{Int}, model)
# vertices = X
# faces = T
coordinates = [vertices[i][j] for i = 1:length(vertices), j = 1:3]
coordinates = coordinates * RotX(pi/2)
connectivity = [faces[i][j] for i = 1:length(faces), j = 1:3]

function plot_moomoo()
  mesh(
      coordinates, connectivity,
      shading=true,
      transparency=true,
      figure=(resolution=(700, 1000),),
      color = (:red, 0.1)
  )
end

function plot_vertex_normals()
  plot_moomoo()
  VN = vertex_normals(X, T)
  arrows!(
    Point3.(eachrow(X)),
    Point3.(eachrow(VN)),
    arrowsize = 0.1,
    linewidth = 20,
    linecolor = :green,
    arrowcolor = :darkblue)
end

function plot_spheres(idxs, color)
  for X_idx in idxs
    mesh!(Sphere(Point3f0(X[X_idx,:]), 0.5f0), transparency=false, color=color)
  end
end


function plot_crit_points()
  # plot_moomoo()
  bounds = Node([
    Point3f0(0., 0., 0.),
    Point3f0(15., 0., 0.),
    Point3f0(-15., 0., 0.),
    Point3f0(0., 0., 20.),
    Point3f0(0., 0., -30.),
    Point3f0(0., 20., 0.),
    Point3f0(0., -15., 0.),
  ])
  lines(bounds, linewidth=0, transparency=true, color=(:black, 0.0),
      shading=true,
      figure=(resolution=(700, 1000),),
  )
  crits = crit_points(X, T)
  plot_spheres(crit_idxs(crit_mins(crits)), :yellow)
  plot_spheres(crit_idxs(crit_maxs(crits)), :green)
  plot_spheres(crit_idxs(crit_saddles(crits)), :blue)
end


function plot_crit_sets()
  crits = crit_points(X, T)
  crit_idx = 2
  my_crit_set = crit_set(crit_saddles(crits)[crit_idx], X, T)
  MeshUtils.showdescriptor(X, T, map(x -> (x[1] in my_crit_set) ? 1.0 : 0.0, 1:size(X, 1)))
  plot_spheres([crit_saddles(crits)[crit_idx][2][1]], :red)
end


plot_crit_points()

rg = reeb_graph(X, T)
crits = crit_points(X, T)
for e in edges(rg)
  println(e)
  src_crits_idx = e.src
  dst_crits_idx = e.dst
  src_vert_idx = crits[src_crits_idx][2][1]
  dst_vert_idx = crits[dst_crits_idx][2][1]
  line_seg = Node([Point3f0(X[src_vert_idx, :]), Point3f0(X[dst_vert_idx, :])])
  lines!(line_seg, linewidth=400)
end
println(rg)

# linesegments!(scene, [Point3f0(0., 0., 0.), Point3f0(0., 0., 1.)], color = :red, fxaa = true)
# line_seg = Node([Point3f0(10., 0., 0.), Point3f0(10., 0., 10.)])
# lines!(line_seg, linewidth=400)


# using Plots
# using GraphRecipes
# am = Matrix(adjacency_matrix(rg))
# n = size(am, 1)
# graphplot(am,
#           markersize = 0.2,
#           node_weights = 1:n,
#           markercolor = range(colorant"yellow", stop=colorant"red", length=n),
#           names = 1:n,
#           fontsize = 10,
#           linecolor = :darkgrey
#           )
