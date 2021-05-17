using FileIO
using MeshIO
using WGLMakie
using Rotations

include("mesh_utils.jl")
include("reeb_graph.jl")
include("mesh_io.jl")

# X, T = MeshUtils.readoff("models/sphere.off")
# rot_mat = I

# X, T = MeshUtils.readoff("models/moomoo.off")
# rot_mat = RotZ(-0.5) * RotX(pi/2)

# X, T = MeshUtils.readoff("models/torusish.off")
# rot_mat = RotX(-0.3)


# model = load("models/cube.off")
# model = load("models/two_cubes.off")
model = load("models/fork.off")
X = decompose(Point3f0, model)
T = decompose(TriangleFace{Int}, model)

# rot_mat = RotX(-0.3) * RotY(-0.3)
# X = map(x -> Point3f0(rot_mat * x), X)


X = [X[i][j] for i = 1:length(X), j = 1:3]
T = [T[i][j] for i = 1:length(T), j = 1:3]
# println(size(X))
# println(size(T))


# vertices = X
# faces = T
# coordinates = [vertices[i][j] for i = 1:length(vertices), j = 1:3]
# coordinates = coordinates * rot_mat
# connectivity = [faces[i][j] for i = 1:length(faces), j = 1:3]


# X = X * rot_mat


function plot_current_mesh()
  fig, ax, msh = mesh(
      X, T,
      shading=true,
      transparency=true,
      figure=(resolution=(700, 1000),),
      color = (:red, 0.1),
  )
  # wireframe!(msh[1], color=(:black, 0.6), linewidth=50)
  # display(fig)

  display(wireframe(msh[1], color=(:black, 0.6), linewidth=50))

  # display(meshscatter(
  #     X,
  #     shading=true,
  #     transparency=true,
  #     figure=(resolution=(700, 1000),),
  #     color = (:black, 0.5),
  #     markersize = 0.1
  # ))
end

function plot_vertex_normals()
  plot_current_mesh()
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
    # mesh!(Sphere(Point3f0(X[X_idx,:]), 0.5f0), transparency=false, color=color)
    plot_sphere_by_pos(Point3f0(X[X_idx,:]), color)
  end
end

function plot_sphere_by_pos(pos, color)
  mesh!(Sphere(pos, 0.5f0), transparency=false, color=color)
end

function plot_crit_areas()
  plot_current_mesh()
  # bounds = Node([
  #   Point3f0(0., 0., 0.),
  #   Point3f0(15., 0., 0.),
  #   Point3f0(-15., 0., 0.),
  #   Point3f0(0., 0., 20.),
  #   Point3f0(0., 0., -30.),
  #   Point3f0(0., 20., 0.),
  #   Point3f0(0., -15., 0.),
  # ])
  # display(lines(bounds, linewidth=0, transparency=true, color=(:black, 0.0),
  #     shading=true,
  #     figure=(resolution=(700, 1000),),
  # ))
  crits = crit_areas(X, T)
  for crit_area in crits
    if crit_area.type == crit_min
      color = :yellow
    elseif crit_area.type == crit_max
      color = :green
    elseif crit_area.type == crit_saddle
      color = :blue
    end
    plot_sphere_by_pos(avg_of(crit_area.idxs, X), color)
  end
end


function plot_crit_sets()
  crits = crit_points(X, T)
  crit_idx = 2
  my_crit_set = crit_set(crit_saddles(crits)[crit_idx], X, T)
  MeshUtils.showdescriptor(X, T, map(x -> (x[1] in my_crit_set) ? 1.0 : 0.0, 1:size(X, 1)))
  plot_spheres([crit_saddles(crits)[crit_idx][2][1]], :red)
end


function avg_of(idxs, X)
  res = Point3f0(0.)
  n = size(idxs, 1)
  for i = 1:n
    res += (1. / n) * X[idxs[i],:]
  end
  return res
end


plot_crit_areas()
# plot_current_mesh()

rg = reeb_graph(X, T)
crits = crit_areas(X, T)
println(crits)
for e in edges(rg)
  println(e)
  src_crits_idx = e.src
  dst_crits_idx = e.dst
  src_vert_idxs = crits[src_crits_idx].idxs
  dst_vert_idxs = crits[dst_crits_idx].idxs
  src_pos = avg_of(src_vert_idxs, X)
  dst_pos = avg_of(dst_vert_idxs, X)
  line_seg = Node([src_pos, dst_pos])
  lines!(line_seg, linewidth=800)
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
