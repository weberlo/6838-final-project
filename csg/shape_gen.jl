include("./csg.jl")

using GeometryBasics


function gen_square_at(s :: Number, center :: Point3)
  return CSG.CSquare(
    center .- (s / 2.),
    center .+ (s / 2.))
end

function treeify(ps)
  res = ps[1]
  for i=2:length(ps)
    res = CSG.CUnion(res, ps[i])
  end
  return res
end

function gen_torus()
  # 2D version of torus-ish shape from GoodNotes
  size = 4.
  c1 = gen_square_at(size, Point3(0., 2., 0.))
  c2 = gen_square_at(size, Point3(-3., 4., 0.))
  c3 = gen_square_at(size, Point3(3., 5., 0.))
  c4 = gen_square_at(size, Point3(0., 7., 0.))
  return CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))
end

function gen_vert_line()
  size = 4.
  c1 = gen_square_at(size, Point3(0., 0., 0.))
  c2 = gen_square_at(size, Point3(3., 3., 0.))
  c3 = gen_square_at(size, Point3(0., 6., 0.))
  c4 = gen_square_at(size, Point3(-3., 9., 0.))
  return CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))
end

function gen_fork_join()
  size = 4.
  c1 = gen_square_at(size, Point3(-3., 0., 0.))
  c2 = gen_square_at(size, Point3(3., 0., 0.))
  c3 = gen_square_at(size, Point3(0., 3., 0.))
  c4 = gen_square_at(size, Point3(6., 3., 0.))
  return CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))
end

function gen_fork()
  size = 4.
  # c1 = gen_square_at(size, Point3(-3., 0.))
  c2 = gen_square_at(size, Point3(3., 0., 0.))
  c3 = gen_square_at(size, Point3(0., 3., 0.))
  c4 = gen_square_at(size, Point3(6., 3., 0.))
  return CSG.CUnion(c2, CSG.CUnion(c3, c4))
end

function gen_join()
  size = 4.
  c1 = gen_square_at(size, Point3(-3., 0., 0.))
  c2 = gen_square_at(size, Point3(3., 0., 0.))
  c3 = gen_square_at(size, Point3(0., 3., 0.))
  # c4 = gen_square_at(size, Point3(6., 3.))
  return CSG.CUnion(CSG.CUnion(c1, c2), c3)
end

function gen_crown()
  size = 4.
  treeify([
    gen_square_at(size, Point3(0., 0., 0.)),
    gen_square_at(size, Point3(-3., 3., -3.)),
    gen_square_at(size, Point3(3., 3., -3.)),
    gen_square_at(size, Point3(-3., 3., 3.)),
    gen_square_at(size, Point3(3., 3., 3.))
  ])
end


function gen_4_torus()
  size = 4.
  treeify([
    gen_square_at(size, Point3(0., 0., 0.)),
    gen_square_at(size, Point3(-3., 3., -3.)),
    gen_square_at(size, Point3(3., 3., -3.)),
    gen_square_at(size, Point3(-3., 3., 3.)),
    gen_square_at(size, Point3(3., 3., 3.)),
    # just a crown with a cube at the top to tie it back
    gen_square_at(size, Point3(0., 6., 0.)),
  ])
end


function gen_flat_line()
  size = 4.
  treeify([
    gen_square_at(size, Point3(-6., 0., -6.)),
    gen_square_at(size, Point3(-3., 0., -3.)),
    gen_square_at(size, Point3(0., 0., 0.)),
    gen_square_at(size, Point3(3., 0., 3.)),
    gen_square_at(size, Point3(6., 0., 6.)),
  ])
end

# TODO add support for this one!
function gen_kissing_torus()
  # where the components only barely touch
  size = 4.
  treeify([
    gen_square_at(size, Point3(0., 0., 0.)),
    gen_square_at(size, Point3(-3., 4., 0.)),
    gen_square_at(size, Point3(3., 4., 0.)),
    gen_square_at(size, Point3(0., 8., 0.)),
  ])
end


function gen_lattice()
  # where the components only barely touch
  size = 4.
  treeify([
    gen_square_at(size, Point3(-6., 0., 0.)),
    gen_square_at(size, Point3(6., 0., 0.)),
    gen_square_at(size, Point3(0., 0., 0.)),
    gen_square_at(size, Point3(-3., 3., 0.)),
    gen_square_at(size, Point3(3., 3., 0.)),
    gen_square_at(size, Point3(0., 6., 0.)),
    gen_square_at(size, Point3(-6., 6., 0.)),
    gen_square_at(size, Point3(6., 6., 0.))
  ])
end

# function gen_double_torus()
#   size = 4.
#   treeify([
#     gen_square_at(size, Point3(0., 2., 0.)),
#     gen_square_at(size, Point3(-3., 4., 0.)),
#     gen_square_at(size, Point3(3., 5., 0.)),
#     gen_square_at(size, Point3(0., 7., 0.)),
#     gen_square_at(size, Point3(3., 11., 0.)),
#     gen_square_at(size, Point3(-3., 11., 0.)),
#     gen_square_at(size, Point3(-3., 11., 0.)),
#     #TODO
#   ])
# end
