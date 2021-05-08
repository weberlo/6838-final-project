include("./csg.jl")

using GeometryBasics


function gen_square_at(s :: Number, center :: Point2)
  return CSG.CSquare(
    Point2(center[1] - s / 2., center[2] - s / 2.),
    Point2(center[1] + s / 2., center[2] + s / 2.))
end

function gen_torus()
  # 2D version of torus-ish shape from GoodNotes
  size = 4.
  c1 = gen_square_at(size, Point2(0., 2.))
  c2 = gen_square_at(size, Point2(-3., 4.))
  c3 = gen_square_at(size, Point2(3., 5.))
  c4 = gen_square_at(size, Point2(0., 7.))
  return CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))
end

function gen_vert_line()
  size = 4.
  c1 = gen_square_at(size, Point2(0., 0.))
  c2 = gen_square_at(size, Point2(3., 3.))
  c3 = gen_square_at(size, Point2(0., 6.))
  c4 = gen_square_at(size, Point2(-3., 9.))
  return CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))
end

function gen_fork_join()
  size = 4.
  c1 = gen_square_at(size, Point2(-3., 0.))
  c2 = gen_square_at(size, Point2(3., 0.))
  c3 = gen_square_at(size, Point2(0., 3.))
  c4 = gen_square_at(size, Point2(6., 3.))
  return CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))
end

function gen_fork()
  size = 4.
  # c1 = gen_square_at(size, Point2(-3., 0.))
  c2 = gen_square_at(size, Point2(3., 0.))
  c3 = gen_square_at(size, Point2(0., 3.))
  c4 = gen_square_at(size, Point2(6., 3.))
  return CSG.CUnion(c2, CSG.CUnion(c3, c4))
end

function gen_join()
  size = 4.
  c1 = gen_square_at(size, Point2(-3., 0.))
  c2 = gen_square_at(size, Point2(3., 0.))
  c3 = gen_square_at(size, Point2(0., 3.))
  # c4 = gen_square_at(size, Point2(6., 3.))
  return CSG.CUnion(CSG.CUnion(c1, c2), c3)
end
