include("./csg.jl")

# using .CSG
using LinearAlgebra
using GeometryBasics
using MLStyle
using LightGraphs

# @data Shape begin
#   Circle(p :: Point2, r :: Int)
# end

extents(s :: CSG.CSquare) = @match s begin
  CSG.CSquare(bot_left, top_right) =>
    (x1=bot_left[1], y1=bot_left[2], x2=top_right[1], y2=top_right[2])
end

function reeb_graph_of(expr :: CSG.Expr)
  @match expr begin
    CSG.CSquare(bot_left, top_right) => begin
      (x1, y1, x2, y2) = extents(expr)
      # min_pt_geom = CSG.CGLine(bot_left, bot_left + Vec2(top_right[1] - bot_left[1], 0.))
      # max_pt_geom = CSG.CGLine(top_right - Vec2(top_right[1] - bot_left[1], 0.), top_right)
      min_pt_geom = CSG.CGLine(Point2(x1, y1), Point2(x2, y1))
      max_pt_geom = CSG.CGLine(Point2(x1, y2), Point2(x2, y2))
      graph = SimpleGraph(2)
      add_edge!(graph, 1, 2)
      (labels=Dict(
        1 => (kind=Min, geom=min_pt_geom),
        2 => (kind=Max, geom=max_pt_geom)),
       graph=graph)
    end
    CSG.CUnion(e1, e2) => begin
      left_res = reeb_graph_of(e1)
      right_res = reeb_graph_of(e2)
      # println(left_res)
      # println(right_res)
      for p in values(left_res.labels)
        intsct = intersect_with(p.geom, e2)
        if !isnothing(intsct)
          println(intsct)
        end
      end
      for p in values(right_res.labels)
        intsct = intersect_with(p.geom, e1)
        if !isnothing(intsct)
          println(intsct)
        end
      end

      @assert false "why are there no intersections?"
      left_res
    end
    _ => @assert false
  end
end

# TODO will need to generalize to allow arbitrary LHS geometry (e.g., arc,
# circle, sphere) and to return the resulting geometry.
function intersect_with(g :: CSG.CritGeom, e :: CSG.Expr)
  @match g begin
    CSG.CGLine(a, b) => @match e begin
      CSG.CSquare(bot_left, top_right) => begin
        # TODO don't assume lines are axis-aligned
        (x1, y1, x2, y2) = extents(e)
        if b[1] < x1
          nothing
        elseif a[1] > x2
          nothing
        elseif b[2] < y1
          nothing
        elseif a[2] > y2
          nothing
        else
          CSG.CGLine(
            Point2(max(a[1], x1), max(a[2], y1)),
            Point2(min(b[1], x2), min(b[2], y2)))
        end
      end
      _ => @assert false
    end
    _ => @assert false
  end
end


# using LinearAlgebra
# using GeometryBasics

c1 = CSG.CSquare(Point2(0., 0.), Point2(1., 1.))
c2 = CSG.CSquare(Point2(0.5, 0.5), Point2(1.5, 1.5))
# println(c1)
# println(c2)
# println(typeof(CSG.CUnion(c1, c2)))
reeb_graph_of(CSG.CUnion(c1, c2))
