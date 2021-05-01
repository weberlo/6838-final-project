include("csg.jl")

using LinearAlgebra
using GeometryBasics
using MLStyle
using LightGraphs

# @data Shape begin
#   Circle(p :: Point2, r :: Int)
# end

function reeb_graph_of(expr :: Expr)
  @match expr begin
    Circle(c, r) => begin
      min_pt_geom = Pt(c + Vec2(0, -r))
      max_pt_geom = Pt(c + Vec2(0, r))
      (labels=Dict(1 => (kind=Min, geom=min_pt_geom), 2 => (kind=Max, geom=max_pt_geom)),
       graph=SimpleGraph(2))
    end
    Diff(e1, e2) => begin
      left_res = reeb_graph_of(e1)
      right_res = reeb_graph_of(e2)
      for p in values(left_res.labels)
        println(intersect_with(p.geom, e2))
      end
      for p in values(right_res.labels)
        println(intersect_with(p.geom, e1))
      end
      @assert false "why are there no intersections?"
      left_res
    end
  end
end

# TODO will need to generalize to allow arbitrary LHS geometry (e.g., arc,
# circle, sphere) and to return the resulting geometry.
function intersect_with(p :: CritGeom, e :: Expr)
  @match p begin
    Pt(p) => @match e begin
      Circle(c, r) => norm(p - c) < r
      Diff(e1, e2) => @assert false
    end
  end
end


using .CSG
using LinearAlgebra
using GeometryBasics

c1 = CSG.Circle(Point2(0, 0), 10)
c2 = CSG.Circle(Point2(5, 0), 10)
println(c1)
println(c2)
CSG.reeb_graph_of(CSG.Diff(c1, c2))
