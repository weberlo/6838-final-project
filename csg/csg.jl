module CSG

# using LinearAlgebra
using GeometryBasics
using MLStyle
# using LightGraphs

# @data Shape begin
#   Circle(p :: Point2, r :: Int)
# end

@data CritKind begin
  Max
  Min
  Saddle
end

@data CritGeom begin
  CGPt(Point2)
  CGLine(a :: Point3, b :: Point3)
end

@data Expr begin
  CCircle(c :: Point3, r :: Int)
  CSquare(bot_left :: Point3, top_right :: Point3)
  CDiff(Expr, Expr)
  CUnion(Expr, Expr)
end

end
