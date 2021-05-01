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
  CGLine(a :: Point2, b :: Point2)
end

@data Expr begin
  CCircle(c :: Point2, r :: Int)
  CSquare(bot_left :: Point2, top_right :: Point2)
  CDiff(Expr, Expr)
  CUnion(Expr, Expr)
end

end
