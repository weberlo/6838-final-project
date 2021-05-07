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

y_coord(g :: CSG.CritGeom) = @match g begin
  CSG.CGLine(a, b) => begin
    @assert a[2] == b[2]
    a[2]
  end
  CSG.CGPt(p) => p[2]
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
      left_tagged = collect(map(collect(left_res.labels)) do pr
        (i, pt) = pr
        (side=:left, graph_idx=i, crit_pt=pt, intsct=nothing)
      end)
      right_tagged = collect(map(collect(right_res.labels)) do pr
        (i, pt) = pr
        (side=:right, graph_idx=i, crit_pt=pt, intsct=nothing)
      end)
      sort!(left_tagged, by = p -> y_coord(p.crit_pt.geom))
      sort!(right_tagged, by = p -> y_coord(p.crit_pt.geom))
      both_crits = vcat(left_tagged, right_tagged)
      sort!(both_crits, by = p -> y_coord(p.crit_pt.geom))

      res_graph = SimpleGraph(length(both_crits))

      for (i, tup) in collect(enumerate(both_crits))[2:end]
        # if tup.side == :left
        #   intsct = intersect_with(tup.crit_pt.geom, e2)
        # else # tup.side == :right
        #   intsct = intersect_with(tup.crit_pt.geom, e1)
        # end
        # println(intsct)

        # if isnothing(intsct)
          # construct graph normally
          if tup.side == :left
            my_geom = e1
          else # tup.side == :right
            my_geom = e2
          end
          if both_crits[i-1].side == :left
            prev_geom = e1
          else # tup.side == :right
            prev_geom = e2
          end
          # println()
          # println(both_crits[i-1].crit_pt.geom)
          # println(intersect_with(both_crits[i-1].crit_pt.geom, my_geom))
          if !isnothing(intersect_with(prev_geom, my_geom))
            add_edge!(res_graph, i-1, i)
          end
        # else
          # WE HAVE SOMETHING
        # end
      end
      println(res_graph)
      println(collect(edges(res_graph)))

      # for (i, p) in left_res.labels
      #   intsct = intersect_with(p.geom, e2)
      #   if !isnothing(intsct)
      #     push!(intscts, (:left, i, intsct))
      #   end
      # end
      # for (i, p) in right_res.labels
      #   intsct = intersect_with(p.geom, e1)
      #   if !isnothing(intsct)
      #     push!(intscts, (:right, i, intsct))
      #   end
      # end
      # sort!(intscts, by=y_coord)
      # println(intscts)

      # for p in values(right_res.labels)
      #   intsct = intersect_with(p.geom, e1)
      #   if !isnothing(intsct)
      #     println(intsct)
      #   end
      # end

      res_graph
    end
    _ => @assert false
  end
end

# TODO will need to generalize to allow arbitrary LHS geometry (e.g., arc,
# circle, sphere) and to return the resulting geometry.
function intersect_with(g :: CSG.CritGeom, e :: CSG.Expr)
  # TODO might need to return both the intersecting geometry *and* the remaining
  # geometry.
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
