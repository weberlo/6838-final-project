include("./csg.jl")

using LinearAlgebra
using GeometryBasics
using DataStructures
using MLStyle
using LightGraphs
using Combinatorics

# println(c1)
# println(c2)
# println(typeof(CSG.CUnion(c1, c2)))
# reeb_graph_of()

extents(s :: CSG.CSquare) = @match s begin
  CSG.CSquare(bot_left, top_right) =>
    (x1=bot_left[1], y1=bot_left[2], x2=top_right[1], y2=top_right[2])
end

function crit_pt_heights(e :: CSG.Expr)
  aux(e :: CSG.Expr) = @match e begin
    CSG.CSquare(bot_left, top_right) => [bot_left[2], top_right[2]]
    CSG.CUnion(e1, e2) => vcat(aux(e1), aux(e2))
  end
  sort(aux(e))
end

all_prims(e :: CSG.Expr) = @match e begin
  CSG.CSquare(bot_left, top_right) => [CSG.CSquare(bot_left, top_right)]
  CSG.CUnion(e1, e2) => vcat(all_prims(e1), all_prims(e2))
end

# num_crit_pts(e :: CSG.Expr) = @match e begin
#   CSG.CSquare(_) => 2
#   CSG.CUnion(e1, e2) => num_crit_pts(e1) + num_crit_pts(e2)
# end
num_crit_pts(e :: CSG.Expr) = length(crit_pt_heights(e))


function reeb_graph_of(e :: CSG.Expr)
  function collisions_at(y :: Number, prims)
    collisions = []
    for (i, prim) in enumerate(prims)
      if prim.bot_left[2] > y || prim.top_right[2] < y
        continue
      else
        push!(collisions, (i, (left=prim.bot_left[1], right=prim.top_right[1])))
      end
    end
    return collisions
  end

  function comp_union_find(intervals)
    # println(intervals)
    res = DisjointSets()
    for (i, _) in intervals
      push!(res, i)
    end
    for ((i1, s1), (i2, s2)) in combinations(intervals, 2)
      if s1.right < s2.left || s1.left > s2.right
        continue
      else
        union!(res, i1, i2)
      end
    end
    return res
  end

  # function all_comps(uf)
  #   res = Set()
  #   for elt in uf
  #     res = res ∪ find_root!(uf, elt)
  #   end
  #   return res
  # end

  function calc_root_to_members(uf)
    res = DefaultDict(Vector{Int})
    for elt in uf
      push!(res[find_root!(uf, elt)], elt)
    end
    return res
  end

  function incr_edge_mult(g, edge_mult, i, j)
    if add_edge!(g, i, j)
      @assert edge_mult[(i, j)] == 0
    end
    edge_mult[(i, j)] += 1
  end

  function insert_vert(g)
    @assert add_vertex!(g)
    println("created vert $(nv(g))")
    return nv(g)
  end

  # prev_frontier = []
  # prev_frontier_uf = DisjointSets()

  # NOTE we assume the topology doesn't change in [y_{i-1} + eps, y_i - eps]

  heights = crit_pt_heights(e)
  adjacent_diffs = (circshift(heights, -1) - heights)[1:end-1]
  epsilon = minimum(adjacent_diffs) / 2.
  # println(epsilon)

  prims = all_prims(e)

  prev_frontier = Dict()
  graph = SimpleGraph()
  edge_mult = DefaultDict(0)
  for crit_pt_idx = 1:num_crit_pts(e)
    println("[crit_pt_idx $crit_pt_idx]")
    curr_y_pos = heights[crit_pt_idx] + epsilon
    curr_y_neg = heights[crit_pt_idx] - epsilon

    frontier = Dict()

    neg_collns = collisions_at(curr_y_neg, prims)
    pos_collns = collisions_at(curr_y_pos, prims)
    # println("neg_collns: $neg_collns")
    # println("pos_collns: $pos_collns")

    neg_uf = comp_union_find(neg_collns)
    pos_uf = comp_union_find(pos_collns)

    # println(in_same_set(neg_uf, 2, 3))
    # println(in_same_set(pos_uf, 2, 3))
    # vert_to_pos = Dict()

    pos_root_to_membs = calc_root_to_members(pos_uf)
    for (comp_root, membs) in pos_root_to_membs
      src_comps = Set()
      for memb in membs
        if memb in neg_uf
          src_comps = src_comps ∪ find_root!(neg_uf, memb)
        end
      end

      if length(src_comps) == 0
        # at a local min, so spawn new vertex
        println("at local min")
        v = insert_vert(graph)
        # TODO we might need to map *all* members of the component (i.e., not just
        # the root) to `v`, because the component root can potentially change
        # across iterations.
        frontier[comp_root] = v
      elseif length(src_comps) == 1
        # get vertex index of source component
        # println("src_comps: $src_comps")
        # println("comp_root: $comp_root")
        # println("prev_frontier: $prev_frontier")
        v = prev_frontier[first(src_comps)]
        # TODO same here as above
        frontier[comp_root] = v
      else  # length(src_comps) > 1
        # topology changed and we're at a *join* point
        v = insert_vert(graph)
        for src_comp_root in src_comps
          # println("src_comp_root: $(src_comp_root)")
          # println("prev_frontier[src_comp_root]: $(prev_frontier[src_comp_root])")
          # println("v: $(v)")
          incr_edge_mult(graph, edge_mult,
            prev_frontier[src_comp_root], v)
        end
        frontier[comp_root] = v
      end
    end

    neg_root_to_membs = calc_root_to_members(neg_uf)
    for (comp_root, membs) in neg_root_to_membs
      dst_comps = Set()
      for memb in membs
        if memb in pos_uf
          dst_comps = dst_comps ∪ find_root!(pos_uf, memb)
        end
      end

      if length(dst_comps) == 0
        # at a local max, so kill this thread
        v = insert_vert(graph)
        # println(prev_frontier)
        incr_edge_mult(graph, edge_mult, prev_frontier[comp_root], v)
      elseif length(dst_comps) == 1
        # do nothing (see GoodNotes)
        # TODO type up explanation from notes
        nothing
      else  # length(dst_comps) > 1
        # topology changed and we're at a *fork* point
        println("at fork point")
        v = insert_vert(graph)
        incr_edge_mult(graph, edge_mult, prev_frontier[comp_root], v)
        for dst_comp_root in dst_comps
          @assert (frontier[dst_comp_root] == prev_frontier[comp_root])  # see GoodNotes
          frontier[dst_comp_root] = v
        end
      end
    end

    println()
    prev_frontier = frontier
  end

  # println("neg comps: $(all_comps(neg_uf))")
  # println("pos comps: $(all_comps(pos_uf))")
  # for root in neg_uf end
  # println(in_same_set(pos_uf, 1, 5))

  # println(neg_uf)
  # println(pos_uf)
  return graph, edge_mult
end


# 2D version of torus-ish shape from GoodNotes
c1 = CSG.CSquare(Point2(-2., 0.), Point2(2., 4.))
c2 = CSG.CSquare(Point2(-5., 2.), Point2(-1., 6.))
c3 = CSG.CSquare(Point2(1., 3.), Point2(5., 7.))
c4 = CSG.CSquare(Point2(-2., 5.), Point2(2., 9.))
test_expr = CSG.CUnion(CSG.CUnion(c1, c2), CSG.CUnion(c3, c4))

graph, edge_mult = reeb_graph_of(test_expr)
edgelabel = []
for edge in edges(graph)
  # println(edge.src)
  # println(edge.dst)
  push!(edgelabel, edge_mult[(edge.src,edge.dst)])
end

using GraphPlot
nodelabel = 1:nv(graph)
gplot(graph, nodelabel=nodelabel, edgelabel=edgelabel,
  edgelabelc="white",
  edgelabeldistx=1.,
  edgelabeldisty=0.,
  arrowlengthfrac=0.1
  )

# function plot_graph(vert_to_pos, graph)
#   plot_current_mesh()
#   # bounds = Node([
#   #   Point3f0(0., 0., 0.),
#   #   Point3f0(15., 0., 0.),
#   #   Point3f0(-15., 0., 0.),
#   #   Point3f0(0., 0., 20.),
#   #   Point3f0(0., 0., -30.),
#   #   Point3f0(0., 20., 0.),
#   #   Point3f0(0., -15., 0.),
#   # ])
#   # display(lines(bounds, linewidth=0, transparency=true, color=(:black, 0.0),
#   #     shading=true,
#   #     figure=(resolution=(700, 1000),),
#   # ))
#   crits = crit_points(X, T)
#   plot_spheres(crit_idxs(crit_mins(crits)), :yellow)
#   plot_spheres(crit_idxs(crit_maxs(crits)), :green)
#   plot_spheres(crit_idxs(crit_saddles(crits)), :blue)
# end
