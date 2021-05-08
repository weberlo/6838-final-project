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

extents(e :: CSG.Expr) = @match e begin
  CSG.CSquare(bot_left, top_right) =>
    (x1=bot_left[1], y1=bot_left[2], x2=top_right[1], y2=top_right[2])
  CSG.CUnion(e1, e2) => begin
    left_res = extents(e1)
    right_res = extents(e2)
    (x1=min(left_res.x1, right_res.x1),
     x2=max(left_res.x2, right_res.x2),
     y1=min(left_res.y1, right_res.y1),
     y2=max(left_res.y2, right_res.y2))
  end
end

lengths(e :: CSG.Expr) = begin
  expr_extents = extents(e)
  (w=(expr_extents.x2 - expr_extents.x1),
   h=(expr_extents.y2 - expr_extents.y1))
end

function crit_pt_heights(e :: CSG.Expr)
  aux(e :: CSG.Expr) = @match e begin
    CSG.CSquare(bot_left, top_right) => [bot_left[2], top_right[2]]
    CSG.CUnion(e1, e2) => vcat(aux(e1), aux(e2))
  end
  sort(unique(aux(e)))
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

using Random

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
    println("TODO stop shuffling intervals once your algorithm is robust to it")
    # NOTE we need to handle arbitrary orderings of shapes for when we move to 3D.
    intervals = shuffle(intervals)

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

  # NOTE we assume the topology doesn't change in [y_{i-1} + eps, y_i - eps]

  heights = crit_pt_heights(e)
  adjacent_diffs = (circshift(heights, -1) - heights)[1:end-1]
  epsilon = minimum(filter(x -> x != 0., adjacent_diffs)) / 2.
  @assert epsilon != 0.
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
        # we need to map *all* members of the component (i.e., not just
        # the root) to `v`, because the component root can potentially change
        # across iterations.
        for memb in membs
          frontier[memb] = v
        end
      elseif length(src_comps) == 1
        # get vertex index of source component
        println("no downward change")
        v = prev_frontier[first(src_comps)]
        for memb in membs
          frontier[memb] = v
        end
      else  # length(src_comps) > 1
        println("at a join point")
        # topology changed and we're at a *join* point
        v = insert_vert(graph)
        for src_comp_root in src_comps
          # println("src_comp_root: $(src_comp_root)")
          # println("prev_frontier[src_comp_root]: $(prev_frontier[src_comp_root])")
          # println("v: $(v)")
          incr_edge_mult(graph, edge_mult,
            prev_frontier[src_comp_root], v)
        end
        for memb in membs
          frontier[memb] = v
        end
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
        println("at a local max")
        # at a local max, so kill this thread
        v = insert_vert(graph)
        # println(prev_frontier)
        incr_edge_mult(graph, edge_mult, prev_frontier[comp_root], v)
      elseif length(dst_comps) == 1
        println("no upward change")
        # do nothing (see GoodNotes)
        # TODO type up explanation from notes
        nothing
      else  # length(dst_comps) > 1
        # topology changed and we're at a *fork* point
        println("at a fork point")
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

  return graph, edge_mult
end

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

test_expr = gen_vert_line()

using Plots

function show_csg(e :: CSG.Expr)
  rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

  expr_extents = extents(e)
  println(expr_extents)
  x_range = expr_extents.x1:expr_extents.x2
  display(plot(x_range, LinRange(expr_extents.y1, expr_extents.y2, length(x_range)), opacity=0.))

  for prim in all_prims(e)
    lens = lengths(prim)
    display(plot!(rectangle(
      lens.w, lens.h,
      prim.bot_left[1], prim.bot_left[2]),
      opacity=.5))
  end
end

show_csg(test_expr)

graph, edge_mult = reeb_graph_of(test_expr)

using GraphPlot

function show_graph(graph, edge_mult)
  edgelabel = []
  for edge in edges(graph)
    push!(edgelabel, edge_mult[(edge.src,edge.dst)])
  end

  nodelabel = 1:nv(graph)
  gplot(graph, nodelabel=nodelabel, edgelabel=edgelabel,
    edgelabelc="white",
    edgelabeldistx=1.,
    edgelabeldisty=0.,
    arrowlengthfrac=0.1
  )
end

show_graph(graph, edge_mult)

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
