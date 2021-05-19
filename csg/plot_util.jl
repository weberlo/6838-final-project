using Plots
using GraphPlot

function show_csg(e :: CSG.Expr)
  rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

  expr_extents = extents(e)
  println(expr_extents)
  x_range = expr_extents.x1:expr_extents.x2
  display(plot(
    x_range,
    LinRange(expr_extents.y1, expr_extents.y2, length(x_range)),
    opacity=0.
  ))

  for prim in all_prims(e)
    lens = lengths(prim)
    display(plot!(rectangle(
      lens.w, lens.h,
      prim.bot_left[1], prim.bot_left[2]),
      opacity=.5))
  end
end


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
