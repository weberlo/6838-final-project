using GeometryBasics
using LinearAlgebra
using DataStructures
using LightGraphs

include("mesh_utils.jl")

@enum CritType crit_min crit_max crit_saddle

# height direction vector
# crit_H = Vec3f0(0., 1., 0.)
crit_H = Vec3f0(0., 0., 1.)

function are_adjacent(v1_idx, v2_idx, T)
  for i = 1:size(T,1)
    if v1_idx in T[i,:] && v2_idx in T[i,:]
      return true
    end
  end
  return false
end

# TODO refactor into uf_util.jl
function calc_root_to_members(uf)
  res = DefaultDict(Vector)
  for elt in uf
    push!(res[find_root!(uf, elt)], elt)
  end
  return res
end

# TODO I think we can do only one round of classification if we restructure our algorithm
function classify_by_link(p_val, link, X, T, post_processing)
  num_neighbs = size(link, 1)
  sides = zeros(num_neighbs)
  num_degen = 0
  for (j, neighb_X_idx) in enumerate(link)
    np = X[neighb_X_idx,:]
    np_val = dot(np, crit_H)
    if np_val > p_val
      sides[j] = 1
    elseif np_val < p_val
      sides[j] = -1
    else
      # display(mesh(
      #     X, T,
      #     shading=true,
      #     transparency=true,
      #     figure=(resolution=(700, 1000),),
      #     color = (:red, 0.1)
      # ))
      # mesh!(Sphere(Point3f0(p), 0.5f0), transparency=false)
      # mesh!(Sphere(Point3f0(np), 0.5f0), transparency=false)
      # @assert false "multiple vertices with same value!"
      sides[j] = 0
      num_degen += 1
    end
  end
  if sum(sides) + num_degen == num_neighbs
    if num_degen > 0 && post_processing
      return nothing
    else
      # push!(crit_pts, (type=crit_min, idx=i, val=p_val))
      return (type=crit_min, val=p_val)
    end
  elseif sum(sides) - num_degen == -num_neighbs
    if num_degen > 0 && post_processing
      return nothing
    else
      # push!(crit_pts, (type=crit_max, idx=i, val=p_val))
      return (type=crit_max, val=p_val)
    end
  else
    group_count = 0
    last = sides[1]
    for j = 2:(num_neighbs+1)
      curr = sides[mod1(j, num_neighbs)]
      if sides[mod1(j, num_neighbs)] != last && last != 0
        group_count += 1
      end
      last = curr
    end
    @assert (group_count % 2 == 0) "group count was not multiple of 2!"
    if group_count > 2
      multiplicity = (group_count - 2) / 2
      # push!(crit_pts, (type=crit_saddle, idx=i, val=p_val, mult=multiplicity))
      return (type=crit_saddle, val=p_val, mult=multiplicity)
    end
  end
  return nothing
end

function crit_areas(X, T)
  nv = size(X, 1)
  nt = size(T, 1)
  links = vert_links(X, T)

  # form critical areas from critical points
  crit_areas = DisjointSets()
  for i = 1:size(X,1)
    push!(crit_areas, i)
  end

  for i = 1:size(T,1)
    p_vals = []
    for j = 1:3
      push!(p_vals, dot(X[T[i,j],:], crit_H))
    end
    if p_vals[1] == p_vals[2] && p_vals[1] == p_vals[3]
      union!(crit_areas, T[i,1], T[i,2])
      union!(crit_areas, T[i,2], T[i,3])
    end
  end
  # for i = 1:size(X,1)
  #   for j = 1:size(X,1)
  #     if i == j
  #       continue
  #     end
  #     p1 = X[i,:]
  #     p2 = X[j,:]
  #     p1_val = dot(p1, crit_H)
  #     p2_val = dot(p2, crit_H)
  #     if are_adjacent(i, j, T) && p1_val == p2_val
  #       union!(crit_areas, i, j)
  #     end
  #   end
  # end

  root_to_membs = calc_root_to_members(crit_areas)

  # crit_pts = []
  # for i = 1:nv
  #   p = X[i,:]
  #   p_val = dot(p, crit_H)
  #   link = links[i]
  #   classified_pt = classify_by_link(p_val, link, X, T, false)

  #   if classified_pt !== nothing
  #     push!(crit_pts, (classified_pt..., idx=i))
  #   end
  # end


  # TODO it looks like we might be losing the NamedTuple multiplicity field here?
  # res = map(
  #   pts -> (type=pts[1].type, idxs=map(p -> p.idx, pts), val=pts[1].val),
  #   values(root_to_membs))
  res = []
  for crit_area_idxs in values(root_to_membs)
    # println("crit_area: $crit_area")
    link = area_link(crit_area_idxs, X, T)
    # println(link)
    p_val = dot(X[crit_area_idxs[1],:], crit_H)
    # println("link: $link")
    # println("p_val: $p_val")
    classified_area = classify_by_link(p_val, link, X, T, true)
    if classified_area !== nothing
      # println("classified_area: $classified_area")
      push!(res, (classified_area..., idxs=crit_area_idxs))
    end
    # println(classified_area === nothing)
  end
  # sort by the critical values
  sort!(res, by = x -> x.val)
  return res
end

function area_link(crit_area_idxs, X, T)
  println("HERE")
  println(crit_area_idxs)
  star = Set()
  for idx in crit_area_idxs
    for i = 1:size(T, 1)
      if idx in T[i,:]
        star = star ∪ Set(T[i,:])
      end
    end
  end
  inner = Set(crit_area_idxs)
  unordered_link = setdiff(star, inner)
  # println(inner)
  # println(unordered_link)

  # if 3 in crit_area.idxs
  #   for vert_idx in inner
  #     mesh!(Sphere(Point3f0(X[vert_idx,:]), 0.5f0), transparency=false, color=:green)
  #   end
  #   for vert_idx in unordered_link
  #     mesh!(Sphere(Point3f0(X[vert_idx,:]), 0.5f0), transparency=false, color=:blue)
  #   end
  # end

  link = vert_link_general(inner, unordered_link, X, T)

  # NOTE: code below plots the link with redness corresponding to order in link
  # n = size(link,1)
  # println("link: $link")
  # for (i, vert_idx) in enumerate(link)
  #   mesh!(Sphere(Point3f0(X[vert_idx,:]), 0.5f0), transparency=false, color=RGBAf0((i-1)/float(n), 0., 0., 1.))
  # end

  return link
end


"""
Computes the *lower* critical set of a critical point `crit_pt`.
"""
function crit_set(crit_area, X, T)
  if crit_area.type == crit_min || crit_area.type == crit_max
    # critical set of min and max is point itself
    #
    # TODO Hajij and Rosen say we should include the closure of the star of these
    # vertices. should we actually?
    return Set(crit_area.idxs)
  end
  # we know we're working with a saddle point now
  curr_link = area_link(crit_area.idxs, X, T)
  println("curr_link: $curr_link")
  to_visit = Queue{Int32}()
  for idx in curr_link
    if dot(X[idx,:], crit_H) < crit_area.val
      enqueue!(to_visit, idx)
    end
  end

  visited = Set(first(to_visit))
  res = Set()
  while !isempty(to_visit)
    curr_vert_idx = dequeue!(to_visit)
    # curr_link = links[curr_vert_idx]

    # we only add the current vertex to the result if they
    # have a neighbor *above* the tangent plane.
    any_above = false

    # at each iter, we check if we've already visited a neighbor
    # (i.e., a vertex in the link) who's *below* the tangent plane.
    # if not, we add them to the `to_visit` queue.
    for j = 1:size(curr_link, 1)
      curr_neighb_idx = curr_link[j]
      curr_neighb_val = dot(X[curr_neighb_idx, :], crit_H)
      if curr_neighb_val < crit_area.val && !(curr_neighb_idx in visited)
        # unconditionally add the vertex to visited
        push!(visited, curr_neighb_idx)
        enqueue!(to_visit, curr_neighb_idx)
      end
      if curr_neighb_val > crit_area.val
        any_above = true
      end
    end

    if any_above
      push!(res, curr_vert_idx)
    end
  end
  return res
end


function reeb_graph(X, T)
  links = vert_links(X, T)
  crits = crit_areas(X, T)
  num_crits = size(crits, 1)
  crit_sets = Array{Set}(undef, num_crits)
  for i = 1:num_crits
    crit_sets[i] = crit_set(crits[i], X, T)
  end

  res = SimpleGraph(num_crits)

  # TODO we used to have a -1 here. make sure we didn't need it.
  for i = 1:num_crits
    curr_crit = crits[i]
    curr_crit_tag = curr_crit.type

    curr_crit_vert_idx = curr_crit.idxs[1]

    seed_idxs = []

    if curr_crit_tag == crit_min
      # only need one ascending path if at minimum
      push!(seed_idxs, curr_crit_vert_idx)
    elseif curr_crit_tag == crit_max
      # discount all maxima from this *outer* loop (but we'll still visit them
      # as destinations in the inner loops)
      continue
    else
      # TODO need to restructure code to compute *multiple* ascending paths for saddles

      # TODO massage this code into current code to compute seed idxs (an idx
      # for the highest vertex in each group)
      p = X[curr_crit_vert_idx,:]

      # link = links[curr_crit_vert_idx]
      link = area_link(curr_crit.idxs, X, T)

      num_neighbs = size(link, 1)
      link_vals = zeros(num_neighbs)
      p_val = dot(p, crit_H)
      for (j, neighb_X_idx) in enumerate(link)
        np = X[neighb_X_idx,:]
        np_val = dot(np, crit_H)
        link_vals[j] = np_val - p_val
      end

      # shift to start with a negative vertex, so we have all of the positive
      # groups contiguous
      offs = findfirst(map(x -> x < 0, link_vals)) - 1
      link_vals = circshift(link_vals, -offs)
      last_val = link_vals[1]
      group_best_val = -1
      group_best_idx = 0
      for j = 2:(num_neighbs+1)
        curr_idx = mod1(j, num_neighbs)
        curr_val = link_vals[curr_idx]
        # println(sign(curr_val))
        if sign(curr_val) < 0 && sign(last_val) > 0
          # println(link_vals[group_best_idx])
          # println(link)
          push!(seed_idxs, link[group_best_idx])
          group_best_val = -1
          group_best_idx = 0
        elseif sign(curr_val) > 0 && curr_val > group_best_val
          group_best_val = curr_val
          group_best_idx = curr_idx
        end
        last_val = curr_val
      end
    end

    # follow each ascending path (each given by a seed idx) upwards to the next
    # critical set
    for seed_idx in seed_idxs
      curr_vert_idx = seed_idx
      # println(seed_idxs)
      curr_val = dot(X[curr_vert_idx, :], crit_H)
      visited = Set()

      for j = i+1:num_crits
        next_crit_set = crit_sets[j]
        next_crit_val = crits[j].val
        println("next: $next_crit_val")
        println("next_crit_set: $next_crit_set")
        while curr_val < next_crit_val && !(curr_vert_idx in next_crit_set)
          was_update = false
          # greedily climb to highest neighbor in link
          curr_link = links[curr_vert_idx]
          for k = 1:size(curr_link, 1)
            # in here, we're mutably updating the current vertex + value to find
            # the highest vertex.
            neighb_idx = curr_link[k]
            neighb_val = dot(X[neighb_idx, :], crit_H)
            if neighb_val >= curr_val # && !(neighb_idx in visited)
              # println("neighb_val: $neighb_val")
              # println("neighb_idx: $neighb_idx")
              # visited = visited ∪ neighb_idx
              curr_val = neighb_val
              curr_vert_idx = neighb_idx
              was_update = true
            end
          end
          # println(curr_val)
          if !was_update
            # for vert_idx in inner
              mesh!(Sphere(Point3f0(X[curr_vert_idx,:]), 0.5f0), transparency=false, color=:green)
            # end
          end
          @assert was_update "no higher vertex found in link!"
        end
        if curr_vert_idx in next_crit_set
          if curr_crit_tag == crit_saddle
            println("FOUND SADDDLELEEE")
          end
          add_edge!(res, i, j)
          break
        end
      end
    end
  end
  return res
end


function filter_crit_tag(crits, tag)
  return filter(x -> x[1] == tag, crits)
end

function crit_mins(crits)
  return filter_crit_tag(crits, crit_min)
end

function crit_maxs(crits)
  return filter_crit_tag(crits, crit_max)
end

function crit_saddles(crits)
  return filter_crit_tag(crits, crit_saddle)
end

function crit_data(crits)
  # discard tag
  return map(x -> x[2], crits)
end

function crit_idxs(crits)
  return map(x -> x[1], crit_data(crits))
end

function crit_vals(crits)
  return map(x -> x[2], crit_data(crits))
end

function crit_mults(crits)
  @assert all(map(x -> x[1] == crit_saddle, crits)) "attempt to extract multiplicity from non-saddle point"
  return map(x -> x[3], crit_data(crits))
end
