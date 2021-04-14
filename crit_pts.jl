using GeometryBasics
using LinearAlgebra
using DataStructures
using LightGraphs

include("mesh_utils.jl")

@enum CritType crit_min crit_max crit_saddle

# height direction vector
# crit_H = Vec3f0(0., 1., 0.)
crit_H = Vec3f0(0., 0., 1.)

function crit_points(X, T)
  nv = size(X, 1)
  nt = size(T, 1)
  links = vert_links(X, T)

  res = []
  for i = 1:nv
    p = X[i,:]
    num_neighbs = size(links[i], 1)
    sides = zeros(num_neighbs)
    p_val = dot(p, crit_H)
    for (j, neighb_X_idx) in enumerate(links[i])
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
        sides[j] = -1
      end
    end
    if sum(sides) == num_neighbs
      push!(res, (crit_min, (i, p_val)))
    elseif sum(sides) == -num_neighbs
      push!(res, (crit_max, (i, p_val)))
    else
      group_count = 0
      last = sides[1]
      for j = 2:(num_neighbs+1)
        curr = sides[mod1(j, num_neighbs)]
        if sides[mod1(j, num_neighbs)] != last
          group_count += 1
        end
        last = curr
      end
      @assert (group_count % 2 == 0) "group count was not multiple of 2!"
      if group_count > 2
        multiplicity = (group_count - 2) / 2
        push!(res, (crit_saddle, (i, p_val, multiplicity)))
      end
    end
  end
  # sort by the critical values
  sort!(res, by = x -> x[2][2])
  return res
end


"""
Computes the *lower* critical set of a critical point `crit_pt`.
"""
function crit_set(crit_pt, X, T)
  if crit_pt[1] == crit_min || crit_pt[1] == crit_max
    # critical set of min and max is point itself
    #
    # TODO they say we should include the closure of the star of these
    # vertices. should we actually?
    return Set(crit_pt[2][1])
  end
  # we know we're working with a saddle point now
  crit_data = crit_pt[2]
  crit_idx = crit_data[1]
  crit_val = crit_data[2]
  links = vert_links(X, T)
  curr_link = links[crit_idx]
  to_visit = Queue{Int32}()
  for idx in curr_link
    if dot(X[idx,:], crit_H) < crit_val
      enqueue!(to_visit, idx)
    end
  end

  visited = Set(first(to_visit))
  res = Set()
  # println(to_visit)
  # println(isempty(to_visit))
  while !isempty(to_visit)
    curr_vert_idx = dequeue!(to_visit)
    curr_link = links[curr_vert_idx]
    # we only add the current vertex to the result if they
    # have a neighbor *above* the tangent plane.
    any_above = false

    # at each iter, we check if we've already visited a neighbor
    # (i.e., a vertex in the link) who's *below* the tangent plane.
    # if not, we add them to the `to_visit` queue.
    for j = 1:size(curr_link, 1)
      curr_neighb_idx = curr_link[j]
      curr_neighb_val = dot(X[curr_neighb_idx, :], crit_H)
      if curr_neighb_val < crit_val && !(curr_neighb_idx in visited)
        # unconditionally add the vertex to visited
        push!(visited, curr_neighb_idx)
        enqueue!(to_visit, curr_neighb_idx)
      end
      if curr_neighb_val > crit_val
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
  crits = crit_points(X, T)
  num_crits = size(crits, 1)
  crit_sets = Array{Set}(undef, num_crits)
  # all_crit_sets = Array{Array{Int64}}(undef, )
  for i = 1:num_crits
    crit_sets[i] = crit_set(crits[i], X, T)
  end
  # println(crit_sets)

  res = SimpleGraph(num_crits)

  # TODO I don't think we need the -1
  for i = 1:num_crits-1
    curr_crit = crits[i]
    curr_crit_tag = curr_crit[1]

    curr_crit_data = curr_crit[2]
    curr_crit_vert_idx = curr_crit_data[1]

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
      link = links[curr_crit_vert_idx]
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
      println(link_vals)
      link_vals = circshift(link_vals, -offs)
      println(link_vals)
      last_val = link_vals[1]
      group_best_val = -1
      group_best_idx = 0
      for j = 2:(num_neighbs+1)
        curr_idx = mod1(j, num_neighbs)
        curr_val = link_vals[curr_idx]
        # println(sign(curr_val))
        if sign(curr_val) < 0 && sign(last_val) > 0
          println(link_vals[group_best_idx])
          push!(seed_idxs, link[group_best_idx])
          group_best_val = -1
          group_best_idx = 0
        elseif sign(curr_val) > 0 && curr_val > group_best_val
          group_best_val = curr_val
          group_best_idx = curr_idx
        end
        last_val = curr_val
      end
      # println(size(seed_idxs))
      # println(seed_idxs)
      # println(seed_idxs)
      # println()
    end

    # follow each ascending path (each given by a seed idx) upwards to the next
    # critical set
    for seed_idx in seed_idxs
      curr_vert_idx = seed_idx
      curr_val = dot(X[curr_vert_idx, :], crit_H)

      for j = i+1:num_crits
        next_crit_set = crit_sets[j]
        next_crit_val = crits[j][2][2]
        while curr_val < next_crit_val && !(curr_vert_idx in next_crit_set)
          was_update = false
          # greedily climb to highest neighbor in link
          curr_link = links[curr_vert_idx]
          for k = 1:size(curr_link, 1)
            # in here, we're mutably updating the current vertex + value to find
            # the highest vertex.
            neighb_idx = curr_link[k]
            neighb_val = dot(X[neighb_idx, :], crit_H)
            if neighb_val > curr_val
              curr_val = neighb_val
              curr_vert_idx = neighb_idx
              was_update = true
            end
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
