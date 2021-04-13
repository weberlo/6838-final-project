using GeometryBasics
using LinearAlgebra


function one_rings(X, T)
  nv = size(X, 1)
  nt = size(T, 1)
  res = Array{Array{Int64}}(undef, nv)
  for i = 1:nv
    res[i] = []
    for j = 1:nt
      if i in T[j,:]
        push!(res[i], j)
      end
    end
  end
  return res
end


function face_normals(X, T)
  nt = size(T, 1)
  FN = zeros(size(T))
  for j = 1:nt
    tri_idxs = T[j,:]
    tri = [X[tri_idxs[1],:], X[tri_idxs[2],:], X[tri_idxs[3],:]]
    A = tri[2] - tri[1]
    B = tri[3] - tri[1]
    n = cross(A, B)
    n = n ./ norm(n)
    FN[j,:] = n
  end
  return FN
end


function face_areas(X, T)
  nt = size(T, 1)
  res = zeros(size(T, 1))
  for j = 1:nt
    tri_idxs = T[j,:]
    tri = [X[tri_idxs[1],:], X[tri_idxs[2],:], X[tri_idxs[3],:]]
    A = tri[2] - tri[1]
    B = tri[3] - tri[1]
    n = norm(cross(A, B)) / 2
    res[j] = n
  end
  return res
end


function barycentricArea(X, T)
  nv = size(X, 1)
  nt = size(T, 1)
  ba = zeros(nv)
  for i = 1:nt
    t = T[i,:]
    p_idx = t[1]
    q_idx = t[2]
    r_idx = t[3]
    p = X[p_idx,:]
    q = X[q_idx,:]
    r = X[r_idx,:]
    area = norm(cross(q - p, r - p)) / 2.0
    ba[p_idx] += area / 3.0
    ba[q_idx] += area / 3.0
    ba[r_idx] += area / 3.0
  end
  return ba
end


function vertex_normals(X, T)
  nv = size(X, 1)
  nt = size(T, 1)
  FN = face_normals(X, T)
  FA = face_areas(X, T)
  rings = one_rings(X, T)

  res = zeros((nv, 3))
  for i = 1:nv
    ring = rings[i]
    total_area = 0.
    total_norm = Vec3f0(0., 0, 0)
    for neighb_idx in ring
      curr_area = FA[neighb_idx]
      curr_norm = FN[neighb_idx,:]
      total_norm += (curr_area .* curr_norm)
      total_area += curr_area
    end
    res[i,:] = total_norm ./ total_area
  end
  return res
end


"""
Computes link (all vertices in 1-ring except for the center vertex) of each
  vertex.
"""
function vert_links(X, T)
  nv = size(X, 1)
  rings = one_rings(X, T)
  res = Array{Array{Int64}}(undef, nv)
  for i = 1:nv
    ring = rings[i]
    tris = [T[idx,:] for idx in ring]
    link_dict = Dict()
    for tri in tris
      shift_amt = - (findall(x -> x == i, tri)[1] - 1)
      tri = circshift(tri, shift_amt)
      link_dict[tri[2]] = tri[3]
    end
    curr = first(link_dict)[1]
    link = [curr]
    for j = 1:size(tris, 1)
      curr = link_dict[curr]
      push!(link, curr)
    end
    res[i] = link
  end
  return res
end
