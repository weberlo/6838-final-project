function plot_spheres(X, idxs, color, radius=0.5f0)
  for X_idx in idxs
    # mesh!(Sphere(Point3f0(X[X_idx,:]), 0.5f0), transparency=false, color=color)
    plot_sphere_by_pos(Point3f0(X[X_idx,:]), color, radius)
  end
end

function plot_sphere_by_pos(pos, color, radius=0.5f0)
  mesh!(Sphere(pos, radius), transparency=false, color=color)
end
