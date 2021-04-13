using MeshCat
# vis = Visualizer()
# open(vis)

using GeometryBasics
using CoordinateTransformations

using FileIO

# pointcloud = PointCloud([[x, 0 + 0.01 * randn(), 0.5] for x in range(-1, stop=1, length=1000)])
# setobject!(vis[:pointcloud], pointcloud)

# setobject!(vis, HyperSphere(Point(0., 0, 0), 0.5))
# setobject!(vis, GeometryBasics.Simplex(Vec(0., 0, 0), Vec(0., 1, 0), Vec(0., 1, 1)))
# settransform!(vis, Translation(-0.5, -0.5, 0))
# settransform!(vis, Translation(0., 0., 0.))

# show(HyperRectangle(Vec(0., 0, 0), Vec(1., 1, 1)))
# show(HyperSphere(Point(0., 0, 0), ))
# show(GeometryBasics.Simplex(Vec(0., 0, 0), Vec(0., 1, 0), Vec(0., 1, 1)))

# path = joinpath(@__DIR__, "..", "test", "data", "meshes", "cube.dae")
# setobject!(vis["meshes", "dae_file_geometry"],
#     MeshFileGeometry("semicone.stl"))

# settransform!(vis["meshes", "dae_file_geometry"],
#     Translation(0.0, 1.25, 0.0))
# settransform!(vis["meshes", "dae_file_geometry"],
#     Scale(0.2, 0.2, 0.2))

