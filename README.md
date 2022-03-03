# CS272: Geometric Modeling
## Requirements and build:
This is a viewer application, using [GLFW](https://www.glfw.org), [Dear ImGui](https://github.com/ocornut/imgui), [Eigen](https://eigen.tuxfamily.org), and [OpenMesh](https://www.openmesh.org), and based on the [libigl](https://libigl.github.io) viewer.
All dependencies are included in the `external` directory.
You can build this application with `cmake`. It has been tested under **Windows** and **Linux**.

## functionality:
1. B-spline/Bezier Curve approximation/interpolation 
2. Discrete Differential Geometery: estimate surface curvatures using Jit-fitting/ Discrete Gauss-Bonet
3. Mesh-subdivision.
4. Laplacian Smoothing of mesh.
5. Constrained Quad mesh   Optimization: co-planarity of vertex-star;co-planarity of quad;supporting structure; othorgonal quad constrain.


