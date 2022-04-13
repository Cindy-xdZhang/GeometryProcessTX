# CS272: Geometric Modeling
## Requirements and build:
This is a viewer application, using [GLFW](https://www.glfw.org), [Dear ImGui](https://github.com/ocornut/imgui), [Eigen](https://eigen.tuxfamily.org), and [OpenMesh](https://www.openmesh.org).
All dependencies are included in the `external` directory.
You can build this application with `cmake`. It has been tested under **Windows** and **Linux**.

## functionality:
0.Discrete Differential Geometery: estimate surface curvatures using Jit-fitting. 
1. Bezier Curve  &   Approximation of B-spline/interpolation of B-spline.
2. Discrete Differential Geometery: estimate surface curvatures using  Discrete Gauss-Bonet (Normal cycle).
3. Simple remesh.
4. Laplacian Smoothing of mesh.
5. Constrained Quad mesh   Optimization: co-planarity of vertex-star; co-planarity of quad;supporting structure; othorgonal quad constrain.

## Some Sample Result: 
Discrete Gausssian Curvature from jit-fitting![pic](showdemo/1.png).
B-spline interpolation ![pic](showdemo/2.png).
Laplacian smoothing ![pic](showdemo/analysis.png).
Laplacian smoothing ![pic](showdemo/iter50_uniform_0.01.png).
quad mesh optimization ![pic](showdemo/remeshing.png).
quad mesh optimization ![pic](showdemo/combined-ondition_png.png).


