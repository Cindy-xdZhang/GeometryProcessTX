#include "utils/OpenMesh.h"


void OpenMesh::toRenderData(
        OpenMesh::Mesh mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F,
        Eigen::MatrixXi& FtoF, Eigen::MatrixXd& N, Eigen::MatrixXd& UV,
        Eigen::MatrixXd& P1, Eigen::MatrixXd& P2)
{
    // Triangulate
    std::map<int, int> face_map;
    std::map<int, OpenMesh::Mesh::Normal> normals_orig;
    OpenMesh::Mesh triMesh(mesh);
    for (OpenMesh::Mesh::FaceIter f_it = triMesh.faces_begin(); f_it != triMesh.faces_end(); ++f_it)
    {
        OpenMesh::Mesh::FaceHandle fh = *f_it;
        if (face_map.find(fh.idx()) == face_map.end())
        {
            face_map[fh.idx()] = fh.idx();
        }

        OpenMesh::Mesh::Normal n = triMesh.normal(fh);
        if (normals_orig.find(fh.idx()) == normals_orig.end())
        {
            normals_orig[fh.idx()] = n;
        }

        OpenMesh::Mesh::HalfedgeHandle base_heh(triMesh.halfedge_handle(fh));
        OpenMesh::Mesh::VertexHandle start_vh = triMesh.from_vertex_handle(base_heh);
        OpenMesh::Mesh::HalfedgeHandle prev_heh(triMesh.prev_halfedge_handle(base_heh));
        OpenMesh::Mesh::HalfedgeHandle next_heh(triMesh.next_halfedge_handle(base_heh));

        while (triMesh.to_vertex_handle(triMesh.next_halfedge_handle(next_heh)) != start_vh)
        {
            OpenMesh::Mesh::HalfedgeHandle next_next_heh(triMesh.next_halfedge_handle(next_heh));

            OpenMesh::Mesh::FaceHandle new_fh = triMesh.new_face();
            triMesh.set_halfedge_handle(new_fh, base_heh);

            face_map[new_fh.idx()] = fh.idx();

            normals_orig[new_fh.idx()] = n;

            OpenMesh::Mesh::HalfedgeHandle new_heh = triMesh.new_edge(triMesh.to_vertex_handle(next_heh), start_vh);

            triMesh.set_next_halfedge_handle(base_heh, next_heh);
            triMesh.set_next_halfedge_handle(next_heh, new_heh);
            triMesh.set_next_halfedge_handle(new_heh, base_heh);

            triMesh.set_face_handle(base_heh, new_fh);
            triMesh.set_face_handle(next_heh, new_fh);
            triMesh.set_face_handle(new_heh, new_fh);

            triMesh.copy_all_properties(prev_heh, new_heh, true);
            triMesh.copy_all_properties(prev_heh, triMesh.opposite_halfedge_handle(new_heh), true);
            triMesh.copy_all_properties(fh, new_fh, true);

            base_heh = triMesh.opposite_halfedge_handle(new_heh);
            next_heh = next_next_heh;
        }

        triMesh.set_halfedge_handle(fh, base_heh);  //the last face takes the handle _fh

        triMesh.set_next_halfedge_handle(base_heh, next_heh);
        triMesh.set_next_halfedge_handle(triMesh.next_halfedge_handle(next_heh), base_heh);

        triMesh.set_face_handle(base_heh, fh);
    }

    // Resize arrays
    V.resize(triMesh.n_vertices(), 3);
    F.resize(triMesh.n_faces(), 3);
    FtoF.resize(triMesh.n_faces(), 1);
    N.resize(triMesh.n_faces(), 3);
    if (mesh.has_vertex_texcoords2D())
    {
        UV.resize(triMesh.n_vertices(), 2);
    }

    // Vertices
    for (OpenMesh::Mesh::VertexIter v_it = triMesh.vertices_begin();
         v_it != triMesh.vertices_end(); v_it++)
    {
        OpenMesh::Mesh::Point p = triMesh.point(*v_it);
        V(v_it->idx(), 0) = p[0];
        V(v_it->idx(), 1) = p[1];
        V(v_it->idx(), 2) = p[2];
    }

    // Faces
    for (OpenMesh::Mesh::FaceIter f_it = triMesh.faces_begin();
         f_it != triMesh.faces_end(); ++f_it)
    {
        int vi = 0;
        for (OpenMesh::Mesh::ConstFaceVertexCCWIter fvi = triMesh.cfv_ccwbegin(*f_it);
             vi < 3 && fvi != triMesh.cfv_ccwend(*f_it); ++fvi)
        {
            F(f_it->idx(), vi) = fvi->idx();
            vi++;
        }
    }

    // Face map
    for (auto& it : face_map)
    {
        FtoF(it.first, 0) = it.second;
    }

    // Normals
    for (auto& it : normals_orig)
    {
        OpenMesh::Mesh::Normal n = it.second;
        N(it.first, 0) = n[0];
        N(it.first, 1) = n[1];
        N(it.first, 2) = n[2];
    }

    // TexCoords
    if (mesh.has_vertex_texcoords2D())
    {
        for (OpenMesh::Mesh::VertexIter v_it = triMesh.vertices_begin();
             v_it != triMesh.vertices_end(); v_it++)
        {
            OpenMesh::Mesh::TexCoord2D tex = triMesh.texcoord2D(*v_it);
            UV(v_it->idx(), 0) = tex[0];
            UV(v_it->idx(), 1) = tex[1];
        }
    }

    // Edges
    P1.resize(mesh.n_edges(), 3);
    P2.resize(mesh.n_edges(), 3);
    for (OpenMesh::Mesh::EdgeIter e_it = mesh.edges_begin();
         e_it != mesh.edges_end(); ++e_it)
    {
        OpenMesh::SmartVertexHandle vh1 = mesh.halfedge_handle(*e_it, 0).to();
        OpenMesh::SmartVertexHandle vh2 = mesh.halfedge_handle(*e_it, 0).from();
        OpenMesh::Mesh::Point v1 = mesh.point(mesh.vertex_handle(vh1.idx()));
        OpenMesh::Mesh::Point v2 = mesh.point(mesh.vertex_handle(vh2.idx()));

        P1(e_it->idx(), 0) = v1[0];
        P1(e_it->idx(), 1) = v1[1];
        P1(e_it->idx(), 2) = v1[2];
        P2(e_it->idx(), 0) = v2[0];
        P2(e_it->idx(), 1) = v2[1];
        P2(e_it->idx(), 2) = v2[2];
    }
}
