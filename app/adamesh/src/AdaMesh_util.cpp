#include <array>
#include "AdaMesh.hpp"

void wmtk::AdaMesh::serialize(std::string filename)
{
    // topology, without consolidate
    Eigen::MatrixXi T = -Eigen::MatrixXi::Ones(tet_capacity(), 4);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(vert_capacity(), 3);
    for_each_tetra([&](auto& t) {
        auto ti = (t.tid(*this));
        auto vs = oriented_tet_vids(t);
        for (auto j = 0; j < 4; j++) T(ti, j) = vs[j];
    });

    for_each_vertex([&](auto& v) {
        auto vi = v.vid(*this);
        V.row(vi) = m_vertex_attribute[vi].m_posf;
    });


    // tags, Face, is_surface_fs, is_bbox_fs
    
    Eigen::MatrixXi FaceTag = -Eigen::MatrixXi::Ones(tet_capacity()*6, 2);
    for (auto i=0; i<m_face_attribute.size(); i++) {
      auto& att = m_face_attribute[i];
      FaceTag.row(i) << att.m_is_surface_fs, att.m_is_bbox_fs ;
    }
}


void wmtk::AdaMesh::deserialize(std::string filename)
{
    // topology, without consolidate
    Eigen::MatrixXi T;
    Eigen::MatrixXd V;
    Eigen::MatrixXi FaceTag;

    std::vector<std::array<size_t, 4>> tets(T.rows());
    for (auto i=0 ;i<T.rows(); i++) {
      tets.
    }
    init(V.rows(), tets);


    // for_each_tetra([&](auto& t) {
    //     auto ti = (t.tid(*this));
    //     auto vs = oriented_tet_vids(t);
    //     for (auto j = 0; j < 4; j++) T(ti, j) = vs[j];
    // });

    // for_each_vertex([&](auto& v) {
    //     auto vi = v.vid(*this);
    //     V.row(vi) = m_vertex_attribute[vi].m_posf;
    // });


    // tags, Face, is_surface_fs, is_bbox_fs
    
    // for (auto i=0; i<m_face_attribute.size(); i++) {
    //   auto& att = m_face_attribute[i];
    //   FaceTag.row(i) << att.m_is_surface_fs, att.m_is_bbox_fs ;
    // }
}