#include <igl/barycenter.h>
#include <igl/read_triangle_mesh.h>
#include <cassert>
#include <catch2/catch.hpp>

#include "AdaMesh.hpp"
#include "spdlog/spdlog.h"
#include "tetwild/TetWild.h"
#include "wmtk/utils/InsertTriangleUtils.hpp"
#include "wmtk/utils/Logger.hpp"


template <typename T, typename EigenT>
auto int_data_convert(const EigenT& F)
{
    std::vector<T> faces(F.rows());
    for (int i = 0; i < F.rows(); i++)
        for (int j = 0; j < 3; j++) faces[i][j] = F(i, j);
    return faces;
}

auto envelope_from_v_f = [](auto& vertices, auto& F, auto eps) {
    auto env_faces = int_data_convert<fastEnvelope::Vector3i>(F);

    fastEnvelope::FastEnvelope envelope;

    envelope.init(vertices, env_faces, eps);
    return envelope;
};

auto param_from_v_f = [](auto& vertices, auto& F) {
    auto faces = int_data_convert<std::array<size_t, 3>>(F);

    tetwild::Parameters params;
    params.lr = 1 / 15.0;
    params.init(vertices, faces);
    return params;
};

TEST_CASE("insert-exist-points", "[adamesh]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXd F;
    std::string input_path = WMT_DATA_DIR "/bunny.off";
    igl::read_triangle_mesh(input_path, V, F);

    std::vector<Eigen::Vector3d> vertices(V.rows());
    for (int i = 0; i < V.rows(); i++) vertices[i] = V.row(i);

    auto params = param_from_v_f(vertices, F);
    auto envelope = envelope_from_v_f(vertices, F, params.eps);
    wmtk::AdaMesh mesh(params, envelope);

    auto faces = int_data_convert<std::array<size_t, 3>>(F);
    wmtk::remove_duplicates(vertices, faces, params.diag_l);

    std::vector<size_t> partition_id(vertices.size(), 0); // serial
    mesh.init_from_input_surface(vertices, faces, partition_id);
    std::vector<int> vids;
    mesh.insert_all_points(vertices, vids);
    REQUIRE([&]() -> bool {
        for (auto i = 0; i < vids.size(); i++)
            if (i != vids[i]) return false;
        return true;
    }());

    Eigen::MatrixXd F_BC;
    igl::barycenter(V, F, F_BC);
    std::vector<Eigen::Vector3d> face_bc(F_BC.rows());
    for (auto i = 0; i < F_BC.rows(); i++) face_bc[i] = F_BC.row(i);

    wmtk::logger().enable_backtrace(100);
    std::vector<int> bc_vids;
    mesh.insert_all_points(face_bc, bc_vids);
    REQUIRE([&]() -> bool { // appending all faces
        auto num = vertices.size();
        for (auto i = 0; i < bc_vids.size(); i++)
            if (mesh.m_vertex_attribute[bc_vids[i]].m_posf != face_bc[i]) return false;
        return true;
    }());
}


TEST_CASE("insert-check-tag", "[adamesh]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXd F;
    std::string input_path = WMT_DATA_DIR "/bunny.off";
    igl::read_triangle_mesh(input_path, V, F);

    std::vector<Eigen::Vector3d> vertices(V.rows());
    for (int i = 0; i < V.rows(); i++) vertices[i] = V.row(i);

    auto params = param_from_v_f(vertices, F);
    auto envelope = envelope_from_v_f(vertices, F, params.eps);
    wmtk::AdaMesh mesh(params, envelope);

    auto faces = int_data_convert<std::array<size_t, 3>>(F);
    wmtk::remove_duplicates(vertices, faces, params.diag_l);

    std::vector<size_t> partition_id(vertices.size(), 0); // serial
    mesh.init_from_input_surface(vertices, faces, partition_id);

    // output surface
    mesh.output_faces("initsurf.off", [](auto& at) { return at.m_is_surface_fs; });
    mesh.output_faces("initbbox.off", [](auto& at) { return at.m_is_bbox_fs >= 0; });
    std::vector<int> vids;
    mesh.insert_all_points(vertices, vids);
    REQUIRE([&]() -> bool {
        for (auto i = 0; i < vids.size(); i++)
            if (i != vids[i]) return false;
        return true;
    }());

    Eigen::MatrixXd F_BC;
    igl::barycenter(V, F, F_BC);
    std::vector<Eigen::Vector3d> face_bc(F_BC.rows());
    for (auto i = 0; i < F_BC.rows(); i++) face_bc[i] = F_BC.row(i);

    wmtk::logger().enable_backtrace(100);
    std::vector<int> bc_vids;
    mesh.insert_all_points(face_bc, bc_vids);
    REQUIRE([&]() -> bool { // appending all faces
        auto num = vertices.size();
        for (auto i = 0; i < bc_vids.size(); i++)
            if (mesh.m_vertex_attribute[bc_vids[i]].m_posf != face_bc[i]) return false;
        return true;
    }());
    mesh.output_faces("pointsurf.off", [](auto& at) { return at.m_is_surface_fs; });
    mesh.output_faces("pointbbox.off", [](auto& at) { return at.m_is_bbox_fs >= 0; });
}


TEST_CASE("offset-bunny", "[adamesh]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXd F;
    std::string input_path = WMT_DATA_DIR "/bunny.off";
    igl::read_triangle_mesh(input_path, V, F);

    std::vector<Eigen::Vector3d> vertices(V.rows());
    for (int i = 0; i < V.rows(); i++) vertices[i] = V.row(i);

    auto params = param_from_v_f(vertices, F);
    auto envelope = envelope_from_v_f(vertices, F, params.eps);
    wmtk::AdaMesh mesh(params, envelope);

    auto faces = int_data_convert<std::array<size_t, 3>>(F);
    wmtk::remove_duplicates(vertices, faces, params.diag_l);
    
    {
        mesh.init_from_delaunay_box_mesh(vertices, 1.0);
        // std::vector<size_t> partition_id(vertices.size(), 0); // serial
        // mesh.insert_triangles_to_mesh(faces, partition_id);
        // mesh.finalize_triangle_insertion(faces);
    }

    auto shifted_points = vertices;
    for (auto& s : shifted_points) s[0] += 0.1;

    std::vector<int> new_vid;
    mesh.insert_all_points(shifted_points, new_vid);

    auto mapped_faces = faces;
    for (auto& f: mapped_faces) {
        for (auto& v: f) v = new_vid[v];
    }
    std::vector<size_t> partition_id(mesh.vert_capacity(), 0); // serial
    mesh.insert_triangles_to_mesh(mapped_faces, partition_id);
    mesh.finalize_triangle_insertion(mapped_faces);
    mesh.output_faces("pointsurf.off", [](auto& at) { return at.m_is_surface_fs; });
    mesh.output_faces("pointbbox.off", [](auto& at) { return at.m_is_bbox_fs >= 0; });
}