#include <igl/read_triangle_mesh.h>
#include <cassert>
#include <catch2/catch.hpp>

#include "src/AdaMesh.hpp"
#include "tetwild/TetWild.h"
#include "wmtk/utils/InsertTriangleUtils.hpp"


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
	auto faces = int_data_convert<std::array<size_t,3>>(F);

    tetwild::Parameters params;
    params.lr = 1 / 15.0;
    params.init(vertices, faces);
    return params;
};


TEST_CASE("insertion", "[adamesh]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXd F;
    std::string input_path = WMT_DATA_DIR "/bunny.off";
    igl::read_triangle_mesh(input_path, V, F);
    wmtk::logger().info("Read Mesh V={}, F={}", V.rows(), F.rows());

    std::vector<Eigen::Vector3d> vertices(V.rows());
    for (int i = 0; i < V.rows(); i++) vertices[i] = V.row(i);

    auto params = param_from_v_f(vertices, F);
    auto envelope = envelope_from_v_f(vertices, F, params.eps);

    wmtk::logger().info("input_surface.params.eps {}", params.eps);

    wmtk::AdaMesh mesh(params, envelope);

	auto faces = int_data_convert<std::array<size_t,3>>(F);
    wmtk::remove_duplicates(vertices, faces, params.diag_l);

    std::vector<size_t> partition_id(vertices.size(), 0); // serial
    mesh.init_from_input_surface(vertices, faces, partition_id);

    REQUIRE(mesh.check_attributes());
}