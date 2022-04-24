#pragma once

#include <tetwild/TetWild.h>

#include <tetwild/Rational.h>
#include <Eigen/Core>
#include <stdexcept>

namespace wmtk {
using RowMatd = Eigen::Matrix<double, -1,-1,Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1,-1,Eigen::RowMajor>;
using Vec3d = Eigen::RowVector3d;

using Vector3r = Eigen::Matrix<apps::Rational, 3, 1>;
using Vector2r = Eigen::Matrix<apps::Rational, 2, 1>;
struct AdaMesh : public tetwild::TetWild
{
public:
    AdaMesh(const RowMatd& V, const RowMati& T); // initialize topology and attrs
    void insert_all_points(
        const std::vector<Vec3d>& points,
        const std::vector<int>& hint_tid); // insert points
    void insert_all_triangles(const std::vector<std::array<size_t, 3>>& tris);

public: // callbacks
    struct TriangleInsertionLocalInfoCache
    {
        // local info: for each face insertion
        int face_id;
        std::vector<std::array<size_t, 3>> old_face_vids;
    };
    TriangleInsertionLocalInfoCache triangle_insertion_local_cache;
    std::map<std::array<size_t, 3>, std::vector<int>> tet_face_tags;
    bool triangle_insertion_before(const std::vector<Tuple>& tup) override;
    bool triangle_insertion_after(const std::vector<std::vector<Tuple>>& new_faces) override;
    void finalize_triangle_insertion(const std::vector<std::array<size_t, 3>>& tris);

};
} // namespace wmtk::prism