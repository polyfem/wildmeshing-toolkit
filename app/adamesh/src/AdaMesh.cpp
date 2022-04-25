#include "AdaMesh.hpp"

#include <wmtk/TetMeshOperations.hpp>
#include "prism/geogram/AABB.hpp"
#include "prism/geogram/AABB_tet.hpp"
#include "prism/predicates/inside_prism_tetra.hpp"
#include "tetwild/TetWild.h"
#include "wmtk/ConcurrentTetMesh.h"

#include <geogram/numerics/predicates.h>

auto replace = [](auto& arr, size_t a, size_t b) {
    for (auto i = 0; i < arr.size(); i++) {
        if (arr[i] == a) {
            arr[i] = b;
            return i;
        }
    }
    assert(false);
    return -1;
};


auto degenerate_config =
    [](const auto& tetv, const auto& tet, const Eigen::Vector3d& pt) -> std::array<int, 3> {
    using GEO::PCK::orient_3d;
    using GEO::PCK::points_are_colinear_3d;
    std::array<bool, 4> colinear{false, false, false, false};
    for (auto i = 0; i < 4; i++) {
        if (pt == tetv[tet[i]].m_posf) return {int(tet[i]), -1, -1}; // vertex
    }
    for (auto i = 0; i < 4; i++) {
        if (orient_3d(
                pt.data(),
                tetv[tet[(i + 1) % 4]].m_posf.data(),
                tetv[tet[(i + 2) % 4]].m_posf.data(),
                tetv[tet[(i + 3) % 4]].m_posf.data()) == 0) {
            for (auto j = 0; j < 3; j++) {
                if (points_are_colinear_3d(
                        pt.data(),
                        tetv[tet[(i + 1 + j) % 4]].m_posf.data(),
                        tetv[tet[(i + 1 + (j + 1) % 3) % 4]].m_posf.data()))
                    return {
                        (int)tet[(i + 1 + j) % 4],
                        (int)tet[(i + 1 + (j + 1) % 3) % 4],
                        -1}; // edge
            }
            return {(int)tet[(i + 1) % 4], (int)tet[(i + 2) % 4], (int)tet[(i + 3) % 4]}; // face
        }
    }
    return {-1, -1, -1}; // general
};


namespace wmtk {
AdaMesh::AdaMesh(
    tetwild::Parameters& params,
    fastEnvelope::FastEnvelope& envelope,
    const RowMatd& V,
    const RowMati& T)
    : tetwild::TetWild(params, envelope)
{
    // TODO: manually set tags.
    assert(false);
    assert(false);
    m_vertex_attribute.resize(V.rows());
    tet_attrs.resize(T.rows());

    std::vector<std::array<size_t, 4>> tets(T.rows());
    for (auto i = 0; i < tets.size(); i++) {
        for (auto j = 0; j < 4; j++) tets[i][j] = (size_t)T(i, j);
    }
    init(V.rows(), tets);

    // attrs
    for (auto i = 0; i < V.rows(); i++) {
        m_vertex_attribute[i] = tetwild::VertexAttributes(Vector3d(V(i, 0), V(i, 1), V(i, 2)));
    }
}

AdaMesh::AdaMesh(tetwild::Parameters& params,
    fastEnvelope::FastEnvelope& envelope):
tetwild::TetWild(params, envelope) {
}


struct DivideTet : public TetMesh::OperationBuilder
{
    const AdaMesh& m;
    TetMesh::Tuple tet;
    size_t ux;

    DivideTet(const AdaMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        tet = t;
        return {t.tid(m)};
    }
    int request_vert_slots() { return 1; }
    std::vector<std::array<size_t, 4>> replacing_tets(const std::vector<size_t>& slots)
    {
        assert(slots.size() == 1);
        ux = slots.front();

        std::array<size_t, 4> t_conn;
        auto vs = m.oriented_tet_vertices(tet);
        for (auto i = 0; i < 4; i++) t_conn[i] = vs[i].vid(m);
        auto new_tets = std::vector<std::array<size_t, 4>>(4, t_conn);
        for (auto i = 0; i < 4; i++) new_tets[i][i] = ux;

        return new_tets;
    }
};


struct SplitFace : public TetMesh::OperationBuilder
{
    const AdaMesh& m;
    std::vector<size_t> affected;
    std::array<size_t, 3> tri_verts;
    size_t ux;

    SplitFace(const AdaMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        auto oppo = t.switch_face(m);
        affected = {t.tid(m), oppo.tid(m)};

        return affected;
    }
    int request_vert_slots() { return 1; }
    std::vector<std::array<size_t, 4>> replacing_tets(const std::vector<size_t>& slots)
    {
        assert(slots.size() == 1);
        ux = slots.front();

        auto new_tets = std::vector<std::array<size_t, 4>>();

        new_tets.reserve(2 * 3);
        for (auto i = 0; i < 2; i++) {
            auto t_conn = m.oriented_tet_vids(m.tuple_from_tet(affected[i]));
            for (auto j = 0; j < 3; j++) {
                new_tets.push_back(t_conn);
                replace(new_tets.back(), tri_verts[j], ux);
            }
        }
        return new_tets;
    }
};

void AdaMesh::insert_all_points(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<int>& hint_tid)
{
    if (hint_tid.empty()) { // find hint
        auto tree = prism::geogram::AABB_tet(tetV, tetT);
    }
    std::function<int(size_t, const Vector3d&)> find_containing_tet;
    std::map<int, std::set<int>> split_maps;
    auto& m = *this;
    find_containing_tet = [&m, &split_maps, &tetv = m_vertex_attribute, &find_containing_tet](
                              size_t tid,
                              const Vector3d& pt) -> int {
        auto it = split_maps.find(tid);
        if (it == split_maps.end()) { // leaf
            auto vs = m.oriented_tet_vertices(m.tuple_from_tet(tid));
            if (::prism::predicates::point_in_tetrahedron(
                    pt,
                    tetv[vs[0].vid(m)].m_posf,
                    tetv[vs[1].vid(m)].m_posf,
                    tetv[vs[2].vid(m)].m_posf,
                    tetv[vs[3].vid(m)].m_posf))
                return tid;
        } else {
            for (auto v : it->second) {
                auto res = find_containing_tet(v, pt);
                if (res != -1) return res;
            }
        }
        return -1;
    };


    std::vector<int> new_vid(points.size());
    for (auto i = 0; i < points.size(); i++) {
        auto pt = points[i];
        auto tid = find_containing_tet(hint_tid[i], pt); // the final tid

        auto config =
            degenerate_config(m_vertex_attribute, m.oriented_tet_vids(m.tuple_from_tet(tid)), pt);
        wmtk::logger().trace("insert {} with config {}", i, config);

        std::vector<Tuple> new_tets;
        if (config[0] != -1) {
            if (config[1] == -1) { // point degenerate
                // vertex
                new_vid[i] = config[0];
                continue;
            } else if (config[2] == -1) {
                // edge
                auto tup = m.tuple_from_edge({(size_t)config[0], (size_t)config[1]});
                auto spl_edge = wmtk::SplitEdge(m);
                auto suc = m.customized_operation(spl_edge, tup, new_tets);

                for (auto j = 0; j < spl_edge.affected.size();
                     j++) { // this follows from the convention inside
                    auto& s = split_maps[spl_edge.affected[j]];
                    for (auto k = 0; k < 2; k++) s.insert(new_tets[j * 2 + k].tid(m));
                }
                new_vid[i] = spl_edge.ux;
            } else {
                // face
                SplitFace spl_face(m);
                spl_face.tri_verts = {(size_t)config[0], (size_t)config[1], (size_t)config[2]};
                auto [tup, fid] = m.tuple_from_face(spl_face.tri_verts);
                auto suc = m.customized_operation(spl_face, tup, new_tets);

                assert(suc);
                for (auto j = 0; j < 2; j++) { // this follows from the convention inside splitface
                    auto& s = split_maps[spl_face.affected[j]];
                    for (auto k = 0; k < 3; k++) s.insert(new_tets[j * 3 + k].tid(m));
                }
                new_vid[i] = spl_face.ux;
            }
        } else {
            // general position, insert the single point
            DivideTet divide_tet(m);
            auto tup = m.tuple_from_tet(tid);

            auto suc = m.customized_operation(divide_tet, tup, new_tets);

            assert(suc);
            auto& s = split_maps[tid];
            for (auto t : new_tets) {
                s.insert(t.tid(m));
            }
            new_vid[i] = divide_tet.ux;
        }
        m.m_vertex_attribute[new_vid[i]] = tetwild::VertexAttributes(pt);
    }
}
} // namespace wmtk