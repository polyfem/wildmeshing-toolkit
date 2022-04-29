#include "AdaMesh.hpp"

#include <array>
#include <cassert>
#include <vector>
#include <wmtk/TetMeshOperations.hpp>
#include "prism/geogram/AABB.hpp"
#include "prism/geogram/AABB_tet.hpp"
#include "prism/predicates/inside_prism_tetra.hpp"
#include "spdlog/spdlog.h"
#include "tetwild/TetWild.h"
#include "wmtk/ConcurrentTetMesh.h"
#include "wmtk/TetMesh.h"
#include "wmtk/utils/InsertTriangleUtils.hpp"
#include "wmtk/utils/Logger.hpp"

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
struct SplitFace : public TetMesh::OperationBuilder
{
    const TetMesh& m;
    std::vector<size_t> affected;
    std::array<size_t, 3> tri_verts;
    size_t ux;

    SplitFace(const TetMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        auto oppo = m.switch_tetrahedron(t).value();
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
} // namespace wmtk

namespace adamesh {
struct SplitFace : public wmtk::SplitFace
{
    struct
    {
        Eigen::Vector3d pos;
        tetwild::FaceAttributes face_tag;
        std::map<std::array<size_t, 3>, tetwild::FaceAttributes> f_attrs;
    } cache;
    wmtk::AdaMesh& m_app;
    SplitFace(wmtk::AdaMesh& m_)
        : wmtk::SplitFace(m_)
        , m_app(m_)
    {}

    bool before(const wmtk::TetMesh::Tuple&);
    bool after(const std::vector<wmtk::TetMesh::Tuple>&);
};

struct DivideTet : public wmtk::DivideTet
{
    struct
    {
        Eigen::Vector3d pos;
        std::vector<std::pair<std::array<size_t, 3>, tetwild::FaceAttributes>> f_attrs;
        std::vector<int> old_fids;
    } cache;

    wmtk::AdaMesh& m_app;
    DivideTet(wmtk::AdaMesh& m_)
        : wmtk::DivideTet(m_)
        , m_app(m_)
    {}
    bool before(const wmtk::TetMesh::Tuple&);
    bool after(const std::vector<wmtk::TetMesh::Tuple>&);
};
} // namespace adamesh

bool adamesh::SplitFace::before(const wmtk::TetMesh::Tuple& loc)
{
    cache.face_tag = m_app.m_face_attribute[loc.fid(m)];

    for (auto& t : {loc, loc.switch_tetrahedron(m).value()}) {
        auto vs = m.oriented_tet_vids(t);
        for (int j = 0; j < 4; j++) {
            std::array<size_t, 3> f_vids = {{vs[(j + 1) % 4], vs[(j + 2) % 4], vs[(j + 3) % 4]}};
            std::sort(f_vids.begin(), f_vids.end());
            auto [_, global_fid] = m.tuple_from_face(f_vids);
            cache.f_attrs.emplace(f_vids, m_app.m_face_attribute[global_fid]);
        }
    }
    return true;
}

bool adamesh::SplitFace::after(const std::vector<wmtk::TetMesh::Tuple>& locs)
{
    // vertex
    m_app.m_vertex_attribute[ux] = tetwild::VertexAttributes(cache.pos);
    if (cache.face_tag.m_is_surface_fs) m_app.m_vertex_attribute[ux].m_is_on_surface = true;
    if (cache.face_tag.m_is_bbox_fs >= 0)
        m_app.m_vertex_attribute[ux].on_bbox_faces.push_back(cache.face_tag.m_is_bbox_fs);

    // face: 2 tet -> 6 tet, bounding faces directly inherit,
    for (auto& [old_vids, attr] : cache.f_attrs) {
        std::vector<int> j_vn; // vertex apart from  tri
        for (int j = 0; j < 3; j++) {
            if (old_vids[j] != tri_verts[0] && old_vids[j] != tri_verts[1] &&
                old_vids[j] != tri_verts[2]) {
                j_vn.push_back(j);
            }
        }

        if (j_vn.size() == 0) { // tri 0-1-2 of interest, split to 3.
            for (auto j = 0; j < 3; j++) {
                auto vids = old_vids;
                vids[j] = ux;
                auto [_, global_fid] = m.tuple_from_face(vids);
                m_app.m_face_attribute[global_fid] = attr;
            }
        } else if (j_vn.size() == 1) { // tri 0-1-T/B
            auto [_, global_fid] = m.tuple_from_face(old_vids); // inherit boundary
            m_app.m_face_attribute[global_fid] = attr;

            for (auto j = 0; j < 3; j++) { // internal (0, x, T), duplicate appear here.
                auto [_, global_fid] =
                    m.tuple_from_face({{ux, tri_verts[j], old_vids[j_vn.front()]}});
                m_app.m_face_attribute[global_fid].reset();
            }
        } else { // j_vn.size() == 2
            assert(false);
        }
    }


    // tet
    for (auto& loc : locs) {
        m_app.m_tet_attribute[loc.tid(m)].m_quality = m_app.get_quality(loc);
    }

    return true;
}


bool adamesh::DivideTet::before(const wmtk::TetMesh::Tuple& t)
{
    auto vs = m.oriented_tet_vids(t);
    for (int j = 0; j < 4; j++) {
        std::array<size_t, 3> f_vids = {{vs[(j + 1) % 4], vs[(j + 2) % 4], vs[(j + 3) % 4]}};
        std::sort(f_vids.begin(), f_vids.end());
        auto [_, global_fid] = m.tuple_from_face(f_vids);
        cache.old_fids.push_back(global_fid);
        cache.f_attrs.emplace_back(f_vids, m_app.m_face_attribute[global_fid]);
        cache.old_fids.push_back(global_fid);
    }
    return true;
}

bool adamesh::DivideTet::after(const std::vector<wmtk::TetMesh::Tuple>& locs)
{
    // vertex
    m_app.m_vertex_attribute[ux] = tetwild::VertexAttributes(cache.pos);

    // face: 2 tet -> 6 tet, bounding faces directly inherit,
    for (auto f : cache.old_fids) {
        m_app.m_face_attribute[f].reset();
    }
    for (auto& [old_vids, attr] : cache.f_attrs) {
        auto [_, global_fid] = m.tuple_from_face(old_vids);
        m_app.m_face_attribute[global_fid] = attr;
    }

    // tet
    for (auto& loc : locs) {
        m_app.m_tet_attribute[loc.tid(m)].m_quality = m_app.get_quality(loc);
    }

    return true;
}


namespace wmtk {
// AdaMesh::AdaMesh(
//     tetwild::Parameters& params,
//     fastEnvelope::FastEnvelope& envelope,
//     const RowMatd& V,
//     const RowMati& T)
//     : tetwild::TetWild(params, envelope)
// {
//     // TODO: manually set tags.
//     assert(false);
//     assert(false);
//     m_vertex_attribute.resize(V.rows());
//     tet_attrs.resize(T.rows());

//     std::vector<std::array<size_t, 4>> tets(T.rows());
//     for (auto i = 0; i < tets.size(); i++) {
//         for (auto j = 0; j < 4; j++) tets[i][j] = (size_t)T(i, j);
//     }
//     init(V.rows(), tets);

//     // attrs
//     for (auto i = 0; i < V.rows(); i++) {
//         m_vertex_attribute[i] = tetwild::VertexAttributes(Vector3d(V(i, 0), V(i, 1), V(i, 2)));
//     }
// }

auto aabb_tree_tid(AdaMesh& m, const std::vector<Eigen::Vector3d>& points)
{
    using RowMatd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
    using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;

    RowMatd V = RowMatd::Zero(m.vert_capacity(), 3);
    m.for_each_vertex([&](auto& v) {
        auto i = v.vid(m);
        assert(m.m_vertex_attribute[i].m_is_rounded);
        V.row(i) = m.m_vertex_attribute[i].m_posf;
    });

    std::vector<std::array<size_t, 4>> tets;
    static_cast<wmtk::TetMesh>(m).for_each_tetra(
        [&](auto& t) { tets.push_back(m.oriented_tet_vids(t)); });


    RowMati T = RowMati(tets.size(), 4);
    for (auto i = 0; i < tets.size(); i++)
        for (auto j = 0; j < 4; j++) T(i, j) = tets[i][j];

    auto tree = prism::geogram::AABB_tet(V, T);

    std::vector<int> hint_tids;
    for (auto p : points) {
        auto [tid, bc] = tree.point_query(p);
        if (tid < 0) {
            wmtk::logger().critical("Outside, need expansion.");
            assert(false);
        }
        hint_tids.push_back(tid);
    }
    return hint_tids;
}

void AdaMesh::insert_all_points(
    const std::vector<Eigen::Vector3d>& points,
    std::vector<int>& new_vid)
{
    new_vid.resize(points.size());
    auto& m = *this;

    auto hint_tid = aabb_tree_tid(*this, points);

    using tet_ts = std::pair<size_t, size_t>;
    std::map<tet_ts, std::vector<tet_ts>> split_maps;
    std::vector<size_t> split_map_hash(m.tet_capacity(), 0);

    std::function<int(size_t, size_t, const Vector3d&)> find_containing_tet;
    find_containing_tet = [&m,
                           &split_maps,
                           &split_map_hash,
                           &tetv = m_vertex_attribute,
                           &find_containing_tet](size_t tid, size_t ts, const Vector3d& pt) -> int {
        auto it = split_maps.find({tid, ts});
        if (it == split_maps.end()) { // leaf
            assert(ts == split_map_hash[tid]); // up to date.
            auto vs = m.oriented_tet_vids(m.tuple_from_tet(tid));
            wmtk::logger().trace("Checking tri {} with {}", tid, ts);
            for (auto i = 0; i < 4; i++) wmtk::logger().trace("{}", tetv[vs[i]].m_posf.transpose());
            if (::prism::predicates::point_in_tetrahedron(
                    pt,
                    tetv[vs[0]].m_posf,
                    tetv[vs[1]].m_posf,
                    tetv[vs[2]].m_posf,
                    tetv[vs[3]].m_posf))
                return tid;
        } else {
            for (auto [tid, ts] : it->second) {
                auto res = find_containing_tet(tid, ts, pt);
                if (res != -1) return res;
            }
            // here is reachable, when a same-id tet is split twice.
        }
        return -1;
    };

    auto record_split = [&split_map_hash, &split_maps](
                            const std::vector<size_t>& remove,
                            const std::vector<std::vector<size_t>>& new_list) {
        for (auto i = 0; i < remove.size(); i++) {
            split_map_hash[remove[i]]++;
        } // the order of increment vs emplace is important here, similar to aliasling.
        for (auto i = 0; i < remove.size(); i++) {
            auto ts = split_map_hash[remove[i]] - 1; // just split
            auto& s = split_maps[{remove[i], ts}];
            for (auto tt : new_list[i]) {
                if (split_map_hash.size() < tt + 1) split_map_hash.resize(tt + 1, 0);
                s.emplace_back(tt, split_map_hash[tt]);
            }
        }

        // wmtk::logger().critical(split_maps);
    };

    ///
    /// Duplicate points: nothing, keep vid
    /// Split Edge: Tet Quality same.
    ////   Face Perserve taggings same as tetwild::edge-splits.
    ///    Vertex Same.
    /// Split Face: Tet Quality easy. Face tagging: inherit directly. Vertex Tagging, only if it is face
    /// Divide Tet: Face Tagging: all internal, but keep previously. Vertex Tagging, no new tag.
    ///
    for (auto i = 0; i < points.size(); i++) {
        auto pt = points[i];
        wmtk::logger().trace("point {} with hint {}", pt.transpose(), hint_tid[i]);
        auto tid = find_containing_tet(hint_tid[i], 0, pt);
        assert(tid != -1);

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
                wmtk::logger().trace("Splitting Edge... {}", config);
                auto tup = m.tuple_from_edge({(size_t)config[0], (size_t)config[1]});
                auto spl_edge = tetwild::SplitEdge(m);

                auto suc = m.customized_operation(spl_edge, tup, new_tets);
                m_vertex_attribute[spl_edge.ux] = tetwild::VertexAttributes(pt);

                std::vector<std::vector<size_t>> new_list;
                for (auto j = 0; j < spl_edge.affected.size();
                     j++) { // this follows from the convention inside
                    new_list.emplace_back(
                        std::vector<size_t>{new_tets[j * 2].tid(m), new_tets[j * 2 + 1].tid(m)});
                }
                record_split(spl_edge.affected, new_list);

                new_vid[i] = spl_edge.ux;
            } else {
                // face
                wmtk::logger().trace("Splitting Face... {}", config);
                adamesh::SplitFace spl_face(m);
                spl_face.tri_verts = {(size_t)config[0], (size_t)config[1], (size_t)config[2]};
                spl_face.cache.pos = pt;
                auto [tup, fid] = m.tuple_from_face(spl_face.tri_verts);
                auto suc = m.customized_operation(spl_face, tup, new_tets);

                assert(suc);
                std::vector<std::vector<size_t>> new_list;
                for (auto j = 0; j < 2; j++) { // this follows from the convention inside splitface
                    new_list.emplace_back(std::vector<size_t>{
                        new_tets[j * 2 + 0].tid(m),
                        new_tets[j * 2 + 1].tid(m),
                        new_tets[j * 2 + 2].tid(m)});
                }
                record_split(spl_face.affected, new_list);
                new_vid[i] = spl_face.ux;
            }
        } else {
            // general position, insert the single point
            wmtk::logger().trace("Divide Tet... {}", config);
            adamesh::DivideTet divide_tet(m);
            divide_tet.cache.pos = pt;
            auto tup = m.tuple_from_tet(tid);

            auto suc = m.customized_operation(divide_tet, tup, new_tets);
            if (!suc) wmtk::logger().dump_backtrace();
            assert(suc);
            {
                std::vector<std::vector<size_t>> sub_tids(1);
                for (const auto& t : new_tets) {
                    sub_tids.front().push_back(t.tid(m));
                }
                wmtk::logger().trace("recording (tid {})::{}", tid, sub_tids.front());
                record_split({(size_t)tid}, sub_tids);
            }

            new_vid[i] = divide_tet.ux;
        }
        m.m_vertex_attribute[new_vid[i]] = tetwild::VertexAttributes(pt);
    }
}


auto internal_insert_single_triangle(
    wmtk::TetMesh& m,
    tetwild::TetWild::VertAttCol& m_vertex_attribute,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::array<size_t, 3>& face,
    std::vector<std::array<size_t, 3>>& marked_tet_faces,
    const std::function<bool(const std::vector<wmtk::TetMesh::Tuple>&)>& try_acquire_edge,
    const std::function<bool(const std::vector<wmtk::TetMesh::Tuple>&)>& try_acquire_tetra)
{
    auto vertex_pos_r = [&m_vertex_attribute](size_t i) -> tetwild::Vector3r {
        return m_vertex_attribute[i].m_pos;
    };

    auto face_verts = std::array<Eigen::Vector3d, 3>();
    for (auto j = 0; j < 3; j++) face_verts[j] = vertices[face[j]];

    std::vector<wmtk::TetMesh::Tuple> initial_tets;
    const auto& [flag, intersected_tets, intersected_edges, intersected_pos] =
        wmtk::triangle_insert_prepare_info_non_point<apps::Rational>(
            m,
            face_verts,
            marked_tet_faces, // output
            try_acquire_tetra,
            vertex_pos_r,
            initial_tets);

    if (!flag) {
        return false;
    }

    if (try_acquire_edge(intersected_edges) == false ||
        try_acquire_tetra(intersected_tets) == false) {
        return false;
    }

    // these are only those on edges.
    std::vector<size_t> new_edge_vids;
    std::vector<size_t> new_center_vids;
    std::vector<std::array<size_t, 4>> center_split_tets;

    ///inert a triangle
    m.triangle_insertion(
        intersected_tets,
        intersected_edges,
        new_edge_vids,
        new_center_vids,
        center_split_tets);

    assert(new_center_vids.size() == center_split_tets.size());
    for (auto i = 0; i < new_center_vids.size(); i++) {
        auto vid = new_center_vids[i];
        auto& vs = center_split_tets[i];
        m_vertex_attribute[vid] = tetwild::VertexAttributes(tetwild::Vector3r(
            (m_vertex_attribute[vs[0]].m_pos + m_vertex_attribute[vs[1]].m_pos +
             m_vertex_attribute[vs[2]].m_pos + m_vertex_attribute[vs[3]].m_pos) /
            4));
    }
    assert(new_edge_vids.size() == intersected_pos.size());

    for (auto i = 0; i < intersected_pos.size(); i++) {
        m_vertex_attribute[new_edge_vids[i]] = tetwild::VertexAttributes(intersected_pos[i]);
    }

    return true;
}

void wmtk::AdaMesh::insert_triangles_to_mesh(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::array<size_t, 3>>& faces)
{
    // match faces preserved in delaunay
    tbb::concurrent_vector<bool> is_matched;
    // wmtk::match_tet_faces_to_triangles(*this, faces, is_matched, tet_face_tags);
    // wmtk::logger().info("is_matched: {}", std::count(is_matched.begin(), is_matched.end(),
    // true));
    std::vector<tbb::concurrent_priority_queue<std::tuple<double, int, size_t>>> insertion_queues(
        NUM_THREADS);
    tbb::concurrent_priority_queue<std::tuple<double, int, size_t>> expired_queue;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 100.0);
    for (size_t face_id = 0; face_id < faces.size(); face_id++) {
        if (is_matched[face_id]) continue;
        double rand = distribution(generator);
        insertion_queues[0].emplace(rand, 0, face_id);
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        wmtk::logger().debug("{}: {}", i, insertion_queues[i].size());
    }

    tbb::task_arena arena(NUM_THREADS);
    tbb::task_group tg;

    arena.execute([&, &m = *this, &tet_face_tags = this->tet_face_tags]() {
        for (int task_id = 0; task_id < m.NUM_THREADS; task_id++) {
            tg.run([&insertion_queues, &expired_queue, &tet_face_tags, &m, &faces, task_id] {
                auto try_acquire_tetra = [&m, task_id](const auto& intersected_tets) {
                    for (auto t_int : intersected_tets) {
                        for (auto v_int : m.oriented_tet_vertices(t_int)) {
                            if (!m.try_set_vertex_mutex_one_ring(v_int, task_id)) {
                                return false;
                            }
                        }
                    }
                    return true;
                };

                auto try_acquire_edge = [&m, task_id](const auto& intersected_edges) {
                    for (auto e_int : intersected_edges) {
                        if (!m.try_set_vertex_mutex_one_ring(e_int, task_id)) {
                            return false;
                        }
                        if (!m.try_set_vertex_mutex_one_ring(e_int.switch_vertex(m), task_id)) {
                            return false;
                        }
                    }
                    return true;
                };

                std::default_random_engine generator;
                std::uniform_real_distribution<double> distribution(0.0, 100.0);
                auto retry_processing = [&,
                                         &Q = insertion_queues[task_id]](auto id, auto retry_time) {
                    double rand = distribution(generator);
                    if (retry_time < 5) {
                        Q.push(std::make_tuple(rand, retry_time + 1, id));
                    } else {
                        expired_queue.push(std::make_tuple(rand, 0, id));
                    }
                };
                auto try_acquire_triangle = [&m, task_id](const auto& f) {
                    return m.try_set_face_mutex_two_ring(f[0], f[1], f[2], task_id);
                };

                auto supply_element = [&m,
                                       &Q = insertion_queues[task_id],
                                       &face_id_cache =
                                           m.triangle_insertion_local_cache.local().face_id,
                                       &retry_processing](const auto& func) {
                    std::tuple<double, int, size_t> eiq;
                    while (Q.try_pop(eiq)) {
                        const auto& [_, retry_time, face_id] = eiq;

                        face_id_cache = face_id;
                        if (func(face_id) == false) {
                            retry_processing(face_id, retry_time);
                            continue;
                        };
                        m.release_vertex_mutex_in_stack();
                    }
                };

                supply_element([&](auto face_id) {
                    std::vector<std::array<size_t, 3>> marked_tet_faces;
                    auto success = internal_insert_single_triangle(
                        m,
                        m.m_vertex_attribute,
                        vertices,
                        faces[face_id],
                        marked_tet_faces,
                        try_acquire_edge,
                        try_acquire_tetra);
                    if (!success) return false;
                    for (auto& f : marked_tet_faces) tet_face_tags[f].push_back(face_id);
                    return true;
                });
            }); // tg.run
        } // parallel for loop
    });
    arena.execute([&] { tg.wait(); });

    wmtk::logger().info("expired size: {}", expired_queue.size());

    auto check_acquire = [](const auto&) { return true; };

    auto supply_element = [&Q = expired_queue,
                           &face_id_cache =
                               triangle_insertion_local_cache.local().face_id](const auto& func) {
        std::tuple<double, int, size_t> eiq;
        while (Q.try_pop(eiq)) {
            const auto& [_, retry_time, face_id] = eiq;

            face_id_cache = face_id;
            if (func(face_id) == false) {
                continue;
            };
        }
    };

    supply_element([&](auto face_id) {
        std::vector<std::array<size_t, 3>> marked_tet_faces;
        auto success = internal_insert_single_triangle(
            *this,
            m_vertex_attribute,
            vertices,
            faces[face_id],
            marked_tet_faces,
            check_acquire,
            check_acquire);
        if (!success) return false;
        for (auto& f : marked_tet_faces) tet_face_tags[f].push_back(face_id);
        return true;
    });
}
} // namespace wmtk