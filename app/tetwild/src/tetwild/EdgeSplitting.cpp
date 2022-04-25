#include "TetWild.h"

#include <igl/Timer.h>
#include <wmtk/ExecutionScheduler.hpp>
#include <wmtk/utils/ExecutorUtils.hpp>
#include <wmtk/utils/Logger.hpp>

namespace tetwild {
bool SplitEdge::before(const wmtk::TetMesh::Tuple& loc0)
{
    auto& m = static_cast<const tetwild::TetWild&>(this->m);
    split_cache.changed_faces.clear();

    split_cache.v1_id = loc0.vid(m);
    auto loc1 = loc0.switch_vertex(m);
    split_cache.v2_id = loc1.vid(m);
    //
    size_t v1_id = split_cache.v1_id;
    size_t v2_id = split_cache.v2_id;

    split_cache.is_edge_on_surface = m.is_edge_on_surface(loc0);

    /// save face track info
    auto comp = [](const std::pair<FaceAttributes, std::array<size_t, 3>>& v1,
                   const std::pair<FaceAttributes, std::array<size_t, 3>>& v2) {
        return v1.second < v2.second;
    };
    auto is_equal = [](const std::pair<FaceAttributes, std::array<size_t, 3>>& v1,
                       const std::pair<FaceAttributes, std::array<size_t, 3>>& v2) {
        return v1.second == v2.second;
    };

    auto tets = m.get_incident_tets_for_edge(loc0);
    for (auto& t : tets) {
        auto vs = m.oriented_tet_vertices(t);
        for (int j = 0; j < 4; j++) {
            std::array<size_t, 3> f_vids = {{
                vs[(j + 1) % 4].vid(m),
                vs[(j + 2) % 4].vid(m),
                vs[(j + 3) % 4].vid(m),
            }}; // todo: speedup
            std::sort(f_vids.begin(), f_vids.end());
            auto [_, global_fid] = m.tuple_from_face(f_vids);
            split_cache.changed_faces.push_back(
                std::make_pair(m.m_face_attribute[global_fid], f_vids));
        }
    }
    wmtk::vector_unique(split_cache.changed_faces, comp, is_equal);
    return true;
}
bool SplitEdge::after(const std::vector<wmtk::TetMesh::Tuple>& locs)
{
    size_t v_id = ux;

    size_t v1_id = split_cache.v1_id;
    size_t v2_id = split_cache.v2_id;

    /// check inversion & rounding
    m_app.m_vertex_attribute[v_id].m_posf =
        (m_app.m_vertex_attribute[v1_id].m_posf + m_app.m_vertex_attribute[v2_id].m_posf) / 2;
    m_app.m_vertex_attribute[v_id].m_is_rounded = true;

    for (auto& loc : locs) {
        if (m_app.is_inverted(loc)) {
            m_app.m_vertex_attribute[v_id].m_is_rounded = false;
            break;
        }
    }
    if (!m_app.m_vertex_attribute[v_id].m_is_rounded) {
        m_app.m_vertex_attribute[v_id].m_pos =
            (m_app.m_vertex_attribute[v1_id].m_pos + m_app.m_vertex_attribute[v2_id].m_pos) / 2;
        m_app.m_vertex_attribute[v_id].m_posf = to_double(m_app.m_vertex_attribute[v_id].m_pos);
    } else
        m_app.m_vertex_attribute[v_id].m_pos = to_rational(m_app.m_vertex_attribute[v_id].m_posf);

    /// update quality
    for (auto& loc : locs) {
        m_app.m_tet_attribute[loc.tid(m)].m_quality = m_app.get_quality(loc);
    }

    /// update vertex attribute
    // bbox
    m_app.m_vertex_attribute[v_id].on_bbox_faces = wmtk::set_intersection(
        m_app.m_vertex_attribute[v1_id].on_bbox_faces,
        m_app.m_vertex_attribute[v2_id].on_bbox_faces);
    // surface
    m_app.m_vertex_attribute[v_id].m_is_on_surface = split_cache.is_edge_on_surface;

    /// update face attribute
    // add new and erase old
    for (auto& info : split_cache.changed_faces) {
        auto& f_attr = info.first;
        auto& old_vids = info.second;
        std::vector<int> j_vn;
        for (int j = 0; j < 3; j++) {
            if (old_vids[j] != v1_id && old_vids[j] != v2_id) {
                j_vn.push_back(j);
            }
        }
        if (j_vn.size() == 1) {
            auto [_1, global_fid1] = m.tuple_from_face({{v1_id, v_id, old_vids[j_vn[0]]}});
            m_app.m_face_attribute[global_fid1] = f_attr;
            auto [_2, global_fid2] = m.tuple_from_face({{v2_id, v_id, old_vids[j_vn[0]]}});
            m_app.m_face_attribute[global_fid2] = f_attr;
        } else { // j_vn.size() == 2
            auto [_, global_fid] = m.tuple_from_face(old_vids);
            m_app.m_face_attribute[global_fid] = f_attr;
            //
            auto [_2, global_fid2] = m.tuple_from_face(
                {{old_vids[j_vn[0]], old_vids[j_vn[1]], v_id}}); // todo: avoid dup comp
            m_app.m_face_attribute[global_fid2].reset();
        }
    }

    m_app.m_vertex_attribute[v_id].partition_id = m_app.m_vertex_attribute[v1_id].partition_id;
    m_app.m_vertex_attribute[v_id].m_sizing_scalar =
        (m_app.m_vertex_attribute[v1_id].m_sizing_scalar +
         m_app.m_vertex_attribute[v2_id].m_sizing_scalar) /
        2;

    m_app.cnt_split++;

    return true;
}
} // namespace tetwild

void tetwild::TetWild::split_all_edges()
{
    igl::Timer timer;
    double time;
    timer.start();
    auto collect_all_ops = std::vector<std::pair<std::string, Tuple>>();
    for (auto& loc : get_edges()) collect_all_ops.emplace_back("edge_split", loc);
    time = timer.getElapsedTime();
    wmtk::logger().info("edge split prepare time: {}s", time);
    auto setup_and_execute = [&](auto& executor) {
        executor.edit_operation_maps["edge_split"] =
            [](auto& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
            std::vector<Tuple> ret;
            tetwild::SplitEdge sp(m);
            if (m.customized_operation(sp, t, ret))
                return ret;
            else
                return {};
        };
        executor.renew_neighbor_tuples = wmtk::renewal_simple;

        executor.priority = [&](auto& m, auto op, auto& t) { return m.get_length2(t); };
        executor.num_threads = NUM_THREADS;
        executor.should_process = [&](const auto& m, const auto& ele) {
            auto [weight, op, tup] = ele;
            auto length = m.get_length2(tup);
            if (length != weight) return false;
            //
            size_t v1_id = tup.vid(*this);
            size_t v2_id = tup.switch_vertex(*this).vid(*this);
            double sizing_ratio = (m_vertex_attribute[v1_id].m_sizing_scalar +
                                   m_vertex_attribute[v2_id].m_sizing_scalar) /
                                  2;
            if (length < m_params.splitting_l2 * sizing_ratio * sizing_ratio) return false;
            return true;
        };
        executor(*this, collect_all_ops);
    };
    if (NUM_THREADS > 0) {
        timer.start();
        auto executor = wmtk::ExecutePass<TetWild, wmtk::ExecutionPolicy::kPartition>();
        executor.lock_vertices = [&](auto& m, const auto& e, int task_id) -> bool {
            return m.try_set_edge_mutex_two_ring(e, task_id);
        };
        setup_and_execute(executor);
        time = timer.getElapsedTime();
        wmtk::logger().info("edge split operation time parallel: {}s", time);
    } else {
        timer.start();
        auto executor = wmtk::ExecutePass<TetWild, wmtk::ExecutionPolicy::kSeq>();
        setup_and_execute(executor);
        time = timer.getElapsedTime();
        wmtk::logger().info("edge split operation time serial: {}s", time);
    }
}
