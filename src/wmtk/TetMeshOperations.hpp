#pragma once
#include "TetMesh.h"

// Customizable Operations
namespace wmtk {
struct SplitEdge : public TetMesh::OperationBuilder
{
    const TetMesh& m;
    std::vector<size_t> affected;
    std::array<size_t, 2> edge_verts;
    size_t ux;

    SplitEdge(const TetMesh& _m)
        : m(_m)
    {}
    std::vector<size_t> removed_tids(const TetMesh::Tuple& t)
    {
        auto incidents = m.get_incident_tets_for_edge(t);
        for (auto i : incidents) {
            affected.push_back(i.tid(m));
        }
        edge_verts = {t.vid(m), t.switch_vertex(m).vid(m)};
        logger().trace("{}", edge_verts);

        return affected;
    }
    int request_vert_slots() { return 1; }
    std::vector<std::array<size_t, 4>> replacing_tets(const std::vector<size_t>& slots)
    {
        assert(slots.size() == 1);
        ux = slots.front();

        auto num = affected.size();
        std::vector<std::array<size_t, 4>> new_tet_conn(2 * num);

        for (auto i = 0; i < affected.size(); i++) {
            auto t_id = affected[i];
            auto t_conn = m.oriented_tet_vids(m.tuple_from_tet(affected[i]));
            for (auto j = 0; j < 2; j++) {
                auto& curent_t = new_tet_conn[i + num * j];
                curent_t = t_conn;
                // 1-j is because of conforming to existing tests.
                std::replace(curent_t.begin(), curent_t.end(), edge_verts[1 - j], ux);
            }
        }
        return new_tet_conn;
    }
};
} // namespace wmtk