#include "AdaMesh.hpp"
#include "prism/predicates/inside_prism_tetra.hpp"
#include "wmtk/utils/InsertTriangleUtils.hpp"

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

    m.for_each_tetra([&](auto& t) {
        auto vs = m.oriented_tet_vids(t);
        for (auto j = 0; j < 3; j++) {
            if (::prism::predicates::point_in_tetrahedron(
                    face_verts[j],
                    m_vertex_attribute[vs[0]].m_posf,
                    m_vertex_attribute[vs[1]].m_posf,
                    m_vertex_attribute[vs[2]].m_posf,
                    m_vertex_attribute[vs[3]].m_posf)) {
                initial_tets.push_back(t);
                return;
            }
        }
    });
    // wmtk::logger().info("size tets {}", initial_tets.size());


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
    // tbb::concurrent_vector<bool> is_matched;
    // wmtk::match_tet_faces_to_triangles(*this, faces, is_matched, tet_face_tags);
    // wmtk::logger().info("is_matched: {}", std::count(is_matched.begin(), is_matched.end(),
    // true));
    std::vector<tbb::concurrent_priority_queue<std::tuple<double, int, size_t>>> insertion_queues(
        NUM_THREADS);
    tbb::concurrent_priority_queue<std::tuple<double, int, size_t>> expired_queue;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 100.0);
    for (size_t face_id = 0; face_id < faces.size(); face_id++) {
        // if (is_matched[face_id]) continue;
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
            tg.run([&insertion_queues,
                    &expired_queue,
                    &tet_face_tags,
                    &m,
                    &faces,
                    &vertices,
                    task_id] {
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
                    assert(success);
                    if (!success) return false;
                    for (auto& f : marked_tet_faces) tet_face_tags[f].push_back(face_id);
                    // wmtk::logger().critical("size {}", tet_face_tags.size());
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
    wmtk::logger().critical("T F T {}", tet_face_tags.size());

    finalize_triangle_insertion(vertices, faces);
}