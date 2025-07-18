#pragma once

#include "wmtk/TetMesh.h"
#include "wmtk/TriMesh.h"
#include "wmtk/utils/Logger.hpp"

// clang-format off
#include <functional>
#include <limits>
#include <wmtk/utils/DisableWarnings.hpp>
#include <tbb/concurrent_priority_queue.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#include <wmtk/utils/EnableWarnings.hpp>
// clang-format on


#include <atomic>
#include <cassert>
#include <cstddef>
#include <queue>
#include <stdexcept>
#include <type_traits>

namespace wmtk {
enum class ExecutionPolicy { kSeq, kUnSeq, kPartition, kColor, kMax };

using Op = std::string;

template <class AppMesh, ExecutionPolicy policy = ExecutionPolicy::kSeq>
struct ExecutePass
{
    using Tuple = typename AppMesh::Tuple;
    /**
     * @brief A dictionary that registers names with operations.
     *
     */
    std::map<
        Op, // strings
        std::function<std::optional<std::vector<Tuple>>(AppMesh&, const Tuple&)>>
        edit_operation_maps;
    /**
     * @brief Priority function (default to edge length)
     *
     */
    std::function<double(const AppMesh&, Op op, const Tuple&)> priority = [](auto&, auto, auto&) {
        return 0.;
    };
    /**
     * @brief check on wheather new operations should be added to the prioirity queue
     *
     */
    std::function<bool(double)> should_renew = [](auto) { return true; };
    /**
     * @brief renew neighboring Tuples after each operation depends on the operation
     *
     */
    std::function<std::vector<std::pair<Op, Tuple>>(const AppMesh&, Op, const std::vector<Tuple>&)>
        renew_neighbor_tuples =
            [](auto&, auto, auto&) -> std::vector<std::pair<Op, Tuple>> { return {}; };
    /**
     * @brief lock the vertices concerned depends on the operation
     *
     */
    std::function<bool(AppMesh&, const Tuple&, int task_id)> lock_vertices =
        [](const AppMesh&, const Tuple&, int task_id) { return true; };
    /**
     * @brief Stopping Criterion based on the whole mesh
        For efficiency, not every time is checked.
        In serial, this may go over all the elements. For parallel, this involves synchronization.
        So there is a checking frequency.
     *
     */
    std::function<bool(const AppMesh&)> stopping_criterion = [](const AppMesh&) {
        return false; // non-stop, process everything
    };
    /**
     * @brief checking frequency to decide whether to stop execution given the stopping criterion
     *
     */
    size_t stopping_criterion_checking_frequency = std::numeric_limits<size_t>::max();
    /**
     * @brief Should Process drops some Tuple from being processed.
         For example, if the energy is out-dated.
         This is in addition to calling tuple valid.
     *
     */
    std::function<bool(const AppMesh&, const std::tuple<double, Op, Tuple>& t)>
        is_weight_up_to_date = [](const AppMesh& m, const std::tuple<double, Op, Tuple>& t) {
            // always do.
            assert(std::get<2>(t).is_valid(m));
            return true;
        };
    /**
     * @brief
     *used to collect operations that are not finished and used for later re-execution
     */
    std::function<void(const AppMesh&, Op, const Tuple& t)> on_fail = [](auto&, auto, auto&) {};


    size_t num_threads = 1;

    /**
     * To Avoid mutual locking, retry limit is set, and then put in a serial queue in the end.
     *
     */
    size_t max_retry_limit = 10;
    /**
     * @brief Construct a new Execute Pass object. It contains the name-to-operation map and the
     *functions that define the rules for operations
     *@note the constructor is differentiated by the type of mesh, namingly wmtk::TetMesh or
     *wmtk::TriMesh
     */
    ExecutePass()
    {
        if constexpr (std::is_base_of<wmtk::TetMesh, AppMesh>::value) {
            edit_operation_maps = {
                {"edge_collapse",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.collapse_edge(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"edge_swap",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.swap_edge(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"edge_swap_44",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.swap_edge_44(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"edge_split",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.split_edge(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"face_swap",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.swap_face(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"vertex_smooth",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     if (m.smooth_vertex(t))
                         return std::vector<Tuple>{};
                     else
                         return {};
                 }}};
        }
        if constexpr (std::is_base_of<wmtk::TriMesh, AppMesh>::value) {
            edit_operation_maps = {
                {"edge_collapse",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.collapse_edge(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"edge_swap",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.swap_edge(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"edge_split",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     std::vector<Tuple> ret;
                     if (m.split_edge(t, ret))
                         return ret;
                     else
                         return {};
                 }},
                {"vertex_smooth",
                 [](AppMesh& m, const Tuple& t) -> std::optional<std::vector<Tuple>> {
                     if (m.smooth_vertex(t))
                         return std::vector<Tuple>{};
                     else
                         return {};
                 }}};
        }
    };

private:
    void operation_cleanup(AppMesh& m)
    { //
        // class ResourceManger
        // what about RAII mesh edit locking?
        // release mutex, but this should be implemented in TetMesh class.
        if constexpr (policy == ExecutionPolicy::kSeq)
            return;
        else {
            m.release_vertex_mutex_in_stack();
        }
    }

    size_t get_partition_id(const AppMesh& m, const Tuple& e)
    {
        if constexpr (policy == ExecutionPolicy::kSeq) return 0;
        if constexpr (std::is_base_of<wmtk::TetMesh, AppMesh>::value)
            return m.get_partition_id(e);
        else if constexpr (std::is_base_of<wmtk::TriMesh, AppMesh>::value) // TODO: make same
                                                                           // interface.
            return m.vertex_attrs[e.vid(m)].partition_id; // TODO: this is temporary.
        return 0;
    }

public:
    /**
     * @brief Executes the operations for an application when the lambda function is invoked. The
     * rules that are customizly defined for applications are applied.
     *
     * @param m
     * @param operation_tuples a vector of pairs of operation's name and the Tuple to be operated on
     * @returns true if finished successfully
     */
    bool operator()(AppMesh& m, const std::vector<std::pair<Op, Tuple>>& operation_tuples)
    {
        auto cnt_update = std::atomic<int>(0);
        m_cnt_success = 0;
        m_cnt_fail = 0;
        auto stop = std::atomic<bool>(false);
        m_parallel_queues = std::vector<tbb::concurrent_priority_queue<Elem>>(num_threads);
        m_serial_queue = tbb::concurrent_priority_queue<Elem>();

        auto run_single_queue = [&](auto& Q, int task_id) {
            auto ele_in_queue = Elem();
            while ([&]() { return Q.try_pop(ele_in_queue); }()) {
                auto& [weight, op, tup, retry] = ele_in_queue;
                if (!tup.is_valid(m)) continue;

                std::vector<Elem> renewed_elements;
                {
                    // Note that returning `Tuples` would be invalid.
                    auto locked_vid = lock_vertices(m, tup, task_id);
                    if (!locked_vid) {
                        retry++;
                        if (retry < max_retry_limit) {
                            Q.emplace(ele_in_queue);
                        } else {
                            retry = 0;
                            m_serial_queue.emplace(ele_in_queue);
                        }
                        continue;
                    }
                    if (tup.is_valid(m)) {
                        if (!is_weight_up_to_date(
                                m,
                                std::tuple<double, Op, Tuple>(
                                    std::get<0>(ele_in_queue),
                                    std::get<1>(ele_in_queue),
                                    std::get<2>(ele_in_queue)))) {
                            operation_cleanup(m);
                            continue;
                        } // this can encode, in qslim, recompute(energy) == weight.
                        auto newtup = edit_operation_maps[op](m, tup);
                        std::vector<std::pair<Op, Tuple>> renewed_tuples;
                        if (newtup) {
                            renewed_tuples = renew_neighbor_tuples(m, op, newtup.value());
                            m_cnt_success++;
                            cnt_update++;
                        } else {
                            on_fail(m, op, tup);
                            m_cnt_fail++;
                        }
                        for (auto& [o, e] : renewed_tuples) {
                            auto val = priority(m, o, e);
                            if (should_renew(val)) renewed_elements.emplace_back(val, o, e, 0);
                        }
                    }
                    operation_cleanup(m); // Maybe use RAII
                }
                for (auto& e : renewed_elements) {
                    Q.emplace(e);
                }

                if (stop.load(std::memory_order_acquire)) return;
                if (m_cnt_success > stopping_criterion_checking_frequency) {
                    if (stopping_criterion(m)) {
                        stop.store(true);
                        return;
                    }
                    cnt_update.store(0, std::memory_order_release);
                }
            }
        };

        if constexpr (policy == ExecutionPolicy::kSeq) {
            for (auto& [op, e] : operation_tuples) {
                if (!e.is_valid(m)) continue;
                m_serial_queue.emplace(priority(m, op, e), op, e, 0);
            }
            run_single_queue(m_serial_queue, 0);
        } else {
            for (auto& [op, e] : operation_tuples) {
                if (!e.is_valid(m)) continue;
                m_parallel_queues[get_partition_id(m, e)].emplace(priority(m, op, e), op, e, 0);
            }
            // Comment out parallel: work on serial first.
            tbb::task_arena arena(num_threads);
            tbb::task_group tg;
            arena.execute([this, &run_single_queue, &tg]() {
                for (int task_id = 0; task_id < this->m_parallel_queues.size(); task_id++) {
                    tg.run([&run_single_queue, this, task_id] {
                        run_single_queue(this->m_parallel_queues[task_id], task_id);
                    });
                }
                tg.wait();
            });
            logger().debug("Parallel Complete, remains element {}", m_serial_queue.size());
            run_single_queue(m_serial_queue, 0);
        }

        logger().info("cnt_success {} cnt_fail {}", cnt_success(), cnt_fail());
        return true;
    }

public:
    int cnt_success() const { return m_cnt_success; }
    int cnt_fail() const { return m_cnt_fail; }

    using Elem = std::tuple<double, Op, Tuple, size_t>;

    const tbb::concurrent_priority_queue<Elem>& serial_queue() const { return m_serial_queue; }

    const std::vector<tbb::concurrent_priority_queue<Elem>>& parallel_queues() const
    {
        return m_parallel_queues;
    }

    std::atomic<int> m_cnt_success = std::atomic<int>(0);
    std::atomic<int> m_cnt_fail = std::atomic<int>(0);

protected:
    tbb::concurrent_priority_queue<Elem> m_serial_queue;
    std::vector<tbb::concurrent_priority_queue<Elem>> m_parallel_queues;
};
} // namespace wmtk
