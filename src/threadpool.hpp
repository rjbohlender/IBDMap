//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_THREADPOOL_HPP
#define CARVAIBD_THREADPOOL_HPP

#include "jointhreads.hpp"
#include "parameters.hpp"
#include "threadsafequeue.hpp"
#include <atomic>
#include <future>
#include <thread>
#include <vector>

/**
 * @brief  A generic threadpool that calls a run member function.
 * @tparam Stat The statistic object that will have it's run function called.
 *
 * @note Origiinally we used two template parameters, Res and Args for a packaged
 * task. I'm not directly using the Statistic objects to avoid the memory accumulation
 * problem that results in excessive memory usage. We should now be able to analyze
 * hapibd data easily, with little memory usage, despite the large number of breakpoints.
 *
 * The original approach was a result of wanting to preserve all the data in memory
 * so that we could combine and calculate multiple testing corrected p-values in
 * carvaIBD. However, we now use a script to post-process carvaIBD's output. Thus,
 * we no longer need to keep the objects in memory, and we can free up that burden.
 */
template<typename Stat>
class ThreadPool {
    std::atomic<int> nsubmitted;
    mutable std::mutex mut;
    std::condition_variable not_empty;
    std::condition_variable not_full;
    std::condition_variable all_tasks_done;
    std::queue<Stat> work_queue;
    std::vector<std::thread> threads;
    JoinThreads joiner;
    Parameters params;
    void worker_thread() {
        std::unique_lock<std::mutex> lk(mut);
        while (!done || ntasks > 0) {
            not_empty.wait_for(lk, std::chrono::seconds(1), [this] { return !work_queue.empty() || done; });
            if (!work_queue.empty()) {
                Stat task = std::move(work_queue.front());
                work_queue.pop();

                // Done with queue; unlock
                lk.unlock();
                not_full.notify_one();
                task.run();

                // Re-acquire lock before modifying ntasks and notifying
                lk.lock();
                ntasks--;

                // Notify if all tasks are complete
                if (ntasks == 0) {
                    all_tasks_done.notify_all();
                }
                // Lock will be reacquired at top of loop
            }
        }
    }

public:
    std::atomic<bool> done;
    std::atomic<int> ntasks;
    explicit ThreadPool(Parameters params_) : done(false), ntasks(0), nsubmitted(0), joiner(threads), params(std::move(params_)) {
        unsigned const thread_count = params.nthreads > 2 ? params.nthreads - 2 : 1;
        try {
            for (unsigned i = 0; i < thread_count; ++i) {
                threads.emplace_back(
                        std::thread(&ThreadPool::worker_thread, this));
            }
        } catch (...) {
            done = true;
            throw;
        }
    }
    ~ThreadPool() {
        done = true;
        not_empty.notify_all();
        wait_for_completion();
    }

    void submit(Stat &&f) {
        std::unique_lock lk(mut);
        not_full.wait(lk, [this]() { return ntasks < params.nthreads + 5; });
        work_queue.push(std::move(f));
        ntasks++;
        nsubmitted++;
        lk.unlock();
        not_empty.notify_one();
    }

    void wait_for_completion() {
        std::unique_lock<std::mutex> lk(mut);
        all_tasks_done.wait(lk, [this]() { return work_queue.empty() && ntasks == 0; });
    }
};

#endif//CARVAIBD_THREADPOOL_HPP
