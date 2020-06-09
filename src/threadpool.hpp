//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_THREADPOOL_HPP
#define CARVAIBD_THREADPOOL_HPP

#include "parameters.hpp"
#include "threadsafequeue.hpp"
#include "jointhreads.hpp"
#include <atomic>
#include <vector>
#include <thread>
#include <future>

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
  std::atomic<int> ntasks;
  std::atomic<int> nsubmitted;
  mutable std::mutex mut;
  std::condition_variable cv;
  std::queue<Stat> work_queue;
  std::vector<std::thread> threads;
  JoinThreads joiner;
  Parameters params;
  void worker_thread()
  {
    std::unique_lock<std::mutex> lk(mut);
    while(!done || ntasks > 0)
	{
      cv.wait_for(lk, std::chrono::seconds(1), [this] { return !work_queue.empty() || done; });
      if(!work_queue.empty()) {
        Stat task = std::move(work_queue.front());
        work_queue.pop();

        // Done with queue; unlock
        lk.unlock();
        task.run();
        ntasks--;

        lk.lock();
      }
	}
  }
public:
  std::atomic<bool> done;
  explicit ThreadPool(Parameters params_) :
  	done(false), ntasks(0), nsubmitted(0), joiner(threads), params(std::move(params_))
  {
    unsigned const thread_count=params.nthreads - 2;
	try
	{
	  for(unsigned i = 0; i < thread_count; ++i)
	  {
	    threads.emplace_back(
	    	std::thread(&ThreadPool::worker_thread, this));
	  }
	}
	catch (...)
	{
	  done=true;
	  throw;
	}
  }
  ~ThreadPool()
  {
    done=true;
    cv.notify_all();
	while (!work_queue.empty() || ntasks > 0) {
	  std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
	}
  }

  void submit(Stat &&f)
  {
    std::unique_lock lk(mut);
    while(ntasks > params.nthreads + 5) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
	}
	work_queue.push(std::move(f));
	ntasks++;
    nsubmitted++;
    lk.unlock();
    cv.notify_all();
  }
};

#endif //CARVAIBD_THREADPOOL_HPP
