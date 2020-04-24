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

template<typename Res, typename... Args>
class ThreadPool {
  std::atomic<bool> done;
  std::atomic<int> ntasks;
  std::atomic<int> nsubmitted;
  mutable std::mutex mut;
  std::condition_variable cv;
  std::queue<std::packaged_task<Res(Args...)>> work_queue;
  std::vector<std::thread> threads;
  JoinThreads joiner;
  Parameters params;
  void worker_thread()
  {
    std::unique_lock<std::mutex> lk(mut);
    while(!done || ntasks > 0)
	{
      std::packaged_task<Res(Args...)> task;
      cv.wait_for(lk, std::chrono::seconds(1), [this] { return !work_queue.empty() || done; });
      if(!work_queue.empty()) {
        task = std::move(work_queue.front());
        work_queue.pop();

        auto f = task.get_future();
        // Done with queue; unlock
        lk.unlock();
        task();
        ntasks--;

        lk.lock();

        f.get(); // Throw any errors -- holy fuck finally this is it
      }
	}
  }
public:
  explicit ThreadPool(const Parameters &params_) :
  	done(false), ntasks(0), nsubmitted(0), joiner(threads), params(params_)
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

  void submit(std::packaged_task<Res(Args...)> &&f)
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
