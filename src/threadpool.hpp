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
  mutable std::mutex mut;
  std::condition_variable cv;
  ThreadSafeQueue<std::packaged_task<Res(Args...)>> work_queue;
  std::vector<std::thread> threads;
  JoinThreads joiner;
  void worker_thread()
  {
    std::unique_lock<std::mutex> lk(mut);
    while(!done || ntasks > 0)
	{
      std::packaged_task<Res(Args...)> task;
      cv.wait(lk, [this] { return !work_queue.empty() || done; });
      if (work_queue.try_pop(task))
      {
        lk.unlock();
        task();
        ntasks--;
        lk.lock();
      }
      else
	  {
        std::this_thread::yield();
	  }
	}
  }
public:
  explicit ThreadPool(const Parameters &params) :
  	done(false), ntasks(0), joiner(threads)
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
    work_queue.push(std::move(f));
    ntasks++;
    cv.notify_all();
  }

  int pending() {
    return ntasks;
  }
};

#endif //CARVAIBD_THREADPOOL_HPP
