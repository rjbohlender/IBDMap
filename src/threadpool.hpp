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
  std::atomic_bool done;
  std::atomic<int> ntasks;
  ThreadSafeQueue<std::packaged_task<Res(Args...)>> work_queue;
  std::vector<std::thread> threads;
  JoinThreads joiner;
  void worker_thread()
  {
    while(!done || ntasks > 0)
	{
      std::packaged_task<Res(Args...)> task;
      if (work_queue.try_pop(task))
      {
        task();
        ntasks--;
      }
      else
	  {
        std::this_thread::yield();
	  }
	}
  }
public:
  explicit ThreadPool(Parameters params) :
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
	while (!work_queue.empty() || ntasks > 0) {
	  std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
	}
    done=true;
  }

  void submit(std::packaged_task<Res(Args...)> f)
  {
    work_queue.push(std::move(f));
    ntasks++;
  }
};

#endif //CARVAIBD_THREADPOOL_HPP
