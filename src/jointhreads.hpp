//
// Created by Bohlender,Ryan James on 9/4/19.
//

#ifndef CARVAIBD_JOINTHREADS_HPP
#define CARVAIBD_JOINTHREADS_HPP

#include <vector>
#include <thread>

class JoinThreads {
  std::vector<std::thread>& threads;
public:
  JoinThreads(std::vector<std::thread>& threads_)
  : threads(threads_) {}
  ~JoinThreads() {
    for(unsigned long i = 0; i < threads.size(); ++i) {
      if(threads[i].joinable()) {
        threads[i].join();
      }
    }
  }
};

#endif //CARVAIBD_JOINTHREADS_HPP
