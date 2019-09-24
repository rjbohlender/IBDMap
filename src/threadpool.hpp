//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_THREADPOOL_HPP
#define CARVAIBD_THREADPOOL_HPP

#include "parameters.hpp"
#include <atomic>

class ThreadPool {
  std::atomic_bool done;
public:
  ThreadPool(Parameters params);

};

#endif //CARVAIBD_THREADPOOL_HPP
