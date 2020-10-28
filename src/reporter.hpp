//
// Created by Bohlender,Ryan James on 10/2/19.
//

#ifndef CARVAIBD_REPORTER_HPP
#define CARVAIBD_REPORTER_HPP

#include <mutex>
#include <iostream>
#include <thread>
#include <queue>
#include <fmt-7.0.3/include/fmt/ostream.h>
#include "threadsafequeue.hpp"

/**
 * @brief Thread safe shared writing class
 */
class Reporter {
  std::atomic_bool done;
  std::atomic<int> nstrings;
  std::atomic<int> nsubmitted;
  std::atomic<int> nwritten;
  std::string out_path;
  std::thread print_thread;
  std::queue<std::string> string_queue;
  mutable std::mutex mut;
  std::condition_variable data_cond;

  void print() {
	std::unique_lock<std::mutex> lk(mut);
	std::string s;
	std::ofstream ofs;
	if (!out_path.empty()) {
	  ofs = std::ofstream(out_path);
	}
	while (!done || nstrings > 0) {
	  data_cond.wait_for(lk, std::chrono::seconds(1), [this] { return !string_queue.empty() || done; });
	  if (!string_queue.empty()) {
		s = std::move(string_queue.front());
		string_queue.pop();
		lk.unlock();
		if (!out_path.empty()) {
		  fmt::print(ofs, s);
		} else {
		  fmt::print(s);
		}
		nstrings--;
		nwritten++;
		lk.lock();
	  }
	}
  }

public:
  explicit Reporter(std::string output) :
	  done(false), nstrings(0), nsubmitted(0), nwritten(0), out_path(std::move(output)) {
	print_thread = std::thread(&Reporter::print, std::ref(*this));
  }

  ~Reporter() {
	while (!string_queue.empty() || nstrings > 0) {
	  fmt::print(std::cerr, "queue size: {} nstrings: {}\n", string_queue.size(), nstrings);
	  std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
	}
	done = true;
	data_cond.notify_all();
	if (print_thread.joinable()) {
	  fmt::print(std::cerr, "In reporter:\nnwritten: {}\nnsubmitted: {}\n", nwritten, nsubmitted);
	  print_thread.join();
	}
  }

  void submit(const std::string &s) {
	std::unique_lock<std::mutex> lk(mut);
	string_queue.push(s);
	nstrings++;
	nsubmitted++;
	lk.unlock();
	data_cond.notify_all();
  }
};

#endif //CARVAIBD_REPORTER_HPP
