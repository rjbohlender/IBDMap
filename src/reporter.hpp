//
// Created by Bohlender,Ryan James on 10/2/19.
//

#ifndef CARVAIBD_REPORTER_HPP
#define CARVAIBD_REPORTER_HPP

#include <mutex>
#include <iostream>
#include <thread>
#include <queue>
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
            data_cond.wait(lk, [this] { return !string_queue.empty() || done; });
            if (!string_queue.empty()) {
                s = string_queue.front();
                string_queue.pop();
                lk.unlock();
                if (!out_path.empty()) {
		            ofs << s << std::flush;
                } else {
                    std::cout << s << std::flush;
                }
                nstrings--;
                nwritten++;
                lk.lock();
            }
        }
    }

public:
    Reporter(std::string &output) :
            done(false), nstrings(0), nsubmitted(0), nwritten(0), out_path(output) {
        print_thread = std::thread(&Reporter::print, std::ref(*this));
    }

    ~Reporter() {
        while (!string_queue.empty() || nstrings > 0) {
            std::cerr << "queue size: " << string_queue.size() << " nstrings: " << nstrings << std::endl;
            std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
        }
        done = true;
        data_cond.notify_all();
        if (print_thread.joinable()) {
            std::cerr << "In reporter:" << std::endl;
            std::cerr << "nwritten: " << nwritten << std::endl;
            std::cerr << "nsubmitted: " << nsubmitted << std::endl;
            print_thread.join();
        }
    }

    void submit(const std::string &s) {
        std::unique_lock<std::mutex> lk(mut);
        string_queue.push(s);
        nstrings++;
        nsubmitted++;
        data_cond.notify_all();
    }
};

#endif //CARVAIBD_REPORTER_HPP
