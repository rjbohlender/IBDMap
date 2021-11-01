/**
 * @file reporter.hpp
 * @brief Results output classes
 *
 * Created by Bohlender,Ryan James on 10/2/19.
 */

#ifndef CARVAIBD_REPORTER_HPP
#define CARVAIBD_REPORTER_HPP

#include "split.hpp"
#include "threadsafequeue.hpp"
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <fmt/include/fmt/ostream.h>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <zstr/src/zstr.hpp>

/**
 * @brief Container for the results so they can be sorted.
 *
 * Our program is multithreaded and as a result the output is not necessarily
 * ordered the same as it was input. The final step of analysis requires the
 * output be sorted, so output is sorted in a final sorting step after writing
 * out the results of every thread.
 */
struct OutContainer {
    std::string chrom;
    int pos;
    std::vector<std::vector<std::string>> data;// Two dimensional to handle multiple phenotypes in a single file
};

/**
 * @brief Thread safe shared writing class
 *
 * Our program is multithreaded and thus can run into conflicts when writing
 * with more than one thread to the same output. To solve this problem we
 * implement a writing object that has a single shared state distributed amongst
 * the threads. Using a lock and a queue, each thread can acquire the lock,
 * submit a string for printing, and return. A condition variable is used to
 * notify the print thread that a string has been submitted and the queue
 * should be checked for additional work. Once the work is complete the lock is
 * returned and the thread waits for additional work.
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
        zstr::ofstream ofs(out_path);

        while (!done || nstrings > 0) {
            data_cond.wait_for(lk, std::chrono::seconds(1), [this] { return !string_queue.empty() || done; });
            if (!string_queue.empty()) {
                s = std::move(string_queue.front());
                string_queue.pop();
                lk.unlock();
                if (!out_path.empty()) {
                    fmt::print(ofs, s);
                    ofs.flush();
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
    /**
   * @brief Contructor for the reporter
   * @param output Path to the output file.
   *
   * Should always be instantiated inside a std::make_shared call.
   */
    explicit Reporter(std::string output) : done(false), nstrings(0), nsubmitted(0), nwritten(0), out_path(std::move(output)) {
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

    /**
   * @brief Locking submission function for submitting work
   * @param s The string to be printed.
   */
    void submit(const std::string &s) {
        std::unique_lock<std::mutex> lk(mut);
        string_queue.push(s);
        nstrings++;
        nsubmitted++;
        lk.unlock();
        data_cond.notify_all();
    }

    /**
   * @brief Ensure that the output is in sorted order.
   *
   * Read through all the output, collect it, sort it by position, then
   * write to a final file. Output is initially written to the output location
   * with a .tmp extension added. All output is zipped to save space as the
   * output can be quite large for large jobs.
   */
    void sort() {
        // Ensure printing is completed and print thread is terminated before sorting.
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

        // Begin sorting
        if (!out_path.empty()) {
            zstr::ifstream reinput(out_path);
            std::string line;
            std::vector<OutContainer> sortable_output;

            int line_no = 0;
            while (std::getline(reinput, line)) {
                line_no++;
                RJBUtil::Splitter<std::string> splitter(line, " \t");
                std::vector<std::string> stats;
                if (splitter.size() > 0) {// Either the first phenotype, or a new position
                    if (sortable_output.empty() || sortable_output.back().pos != std::stoi(splitter[1])) {
                        // Same position output is guaranteed to be adjacent so we only need to check the back.
                        // We need to allocate a new output container if we've reached a position we haven't seen
                        std::copy(splitter.begin() + 2, splitter.end(), std::back_inserter(stats));
                        OutContainer oc{
                                splitter[0],
                                std::stoi(splitter[1]),
                                std::vector<std::vector<std::string>>()};
                        oc.data.emplace_back(std::move(stats));
                        sortable_output.emplace_back(std::move(oc));
                    } else {// Continuing the output
                        std::copy(splitter.begin() + 2, splitter.end(), std::back_inserter(stats));
                        sortable_output.back().data.emplace_back(std::move(stats));
                    }
                }
            }
            std::sort(sortable_output.begin(),
                      sortable_output.end(),
                      [](OutContainer &a, OutContainer &b) { return a.pos < b.pos; });
            zstr::ofstream reoutput(out_path);
            for (const auto &v : sortable_output) {
                for (const auto &w : v.data) {
                    fmt::print(reoutput, "{}\t{}\t{}\n", v.chrom, v.pos, fmt::join(w.begin(), w.end(), "\t"));
                    reoutput.flush();
                }
            }
        }
    }
};

#endif//CARVAIBD_REPORTER_HPP
