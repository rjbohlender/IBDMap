// -*- mode: C++ -*-
/**
 * @file split.hpp
 * @brief Implements a simple (read amateur) Splitter.
 *
 */
#ifndef PGEN_SPLIT_HPP
#define PGEN_SPLIT_HPP 1

#include <algorithm>
#include <string>
#include <vector>
#include <sstream>

namespace RJBUtil {
    /**
 * @class Splitter split.hpp "src/util/split.hpp"
 * @brief A utility for quickly splitting strings into readable substring segments.
 *
 * @remark  Removed string_view, as string_views can't be built from iterators.
 */
    template<typename String_t>
    class Splitter {
    public:
        /** Aliases */
        using string_type = String_t;
        using const_iter = typename string_type::const_iterator;
        using size_type = typename String_t::size_type;

        /** Constructors */
        Splitter() = default;
        Splitter(const string_type &str, const string_type &delim, bool allow_adjacent = false)
            : data_(str), delim_(delim), allow_adjacent_(allow_adjacent) {
            split();
        }

        Splitter(const string_type &&str, const string_type &&delim, bool allow_adjacent = false)
            : data_(str), delim_(delim), allow_adjacent_(allow_adjacent) {
            split();
        }

        Splitter(const string_type &str, const string_type &&delim, bool allow_adjacent = false)
            : data_(str), delim_(delim), allow_adjacent_(allow_adjacent) {
            split();
        }

        Splitter(const string_type &&str, const string_type &delim, bool allow_adjacent = false)
            : data_(str), delim_(delim), allow_adjacent_(allow_adjacent) {
            split();
        }

        Splitter &operator=(const Splitter &rhs) {
            data_ = rhs.data_;
            delim_ = rhs.delim_;
            tokens_ = rhs.tokens_;
            allow_adjacent_ = rhs.allow_adjacent_;

            return *this;
        }

        Splitter &operator=(Splitter &&rhs) noexcept {
            data_ = std::move(rhs.data_);
            delim_ = std::move(rhs.delim_);
            tokens_ = std::move(rhs.tokens_);
            allow_adjacent_ = std::move(rhs.allow_adjacent_);

            return *this;
        }

        auto join(const string_type &delim) {
            std::stringstream os;
            os << tokens_[0];
            for (int i = 1; i < tokens_.size(); i++) {
                os << delim << tokens_[i];
            }
            return os.str();
        }

        auto begin() {
            return tokens_.begin();
        }

        auto cbegin() {
            return tokens_.cbegin();
        }

        auto end() {
            return tokens_.end();
        }

        auto cend() {
            return tokens_.cend();
        }

        auto size() {
            return tokens_.size();
        }

        auto empty() {
            return tokens_.empty();
        }

        auto erase(std::vector<std::string>::iterator &pos) {
            tokens_.erase(pos);
        }

        auto erase(std::vector<std::string>::const_iterator &pos) {
            tokens_.erase(pos);
        }

        auto front() {
            return tokens_.front();
        }

        auto back() {
            return tokens_.back();
        }

        template<typename _Integer>
        auto operator[](_Integer i) {
            return tokens_[i];
        }

        template<typename _Integer>
        auto at(_Integer i) {
            return tokens_.at(i);
        }

    private:
        /**
   * @brief A split member function to simplify multiple constructors should we need them.
   */
        void split() {
#if 0
	const_iter cit = data_.cbegin();
	size_type cur = 0;
	size_type next = 0;
	size_type npos = std::string::npos;
	while ((cit + cur) != data_.cend()) {
	  next = data_.find_first_of(delim_, cur);
	  // Skip adjacent.
	  if (next - cur == 0) {
		cur++;
		next = cur;
		continue;
	  }
	  if (next == npos) {
		if (data_.cend() - (cit + cur) > 0) {
		  tokens_.emplace_back(std::string(cit + cur, data_.cend()));
		}
		break;
	  }

	  tokens_.emplace_back(std::string(cit + cur, cit + next));

	  // Prepare for next loop
	  next++;
	  cur = next;

	}
#else
            // Faster version described at: https://www.bfilipek.com/2018/07/string-view-perf-followup.html
            for (auto first = data_.data(), second = data_.data(), last = first + data_.size(); second != last && first != last;
                 first = second + 1) {
                second = std::find_first_of(first, last, std::cbegin(delim_), std::cend(delim_));

                if (first != second || allow_adjacent_)
                    tokens_.emplace_back(first, second - first);
            }
#endif
        }

        /** Private members */
        string_type data_;
        string_type delim_;
        std::vector<std::string> tokens_;
        bool allow_adjacent_ = false;
    };
}// namespace RJBUtil
#endif// PGEN_SPLIT_HPP