#include "ibdlib.hpp"

#include <SugarPP/include/sugarpp/range/enumerate.hpp>
#include <fmt/format.h>
#include <future>
#include <iostream>
#include <iomanip>
#include <zstr.hpp>

#define MULTI 1

auto parse_line(const std::string &line, bool new_) {
    Splitter<std::string> tok(line, " \t");
    std::string chrom_str;

    int orig_idx = 2;
    int line_start = 3;
    if (new_) {
        orig_idx = 5;
        line_start = 6;
    }
    if (boost::starts_with(tok[0], "chr")) {
        chrom_str = tok[0].substr(3, tok[0].size());
    } else {
        chrom_str = tok[0];
    }
    int chrom = std::stoi(chrom_str);
    int pos = std::stoi(tok[1]);
    double orig = std::stod(tok[orig_idx]);

    arma::rowvec vals(tok.size() - line_start, arma::fill::zeros);
    for (int i = line_start; i < tok.size(); i++) {
        vals(i - line_start) = std::stod(tok[i]);
    }
    return std::make_tuple(chrom, orig, pos, vals);
}

auto run_ibdlen(const IBDRArgs &args, GeneticMap &gmap) {
    ibd_t ibdlen;
    ibd_t ibdfrac;
    std::unordered_map<int, int> breakpoints;
    double sum_ibdlen = 0;
    std::vector<std::future<void>> futures;
    std::vector<int> chroms;
    {
        // Release GIL so we can do work multithreaded
        ScopedGILRelease scoped = ScopedGILRelease();
        build_chroms(args, chroms);

        for (auto chrom : chroms) {
            ibdlen[chrom] = {};
            breakpoints[chrom] = 0;
        }

        futures.reserve(chroms.size());
        for (auto chrom : chroms) {
#if MULTI
            futures.push_back(std::async(std::launch::async, parse_ibdlen, args, gmap, std::reference_wrapper(ibdlen), std::reference_wrapper(breakpoints), chrom));
#else
            parse_ibdlen(args, gmap, std::reference_wrapper(ibdlen), std::reference_wrapper(breakpoints), chrom);
#endif
        }

        for (auto &f : futures) {
            f.wait();
            f.get();
        }
        for (const auto &[k, v] : ibdlen) {
            for (const auto &[p, d] : v) {
                sum_ibdlen += d;
            }
        }
    }
    std::cerr << fmt::format("sum_ibdlen: {}\n", sum_ibdlen);
    for (const auto &[chrom, pos_data] : ibdlen) {
        for (const auto &[pos, len] : pos_data) {
            ibdfrac[chrom][pos] = len / sum_ibdlen;
        }
    }

    return std::make_tuple(ibdlen, ibdfrac, breakpoints);
}

auto run_avg(const IBDRArgs &args, ibd_t &ibdfrac) {
    // Release GIL so we can do work multithreaded
    ScopedGILRelease scoped = ScopedGILRelease();
    size_t arr_size = args.nruns * args.nperm;
    std::vector<ibd_t> original;
    perm_t permuted;
    std::vector<int> chroms;
    std::vector<std::future<void>> futures;

    build_chroms(args, chroms);

    for (int i = 0; i < args.phenotypes; i++) {
        original.emplace_back();
        permuted.push_back(arma::mat(chroms.size(), arr_size, arma::fill::zeros));
        for (auto chrom : chroms) {
            original.back()[chrom] = {};
        }
    }
    futures.reserve(chroms.size());
    for (auto chrom : chroms) {
#if MULTI
        futures.push_back(std::async(std::launch::async, parse_avg, args, std::reference_wrapper(ibdfrac), chrom, std::reference_wrapper(original), std::reference_wrapper(permuted)));
#else
        parse_avg(args, std::reference_wrapper(ibdfrac), chrom, std::reference_wrapper(original), std::reference_wrapper(permuted));
#endif
    }
    for (auto &f : futures) {
        f.wait();
        f.get();
    }

    std::vector<double> orig_avg;
    perm_t perm_avg;
    for (int j = 0; j < args.phenotypes; j++) {
        orig_avg.push_back(0);
        perm_avg.emplace_back(1, arr_size);
        perm_avg.back() = arma::sum(permuted[j], 0);
        // Sum over chromosomes
        for (auto i : chroms) {
            for (const auto &[k, v] : original[j][i]) {
                orig_avg.back() += v;
            }
        }
    }

    return std::make_tuple(orig_avg, perm_avg);
}

void build_chroms(const IBDRArgs &args, std::vector<int> &chroms) {
    if (args.is_single) {
        chroms.push_back(args.single);
        if (args.is_null) {
            chroms.push_back(args.null);
        }
    } else {
        for (int i = 1; i < 23; i++) {
            chroms.push_back(i);
        }
    }
}

auto run_parse(const IBDRArgs &args, const std::vector<double> &orig_avg, const perm_t &perm_avg, const std::unordered_map<int, int> &breakpoints) {
    // Release GIL so we can do work multithreaded
    // ScopedGILRelease scoped = ScopedGILRelease();
    size_t arr_size = args.nruns * args.nperm;
    std::vector<orig_t> original;
    std::vector<int> chroms;
    std::vector<arma::mat> evd;
    std::vector<arma::mat> fdr;
    std::vector<std::future<void>> futures;
    std::vector<ibd_t> deltas;

    build_chroms(args, chroms);

    int total_breakpoints = 0;
    for (const auto &[k, v] : breakpoints) {
        total_breakpoints += v;
    }

    for (int i = 0; i < args.phenotypes; i++) {
        original.emplace_back();
        evd.emplace_back(chroms.size(), arr_size);
        evd.back().fill(arma::datum::inf);
        deltas.emplace_back();
        // Allocate memory for matrix
        if (args.fdr) {
            auto fdr_ptr = new double[total_breakpoints * arr_size];
            fdr.emplace_back(fdr_ptr, total_breakpoints, arr_size, false, true);
        }
        for (auto chrom : chroms) {
            original.back()[chrom] = {};
            deltas.back()[chrom] = {};
        }
    }
    p::list fdr_list;
    if (args.fdr) {
        for (auto &v : fdr) {
            p::tuple shape = p::make_tuple(total_breakpoints, arr_size);
            np::dtype dtype = np::dtype::get_builtin<double>();
            p::tuple stride = p::make_tuple(sizeof(double), sizeof(double));
            np::ndarray fdr_arr = np::from_data(v.memptr(), dtype, shape, stride, p::object());
            fdr_list.append(fdr_arr);
        }
    }

    futures.reserve(chroms.size());
    for (auto chrom : chroms) {
#if MULTI
        futures.push_back(std::async(std::launch::async, parse, args, chrom, orig_avg, std::reference_wrapper(original), perm_avg, std::reference_wrapper(evd), std::reference_wrapper(deltas), std::reference_wrapper(fdr), std::reference_wrapper(breakpoints)));
#else
        parse(args, chrom, orig_avg, std::reference_wrapper(original), perm_avg, std::reference_wrapper(evd), std::reference_wrapper(deltas), std::reference_wrapper(fdr), std::reference_wrapper(breakpoints));
#endif
    }
    for (auto &fut : futures) {
        fut.wait();
        fut.get();
    }

    return std::make_tuple(original, evd, deltas, fdr_list);
}


p::tuple run_stages(const IBDRArgs &args, GeneticMap &gmap) {
    p::list orig_list;
    p::list delta_list;
    p::list evd_list;
    p::tuple ibdlen_ret;
    p::tuple avg_ret;
    p::tuple parse_ret;
    p::dict ibdlen_dict;

    arma::wall_clock timer;
#if 1
    timer.tic();
    auto [ibdlen, ibdfrac, breakpoints] = run_ibdlen(args, gmap);
    fmt::print(stderr, "run_ibdlen time: {}\n", timer.toc());
    timer.tic();
    auto [orig_avg, p_avg] = run_avg(args, ibdfrac);
    fmt::print(stderr, "run_avg time: {}\n", timer.toc());
    std::cerr << fmt::format("Observed averages: {}", fmt::join(orig_avg.begin(), orig_avg.end(), " ")) << std::endl;
    timer.tic();
    auto [original, evd, deltas, fdr] = run_parse(args, orig_avg, p_avg, breakpoints);
    fmt::print(stderr, "run_parse time: {}\n", timer.toc());
#else
    ibdlen_ret = run_ibdlen(args, gmap);
    ibd_t ibdfrac = p::extract<ibd_t>(ibdlen_ret[1]);
    avg_ret = run_avg(args, ibdfrac);
    std::vector<double> orig_avg = p::extract<std::vector<double>>(avg_ret[0]);
    perm_t p_avg = p::extract<perm_t>(avg_ret[1]);
    parse_ret = run_parse(args, orig_avg, p_avg);
    std::vector<orig_t> original = p::extract<std::vector<orig_t>>(parse_ret[0]);
    std::vector<arma::mat> evd = p::extract<std::vector<arma::mat>>(parse_ret[1]);
    std::vector<ibd_t> deltas = p::extract<std::vector<ibd_t>>(parse_ret[2]);
#endif

    // Convert to python objects
    for (const auto &[lno, orig] : Enumerate(original)) {
        p::dict chrom_dict;
        p::dict delta_chrom;
        for (const auto &[chrom, chrom_data] : orig) {
            p::dict pos_dict;
            p::dict delta_pos;
            p::dict ibdlen_pos;
            for (const auto &[pos, pos_data] : chrom_data) {
                p::list data;
                data.append(pos_data.first);
                data.append(pos_data.second);
                pos_dict[pos] = data;
                delta_pos[pos] = deltas[lno][chrom][pos];
                if (lno == 0)
                    ibdlen_pos[pos] = ibdlen[chrom][pos];
            }
            chrom_dict[chrom] = pos_dict;
            delta_chrom[chrom] = delta_pos;
            if (lno == 0)
                ibdlen_dict[chrom] = ibdlen_pos;
        }
        orig_list.append(chrom_dict);
        delta_list.append(delta_chrom);
        arma::rowvec evd_tmp = arma::min(evd[lno], 0);
        p::list evd_data;
        for (auto v : evd_tmp) {
            evd_data.append(v);
        }
        evd_list.append(evd_data);
    }

    std::string ret = fmt::format("p_avg size: {}", p_avg[0].size());

    return p::make_tuple(orig_list, evd_list, delta_list, ibdlen_dict, fdr);
}

void parse_ibdlen(const IBDRArgs &args, const GeneticMap &gmap, std::unordered_map<int, std::unordered_map<int, double>> &ibdlen, std::unordered_map<int, int> &breakpoints, int i) {
    breakpoints[i] = 0;

    std::string path = fmt::format(args.prefix + args.suffix, fmt::arg("i", i), fmt::arg("j", args.at));

    zstr::ifstream ifs(path);
    if (!ifs.good()) {
        throw std::runtime_error("Failed to open file path in parse_ibdlen.");
    }

    double predis = 0;
    int prechr = 0;
    int prepos = 0;
    int lineno = 0;
    double dis = 0;
    std::string line;
    std::string chr;
    while (std::getline(ifs, line)) {
        if (args.phenotypes > 1) {
            if (lineno % args.phenotypes != 0) {
                lineno++;
                continue;
            }
        }
        breakpoints[i]++;
        auto [chrom, orig, pos, vals] = parse_line(line, args.new_);
        if (chrom != prechr) {
            chr = fmt::format("chr{}", chrom);
        }
        dis = calculate_dis(gmap.find_nearest(fmt::format("chr{}", chrom), pos), pos);
        if (chrom == prechr) {
            ibdlen[prechr][prepos] = dis - predis;
        } else if (prechr != 0) {
            ibdlen[prechr][prepos] = predis + 1;
        }
        prechr = chrom;
        prepos = pos;
        predis = dis;
    }
    // For final variants
    ibdlen[prechr][prepos] = 1;
}

void parse_avg(const IBDRArgs &args, ibd_t &ibdfrac, int i, std::vector<ibd_t> &original, perm_t &permuted) {
    std::map<int, int> chrom_idx_map;
    if (args.is_single) {
        if (args.is_null) {
            if (i == args.single) {
                chrom_idx_map[i] = 0;
            } else {
                chrom_idx_map[i] = 1;
            }
        }
        chrom_idx_map[i] = 0;
    } else {
        for(int k = 0; k < 22; k++) {
            chrom_idx_map[k + 1] = k;
        }
    }
    for (int j = args.at; j < args.at + args.nruns; j++) {
        std::string fpath = fmt::format(args.prefix + args.suffix, fmt::arg("i", i), fmt::arg("j", j));

        zstr::ifstream ifs(fpath);
        int lineno = 0;
        std::string line;
        while (std::getline(ifs, line)) {

            auto [chrom, orig, pos, vals] = parse_line(line, args.new_);
            orig *= ibdfrac[chrom][pos];
            vals *= ibdfrac[chrom][pos];

            if (std::all_of(vals.begin(), vals.end(), [](const auto &val) { return val == 0; })) {
                lineno++;
                continue;
            }
            if (j == args.at) {
                original[lineno % args.phenotypes][chrom][pos] = orig;
                arma::span col_span(0, args.nperm - 1);
                permuted[lineno % args.phenotypes](chrom_idx_map[i], col_span) += vals;
            } else {
                assert(args.nperm * (j - args.at) - args.nperm * (j - args.at) == vals.n_elem);
                arma::span col_span(args.nperm * (j - args.at), args.nperm * (j - args.at));
                permuted[lineno % args.phenotypes](chrom_idx_map[i], col_span) += vals;
            }
        }
    }
}

void parse(const IBDRArgs &args, int i, const std::vector<double> &orig_avg, std::vector<orig_t> &original, const perm_t &p_avg, std::vector<arma::mat> &evd, std::vector<ibd_t> &deltas, std::vector<arma::mat> &fdr, const std::unordered_map<int, int> &breakpoints) {
    size_t arr_size = args.nruns * args.nperm;
    std::vector<std::unordered_map<int, std::unordered_map<int, std::pair<double, double>>>> orig_buffer;
    arma::rowvec evd_buffer(arr_size, arma::fill::zeros);
    std::string line;

    int row_offset = 0;
    std::map<int, int> chrom_idx_map;
    if (args.is_single) {
        if (args.is_null) {
            if (i == args.single) {
                chrom_idx_map[i] = 0;
            } else {
                row_offset += breakpoints.at(args.single);
                chrom_idx_map[i] = 1;
            }
        }
        chrom_idx_map[i] = 0;
    } else {
        for(int k = 0; k < 22; k++) {
            // Offset for FDR matrix start row
            if (k + 1 < i) {
                row_offset += breakpoints.at(k + 1);
            }
            chrom_idx_map[k + 1] = k;
        }
    }

    orig_buffer.reserve(args.phenotypes);
    for (int k = 0; k < args.phenotypes; k++) {
        orig_buffer.emplace_back();
        orig_buffer.back()[i] = {};
    }

    std::vector<std::unique_ptr<zstr::ifstream>> ifs_buffer;
    for (int j = args.at; j < args.at + args.nruns; j++) {
        ifs_buffer.emplace_back(std::make_unique<zstr::ifstream>(fmt::format(args.prefix + args.suffix, fmt::arg("i", i), fmt::arg("j", j))));
    }

    int lineno = -1;
    int breakpointno = -1;
    while(true) {
        lineno++;
        if (lineno % args.phenotypes == 0) {
            breakpointno++;
        }
        for (auto [j, f] : Enumerate(ifs_buffer)) {
            if (!std::getline(*f, line)) {
                continue;
            }
            auto [chrom, orig, pos, vals] = parse_line(line, args.new_);
            if (std::all_of(vals.begin(), vals.end(), [](auto v) { return v == 0; })) {
                continue;
            }
            orig -= orig_avg[lineno % args.phenotypes];
            arma::span col_span(args.nperm * j, args.nperm * (j + 1) - 1);
            arma::rowvec avg_vec = p_avg[lineno % args.phenotypes](arma::span(), col_span);
            vals -= p_avg[lineno % args.phenotypes](arma::span(), col_span);

            deltas[lineno % args.phenotypes][chrom][pos] = orig;
            if (orig_buffer[lineno % args.phenotypes][chrom].count(pos) == 0) {
                orig_buffer[lineno % args.phenotypes][chrom][pos] = {0, 0};
            }
            double count = static_cast<double>(arma::accu(vals >= orig));
            orig_buffer[lineno % args.phenotypes][chrom][pos].first += count;
            orig_buffer[lineno % args.phenotypes][chrom][pos].second += static_cast<double>(vals.size());
            evd_buffer(col_span) = vals;
            if (args.fdr) {
                fdr[lineno % args.phenotypes](row_offset + breakpointno, col_span) = vals;
            }
        }
        if (line.empty()) {
            break;
        }
        evd_buffer = rank(evd_buffer, "descend") / (static_cast<double>(arr_size) + 1.);
        evd[lineno % args.phenotypes].row(i - 1) = minimum(evd[lineno % args.phenotypes].row(i - 1), evd_buffer);
    }

    for (const auto &[lno, phenotype] : Enumerate(orig_buffer)) {
        for (const auto &[chrom, chrom_data] : phenotype) {
            for (const auto &[pos, data] : chrom_data) {
                original[lno][chrom][pos] = std::make_pair((data.first + 1.) / (data.second + 1.), data.first);
            }
        }
    }
}

double
calculate_dis(std::pair<std::pair<int, double>, std::pair<int, double>> bounds,
              int pos) {
    if (bounds.first == bounds.second) {
        return bounds.first.second;
    } else {
        int pos1 = bounds.first.first;
        int pos2 = bounds.second.first;
        double cm1 = bounds.first.second;
        double cm2 = bounds.second.second;
        double dis = cm1 + (double(pos) - double(pos1)) / (double(pos2) - double(pos1)) * (cm2 - cm1);
        return dis;
    }
}

arma::rowvec rank(arma::rowvec &v, const char *direction) {
    if (strcmp(direction, "ascend") != 0 && strcmp(direction, "descend") != 0)
        throw (std::logic_error("Order argument for rank() must be either 'ascend' or 'descend'"));

    arma::uvec sort_indices;
    try {
        sort_indices = arma::sort_index(v, direction);
    } catch (const std::logic_error &e) {
        std::cerr << "NANs among ranked values. Replacing with 0.\n";
        v.replace(arma::datum::nan, 0);
        sort_indices = arma::sort_index(v, direction);
    }
    arma::vec sorted = v(sort_indices);

    arma::rowvec ranks = arma::rowvec(v.n_elem, arma::fill::zeros);
    arma::sword i = 0, j = 0;

    while (i < v.n_elem) {
        j = i + 1;
        // Find the next different value
        while (j < v.n_elem) {
            if (sorted(i) != sorted(j))
                break;
            j++;
        }
        // Adjusted rank
        for (arma::uword k = i; k < j; k++) {
            ranks(sort_indices(k)) = 1. + (i + j - 1.) / 2.0f;
        }
        // Update i
        i = j;
    }

    return ranks;
}

arma::rowvec minimum(const arma::rowvec &first, const arma::rowvec &second) {
    assert(first.n_elem == second.n_elem);
    arma::rowvec ret(first.n_elem, arma::fill::zeros);

    for (arma::uword i = 0; i < ret.n_elem; i++) {
        ret(i) = std::min(first(i), second(i));
    }
    return ret;
}

void print_vec(arma::rowvec &vec) {
    std::ostringstream oss;
    for (auto [i, v] : Enumerate(vec)) {
        oss << std::setprecision(4) << v;
        if (i == vec.n_elem - 1) {
            oss << "\n\n";
        } else if(i > 0 && (i + 1) % 4 == 0) {
            oss << "\n";
        } else {
            oss << " ";
        }
    }
    std::cerr << oss.str();
}

void print_vec(arma::rowvec &&vec) {
    std::ostringstream oss;
    for (auto [i, v] : Enumerate(vec)) {
        oss << std::setprecision(4) << v;
        if (i == vec.n_elem - 1) {
            oss << "\n\n";
        } else if(i > 0 && (i + 1) % 4 == 0) {
            oss << "\n";
        } else {
            oss << " ";
        }
    }
    std::cerr << oss.str();
}
