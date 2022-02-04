#ifndef IBDLIB_LIBRARY_H
#define IBDLIB_LIBRARY_H

#include "ibdlibclasses.hpp"
#include "geneticmap.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <tuple>
#include <unordered_map>
#include <armadillo>

namespace p = boost::python;
namespace np = boost::python::numpy;

using ibd_t = std::unordered_map<int, std::unordered_map<int, double>>;
using orig_t = std::unordered_map<int, std::unordered_map<int, std::pair<double, double>>>;
using bp_t = std::unordered_map<int, int>;
using perm_t = std::vector<arma::mat>;

arma::rowvec rank(arma::rowvec &v, const char *direction);

arma::rowvec minimum(const arma::rowvec &first, const arma::rowvec &second);


/** Encapsulate stages so that we don't have to keep converting to python objects.
 *
 * @param args Program arguments
 * @param gmap Genetic map
 * @return Final stage results.
 */
p::tuple run_stages(const IBDRArgs &args, GeneticMap &gmap);

auto parse_line(const std::string &line, bool new_ = false);
void parse_ibdlen(const IBDRArgs &args, const GeneticMap &gmap, ibd_t &ibdlen, bp_t &breakpoints, int i);
void parse_avg(const IBDRArgs &val, ibd_t &ibdfrac, int i, std::vector<ibd_t> &original, perm_t &permuted);
void parse(const IBDRArgs &args, int i, const std::vector<double> &orig_avg, std::vector<orig_t> &original, const perm_t &p_avg, std::vector<arma::mat> &evd, std::vector<ibd_t> &deltas, std::vector<arma::mat> &fdr, const std::unordered_map<int, int> &breakpoints);
auto run_ibdlen(const IBDRArgs &args, GeneticMap &gmap);
auto run_avg(const IBDRArgs &args, ibd_t &ibdfrac);
auto run_parse(const IBDRArgs &args, const std::vector<double> &orig_avg, const perm_t &perm_avg, const std::unordered_map<int, int> &breakpoints);
void build_chroms(const IBDRArgs &args, std::vector<int> &chroms);

void print_vec(arma::rowvec &vec);
void print_vec(arma::rowvec &&vec);

double
calculate_dis(std::pair<std::pair<int, double>, std::pair<int, double>> bounds,
              int pos);

typedef std::shared_ptr<GeneticMap> gmap_ptr;
typedef std::shared_ptr<IBDRArgs> ibdr_ptr;

BOOST_PYTHON_MODULE(ibdlib) {
    np::initialize();
    p::class_<GeneticMap, gmap_ptr>("GeneticMap", p::init<std::string>())
            .def(p::init<p::list>())
            .def("find_nearest", &GeneticMap::find_nearest);
    p::class_<IBDRArgs, ibdr_ptr>("IBDRArgs")
            .def_readwrite("at", &IBDRArgs::at)
            .def_readwrite("nperm", &IBDRArgs::nperm)
            .def_readwrite("nruns", &IBDRArgs::nruns)
            .def_readwrite("is_single", &IBDRArgs::is_single)
            .def_readwrite("is_null", &IBDRArgs::is_null)
            .def_readwrite("single", &IBDRArgs::single)
            .def_readwrite("null", &IBDRArgs::null)
            .def_readwrite("new_", &IBDRArgs::new_)
            .def_readwrite("fdr", &IBDRArgs::fdr)
            .def_readwrite("prefix", &IBDRArgs::prefix)
            .def_readwrite("suffix", &IBDRArgs::suffix)
            .def_readwrite("gmap_path", &IBDRArgs::gmap_path)
            .def_readwrite("phenotypes", &IBDRArgs::phenotypes);
    p::def("run_stages", run_stages);
}
#endif//IBDLIB_LIBRARY_H
