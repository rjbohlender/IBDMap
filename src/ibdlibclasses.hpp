//
// Created by Bohlender,Ryan J on 11/23/21.
//

#ifndef CARVAIBD_IBDLIBCLASSES_HPP
#define CARVAIBD_IBDLIBCLASSES_HPP

#include <boost/python.hpp>
#include <string>
#include <optional>

struct IBDRArgs {
    int at;
    int phenotypes;
    int nperm;
    int nruns;
    bool is_single;
    bool is_null;
    int single;
    int null;
    bool new_;
    bool fdr;
    std::string prefix;
    std::string suffix;
    std::string gmap_path;
};

class ScopedGILRelease {
public:
    inline ScopedGILRelease() {
        m_thread_state = PyEval_SaveThread();
    }

    inline ~ScopedGILRelease() {
        PyEval_RestoreThread(m_thread_state);
    }

private:
    PyThreadState *m_thread_state;
};

#endif //CARVAIBD_IBDLIBCLASSES_HPP
