#include "../binary_format/ibd_reader.hpp"
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <file.ibdf>\n", argv[0]);
        return 1;
    }

    IbdfReader reader(argv[1]);
    uint64_t n = reader.n_positions();

    fprintf(stderr, "Positions: %lu, Samples: %u, Checkpoint interval: %u\n",
            n, reader.n_samples(), reader.checkpoint_interval());

    std::vector<size_t> sizes;

    for (uint64_t i = 0; i < n; ++i) {
        if (!reader.index()[i].is_checkpoint()) continue;

        auto cp = ibdf::decode_checkpoint(reader.block_ptr(i));
        size_t sz = cp.n_pairs();
        sizes.push_back(sz);

        printf("%llu\t%llu\t%zu\n", i, reader.index()[i].bp_pos, sz);
    }

    // Summary stats
    auto [mn, mx] = std::minmax_element(sizes.begin(), sizes.end());
    double mean = std::accumulate(sizes.begin(), sizes.end(), 0.0) / sizes.size();

    std::vector<size_t> sorted = sizes;
    std::sort(sorted.begin(), sorted.end());
    size_t median = sorted[sorted.size() / 2];
    size_t p95 = sorted[static_cast<size_t>(sorted.size() * 0.95)];
    size_t p99 = sorted[static_cast<size_t>(sorted.size() * 0.99)];

    fprintf(stderr, "\n--- Active set size summary ---\n");
    fprintf(stderr, "  Min:    %zu\n", *mn);
    fprintf(stderr, "  Median: %zu\n", median);
    fprintf(stderr, "  Mean:   %.1f\n", mean);
    fprintf(stderr, "  P95:    %zu\n", p95);
    fprintf(stderr, "  P99:    %zu\n", p99);
    fprintf(stderr, "  Max:    %zu\n", *mx);

    return 0;
}
