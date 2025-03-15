#include <benchmark/benchmark.h>
#include <vector>
#include <memory>
#include <random>
#include "neighborlist.h"  // Contains the refactored NeighborList and vec definitions

// Fixture for benchmarking NeighborList
class NeighborListBenchmark : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State& state) override {
        // Initialize NeighborList with typical parameters.
        double cutoff_radius = 1.0;
        std::vector<double> box = {10.0, 10.0, 10.0};
        double threshold = 0.5;
        neighborlist = std::make_unique<NeighborList>(cutoff_radius, box, threshold);

        // Generate random positions for actins and myosins based on provided arguments.
        size_t num_actins = static_cast<size_t>(state.range(0));
        size_t num_myosins = static_cast<size_t>(state.range(1));
        actin_positions = GenerateRandomPositions(num_actins);
        myosin_positions = GenerateRandomPositions(num_myosins);

        // Initialize the neighbor list with these positions.
        neighborlist->initialize(actin_positions, myosin_positions);
    }

    void TearDown(const ::benchmark::State& state) override {
        // Any cleanup if necessary.
    }

protected:
    std::unique_ptr<NeighborList> neighborlist;
    std::vector<vec> actin_positions;
    std::vector<vec> myosin_positions;

private:
    std::vector<vec> GenerateRandomPositions(size_t count) {
        std::vector<vec> positions;
        positions.reserve(count);
        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> dist(0.0, 10.0);
        for (size_t i = 0; i < count; ++i) {
            positions.push_back({dist(rng), dist(rng), dist(rng)});
        }
        return positions;
    }
};

// Benchmark for Step 0: Clearing and resizing internal containers.
BENCHMARK_DEFINE_F(NeighborListBenchmark, ClearAndResizeContainers)(benchmark::State& state) {
    for (auto _ : state) {
        neighborlist->clearAndResizeContainers();
    }
}
BENCHMARK_REGISTER_F(NeighborListBenchmark, ClearAndResizeContainers)
    ->Args({100, 50})
    ->Args({200, 100});

// Benchmark for Step 1: Populating the cell list.
BENCHMARK_DEFINE_F(NeighborListBenchmark, PopulateCellList)(benchmark::State& state) {
    // Ensure containers are cleared first.
    neighborlist->clearAndResizeContainers();
    for (auto _ : state) {
        neighborlist->populateCellList();
    }
}
BENCHMARK_REGISTER_F(NeighborListBenchmark, PopulateCellList)
    ->Args({100, 50})
    ->Args({200, 100});

// Benchmark for Step 2: Creating thread-local storage.
BENCHMARK_DEFINE_F(NeighborListBenchmark, CreateThreadLocalStorage)(benchmark::State& state) {
    for (auto _ : state) {
        std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> tls;
        neighborlist->createThreadLocalStorage(tls);
        benchmark::DoNotOptimize(tls);
    }
}
BENCHMARK_REGISTER_F(NeighborListBenchmark, CreateThreadLocalStorage)
    ->Args({100, 50})
    ->Args({200, 100});

// Benchmark for Step 3: Computing neighbors in parallel.
BENCHMARK_DEFINE_F(NeighborListBenchmark, ComputeNeighbors)(benchmark::State& state) {
    std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> tls;
    neighborlist->createThreadLocalStorage(tls);
    for (auto _ : state) {
        neighborlist->computeNeighbors(tls);
    }
}
BENCHMARK_REGISTER_F(NeighborListBenchmark, ComputeNeighbors)
    ->Args({100, 50})
    ->Args({200, 100})
    ->Args({400, 200})
    ->Args({800, 400});


// Benchmark for Step 4: Merging thread-local neighbor lists.
BENCHMARK_DEFINE_F(NeighborListBenchmark, MergeThreadLocalLists)(benchmark::State& state) {
    std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> tls;
    neighborlist->createThreadLocalStorage(tls);
    neighborlist->computeNeighbors(tls);
    // Clear the main neighbor list to simulate merging.
    neighborlist->clearAndResizeContainers();
    for (auto _ : state) {
        neighborlist->mergeThreadLocalLists(tls);
    }
}
BENCHMARK_REGISTER_F(NeighborListBenchmark, MergeThreadLocalLists)
    ->Args({100, 50})
    ->Args({200, 100});

// Benchmark for Step 5: Updating last known positions.
BENCHMARK_DEFINE_F(NeighborListBenchmark, UpdateLastKnownPositions)(benchmark::State& state) {
    for (auto _ : state) {
        neighborlist->updateLastKnownPositions();
    }
}
BENCHMARK_REGISTER_F(NeighborListBenchmark, UpdateLastKnownPositions)
    ->Args({100, 50})
    ->Args({200, 100});

BENCHMARK_MAIN();
