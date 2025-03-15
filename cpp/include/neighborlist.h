#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>
#include <unordered_map>
#include <cmath>
#include <utility>
#include <set>
#include <algorithm>
#include <omp.h>

#ifndef UTILS_H
#include "utils.h"
#endif
using vec = utils::vec;



// Enum class for particle species.
enum class ParticleType {
    Actin,
    Myosin
};

class NeighborList {
public:
    // Constructors.
    NeighborList(double cutoff_radius, const std::vector<double>& box, double threshold);
    NeighborList();

    // Initialize the neighbor list with positions for actin and myosin.
    void initialize(const std::vector<vec>& actin_positions, const std::vector<vec>& myosin_positions);

    // Rebuild the neighbor list using a cell list approach.
    void rebuild_neighbor_list();

    // Set current positions for each species.
    void set_species_positions(const std::vector<vec>& actin_positions, const std::vector<vec>& myosin_positions);

    // Check if the neighbor list needs to be rebuilt (based on displacement).
    bool needs_rebuild() const;

      // Step 0: Clear and resize internal containers.
    void clearAndResizeContainers();

    // Step 1: Populate the cell list with particle indices.
    void populateCellList();

    // Step 2: Create thread-local storage for neighbor lists.
    // The TLS container is passed by reference and resized appropriately.
    void createThreadLocalStorage(std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> &tls) const;

    // Step 3: Compute neighbors in parallel.
    void computeNeighbors(std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> &tls);

    // Step 4: Merge thread-local neighbor lists into the main neighbor list.
    void mergeThreadLocalLists(const std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> &tls);

    // Step 5: Update the last known positions.
    void updateLastKnownPositions();
    // Return a pair of vectors (first: actin neighbors, second: myosin neighbors) for a given particle.
    std::pair<std::vector<int>, std::vector<int>> get_neighbors_by_type(int index) const;

    // Access the neighbor list for a specific particle.
    const std::vector<std::pair<int, ParticleType>>& get_neighbors(int index) const;

private:
    // Helper functions.
    std::tuple<int, int, int> get_cell_index(const vec& position) const;
    double displacement(const vec& current, const vec& last) const;
    void concatenate_positions();
    void track_species_types();
    // Data members.
    double cutoff_radius_;
    double threshold_;
    std::vector<double> box_;
    int num_cells_x_, num_cells_y_, num_cells_z_;
    double cell_size_x_, cell_size_y_, cell_size_z_;
    std::vector<std::tuple<int, int, int>> neighboring_cells;

    std::vector<vec> actin_positions_;
    std::vector<vec> myosin_positions_;
    size_t n_actins_;
    std::vector<vec> last_actin_positions_;
    std::vector<vec> last_myosin_positions_;

    std::vector<vec> all_positions_;
    std::vector<ParticleType> species_types_;
    std::vector<std::vector<std::pair<int, ParticleType>>> neighbor_list_;

    struct pair_hash {
    // Hash function for std::pair
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        return std::hash<T1>{}(pair.first) ^ (std::hash<T2>{}(pair.second) << 1);
    }

    // Hash function for std::tuple<int, int, int>
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        std::size_t h1 = std::hash<int>{}(std::get<0>(t));
        std::size_t h2 = std::hash<int>{}(std::get<1>(t));
        std::size_t h3 = std::hash<int>{}(std::get<2>(t));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
    };

    std::unordered_map<std::tuple<int, int, int>, std::vector<int>, pair_hash> cell_list_;
};

#endif // NEIGHBORLIST_H