#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
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
    void initialize(const std::vector<double>& actin_x, const std::vector<double>& actin_y, const std::vector<double>& actin_z,
                    const std::vector<double>& myosin_x, const std::vector<double>& myosin_y, const std::vector<double>& myosin_z);

    // Rebuild the neighbor list using a cell list approach.
    void rebuild_neighbor_list();

    // Set current positions for each species.
    void set_species_positions(const std::vector<double>& actin_x, const std::vector<double>& actin_y, const std::vector<double>& actin_z,
                               const std::vector<double>& myosin_x, const std::vector<double>& myosin_y, const std::vector<double>& myosin_z);

    // Check if the neighbor list needs to be rebuilt (based on displacement).
    bool needs_rebuild() const;

      // Step 0: Clear and resize internal containers.
    void clearAndResizeContainers();

    // Step 1: Populate the cell list with particle indices.
    void populateCellList();

    // Step 2: Compute neighbor pairs in parallel.
    void computeNeighbors(std::vector<std::vector<std::pair<size_t, size_t>>> &tls);

    // Step 3: Merge thread-local pairs into the main neighbor list.
    void mergeThreadLocalPairs(const std::vector<std::vector<std::pair<size_t, size_t>>> &tls);

    // Step 4: Update the last known positions.
    void updateLastKnownPositions();
    // Return a pair of vectors (first: actin neighbors, second: myosin neighbors) for a given particle.
    std::pair<std::vector<int>, std::vector<int>> get_neighbors_by_type(int index) const;

    // Access the neighbor list for a specific particle.
    const std::vector<std::pair<int, ParticleType>>& get_neighbors(int index) const;

private:
    // Helper functions.
    std::tuple<int, int, int> get_cell_index(double x, double y, double z) const;
    double displacement(double cx, double cy, double cz,
                        double lx, double ly, double lz) const;
    void concatenate_positions();
    void track_species_types();
    // Data members.
    double cutoff_radius_;
    double threshold_;
    std::vector<double> box_;
    int num_cells_x_, num_cells_y_, num_cells_z_;
    double cell_size_x_, cell_size_y_, cell_size_z_;
    std::vector<std::tuple<int, int, int>> neighboring_cells;

    // Structure-of-arrays storage for positions.
    std::vector<double> actin_x_, actin_y_, actin_z_;
    std::vector<double> myosin_x_, myosin_y_, myosin_z_;
    size_t n_actins_;
    std::vector<double> last_actin_x_, last_actin_y_, last_actin_z_;
    std::vector<double> last_myosin_x_, last_myosin_y_, last_myosin_z_;

    std::vector<double> all_x_, all_y_, all_z_;
    std::vector<ParticleType> species_types_;
    std::vector<std::vector<std::pair<int, ParticleType>>> neighbor_list_;

    std::vector<std::vector<int>> cell_list_;
};

#endif // NEIGHBORLIST_H
