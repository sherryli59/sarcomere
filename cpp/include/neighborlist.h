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



enum class ParticleType {
    Actin,
    Myosin
};

// Define a NeighborList class that manages neighbors with a unified structure
class NeighborList {
public:
    NeighborList(double cutoff_radius, const std::vector<double>& box, double threshold)
        : cutoff_radius_(cutoff_radius), threshold_(threshold), box_(box) {
        num_cells_x_ = static_cast<int>(std::floor(box[0] / cutoff_radius));
        num_cells_y_ = static_cast<int>(std::floor(box[1] / cutoff_radius));

        // Adjust cell size to ensure full coverage
        cell_size_x_ = box[0] / num_cells_x_;
        cell_size_y_ = box[1] / num_cells_y_;

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                if (num_cells_x_ > 2 && num_cells_y_ > 2) {
                    neighboring_cells.emplace_back(dx, dy);
                    continue;
                }
                if (dx == 0 && dy == 0) {
                    neighboring_cells.emplace_back(dx, dy);
                    continue;
                }
                // If num_cells_x_ == 1, don't add any neighbors in x-direction
                if (num_cells_x_ == 1 && dx != 0) continue;
                // If num_cells_y_ == 1, don't add any neighbors in y-direction
                if (num_cells_y_ == 1 && dy != 0) continue;
                // If num_cells_x_ == 2, only include dx = 1
                if (num_cells_x_ == 2 && dx == -1) continue;
                // If num_cells_y_ == 2, only include dy = 1
                if (num_cells_y_ == 2 && dy == -1) continue;
                // Default case: Include all normal neighbors
                neighboring_cells.emplace_back(dx, dy);
            }
        }
    }

    NeighborList() {}

    // Initialize with separate positions for each species
    void initialize(const std::vector<vec>& actin_positions, const std::vector<vec>& myosin_positions) {

        actin_positions_ = actin_positions;
        n_actins_ = actin_positions.size();
        myosin_positions_ = myosin_positions;

        // Initialize last known positions for displacement tracking
        last_actin_positions_ = actin_positions;
        last_myosin_positions_ = myosin_positions;

        // Concatenate positions into a single vector for the neighbor list
        concatenate_positions();

        // Track species types for each concatenated particle
        track_species_types();

        // Resize neighbor list to match all_positions_ size
        neighbor_list_.resize(all_positions_.size());

        rebuild_neighbor_list();
    }

    // Rebuild the neighbor list using a cell list approach
    void rebuild_neighbor_list() {
        printf("Rebuilding neighbor list\n");
        // Clear previous neighbor lists and cell list
        neighbor_list_.clear();
        neighbor_list_.resize(all_positions_.size());
        cell_list_.clear();
        // Populate the cell list with particle indices
        for (size_t i = 0; i < all_positions_.size(); ++i) {
            std::pair<int, int> cell = get_cell_index(all_positions_[i]);
            cell_list_[cell].push_back(i);
        }

        // Thread-local storage for neighbor lists
        std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> thread_local_neighbor_list(omp_get_max_threads());
        for (auto& local_list : thread_local_neighbor_list) {
            local_list.resize(all_positions_.size());
        }

        // Parallel loop to compute neighbors
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();

            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < all_positions_.size(); ++i) {
                vec position = all_positions_[i];
                std::pair<int, int> cell = get_cell_index(position);
                
                // Check current cell and precomputed neighboring cells
                for (const auto& offset : neighboring_cells) {
                    int neighbor_x = (cell.first + offset.first + num_cells_x_) % num_cells_x_;
                    int neighbor_y = (cell.second + offset.second + num_cells_y_) % num_cells_y_;
                    std::pair<int, int> neighbor_cell = {neighbor_x, neighbor_y};
                    // Skip if the neighbor cell is empty
                    if (cell_list_.count(neighbor_cell) == 0) continue;

                    // Iterate over particles in the neighboring cell
                    for (int j : cell_list_[neighbor_cell]) {
                        if (i >= j) continue;
                        // Compute the distance and check if it's within the cutoff radius
                        double distance = all_positions_[i].distance(all_positions_[j], box_);
                        if (distance < cutoff_radius_) {
                            // Update thread-local neighbor list
                            thread_local_neighbor_list[thread_id][i].emplace_back(j, species_types_[j]);
                            thread_local_neighbor_list[thread_id][j].emplace_back(i, species_types_[i]);
                        }
                    }
                }
            }
        }

        // Merge thread-local neighbor lists into the global neighbor list
        for (const auto& local_list : thread_local_neighbor_list) {
            for (size_t i = 0; i < all_positions_.size(); ++i) {
                neighbor_list_[i].insert(
                    neighbor_list_[i].end(),
                    local_list[i].begin(),
                    local_list[i].end()
                );
            }
        }


        
        // Update last known positions after rebuilding the neighbor list
        last_actin_positions_ = actin_positions_;
        last_myosin_positions_ = myosin_positions_;
    }

    // Set current particle positions separately for each species
    void set_species_positions(const std::vector<vec>& actin_positions, const std::vector<vec>& myosin_positions) {
        actin_positions_ = actin_positions;
        n_actins_ = actin_positions.size();
        myosin_positions_ = myosin_positions;
        concatenate_positions();
        track_species_types();
    }

    // Check if neighbor list needs to be rebuilt based on particle displacement
    bool needs_rebuild() const {
        bool needs_rebuild = false;
        // Check each actin particle’s displacement
        #pragma omp parallel for
        for (size_t i = 0; i < actin_positions_.size(); ++i) {
            if (displacement(actin_positions_[i], last_actin_positions_[i]) > threshold_) {
                needs_rebuild = true;
                #pragma omp cancel for
            }
        }
        #pragma omp barrier 
        if (needs_rebuild) return true;
        // Check each myosin particle’s displacement
        #pragma omp parallel for
        for (size_t i = 0; i < myosin_positions_.size(); ++i) {
            if (displacement(myosin_positions_[i], last_myosin_positions_[i]) > threshold_) {
                needs_rebuild = true;
                #pragma omp cancel for
            }
        }
        #pragma omp barrier 
        return false;
    }
    // Access the neighbor list for a specific particle
    const std::vector<std::pair<int, ParticleType>>& get_neighbors(int index) const {
        return neighbor_list_[index];
    }

    // Access neighbors of a specific particle, separated by type
    std::pair<std::vector<int>, std::vector<int>> get_neighbors_by_type(int index) const {
        std::vector<int> actin_neighbors;
        std::vector<int> myosin_neighbors;

        for (const auto& [neighbor_index, type] : neighbor_list_[index]) {

            if (type == ParticleType::Actin) {
                actin_neighbors.push_back(neighbor_index);
            } else if (type == ParticleType::Myosin) {
                myosin_neighbors.push_back(neighbor_index-n_actins_);
            }
        }
        return {actin_neighbors, myosin_neighbors};
    }

private:
    double cutoff_radius_;
    double threshold_;  // Skin distance threshold to determine neighbor list rebuild
    double cell_size_x_, cell_size_y_;
    int num_cells_x_, num_cells_y_, n_actins_;
    std::vector<double> box_;
    std::vector<std::pair<int, int>> neighboring_cells; 

    std::vector<vec> actin_positions_;
    std::vector<vec> myosin_positions_;
    std::vector<vec> last_actin_positions_; // Previous actin positions for displacement tracking
    std::vector<vec> last_myosin_positions_; // Previous myosin positions for displacement tracking
    std::vector<vec> all_positions_;
    std::vector<ParticleType> species_types_;

    // Unified neighbor list storing pairs of (index, type) for each particle
    std::vector<std::vector<std::pair<int, ParticleType>>> neighbor_list_;

    // Hash function for std::pair<int, int> to use with unordered_map
    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator()(const std::pair<T1, T2>& pair) const {
            return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
        }
    };

    // Cell list to hold particle indices by cell
    std::unordered_map<std::pair<int, int>, std::vector<int>, pair_hash> cell_list_;

    std::pair<int, int> get_cell_index(const vec& position) const {
        int cell_x = static_cast<int>(std::floor((position.x+box_[0]) / cell_size_x_)) % num_cells_x_;
        int cell_y = static_cast<int>(std::floor((position.y+box_[1]) / cell_size_y_))% num_cells_y_;

        // Correct modulo operation to handle negative values
        cell_x = (cell_x % num_cells_x_ + num_cells_x_) % num_cells_x_;
        cell_y = (cell_y % num_cells_y_ + num_cells_y_) % num_cells_y_;

        return {cell_x, cell_y};
    }
    

    // Calculate displacement between two points with periodic boundary conditions
    double displacement(const vec& current, const vec& last) const {
        return current.distance(last, box_);
    }

    // Concatenate positions from both species into a single vector
    void concatenate_positions() {
        all_positions_.clear();
        all_positions_.insert(all_positions_.end(), actin_positions_.begin(), actin_positions_.end());
        all_positions_.insert(all_positions_.end(), myosin_positions_.begin(), myosin_positions_.end());
    }

    // Track species types for concatenated particles in all_positions_
    void track_species_types() {
        species_types_.clear();
        species_types_.insert(species_types_.end(), actin_positions_.size(), ParticleType::Actin);
        species_types_.insert(species_types_.end(), myosin_positions_.size(), ParticleType::Myosin);
    }
};

#endif // NEIGHBORLIST_H
