#include "neighborlist.h"
#include <cstdio>
#include <omp.h>

// Constructor with parameters.
NeighborList::NeighborList(double cutoff_radius, const std::vector<double>& box, double threshold)
    : cutoff_radius_(cutoff_radius), threshold_(threshold), box_(box)
{
    num_cells_x_ = static_cast<int>(std::floor(box[0] / cutoff_radius));
    num_cells_y_ = static_cast<int>(std::floor(box[1] / cutoff_radius));

    // Adjust cell size to ensure full coverage.
    cell_size_x_ = box[0] / num_cells_x_;
    cell_size_y_ = box[1] / num_cells_y_;

    // Precompute offsets for neighboring cells.
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            // For larger grids include all neighbors.
            if (num_cells_x_ > 2 && num_cells_y_ > 2) {
                neighboring_cells.emplace_back(dx, dy);
                continue;
            }
            // Always include the current cell.
            if (dx == 0 && dy == 0) {
                neighboring_cells.emplace_back(dx, dy);
                continue;
            }
            // If only one cell exists in a given direction, skip off‐grid neighbors.
            if (num_cells_x_ == 1 && dx != 0) continue;
            if (num_cells_y_ == 1 && dy != 0) continue;
            // For two cells, only include positive offsets.
            if (num_cells_x_ == 2 && dx == -1) continue;
            if (num_cells_y_ == 2 && dy == -1) continue;
            neighboring_cells.emplace_back(dx, dy);
        }
    }
}

// Default constructor.
NeighborList::NeighborList() {}

// Initialize the neighbor list using separate actin and myosin positions.
void NeighborList::initialize(const std::vector<vec>& actin_positions,
                              const std::vector<vec>& myosin_positions) {
    actin_positions_ = actin_positions;
    n_actins_ = actin_positions.size();
    myosin_positions_ = myosin_positions;

    // Save current positions for later displacement checking.
    last_actin_positions_ = actin_positions;
    last_myosin_positions_ = myosin_positions;

    // Concatenate positions and record species.
    concatenate_positions();
    track_species_types();

    // Resize the neighbor list container.
    neighbor_list_.resize(all_positions_.size());

    // Build the neighbor list.
    rebuild_neighbor_list();
}

// Rebuild the neighbor list using a cell list approach.
void NeighborList::rebuild_neighbor_list() {
    printf("Rebuilding neighbor list\n");

    // Clear any previous neighbor list.
    neighbor_list_.clear();
    neighbor_list_.resize(all_positions_.size());

    // --- Build the cell list as a flat vector ---
    // Resize the cell list to cover every cell.
    cell_list_.clear();
    cell_list_.resize(num_cells_x_ * num_cells_y_);

    // A lambda to compute the flat index from a (cell_x, cell_y) pair.
    auto get_flat_index = [this](const std::pair<int, int>& cell) -> int {
        return cell.first + cell.second * num_cells_x_;
    };

    // Populate the cell list with particle indices.
    for (size_t i = 0; i < all_positions_.size(); ++i) {
        std::pair<int, int> cell = get_cell_index(all_positions_[i]);
        int idx = get_flat_index(cell);
        cell_list_[idx].push_back(i);
    }

    // --- Create thread-local storage for neighbor lists ---
    int max_threads = omp_get_max_threads();
    std::vector<std::vector<std::vector<std::pair<int, ParticleType>>>> thread_local_neighbor_list(max_threads);
    #pragma omp parallel for schedule(static)
    for (int t = 0; t < max_threads; ++t) {
        thread_local_neighbor_list[t].resize(all_positions_.size());
    }

    // Precompute the cutoff squared to avoid computing square roots.
    double cutoff2 = cutoff_radius_ * cutoff_radius_;

    // --- Compute neighbors in parallel ---
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static)
        for (size_t i = 0; i < all_positions_.size(); ++i) {
            const vec& pos_i = all_positions_[i];
            std::pair<int, int> cell = get_cell_index(pos_i);

            // Loop over the current cell and all neighboring cells.
            for (const auto& offset : neighboring_cells) {
                int neighbor_x = (cell.first + offset.first + num_cells_x_) % num_cells_x_;
                int neighbor_y = (cell.second + offset.second + num_cells_y_) % num_cells_y_;
                std::pair<int, int> neighbor_cell = {neighbor_x, neighbor_y};
                int cell_idx = get_flat_index(neighbor_cell);

                // Skip if there are no particles in this cell.
                if (cell_idx < 0 || cell_idx >= static_cast<int>(cell_list_.size()) ||
                    cell_list_[cell_idx].empty())
                    continue;

                for (int j : cell_list_[cell_idx]) {
                    if (i >= static_cast<size_t>(j))
                        continue;
                    // Use squared distance to avoid the expensive sqrt.
                    double dist2 = pos_i.distance_squared(all_positions_[j], box_);
                    if (dist2 < cutoff2) {
                        thread_local_neighbor_list[thread_id][i].emplace_back(j, species_types_[j]);
                        thread_local_neighbor_list[thread_id][j].emplace_back(i, species_types_[i]);
                    }
                }
            }
        }
    }

    // --- Merge the thread-local neighbor lists into the main neighbor list ---
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < all_positions_.size(); ++i) {
        for (int t = 0; t < max_threads; ++t) {
            neighbor_list_[i].insert(
                neighbor_list_[i].end(),
                thread_local_neighbor_list[t][i].begin(),
                thread_local_neighbor_list[t][i].end()
            );
        }
    }

    // Update the last known positions.
    last_actin_positions_ = actin_positions_;
    last_myosin_positions_ = myosin_positions_;
}

// Set new positions for the two species.
void NeighborList::set_species_positions(const std::vector<vec>& actin_positions,
                                         const std::vector<vec>& myosin_positions) {
    actin_positions_ = actin_positions;
    // (If you store n_actins_ as a member, update it accordingly.)
    n_actins_ = actin_positions.size();
    myosin_positions_ = myosin_positions;
    concatenate_positions();
    track_species_types();
}

// Check if the neighbor list needs rebuilding based on displacement.
bool NeighborList::needs_rebuild() const {
    bool rebuild = false;
    // Check actin displacements.
    #pragma omp parallel for schedule(static) shared(rebuild)
    for (size_t i = 0; i < actin_positions_.size(); ++i) {
        if (displacement(actin_positions_[i], last_actin_positions_[i]) > threshold_) {
            #pragma omp atomic write
            rebuild = true;
        }
    }
    if (rebuild)
        return true;
    // Check myosin displacements.
    #pragma omp parallel for schedule(static) shared(rebuild)
    for (size_t i = 0; i < myosin_positions_.size(); ++i) {
        if (displacement(myosin_positions_[i], last_myosin_positions_[i]) > threshold_) {
            #pragma omp atomic write
            rebuild = true;
        }
    }
    return rebuild;
}

// Return the neighbor list for a specific particle.
const std::vector<std::pair<int, ParticleType>>& NeighborList::get_neighbors(int index) const {
    return neighbor_list_[index];
}

// Return the neighbors of a particle separated by species.
std::pair<std::vector<int>, std::vector<int>> NeighborList::get_neighbors_by_type(int index) const {
    std::vector<int> actin_neighbors;
    std::vector<int> myosin_neighbors;

    for (const auto& entry : neighbor_list_[index]) {
        int neighbor_index = entry.first;
        ParticleType type = entry.second;
        if (type == ParticleType::Actin) {
            actin_neighbors.push_back(neighbor_index);
        } else if (type == ParticleType::Myosin) {
            // Adjust index for myosin particles (assuming actin particles come first).
            actin_neighbors; // (unused here)
            myosin_neighbors.push_back(neighbor_index - static_cast<int>(n_actins_));
        }
    }
    return {actin_neighbors, myosin_neighbors};
}

// Get the cell index corresponding to a position.
std::pair<int, int> NeighborList::get_cell_index(const vec& position) const {
    int cell_x = static_cast<int>(std::floor((position.x + box_[0]) / cell_size_x_)) % num_cells_x_;
    int cell_y = static_cast<int>(std::floor((position.y + box_[1]) / cell_size_y_)) % num_cells_y_;

    // Correct for negative values.
    cell_x = (cell_x % num_cells_x_ + num_cells_x_) % num_cells_x_;
    cell_y = (cell_y % num_cells_y_ + num_cells_y_) % num_cells_y_;

    return {cell_x, cell_y};
}

// Compute the displacement between two positions (using the distance method).
double NeighborList::displacement(const vec& current, const vec& last) const {
    return current.distance(last, box_);
}

// Concatenate actin and myosin positions into one vector.
void NeighborList::concatenate_positions() {
    all_positions_.clear();
    all_positions_.insert(all_positions_.end(), actin_positions_.begin(), actin_positions_.end());
    all_positions_.insert(all_positions_.end(), myosin_positions_.begin(), myosin_positions_.end());
}

// Record the species type for each particle in the concatenated list.
void NeighborList::track_species_types() {
    species_types_.clear();
    species_types_.insert(species_types_.end(), actin_positions_.size(), ParticleType::Actin);
    species_types_.insert(species_types_.end(), myosin_positions_.size(), ParticleType::Myosin);
}
