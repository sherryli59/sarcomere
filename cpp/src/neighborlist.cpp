#include "neighborlist.h"

// Constructor with parameters for 3D.
NeighborList::NeighborList(double cutoff_radius, const std::vector<double>& box, double threshold)
    : cutoff_radius_(cutoff_radius), threshold_(threshold), box_(box)
{
    // Assume box has three dimensions.
    num_cells_x_ = static_cast<int>(std::floor(box[0] / cutoff_radius));
    num_cells_y_ = static_cast<int>(std::floor(box[1] / cutoff_radius));
    num_cells_z_ = static_cast<int>(std::floor(box[2] / cutoff_radius));

    cell_size_x_ = box[0] / num_cells_x_;
    cell_size_y_ = box[1] / num_cells_y_;
    cell_size_z_ = box[2] / num_cells_z_;

    // Precompute offsets for neighboring cells in 3D.
    // Using a triple nested loop for dx, dy, dz from -1 to 1.
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                // If we have more than one cell in every dimension, include all neighbors.
                if (num_cells_x_ > 2 && num_cells_y_ > 2 && num_cells_z_ > 2) {
                    neighboring_cells.emplace_back(std::make_tuple(dx, dy, dz));
                    continue;
                }
                // Always include (0,0,0).
                if (dx == 0 && dy == 0 && dz == 0) {
                    neighboring_cells.emplace_back(std::make_tuple(dx, dy, dz));
                    continue;
                }
                // If there is only one cell in a given direction, skip offsets in that direction.
                if (num_cells_x_ == 1 && dx != 0) continue;
                if (num_cells_y_ == 1 && dy != 0) continue;
                if (num_cells_z_ == 1 && dz != 0) continue;
                // For two cells, only include positive offsets.
                if (num_cells_x_ == 2 && dx == -1) continue;
                if (num_cells_y_ == 2 && dy == -1) continue;
                if (num_cells_z_ == 2 && dz == -1) continue;
                // Otherwise include this neighbor.
                neighboring_cells.emplace_back(std::make_tuple(dx, dy, dz));
            }
        }
    }
}

// Default constructor.
NeighborList::NeighborList() {}

// Compute the 3D cell index for a position using PBC.
// Returns a tuple of (cell_x, cell_y, cell_z).
std::tuple<int, int, int> NeighborList::get_cell_index(const vec& position) const {
    int cell_x = static_cast<int>(std::floor((position.x + box_[0]) / cell_size_x_)) % num_cells_x_;
    int cell_y = static_cast<int>(std::floor((position.y + box_[1]) / cell_size_y_)) % num_cells_y_;
    int cell_z = static_cast<int>(std::floor((position.z + box_[2]) / cell_size_z_)) % num_cells_z_;

    // Correct negative indices.
    cell_x = (cell_x % num_cells_x_ + num_cells_x_) % num_cells_x_;
    cell_y = (cell_y % num_cells_y_ + num_cells_y_) % num_cells_y_;
    cell_z = (cell_z % num_cells_z_ + num_cells_z_) % num_cells_z_;

    return std::make_tuple(cell_x, cell_y, cell_z);
}

// Initialize the neighbor list using separate actin and myosin positions.
void NeighborList::initialize(const std::vector<vec>& actin_positions, const std::vector<vec>& myosin_positions) {
    actin_positions_ = actin_positions;
    n_actins_ = actin_positions.size();
    myosin_positions_ = myosin_positions;

    // Save current positions for later displacement checking.
    last_actin_positions_ = actin_positions;
    last_myosin_positions_ = myosin_positions;

    // Concatenate positions into a single vector.
    concatenate_positions();

    // Record species type for each particle.
    track_species_types();

    // Resize the neighbor list container.
    neighbor_list_.resize(all_positions_.size());

    // Build the neighbor list.
    rebuild_neighbor_list();
}

// Rebuild the neighbor list using a cell list approach in 3D.
void NeighborList::rebuild_neighbor_list() {
    printf("Rebuilding neighbor list\n");

    neighbor_list_.clear();
    neighbor_list_.resize(all_positions_.size());
    cell_list_.clear();  // cell_list_ is now a map from tuple<int, int, int> to vector<int>

    // Populate the cell list with particle indices.
    for (size_t i = 0; i < all_positions_.size(); ++i) {
        auto cell = get_cell_index(all_positions_[i]);
        cell_list_[cell].push_back(i);
    }

    // Create thread-local storage for neighbor lists.
    std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> thread_local_neighbor_list(omp_get_max_threads());
    #pragma omp parallel for
    for (auto& local_list : thread_local_neighbor_list) {
        local_list.resize(all_positions_.size());
    }

    // Compute neighbors in parallel.
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < all_positions_.size(); ++i) {
            vec position = all_positions_[i];
            auto cell = get_cell_index(position);
            
            // Check particles in the current cell and its neighbors in 3D.
            for (const auto& offset_tuple : neighboring_cells) {
                int dx = std::get<0>(offset_tuple);
                int dy = std::get<1>(offset_tuple);
                int dz = std::get<2>(offset_tuple);
                int neighbor_x = (std::get<0>(cell) + dx + num_cells_x_) % num_cells_x_;
                int neighbor_y = (std::get<1>(cell) + dy + num_cells_y_) % num_cells_y_;
                int neighbor_z = (std::get<2>(cell) + dz + num_cells_z_) % num_cells_z_;
                auto neighbor_cell = std::make_tuple(neighbor_x, neighbor_y, neighbor_z);
                if (cell_list_.count(neighbor_cell) == 0)
                    continue;
                for (int j : cell_list_[neighbor_cell]) {
                    if (i >= static_cast<size_t>(j))
                        continue;
                    double distance = all_positions_[i].distance(all_positions_[j], box_);
                    if (distance < cutoff_radius_) {
                        thread_local_neighbor_list[thread_id][i].emplace_back(j, species_types_[j]);
                        thread_local_neighbor_list[thread_id][j].emplace_back(i, species_types_[i]);
                    }
                }
            }
        }
    }

    // Merge thread-local neighbor lists into the main neighbor list.
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < all_positions_.size(); ++i) {
        for (const auto& local_list : thread_local_neighbor_list) {
            neighbor_list_[i].insert(
                neighbor_list_[i].end(),
                local_list[i].begin(),
                local_list[i].end()
            );
        }
    }

    // Update the last known positions.
    last_actin_positions_ = actin_positions_;
    last_myosin_positions_ = myosin_positions_;
}

// Set new positions for the two species.
void NeighborList::set_species_positions(const std::vector<vec>& actin_positions, const std::vector<vec>& myosin_positions) {
    actin_positions_ = actin_positions;
    n_actins_ = actin_positions.size();
    myosin_positions_ = myosin_positions;
    concatenate_positions();
    track_species_types();
}

// Check if the neighbor list needs rebuilding based on displacement.
bool NeighborList::needs_rebuild() const {
    bool needs_rebuild = false;
    // Check actin displacements.
    #pragma omp parallel for
    for (size_t i = 0; i < actin_positions_.size(); ++i) {
        if (displacement(actin_positions_[i], last_actin_positions_[i]) > threshold_) {
            needs_rebuild = true;
            #pragma omp cancel for
        }
    }
    #pragma omp barrier 
    if (needs_rebuild)
        return true;
    // Check myosin displacements.
    #pragma omp parallel for
    for (size_t i = 0; i < myosin_positions_.size(); ++i) {
        if (displacement(myosin_positions_[i], last_myosin_positions_[i]) > threshold_) {
            needs_rebuild = true;
            #pragma omp cancel for
        }
    }
    #pragma omp barrier 
    return needs_rebuild;
}

// Return the neighbor list for a specific particle.
const std::vector<std::pair<int, ParticleType>>& NeighborList::get_neighbors(int index) const {
    return neighbor_list_[index];
}

// Return the neighbors of a specific particle, separated by species.
std::pair<std::vector<int>, std::vector<int>> NeighborList::get_neighbors_by_type(int index) const {
    std::vector<int> actin_neighbors;
    std::vector<int> myosin_neighbors;

    for (const auto& entry : neighbor_list_[index]) {
        int neighbor_index = entry.first;
        ParticleType type = entry.second;
        if (type == ParticleType::Actin) {
            actin_neighbors.push_back(neighbor_index);
        } else if (type == ParticleType::Myosin) {
            myosin_neighbors.push_back(neighbor_index - static_cast<int>(n_actins_));
        }
    }
    return std::make_pair(actin_neighbors, myosin_neighbors);
}


// Compute the displacement between two positions (with periodic boundaries).
double NeighborList::displacement(const vec& current, const vec& last) const {
    return current.distance(last, box_);
}

// Concatenate the positions from actin and myosin into one vector.
void NeighborList::concatenate_positions() {
    all_positions_.clear();
    all_positions_.insert(all_positions_.end(), actin_positions_.begin(), actin_positions_.end());
    all_positions_.insert(all_positions_.end(), myosin_positions_.begin(), myosin_positions_.end());
}

// Track the species type for each particle in the concatenated list.
void NeighborList::track_species_types() {
    species_types_.clear();
    species_types_.insert(species_types_.end(), actin_positions_.size(), ParticleType::Actin);
    species_types_.insert(species_types_.end(), myosin_positions_.size(), ParticleType::Myosin);
}
