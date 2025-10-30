#include "neighborlist.h"

// Constructor with parameters for 3D.
NeighborList::NeighborList(double cutoff_radius, const std::vector<double>& box, double threshold,
                           const std::array<bool,3>& periodic_axes)
    : cutoff_radius_(cutoff_radius), threshold_(threshold), box_(box), periodic_axes_(periodic_axes)
{
    // Assume box has three dimensions.
    double Lx = box.size() > 0 ? box[0] : cutoff_radius;
    double Ly = box.size() > 1 ? box[1] : cutoff_radius;
    double Lz = box.size() > 2 ? box[2] : cutoff_radius;
    num_cells_x_ = std::max(1, static_cast<int>(std::floor(Lx / cutoff_radius)));
    num_cells_y_ = std::max(1, static_cast<int>(std::floor(Ly / cutoff_radius)));
    num_cells_z_ = std::max(1, static_cast<int>(std::floor(Lz / cutoff_radius)));
    cell_size_x_ = (num_cells_x_ > 0) ? Lx / num_cells_x_ : Lx;
    cell_size_y_ = (num_cells_y_ > 0) ? Ly / num_cells_y_ : Ly;
    cell_size_z_ = (num_cells_z_ > 0) ? Lz / num_cells_z_ : Lz;

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

void NeighborList::set_periodic_axes(const std::array<bool,3>& periodic_axes) {
    periodic_axes_ = periodic_axes;
}

// Compute the 3D cell index for a position using PBC.
// Returns a tuple of (cell_x, cell_y, cell_z).
std::tuple<int, int, int> NeighborList::get_cell_index(double x, double y, double z) const {
    auto compute_index = [](double coord, double L, double cell_size, int num_cells, bool periodic) {
        if (cell_size <= 0 || num_cells <= 0) {
            return 0;
        }
        double shifted = coord + 0.5 * L;
        int idx = static_cast<int>(std::floor(shifted / cell_size));
        if (periodic) {
            int mod = idx % num_cells;
            if (mod < 0) {
                mod += num_cells;
            }
            return mod;
        }
        return std::clamp(idx, 0, num_cells - 1);
    };
    double Lx = box_.size() > 0 ? box_[0] : cell_size_x_;
    double Ly = box_.size() > 1 ? box_[1] : cell_size_y_;
    double Lz = box_.size() > 2 ? box_[2] : cell_size_z_;
    int cell_x = compute_index(x, Lx, cell_size_x_, num_cells_x_, periodic_axes_[0]);
    int cell_y = compute_index(y, Ly, cell_size_y_, num_cells_y_, periodic_axes_[1]);
    int cell_z = compute_index(z, Lz, cell_size_z_, num_cells_z_, periodic_axes_[2]);

    return std::make_tuple(cell_x, cell_y, cell_z);
}

// Initialize the neighbor list using separate actin and myosin positions.
void NeighborList::initialize(const std::vector<double>& actin_x, const std::vector<double>& actin_y, const std::vector<double>& actin_z,
                              const std::vector<double>& myosin_x, const std::vector<double>& myosin_y, const std::vector<double>& myosin_z) {
    // Store positions in structure-of-arrays form.
    set_species_positions(actin_x, actin_y, actin_z, myosin_x, myosin_y, myosin_z);

    // Save current positions for later displacement checking.
    last_actin_x_ = actin_x_;
    last_actin_y_ = actin_y_;
    last_actin_z_ = actin_z_;
    last_myosin_x_ = myosin_x_;
    last_myosin_y_ = myosin_y_;
    last_myosin_z_ = myosin_z_;

    // Resize the neighbor list container.
    neighbor_list_.resize(all_x_.size());

    // Build the neighbor list.
    rebuild_neighbor_list();
}

// In neighborlist.cpp

// Helper function to clear and resize internal containers.
void NeighborList::clearAndResizeContainers() {
    neighbor_list_.clear();
    neighbor_list_.resize(all_x_.size());
    cell_list_.clear();
    cell_list_.resize(num_cells_x_ * num_cells_y_ * num_cells_z_);
}

// Step 1: Populate the cell list with particle indices.
void NeighborList::populateCellList() {
    for (size_t i = 0; i < all_x_.size(); ++i) {
        auto cell = get_cell_index(all_x_[i], all_y_[i], all_z_[i]);
        int idx = std::get<0>(cell) + num_cells_x_ * (std::get<1>(cell) + num_cells_y_ * std::get<2>(cell));
        cell_list_[idx].push_back(i);
    }
}

// Step 2: Compute neighbor pairs in parallel.
void NeighborList::computeNeighbors(std::vector<std::vector<std::pair<size_t, size_t>>> &tls)
{
    tls.resize(omp_get_max_threads());
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto &local_pairs = tls[thread_id];
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < all_x_.size(); ++i) {
            auto cell = get_cell_index(all_x_[i], all_y_[i], all_z_[i]);

            // Check particles in the current cell and its neighbors in 3D.
            for (const auto& offset_tuple : neighboring_cells) {
                int dx = std::get<0>(offset_tuple);
                int dy = std::get<1>(offset_tuple);
                int dz = std::get<2>(offset_tuple);
                int neighbor_x = std::get<0>(cell) + dx;
                int neighbor_y = std::get<1>(cell) + dy;
                int neighbor_z = std::get<2>(cell) + dz;
                if (periodic_axes_[0]) {
                    neighbor_x = (neighbor_x % num_cells_x_ + num_cells_x_) % num_cells_x_;
                } else if (neighbor_x < 0 || neighbor_x >= num_cells_x_) {
                    continue;
                }
                if (periodic_axes_[1]) {
                    neighbor_y = (neighbor_y % num_cells_y_ + num_cells_y_) % num_cells_y_;
                } else if (neighbor_y < 0 || neighbor_y >= num_cells_y_) {
                    continue;
                }
                if (periodic_axes_[2]) {
                    neighbor_z = (neighbor_z % num_cells_z_ + num_cells_z_) % num_cells_z_;
                } else if (neighbor_z < 0 || neighbor_z >= num_cells_z_) {
                    continue;
                }
                int neighbor_index = neighbor_x + num_cells_x_ * (neighbor_y + num_cells_y_ * neighbor_z);
                if (cell_list_[neighbor_index].empty())
                    continue;
                for (int j : cell_list_[neighbor_index]) {
                    if (i >= static_cast<size_t>(j))
                        continue;
                    double dx = utils::wrap_axis(all_x_[i] - all_x_[j], box_.size() > 0 ? box_[0] : 0.0, periodic_axes_[0]);
                    double dy = utils::wrap_axis(all_y_[i] - all_y_[j], box_.size() > 1 ? box_[1] : 0.0, periodic_axes_[1]);
                    double dz = utils::wrap_axis(all_z_[i] - all_z_[j], box_.size() > 2 ? box_[2] : 0.0, periodic_axes_[2]);
                    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
                    if (distance < cutoff_radius_) {
                        local_pairs.emplace_back(i, j);
                    }
                }
            }
        }
    }
}

// Step 3: Merge thread-local pairs into the main neighbor list.
void NeighborList::mergeThreadLocalPairs(const std::vector<std::vector<std::pair<size_t, size_t>>> &tls)
{
    std::vector<size_t> counts(all_x_.size(), 0);
    for (const auto &local_pairs : tls) {
        for (const auto &p : local_pairs) {
            counts[p.first]++;
            counts[p.second]++;
        }
    }

    #pragma omp parallel for
    for (size_t i = 0; i < all_x_.size(); ++i) {
        neighbor_list_[i].clear();
        neighbor_list_[i].reserve(counts[i]);
    }

    for (const auto &local_pairs : tls) {
        for (const auto &p : local_pairs) {
            size_t i = p.first;
            size_t j = p.second;
            neighbor_list_[i].emplace_back(j, species_types_[j]);
            neighbor_list_[j].emplace_back(i, species_types_[i]);
        }
    }
}

// Step 5: Update last known positions.
void NeighborList::updateLastKnownPositions() {
    last_actin_x_ = actin_x_;
    last_actin_y_ = actin_y_;
    last_actin_z_ = actin_z_;
    last_myosin_x_ = myosin_x_;
    last_myosin_y_ = myosin_y_;
    last_myosin_z_ = myosin_z_;
}

// Main rebuild function that calls each step.
void NeighborList::rebuild_neighbor_list() {
    // Step 0: Clear and prepare containers.
    clearAndResizeContainers();

    // Step 1: Populate cell list.
    populateCellList();

    // Step 2: Compute neighbor pairs in parallel.
    std::vector<std::vector<std::pair<size_t, size_t>>> tls;
    computeNeighbors(tls);

    // Step 3: Merge thread-local pairs.
    mergeThreadLocalPairs(tls);

    // Step 4: Update last known positions.
    updateLastKnownPositions();
}

// // Rebuild the neighbor list using a cell list approach in 3D.
// void NeighborList::rebuild_neighbor_list() {
//     //printf("Rebuilding neighbor list\n");

//     neighbor_list_.clear();
//     neighbor_list_.resize(all_positions_.size());
//     cell_list_.clear();  // cell_list_ is now a map from tuple<int, int, int> to vector<int>

//     // Populate the cell list with particle indices.
//     for (size_t i = 0; i < all_positions_.size(); ++i) {
//         auto cell = get_cell_index(all_positions_[i]);
//         cell_list_[cell].push_back(i);
//     }

//     // Create thread-local storage for neighbor lists.
//     std::vector<std::vector<std::vector<std::pair<size_t, ParticleType>>>> thread_local_neighbor_list(omp_get_max_threads());
//     #pragma omp parallel for
//     for (auto& local_list : thread_local_neighbor_list) {
//         local_list.resize(all_positions_.size());
//     }

//     // Compute neighbors in parallel.
//     #pragma omp parallel
//     {
//         int thread_id = omp_get_thread_num();
//         #pragma omp for schedule(dynamic)
//         for (size_t i = 0; i < all_positions_.size(); ++i) {
//             vec position = all_positions_[i];
//             auto cell = get_cell_index(position);
            
//             // Check particles in the current cell and its neighbors in 3D.
//             for (const auto& offset_tuple : neighboring_cells) {
//                 int dx = std::get<0>(offset_tuple);
//                 int dy = std::get<1>(offset_tuple);
//                 int dz = std::get<2>(offset_tuple);
//                 int neighbor_x = (std::get<0>(cell) + dx + num_cells_x_) % num_cells_x_;
//                 int neighbor_y = (std::get<1>(cell) + dy + num_cells_y_) % num_cells_y_;
//                 int neighbor_z = (std::get<2>(cell) + dz + num_cells_z_) % num_cells_z_;
//                 auto neighbor_cell = std::make_tuple(neighbor_x, neighbor_y, neighbor_z);
//                 if (cell_list_.count(neighbor_cell) == 0)
//                     continue;
//                 for (int j : cell_list_[neighbor_cell]) {
//                     if (i >= static_cast<size_t>(j))
//                         continue;
//                     double distance = all_positions_[i].distance(all_positions_[j], box_);
//                     if (distance < cutoff_radius_) {
//                         thread_local_neighbor_list[thread_id][i].emplace_back(j, species_types_[j]);
//                         thread_local_neighbor_list[thread_id][j].emplace_back(i, species_types_[i]);
//                     }
//                 }
//             }
//         }
//     }

//     // Merge thread-local neighbor lists into the main neighbor list.
//     #pragma omp parallel for schedule(dynamic)
//     for (size_t i = 0; i < all_positions_.size(); ++i) {
//         for (const auto& local_list : thread_local_neighbor_list) {
//             neighbor_list_[i].insert(
//                 neighbor_list_[i].end(),
//                 local_list[i].begin(),
//                 local_list[i].end()
//             );
//         }
//     }

//     // Update the last known positions.
//     last_actin_positions_ = actin_positions_;
//     last_myosin_positions_ = myosin_positions_;
// }

// Set new positions for the two species.
void NeighborList::set_species_positions(const std::vector<double>& actin_x, const std::vector<double>& actin_y, const std::vector<double>& actin_z,
                                         const std::vector<double>& myosin_x, const std::vector<double>& myosin_y, const std::vector<double>& myosin_z) {
    n_actins_ = actin_x.size();
    actin_x_ = actin_x;
    actin_y_ = actin_y;
    actin_z_ = actin_z;

    myosin_x_ = myosin_x;
    myosin_y_ = myosin_y;
    myosin_z_ = myosin_z;

    concatenate_positions();
    track_species_types();
}

// Check if the neighbor list needs rebuilding based on displacement.
bool NeighborList::needs_rebuild() const {
    bool needs_rebuild = false;
    // Check actin displacements.
    #pragma omp parallel for
    for (size_t i = 0; i < actin_x_.size(); ++i) {
        if (displacement(actin_x_[i], actin_y_[i], actin_z_[i],
                         last_actin_x_[i], last_actin_y_[i], last_actin_z_[i]) > threshold_) {
            needs_rebuild = true;
            #pragma omp cancel for
        }
    }
    #pragma omp barrier
    if (needs_rebuild)
        return true;
    // Check myosin displacements.
    #pragma omp parallel for
    for (size_t i = 0; i < myosin_x_.size(); ++i) {
        if (displacement(myosin_x_[i], myosin_y_[i], myosin_z_[i],
                         last_myosin_x_[i], last_myosin_y_[i], last_myosin_z_[i]) > threshold_) {
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
double NeighborList::displacement(double cx, double cy, double cz,
                                  double lx, double ly, double lz) const {
    double dx = utils::wrap_axis(cx - lx, box_.size() > 0 ? box_[0] : 0.0, periodic_axes_[0]);
    double dy = utils::wrap_axis(cy - ly, box_.size() > 1 ? box_[1] : 0.0, periodic_axes_[1]);
    double dz = utils::wrap_axis(cz - lz, box_.size() > 2 ? box_[2] : 0.0, periodic_axes_[2]);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Concatenate the positions from actin and myosin into one vector.
void NeighborList::concatenate_positions() {
    all_x_.clear();
    all_y_.clear();
    all_z_.clear();
    all_x_.insert(all_x_.end(), actin_x_.begin(), actin_x_.end());
    all_x_.insert(all_x_.end(), myosin_x_.begin(), myosin_x_.end());
    all_y_.insert(all_y_.end(), actin_y_.begin(), actin_y_.end());
    all_y_.insert(all_y_.end(), myosin_y_.begin(), myosin_y_.end());
    all_z_.insert(all_z_.end(), actin_z_.begin(), actin_z_.end());
    all_z_.insert(all_z_.end(), myosin_z_.begin(), myosin_z_.end());
}

// Track the species type for each particle in the concatenated list.
void NeighborList::track_species_types() {
    species_types_.clear();
    species_types_.insert(species_types_.end(), actin_x_.size(), ParticleType::Actin);
    species_types_.insert(species_types_.end(), myosin_x_.size(), ParticleType::Myosin);
}
