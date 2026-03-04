#include "h5_utils.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include "omp.h"
#include "geometry.h"

//---------------------------------------------------------------------
// Function Definitions
//---------------------------------------------------------------------

std::vector<std::vector<int>> serializeActinIndicesPerActin(
    utils::MoleculeConnection& actinIndicesPerActin, int n_actins, int& max_bonds)
{
    std::vector<std::vector<int>> serialized(n_actins);
    #pragma omp parallel for
    for (int i = 0; i < n_actins; ++i) {
        serialized[i] = actinIndicesPerActin.getConnections(i);
    }
    // Pad each list to have the same length (max_bonds).
    for (auto& bonds : serialized) {
        while (bonds.size() < static_cast<size_t>(max_bonds)) {
            bonds.push_back(-1);  // Use -1 as a placeholder for no bond.
        }
    }
    return serialized;
}


std::vector<double> flatten_3d_array(const Filament::VecArray& array)
{
    size_t rows = array.size();
    std::vector<double> flattened(rows * 3);
    #pragma omp parallel for
    for (size_t i = 0; i < rows; ++i) {
        flattened[i * 3]     = array[i].x;
        flattened[i * 3 + 1] = array[i].y;
        flattened[i * 3 + 2] = array[i].z;
    }
    return flattened;
}

std::vector<double> flatten_2d_array(std::vector<std::vector<double>> array)
{
    size_t rows = array.size();
    size_t cols = array[0].size();
    std::vector<double> flattened(rows * cols);
    #pragma omp parallel for
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            flattened[i * cols + j] = array[i][j];
        }
    }
    return flattened;
}

void create_empty_dataset(H5::H5File& file, const std::string& groupName,
                          const std::string& datasetName,
                          const std::vector<hsize_t>& initialDims,
                          const std::vector<hsize_t>& maxDims,
                          const std::vector<hsize_t>& chunkDims)
{
    try {
        // Open the group if it exists; otherwise create it.
        H5::Group group;
        if (!file.nameExists(groupName)) {
            group = file.createGroup(groupName);
        } else {
            group = file.openGroup(groupName);
        }
        // Create a dataspace with unlimited maximum size.
        H5::DataSpace dataspace(initialDims.size(), initialDims.data(), maxDims.data());
        // Create a dataset creation property list and enable chunking.
        H5::DSetCreatPropList prop;
        prop.setChunk(chunkDims.size(), chunkDims.data());
        // Create the dataset.
        H5::FloatType datatype(H5::PredType::IEEE_F64LE);
        H5::DataSet dataset = group.createDataSet(datasetName, datatype, dataspace, prop);
    }
    catch (H5::FileIException& error) {
        error.printErrorStack();
        std::cerr << "Error creating empty dataset: " << datasetName << std::endl;
    }
    catch (H5::DataSetIException& error) {
        error.printErrorStack();
        std::cerr << "Error creating empty dataset: " << datasetName << std::endl;
    }
}

void create_empty_dataset_int(H5::H5File& file, const std::string& groupName,
                          const std::string& datasetName,
                          const std::vector<hsize_t>& initialDims,
                          const std::vector<hsize_t>& maxDims,
                          const std::vector<hsize_t>& chunkDims)
{
    try {
        H5::Group group;
        if (!file.nameExists(groupName)) {
            group = file.createGroup(groupName);
        } else {
            group = file.openGroup(groupName);
        }
        H5::DataSpace dataspace(initialDims.size(), initialDims.data(), maxDims.data());
        H5::DSetCreatPropList prop;
        prop.setChunk(chunkDims.size(), chunkDims.data());
        H5::IntType datatype(H5::PredType::STD_I32LE);
        H5::DataSet dataset = group.createDataSet(datasetName, datatype, dataspace, prop);
    }
    catch (H5::Exception& error) {
        error.printErrorStack();
        std::cerr << "Error creating empty int dataset: " << datasetName << std::endl;
    }
}

void append_to_dataset(H5::Group& group, const std::string& datasetName,
                       const std::vector<double> newData,
                       const std::vector<hsize_t>& newDims)
{
    try {
        // Open the dataset.
        H5::DataSet dataset = group.openDataSet(datasetName);

        // Get the current dataspace and size of the dataset.
        H5::DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();
        std::vector<hsize_t> currentSize(rank);
        filespace.getSimpleExtentDims(currentSize.data(), NULL);

        // Diagnostic: if shapes don't match, print both and throw.
        if (newDims.size() != currentSize.size()) {
            std::cerr << "Dimension mismatch when appending to dataset '" << datasetName << "'.\n";
            std::cerr << "  dataset rank=" << rank << " currentSize=[";
            for (size_t ii = 0; ii < currentSize.size(); ++ii) {
                if (ii) std::cerr << ", ";
                std::cerr << currentSize[ii];
            }
            std::cerr << "] newDims=[";
            for (size_t ii = 0; ii < newDims.size(); ++ii) {
                if (ii) std::cerr << ", ";
                std::cerr << newDims[ii];
            }
            std::cerr << "]\n";
            throw std::runtime_error("Dimensions of new data do not match dataset.");
        }
        for (size_t i = 1; i < newDims.size(); ++i) {
            if (newDims[i] != currentSize[i]) {
                std::cerr << "Non-append dimension mismatch for dataset '" << datasetName << "': index " << i
                          << " new=" << newDims[i] << " existing=" << currentSize[i] << "\n";
                throw std::runtime_error("Non-append dimensions do not match dataset.");
            }
        }

        // Extend the dataset along the first dimension.
        std::vector<hsize_t> newSize = currentSize;
        newSize[0] += newDims[0];
        dataset.extend(newSize.data());

        // Select the hyperslab where the new data will be written.
        filespace = dataset.getSpace();
        std::vector<hsize_t> offset(newSize.size(), 0);
        offset[0] = currentSize[0];
        filespace.selectHyperslab(H5S_SELECT_SET, newDims.data(), offset.data());

        // Define the memory space for the new data.
        H5::DataSpace memspace(newDims.size(), newDims.data());

        // Write the new data to the hyperslab.
        dataset.write(newData.data(), H5::PredType::IEEE_F64LE, memspace, filespace);
    }
    catch (H5::FileIException& error) {
        error.printErrorStack();
        std::cerr << "Error appending data to dataset in group: " << datasetName << std::endl;
    }
    catch (H5::DataSetIException& error) {
        error.printErrorStack();
        std::cerr << "Error appending data to dataset in group: " << datasetName << std::endl;
    }
    catch (H5::DataSpaceIException& error) {
        error.printErrorStack();
        std::cerr << "Error appending data to dataset in group: " << datasetName << std::endl;
    }
    catch (std::runtime_error& error) {
        std::cerr << "Runtime error: " << error.what() << std::endl;
    }
}

void append_to_dataset_int(H5::Group& group, const std::string& datasetName,
                       const std::vector<int> newData,
                       const std::vector<hsize_t>& newDims)
{
    try {
        H5::DataSet dataset = group.openDataSet(datasetName);
        H5::DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();
        std::vector<hsize_t> currentSize(rank);
        filespace.getSimpleExtentDims(currentSize.data(), NULL);

        if (newDims.size() != currentSize.size()) {
            std::cerr << "Dimension mismatch when appending to int dataset '" << datasetName << "'.\n";
            std::cerr << "  dataset rank=" << rank << " currentSize=[";
            for (size_t ii = 0; ii < currentSize.size(); ++ii) {
                if (ii) std::cerr << ", ";
                std::cerr << currentSize[ii];
            }
            std::cerr << "] newDims=[";
            for (size_t ii = 0; ii < newDims.size(); ++ii) {
                if (ii) std::cerr << ", ";
                std::cerr << newDims[ii];
            }
            std::cerr << "]\n";
            throw std::runtime_error("Dimensions of new data do not match dataset.");
        }

        std::vector<hsize_t> newSize = currentSize;
        newSize[0] += newDims[0];
        dataset.extend(newSize.data());
        H5::DataSpace newFilespace = dataset.getSpace();
        std::vector<hsize_t> offset(currentSize.size(), 0);
        offset[0] = currentSize[0];
        newFilespace.selectHyperslab(H5S_SELECT_SET, newDims.data(), offset.data());
        H5::DataSpace memspace(newDims.size(), newDims.data());
        dataset.write(newData.data(), H5::PredType::STD_I32LE, memspace, newFilespace);
    }
    catch (H5::Exception& error) {
        error.printErrorStack();
        std::cerr << "Error appending int dataset: " << datasetName << std::endl;
    }
}

void create_file(std::string& filename, Filament& actin, Myosin& myosin,
                 int max_myosin_bonds)
{
    // Remove any existing file with the same name.
    std::remove(filename.c_str());
    H5::H5File file(filename, H5F_ACC_TRUNC);

    hsize_t n_actins  = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);

    // Create datasets for  actin.
    std::vector<hsize_t> initialDims = {0, n_actins, 3};
    std::vector<hsize_t> maxDims     = {H5S_UNLIMITED, n_actins, 3};
    std::vector<hsize_t> chunkDims   = {10, n_actins, 3};
    create_empty_dataset(file, "/actin", "center", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "velocity", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "force", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "torque", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "direction", initialDims, maxDims, chunkDims);

    initialDims = {0, n_actins, 1};
    maxDims     = {H5S_UNLIMITED, n_actins, 1};
    chunkDims   = {10, n_actins, 1};
    create_empty_dataset_int(file, "/actin", "cb_status", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "f_load", initialDims, maxDims, chunkDims);


    // Create datasets for any custom actin features.
    for (const auto& feature : actin.custom_features) {
        create_empty_dataset(file, "/actin", feature.first, initialDims, maxDims, chunkDims);
    }


    // Create datasets for myosin.
    initialDims = {0, n_myosins, 3};
    maxDims     = {H5S_UNLIMITED, n_myosins, 3};
    chunkDims   = {10, n_myosins, 3};
    create_empty_dataset(file, "/myosin", "center", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "velocity", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "force", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "torque", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "direction", initialDims, maxDims, chunkDims);

    initialDims = {0, n_myosins, 1};
    maxDims     = {H5S_UNLIMITED, n_myosins, 1};
    chunkDims   = {10, n_myosins, 1};
    for (const auto& feature : myosin.custom_features) {
        create_empty_dataset(file, "/myosin", feature.first, initialDims, maxDims, chunkDims);
    }
    

    // Create dataset for actin indices (connections) with a fixed maximum number of bonds.
    initialDims = {0, n_actins, 10};
    maxDims     = {H5S_UNLIMITED, n_actins, 10};
    chunkDims   = {10, n_actins, 10};
    create_empty_dataset(file, "/actin", "indices_per_actin", initialDims, maxDims, chunkDims);


    hsize_t max_n_actin_bonds = n_actins * 5;

    initialDims = {0, max_n_actin_bonds, 2};
    maxDims     = {H5S_UNLIMITED, max_n_actin_bonds, 2};
    chunkDims   = {10, max_n_actin_bonds, 2};
    create_empty_dataset(file, "/actin", "bonds", initialDims, maxDims, chunkDims);

    initialDims = {0, max_n_actin_bonds, 1};
    maxDims     = {H5S_UNLIMITED, max_n_actin_bonds, 1};
    chunkDims   = {10, max_n_actin_bonds, 1};
    create_empty_dataset(file, "/actin", "bond_pair_load", initialDims, maxDims, chunkDims);

    // Each catch-bonded actin pair (i, j) gets an entry containing their current distance.
    create_empty_dataset(file, "/actin", "cb_distance", initialDims, maxDims, chunkDims);

    hsize_t max_n_myosin_bonds = n_myosins * 4;
    initialDims = {0, max_n_myosin_bonds, 2};
    maxDims     = {H5S_UNLIMITED, max_n_myosin_bonds, 2};
    chunkDims   = {10, max_n_myosin_bonds, 2};
    create_empty_dataset(file, "/myosin", "bonds", initialDims, maxDims, chunkDims);

    hsize_t max_n_am_bonds = n_myosins * max_myosin_bonds;
    initialDims = {0, max_n_am_bonds, 2};
    maxDims     = {H5S_UNLIMITED, max_n_am_bonds, 2};
    chunkDims   = {10, max_n_am_bonds, 2};
    create_empty_dataset(file, "/actin_myo", "bonds", initialDims, maxDims, chunkDims);

    initialDims = {0, max_n_am_bonds, 1};
    maxDims     = {H5S_UNLIMITED, max_n_am_bonds, 1};
    chunkDims   = {10, max_n_am_bonds, 1};
    create_empty_dataset(file, "/actin_myo", "distance", initialDims, maxDims, chunkDims);

    // Dataset to record catch-bond breakage events:
    // columns: i, j, step, distance, cos_angle,
    //          tension_i, tension_j, crosslink_ratio_i, crosslink_ratio_j,
    //          myosin_count_i, myosin_count_j,
    //          myosin_ids_i[0:max_myosin_bonds-1], myosin_ids_j[0:max_myosin_bonds-1]
    hsize_t event_width = 11 + 2 * max_myosin_bonds;
    initialDims = {0, event_width};
    maxDims     = {H5S_UNLIMITED, event_width};
    chunkDims   = {10, event_width};
    create_empty_dataset(file, "/catch_bond", "breakage", initialDims, maxDims, chunkDims);

    // Dataset to record removals triggered by max_strong_actin_bonds
    // columns: i, j, step, bond_count_i, bond_count_j
    hsize_t limit_width = 5;
    initialDims = {0, limit_width};
    maxDims     = {H5S_UNLIMITED, limit_width};
    chunkDims   = {10, limit_width};
    create_empty_dataset(file, "/catch_bond", "limit_removal", initialDims, maxDims, chunkDims);

    // Dataset to record completed lifetimes for detached actin–actin bonds (seconds)
    // 1D array where each entry is a single completed lifetime value
    initialDims = {0};
    maxDims     = {H5S_UNLIMITED};
    chunkDims   = {100};
    create_empty_dataset(file, "/catch_bond", "completed_lifetimes", initialDims, maxDims, chunkDims);

    // Initialize /state group and datasets so newly created files contain
    // the state arrays expected by code paths that open /state early.

    // current_step: 2D dataset with shape (n_frames, 1) to match save_state append {1,1}
    initialDims = {0, 1};
    maxDims     = {H5S_UNLIMITED, 1};
    chunkDims   = {10, 1};
    create_empty_dataset_int(file, "/state", "current_step", initialDims, maxDims, chunkDims);

    // Flattened actin-actin arrays: store as 2D datasets with shape (n_frames, n_actins*n_actins)
    initialDims = {0, static_cast<hsize_t>(n_actins * n_actins)};
    maxDims     = {H5S_UNLIMITED, static_cast<hsize_t>(n_actins * n_actins)};
    chunkDims   = {1, static_cast<hsize_t>(n_actins * n_actins)};
    create_empty_dataset_int(file, "/state", "actin_actin_bonds_prev", initialDims, maxDims, chunkDims);
    create_empty_dataset_int(file, "/state", "actin_actin_status_prev", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/state", "actin_actin_lifetime_prev", initialDims, maxDims, chunkDims);

    // am_bonds_prev: flattened per-frame array of size actin.n * myosin.n
    initialDims = {0, static_cast<hsize_t>(n_actins * n_myosins)};
    maxDims     = {H5S_UNLIMITED, static_cast<hsize_t>(n_actins * n_myosins)};
    chunkDims   = {1, static_cast<hsize_t>(n_actins * n_myosins)};
    create_empty_dataset_int(file, "/state", "am_bonds_prev", initialDims, maxDims, chunkDims);

    // actin_recovery_until: flattened per-frame array size n_actins * n_actins
    initialDims = {0, static_cast<hsize_t>(n_actins * n_actins)};
    maxDims     = {H5S_UNLIMITED, static_cast<hsize_t>(n_actins * n_actins)};
    chunkDims   = {1, static_cast<hsize_t>(n_actins * n_actins)};
    create_empty_dataset_int(file, "/state", "actin_recovery_until", initialDims, maxDims, chunkDims);
}


void append_to_file(std::string& filename, Filament& actin, Myosin& myosin,
    std::vector<double>& flatActinBonds,
    std::vector<double>& flatMyosinBonds,
    std::vector<double>& flatActinMyosinBonds,
    int max_myosin_bonds)

{
    H5::H5File file(filename, H5F_ACC_RDWR);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));
    H5::Group group_am(file.openGroup("/actin_myo"));

    hsize_t n_actins  = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);

    std::vector<double> flattened_actin_center = flatten_3d_array(actin.center);
    append_to_dataset(group_actin, "center", flattened_actin_center, {1, n_actins, 3});

    std::vector<double> flattened_actin_velocity = flatten_3d_array(actin.velocity);
    append_to_dataset(group_actin, "velocity", flattened_actin_velocity, {1, n_actins, 3});

    std::vector<double> flattened_actin_force = flatten_3d_array(actin.force);
    append_to_dataset(group_actin, "force", flattened_actin_force, {1, n_actins, 3});

    std::vector<double> flattened_actin_torque = flatten_3d_array(actin.torque);
    append_to_dataset(group_actin, "torque", flattened_actin_torque, {1, n_actins, 3});

    
    std::vector<double> flattened_actin_direction = flatten_3d_array(actin.direction);
    append_to_dataset(group_actin, "direction", flattened_actin_direction, {1, n_actins, 3});

    append_to_dataset_int(group_actin, "cb_status", actin.cb_status, {1, n_actins, 1});
    append_to_dataset(group_actin, "f_load", actin.f_load, {1, n_actins, 1});


    for (const auto& feature : actin.custom_features) {
        append_to_dataset(group_actin, feature.first, feature.second, {1, n_actins, 1});
    }

    std::vector<double> flattened_myosin_center = flatten_3d_array(myosin.center);
    append_to_dataset(group_myosin, "center", flattened_myosin_center, {1, n_myosins, 3});

    std::vector<double> flattened_myosin_velocity = flatten_3d_array(myosin.velocity);
    append_to_dataset(group_myosin, "velocity", flattened_myosin_velocity, {1, n_myosins, 3});

    std::vector<double> flattened_myosin_force = flatten_3d_array(myosin.force);
    append_to_dataset(group_myosin, "force", flattened_myosin_force, {1, n_myosins, 3});

    std::vector<double> flattened_myosin_torque = flatten_3d_array(myosin.torque);
    append_to_dataset(group_myosin, "torque",flattened_myosin_torque, {1, n_myosins, 3});

    std::vector<double> flattened_myosin_direction = flatten_3d_array(myosin.direction);
    append_to_dataset(group_myosin, "direction", flattened_myosin_direction, {1, n_myosins, 3});

    for (const auto& feature : myosin.custom_features) {
        append_to_dataset(group_myosin, feature.first, feature.second, {1, n_myosins, 1});
    }

    // Compute derived metrics before padding.
    const std::vector<double>* actin_f_load_cb = nullptr;
    auto cb_feature_it = actin.custom_features.find("f_load_cb");
    if (cb_feature_it != actin.custom_features.end()) {
        actin_f_load_cb = &cb_feature_it->second;
    }

    std::vector<double> actin_pair_loads;
    std::vector<double> actin_cb_distances;
    actin_pair_loads.reserve(flatActinBonds.size() / 2);
    actin_cb_distances.reserve(flatActinBonds.size() / 2);
    for (size_t idx = 0; idx + 1 < flatActinBonds.size(); idx += 2) {
        int a = static_cast<int>(flatActinBonds[idx]);
        int b = static_cast<int>(flatActinBonds[idx + 1]);
        if (a < 0 || b < 0) {
            actin_pair_loads.push_back(-1.0);
            actin_cb_distances.push_back(-1.0);
            continue;
        }
        vec dir_a = actin.direction[a];
        vec dir_b = actin.direction[b];
        double norm_a = dir_a.norm();
        double norm_b = dir_b.norm();
        if (norm_a > 1e-12) {
            dir_a = dir_a / norm_a;
        }
        if (norm_b > 1e-12) {
            dir_b = dir_b / norm_b;
        }
        double cos_val = dir_a.dot(dir_b);
        double load_a = actin_f_load_cb ? (*actin_f_load_cb)[a] : actin.f_load[a];
        double load_b = actin_f_load_cb ? (*actin_f_load_cb)[b] : actin.f_load[b];
        double pair_load = (cos_val < 0.0)
            ? std::abs(cos_val) * std::min(load_a, load_b)
            : 0.0;
        actin_pair_loads.push_back(pair_load);
        double aa_dist = geometry::segment_segment_distance(
            actin.left_end[a], actin.right_end[a],
            actin.left_end[b], actin.right_end[b],
            actin.box, actin.periodic_axes);
        actin_cb_distances.push_back(aa_dist);
    }

    std::vector<double> actin_myo_distances;
    actin_myo_distances.reserve(flatActinMyosinBonds.size() / 2);
    for (size_t idx = 0; idx + 1 < flatActinMyosinBonds.size(); idx += 2) {
        int a = static_cast<int>(flatActinMyosinBonds[idx]);
        int m = static_cast<int>(flatActinMyosinBonds[idx + 1]);
        if (a < 0 || m < 0) {
            actin_myo_distances.push_back(-1.0);
            continue;
        }
        double dist = geometry::segment_segment_distance(
            actin.left_end[a], actin.right_end[a],
            myosin.left_end[m], myosin.right_end[m],
            actin.box, actin.periodic_axes);
        actin_myo_distances.push_back(dist);
    }

    // Compute maximum possible bonds per frame.
    int max_n_actin_bonds = actin.n * 5;
    int max_n_myosin_bonds = myosin.n * 4;
    int max_n_am_bonds = myosin.n * max_myosin_bonds;

    // Compute current number of bonds (each bond occupies two entries).
    int current_n_actinBonds = static_cast<int>(flatActinBonds.size()) / 2;
    int current_n_myosinBonds = static_cast<int>(flatMyosinBonds.size()) / 2;
    int current_n_amBonds = static_cast<int>(flatActinMyosinBonds.size()) / 2;

    // Pad flatActinBonds if needed.
    if (current_n_actinBonds < max_n_actin_bonds) {
        int padRows = max_n_actin_bonds - current_n_actinBonds;
        for (int i = 0; i < padRows; i++) {
            flatActinBonds.push_back(-1);
            flatActinBonds.push_back(-1);
            actin_pair_loads.push_back(-1.0);
            actin_cb_distances.push_back(-1.0);
        }
    }
    // Pad flatMyosinBonds if needed.
    if (current_n_myosinBonds < max_n_myosin_bonds) {
        int padRows = max_n_myosin_bonds - current_n_myosinBonds;
        for (int i = 0; i < padRows; i++) {
            flatMyosinBonds.push_back(-1);
            flatMyosinBonds.push_back(-1);
        }
    }
    // Pad actin-myosin bonds if needed.
    if (current_n_amBonds < max_n_am_bonds) {
        int padRows = max_n_am_bonds - current_n_amBonds;
        for (int i = 0; i < padRows; i++) {
            flatActinMyosinBonds.push_back(-1);
            flatActinMyosinBonds.push_back(-1);
            actin_myo_distances.push_back(-1.0);
        }
    }

    // Now that each flattened vector has a size corresponding to (max_n * 2),
    // the new dimensions for appending one frame are {1, max_n, 2}.
    hsize_t newActinDims[3] = {1, static_cast<hsize_t>(max_n_actin_bonds), 2};
    hsize_t newMyosinDims[3] = {1, static_cast<hsize_t>(max_n_myosin_bonds), 2};
    hsize_t newAmDims[3] = {1, static_cast<hsize_t>(max_n_am_bonds), 2};

    append_to_dataset(group_actin, "bonds", flatActinBonds, { newActinDims[0], newActinDims[1], newActinDims[2] });
    append_to_dataset(group_actin, "bond_pair_load", actin_pair_loads, {1, static_cast<hsize_t>(max_n_actin_bonds), 1});
    append_to_dataset(group_actin, "cb_distance", actin_cb_distances, {1, static_cast<hsize_t>(max_n_actin_bonds), 1});
    append_to_dataset(group_myosin, "bonds", flatMyosinBonds, { newMyosinDims[0], newMyosinDims[1], newMyosinDims[2] });
    append_to_dataset(group_am, "bonds", flatActinMyosinBonds, { newAmDims[0], newAmDims[1], newAmDims[2] });
    append_to_dataset(group_am, "distance", actin_myo_distances, {1, static_cast<hsize_t>(max_n_am_bonds), 1});
    // int max_bonds = 10;
    // // Serialize actinIndicesPerActin.
    // auto serialized_indices = serializeActinIndicesPerActin(actinIndicesPerActin, actin.n, max_bonds);

    // // Flatten the serialized indices into a 1D array.
    // std::vector<double> flattened_indices;
    // for (const auto& indices : serialized_indices) {
    //     for (int index : indices) {
    //         flattened_indices.push_back(static_cast<double>(index));
    //     }
    // }
    // // Append the serialized indices to the dataset.
    // append_to_dataset(group_actin, "indices_per_actin", flattened_indices,
    //                     {1, n_actins, static_cast<hsize_t>(max_bonds)});
}

std::vector<double> load_from_dataset(H5::Group& group, const std::string& datasetName,
                                      std::vector<hsize_t>& dims)
{
    std::vector<double> data;
    try {
        H5::DataSet dataset = group.openDataSet(datasetName);
        H5::DataSpace dataspace = dataset.getSpace();
        int rank = dataspace.getSimpleExtentNdims();
        dims.resize(rank);
        dataspace.getSimpleExtentDims(dims.data(), NULL);

        hsize_t numElements = 1;
        for (hsize_t dim : dims) {
            numElements *= dim;
        }
        data.resize(numElements);
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    }
    catch (H5::Exception& error) {
        std::cerr << "Error reading dataset " << datasetName << " from file" << std::endl;
        throw error;
    }
    return data;
}

int load_from_file(std::string& filename, Filament& actin, Myosin& myosin,
                    std::vector<std::vector<int>>& actin_actin_bonds, int& n_frames, int frame_index)
{
    H5::H5File file(filename, H5F_ACC_RDONLY);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));

    hsize_t n_actins  = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);

    std::vector<hsize_t> dims;
    std::vector<double> actin_center_all = load_from_dataset(group_actin, "center", dims);
    n_frames = static_cast<int>(dims[0]);
    if (n_frames <= 0) {
        throw std::runtime_error("No frames available in actin center dataset for resume.");
    }
    int target_frame = frame_index;
    if (target_frame < 0 || target_frame >= n_frames) {
        target_frame = n_frames - 1;
    }

    size_t actin_frame_stride = static_cast<size_t>(n_actins) * 3;
    size_t actin_start = static_cast<size_t>(target_frame) * actin_frame_stride;
    std::vector<double> actin_center(
        actin_center_all.begin() + actin_start,
        actin_center_all.begin() + actin_start + actin_frame_stride
    );
    for (int i = 0; i < actin.n; i++) {
        actin.center[i].x = actin_center[3 * i];
        actin.center[i].y = actin_center[3 * i + 1];
        actin.center[i].z = actin_center[3 * i + 2];
    }
    std::vector<double> actin_direction_all = load_from_dataset(group_actin, "direction", dims);
    std::vector<double> actin_direction(
        actin_direction_all.begin() + actin_start,
        actin_direction_all.begin() + actin_start + actin_frame_stride
    );
    for (int i = 0; i < actin.n; i++) {
        actin.direction[i].x = actin_direction[3 * i];
        actin.direction[i].y = actin_direction[3 * i + 1];
        actin.direction[i].z = actin_direction[3 * i + 2];
    }

    // Optional per-actin load arrays
    try {
        std::vector<double> actin_f_load_all = load_from_dataset(group_actin, "f_load", dims);
        size_t load_stride = static_cast<size_t>(n_actins);
        size_t load_start = static_cast<size_t>(target_frame) * load_stride;
        for (int i = 0; i < actin.n; ++i) {
            actin.f_load[i] = actin_f_load_all[load_start + i];
        }
    } catch (const H5::Exception&) {
        for (int i = 0; i < actin.n; ++i) {
            actin.f_load[i] = 0.0;
        }
    }
    auto feature_cb_it = actin.custom_features.find("f_load_cb");
    if (feature_cb_it != actin.custom_features.end()) {
        try {
            std::vector<double> actin_f_load_cb_all = load_from_dataset(group_actin, "f_load_cb", dims);
            size_t load_stride = static_cast<size_t>(n_actins);
            size_t load_start = static_cast<size_t>(target_frame) * load_stride;
            for (int i = 0; i < actin.n; ++i) {
                feature_cb_it->second[i] = actin_f_load_cb_all[load_start + i];
            }
        } catch (const H5::Exception&) {
            std::fill(feature_cb_it->second.begin(), feature_cb_it->second.end(), 0.0);
        }
    }
    actin.update_endpoints();
    for (auto& row : actin_actin_bonds) {
        std::fill(row.begin(), row.end(), 0);
    }

    // The bonds dataset is assumed to have dimensions: (n_frames, max_n_actin_bonds, 2)
    std::vector<double> flatActinBondsAll = load_from_dataset(group_actin, "bonds", dims);
    // dims[0] is n_frames, dims[1] is max_n_actin_bonds, dims[2] should be 2.
    size_t bonds_per_frame = dims[1] * dims[2];
    size_t bonds_start = static_cast<size_t>(target_frame) * bonds_per_frame;
    std::vector<double> flatActinBonds(
        flatActinBondsAll.begin() + bonds_start,
        flatActinBondsAll.begin() + bonds_start + bonds_per_frame
    );
    // For each bonded pair in the flat array, update the matrix.
    for (size_t i = 0; i < flatActinBonds.size(); i += 2) {
        int a = static_cast<int>(flatActinBonds[i]);
        int b = static_cast<int>(flatActinBonds[i+1]);
        if (a != -1 && b != -1) {
            actin_actin_bonds[a][b] = 1;
            actin_actin_bonds[b][a] = 1;
        }
    }
    
    std::vector<double> myosin_center_all = load_from_dataset(group_myosin, "center", dims);
    size_t myosin_frame_stride = static_cast<size_t>(n_myosins) * 3;
    size_t myosin_start = static_cast<size_t>(target_frame) * myosin_frame_stride;
    std::vector<double> myosin_center(
        myosin_center_all.begin() + myosin_start,
        myosin_center_all.begin() + myosin_start + myosin_frame_stride
    );
    for (int i = 0; i < myosin.n; i++) {
        myosin.center[i].x = myosin_center[3 * i];
        myosin.center[i].y = myosin_center[3 * i + 1];
        myosin.center[i].z = myosin_center[3 * i + 2];
    }
    std::vector<double> myosin_direction_all = load_from_dataset(group_myosin, "direction", dims);
    std::vector<double> myosin_direction(
        myosin_direction_all.begin() + myosin_start,
        myosin_direction_all.begin() + myosin_start + myosin_frame_stride
    );
    for (int i = 0; i < myosin.n; i++) {
        myosin.direction[i].x = myosin_direction[3 * i];
        myosin.direction[i].y = myosin_direction[3 * i + 1];
        myosin.direction[i].z = myosin_direction[3 * i + 2];
    }
    myosin.update_endpoints();
    return target_frame;
}
