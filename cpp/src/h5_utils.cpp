#include "h5_utils.h"
#include <stdexcept>
#include <cmath>
#include "omp.h"

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

void append_to_dataset(H5::Group& group, const std::string& datasetName,
                       const std::vector<double> newData,
                       const std::vector<hsize_t>& newDims)
{
    try {
        // Open the dataset.
        H5::DataSet dataset = group.openDataSet(datasetName);

        // Get the current dataspace and size of the dataset.
        H5::DataSpace filespace = dataset.getSpace();
        std::vector<hsize_t> currentSize(filespace.getSimpleExtentNdims());
        filespace.getSimpleExtentDims(currentSize.data(), NULL);

        // Check that the new data dimensions match (except for the first dimension).
        if (newDims.size() != currentSize.size()) {
            throw std::runtime_error("Dimensions of new data do not match dataset.");
        }
        for (size_t i = 1; i < newDims.size(); ++i) {
            if (newDims[i] != currentSize[i]) {
                printf("datasetName: %s\n", datasetName.c_str());
                printf("newDims[%zu]: %zu, currentSize[%zu]: %zu\n", i, newDims[i], i, currentSize[i]);
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
    create_empty_dataset(file, "/actin", "cb_strength", initialDims, maxDims, chunkDims);
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

    append_to_dataset(group_actin, "cb_strength", actin.cb_strength, {1, n_actins, 1});
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
        }
    }

    // Now that each flattened vector has a size corresponding to (max_n * 2),
    // the new dimensions for appending one frame are {1, max_n, 2}.
    hsize_t newActinDims[3] = {1, static_cast<hsize_t>(max_n_actin_bonds), 2};
    hsize_t newMyosinDims[3] = {1, static_cast<hsize_t>(max_n_myosin_bonds), 2};
    hsize_t newAmDims[3] = {1, static_cast<hsize_t>(max_n_am_bonds), 2};

    append_to_dataset(group_actin, "bonds", flatActinBonds, { newActinDims[0], newActinDims[1], newActinDims[2] });
    append_to_dataset(group_myosin, "bonds", flatMyosinBonds, { newMyosinDims[0], newMyosinDims[1], newMyosinDims[2] });
    append_to_dataset(group_am, "bonds", flatActinMyosinBonds, { newAmDims[0], newAmDims[1], newAmDims[2] });
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

void load_from_file(std::string& filename, Filament& actin, Myosin& myosin,
                    std::vector<std::vector<int>>& actin_actin_bonds, int& n_frames)
{
    H5::H5File file(filename, H5F_ACC_RDONLY);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));

    hsize_t n_actins  = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);

    std::vector<hsize_t> dims;
    std::vector<double> actin_center = load_from_dataset(group_actin, "center", dims);
    n_frames = static_cast<int>(dims[0]);

    // Extract the most recent actin center data.
    actin_center = std::vector<double>(actin_center.end() - 3 * n_actins, actin_center.end());
    for (int i = 0; i < actin.n; i++) {
        actin.center[i].x = actin_center[3 * i];
        actin.center[i].y = actin_center[3 * i + 1];
        actin.center[i].z = actin_center[3 * i + 2];
    }
    std::vector<double> actin_direction = load_from_dataset(group_actin, "direction", dims);
    actin_direction = std::vector<double>(actin_direction.end() - 3 * n_actins, actin_direction.end());
    for (int i = 0; i < actin.n; i++) {
        actin.direction[i].x = actin_direction[3 * i];
        actin.direction[i].y = actin_direction[3 * i + 1];
        actin.direction[i].z = actin_direction[3 * i + 2];
    }
    actin.update_endpoints();

    // The bonds dataset is assumed to have dimensions: (n_frames, max_n_actin_bonds, 2)
    std::vector<double> flatActinBonds = load_from_dataset(group_actin, "bonds", dims);
    // dims[0] is n_frames, dims[1] is max_n_actin_bonds, dims[2] should be 2.
    // Extract the most recent frame:
    size_t bonds_per_frame = dims[1] * dims[2];
    flatActinBonds = std::vector<double>(flatActinBonds.end() - bonds_per_frame, flatActinBonds.end());
    // For each bonded pair in the flat array, update the matrix.
    for (size_t i = 0; i < flatActinBonds.size(); i += 2) {
        int a = static_cast<int>(flatActinBonds[i]);
        int b = static_cast<int>(flatActinBonds[i+1]);
        if (a != -1 && b != -1) {
            actin_actin_bonds[a][b] = 1;
            actin_actin_bonds[b][a] = 1;
        }
    }
    
    std::vector<double> myosin_center = load_from_dataset(group_myosin, "center", dims);
    myosin_center = std::vector<double>(myosin_center.end() - 3 * n_myosins, myosin_center.end());
    for (int i = 0; i < myosin.n; i++) {
        myosin.center[i].x = myosin_center[3 * i];
        myosin.center[i].y = myosin_center[3 * i + 1];
        myosin.center[i].z = myosin_center[3 * i + 2];
    }
    std::vector<double> myosin_direction = load_from_dataset(group_myosin, "direction", dims);
    myosin_direction = std::vector<double>(myosin_direction.end() - 3 * n_myosins, myosin_direction.end());
    for (int i = 0; i < myosin.n; i++) {
        myosin.direction[i].x = myosin_direction[3 * i];
        myosin.direction[i].y = myosin_direction[3 * i + 1];
        myosin.direction[i].z = myosin_direction[3 * i + 2];
    }
    myosin.update_endpoints();
}
