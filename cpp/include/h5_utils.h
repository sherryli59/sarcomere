
#ifndef H5_UTILS_H
#define H5_UTILS_H
#include <iostream>
#include <H5Cpp.h>
#include <vector>

#ifndef COMPONENTS_H
#include "components.h"
#endif
using vec = utils::vec;

std::vector<double> flatten_2d_array(double** array, size_t rows, size_t cols) {
    if (array == nullptr) {
        throw std::invalid_argument("Input array is null.");
    }

    std::vector<double> flattened(rows * cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            flattened[i * cols + j] = array[i][j];
        }
    }
    return flattened;
}

std::vector<double> flatten_2d_array(std::vector<vec> array) {
    size_t rows = array.size();
    std::vector<double> flattened(rows * 2);
    for (size_t i = 0; i < rows; ++i) {
            flattened[i * 2] = array[i].x;
            flattened[i * 2 + 1] = array[i].y;
    }
    return flattened;
}

void create_empty_dataset(H5::H5File& file, const std::string& groupName,
     const std::string& datasetName, const std::vector<hsize_t>& initialDims, const std::vector<hsize_t>& maxDims, const std::vector<hsize_t>& chunkDims) {
    try {
        // Open the file or create it if it doesn't exist.
        H5::Group group;
        if (!file.nameExists(groupName)) {
            group = file.createGroup(groupName);
        } else {
            group = file.openGroup(groupName);
        }
        // Create a dataspace for the dataset with unlimited maximum size.
        H5::DataSpace dataspace(initialDims.size(), initialDims.data(), maxDims.data());
        // Create the dataset creation property list and set it to enable chunking.
        H5::DSetCreatPropList prop;
        
        prop.setChunk(chunkDims.size(), chunkDims.data());
        // Create the dataset.
        H5::FloatType datatype(H5::PredType::IEEE_F64LE);
        H5::DataSet dataset = group.createDataSet(datasetName, datatype, dataspace, prop);
    } catch (H5::FileIException& error) {
        error.printErrorStack();
        std::cerr << "Error creating empty dataset: " << datasetName << std::endl;
    } catch (H5::DataSetIException& error) {
        error.printErrorStack();
        std::cerr << "Error creating empty dataset: " << datasetName << std::endl;
    }
}

void append_to_dataset(H5::Group& group, const std::string& datasetName, const std::vector<double> newData, const std::vector<hsize_t>& newDims) {
    try {

        // Open the dataset.
        H5::DataSet dataset = group.openDataSet(datasetName);

        // Get the current dataspace of the dataset.
        H5::DataSpace filespace = dataset.getSpace();

        // Get the current size of the dataset.
        std::vector<hsize_t> currentSize(filespace.getSimpleExtentNdims());
        filespace.getSimpleExtentDims(currentSize.data(), NULL);

        // Ensure new data dimensions match except for the dimension being appended to.
        if (newDims.size() != currentSize.size()) {
            throw std::runtime_error("Dimensions of new data do not match dataset.");
        }

        for (size_t i = 1; i < newDims.size(); ++i) {
            if (newDims[i] != currentSize[i]) {
                throw std::runtime_error("Non-append dimensions do not match dataset.");
            }
        }

        // Extend the dataset to accommodate the new data.
        std::vector<hsize_t> newSize = currentSize;
        newSize[0] += newDims[0];
        dataset.extend(newSize.data());

        // Select the hyperslab in the file where the new data will be written.
        filespace = dataset.getSpace();
        std::vector<hsize_t> offset(newSize.size(), 0);
        offset[0] = currentSize[0];
        filespace.selectHyperslab(H5S_SELECT_SET, newDims.data(), offset.data());

        // Define the memory space.
        H5::DataSpace memspace(newDims.size(), newDims.data());

        // Write the data to the hyperslab.
        dataset.write(newData.data(), H5::PredType::IEEE_F64LE, memspace, filespace);
    
    } catch (H5::FileIException& error) {
        error.printErrorStack();
        std::cerr << "Error appending data to dataset in group: " << datasetName << std::endl;
    } catch (H5::DataSetIException& error) {
        error.printErrorStack();
        std::cerr << "Error appending data to dataset in group: " << datasetName << std::endl;
    } catch (H5::DataSpaceIException& error) {
        error.printErrorStack();
        std::cerr << "Error appending data to dataset in group: " << datasetName << std::endl;
    } catch (std::runtime_error& error) {
        std::cerr << "Runtime error: " << error.what() << std::endl;
    }
}


void create_file(std::string& filename,  Filament& actin, Myosin& myosin){
    std::remove(filename.c_str());
    H5::H5File file(filename, H5F_ACC_TRUNC);
    hsize_t n_actins = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    std::vector<hsize_t> initialDims = {0, n_actins, 2};
    std::vector<hsize_t> maxDims = {H5S_UNLIMITED, n_actins, 2};
    std::vector<hsize_t> chunkDims = {10, n_actins, 2};
    create_empty_dataset(file, "/actin", "center", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "velocity", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "force", initialDims, maxDims, chunkDims);
    initialDims = {0, n_actins, 1};
    maxDims = {H5S_UNLIMITED, n_actins, 1};
    chunkDims = {10, n_actins, 1};
    create_empty_dataset(file, "/actin", "theta", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "angular_force", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/actin", "cb_strength", initialDims, maxDims, chunkDims);

    // Loop through custom features and create datasets for each
    for (const auto& feature : actin.custom_features) {
        create_empty_dataset(file, "/actin", feature.first, initialDims, maxDims, chunkDims);
    }

    initialDims = {0, n_myosins, 2};
    maxDims = {H5S_UNLIMITED, n_myosins, 2};
    chunkDims = {10, n_myosins, 2};
    create_empty_dataset(file, "/myosin", "center", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "velocity", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "force", initialDims, maxDims, chunkDims);
    initialDims = {0, n_myosins, 1};
    maxDims = {H5S_UNLIMITED, n_myosins, 1};
    chunkDims = {10, n_myosins, 1};
    create_empty_dataset(file, "/myosin", "theta", initialDims, maxDims, chunkDims);
    create_empty_dataset(file, "/myosin", "angular_force", initialDims, maxDims, chunkDims);
    for (const auto& feature : myosin.custom_features) {
        create_empty_dataset(file, "/myosin", feature.first, initialDims, maxDims, chunkDims);
    }


}

void append_to_file(std::string& filename, Filament& actin, Myosin& myosin, double& total_energy){
    H5::H5File file(filename, H5F_ACC_RDWR);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));
    hsize_t n_actins = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    std::vector<double> flattened_actin_center = flatten_2d_array(actin.center);
    append_to_dataset(group_actin, "center", flattened_actin_center, {1, n_actins, 2});

    std::vector<double> flattened_actin_velocity = flatten_2d_array(actin.velocity);
    append_to_dataset(group_actin, "velocity", flattened_actin_velocity, {1, n_actins, 2});

    std::vector<double> flattened_actin_force = flatten_2d_array(actin.force);
    append_to_dataset(group_actin, "force", flattened_actin_force, {1, n_actins, 2});

    append_to_dataset(group_actin, "theta", actin.theta, {1, n_actins, 1});
    append_to_dataset(group_actin, "angular_force", actin.angular_force, {1, n_actins, 1});
    
    append_to_dataset(group_actin, "cb_strength", actin.cb_strength, {1, n_actins, 1});

    // Loop through custom features and append data for each
    for (const auto& feature : actin.custom_features) {
        append_to_dataset(group_actin, feature.first, feature.second, {1, n_actins, 1});
    }


    std::vector<double> flattened_myosin_center = flatten_2d_array(myosin.center);
    append_to_dataset(group_myosin, "center", flattened_myosin_center, {1, n_myosins, 2});

    std::vector<double> flattened_myosin_velocity = flatten_2d_array(myosin.velocity);
    append_to_dataset(group_myosin, "velocity", flattened_myosin_velocity, {1, n_myosins, 2});

    std::vector<double> flattened_myosin_force = flatten_2d_array(myosin.force);
    append_to_dataset(group_myosin, "force", flattened_myosin_force, {1, n_myosins, 2});

    append_to_dataset(group_myosin, "theta", myosin.theta, {1, n_myosins, 1});
    append_to_dataset(group_myosin, "angular_force", myosin.angular_force, {1, n_myosins, 1});

    for (const auto& feature : myosin.custom_features) {
        append_to_dataset(group_myosin, feature.first, feature.second, {1, n_myosins, 1});
    }

}

std::vector<double> load_from_dataset(H5::Group& group, const std::string& datasetName, std::vector<hsize_t>& dims){
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


void load_from_file(std::string& filename, Filament& actin, Myosin& myosin,  double& total_energy){
    H5::H5File file(filename, H5F_ACC_RDONLY);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));
    hsize_t n_actins = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    std::vector<hsize_t> dims;
    std::vector<double> actin_center = load_from_dataset(group_actin, "center", dims);
    //load actin
    actin_center = std::vector<double>(actin_center.end() - 2*n_actins, actin_center.end());
    for (int i = 0; i < actin.n; i++){
        actin.center[i].x = actin_center[2*i];
        actin.center[i].y = actin_center[2*i+1];
    }
    std::vector<double> actin_theta = load_from_dataset(group_actin, "theta", dims);
    actin_theta = std::vector<double>(actin_theta.end() - n_actins, actin_theta.end());
    for (int i = 0; i < actin.n; i++){
        actin.theta[i] = actin_theta[i];
    }
    actin.update_endpoints();
    //load myosin
    std::vector<double> myosin_center = load_from_dataset(group_myosin, "center", dims);
    myosin_center = std::vector<double>(myosin_center.end() - 2*n_myosins, myosin_center.end());
    for (int i = 0; i < myosin.n; i++){
        myosin.center[i].x = myosin_center[2*i];
        myosin.center[i].y = myosin_center[2*i + 1];
    }
    std::vector<double> myosin_theta = load_from_dataset(group_myosin, "theta", dims);
    myosin_theta = std::vector<double>(myosin_theta.end() - n_myosins, myosin_theta.end());
    for (int i = 0; i < myosin.n; i++){
        myosin.theta[i] = myosin_theta[i];
    }
    myosin.update_endpoints();
    // //load energy
    // std::vector<double> energy = load_from_dataset(group_energy, "total_energy", dims);
    // //last element of energy is the total energy
    // total_energy = energy.back();
}
#endif