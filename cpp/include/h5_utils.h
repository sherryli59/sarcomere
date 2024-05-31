
#ifndef H5_UTILS_H
#define H5_UTILS_H
#include <iostream>
#include <H5Cpp.h>
#include <vector>



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
        prop.setChunk(initialDims.size(), chunkDims.data());

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

void append_to_dataset(H5::Group& group, const std::string& datasetName, const double* newData, const std::vector<hsize_t>& newDims) {
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
        dataset.write(newData, H5::PredType::NATIVE_DOUBLE, memspace, filespace);
    
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


void create_file(std::string& filename, int& n_actins0, int& n_myosins0, int& n_alpha_actinins0){
    std::remove(filename.c_str());

    H5::H5File file(filename, H5F_ACC_TRUNC);
    hsize_t n_actins = static_cast<hsize_t>(n_actins0);
    hsize_t n_myosins = static_cast<hsize_t>(n_myosins0);
    hsize_t n_alpha_actinins = static_cast<hsize_t>(n_alpha_actinins0);

    std::vector<hsize_t> initialDims = {0, n_actins, 2};
    std::vector<hsize_t> maxDims = {H5S_UNLIMITED, n_actins, 2};
    std::vector<hsize_t> chunkDims = {10, n_actins, 2};
    create_empty_dataset(file, "/actin", "xs", initialDims, maxDims, chunkDims);
    initialDims = {0, n_actins, 1};
    maxDims = {H5S_UNLIMITED, n_actins, 1};
    chunkDims = {10, n_actins, 1};
    create_empty_dataset(file, "/actin", "thetas", initialDims, maxDims, chunkDims);
    initialDims = {0, n_myosins, 2};
    maxDims = {H5S_UNLIMITED, n_myosins, 2};
    chunkDims = {10, n_myosins, 2};
    create_empty_dataset(file, "/myosin", "xs", initialDims, maxDims, chunkDims);
    initialDims = {0, n_myosins, 1};
    maxDims = {H5S_UNLIMITED, n_myosins, 1};
    chunkDims = {10, n_myosins, 1};
    create_empty_dataset(file, "/myosin", "thetas", initialDims, maxDims, chunkDims);
    initialDims = {0, n_alpha_actinins, 2};
    maxDims = {H5S_UNLIMITED, n_alpha_actinins, 2};
    chunkDims = {10, n_alpha_actinins, 2};
    create_empty_dataset(file, "/alpha_actinin", "xs", initialDims, maxDims, chunkDims);
    initialDims = {0, 1};
    maxDims = {H5S_UNLIMITED, 1};
    chunkDims = {10, 1};
    create_empty_dataset(file, "/energy", "total_energy", initialDims, maxDims, chunkDims);
}

void append_to_file(std::string& filename, Filament& actin, Myosin& myosin, AlphaActinin& alpha_actinin, double& total_energy){
    H5::H5File file(filename, H5F_ACC_RDWR);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));
    H5::Group group_alpha_actinin(file.openGroup("/alpha_actinin"));
    H5::Group group_energy(file.openGroup("/energy"));
    hsize_t n_actins = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    hsize_t n_alpha_actinins = static_cast<hsize_t>(alpha_actinin.n);
    std::vector<double> flattened_actin_xs = flatten_2d_array(actin.xs, n_actins, 2);

    append_to_dataset(group_actin, "xs", flattened_actin_xs.data(), {1, n_actins, 2});
    append_to_dataset(group_actin, "thetas", actin.thetas, {1, n_actins, 1});
    std::vector<double> flattened_myosin_xs = flatten_2d_array(myosin.xs, n_myosins, 2);
    append_to_dataset(group_myosin, "xs", flattened_myosin_xs.data(), {1, n_myosins, 2});
    append_to_dataset(group_myosin, "thetas", myosin.thetas, {1, n_myosins, 1});
    std::vector<double> flattened_alpha_actinin_xs = flatten_2d_array(alpha_actinin.xs, n_alpha_actinins, 2);
    append_to_dataset(group_alpha_actinin, "xs", flattened_alpha_actinin_xs.data(), {1, n_alpha_actinins, 2});
    append_to_dataset(group_energy, "total_energy", &total_energy, {1, 1});
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


void load_from_file(std::string& filename, Filament& actin, Myosin& myosin, AlphaActinin& alpha_actinin, double& total_energy){
    H5::H5File file(filename, H5F_ACC_RDONLY);
    H5::Group group_actin(file.openGroup("/actin"));
    H5::Group group_myosin(file.openGroup("/myosin"));
    H5::Group group_alpha_actinin(file.openGroup("/alpha_actinin"));
    H5::Group group_energy(file.openGroup("/energy"));
    hsize_t n_actins = static_cast<hsize_t>(actin.n);
    hsize_t n_myosins = static_cast<hsize_t>(myosin.n);
    hsize_t n_alpha_actinins = static_cast<hsize_t>(alpha_actinin.n);
    std::vector<hsize_t> dims;
    std::vector<double> actin_xs = load_from_dataset(group_actin, "xs", dims);
    //load actin
    actin_xs = std::vector<double>(actin_xs.end() - 2*n_actins, actin_xs.end());
    for (int i = 0; i < actin.n; i++){
        actin.xs[i][0] = actin_xs[2*i];
        actin.xs[i][1] = actin_xs[2*i + 1];
    }
    std::vector<double> actin_thetas = load_from_dataset(group_actin, "thetas", dims);
    actin_thetas = std::vector<double>(actin_thetas.end() - n_actins, actin_thetas.end());
    for (int i = 0; i < actin.n; i++){
        actin.thetas[i] = actin_thetas[i];
    }
    //load myosin
    std::vector<double> myosin_xs = load_from_dataset(group_myosin, "xs", dims);
    myosin_xs = std::vector<double>(myosin_xs.end() - 2*n_myosins, myosin_xs.end());
    for (int i = 0; i < myosin.n; i++){
        myosin.xs[i][0] = myosin_xs[2*i];
        myosin.xs[i][1] = myosin_xs[2*i + 1];
    }
    std::vector<double> myosin_thetas = load_from_dataset(group_myosin, "thetas", dims);
    myosin_thetas = std::vector<double>(myosin_thetas.end() - n_myosins, myosin_thetas.end());
    for (int i = 0; i < myosin.n; i++){
        myosin.thetas[i] = myosin_thetas[i];
    }
    //load alpha_actinin
    std::vector<double> alpha_actinin_xs = load_from_dataset(group_alpha_actinin, "xs", dims);
    alpha_actinin_xs = std::vector<double>(alpha_actinin_xs.end() - 2*n_alpha_actinins, alpha_actinin_xs.end());
    for (int i = 0; i < alpha_actinin.n; i++){
        alpha_actinin.xs[i][0] = alpha_actinin_xs[2*i];
        alpha_actinin.xs[i][1] = alpha_actinin_xs[2*i + 1];
    }
    //load energy
    std::vector<double> energy = load_from_dataset(group_energy, "total_energy", dims);
    //last element of energy is the total energy
    total_energy = energy.back();
}
#endif