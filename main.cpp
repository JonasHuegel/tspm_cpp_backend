

#include "utils/sequencing.h"


#pragma clang diagnostic get push
#pragma ide diagnostic ignored "openmp-use-default-none"

int main(int argc, char *argv[]) {
    omp_set_num_threads(16);
    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_fourtimes_processed.csv";
//    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_debug.csv";
    std::string description = "tspm_test_";
    std::string outputDir = "/home/jonas/CLionProjects/tspm_cpp_backend/out/data/";
    int patIdColumn = 0;
    int phenotypeIdColumn = 1;
    int dateColumn = 3;
    int patientCount = 1168; //TODO add function to determine (extract from last line in file)
    bool createDuration = true;
    bool removeSparseBuckets = true;
    double sparsity = 0.005;

    std::vector<std::string> inputFilePaths;
    inputFilePaths.push_back(fileName);

    char inputFileDelimiter = ',';
    int patIdColumns[1] = {0};
    int phenxIDColumns[1] = {1};
    int dateColumns[1] = {3};
    sequenceWorkflow(createDuration, removeSparseBuckets, inputFilePaths, inputFileDelimiter,
                     outputDir, description, patIdColumns, phenxIDColumns,
                     dateColumns, patientCount, sparsity);


}


