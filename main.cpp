#include "utils/sequencing.h"
#include "utils/workflows.h"


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
    bool removeSparseSequences = true;
    bool createTemporalBuckets = true;
    bool durationInWeeks = false;
    bool durationInMonths = false;


    std::vector<std::string> inputFilePaths;
    inputFilePaths.push_back(fileName);

    char inputFileDelimiter = ',';
    int patIdColumns[1] = {0};
    int phenxIDColumns[1] = {1};
    int dateColumns[1] = {3};
    std::vector<temporalSequence> seq = sequenceWorkflowFromCsVFiles( inputFilePaths, inputFileDelimiter, patIdColumns, phenxIDColumns,
                                  dateColumns,outputDir, description, removeSparseSequences,
                                  sparsity,createTemporalBuckets,durationInWeeks, durationInMonths, removeSparseBuckets,7,16 );
    std::cout<< "Number of sequences: " << seq.size();

//    std::vector<dbMartEntry> dbMart;
//    for (int i = 0; i < 20; ++i) {
//        dbMartEntry entry;
//        entry.patID = i%3;
//        entry.phenID = i%3 + i;
//        entry.date = INT64_MAX/(i+1);
//        dbMart.emplace_back(entry);
//    }
//
//    std::sort(dbMart.begin(), dbMart.end(), dbMartSorter);
//    size_t startPositions[3];
//    startPositions[0] = 0;
//    for (int i = 1; i < 20; ++i) {
//        if(dbMart[i].patID != dbMart[i-1].patID){
//            startPositions[dbMart[i].patID] = i;
//        }
//    }
//    std::cout << extractSequencesFromArray(dbMart, 3, startPositions,outputDir,description,7, 1);
//    std::cout << std::endl;
//    std::map<long, size_t> sequences = summarizeSequences(3, false, outputDir,description);
//    size_t sparsityThreshold = 3 * sparsity;
//    std::cout << "sparsity= " << sparsity << " sparsity threshold: " << sparsityThreshold <<std::endl;
//    for(auto it = sequences.begin(); it != sequences.end();){
//        if(it->second < sparsityThreshold){
//            it = sequences.erase(it);
//        } else{
//            ++it;
//        }
//    }
//    std::vector<temporalSequence> sparseSequences = extractNonSparseSequences(dbMart, 3, startPositions,
//                                                                              sequences, 16, false, false);
    return 0;


}


