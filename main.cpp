#include "utils/sequencing.h"
#include "utils/workflows.h"




int main(int argc, char *argv[]) {
    omp_set_num_threads(16);
//    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_fourtimes_processed.csv";
    std::string fileName = "/home/jonas/dbmart_num.test.csv";
    std::string description = "tspm_test_";
    std::string outputDir = "/home/jonas/CLionProjects/tspm_cpp_backend/out/data/";

    bool removeSparseBuckets = true;
    double sparsity = 0.05;
    bool removeSparseSequences = true;
    bool createTemporalBuckets = false;
    double durationPeriods = DURATION_IN_MONTHS;
    unsigned int daysForCoOoccurence = 14;
    bool durationSparsity = true;
    double durationSparsityValue = 0.05;
    bool storeSeqDuringCreation = false;

    std::vector<std::string> inputFilePaths;
    inputFilePaths.push_back(fileName);

    char inputFileDelimiter = ',';
    int patIdColumns[1] = {0};
    int phenxIDColumns[1] = {1};
    int dateColumns[1] = {2};
//    std::vector<temporalSequence> seq = sequenceWorkflowFromCsVFiles( inputFilePaths, inputFileDelimiter, patIdColumns, phenxIDColumns,
//                                  dateColumns, storeSeqDuringCreation, outputDir, description, removeSparseSequences,
//                                  sparsity,createTemporalBuckets, durationPeriods, daysForCoOoccurence,
//                                  durationSparsity, durationSparsityValue, removeSparseBuckets,7,1);
//    std::cout<< "Number of sequences: " << seq.size() << std::endl;

    std::vector<dbMartEntry> dbMart;
    for (int i = 0; i < 20; ++i) {
        dbMartEntry entry;
        entry.patID = i%3;
        entry.phenID = i%3 + i;
        entry.date = INT64_MAX/(i+1);
        dbMart.emplace_back(entry);
    }

    std::sort(dbMart.begin(), dbMart.end(), dbMartSorter);
    size_t startPositions[3];
    startPositions[0] = 0;
    for (int i = 1; i < 20; ++i) {
        if(dbMart[i].patID != dbMart[i-1].patID){
            startPositions[dbMart[i].patID] = i;
        }
    }
    durationSparsityValue = 0;
    durationSparsity = false;
    std::vector<temporalSequence> nonSparseSequences = sequenceWorkflow(dbMart,
                     storeSeqDuringCreation,
                     outputDir,
                     description,
                     removeSparseSequences,
                     sparsity,
                     createTemporalBuckets,
                     durationPeriods,
                     daysForCoOoccurence,
                     durationSparsity,
                     durationSparsityValue,
                     removeSparseBuckets,
                     7,
                     16);

    std::vector<unsigned int> phenx;
    phenx.emplace_back(3);
    int t = 16;
    auto  a = extractEndPhenxWithGivenStartPhenx(nonSparseSequences,3,0,7,phenx,t);
    extractSequencesWithEnd(nonSparseSequences,0,7, a, t);
//    std::cout << extractSequencesFromArray(dbMart, 3, startPositions,outputDir,description,7, 1);
//    std::cout << std::endl;
//    std::map<std::uint64_t, size_t> sequences = summarizeSequencesFromFiles(3, false, outputDir,description);
//    size_t sparsityThreshold = 3 * sparsity;
//    std::cout << "sparsity= " << sparsity << " sparsity threshold: " << sparsityThreshold <<std::endl;
//    for(auto it = sequences.begin(); it != sequences.end();){
//        if(it->second < sparsityThreshold){
//            it = sequences.erase(it);
//        } else{
//            ++it;
//        }
//    }
//    std::vector<temporalSequence> nonSparseSequences = extractNonSparseSequences(dbMart, 3, startPositions,
//                                                                              sequences, 16, false, false);
    return 0;


}


