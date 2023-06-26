#include "../tspmPlus.h"


int main(int argc, char *argv[]) {

    std::cout << sizeof (unsigned int) << " " << sizeof (std::uint64_t)<< std::endl;
    omp_set_num_threads(11);
//    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_fourtimes_processed.csv";
    std::string fileName = "/home/jonas/dbmart_num.test.csv";
    std::string description = "tspm_test_";
    std::string outputDir = "/home/jonas/CLionProjects/tspm_cpp_backend/out/data/";

    bool removeSparseBuckets = true;
    double sparsity = 0.05;
    bool removeSparseSequences = true;
    bool createTemporalBuckets = false;
    double durationPeriods = tspm::DURATION_IN_MONTHS;
    unsigned int coOccurrence = 14;
    bool durationSparsity = false;
    double durationSparsityValue = 0.05;
    bool storeSeqDuringCreation = false;
    unsigned int numOfThreads = 16;
    std::vector<unsigned int> thresholds = {0,1,3};

    std::vector<std::string> inputFilePaths;
    inputFilePaths.push_back(fileName);

    char inputFileDelimiter = ',';
    int patIdColumns[1] = {0};
    int phenxIDColumns[1] = {1};
    int dateColumns[1] = {2};
    std::vector<tspm::temporalSequence> seq = tspm::sequenceWorkflowFromCsVFiles( inputFilePaths, inputFileDelimiter, patIdColumns, phenxIDColumns,
                                  dateColumns, storeSeqDuringCreation, outputDir, description, removeSparseSequences,
                                  sparsity,createTemporalBuckets, durationPeriods, coOccurrence,
                                  durationSparsity, durationSparsityValue, removeSparseBuckets,7,16);
    std::cout<< "Number of sequences: " << seq.size() << std::endl;

    auto sum = tspm::summarizeSequencesAsVector(seq,true,thresholds,numOfThreads);

    std::vector<tspm::dbMartEntry> dbMart;
    for (int i = 0; i < 20; ++i) {
        tspm::dbMartEntry entry;
        entry.patID = i%3;
        entry.phenID = i%3 + i;
        entry.date = INT64_MAX/(i+1);
        dbMart.emplace_back(entry);
    }

    std::sort(dbMart.begin(), dbMart.end(), tspm::dbMartSorter);
    size_t startPositions[3];
    startPositions[0] = 0;
    for (int i = 1; i < 20; ++i) {
        if(dbMart[i].patID != dbMart[i-1].patID){
            startPositions[dbMart[i].patID] = i;
        }
    }
    durationSparsityValue = 0;
    durationSparsity = false;
    std::vector<tspm::temporalSequence> nonSparseSequences = tspm::sequenceWorkflow(dbMart,
                                                                        storeSeqDuringCreation,
                                                                        outputDir,
                                                                        description,
                                                                        removeSparseSequences,
                                                                        sparsity,
                                                                        createTemporalBuckets,
                                                                        durationPeriods,
                                                                        coOccurrence,
                                                                        durationSparsity,
                                                                        durationSparsityValue,
                                                                        removeSparseBuckets,
                                                                        7,
                                                                        16);
    std::cout << nonSparseSequences.size() << std::endl;

    std::vector<unsigned int> phenx;
    phenx.emplace_back(3);
    int t = 16;
    auto  endPhenx = tspm::extractEndPhenxWithGivenStartPhenx(nonSparseSequences,3,0,7,phenx,t);
    auto sequencesOfInterest = tspm::extractSequencesWithEnd(nonSparseSequences,0,7, endPhenx, t);
//    std::cout << tspm::extractSequencesFromArray(dbMart, 3, startPositions,outputDir,description,7, 1);
//    std::cout << std::endl;
//    std::map<std::uint64_t, size_t> sequences = tspm::summarizeSequencesFromFiles(3, false, outputDir,description);
//    size_t sparsityThreshold = 3 * sparsity;
//    std::cout << "sparsity= " << sparsity << " sparsity threshold: " << sparsityThreshold <<std::endl;
//    for(auto it = sequences.begin(); it != sequences.end();){
//        if(it->second < sparsityThreshold){
//            it = sequences.erase(it);
//        } else{
//            ++it;
//        }
//    }
//    std::vector<tspm::temporalSequence> nonSparseSequences = tspm::extractNonSparseSequences(dbMart, 3, startPositions,
//                                                                              sequences, 16, false, false);
    return 0;


}


