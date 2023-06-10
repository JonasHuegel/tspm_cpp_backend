//
// Created by jonas on 18.04.23.
//
#include "workflows.h"


std::filesystem::path createOutputFilePath(const std::string &basicString);




std::vector<temporalSequence> sequenceWorkflow(std::vector<dbMartEntry> &dbMart, bool storeSeqDuringCreation, const std::string &outPutDirectory,
                                               const std::string &outputFilePrefix, bool removeSparseSequences,
                                               double sparsity_value, bool createTemporalBuckets, double durationPeriods,
                                               unsigned int daysForCoOccurrence, bool durationSparsity,
                                               double durationSparsityValue, bool removeSparseTemporalBuckets,
                                               unsigned int patIdLength, unsigned int numOfThreads) {
    std::vector<size_t> startPositions = extractStartPositions(dbMart);
    size_t numOfPatients = startPositions.size();

    std::map<std::int64_t, size_t> sequenceCount;
    //===== extract sequence
    if(storeSeqDuringCreation) {
        std::filesystem::path outputPath = createOutputFilePath(outPutDirectory);
        std::cout << "extracting transitive sequences" << std::endl;
        size_t numOfSequences = extractSequencesFromArray(dbMart, numOfPatients, startPositions.data(),
                                                          outputPath.string(), outputFilePrefix, patIdLength,
                                                          numOfThreads);
        std::cout << "Number of extracted sequences: " << numOfSequences << std::endl;
        sequenceCount = summarizeSequencesFromFiles(outputPath.string(), outputFilePrefix,
                                                    numOfPatients, false, patIdLength, 0, numOfThreads);
    }else{
        sequenceCount = summarizeSequencesFromDbMart(dbMart, startPositions, numOfThreads);
    }
    std::cout << "Number of overall unique sequences: " << sequenceCount.size() <<std::endl;
    size_t sum =0;
    for(auto entry : sequenceCount){
        sum += entry.second;
    }
    std::cout << "Sum of overall unique sequences: " << sum << std::endl;

    if(removeSparseSequences) {
        //===== remove sparse sequences
        std::cout << "determine sparse sequences" << std::endl;
        size_t sparsityThreshold = numOfPatients * sparsity_value;
        std::cout << "sparsity= " << sparsity_value << " sparsity threshold: " << sparsityThreshold << std::endl;
        for (auto it = sequenceCount.begin(); it != sequenceCount.end();) {
            if (it->second < sparsityThreshold) {
                it = sequenceCount.erase(it);
            } else {
                ++it;
            }
        }
        size_t numOfUniqueSequences = sequenceCount.size();
        std::cout << "Number of unique non-sparse sequences: " << numOfUniqueSequences << std::endl;
    }

    std::vector<temporalSequence> sequences;
    if(!createTemporalBuckets){
        if(durationSparsity){
            std::cout << "extracting sequences with non-sparse duration" <<std::endl;
            std::vector<temporalSequence> nonSparseSequences;
            nonSparseSequences = extractNonSparseSequences(dbMart, numOfPatients, startPositions.data(),
                                                           sequenceCount, numOfThreads,
                                                           DURATION_IN_MONTHS, 14);
            std::cout << "extracted non-sparse sequences: "<<nonSparseSequences.size()  << "! Removing sparse durations" <<std::endl;
            sequences = extractMonthlySequences(nonSparseSequences, durationSparsity,
                                                durationSparsityValue,numOfPatients,numOfThreads);

        }else {
            std::cout << "extracting (non-sparse) sequences" << std::endl;
            sequences = extractNonSparseSequences(dbMart, numOfPatients, startPositions.data(), sequenceCount,
                                                  numOfThreads, DURATION_IN_MONTHS,14);
        }


    }else {
        std::cout << "creating sequences with temporal buckets" << std::endl;

        sequences = extractTemporalBuckets(dbMart, numOfPatients, startPositions.data(), sequenceCount, numOfThreads,
                                           durationPeriods, daysForCoOccurrence, sparsity_value,
                                           removeSparseTemporalBuckets);
    }
    return sequences;

}

std::vector<temporalSequence> sequenceWorkflowFromCsVFiles(const std::vector<std::string>& inputFilePaths, char inputFileDelimiter,
                                                           int patIDColumns[], int phenxColumns[], int dateColumns[], bool storeSeqDuringCreation, const std::string& outPutDirectory,
                                                           const std::string& outputFilePrefix, bool removeSparseSequences, double sparsity_value,
                                                           bool createTemporalBuckets, double durationPeriods,
                                                           unsigned int daysForCoOccurrence, bool durationSparsity,
                                                           double durationSparsityValue, bool removeSparseTemporalBuckets, unsigned int patIdLength, unsigned int numOfThreads){
    std::vector<dbMartEntry> dbMart;
    for(int i = 0; i < inputFilePaths.size();++i) {
        FILE *csvFilePointer = fopen(inputFilePaths[i].c_str(), "r");
        if (csvFilePointer == nullptr) {
            return {};
        }
        std::vector<dbMartEntry> localDBMart = extractDBMartFromCsv(csvFilePointer, patIDColumns[i], phenxColumns[i],
                                                                    dateColumns[i], inputFileDelimiter);
        dbMart.insert(dbMart.end(), localDBMart.begin(), localDBMart.end());
    }

    return sequenceWorkflow(dbMart, storeSeqDuringCreation, outPutDirectory, outputFilePrefix, removeSparseSequences, sparsity_value,
                            createTemporalBuckets, durationPeriods, daysForCoOccurrence, durationSparsity, durationSparsityValue,
                            removeSparseTemporalBuckets, patIdLength, numOfThreads);


}
