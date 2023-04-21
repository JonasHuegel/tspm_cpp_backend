//
// Created by jonas on 18.04.23.
//
#include "workflows.h"


std::vector<temporalSequence> sequenceWorkflow(std::vector<dbMartEntry> &dbMart, const std::string& outPutDirectory,
                                               const std::string& outputFilePrefix, bool removeSparseSequences,
                                               double sparsity_value, bool createTemporalBuckets, bool durationInWeeks,
                                               bool durationInMonths, bool removeSparseTemporalBuckets, int patIdLength,
                                               int numOfThreads){
    std::vector<size_t> startPositions = extractStartPositions(dbMart);
    size_t numOfPatients = startPositions.size();

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string stringPath = outPutDirectory;
    stringPath.append(oss.str()).append("/");
    std::filesystem::path outputPath = std::filesystem::u8path(stringPath);
    std::filesystem::create_directory(outputPath);

    //===== extract sequence
    std::cout << "extracting transitive sequences" << std::endl;
    size_t numOfSequences  = extractSequencesFromArray(dbMart, numOfPatients,startPositions.data(), outputPath.string(), outputFilePrefix,patIdLength, numOfThreads);
    std::cout << "Number of extracted sequences: " << numOfSequences << std::endl;
    std::map<long, size_t> sequenceCount = summarizeSequences(numOfPatients, false,
                                                              outputPath.string(),outputFilePrefix);
    if(removeSparseSequences) {
        //===== remove sparse sequences
        std::cout << "removing sparse sequences" << std::endl;
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
        std::cout << "Number of unique sparse sequences: " << numOfUniqueSequences << std::endl;
    }

    std::vector<temporalSequence> sequences;
    if(!createTemporalBuckets){
        std::cout << "extracting sequences" << std::endl;
        sequences = extractNonSparseSequences(dbMart, numOfPatients,startPositions.data(), sequenceCount, numOfThreads, durationInWeeks, durationInMonths);

    }else {
        std::cout << "creating temporal sequences" << std::endl;

        sequences = extractTemporalSequences(dbMart, numOfPatients, startPositions.data(), sequenceCount, numOfThreads,
                                             durationInWeeks, durationInMonths, sparsity_value,
                                             removeSparseTemporalBuckets);
    }
    return sequences;

}






std::vector<temporalSequence> sequenceWorkflowFromCsVFiles(const std::vector<std::string>& inputFilePaths, char inputFileDelimiter,
                     int patIDColumns[], int phenxColumns[], int dateColumns[], const std::string& outPutDirectory,
                     const std::string& outputFilePrefix, bool removeSparseSequences, double sparsity_value,
                     bool createTemporalBuckets, bool durationInWeeks, bool durationInMonths,
                     bool removeSparseTemporalBuckets, int patIdLength, int numOfThreads){
    std::vector<dbMartEntry> dbMart;
    for(int i = 0; i < inputFilePaths.size();++i) {
        FILE *csvFilePointer = fopen(inputFilePaths[i].c_str(), "r");
        if (csvFilePointer == nullptr) {
            exit(EXIT_FAILURE);
        }
        std::vector<dbMartEntry> localDBMart = extractDBMartFromCsv(csvFilePointer, patIDColumns[i], phenxColumns[i],
                                                                    dateColumns[i], inputFileDelimiter);
        dbMart.insert(dbMart.end(), localDBMart.begin(), localDBMart.end());
    }

    return  sequenceWorkflow(dbMart, outPutDirectory, outputFilePrefix, removeSparseSequences, sparsity_value,
                             createTemporalBuckets, durationInWeeks, durationInMonths, removeSparseTemporalBuckets,
                             patIdLength, numOfThreads);


}
