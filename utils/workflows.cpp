//
// Created by jonas on 18.04.23.
//
#include "workflows.h"

namespace tspm{

    std::vector<temporalSequence> sequenceWorkflow(std::vector<dbMartEntry> &dbMart, bool storeSeqDuringCreation, const std::string &outPutDirectory,
                                                   const std::string &outputFilePrefix, bool removeSparseSequences,
                                                   double sparsity_value, bool createTemporalBuckets, double durationPeriods,
                                                   unsigned int daysForCoOccurrence, bool durationSparsity,
                                                   double durationSparsityValue, bool removeSparseTemporalBuckets,
                                                   unsigned int patIdLength, unsigned int numOfThreads, unsigned int phenxIdLength) {
        std::vector<size_t> startPositions = extractStartPositions(dbMart);
        std::vector<temporalSequence> sequences;
        //===== extract sequence
        if(storeSeqDuringCreation) {
            std::vector<temporalSequence> allSequences;
            std::filesystem::path outputPath = createOutputFilePath(outPutDirectory);
            std::cout << "extracting transitive sequences" << std::endl;
            size_t numOfSequences = writeSequencesFromArrayToFile(dbMart, startPositions,
                                                              outputPath.string(), outputFilePrefix, patIdLength,
                                                              numOfThreads);
            std::cout << "Number of extracted sequences: " << numOfSequences << std::endl;
            allSequences = readSequencesFromFiles(outputPath.string(), outputFilePrefix,
                                                  startPositions.size(), false, patIdLength, 0, numOfThreads);
            if(removeSparseSequences) {
                allSequences = tspm::removeSparseSequences(allSequences,startPositions.size(), sparsity_value, numOfThreads);
                std::cout << "extracted non-sparse sequences: "<<allSequences.size()  << "!" << std::endl;

            }
            sequences = allSequences;
        }else if(removeSparseSequences){
                sequences = extractNonSparseSequences(dbMart,startPositions,sparsity_value,numOfThreads,
                                                      durationPeriods, daysForCoOccurrence,phenxIdLength);
             std::cout << "extracted non-sparse sequences: "<<sequences.size()  << "!" << std::endl;
        }else{
            sequences = extractSparseSequences(dbMart,startPositions, numOfThreads, durationPeriods, daysForCoOccurrence);
            std::cout << "extracted sparse sequences: "<<sequences.size()  << "!" << std::endl;
        }
        if(!createTemporalBuckets){
            if(durationSparsity){
                std::cout << "! Removing sparse durations" <<std::endl;
                sequences = extractMonthlySequences(sequences, durationSparsity,
                                                    durationSparsityValue,startPositions.size(),numOfThreads);
            }
        }else {
            std::cout << "creating sequences with temporal buckets" << std::endl;
            sequences = extractTemporalBuckets(sequences, startPositions.size(), numOfThreads,
                                               durationPeriods, daysForCoOccurrence, sparsity_value,
                                               removeSparseTemporalBuckets);
        }
        return sequences;
    }


//    std::vector<temporalSequence> sequenceWorkflow(std::vector<dbMartEntry> &dbMart, bool storeSeqDuringCreation, const std::string &outPutDirectory,
//                                                   const std::string &outputFilePrefix, bool removeSparseSequences,
//                                                   double sparsity_value, bool createTemporalBuckets, double durationPeriods,
//                                                   unsigned int daysForCoOccurrence, bool durationSparsity,
//                                                   double durationSparsityValue, bool removeSparseTemporalBuckets,
//                                                   unsigned int patIdLength, unsigned int numOfThreads) {
//        std::vector<size_t> startPositions = extractStartPositions(dbMart);
//        std::map<std::int64_t, size_t> sequenceCount;
//        //===== extract sequence
//        if(storeSeqDuringCreation) {
//            std::filesystem::path outputPath = createOutputFilePath(outPutDirectory);
//            std::cout << "extracting transitive sequences" << std::endl;
//            size_t numOfSequences = extractSequencesFromArray(dbMart, startPositions,
//                                                              outputPath.string(), outputFilePrefix, patIdLength,
//                                                              numOfThreads);
//            std::cout << "Number of extracted sequences: " << numOfSequences << std::endl;
//            if(removeSparseSequences) {
//                sequenceCount = summarizeSequencesFromFiles(outputPath.string(), outputFilePrefix,
//                                                            startPositions.size(), false, patIdLength, 0, numOfThreads);
//
//            }
//            return {};
//        }else if(removeSparseSequences){
//            sequenceCount = summarizeSequencesFromDbMart(dbMart, startPositions, numOfThreads);
//        }else{
//            return extractSparseSequences(dbMart,startPositions, numOfThreads, durationPeriods, daysForCoOccurrence);
//        }
//
//
//        std::vector<temporalSequence> sequences;
//        if(!createTemporalBuckets){
//            if(durationSparsity){
//                std::cout << "extracting sequences with non-sparse duration" <<std::endl;
//                std::vector<temporalSequence> nonSparseSequences;
//                nonSparseSequences = extractNonSparseSequences(dbMart, startPositions,
//                                                               sequenceCount, numOfThreads,
//                                                               DURATION_IN_MONTHS, 14);
//                std::cout << "extracted non-sparse sequences: "<<nonSparseSequences.size()  << "! Removing sparse durations" <<std::endl;
//                sequences = extractMonthlySequences(nonSparseSequences, durationSparsity,
//                                                    durationSparsityValue,startPositions.size(),numOfThreads);
//
//            }else {
//                std::cout << "extracting (non-sparse) sequences" << std::endl;
//                sequences = extractNonSparseSequences(dbMart, startPositions, sequenceCount,
//                                                      numOfThreads, DURATION_IN_MONTHS,14);
//            }
//
//
//        }else {
//            std::cout << "creating sequences with temporal buckets" << std::endl;
//
//            sequences = extractTemporalBuckets(dbMart, startPositions, sequenceCount, numOfThreads,
//                                               durationPeriods, daysForCoOccurrence, sparsity_value,
//                                               removeSparseTemporalBuckets);
//        }
//        return sequences;
//
//    }

    std::vector<temporalSequence> sequenceWorkflowFromCsVFiles(const std::vector<std::string>& inputFilePaths, char inputFileDelimiter,
                                                               int patIDColumns[], int phenxColumns[], int dateColumns[], bool storeSeqDuringCreation, const std::string& outPutDirectory,
                                                               const std::string& outputFilePrefix, bool removeSparseSequences, double sparsity_value,
                                                               bool createTemporalBuckets, double durationPeriods,
                                                               unsigned int daysForCoOccurrence, bool durationSparsity,
                                                               double durationSparsityValue, bool removeSparseTemporalBuckets, unsigned int patIdLength, unsigned int numOfThreads, unsigned int phenxIdLength){
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
                                removeSparseTemporalBuckets, patIdLength, numOfThreads, phenxIdLength);
    }
}//tspm
