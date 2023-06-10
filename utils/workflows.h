//
// Created by jonas on 18.04.23.
//
#ifndef TSPM_CPP_BACKEND_WORKFLOWS_H
#define TSPM_CPP_BACKEND_WORKFLOWS_H
#include "sequencing.h"
#include "sorter.h"
#include <ctime>

std::vector<temporalSequence> sequenceWorkflow(std::vector<dbMartEntry> &dbMart, bool storeSeqDuringCreation, const std::string &outPutDirectory,
                                               const std::string &outputFilePrefix, bool removeSparseSequences,
                                               double sparsity_value, bool createTemporalBuckets, double durationPeriods,
                                               unsigned int daysForCoOccurrence, bool durationSparsity,
                                               double durationSparsityValue, bool removeSparseTemporalBuckets,
                                               unsigned int patIdLength, unsigned int numOfThreads);



std::vector<temporalSequence> sequenceWorkflowFromCsVFiles(const std::vector<std::string>& inputFilePaths, char inputFileDelimiter,
                                                           int patIDColumns[], int phenxColumns[], int dateColumns[], bool storeSeqDuringCreation, const std::string& outPutDirectory,
                                                           const std::string& outputFilePrefix, bool removeSparseSequences, double sparsity_value,
                                                           bool createTemporalBuckets, double durationPeriods,
                                                           unsigned int daysForCoOccurrence, bool durationSparsity,
                                                           double durationSparsityValue, bool removeSparseTemporalBuckets, unsigned int patIdLength, unsigned int numOfThreads);
#endif //TSPM_CPP_BACKEND_WORKFLOWS_H
