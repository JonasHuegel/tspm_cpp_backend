//
// Created by jonas on 03.03.23.
//

#ifndef TSPM_CPP_BACKEND_SEQUENCING_H
#define TSPM_CPP_BACKEND_SEQUENCING_H
#include "sorter.h"
#include "utils.h"
#include <algorithm>
#include <set>
#include <mutex>
#include <omp.h>
#include <math.h>

std::vector<size_t> extractStartPositions(std::vector<dbMartEntry> &dbMart);

size_t createSequencesFromFiles (std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                 const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                 int patIDColumns[], int phenxColumns[], int dateColumns[], size_t numOfPatients,
                                 int  patIdLength, int numOfThreads);

std::vector<temporalSequence> extractTemporalSequences(std::vector<dbMartEntry> &dbMart, size_t numOfPatients,
                                                       const size_t *startPositions,
                                                       std::map<long, size_t> &nonSparseSequencesIDs, int numOfThreads,
                                                       bool durationInWeeks, bool durationInMonths,
                                                       size_t sparsityThreshold, bool removeSparseBuckets);

std::vector<temporalSequence> extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, size_t numOfPatients,
                                                        const size_t *startPositions, std::map<long,
                                                        size_t> &nonSparseSequencesIDs, int numOfThreads,
                                                        bool durationInWeeks = false, bool durationInMonths = false);

size_t extractSequencesFromArray(std::vector<dbMartEntry> &dbMart, size_t numOfPatients, const size_t * startPositions,
                                 const std::string& outPutDirectory,const std::string& outputFilePrefix,
                                 int patIDLength = 7, int numOfThreads = 1);


unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration);

#endif //TSPM_CPP_BACKEND_SEQUENCING_H
