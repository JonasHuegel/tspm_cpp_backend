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

const int daysPerWeek = 7;
const double daysPerMonth = 30.437;
const double DURATION_IN_DAYS = 1;
const double DURATION_IN_WEEKS = 7;
const double DURATION_IN_MONTHS = 30.472;
const double DURATION_IN_YEARS = 365.25;

std::vector<temporalSequence> extractMonthlySequences(std::vector<temporalSequence> &sequences, bool durationSparsity,
                                                      double sparsity, size_t numOfPatients, int numOfThreads = 1,
                                                      double durationPeriods = DURATION_IN_MONTHS,
                                                      unsigned int daysForCoOoccurence = 14,  unsigned int bitShift = 52);

unsigned int getDurationPeriod(unsigned int duration, double durationPeriods, unsigned int daysForCoOoccurence);


std::vector<std::vector<temporalSequence>> splitSequenceVectorInChunkes(std::vector<temporalSequence> &sequences, unsigned int chunks, double durationPeriods = DURATION_IN_DAYS,
                                                                        unsigned int daysForCoOoccurence = 0);
std::vector<size_t> extractStartPositions(std::vector<dbMartEntry> &dbMart);

size_t createSequencesFromFiles (std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                 const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                 int patIDColumns[], int phenxColumns[], int dateColumns[], size_t numOfPatients,
                                 int  patIdLength, int numOfThreads);

std::vector<temporalSequence> extractTemporalBuckets(std::vector<dbMartEntry> &dbMart, size_t numOfPatients,
                                                     const size_t *startPositions,
                                                     std::map<std::int64_t, size_t> &nonSparseSequencesIDs, int numOfThreads,
                                                     double durationPeriods, unsigned int daysForCoOoccurence,
                                                     size_t sparsityThreshold, bool removeSparseBuckets);

std::vector<temporalSequence>
extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, size_t numOfPatients, const size_t *startPositions,
                          std::map<std::int64_t, size_t> &nonSparseSequencesIDs, int numOfThreads, double durationPeriod = DURATION_IN_MONTHS,
                          int daysForCoOccurence = 14);

size_t extractSequencesFromArray(std::vector<dbMartEntry> &dbMart, size_t numOfPatients, const size_t * startPositions,
                                 const std::string& outPutDirectory,const std::string& outputFilePrefix,
                                 int patIDLength = 7, int numOfThreads = 1);


unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration);

#endif //TSPM_CPP_BACKEND_SEQUENCING_H
