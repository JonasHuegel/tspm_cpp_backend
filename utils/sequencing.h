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

namespace tspm {
    const int daysPerWeek = 7;
    const double daysPerMonth = 30.437;
    const double DURATION_IN_DAYS = 1;
    const double DURATION_IN_WEEKS = 7;
    const double DURATION_IN_MONTHS = 30.472;
    const double DURATION_IN_YEARS = 365.25;

    std::vector<temporalSequence>
    applyDurationSparsity(std::vector<temporalSequence> &sequences, bool storeDurationInSequence,
                            double sparsity, size_t numOfPatients, int numOfThreads = 1,
                            unsigned int bitShift = 52);

    unsigned int getDurationPeriod(unsigned int duration, double durationPeriods, unsigned int daysForCoOccurrence);


    std::vector<std::vector<temporalSequence>>
    splitSequenceVectorInChunks(std::vector<temporalSequence> &sequences, unsigned int chunks,
                                double durationPeriods = DURATION_IN_DAYS,
                                unsigned int daysForCoOccurrence = 0);

    std::vector<size_t> extractStartPositions(std::vector<dbMartEntry> &dbMart);

    size_t createSequencesFromFiles(std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                    const std::string &outPutDirectory, const std::string &outputFilePrefix,
                                    int patIDColumns[], int phenxColumns[], int dateColumns[], size_t numOfPatients,
                                    int patIdLength, int numOfThreads, unsigned int phenxIdLength);

    std::vector<temporalSequence> extractDynamicTemporalBuckets(std::vector<temporalSequence> &nonSparseSequences, size_t numOfPatients,
                                                         int numOfThreads, double sparsity, bool removeSparseBuckets);

    std::vector<temporalSequence> extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,double sparsity,
                                                            unsigned int &numOfThreads, double durationPeriod = DURATION_IN_MONTHS, int daysForCoOccurrence = 14, unsigned int phenxIdLength = 7);

    std::vector<temporalSequence>
    extractSequencesOfInterest(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              std::map<std::int64_t, size_t> &sequenceOfInterest, int numOfThreads,
                              double durationPeriod = DURATION_IN_MONTHS,
                              int daysForCoOccurrence = 14, unsigned int phenxIdLength = 7);

    std::vector<temporalSequence>
    extractSparseSequences(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              int numOfThreads, double durationPeriod = DURATION_IN_MONTHS,
                              int daysForCoOccurrence = 14, unsigned int phenxIdLength = 7);

    size_t
    writeSequencesFromArrayToFile(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              const std::string &outPutDirectory, const std::string &outputFilePrefix,
                              int patIDLength = 7, int numOfThreads = 1 , unsigned int phenxIdLength = 7);


    unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration);

}//tspm
#endif //TSPM_CPP_BACKEND_SEQUENCING_H
