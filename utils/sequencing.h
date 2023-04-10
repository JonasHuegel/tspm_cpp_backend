//
// Created by jonas on 03.03.23.
//

#ifndef TSPM_CPP_BACKEND_SEQUENCING_H
#define TSPM_CPP_BACKEND_SEQUENCING_H
#include "utils.h"
#include <algorithm>
#include <set>
#include <mutex>
#include <omp.h>
#include <math.h>


struct temporalSequence{
    long seqID;
    unsigned int duration;
    unsigned int patientID;
}__attribute__((packed));

int extractSequencesFromArray(dbMartEntry * dbMart, size_t numOfPatients, size_t * startPositions,
                              size_t numberOfDbMartEntries,  const std::string& outPutDirectory,
                              const std::string& outputFilePrefix, int patIDLength = 7);

unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration);


size_t createSequencesFromFiles (std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                 const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                 int patIDColumns[], int phenxColumns[], int dateColumns[], long maxPatID);

std::vector<temporalSequence> createSequencesWithDuration(std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                                               const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                                               int patIDColumns[], int phenxColumns[], int dateColumns[], long maxPatID, size_t sequenceCount,
                                                               const std::map<long,size_t>& sequenceMap);

std::vector<temporalSequence>
extractTemporalSequences(const std::vector<std::string> &inputFilePaths, char inputFileDelimiter,
                         const std::string &outPutDirectory, const std::string &outputFilePrefix, int *patIdColumns,
                         int *phenxColumn, int *dateColumns, long maxPatID, size_t sequenceCount,
                         const std::map<long, size_t>& sequences, size_t sparsityThreshold, bool removeSparseBuckets);

int sequenceWorkflow(bool temporal, bool removeSparseBuckets, const std::vector<std::string>& inputFilePaths, char inputFileDelimiter,
                     const std::string& outPutDirectory, const std::string& outputFilePrefix,
                     int patIDColumns[], int phenxColumns[], int dateColumns[], long maxPatID, double sparsity_value);


long writeSequencesAsCsV(std::string fileName, std::string filepath, char delimiter, size_t numOfSequences, temporalSequence * temporalSequences, bool debug =false);
#endif //TSPM_CPP_BACKEND_SEQUENCING_H
