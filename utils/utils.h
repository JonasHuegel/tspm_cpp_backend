
#ifndef TSPM_UTILS_H
#define TSPM_UTILS_H
#include "../lib/ips4o/ips4o.hpp"
#include "dbMartEntry.h"
#include "temporalSequence.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstring>
#include <map>
#include <omp.h>
#include <mutex>
#include <set>
#include <filesystem>
#include "sequencing.h"

namespace tspm {
    size_t writeSequencesToFile(std::string patientFilename, std::vector<std::int64_t> &sequences);

    std::vector<temporalSequence>
            readSequencesFromFiles(const std::string &outputDir, const std::string &file_prefix, int numberOfPatients,
                                   bool storesDuration, unsigned int patIdLength, unsigned int bitShift = 24,
                                   int numOfThreads = 1);

    std::vector<temporalSequence>
            removeSparseSequences(std::vector<temporalSequence> &sequences, size_t numOfPatients, double sparsity,
                               unsigned int numOfThreads);

    std::vector<std::pair<temporalSequence, size_t>>
    summarizeSequencesAsVector(std::vector<temporalSequence> &sequences, bool includeDuration,
                               std::vector<size_t> durationBuckets, unsigned int &numOfThreads);

    std::map<std::int64_t, size_t>
            summarizeSequencesAsMap(std::vector<temporalSequence> &sequences, unsigned int &numOfThreads);

    std::vector<std::string> getTokensFromLine(const std::string &line, char delim);

    std::int64_t createSequence(int phenotypeA, int phenotypeB, unsigned int phenotypelenght = 7);

    size_t writeSequencesToBinaryFile(const std::string& patientFilename, std::vector<std::int64_t> &sequences);

    std::int64_t addDurationToSequence(std::int64_t &back, std::int64_t startDate, std::int64_t endDate);

    std::map<std::int64_t, size_t>
    summarizeSequencesFromFiles(const std::string &outputDir, const std::string &file_prefix, int numberOfPatients,
                                bool storesDuration, unsigned int patIdLength, unsigned int bitShift = 24,
                                int numOfThreads = 1);

    std::int64_t getFileSize(const std::string &filename);

    unsigned int getDuration(std::int64_t startDate, std::int64_t endDate);

    std::int64_t getTimeFromString(const char *date_string);

    std::vector<dbMartEntry>
            extractDBMartFromCsv(FILE *csv_file, int patIdColumn, int phenotypeIDColumn,
                                 int dateColumn, char delim = ',');

    std::int64_t
    writeSequencesAsCsV(std::string fileName, std::string filepath, char delimiter, size_t numOfSequences,
                        temporalSequence *temporalSequences, bool debug = false);

    std::filesystem::path createOutputFilePath(const std::string &outPutDirectory);

    std::vector<size_t> getSequenceStartPositions(std::vector<temporalSequence> &sequences);

    unsigned int getStartPhenx(temporalSequence &sequence, unsigned int lengthOfPhenx);

    unsigned int getStartPhenx(int64_t sequence, unsigned int lengthOfPhenx);

    bool isPhenxOfInterest(unsigned int phenx, std::vector<unsigned int> phenxsOfInterest);

    unsigned int getEndPhenx(temporalSequence &sequence, unsigned int lengthOfPhenx);

    unsigned int getEndPhenx(int64_t sequence, unsigned int lengthOfPhenx);

//we assume that the duration is not stored in the sequence id, but instead in the duration field of the struct
    unsigned int getCandidateBucket(unsigned int duration, std::vector<unsigned int> lowerBucketLimits);

    std::set<unsigned int>
    extractEndPhenxWithGivenStartPhenx(std::vector<temporalSequence> &originalSequences, std::uint64_t minDuration,
                                       unsigned int bitShift, unsigned int lengthOfPhenx,
                                       std::vector<unsigned int> &phenxOfInterest, int &numOfThreads);

    std::vector<temporalSequence>
    extractSequencesWithSpecificStart(std::vector<temporalSequence> &originalSequences, std::uint64_t minDuration,
                                      unsigned int bitShift, unsigned int lengthOfPhenx,
                                      std::vector<unsigned int> &phenxOfInterest, int &numOfThreads);

    std::vector<temporalSequence>
            extractSequencesWithEnd(std::vector<temporalSequence> &originalSequences, unsigned int bitShift,
                                    unsigned int lengthOfPhenx, std::set<unsigned int> &allEndPhenx, int &numOfThreads);
}//tspm
#endif //TSPM_UTILS_H
