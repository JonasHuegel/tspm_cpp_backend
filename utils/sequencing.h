/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* sequencing.h
***********************************************************************************
* MIT License
*
* Copyright (c) 2023 Jonas Hügel
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
**********************************************************************************/

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
    /**
     * constants required to calculate the durations of a sequence.
    */
    const int daysPerWeek = 7;
    const double daysPerMonth = 30.437;
    const double DURATION_IN_DAYS = 1;
    const double DURATION_IN_WEEKS = 7;
    const double DURATION_IN_MONTHS = 30.472;
    const double DURATION_IN_YEARS = 365.25;

    /**
     * Function to remove sparse sequences based on  their duration buckets.
    */
    std::vector<temporalSequence>
    applyDurationSparsity(std::vector<temporalSequence> &sequences, bool storeDurationInSequence,
                            double sparsity, size_t numOfPatients, int numOfThreads = 1,
                            unsigned int bitShift = 52);
    /**
     * Function to get the bucket that represents the duration periond, e.g. 0 (co-occurence), <1, 1<x<3, >3 months.
    */
    unsigned int getDurationPeriod(unsigned int duration, double durationPeriods, unsigned int daysForCoOccurrence);

    /**
     * utils function from the sequencing process, allows to split one large vector of mined sequences in multiple chunks.
     * to allow to process them in parallel
    */
    std::vector<std::vector<temporalSequence>>
    splitSequenceVectorInChunks(std::vector<temporalSequence> &sequences, unsigned int chunks,
                                double durationPeriods = DURATION_IN_DAYS,
                                unsigned int daysForCoOccurrence = 0);
    /**
     * requires a sorted sequence vector and extracts the start postions of each sequence block.
    */
    std::vector<size_t> extractStartPositions(std::vector<dbMartEntry> &dbMart);

    size_t createSequencesFromFiles(std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                    const std::string &outPutDirectory, const std::string &outputFilePrefix,
                                    int patIDColumns[], int phenxColumns[], int dateColumns[], size_t numOfPatients,
                                    int patIdLength, int numOfThreads, unsigned int phenxIdLength);
    /**
     * Function that groups all sequences in a number of buckets depending on their duration. The number of buckets is computed for log(min(duration)-max(duration)) for each sequence id.
    */
    std::vector<temporalSequence> extractDynamicTemporalBuckets(std::vector<temporalSequence> &nonSparseSequences, size_t numOfPatients,
                                                         int numOfThreads, double sparsity, bool removeSparseBuckets);
    /**
     * Extracts all transitive sequences from the dbMart and applies the sparsity screening afterwards to remove sparse sequences.
    */
    std::vector<temporalSequence> extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,double sparsity,
                                                            unsigned int &numOfThreads, double durationPeriod = DURATION_IN_MONTHS, int daysForCoOccurrence = 14, unsigned int phenxIdLength = 7);
    /**
     * Function that extracts only a given subset of sequences from the dbMart.
    */
    std::vector<temporalSequence>
    extractSequencesOfInterest(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              std::map<std::int64_t, size_t> &sequenceOfInterest, int numOfThreads,
                              double durationPeriod = DURATION_IN_MONTHS,
                              int daysForCoOccurrence = 14, unsigned int phenxIdLength = 7);
    /**
     * Function to extract all sequencues without the sparsity screening.
    */
    std::vector<temporalSequence>
    extractSparseSequences(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              int numOfThreads, double durationPeriod = DURATION_IN_MONTHS,
                              int daysForCoOccurrence = 14, unsigned int phenxIdLength = 7);
    /**
     * Function that extracts all sequences from the dbmart and stores them in files. The function creates one file per patient.
    */
    size_t
    writeSequencesFromArrayToFile(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              const std::string &outPutDirectory, const std::string &outputFilePrefix,
                              int patIDLength = 7, int numOfThreads = 1 , unsigned int phenxIdLength = 7);

    /**
     * Function to  caluclate the duration bucket base on the given range(threshold) for the sequence id.
    */
    unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration);

}//tspm
#endif //TSPM_CPP_BACKEND_SEQUENCING_H
