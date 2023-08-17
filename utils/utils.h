/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* utils.h
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
    /**
     * Function to store the sequences in a non-binary format in a text file. The durations should be already shifted into the sequence.
     * */    
    size_t writeSequencesToFile(std::string patientFilename, std::vector<std::int64_t> &sequences);

    /**
     * Function to read sequences from a binary file. Provides the option to shift the duration out of the sequence.
    */
    std::vector<temporalSequence>
            readSequencesFromFiles(const std::string &outputDir, const std::string &file_prefix, int numberOfPatients,
                                   bool storesDuration, unsigned int patIdLength, unsigned int bitShift = 24,
                                   int numOfThreads = 1);
    /**
     * Function to remove all sequences that are not present in a given percentage of the patients.
    */
    std::vector<temporalSequence>
            removeSparseSequences(std::vector<temporalSequence> &sequences, size_t numOfPatients, double sparsity,
                               unsigned int numOfThreads);
    
    /**
     * Function to summarize the sequences by counting how many sequences are in each bucket foreach sequence id. Returns the result as vector.
    */
    std::vector<std::pair<temporalSequence, size_t>>
    summarizeSequencesAsVector(std::vector<temporalSequence> &sequences, bool includeDuration,
                               std::vector<unsigned int> durationBuckets, bool patientLevel = false, unsigned int numOfThreads = 1);
    /**
     * Function to summarize the sequences by counting how many sequences are in each bucket foreach sequence id. Returns the result as map.
    */
    std::map<std::int64_t, size_t>
            summarizeSequencesAsMap(std::vector<temporalSequence> &sequences, unsigned int &numOfThreads);

     /**
      * Funtion to split a line in multiple strings, based in a given delimiter, used to extract the values from the dbMart text files.
     */
    std::vector<std::string> getTokensFromLine(const std::string &line, char delim);

     /**
      * Function to merge two phenotypes into a seqeunce.
     */
    std::int64_t createSequence(int phenotypeA, int phenotypeB, unsigned int phenotypelenght = 7);

     /**
      * Function to write the sequences in a binary file. The duration should already be shifted into the file
     */
    size_t writeSequencesToBinaryFile(const std::string& patientFilename, std::vector<std::int64_t> &sequences);

     /**
      * Function to read in the binary files and summarize the files by counting the number of sequences per sequence ID. Returns the results as a map.
     */
    std::map<std::int64_t, size_t>
    summarizeSequencesFromFiles(const std::string &outputDir, const std::string &file_prefix, int numberOfPatients,
                                bool storesDuration, unsigned int patIdLength, unsigned int bitShift = 24,
                                int numOfThreads = 1);

     /**
      * Function that returns the size of the given file. Works on MAC, Linux and Windows, requires C++17.
     */
    std::int64_t getFileSize(const std::string &filename);

    /**
     * Function to calculate the duration of sequence. 
    */
    unsigned int getDuration(std::int64_t startDate, std::int64_t endDate);
    
     /**
      * Function extract the timestamp (date) from the string and transform into int 64 representation
     */
    std::int64_t getTimeFromString(const char *date_string);

     /**
      * Function to read a dbMart from a csc file and store it as an vector of dbMartEntries. 
     */
    std::vector<dbMartEntry>
            extractDBMartFromCsv(FILE *csv_file, int patIdColumn, int phenotypeIDColumn,
                                 int dateColumn, char delim = ',');

    /**
     * Function to store a m array of sequences as a csv file. Can include the duration.
     *
     */                             
    std::int64_t
    writeSequencesAsCsV(std::string fileName, std::string filepath, char delimiter, size_t numOfSequences,
                        temporalSequence *temporalSequences, bool debug = false);

    /**
     * Function to create the path to the outputh direcotory, creates all required directories on the path.
    */
    std::filesystem::path createOutputFilePath(const std::string &outPutDirectory);
    
      /**
       * Function to extract tje startPositons of each sequence ID. Requires the sequences to be sorted after their ID. 
      */
    std::vector<size_t> getSequenceStartPositions(std::vector<temporalSequence> &sequences);

    /**
     * Extracts and returns the startPhenx of a sequence. Using the struct as input.
     */     
    unsigned int getStartPhenx(temporalSequence &sequence, unsigned int lengthOfPhenx);

     /**
     * Extracts and returns the startPhenx of a sequence. Using the numeric sequence representation as input. 
     */
    unsigned int getStartPhenx(int64_t sequence, unsigned int lengthOfPhenx);

     /**
      * Function to check if a phenx is part of given list of phenx.
     */
    bool isPhenxOfInterest(unsigned int phenx, std::vector<unsigned int> phenxsOfInterest);
     /**
      * Function to extract the end from the sequence struct.
     */
    unsigned int getEndPhenx(temporalSequence &sequence, unsigned int lengthOfPhenx);

     /**
      * Function to extract the end from the numeric sequence representation.
     */
    unsigned int getEndPhenx(int64_t sequence, unsigned int lengthOfPhenx);

//
     /**
      * Function to determine the candidate bucket based on the given duration and bucket thresholds.
      * Requires that the duration is not stored in the sequence id, but instead in the duration field of the struct
     */
    unsigned int getCandidateBucket(unsigned int duration, std::vector<unsigned int> lowerBucketLimits);


     /**
      * Function to extract the set of all phenx that are the end phenx of sequences that start with a given phen,x
     */
    std::set<unsigned int>
    extractEndPhenxWithGivenStartPhenx(std::vector<temporalSequence> &originalSequences, std::uint64_t minDuration,
                                       unsigned int bitShift, unsigned int lengthOfPhenx,
                                       std::vector<unsigned int> &phenxOfInterest, int &numOfThreads);
   
     /**
      * Function to extract all sequences that starts with one phenx of given list of phenx. 
     */
    std::vector<temporalSequence>
    extractSequencesWithSpecificStart(std::vector<temporalSequence> &originalSequences, std::uint64_t minDuration,
                                      unsigned int bitShift, unsigned int lengthOfPhenx,
                                      std::vector<unsigned int> &phenxOfInterest, int &numOfThreads);
     /**
      * Function to extract all sequences that end with one phenx of given list of phenx.
     */
    std::vector<temporalSequence>
            extractSequencesWithEnd(std::vector<temporalSequence> &originalSequences, unsigned int bitShift,
                                    unsigned int lengthOfPhenx, std::set<unsigned int> &allEndPhenx, int &numOfThreads);
}//tspm
#endif //TSPM_UTILS_H
