/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* workflows.cpp
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
                std::cout << "removing sparse sequence" << std::endl;
                allSequences = tspm::removeSparseSequences(allSequences,startPositions.size(), sparsity_value, numOfThreads);
                std::cout << "extracted non-sparse sequences: "<<allSequences.size()  << "!" << std::endl;

            }
            sequences = allSequences;
        }else if(removeSparseSequences){
            std::cout << "extracting transitive sequences" << std::endl;
            sequences = extractNonSparseSequences(dbMart,startPositions,sparsity_value,numOfThreads,
                                                      durationPeriods, daysForCoOccurrence,phenxIdLength);
             std::cout << "extracted non-sparse sequences: "<<sequences.size()  << "!" << std::endl;
        }else{
            sequences = extractSparseSequences(dbMart,startPositions, numOfThreads, durationPeriods, daysForCoOccurrence);
            std::cout << "extracted sparse sequences: "<<sequences.size()  << "!" << std::endl;
        }
        if(!createTemporalBuckets){
            if(durationSparsity){
                std::cout << "Removing sparse durations" <<std::endl;
                bool storeDurationInSequence = false;
                sequences = applyDurationSparsity(sequences, storeDurationInSequence,
                                                    durationSparsityValue,startPositions.size(),numOfThreads);
            }
        }else {
            std::cout << "creating sequences with temporal buckets" << std::endl;
            sequences = extractDynamicTemporalBuckets(sequences, startPositions.size(), numOfThreads,
                                               sparsity_value, removeSparseTemporalBuckets);
        }
        return sequences;
    }


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
