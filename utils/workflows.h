/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* workflows.h
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
#ifndef TSPM_CPP_BACKEND_WORKFLOWS_H
#define TSPM_CPP_BACKEND_WORKFLOWS_H
#include "sequencing.h"
#include "sorter.h"
#include <ctime>
namespace tspm {

    /**
     * Main workflow function which orchestarte the sequencing from dbMart that is already stored as a vector of dbMartEntries. 
    */
    std::vector<temporalSequence>
    sequenceWorkflow(std::vector<dbMartEntry> &dbMart, bool storeSeqDuringCreation, const std::string &outPutDirectory,
                     const std::string &outputFilePrefix, bool removeSparseSequences,
                     double sparsity_value, bool createTemporalBuckets, double durationPeriods,
                     unsigned int daysForCoOccurrence, bool durationSparsity,
                     double durationSparsityValue, bool removeSparseTemporalBuckets,
                     unsigned int patIdLength, unsigned int numOfThreads, unsigned int phenxIdLength = 7);

    /**
     * Main workflow function which orchestarte the sequencing from dbMart that is stored as one or multiple csv files.
    */
    std::vector<temporalSequence>
    sequenceWorkflowFromCsVFiles(const std::vector<std::string> &inputFilePaths, char inputFileDelimiter,
                                 int patIDColumns[], int phenxColumns[], int dateColumns[], bool storeSeqDuringCreation,
                                 const std::string &outPutDirectory,
                                 const std::string &outputFilePrefix, bool removeSparseSequences, double sparsity_value,
                                 bool createTemporalBuckets, double durationPeriods,
                                 unsigned int daysForCoOccurrence, bool durationSparsity,
                                 double durationSparsityValue, bool removeSparseTemporalBuckets,
                                 unsigned int patIdLength, unsigned int numOfThreads, unsigned int phenxIdLength =7);
}//tspm
#endif //TSPM_CPP_BACKEND_WORKFLOWS_H
