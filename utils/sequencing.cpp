/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* sequencing.cpp
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

#include "sequencing.h"
#include <string>

namespace tspm {

    std::vector<temporalSequence> extractDynamicTemporalBuckets( std::vector<temporalSequence> &nonSparseSequences, size_t numOfPatients,
                                                          int numOfThreads, double sparsity, bool removeSparseBuckets)  {
        ips4o::parallel::sort(nonSparseSequences.begin(),nonSparseSequences.end(), timedSequencesSorter, numOfThreads);
        std::vector<size_t> startPositions = getSequenceStartPositions(nonSparseSequences);
        std::cout << "computing dynamic temporal buckets" << std::endl;
        size_t sparsityThreshold = sparsity * numOfPatients;
#pragma omp parallel for shared(nonSparseSequences, removeSparseBuckets, sparsityThreshold,std::cout, startPositions) default (none)
        for (size_t i = 0; i < startPositions.size(); ++i) {
            size_t startPos = startPositions[i];
            size_t endPos;
            if (i < startPositions.size() - 1) {
                endPos = startPositions[i + 1];
            } else {
                endPos = nonSparseSequences.size();
            }
            auto startIterator = nonSparseSequences.begin() + startPos;
            auto endIterator = nonSparseSequences.begin() + endPos;

            //       compute buckets
            unsigned int min = UINT32_MAX, max = 0;

            for (auto it = startIterator; it != endIterator; ++it) {
                min = std::min(min, it->duration);
                max = std::max(max, it->duration);
            }
            unsigned int range = max - min;
            int numberOfBuckets = std::max(4, int(std::log(range)));

//      add bucket id to sequence id
            for (auto it = startIterator; it != endIterator; ++it) {
                unsigned int bucket = getBucket(min, max, range / numberOfBuckets, it->duration);
                int bitShift = sizeof(std::int64_t) * 8 - 8;
                it->seqID = bucket << bitShift && it->seqID;
            }
        }
        if(removeSparseBuckets){
            nonSparseSequences = removeSparseSequences(nonSparseSequences,numOfPatients, sparsity, numOfThreads);
        }
        return nonSparseSequences;
    }

    std::vector<size_t> extractStartPositions(std::vector<dbMartEntry> &dbMart) {
        ips4o::parallel::sort(dbMart.begin(), dbMart.end(), dbMartSorter);
        std::vector<size_t> startPositions;
        startPositions.emplace_back(0);
        for (int i = 1; i < dbMart.size(); ++i) {
            if (dbMart[i].patID != dbMart[i - 1].patID) {
                startPositions.emplace_back(i);
            }
        }
        return startPositions;
    }


    size_t
    writeSequencesFromArrayToFile(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                              const std::string &outPutDirectory, const std::string &outputFilePrefix, int patIDLength,
                              int numOfThreads, unsigned int phenxIdLength) {
        omp_set_num_threads(numOfThreads);
        size_t numOfPatients = startPositions.size();
        auto numOfSequences = new size_t[numOfThreads];
        std::fill_n(numOfSequences, numOfThreads, 0);
        size_t numberOfDbMartEntries = dbMart.size();
        omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, patIDLength, outPutDirectory, outputFilePrefix, numOfSequences, phenxIdLength)
        for (size_t i = 0; i < numOfPatients; ++i) {
            size_t startPos = startPositions[i];
            size_t endPos = i < numOfPatients - 1 ? startPositions[i + 1] : numberOfDbMartEntries;
            size_t numOfPatientEntries = endPos - startPos;

            size_t numberOfSequences = (numOfPatientEntries * (numOfPatientEntries - 1)) / 2;
            std::vector<std::int64_t> sequences;
            sequences.reserve(numberOfSequences);
            for (size_t j = startPos; j < endPos - 1; ++j) {
                for (size_t k = j + 1; k < endPos; ++k) {
                    sequences.emplace_back(createSequence(dbMart[j].phenID, dbMart[k].phenID, phenxIdLength));
                }
            }
            numOfSequences[omp_get_thread_num()] += sequences.size();
            std::string patIDString = std::to_string(i);
            patIDString.insert(patIDString.begin(), patIDLength - patIDString.size(), '0');
            std::string patientFileName = std::string(outPutDirectory).append(outputFilePrefix).append(patIDString);
            writeSequencesToBinaryFile(patientFileName, sequences);
//        writeSequencesToFile(patientFileName, sequences);
        }
        size_t sumOfSequences = 0;
        for (int i = 0; i < numOfThreads; ++i) {
            sumOfSequences += numOfSequences[i];
        }
        return sumOfSequences;
    }


    std::vector<size_t> getStartPositionsFromSequenceVector(std::vector<temporalSequence> &sequence) {
        std::vector<size_t> startPositions;
        unsigned int lastPatId = sequence[0].patientID;
        for (int i = 1; i < sequence.size(); ++i) {
            if (lastPatId != sequence[i].patientID) {
                startPositions.emplace_back(i);
                lastPatId = sequence[i].patientID;
            }
        }
        return startPositions;
    }


    std::vector<std::vector<temporalSequence>> splitSequenceVectorInChunks(std::vector<temporalSequence> &sequences,
                                                                           unsigned int chunks, double durationPeriods,
                                                                           unsigned int daysForCoOccurrence) {
        ips4o::parallel::sort(sequences.begin(), sequences.end(), timedSequencesSorter);
        std::vector<std::vector<temporalSequence>> localSequences;
        //    Split Sequences in sub vectors for parallel access
        auto endPos = sequences.begin();
        size_t numOfSequencesPerChunk = sequences.size() / chunks;
        for (size_t i = 0; i < chunks; ++i) {
            if (sequences.empty()) {
                localSequences.emplace_back();
                continue;
            }
            if (numOfSequencesPerChunk >= sequences.size()) {
                endPos = sequences.end();
            } else {
                endPos = sequences.begin() + numOfSequencesPerChunk;
                std::int64_t lastSeq = endPos->seqID;
                std::int64_t lastDur = endPos->duration;
                for (; endPos != sequences.end() && lastSeq == endPos->seqID &&
                       endPos->duration == lastDur; ++endPos);
            }
            localSequences.emplace_back(std::vector<temporalSequence>(sequences.begin(), endPos));
            sequences.erase(sequences.begin(), endPos);
            sequences.shrink_to_fit();
        }
        return localSequences;
    }


    std::vector<temporalSequence>
    applyDurationSparsity(std::vector<temporalSequence> &sequences, bool storeDurationInSequence,
                            double sparsity, size_t numOfPatients, int numOfThreads,
                            unsigned int bitShift) {
        ips4o::parallel::sort(sequences.begin(), sequences.end(), timedSequencesSorter, numOfThreads);
        size_t sparsityThreshold = numOfPatients * sparsity;
        std::vector<size_t> startPositions = getSequenceStartPositions(sequences);
#pragma omp parallel for default (none) shared(startPositions, bitShift, sparsityThreshold, storeDurationInSequence, sequences)
        for (size_t i = 0; i < startPositions.size(); ++i) {
            size_t startPos = startPositions[i];
            size_t endPos;
            if (i < startPositions.size() - 1) {
                endPos = startPositions[i + 1];
            } else {
                endPos = sequences.size();
            }
            auto sparsityIt = sequences.begin() + startPos;
            std::set<unsigned int> sequenceInPatient;
            std::uint64_t lastSequence;
            if(storeDurationInSequence) {
                lastSequence =
                        sparsityIt->duration << bitShift | sparsityIt->seqID;
                sparsityIt->seqID =
                        sparsityIt->duration << bitShift | sparsityIt->seqID;
            }else {
                lastSequence = sparsityIt->duration;
            }
            ++sparsityIt;
            size_t count = 0;
            // since the sequences are order by ID as first and duration as second iterator, we can iterate over them and
            // integrate the duration in month and directly remove check for sparsity if the updated sequenceID change
            auto endIt = sequences.begin() + endPos;
            while (sparsityIt != endIt) {
                if(storeDurationInSequence) {
                    sparsityIt->seqID = sparsityIt->duration << bitShift | sparsityIt->seqID;
                }
                if (sparsityIt->seqID == lastSequence) {
                    sequenceInPatient.insert(sparsityIt->patientID);
                    ++count;
                    ++sparsityIt;
                } else {
                    lastSequence = sparsityIt->seqID;
                    if (sequenceInPatient.size() < sparsityThreshold) {
                        for (auto it = sparsityIt - (count + 1); it != sparsityIt; ++it) {
                            //use UINT32_MAX later after sorting to identify and remove sparse sequences!
                            it->patientID = UINT32_MAX;
                        }
                    }
                    ++sparsityIt;
                    sequenceInPatient.clear();
                    count = 0;
                }
            }
        }

        ips4o::parallel::sort(sequences.begin(), sequences.end(), timedSequenceByPatientIDSorter, numOfThreads);
        auto it = sequences.begin();
        while(it!= sequences.end() && it->patientID < UINT32_MAX ){
            ++it;
        }
        sequences.erase(it, sequences.end());
        sequences.shrink_to_fit();
        return sequences;
    }




    std::vector<temporalSequence> extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions, double sparsity,
                                                            unsigned int &numOfThreads, double durationPeriod, int daysForCoOccurrence, unsigned int phenxIdLength){

        std::vector<temporalSequence> sequences = extractSparseSequences(dbMart, startPositions, numOfThreads,
                                                                            durationPeriod, daysForCoOccurrence);
        std::cout << "Number of overall extracted Sequences: "<< sequences.size() << std::endl;
        std::cout << "Removing sparse sequences:" << std::endl;
        sequences = removeSparseSequences (sequences, startPositions.size(), sparsity, numOfThreads);
        return sequences;
    }



    unsigned int getDurationPeriod(unsigned int duration, double durationPeriods, unsigned int daysForCoOccurrence) {
        //if the duration is less than 2 weeks it is considered as co-occurrence and therefore the distance is 0
        if (daysForCoOccurrence <= 0)
            daysForCoOccurrence = 1;
        if (duration < daysForCoOccurrence) {
            return 0;
        }
        if (durationPeriods == DURATION_IN_DAYS)
            return duration;
        else
            return std::ceil(duration / durationPeriods);
    }

    std::vector<temporalSequence>
    extractSequencesOfInterest(std::vector<dbMartEntry> &dbMart,  std::vector<size_t> &startPositions,
                              std::map<std::int64_t, size_t> &sequenceOfInterest, int numOfThreads,
                              double durationPeriod,
                              int daysForCoOccurrence, unsigned int phenxIdLength) {
        size_t numOfPatients = startPositions.size();
        size_t numberOfDbMartEntries = dbMart.size();
        std::vector<temporalSequence> localSequences[numOfThreads];
        omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, sequenceOfInterest, localSequences, durationPeriod, daysForCoOccurrence, daysPerWeek, daysPerMonth, phenxIdLength)
        for (size_t i = 0; i < numOfPatients; ++i) {
            size_t startPos = startPositions[i];
            size_t endPos = i < numOfPatients - 1 ? startPositions[i + 1] : numberOfDbMartEntries;

            for (size_t j = startPos; j < endPos - 1; ++j) {
                for (size_t k = j + 1; k < endPos; ++k) {
                    std::int64_t sequence = createSequence(dbMart[j].phenID, dbMart[k].phenID, phenxIdLength);
                    if (sequenceOfInterest.find(sequence) != sequenceOfInterest.end()) {
                        unsigned int duration = getDuration(dbMart[j].date, dbMart[k].date);
                        duration = getDurationPeriod(duration, durationPeriod, daysForCoOccurrence);
                        temporalSequence sequenceStruct = {sequence, duration, ((unsigned int) i)};
                        localSequences[omp_get_thread_num()].emplace_back(sequenceStruct);
                    }
                }
            }
        }
        std::cout << "merging sequencing vectors from all threads" << std::endl;
        std::vector<temporalSequence> allSequences;
        size_t sumOfSequences = 0;
        for (int i = 0; i < numOfThreads; ++i) {
            sumOfSequences += localSequences[i].size();
            allSequences.insert(allSequences.end(), localSequences[i].begin(), localSequences[i].end());
            localSequences[i].clear();
            localSequences[i].shrink_to_fit();
        }
        if (allSequences.size() != sumOfSequences) {
            std::cout << "Error during vector merging! Expected " << sumOfSequences
                      << " sequences, but allSequences stores "
                      << allSequences.size() << " sequences!" << std::endl;
        }
        return allSequences;
    }



    std::vector<temporalSequence>
    extractSparseSequences(std::vector<dbMartEntry> &dbMart, std::vector<size_t> &startPositions,
                           int numOfThreads, double durationPeriod, int daysForCoOccurrence, unsigned int phenxIdLength){
        size_t numOfPatients = startPositions.size();
        size_t numberOfDbMartEntries = dbMart.size();
        std::vector<temporalSequence> localSequences[numOfThreads];
        omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, localSequences, durationPeriod, daysForCoOccurrence, daysPerWeek, daysPerMonth,phenxIdLength)
        for (size_t i = 0; i < numOfPatients; ++i) {
            size_t startPos = startPositions[i];
            size_t endPos = i < numOfPatients - 1 ? startPositions[i + 1] : numberOfDbMartEntries;

            for (size_t j = startPos; j < endPos - 1; ++j) {
                for (size_t k = j + 1; k < endPos; ++k) {
                    std::int64_t sequence = createSequence(dbMart[j].phenID, dbMart[k].phenID, phenxIdLength);
                        unsigned int duration = getDuration(dbMart[j].date, dbMart[k].date);
                        duration = getDurationPeriod(duration, durationPeriod, daysForCoOccurrence);
                        temporalSequence sequenceStruct = {sequence, duration, ((unsigned int) i)};
                        localSequences[omp_get_thread_num()].emplace_back(sequenceStruct);
                }
            }
        }
        std::vector<temporalSequence> allSequences;
        size_t sumOfSequences = 0;
        for (int i = 0; i < numOfThreads; ++i) {
            sumOfSequences += localSequences[i].size();
            allSequences.insert(allSequences.end(), localSequences[i].begin(), localSequences[i].end());
            localSequences[i].clear();
            localSequences[i].shrink_to_fit();
        }
        if (allSequences.size() != sumOfSequences) {
            std::cout << "Error during vector merging! Expected " << sumOfSequences
                      << " sequences, but allSequences stores "
                      << allSequences.size() << " sequences!" << std::endl;
        }
        return allSequences;

    }


    unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration) {
        if (threshold == 0)
            return 0;
        if (duration == max)
            return (max - min - 1) / threshold;
        else
            return (duration - min) / threshold;
    }
}//tspm