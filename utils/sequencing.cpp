//
// Created by jonas on 03.03.23.
//

#include "sequencing.h"
#include <parallel/algorithm>
#include "../lib/ips4o/ips4o.hpp"
#include <string>


std::vector<temporalSequence> extractTemporalBuckets(std::vector<dbMartEntry> &dbMart, size_t numOfPatients,
                                                     const size_t *startPositions,
                                                     std::map<long, size_t> &nonSparseSequencesIDs, int numOfThreads,
                                                     bool durationInWeeks, bool durationInMonths,
                                                     size_t sparsityThreshold, bool removeSparseBuckets){

    std::vector<temporalSequence> nonSparseSequences =
            extractNonSparseSequences(dbMart, numOfPatients, startPositions,
                                      nonSparseSequencesIDs,
                                      numOfThreads, durationInWeeks, durationInMonths);


    // split sequence vector in subvectors! //
    std::vector<std::vector<temporalSequence>> globalSequences(numOfThreads);
    globalSequences = splitSequenceVectorInChunkes(nonSparseSequences, numOfThreads);


    std::cout << "creating arrays for start indicies" << std::endl;
    std::vector<std::vector<size_t>> startIndices(omp_get_max_threads());
    for (size_t i = 0; i < omp_get_max_threads(); ++i){
        startIndices[i] =std::vector<size_t>();
        if(globalSequences[i].size() == 0){
            continue;
        }
        startIndices.emplace_back(0);
        unsigned long seq = globalSequences[i][0].seqID;
        for (int j = 0; j <globalSequences[i].size() ; ++j) {
            if(seq != globalSequences[i][j].seqID){
                seq = globalSequences[i][j].seqID;
                startIndices[i].emplace_back(j);
            }
        }
    }
    std::cout << "computing temporal buckets" << std::endl;

#pragma omp parallel for shared(nonSparseSequences, removeSparseBuckets, sparsityThreshold, globalSequences, std::cout, startIndices) default (none)
    for(size_t i = 0; i < omp_get_max_threads(); ++i) {
        for (size_t j = 0; j < startIndices[i].size(); ++j) {
            size_t local_start = startIndices[i][j];
            size_t end;
            if (j < startIndices[i].size() - 1) {
                end = startIndices[i][j + 1];
            } else {
                end = globalSequences[i].size();
            }
            auto startIterator = globalSequences[i].begin() + local_start;
            auto endIterator = globalSequences[i].begin() + end;

//       compute buckets
            unsigned int min = 0, max = 0;

            for (auto it = startIterator; it != endIterator; ++it) {
                min = std::min(min, it->duration);
                max = std::max(max, it->duration);
            }
            unsigned int range = max - min;
            int numberOfBuckets = std::max(4, int(std::log(range)));

//      add bucket id to sequence id
            for (auto it =startIterator; it != endIterator; ++it) {
                unsigned int bucket = getBucket(min, max, range / numberOfBuckets, it->duration);
                int bitShift = sizeof(long) * 8 - 8;
                it->seqID = bucket << bitShift && it->seqID;
            }

            if (removeSparseBuckets) {
                std::set<unsigned int> sequenceInPatient;
                ips4o::parallel::sort(startIterator, endIterator,timedSequencesSorter);
                unsigned long lastSequence = startIterator->seqID;
                size_t count = 0;
                auto it = startIterator;
                while (it != endIterator) {
                    if (it->seqID == lastSequence) {
                        sequenceInPatient.insert(it->patientID);
                        ++count;
                        ++it;
                    } else {
                        lastSequence = it->seqID;
                        if (sequenceInPatient.size() < sparsityThreshold) {
                            it = globalSequences[i].erase(it - count, it);
                        } else {
                            ++it;
                        }
                        sequenceInPatient.clear();
                        count = 0;
                    }
                }
            }
        }
    }

    std::cout << "merging and sorting sequence vectors" << std::endl;
    std::vector<temporalSequence> sortedSequences;
    for(std::vector<temporalSequence> sequences : globalSequences){
        sortedSequences.insert(sortedSequences.end(),sequences.begin(), sequences.end());
        sequences.clear();
        sequences.shrink_to_fit();
    }
    ips4o::parallel::sort(sortedSequences.begin(),sortedSequences.end(),timedSequencesSorter);

    return sortedSequences;

}


size_t createSequencesFromFiles (std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                 const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                 int patIDColumns[], int phenxColumns[], int dateColumns[], size_t numOfPatients, int  patIdLength, int numOfThreads){
    std::vector<dbMartEntry> dbMart;
    for(int i = 0; i < inputFilePaths.size();++i){
        FILE *csvFilePointer = fopen(inputFilePaths[i].c_str(), "r");
        if (csvFilePointer == nullptr) {
            exit(EXIT_FAILURE);
        }
        std::vector<dbMartEntry> localDBMart = extractDBMartFromCsv(csvFilePointer,patIDColumns[i], phenxColumns[i], dateColumns [i], inputFileDelimiter);
        dbMart.insert(dbMart.end(), localDBMart.begin(), localDBMart.end());
    }
    std::vector<size_t> startPositions = extractStartPositions(dbMart);
    return extractSequencesFromArray(dbMart, numOfPatients, startPositions.data(), outPutDirectory, outputFilePrefix, patIdLength, numOfThreads);
}

std::vector<size_t> extractStartPositions(std::vector<dbMartEntry> &dbMart) {
    ips4o::parallel::sort(dbMart.begin(), dbMart.end(), dbMartSorter);
    std::vector<size_t> startPositions;
    startPositions.emplace_back(0);
    for (int i = 1; i < dbMart.size(); ++i) {
        if(dbMart[i].patID != dbMart[i-1].patID){
            startPositions.emplace_back(i);
        }
    }
    return startPositions;
}


size_t extractSequencesFromArray(std::vector<dbMartEntry> &dbMart, size_t numOfPatients, const size_t * startPositions,
                              const std::string& outPutDirectory,const std::string& outputFilePrefix, int patIDLength, int numOfThreads){
    omp_set_num_threads(numOfThreads);
    size_t numOfSequences [numOfThreads]= { 0 };
    size_t numberOfDbMartEntries = dbMart.size();
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, patIDLength, outPutDirectory,outputFilePrefix, numOfSequences)
    for(size_t i = 0; i < numOfPatients; ++i){
        size_t startPos = startPositions[i];
        size_t endPos = i < numOfPatients - 1 ? startPositions[i+1] : numberOfDbMartEntries;
        size_t numOfPatientEntries = endPos - startPos;

        size_t numberOfSequences = (numOfPatientEntries * (numOfPatientEntries - 1)) / 2;
        std::vector<long> sequences;
        sequences.reserve(numberOfSequences);
        for(size_t j = startPos; j < endPos - 1; ++j){
            for (size_t k = j + 1; k < endPos ; ++k) {
                sequences.emplace_back(createSequence(dbMart[j].phenID, dbMart[k].phenID));
            }
        }
        numOfSequences [omp_get_thread_num()] += sequences.size();
        std::string patIDString = std::to_string(i);
        patIDString.insert(patIDString.begin(), patIDLength - patIDString.size(), '0');
        std::string patientFileName = std::string(outPutDirectory).append(outputFilePrefix).append(patIDString);
        writeSequencesToBinaryFile(patientFileName, sequences);
    }
    size_t sumOfSequences = 0;
    for(int i = 0; i< numOfThreads; ++i){
        sumOfSequences += numOfSequences[i];
    }
    return sumOfSequences;
}



std::vector<size_t>getStartPositionsFromSequenceVector(std::vector<temporalSequence> &sequence){
    std::vector<size_t> startPositions;
    unsigned int lastPatId = sequence[0].patientID;
    for (int i = 1; i < sequence.size(); ++i) {
        if(lastPatId!=sequence[i].patientID){
            startPositions.emplace_back(i);
            lastPatId = sequence[i].patientID;
        }
    }
    return startPositions;
}


std::vector<std::vector<temporalSequence>> splitSequenceVectorInChunkes(std::vector<temporalSequence> &sequences, unsigned int chunks){
    ips4o::parallel::sort(sequences.begin(), sequences.end(), timedSequencesSorter);
    std::vector<std::vector<temporalSequence>> localSequences;
    //    Split Sequences in sub vectors for paralle access
    auto endPos = sequences.begin();
    size_t numOfSequencesPerChunk = sequences.size() / chunks;
    for (size_t i = 0; i < chunks; ++i){
        endPos = sequences.begin() + numOfSequencesPerChunk;
        auto it = endPos;
        for (; it != sequences.end() && it->seqID == sequences.end()->seqID; ++it);
        endPos = it;
        localSequences.emplace_back(std::vector<temporalSequence>(sequences.begin(), endPos));
        sequences.erase(sequences.begin(),endPos);
        sequences.shrink_to_fit();
    }
    //for small dataset, we might have some sequences (< chunksize) left at the end
    localSequences[chunks - 1].insert(localSequences[chunks - 1].end(),sequences.begin(),sequences.end());
    sequences.clear();
    sequences.shrink_to_fit();
    return localSequences;
}



std::vector<temporalSequence> extractMonthlySequences(std::vector<temporalSequence> &sequences, bool durationSparsity,
                                                      double sparsity, size_t numOfPatients, int numOfThreads,
                                                      double durationPeriods, unsigned int daysForCoOoccurence, unsigned int bitShift){
    std::vector<std::vector<temporalSequence>> localSequences = splitSequenceVectorInChunkes(sequences, numOfThreads);
    size_t sparsityThreshold = numOfPatients * sparsity;

#pragma omp parallel for default (none) shared(numOfThreads, localSequences, bitShift, sparsityThreshold, durationSparsity, durationPeriods, daysForCoOoccurence)
    for(size_t i = 0; i < numOfThreads; ++i) {
        if(localSequences[i].empty())
            continue;
        auto sparsityIt = localSequences[i].begin();
        std::set<unsigned int> sequenceInPatient;
        unsigned long lastSequence = getDurationPeriod(sparsityIt->duration, durationPeriods, daysForCoOoccurence)<< bitShift & sparsityIt->seqID;;
        size_t count = 0;
        // since the sequences are order by ID as first and duration as second iterator, we can iterate over them and
        // integrate the duration in month and directly remove check for sparsity if the updated sequenceID change
        while (sparsityIt != localSequences[i].end()) {
            sparsityIt->seqID = getDurationPeriod(sparsityIt->duration, durationPeriods, daysForCoOoccurence)<< bitShift | sparsityIt->seqID;
            if (sparsityIt->seqID == lastSequence) {
                sequenceInPatient.insert(sparsityIt->patientID);
                ++count;
                ++sparsityIt;
            } else {
                lastSequence = sparsityIt->seqID;
                if (durationSparsity && sequenceInPatient.size() < sparsityThreshold) {
                    sparsityIt = localSequences[i].erase(sparsityIt - count, sparsityIt);
                } else {
                    ++sparsityIt;
                }
                sequenceInPatient.clear();
                count = 0;
            }
        }
    }
    std::vector<temporalSequence> mergedSequences;
    for(std::vector<temporalSequence> sequences : localSequences){
        mergedSequences.insert(mergedSequences.end(),sequences.begin(), sequences.end());
        sequences.clear();
        sequences.shrink_to_fit();
    }
    return mergedSequences;
}

unsigned int getDurationPeriod(unsigned int duration, double durationPeriods, int daysForCoOoccurence) {
    //if the duration is less than 2 weeks it is considered as co-occurrence and therefore the distance is 0
    if (daysForCoOoccurence <= 0)
        daysForCoOoccurence = 1;
    if(duration <  daysForCoOoccurence){
        return 0;
    }
    if(durationPeriods == DURATION_IN_DAYS)
        return duration;
    else
        return std::ceil(duration/durationPeriods);
}

std::vector<temporalSequence>
extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, size_t numOfPatients, const size_t *startPositions,
                          std::map<long, size_t> &nonSparseSequencesIDs, int numOfThreads, bool durationInWeeks,
                          bool durationInMonths) {
    size_t numberOfDbMartEntries = dbMart.size();
    std::vector<temporalSequence> localSequences[numOfThreads];
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, nonSparseSequencesIDs, localSequences, durationInMonths, durationInWeeks, daysPerWeek, daysPerMonth)
    for(size_t i = 0; i < numOfPatients; ++i){
        size_t startPos = startPositions[i];
        size_t endPos = i < numOfPatients-1 ? startPositions[i+1] : numberOfDbMartEntries;

        for(size_t j = startPos; j < endPos - 1; ++j) {
            for (size_t k = j + 1; k < endPos; ++k) {
                long sequence = createSequence(dbMart[j].phenID, dbMart[k].phenID);
                if (nonSparseSequencesIDs.find(sequence) != nonSparseSequencesIDs.end()) {
                    unsigned int duration = getDuration(dbMart[j].date, dbMart[k].date);
                    if(durationInWeeks && ! durationInMonths){
                        duration  = duration / daysPerWeek;
                    }else if(durationInMonths && !durationInWeeks){
                        duration = duration / daysPerMonth;
                    }
                    temporalSequence sequenceStruct = {sequence, duration, ((unsigned int) i)};
                    localSequences[omp_get_thread_num()].emplace_back(sequenceStruct);
                }
            }
        }
    }
    std::cout << "merging sequencing vectors from all threads" << std::endl;
    std::vector<temporalSequence> allSequences;
    size_t sumOfSequences = 0;
    for(int i = 0; i < numOfThreads; ++i){
        sumOfSequences += localSequences[i].size();
        allSequences.insert(allSequences.end(), localSequences[i].begin(), localSequences[i].end());
        localSequences[i].clear();
        localSequences[i].shrink_to_fit();
    }
    if(allSequences.size() != sumOfSequences){
        std::cout << "Error during vector merging! Expected " << sumOfSequences << " sequences, but allSequences stores "
        << allSequences.size() << " sequences!" <<std::endl;
    }
    return allSequences;
}


unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration) {
    if (threshold == 0)
        return 0;
    if (duration == max)
        return (max-min-1)/threshold;
    else
        return (duration-min)/threshold;
}
