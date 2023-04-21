//
// Created by jonas on 03.03.23.
//

#include "sequencing.h"
#include <parallel/algorithm>
#include "../lib/ips4o/ips4o.hpp"
#include <string>


std::vector<temporalSequence> extractTemporalSequences(std::vector<dbMartEntry> &dbMart, size_t numOfPatients,
                                                       const size_t *startPositions,
                                                       std::map<long, size_t> &nonSparseSequencesIDs, int numOfThreads,
                                                       bool durationInWeeks, bool durationInMonths,
                                                       size_t sparsityThreshold, bool removeSparseBuckets){

    std::vector<temporalSequence> nonSparseSequences =
            extractNonSparseSequences(dbMart, numOfPatients, startPositions,
                                      nonSparseSequencesIDs,
                                      numOfThreads, durationInWeeks, durationInMonths);
    std::vector<std::vector<temporalSequence>> globalSequences(omp_get_max_threads());
    ips4o::parallel::sort(nonSparseSequences.begin(),nonSparseSequences.end(),timedSequencesSorter);

    // split sequence vector in subvectors!
    auto endPos = nonSparseSequences.begin();
    for (size_t i = 0; i < omp_get_max_threads(); ++i){
        endPos = nonSparseSequences.begin() + (nonSparseSequences.size() / omp_get_max_threads() * (i + 1));
        auto it = endPos;
        for (; it != nonSparseSequences.end() && it->seqID == nonSparseSequences.end()->seqID; ++it);
        endPos = it;
        globalSequences[i] = std::vector<temporalSequence>(nonSparseSequences.begin(), endPos);
        nonSparseSequences.erase(nonSparseSequences.begin(),endPos);
        nonSparseSequences.shrink_to_fit();
    }
    nonSparseSequences.clear();
    nonSparseSequences.shrink_to_fit();
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
                startIndices[i].emplace_back(j-1);
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

std::vector<temporalSequence>
extractNonSparseSequences(std::vector<dbMartEntry> &dbMart, size_t numOfPatients, const size_t *startPositions,
                          std::map<long, size_t> &nonSparseSequencesIDs, int numOfThreads, bool durationInWeeks,
                          bool durationInMonths) {
    size_t numberOfDbMartEntries = dbMart.size();
    std::vector<temporalSequence> localSequences[numOfThreads];
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, nonSparseSequencesIDs, localSequences, durationInMonths, durationInWeeks)
    for(size_t i = 0; i < numOfPatients; ++i){
        size_t startPos = startPositions[i];
        size_t endPos = i < numOfPatients-1 ? startPositions[i+1] : numberOfDbMartEntries;

        for(size_t j = startPos; j < endPos - 1; ++j) {
            for (size_t k = j + 1; k < endPos; ++k) {
                long sequence = createSequence(dbMart[j].phenID, dbMart[k].phenID);
                if (nonSparseSequencesIDs.find(sequence) != nonSparseSequencesIDs.end()) {
                    unsigned int duration = getDuration(dbMart[j].date, dbMart[k].date);
                    if(durationInWeeks && ! durationInMonths){
                        duration  = duration / 7;
                    }else if(durationInMonths && !durationInWeeks){
                        duration = duration / 30.437;
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
