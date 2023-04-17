//
// Created by jonas on 03.03.23.
//

#include "sequencing.h"
#include <parallel/algorithm>
#include "../lib/ips4o/ips4o.hpp"



size_t extractSequencesFromArray(dbMartEntry * dbMart, size_t numOfPatients, const size_t * startPositions,
                              size_t numberOfDbMartEntries,  const std::string& outPutDirectory,
                              const std::string& outputFilePrefix, int patIDLength, int numOfThreads){
    omp_set_num_threads(numOfThreads);
    size_t numOfSequences [numOfThreads]= { 0 };

    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, patIDLength, outPutDirectory,outputFilePrefix, numOfSequences)
    for(size_t i = 0; i < numOfPatients; ++i){
        size_t startPos = startPositions[i];
        size_t endPos = i < numOfPatients-1 ? startPositions[i+1] : numberOfDbMartEntries;
        size_t numOfPatientEntries = endPos - startPos;

        size_t numberOfSequences = (numOfPatientEntries * (numOfPatientEntries + 1)) / 2;
        std::vector<long> sequences;
        sequences.reserve(numberOfSequences);
        for(size_t j = startPos; j < endPos-1;++j){
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
createNonSparseTemporalSequences(dbMartEntry *dbMart, size_t numOfPatients, const size_t *startPositions,
                                 size_t numberOfDbMartEntries, std::map<long, size_t> nonSparseSequencesIDs,
                                 int numOfThreads, bool durationInWeeks, bool durationInMonths) {
    std::vector<temporalSequence> localSequences[numOfThreads];
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(numOfPatients, numberOfDbMartEntries, dbMart, startPositions, nonSparseSequencesIDs, localSequences, durationInMonths, durationInWeeks)
    for(size_t i = 0; i < numOfPatients; ++i){
        size_t startPos = startPositions[i];
        size_t endPos = i < numOfPatients-1 ? startPositions[i+1] : numberOfDbMartEntries;
        size_t numOfPatientEntries = endPos - startPos;
        unsigned patientId = dbMart[i].patID;
        size_t numberOfSequences = (numOfPatientEntries * (numOfPatientEntries + 1)) / 2;
        std::vector<long> sequences;
        sequences.reserve(numberOfSequences);
        for(size_t j = startPos; j < endPos -1;++j) {
            for (size_t k = j +1; k < endPos; ++k) {
                long sequence = createSequence(dbMart[j].phenID, dbMart[k].phenID);
                if (nonSparseSequencesIDs.find(sequence) != nonSparseSequencesIDs.end()) {
                    unsigned int duration = getDuration(dbMart[j].date, dbMart[k].date);
                    if(durationInWeeks && ! durationInMonths){
                        duration  = duration / 7;
                    }else if(durationInMonths && !durationInWeeks){
                        duration = duration / 30.437;
                    }
                    temporalSequence sequenceStruct = {sequence, duration, ((unsigned int) patientId)};
                    localSequences[omp_get_thread_num()].emplace_back(sequenceStruct);
                }
            }
        }
    }
    std::cout << "merging sequencing vectors from all threads" << std::endl;
    std::vector<temporalSequence> allSequences;
    size_t sumOfSequences = 0;
    for(int i = 0; i< numOfThreads; ++i){
        sumOfSequences += localSequences[i].size();
        allSequences.insert(allSequences.end(), localSequences[i].begin(), localSequences[i].end());
        localSequences[i].clear();
        localSequences->shrink_to_fit();
    }
    if(allSequences.size() != sumOfSequences){
        std::cout << "Error during vector merging! Expected " << sumOfSequences << " sequences, but allSequences stores "
        << allSequences.size() << " sequences!" <<std::endl;
    }
    return allSequences;
}


int sequenceWorkflow(bool temporal, bool removeSparseBuckets, const std::vector<std::string>& inputFilePaths, char inputFileDelimiter,
                     const std::string& outPutDirectory, const std::string& outputFilePrefix,
                     int patIDColumns[], int phenxColumns[], int dateColumns[], long maxPatID, double sparsity_value){
    //Todo sort files: execute preprocessing.sh script

    //===== extract sequences
    std::cout << "extracting transitive sequences" << std::endl;
    size_t sequenceCount = createSequencesFromFiles( inputFilePaths, inputFileDelimiter, outPutDirectory, outputFilePrefix,
    patIDColumns, phenxColumns, dateColumns, maxPatID);

    //===== remove sparse sequences
    std::cout << "removing sparse sequences" << std::endl;
    std::map<long, size_t> sequences = summarizeSequences(maxPatID, false, outPutDirectory,outputFilePrefix);
    size_t sparsityThreshold = maxPatID * sparsity_value;
    std::cout << "sparsity= " << sparsity_value << " sparsity threshold: " << sparsityThreshold <<std::endl;
    for(auto it = sequences.begin(); it != sequences.end();){
        if(it->second < sparsityThreshold){
            it = sequences.erase(it);
        } else{
            ++it;
        }
    }
    size_t numOfUniqueSequences = sequences.size();
    std::cout << "unique sequences: " << numOfUniqueSequences <<std::endl;
    if(!temporal){
        //TODO: Writeout reduced sequences
        return numOfUniqueSequences;
    }
    std::cout << "creating temporal sequences" << std::endl;
    //====== create buckets for each sequence
    std::vector<temporalSequence> temporalSequences = extractTemporalSequences(inputFilePaths, inputFileDelimiter,
                                                                             outPutDirectory, outputFilePrefix,
                                                                             patIDColumns, phenxColumns, dateColumns,
                                                                             maxPatID,sequenceCount,sequences,
                                                                             sparsityThreshold, removeSparseBuckets);
    writeSequencesAsCsV("test.out", outPutDirectory,'\t',temporalSequences.size(),temporalSequences.data());
    return 0;

}





std::vector<temporalSequence>
extractTemporalSequences(const std::vector<std::string> &inputFilePaths, char inputFileDelimiter,
                         const std::string &outPutDirectory, const std::string &outputFilePrefix, int *patIDColumns,
                         int *phenxColumns, int *dateColumns, long maxPatID, size_t sequenceCount,
                         const std::map<long, size_t>& sequences, size_t sparsityThreshold, bool removeSparseBuckets) {

    std::vector<temporalSequence> timedSequences = createSequencesWithDuration(inputFilePaths, inputFileDelimiter, outPutDirectory, outputFilePrefix,
                                                                               patIDColumns, phenxColumns, dateColumns, maxPatID,sequenceCount,sequences);
    std::cout << "created " <<timedSequences.size() << " sequences with duration. \nSorting:" << std::endl;
    std::vector<std::vector<temporalSequence>> globalSequences(omp_get_max_threads());
    ips4o::parallel::sort(timedSequences.begin(),timedSequences.end(),timedSequencesSorter);

    size_t numOfSeqs = timedSequences.size();
    std::mutex sequenceMutex;
    std::cout << "creating arrays for each thread" << std::endl;

    auto endPos = timedSequences.begin();
    for (size_t i = 0; i < omp_get_max_threads(); ++i){
        endPos = timedSequences.begin() +(timedSequences.size() / omp_get_max_threads() * (i+1));
        auto it = endPos;
        for (; it != timedSequences.end() && it->seqID == timedSequences.end()->seqID; ++it);
        endPos = it;
        globalSequences[i] = std::vector<temporalSequence>(timedSequences.begin(), endPos);
        timedSequences.erase(timedSequences.begin(),endPos);
        timedSequences.shrink_to_fit();
    }
    timedSequences.clear();
    timedSequences.shrink_to_fit();
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

#pragma omp parallel for shared(sequences,sequenceMutex, removeSparseBuckets, sparsityThreshold, globalSequences, std::cout, startIndices, numOfSeqs) default (none)
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

    timedSequences.clear();
    timedSequences.shrink_to_fit();
    std::cout << "merging and sorting sequence vectors" << std::endl;
    std::vector<temporalSequence> sortedSequences;
    for(std::vector<temporalSequence> sequences : globalSequences){
        sortedSequences.insert(sortedSequences.end(),sequences.begin(), sequences.end());
        sequences.clear();
        sequences.shrink_to_fit();
    }
    ips4o::parallel::sort(timedSequences.begin(),timedSequences.end(),timedSequencesSorter);

    return sortedSequences;
}

std::pair<size_t, size_t> findStartAndEnd(std::vector<temporalSequence> sequences, const unsigned long &sequence) {

    size_t start = 0, end = 0;
    bool lookforStart = true;
    bool lookForEnd = false;
    for (size_t i = 0; i < sequences.size(); ++i){
        if(lookforStart && sequences[i].seqID==sequence){
            start = i;
            lookforStart = false;
            lookForEnd = true;
        }
        if(lookForEnd && sequences[i].seqID!=sequence){
            end = i;
            return std::make_pair(start,end);
        }
    }


    return std::make_pair(start,end);
}


long writeSequencesAsCsV(std::string fileName, std::string filepath, char delimiter, size_t numOfSequences, temporalSequence * temporalSequences, bool debug){
    FILE* sequenceFile;
    sequenceFile = fopen((filepath.append(fileName)).c_str(), "w");
    long written = 0;

    if(sequenceFile == nullptr) {
        exit(EXIT_FAILURE);
    }
    if(debug) {
        for (int i = 0; i < numOfSequences; ++i) {
            std::string out = std::to_string(temporalSequences[i].patientID).append(1, delimiter)
                    .append(std::to_string((temporalSequences[i].seqID << 8) >> 8)).append(1, delimiter)
                    .append(std::to_string(temporalSequences[i].seqID >> 63)).append(1, delimiter)
                    .append(std::to_string(temporalSequences[i].duration)).append(1, '\n');
            written += fwrite(out.c_str(), sizeof(char), out.length(), sequenceFile);
        }
    }else {
        for (int i = 0; i < numOfSequences; ++i) {
            std::string out = std::to_string(temporalSequences[i].patientID).append(1, delimiter)
                    .append(std::to_string(temporalSequences[i].seqID)).append(1, '\n');
            written += fwrite(out.c_str(), sizeof(char), out.length(), sequenceFile);
        }
    }

    fclose(sequenceFile);
    return written;


}
unsigned int getBucket(unsigned int min, unsigned int max, int threshold, unsigned int duration) {
    if (threshold == 0)
        return 0;
    if (duration == max)
        return (max-min-1)/threshold;
    else
        return (duration-min)/threshold;
}

std::vector<temporalSequence> createSequencesWithDuration(std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                                          const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                                          int patIDColumns[], int phenxColumns[], int dateColumns[], long maxPatID, size_t sequenceCount,
                                                          const std::map<long,size_t>& sequenceMap){
    std::vector<temporalSequence> sequencesWithDuration;
    int numOfProcs = omp_get_max_threads();
    std::vector<temporalSequence> localSequences[numOfProcs];
    for(std::vector<temporalSequence> vec : localSequences){
        vec.reserve((sequenceCount/numOfProcs*1.1));
    }
    std::vector<size_t> startPositions;
    size_t sequenceNum = 0;
    std::mutex readMutex;
    std::cout << "extracting temporal transitive sequences" << std::endl;
    for (int i = 0; i < inputFilePaths.size(); ++i){
        std::string filePath = inputFilePaths[i];
        std::pair<size_t, size_t>linesAndPatientInFile = determinePatientStartPositionsInFile(filePath,
                                                                                              inputFileDelimiter,
                                                                                              &startPositions);
        size_t patientCount = linesAndPatientInFile.second;

        FILE *csvFilePointer = fopen(filePath.c_str(), "r");
        if (csvFilePointer == nullptr) {
            exit(EXIT_FAILURE);
        }
        // read in header line before iterating over all files
        char line[2048];
        size_t len = 2048;
        if (fgets(line, len, csvFilePointer) == NULL){
            return sequencesWithDuration;
        };
//   ====== Sequence Patients =======
        int patientId = 0;
        int local_patID;

#pragma omp parallel for private(local_patID) default (none) shared(sequenceMap, localSequences, i, patientId, csvFilePointer,outputFilePrefix,outPutDirectory, patIDColumns,phenxColumns, dateColumns, patientCount, readMutex, startPositions)
        for (size_t j = 0; j < patientCount; ++j) {
            std::vector<dbMartEntry> patientEntries;
            readMutex.lock();
            patientEntries = extractPatient(csvFilePointer,
                                            &startPositions, patientId, patIDColumns[i], phenxColumns[i], dateColumns[i]);
            local_patID = patientId;
            ++patientId;
            readMutex.unlock();

            for (int k = 0; k < patientEntries.size() - 1; ++k) {
                for (int l = k + 1; l < patientEntries.size(); ++l) {
                    long sequence = createSequence(patientEntries[k].phenID, patientEntries[l].phenID);
                    if(sequenceMap.find(sequence)!= sequenceMap.end()) {
                        unsigned int duration = getDuration(patientEntries[k].date, patientEntries[l].date);
                        temporalSequence sequenceStruct = {sequence, duration,((unsigned int) patientId)};
                        localSequences[omp_get_thread_num()].emplace_back(sequenceStruct);
                    }
                }
            }
        }
        fclose(csvFilePointer);
    }

    size_t tmp =0;
    for(int i = 0; i < numOfProcs; ++i){
       tmp += localSequences[i].size();
    }

    std::cout << "merging vectors" << std::endl;
    for(int i = 0; i < numOfProcs; ++i){
        sequencesWithDuration.insert(sequencesWithDuration.end(),localSequences[i].begin(),localSequences[i].end());
        localSequences[i].clear();
        localSequences->shrink_to_fit();
    }

    std::cout<<"sequences bevor mergering: " << tmp << " sequences after merging: " << sequencesWithDuration.size()<<std::endl;

    return sequencesWithDuration;
}


size_t createSequencesFromFiles (std::vector<std::string> inputFilePaths, char inputFileDelimiter,
                                 const std::string& outPutDirectory, const std::string& outputFilePrefix,
                                 int patIDColumns[], int phenxColumns[], int dateColumns[], long maxPatID){
    size_t sequenceNum = 0;
    std::mutex readMutex;
    std::mutex countMutex;
    int patIDLength = 7;
    std::vector<size_t> startPositions;
    for (int i = 0; i < inputFilePaths.size(); ++i){
        size_t local_NumOfSequences = 0;
        std::string filePath = inputFilePaths[i];
        std::pair<size_t, size_t>linesAndPatientInFile = determinePatientStartPositionsInFile(filePath,
                                                                                              inputFileDelimiter,
                                                                                              &startPositions);

        size_t patientCount = linesAndPatientInFile.second;
        size_t lineCount = linesAndPatientInFile.first;

        FILE *csvFilePointer = fopen(filePath.c_str(), "r");
        if (csvFilePointer == nullptr) {
            exit(EXIT_FAILURE);
        }
        // read in header line before iterating over all files
        char line[2048];
        size_t len = 2048;
        if(fgets(line, len, csvFilePointer) == NULL){
            return 0;
        }
        std::cout << "extracting sequences" << std::endl;
//   ====== Sequence Patients =======
        int patientId = 0;
        int local_patID;

#pragma omp parallel for private(local_patID) default (none) shared(i, patientId, patIDLength, csvFilePointer,outputFilePrefix,outPutDirectory, patIDColumns,phenxColumns, dateColumns, patientCount, readMutex, countMutex,local_NumOfSequences, startPositions)
        for (size_t j = 0; j < patientCount; ++j) {
            std::vector<dbMartEntry> patientEntries;
            readMutex.lock();
            patientEntries = extractPatient(csvFilePointer,
                                            &startPositions, patientId, patIDColumns[i], phenxColumns[i],
                                            dateColumns[i]);
            local_patID = patientId;
            ++patientId;
            readMutex.unlock();

            //the number sequences is equal to the gaussian sum, initalize vector with size to avoid re
            size_t numberOfSequences = (patientEntries.size() * (patientEntries.size() + 1)) / 2;
            std::vector<long> sequences;
            sequences.reserve(numberOfSequences);

            for (int k = 0; k < patientEntries.size() - 1; ++k) {
                for (int l = k + 1; l < patientEntries.size(); ++l) {
                    sequences.emplace_back(createSequence(patientEntries[k].phenID, patientEntries[l].phenID));
                }
            }

            countMutex.lock();
            local_NumOfSequences += sequences.size();
            countMutex.unlock();

            std::string patIDString = std::to_string(local_patID);
            patIDString.insert(patIDString.begin(), patIDLength - patIDString.size(), '0');
            std::string patientFileName = std::string(outPutDirectory).append(outputFilePrefix).append(patIDString);
            writeSequencesToBinaryFile(patientFileName, sequences);
        }

        sequenceNum+=local_NumOfSequences;
        std::cout << "file: " << inputFilePaths[i] << std::endl;
        std::cout << "number of sequences: " << sequenceNum << std::endl;
        std::cout << "number of patients: " << patientCount<< std::endl;
        std::cout << "lines in file: " << lineCount << std::endl;
        fclose(csvFilePointer);
    }
    return sequenceNum;
}