#include "utils.h"
#include "sorter.h"


std::vector<dbMartEntry> extractDBMartFromCsv(FILE *csv_file, int patIdColumn, int phenotypeIDColumn,
                                              int dateColumn, char delim){
    if (csv_file == nullptr) {
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> lines;
    std::vector<dbMartEntry> dbMartEntries;
    char line[2048];
    size_t len = 2048;
    fgets(line, len, csv_file);
    while(fgets(line, len, csv_file) != nullptr){
        lines = getTokensFromLine(std::string(line), delim);
        dbMartEntry entry;
        entry.patID = atoi(lines[patIdColumn].c_str());
        entry.phenID = atoi(lines[phenotypeIDColumn].c_str());
        entry.date = getTimeFromString(lines[dateColumn].c_str());
        dbMartEntries.emplace_back(entry);
    }
    return  dbMartEntries;


}


size_t writeSequencesToBinaryFile(std::string patientFilename, std::vector<std::int64_t> sequences){
    FILE* patientFile;
    patientFile = fopen(patientFilename.c_str(), "wb");
    size_t written;
    if(patientFile!= nullptr) {
       written = std::fwrite(&sequences[0], 1, sequences.size() * sizeof(std::int64_t), patientFile);
    }else{
        return 0;
    }
    fclose(patientFile);
    return written;
}

size_t writeSequencesToFile(std::string patientFilename, std::vector<std::int64_t> sequences){
    FILE* patientFile;
    patientFile = fopen(patientFilename.append("asChar").c_str(), "w");

    size_t written;
    for(std::int64_t seq : sequences){
        std::string s = std::to_string(seq).append("\n");
        written += std::fwrite(s.c_str(), 1,  s.size(), patientFile);
    }
    fclose(patientFile);
    return written;
}

std::int64_t createSequence(int phenotypeA, int phenotypeB, int phenotypelenght){
    std::string phenA = std::to_string(phenotypeA);
    std::string phenB = std::to_string(phenotypeB);
    phenB.insert(phenB.begin(),phenotypelenght-phenB.size(), '0');
    phenA.append(phenB);
    return atol(phenA.c_str());
}

std::int64_t getTimeFromString(const char * date_string) {
    std::istringstream stringStream(date_string);
    std::tm time = {};
    stringStream >> std::get_time(&time,"%Y-%m-%dT%H:%M:%SZ");

    std::time_t  time_stamp = std::mktime(&time);
    return time_stamp;
}

unsigned int getDuration(std::int64_t startDate, std::int64_t endDate) {
    unsigned int secondsPerDay = 60 * 60 * 24;
    unsigned int duration = std::abs(endDate - startDate) / secondsPerDay;
    return duration;
}

std::vector<std::string> getTokensFromLine(const std::string &line, char delim) {
    std::vector<std::string> vectorizedLine;
    size_t start;
    size_t end = 0;


    while ((start = line.find_first_not_of(delim, end)) != std::string::npos) {
        end = line.find(delim, start);
        vectorizedLine.emplace_back(line.substr(start, end - start));
    }
    return vectorizedLine;
}

std::int64_t getFileSize(const std::string& filename){
    try {
        return std::filesystem::file_size(filename);
    } catch(std::filesystem::filesystem_error& e) {
        std::cout << e.what() << '\n';
        return 0;
    }
}

std::map<std::int64_t, size_t>
summarizeSequencesFromFiles(const std::string &outputDir, const std::string &file_prefix, int numberOfPatients,
                            bool storesDuration, unsigned int patIdLength, unsigned int bitShift) {

    std::map<std::int64_t, size_t> globalSequenceMap;
    const int numOfProcs = omp_get_max_threads();
    std::map<std::int64_t, size_t> localmaps[numOfProcs];
    std::vector<int> patientsPerThread(numOfProcs);
    std::mutex map_mutex;
#pragma omp parallel for default (none) shared(numberOfPatients, file_prefix, outputDir,storesDuration,localmaps,patientsPerThread,map_mutex, globalSequenceMap, patIdLength, bitShift)
    for (size_t i = 0; i < numberOfPatients; ++i) {
        std::string patIDString = std::to_string(i);
        patIDString.insert(patIDString.begin(), patIdLength - patIDString.size(), '0');
        std::string patientFileName = std::string(outputDir).append(file_prefix).append(patIDString);
        std::int64_t numberOfSequences = getFileSize(patientFileName)/sizeof(std::int64_t);
        FILE* patientFile = fopen(patientFileName.c_str(),"rb");

        std::set<std::int64_t> patientSequenceSet;
        for (int j = 0; j < numberOfSequences; ++j) {
            std::int64_t sequence = 0;
            fread(&sequence, sizeof(std::int64_t),1,patientFile);
            if(storesDuration){
                sequence = (sequence<<bitShift)>>bitShift;
            }
            if(patientSequenceSet.insert(sequence).second){
                ++localmaps[omp_get_thread_num()][sequence];
            }

        }

        fclose(patientFile);
        ++patientsPerThread[omp_get_thread_num()];
        if(patientsPerThread[omp_get_thread_num()]>=50) {
            map_mutex.lock();
            for (std::pair<std::int64_t, size_t> entry: localmaps[omp_get_thread_num()]) {
                globalSequenceMap[entry.first] += entry.second;
            }
            map_mutex.unlock();
            localmaps[omp_get_thread_num()].clear();
            patientsPerThread[omp_get_thread_num()] = 0;
        }
    }

    //merge all local maps into final map
    for (int i = 0; i < numOfProcs; ++i) {
        for (std::pair<std::int64_t, size_t> entry: localmaps[i]) {
            globalSequenceMap[entry.first] += entry.second;
        }
    }
    return globalSequenceMap;
}

std::map<std::int64_t, size_t>
summarizeSequencesFromDbMart(std::vector<dbMartEntry> &dbMart, std::vector<size_t> startPositions,
                            unsigned int &numOfThreads, unsigned  int phenxIdLength) {
    std::map<std::int64_t, size_t> globalSequenceMap;
    std::map<std::int64_t, size_t> localmaps[numOfThreads];
    size_t numOfPatients = startPositions.size();
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(startPositions,dbMart,localmaps, globalSequenceMap, phenxIdLength, numOfPatients)
    for (size_t i = 0; i < numOfPatients; ++i) {
        size_t startPos = startPositions[i];
        size_t endPos = i < numOfPatients - 1 ? startPositions[i+1] : dbMart.size();
        std::set<std::int64_t> patientSequenceSet;
        for(size_t j = startPos; j < endPos - 1; ++j){
            for (size_t k = j + 1; k < endPos ; ++k) {
                std::int64_t sequence = createSequence(dbMart[j].phenID, dbMart[k].phenID, phenxIdLength);
                if(patientSequenceSet.insert(sequence).second){
                    if(auto it = localmaps[omp_get_thread_num()].find(sequence); it == localmaps[omp_get_thread_num()].end()){
                        localmaps[omp_get_thread_num()].insert(std::make_pair(sequence,1));
                    }else{
                        it->second += 1;
                    }
                }
            }
        }
    }

    for(int i = 0; i < numOfThreads; ++i){
        if(localmaps[i].empty()){
            continue;
        }
        for(auto mapEntry :localmaps[i]) {
            if(auto it = globalSequenceMap.find(mapEntry.first); it == globalSequenceMap.end()){
                globalSequenceMap.insert(std::make_pair(mapEntry.first,mapEntry.second));
            }else{
                it->second += mapEntry.second;
            }
        }
    }
    return globalSequenceMap;
}


std::int64_t writeSequencesAsCsV(std::string fileName, std::string filepath, char delimiter, size_t numOfSequences, temporalSequence * temporalSequences, bool debug){
    FILE* sequenceFile;
    sequenceFile = fopen((filepath.append(fileName)).c_str(), "w");
    std::int64_t written = 0;

    if(sequenceFile == nullptr) {
        return -1;
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

std::filesystem::path createOutputFilePath(const std::string &outPutDirectory) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string stringPath = outPutDirectory;
    stringPath.append(oss.str()).append("/");
    std::filesystem::path outputPath = std::filesystem::u8path(stringPath);
    std::filesystem::create_directories(outputPath);
    return outputPath;
}


/**
 * Function to identify possible candidate phenx from the sequences that might be of interest given a list of specific phenx
 * @param originalSequences the original sequences, the duration should not be stored in the identifier
 * @param minDuration the minimal duration that at least one sequence of each unique sequence starting with a phenx of interestshould have
 * @param bitShift
 * @param lengthOfPhenx the length defined for each phenx when combining then into a sequence
 * @param phenxOfInterest an array, containing all starting phenx ids
 * @param numOfThreads the number of threads for paralleisation
 */
std::set<unsigned int> extractEndPhenxWithGivenStartPhenx(std::vector<temporalSequence> &originalSequences, std::uint64_t minDuration,
                                                          unsigned int bitShift, unsigned int lengthOfPhenx,
                                                          std::vector<unsigned int> &phenxOfInterest, int &numOfThreads){
    std::set<unsigned int> candidatePhenxs[numOfThreads];

    ips4o::parallel::sort(originalSequences.begin(),originalSequences.end(),
                          timedSequencesSorter, numOfThreads);
    std::vector<size_t>startPositions = getSequenceStartPositions(originalSequences);
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(startPositions,originalSequences,candidatePhenxs, minDuration, lengthOfPhenx, phenxOfInterest )
    for (size_t i = 0; i < startPositions.size(); ++i) {
        size_t startPos = startPositions [i];
        size_t endPos;
        if(i < startPositions.size() - 1){
            endPos = startPositions[i+1];
        }else{
            endPos = originalSequences.size();
        }
        unsigned int startPhenx = getStartPhenx(originalSequences[startPos],lengthOfPhenx);
        if(isPhenxOfInterest( startPhenx, phenxOfInterest)){
            unsigned int mininmalDuration = 0;
            unsigned int endPhenx = getEndPhenx(originalSequences[startPos],lengthOfPhenx);
            for(size_t j = startPos; j < endPos; ++j ){
                temporalSequence seq = originalSequences[j];
                mininmalDuration = std::max(mininmalDuration,  seq.duration);
            }
            if(mininmalDuration >= minDuration){
                candidatePhenxs[omp_get_thread_num()].insert(endPhenx);
            }
        }
    }
    std::set<unsigned int> candidatePhenx;
    for (int i = 0; i < numOfThreads; ++i) {
        candidatePhenx.insert(candidatePhenxs[i].begin(), candidatePhenxs[i].end());
        candidatePhenxs[i].clear();
    }
    return candidatePhenx;
}

unsigned int getEndPhenx(temporalSequence &sequence, unsigned int lengthOfPhenx) {
    std::string sequenceString = std::to_string(sequence.seqID);
    return  std::stoul(sequenceString.substr(sequenceString.length() - lengthOfPhenx, lengthOfPhenx));
}

bool isPhenxOfInterest(unsigned int phenx, std::vector<unsigned int> phenxsOfInterest) {
    return std::find(phenxsOfInterest.begin(), phenxsOfInterest.end(),phenx) != phenxsOfInterest.end();
}

unsigned int getStartPhenx(temporalSequence &sequence, unsigned int lengthOfPhenx){
    //remove possible values stored sequenceID
    std::string sequenceString = std::to_string(sequence.seqID);
    return  std::stoul(sequenceString.substr(0,sequenceString.length() - lengthOfPhenx));
}

std::vector<size_t> getSequenceStartPositions(std::vector<temporalSequence> &sequences) {
    std::vector<size_t> startPositions;
    startPositions.emplace_back(0);
    for (int i = 1; i < sequences.size(); ++i) {
        if(sequences[i-1].seqID != sequences[i].seqID){
            startPositions.emplace_back(i);
        }
    }
    return startPositions;
}


std::vector<temporalSequence> extractSequencesWithSpecificStart(std::vector<temporalSequence> &originalSequences, std::uint64_t minDuration,
                                                                unsigned int bitShift, unsigned int lengthOfPhenx,
                                                                std::vector<unsigned int> &phenxOfInterest, int &numOfThreads){
    std::vector<temporalSequence> candidateSequences[numOfThreads];

    ips4o::parallel::sort(originalSequences.begin(),originalSequences.end(),
                          timedSequencesSorter, numOfThreads);
    std::vector<size_t>startPositions = getSequenceStartPositions(originalSequences);
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared(startPositions,originalSequences,candidateSequences, minDuration, lengthOfPhenx, phenxOfInterest )
    for (size_t i = 0; i < startPositions.size(); ++i) {
        size_t startPos = startPositions [i];
        size_t endPos;
        if(i < startPositions.size() - 1){
            endPos = startPositions[i+1];
        }else{
            endPos = originalSequences.size();
        }

        if(isPhenxOfInterest(getStartPhenx(originalSequences[startPos],lengthOfPhenx), phenxOfInterest)){
            unsigned int mininmalDuration = 0;
            for(size_t j = startPos; j < endPos; ++j ){
                temporalSequence seq = originalSequences[j];
                mininmalDuration = std::max(mininmalDuration,  seq.duration);
            }
            if (mininmalDuration >= minDuration) {
                candidateSequences[omp_get_thread_num()].insert(candidateSequences[omp_get_thread_num()].end(),
                                                                originalSequences.begin()+startPos,
                                                                originalSequences.begin()+endPos);
            }
        }
    }

    std::vector<temporalSequence>  mergedSequenceArrays;
    for (int i = 0; i < numOfThreads; ++i) {
        mergedSequenceArrays.insert(mergedSequenceArrays.end(),candidateSequences[i].begin(), candidateSequences[i].end());
        candidateSequences[i].clear();
    }
    return mergedSequenceArrays;
}

std::vector<temporalSequence> extractSequencesWithEnd(std::vector<temporalSequence> &originalSequences,
                             unsigned int bitShift, unsigned int lengthOfPhenx, std::set<unsigned  int> &allEndPhenx,
                             int &numOfThreads){
    std::vector<temporalSequence> candidateSequences[numOfThreads];
    std::vector<size_t> startPositions = getSequenceStartPositions(originalSequences);
    omp_set_num_threads(numOfThreads);
#pragma omp parallel for default (none) shared (startPositions, originalSequences, lengthOfPhenx, allEndPhenx, candidateSequences)
    for(size_t i = 0; i < startPositions.size(); ++i){
        size_t startPos = startPositions[i];
        unsigned int endPhenx = getEndPhenx(originalSequences[startPos], lengthOfPhenx);
        if(allEndPhenx.find(endPhenx) != allEndPhenx.end()){
            size_t endPos;
            if(i < startPositions.size() - 1){
                endPos = startPositions[i+1];
            }else{
                endPos = originalSequences.size();
            }
            for(size_t j = startPos; j < endPos; ++j) {
                candidateSequences[omp_get_thread_num()].emplace_back(originalSequences[j]);
            }
        }
    }
    std::vector<temporalSequence> allCandidateSequences;
    for(std::vector<temporalSequence> sequences : candidateSequences){
        allCandidateSequences.insert(allCandidateSequences.end(),sequences.begin(), sequences.end());
        sequences.clear();
        sequences.shrink_to_fit();
    }
    return allCandidateSequences;
}

unsigned int getCandidateBucket(unsigned int duration, std::vector<unsigned int> lowerBucketLimits) {
    unsigned int i = 0;
    while (i < lowerBucketLimits.size() && duration > lowerBucketLimits[i]){
        ++i;
    }
    return i;
}

