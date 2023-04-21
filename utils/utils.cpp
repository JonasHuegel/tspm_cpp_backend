#include "utils.h"


std::vector<dbMartEntry> extractDBMartFromCsv(FILE *csv_file, int patIdColumn, int phenotypeIDColumn,
                                              int dateColumn, char delim){
    if (csv_file == nullptr) {
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> lines;
    std::vector<dbMartEntry> dbMartEntries;
    char line[2048];
    size_t len = 2048;
    while(fgets(line, len, csv_file) != NULL){
        lines = getTokensFromLine(std::string(line), delim);
        dbMartEntry entry;
        entry.patID = atoi(lines[patIdColumn].c_str());
        entry.phenID = atoi(lines[phenotypeIDColumn].c_str());
        entry.date = getTimeFromString(lines[dateColumn].c_str());
        dbMartEntries.emplace_back(entry);
    }
    return  dbMartEntries;


}


size_t writeSequencesToBinaryFile(std::string patientFilename, std::vector<long> sequences){
    FILE* patientFile;
    patientFile = fopen(patientFilename.c_str(), "wb");
    size_t written;
    if(patientFile!= nullptr) {
       written = std::fwrite(&sequences[0], 1, sequences.size() * sizeof(long), patientFile);
    }else{
        exit(EXIT_FAILURE);
    }
    fclose(patientFile);
    return written;
}



long createSequence(int phenotypeA, int phenotypeB, int phenotypelenght){
    std::string phenA = std::to_string(phenotypeA);
    std::string phenB = std::to_string(phenotypeB);
    phenB.insert(phenB.begin(),phenotypelenght-phenB.size(), '0');
    phenA.append(phenB);
    return atol(phenA.c_str());
}

long getTimeFromString(const char * date_string) {
    std::istringstream stringStream(date_string);
    std::tm time = {};
    stringStream >> std::get_time(&time,"%Y-%m-%dT%H:%M:%SZ");

    std::time_t  time_stamp = std::mktime(&time);
    return time_stamp;
}


std::vector<dbMartEntry>
extractPatient(FILE *csv_file, std::vector<size_t> *startPositions, int patId, int patIdColumn, int phenotypeIDColumn,
               int dateColumn) {
    if (csv_file == nullptr) {
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> lines;
    std::vector<dbMartEntry> dbMartEntries;
    char line[2048];
    size_t len = 2048;
    for(size_t i = (*startPositions)[patId]; i < (*startPositions)[patId + 1]; ++i){
        if(fgets(line, len, csv_file) == NULL){
            return dbMartEntries;
        }
        lines = getTokensFromLine(std::string(line), ',');
        dbMartEntry entry;
        entry.patID = atoi(lines[patIdColumn].c_str());
        entry.phenID = atoi(lines[phenotypeIDColumn].c_str());
        entry.date = getTimeFromString(lines[dateColumn].c_str());
        dbMartEntries.emplace_back(entry);

    }
    return  dbMartEntries;
}

std::pair<size_t, size_t>
determinePatientStartPositionsInFile(const std::basic_string<char> &filename, char delimiter, std::vector<size_t> *startPositions) {
    FILE *file = fopen(filename.c_str(), "r");
    size_t line_count = 0; //Headerline was skipped
    size_t pat_count = 0;
    int patId = 0;
    char nextPatId[20];
    memset(nextPatId,0, sizeof (nextPatId));
    bool detect_Patid = false;
    std::vector<char> patIDAsString;
    startPositions->emplace_back(1);
    if(file != nullptr) {
        while (!feof(file)) {
            char nextchar = fgetc(file);
            if (nextchar == '\n') {
                ++line_count;
                detect_Patid = true;
            } else if (detect_Patid) {
                if (nextchar != delimiter) {
                    patIDAsString.emplace_back(nextchar);

                } else {
                    patIDAsString.emplace_back('\0');
                    int nextPat = atoi(patIDAsString.data());
                    if (patId != nextPat) {
                        ++pat_count;
                        startPositions->emplace_back(line_count);
                        patId=nextPat;
                    }
                    detect_Patid = false;
                    patIDAsString.clear();
                }
            }
        }
    } else{
        std::cout << "could not open file: " << filename <<std::endl;
        EXIT_FAILURE;
    }
    fclose(file);
//    #TODO check if file closed
startPositions->emplace_back(line_count); // add end
    ++pat_count;//increase the patcount by one for the last patient (no increcement for this patient in the loop!)
    return std::make_pair(line_count,pat_count);
}


unsigned int getDuration(long startDate, long endDate) {
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

long getFileSize(const std::string& filename){
    try {
        return std::filesystem::file_size(filename);
    } catch(std::filesystem::filesystem_error& e) {
        std::cout << e.what() << '\n';
        return 0;
    }
}

std::map<long, size_t> summarizeSequences(int numberOfPatients, bool storesDuration, const std::string& outputDir, const std::string& file_prefix) {

    std::map<long, size_t> globalSequenceMap;
    const int numOfProcs = omp_get_max_threads();
    std::map<long, size_t> localmaps[numOfProcs];
    std::vector<int> patientsPerThread(numOfProcs);
    std::mutex map_mutex;
#pragma omp parallel for default (none) shared(numberOfPatients, file_prefix, outputDir,storesDuration,localmaps,patientsPerThread,map_mutex, globalSequenceMap, std::cout)
    for (size_t i = 0; i < numberOfPatients; ++i) {
        std::string patIDString = std::to_string(i);
        int patIDLength = 7;
        patIDString.insert(patIDString.begin(), patIDLength - patIDString.size(), '0');
        std::string patientFileName = std::string(outputDir).append(file_prefix).append(patIDString);
        long numberOfSequences = getFileSize(patientFileName)/sizeof(long);
        FILE* patientFile = fopen(patientFileName.c_str(),"rb");

        std::set<long> patientSequenceSet;
        for (int j = 0; j < numberOfSequences; ++j) {
            long sequence =0;
            fread(&sequence, sizeof(long ),1,patientFile);
            if(storesDuration){
                sequence = (sequence<<24)>>24;
            }
            if(patientSequenceSet.insert(sequence).second){
                ++localmaps[omp_get_thread_num()][sequence];
            }

        }

        fclose(patientFile);
        ++patientsPerThread[omp_get_thread_num()];
        if(patientsPerThread[omp_get_thread_num()]>=50) {
            map_mutex.lock();
            for (std::pair<long, size_t> entry: localmaps[omp_get_thread_num()]) {
                globalSequenceMap[entry.first] += entry.second;
            }
            map_mutex.unlock();
            localmaps[omp_get_thread_num()] = std::map<long, size_t>();
            patientsPerThread[omp_get_thread_num()] = 0;
        }
    }

    //merge all local maps into final map
    for (int i = 0; i < numOfProcs; ++i) {

        for (std::pair<long, size_t> entry: localmaps[i]) {
            globalSequenceMap[entry.first] += entry.second;
        }
    }
    return globalSequenceMap;
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
