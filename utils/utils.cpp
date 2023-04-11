#include "utils.h"



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
    //TODO define end for cases when a file skips a patient number
    //TODO set end for last patID == largest patID
    for(size_t i = (*startPositions)[patId]; i < (*startPositions)[patId + 1]; ++i){
        fgets(line, len, csv_file);
        lines = getTokensFromLine(std::string(line));
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

//size_t countLinesInFile(const std::basic_string<char>& filename){
//    FILE *file = fopen(filename.c_str(), "r");
//    size_t line_count = 0;
//    int patId = 0;
//    char nextPatId[20];
//    memset(nextPatId,0, sizeof (nextPatId));
//    int posInLine = 0;
//    bool detect_Patid = false;
//    startPositions[0] = 1;
//    if(file != nullptr) {
//        while (!feof(file)) {
//            char nextchar = fgetc(file);
//
//            if (nextchar == '\n') {
//                posInLine = 0;
//                ++line_count;
//                detect_Patid = true;
//            } else if (detect_Patid) {
//                if (nextchar != ',') {
//                    nextPatId[posInLine] = nextchar;
//                    ++posInLine;
//                } else {
//                    nextPatId[posInLine] = '\0';
//                    int nextPat = atoi(nextPatId);
//                    if (patId != nextPat) {
//                        ++patId;
//                        startPositions[patId] = line_count;
//
//                    }
//                    detect_Patid = false;
//                    memset(nextPatId, 0, sizeof(nextPatId));
//                }
//            }
//        }
//    } else{
////        TODO: throw error
//    }
//    fclose(file);
////    #TODO check if file closed
//    return line_count;
//}

//long addDurationToSequence(long &back, long startDate, long endDate) {
//    long milliSecondsPerDay = 1000 * 60 * 60 * 24;
//    long duration = (endDate-startDate)%milliSecondsPerDay;
//    int bitsForDuration = 40;
//    back = (duration << bitsForDuration) | back;
//    return back;
//}


unsigned int getDuration(long startDate, long endDate) {
    unsigned int milliSecondsPerDay = 1000 * 60 * 60 * 24;
    unsigned int duration = (endDate - startDate) / milliSecondsPerDay;
    return duration;
}

std::vector<std::string> getTokensFromLine(const std::string& line) {
    std::vector<std::string> vectorizedLine;
    size_t start;
    size_t end = 0;
    char delim = ',';

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
//        std::cout << "open file: " << patientFileName << std::endl;

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

