

#include "utils.h"
extern size_t* startPositions;



long createSequence(std::string phenotypeA, std::string phenotypeB, int phenotypelenght = 6){
    //    fill up the second string with leading Zeros to avoid collisions with other sequences (creating unique seqs)
    //    e,g. phenA=1,phenB=12 ->seq=100012; without padding seq=112, could be phenA=11&phenB=2 or phenA=1&phenB=12
    phenotypeA.insert(phenotypeB.begin(),phenotypelenght-phenotypeA.size(), '0');
    phenotypeA.append(phenotypeB);
    return atol(phenotypeA.c_str());

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


std::vector<dbMartEntry> extractPatient(FILE* csv_file, int patId, int patIdColumn, int phenotypeIDColumn, int dateColumn) {
    if (csv_file == nullptr) {
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> lines;
    std::vector<dbMartEntry> dbMartEntries;
    char* line = nullptr;
    size_t len = 0;
    for(size_t i = startPositions[patId]; i < startPositions[patId + 1]; ++i){
        getline(&line, &len, csv_file);
        lines = getTokensFromLine(std::string(line));
        dbMartEntry entry;
        entry.patID = atoi(lines[patIdColumn].c_str());
        entry.phenID = atoi(lines[phenotypeIDColumn].c_str());
        entry.date = getTimeFromString(lines[dateColumn].c_str());
        dbMartEntries.push_back(entry);

    }
    return  dbMartEntries;
}


size_t countLinesInFile(const std::basic_string<char>& filename){
    FILE *file = fopen(filename.c_str(), "r");
    size_t line_count = 0;
    int patId = 0;
    char nextPatId[20];
    memset(nextPatId,0, sizeof (nextPatId));
    int posInLine = 0;
    bool detect_Patid = false;
    startPositions[0] = 1;
    if(file != nullptr) {
        while (!feof(file)) {
            char nextchar = fgetc(file);

            if (nextchar == '\n') {
                posInLine = 0;
                ++line_count;
                detect_Patid = true;
            } else if (detect_Patid) {
                if (nextchar != ',') {
                    nextPatId[posInLine] = nextchar;
                    ++posInLine;
                } else {
                    nextPatId[posInLine] = '\0';
                    int nextPat = atoi(nextPatId);
                    if (patId != nextPat) {
                        ++patId;
                        startPositions[patId] = line_count;

                    }
                    detect_Patid = false;
                    memset(nextPatId, 0, sizeof(nextPatId));
                }
            }
        }
    } else{
//        TODO: throw error
    }
    fclose(file);
//    #TODO check if file closed
    return line_count;
}

std::vector<std::string> getTokensFromLine(const std::string& line) {
    std::vector<std::string> vectorizedLine;
    size_t start;
    size_t end = 0;
    char delim = ',';

    while ((start = line.find_first_not_of(delim, end)) != std::string::npos) {
        end = line.find(delim, start);
        vectorizedLine.push_back(line.substr(start, end - start));
    }
    //remove \n from last entry
//    vectorizedLine[vectorizedLine.size()-1] = vectorizedLine.end()->substr(0,vectorizedLine.end()->length()-1);
    return vectorizedLine;
}

