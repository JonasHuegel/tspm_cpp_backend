
#ifndef TSPM_UTILS_H
#define TSPM_UTILS_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstring>

struct dbMartEntry{
    int patID;
    int phenID;
    long date;

};
std::vector<dbMartEntry> extractPatient(FILE* csv_file, int patId, int patIdColumn = 0, int phenotypeIDColumn = 1, int dateColumn = 3);
size_t countLinesInFile(const std::basic_string<char>& filename);
std::vector<std::string> getTokensFromLine(const std::string& line);
long createSequence(int phenotypeA, int phenotypeB, int phenotypelenght = 6);
size_t writeSequencestoBinaryFile(std::string patientFilename, std::vector<long> sequences);
long addDurationToSequence(long &back, long startDate, long endDate);



#endif //TSPM_UTILS_H
