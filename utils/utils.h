
#ifndef TSPM_UTILS_H
#define TSPM_UTILS_H
#include "ips4o.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstring>
#include <map>
#include <omp.h>


struct dbMartEntry{
    int patID;
    int phenID;
    long date;

};

std::vector<dbMartEntry> extractPatient(FILE* csv_file, int patId, int patIdColumn = 0, int phenotypeIDColumn = 1, int dateColumn = 3);
size_t countLinesInFile(const std::basic_string<char>& filename);
std::vector<std::string> getTokensFromLine(const std::string& line);
long createSequence(int phenotypeA, int phenotypeB, int phenotypelenght = 6);
size_t writeSequencesToBinaryFile(std::string patientFilename, std::vector<long> sequences);
long addDurationToSequence(long &back, long startDate, long endDate);
std::map<long, size_t> summarizeSequences(int numberOfPatients, bool storesDuration, const std::string& outputDir, const std::string& file_prefix);



std::pair<size_t,size_t>  countLinesAndPatientsInFile(const std::basic_string<char>& filename, char delimiter);
long getFileSize(const std::string& filename);
unsigned int getDuration(long startDate, long endDate);

#endif //TSPM_UTILS_H
