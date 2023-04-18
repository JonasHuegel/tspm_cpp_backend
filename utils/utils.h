
#ifndef TSPM_UTILS_H
#define TSPM_UTILS_H
#include "../lib/ips4o/ips4o.hpp"
#include "dbMartEntry.h"
#include "temporalSequence.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstring>
#include <map>
#include <omp.h>
#include <mutex>
#include <set>
#include <filesystem>


std::vector<dbMartEntry>
extractPatient(FILE *csv_file, std::vector<size_t> *startPositions, int patId, int patIdColumn = 0, int phenotypeIDColumn = 1,
               int dateColumn = 3);

std::vector<std::string> getTokensFromLine(const std::string &line, char delim);

long createSequence(int phenotypeA, int phenotypeB, int phenotypelenght = 6);

size_t writeSequencesToBinaryFile(std::string patientFilename, std::vector<long> sequences);

long addDurationToSequence(long &back, long startDate, long endDate);

std::map<long, size_t> summarizeSequences(int numberOfPatients, bool storesDuration, const std::string& outputDir, const std::string& file_prefix);

std::pair<size_t, size_t>
determinePatientStartPositionsInFile(const std::basic_string<char> &filename, char delimiter, std::vector<size_t> *startPositions);

long getFileSize(const std::string& filename);

unsigned int getDuration(long startDate, long endDate);

long getTimeFromString(const char * date_string);

std::vector<dbMartEntry> extractDBMartFromCsv(FILE *csv_file, int patIdColumn, int phenotypeIDColumn,
                                              int dateColumn, char delim = ',');

long writeSequencesAsCsV(std::string fileName, std::string filepath, char delimiter, size_t numOfSequences,
                         temporalSequence * temporalSequences, bool debug =false);

#endif //TSPM_UTILS_H
