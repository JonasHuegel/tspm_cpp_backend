#include <fstream>
#include <algorithm>
#include <omp.h>
#include <mutex>
#include <map>
#include "utils/utils.h"

const size_t * startPositions;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

int main(int argc, char *argv[]) {
//  ======Input Variables=====
// TODO: change to real input variables
    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_fourtimes_processed.csv";
    std::string description = "tspm_test_";
    std::string outputDir = "/home/jonas/CLionProjects/tspm_cpp_backend/out/data/";
    int patIdColumn = 0;
    int phenotypeIdColumn = 1;
    int dateColumn = 3;
    int patientCount = 1167; //TODO add function to determine (extract from last line in file)
    bool createDuration = true;
    size_t patientPos[patientCount];
    startPositions = patientPos;
//   ======Initalise variables =====-
    std::cout << "indexing" << std::endl;
    size_t lineCount = countLinesInFile(fileName);

//   ====== Sequence Patients =======
    FILE *csvFilePointer = fopen(fileName.c_str(), "r");
    if (csvFilePointer == nullptr) {
        exit(EXIT_FAILURE);
    }
    char *header = nullptr;
    size_t len = 0;
    getline(&header, &len, csvFilePointer);
    std::cout << "extracting sequences" << std::endl;

    int patientId = 0;
    size_t sequenceNum = 0;
    omp_set_num_threads(16);
    std::mutex readMutex;
    std::mutex countMutex;
    int local_patID;
#pragma omp parallel for private(local_patID)
    for (size_t i = 0; i < patientCount; ++i) {

        std::vector<dbMartEntry> patientEntries;
        readMutex.lock();
        patientEntries = extractPatient(csvFilePointer, patientId, patIdColumn, phenotypeIdColumn, dateColumn);
        ++patientId;
        local_patID = patientId;
        local_patID--;
        readMutex.unlock();


        //the number sequences is equal to the gaussian sum, initalize vector with size to avoid re
        size_t numberOfSequences = (patientEntries.size() * (patientEntries.size() + 1)) / 2;
        std::vector<long> sequences;
        sequences.reserve(numberOfSequences);
        for (int j = 0; j < patientEntries.size() - 1; ++j) {
            for (int k = j + 1; k < patientEntries.size(); ++k) {
                sequences.push_back(createSequence(patientEntries[j].phenID, patientEntries[k].phenID));
                if (createDuration) {
                    sequences.back() = addDurationToSequence(sequences.back(), patientEntries[j].date,
                                                             patientEntries[k].date);
                }
            }
        }
        countMutex.lock();
        sequenceNum += sequences.size();
        std::cout << local_patID << " num of sequences: " << sequences.size() << std::endl;
        countMutex.unlock();


        std::string patIDString = std::to_string(local_patID);
        int patIDLength = 7;
        patIDString.insert(patIDString.begin(), patIDLength - patIDString.size(), '0');
        std::string patientFileName = std::string(outputDir).append(description).append(patIDString);
        writeSequencestoBinaryFile(patientFileName, sequences);
    }
    std::cout << sequenceNum << std::endl;
    std::cout << patientId << std::endl;
    std::cout << lineCount << std::endl;
    fclose(csvFilePointer);

    //TODO sparsity
//    use local hashmaps merge with openmp (map?)
//    - avoid extensive merging, large look ups, cache invalidations




    return 0;
}



