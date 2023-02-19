#include <fstream>
#include <algorithm>
#include <omp.h>
#include <mutex>
#include "utils/utils.h"

const size_t * startPositions;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
int main(int argc, char * argv[]) {
//  ======Input Variables=====
// TODO: change to real input variables
    std::string fileName = "/home/jonas/CLionProjects/tspm_cpp_backend/data/dbmart_processed.csv";
    int patIdColumn = 0;
    int phenotypeIdColumn = 1;
    int dateColumn = 3;
    int patientCount = 1167; //TODO add function to determine (extract from last line in file)
    size_t patientPos[patientCount];
    startPositions = patientPos;
//   ======Initalise variables =====-
    std::cout  << "indexing" << std::endl;
    size_t lineCount = countLinesInFile(fileName);


//   ====== Sequence Patients =======

      // TODO rename the function
    FILE* csvFilePointer = fopen(fileName.c_str(), "r");
    if (csvFilePointer == nullptr) {
        exit(EXIT_FAILURE);
    }
    char* header = nullptr;
    size_t len = 0;
    getline(&header, &len, csvFilePointer);
    std::cout << "extracting sequences" << std::endl;

    int patientId =0;
    size_t sequenceNum=0;
    omp_set_num_threads(16);
    std::mutex readMutex;
    std::mutex countMutex;
    #pragma omp parallel for
    for(size_t i = 0; i < patientCount; ++i) {

        std::vector<dbMartEntry> patientEntries;
        readMutex.lock();
        patientEntries = extractPatient(csvFilePointer, patientId);
        ++patientId;
        int local_patID = patientId;
        readMutex.unlock();


         //TODO create sequences
         //the number sequences is equal to the gaussian sum
//         size_t numberOfSequences = (patientEntries.size() * (patientEntries.size() + 1)) / 2;
         std::vector<long> sequences;//(numberOfSequences);
         for(int j = 0; j < patientEntries.size() - 1; ++j){
             for(int k = j+1; k < patientEntries.size(); ++k){
                 sequences.push_back(createSequence(patientEntries[j].phenID, patientEntries[k].phenID));
             }
         }
         countMutex.lock();
         sequenceNum+= sequences.size();
         std::cout << local_patID << " num of sequences: " <<sequences.size() <<std::endl;
         countMutex.unlock();

         //TODO write sequences
     }
    std::cout <<sequenceNum << std::endl;
    std::cout << patientId << std::endl;
    std::cout << lineCount <<std::endl;
    //TODO sparsity


    return 0;
}
#pragma clang diagnostic pop

