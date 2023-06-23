//
// Created by Jonas on 11.04.23.
//

#ifndef TSPM_CPP_BACKEND_TEMPORALSEQUENCE_H
#define TSPM_CPP_BACKEND_TEMPORALSEQUENCE_H

#include <cstdint>
namespace tspm {
    struct temporalSequence {
        std::int64_t seqID;
        unsigned int duration;
        std::uint32_t patientID;
    }__attribute__((packed));
}//tspm
#endif //TSPM_CPP_BACKEND_TEMPORALSEQUENCE_H
