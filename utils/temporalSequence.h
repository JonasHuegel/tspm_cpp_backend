//
// Created by Jonas on 11.04.23.
//

#ifndef TSPM_CPP_BACKEND_TEMPORALSEQUENCE_H
#define TSPM_CPP_BACKEND_TEMPORALSEQUENCE_H


struct temporalSequence{
    long seqID;
    unsigned int duration;
    unsigned int patientID;
}__attribute__((packed));

#endif //TSPM_CPP_BACKEND_TEMPORALSEQUENCE_H
