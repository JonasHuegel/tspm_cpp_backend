//
// Created by jonas on 11.04.23.
//
#include "sorter.h"
namespace tspm {
/**
 * Comparator for temporalSequence struct. Added all greater as and smaller as cases for all attributes. In the equal
 * cases the next attribute should be used. In the case that all attributes are equal return false. The sequenceID is the
 * most important parameter, when to calculate the min/max values and the buckets for a unique sequence ID. The other
 * attributes are added for completeness.
 * @param first
 * @param second
 * @return
 */
    bool timedSequencesSorter(temporalSequence const &first, temporalSequence const &second) {
        if (first.seqID < second.seqID)
            return true;
        if (first.seqID > second.seqID)
            return false;
        if (first.duration < second.duration)
            return true;
        if (first.duration > second.duration)
            return false;
        if (first.patientID < second.patientID)
            return true;
        if (first.patientID > second.patientID)
            return false;

        return false;

    }

    bool timedSequenceByPatientIDSorter(temporalSequence const &first, temporalSequence const &second) {
        if (first.patientID < second.patientID)
            return true;
        if (first.patientID > second.patientID)
            return false;
        if (first.seqID < second.seqID)
            return true;
        if (first.seqID > second.seqID)
            return false;
        if (first.duration < second.duration)
            return true;
        if (first.duration > second.duration)
            return false;
        return false;

    }


    bool dbMartSorter(const dbMartEntry &first, const dbMartEntry &second) {
        if (first.patID < second.patID)
            return true;
        if (first.patID > second.patID)
            return false;
        if (first.date < second.date)
            return true;
        if (first.date > second.date)
            return false;
        if (first.phenID < second.phenID)
            return true;
        if (first.phenID > second.phenID)
            return false;

        return false;

    }
}//tspm