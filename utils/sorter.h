//
// Created by jonas on 11.04.23.
//
#include "dbMartEntry.h"
#include "temporalSequence.h"

#ifndef TSPM_CPP_BACKEND_SORTER_H
#define TSPM_CPP_BACKEND_SORTER_H
bool timedSequencesSorter(temporalSequence const& first, temporalSequence const& second);
bool dbMartSorter(const dbMartEntry first, const dbMartEntry second);
#endif //TSPM_CPP_BACKEND_SORTER_H
