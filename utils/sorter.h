/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* sort.h
***********************************************************************************
* MIT License
*
* Copyright (c) 2023 Jonas Hügel
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
**********************************************************************************/
#include "dbMartEntry.h"
#include "temporalSequence.h"

#ifndef TSPM_CPP_BACKEND_SORTER_H
#define TSPM_CPP_BACKEND_SORTER_H
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
    bool timedSequencesSorter(temporalSequence const &first, temporalSequence const &second);
    
    /**
     * Comparator for temporalSequence struct. Using the sequence ID as first, and the patient id as second attribute. The sequence duration is the third parameter. Added all greater as and smaller as cases for all attributes. In the equal
     * cases the next attribute should be used. In the case that all attributes are equal return false.
     * @param first
     * @param second
     * @return
     */
    bool timedSequenceByPatientIDSorter(temporalSequence const &first, temporalSequence const &second);
    
    /**
     * Comperator function for temporal sequence structs. It sorts the sequences  by the patient id as the main attribute
    */
    bool timedSequenceByPatientIDSorter(temporalSequence const &first, temporalSequence const &second);

    /**
     * Comperator function to compare the dbMart structs. Sorts by patientID, start_date and phenx.
    */
    bool dbMartSorter(const dbMartEntry &first, const dbMartEntry &second);
}//tspm
#endif //TSPM_CPP_BACKEND_SORTER_H
