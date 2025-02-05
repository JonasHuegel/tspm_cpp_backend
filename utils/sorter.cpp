/**********************************************************************************
* transitive Sequence Pattern Mining Plus algorithm
* sorter.cpp
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
#include "sorter.h"
namespace tspm {
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

    bool timedSequencesBySeqAndIDSorter(temporalSequence const &first, temporalSequence const &second) {
      if (first.seqID < second.seqID)
        return true;
      if (first.seqID > second.seqID)
        return false;
      if (first.patientID < second.patientID)
        return true;
      if (first.patientID > second.patientID)
        return false;
      if (first.duration < second.duration)
        return true;
      if (first.duration > second.duration)
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