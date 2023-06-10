//
// Created by jonas on 11.04.23.
//

#ifndef TSPM_CPP_BACKEND_DBMARTENTRY_H
#define TSPM_CPP_BACKEND_DBMARTENTRY_H

#include <cstdint>
namespace tspm {
    struct dbMartEntry {
        int patID;
        int phenID;
        std::int64_t date;

    };
}//tspm
#endif //TSPM_CPP_BACKEND_DBMARTENTRY_H
