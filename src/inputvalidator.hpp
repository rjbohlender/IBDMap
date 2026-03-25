//
// Created by Bohlender,Ryan James on 2/11/21.
//

#ifndef CARVAIBD_INPUTVALIDATOR_HPP
#define CARVAIBD_INPUTVALIDATOR_HPP

#include <string>

bool exists(const std::string &path);

/**
 * @brief Validator for input files to ensure properly formatted input.
 */
class InputValidator {
    size_t pheno_column_count = 0;

public:
    static void check_gmap(const std::string &line, size_t line_no);
    void check_pheno(const std::string &line, size_t line_no);
};

#endif//CARVAIBD_INPUTVALIDATOR_HPP
