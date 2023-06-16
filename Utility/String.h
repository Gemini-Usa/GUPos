//
// Created by GeminiUsa on 2023/6/13.
//

#ifndef GUPOS_STRING_H
#define GUPOS_STRING_H

#include <string>
#include <vector>

namespace Utility {
    void SplitString(const std::string& str, std::vector<std::string>& substr, const std::string& split);
    void RemoveSpace(std::string &str);
    void SplitValueAndUnit(const std::string& input, double& value, std::string& unit);
}

#endif //GUPOS_STRING_H
