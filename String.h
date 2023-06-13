//
// Created by 崔宇璐 on 2023/6/13.
//

#ifndef GUPOS_STRING_H
#define GUPOS_STRING_H

#include <string>

namespace Utility {
    void SplitString(const std::string& str, std::vector<std::string>& substr, const std::string& split);
    void RemoveSpace(std::string &str);
}

#endif //GUPOS_STRING_H
