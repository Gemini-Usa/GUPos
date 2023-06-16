//
// Created by 崔宇璐 on 2023/6/13.
//

#include <regex>
#include "String.h"

void Utility::SplitString(const std::string& str, std::vector<std::string>& substr, const std::string& split) {
    std::string str0;
    for (auto c : str) {
        if (split.find(c) == std::string::npos) str0.push_back(c);
        else { substr.push_back(str0); str0.clear(); }
    }
    substr.push_back(str0);
    str0.clear();
}

void Utility::RemoveSpace(std::string &str) {
    int start = 0, end = str.size() - 1;
    while (start <= end && str[start] == ' ') start++;
    while (end >= start && str[end] == ' ') end--;
    str = str.substr(start, end - start + 1);
    std::string result = "";
    int i = 0;
    while (i < str.size()) {
        if (str[i] != ' ') {
            result += str[i];
            i++;
        } else {
            result += ' ';
            while (str[i] == ' ') i++;
        }
    }
    str = result;
}

void Utility::SplitValueAndUnit(const std::string &input, double &value, std::string &unit) {
    std::regex pattern(R"((\d+\.?\d*(?:[eE][-+]?\d+)?)(\D*))");
    std::smatch match;
    if (std::regex_search(input, match, pattern)) {
        value = std::stod(match[1]);
        unit = match[2];
    }
}
