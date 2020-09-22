//
// Created by Bohlender,Ryan James on 9/22/20.
//

#ifndef CARVAIBD_INFO_HPP
#define CARVAIBD_INFO_HPP

#include <string>
#include <vector>
#include <map>

class Info {
  std::map<std::string, std::vector<double>> data;
  std::map<std::string, int> field_map;
public:
  explicit Info(const std::string &fpath);

  double get_field(const std::string &segment, const std::string &field);
};

#endif //CARVAIBD_INFO_HPP
