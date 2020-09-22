//
// Created by Bohlender,Ryan James on 9/22/20.
//

#include <fstream>
#include "info.hpp"
#include "split.hpp"

Info::Info(const std::string &fpath) {
  std::string line;
  int lineno = -1;
  std::ifstream ifs(fpath);
  if(!ifs.is_open()) {
    throw(std::runtime_error("ERROR: Failed to open info file."));
  }
  while(std::getline(ifs, line)) {
    lineno++;
	RJBUtil::Splitter<std::string> splitter(line, " \t");
	if (lineno == 0) { // Handle header
	  // First two fields should be fixed as chr and segID
	  if(splitter[0] != "chr") {
	    throw(std::runtime_error("ERROR: Malformed info file. First field should be chr."));
	  }
	  if(splitter[1] != "segID") {
		throw(std::runtime_error("ERROR: Malformed info file. Second field should be segID."));
	  }
	  for(int i = 2; i < splitter.size(); i++) {
	    field_map[splitter[i]] = i - 2;
	  }
	  continue;
	}
	data[splitter[1]] = {};
	for(int i = 2; i < splitter.size(); i++) {
	  try {
		data[splitter[i]].push_back(std::stod(splitter[i]));
	  } catch(std::exception &e) {
		throw(std::runtime_error("ERROR: Error while reading info file. Failed to convert item to double."));
	  }
	}
  }
  if(ifs.bad()) {
    throw(std::runtime_error("ERROR: Error while reading info file."));
  }
}

double Info::get_field(const std::string &segment, const std::string &field) {
  return data[segment][field_map[field]];
}
