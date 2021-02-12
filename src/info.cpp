//
// Created by Bohlender,Ryan James on 9/22/20.
//

#include <fstream>
#include "info.hpp"
#include "split.hpp"
#include "inputvalidator.hpp"

Info::Info(std::istream &ifs, Parameters &params) {
  std::string line;
  long lineno = -1;
  InputValidator iv(true);
  while (std::getline(ifs, line)) {
	lineno++;
	iv.check_info(line, lineno, bool(params.AF), bool(params.cM));
	RJBUtil::Splitter<std::string> splitter(line, " \t");
	if (lineno == 0) { // Handle header
	  // First two fields should be fixed as chr and segID
	  if (splitter[0] != "chr") {
		throw (std::runtime_error("ERROR: Malformed info file. First field should be chr."));
	  }
	  if (splitter[1] != "clusterID") {
		throw (std::runtime_error("ERROR: Malformed info file. Second field should be clusterID."));
	  }
	  for (int i = 2; i < splitter.size(); i++) {
		field_map[splitter[i]] = i - 2;
	  }
	  continue;
	}
	data[splitter[1]] = {};
	for (int i = 2; i < splitter.size(); i++) {
	  try {
		data[splitter[1]].push_back(std::stod(splitter[i]));
	  } catch (std::exception &e) {
		throw (std::runtime_error("ERROR: Error while reading info file. Failed to convert item to double."));
	  }
	}
  }
  if (ifs.bad()) {
	throw (std::runtime_error("ERROR: Error while reading info file."));
  }
}

double Info::get_field(const std::string &segment, const std::string &field) {
  return data[segment][field_map[field]];
}

bool Info::filter_segment(const std::string &segment, const Parameters &params) {
  // If neither filter is used the segment will automatically pass.
  // If only a single is used, passing will depend on the single filter.
  bool AF = true;
  bool cM = true;
  double af_val = get_field(segment, "freq");
  if (params.AF) {
	AF = get_field(segment, "freq") < *params.AF;
  }
  if (params.cM) {
	cM = get_field(segment, "cM") >= *params.cM;
  }
  return !(AF & cM); // If both filters are used the segment must pass both filters.
  // We negate because the logic is, true = filter, false = keep
}
