//
// Created by Bohlender,Ryan James on 9/10/19.
//

#ifndef CARVAIBD_ISGZIPPED_HPP
#define CARVAIBD_ISGZIPPED_HPP

#include <fstream>

template<class StringT>
class IsGzipped {
public:
  bool operator()(StringT file_path) {
	char byte1;
	char byte2;
	std::ifstream isource(file_path, std::ios_base::binary);
	if (check_file_exists(file_path)) {
	  isource.get(byte1);
	  isource.get(byte2);

	  return (byte1 == '\x1F' && byte2 == '\x8B');
	} else {
	  throw (std::logic_error("File doesn't exist."));

	}
  }
  bool check_file_exists(StringT file_path) {
	std::ifstream ifs(file_path);
	return ifs.good();
  }
};

#endif //CARVAIBD_ISGZIPPED_HPP
