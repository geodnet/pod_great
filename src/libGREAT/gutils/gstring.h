/**
* Some tools for string
* @verbatim
	History
	 -1.0 jdhuang	 2019-04-07 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gstring.h
* @brief		for string
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#include <string>
#include <vector>
#include <set>
#include "ExportHeaders/ExportBase.h"
using namespace std;

namespace great
{
	BASE_LIBRARY_EXPORT std::string trim(std::string s);

	BASE_LIBRARY_EXPORT void split(const std::string& s, std::string delim, std::vector< std::string >& ret);
	

	BASE_LIBRARY_EXPORT std::string str2upper(const std::string& str);
	BASE_LIBRARY_EXPORT std::string str_repaceall(const std::string& str,const std::string& substr_old, const std::string& substr_new);
	BASE_LIBRARY_EXPORT std::string str2lower(const std::string& str);

#include <stdarg.h>

	inline BASE_LIBRARY_EXPORT std::string format(const char* fmt, ...)
	{
		va_list ap;
		va_start(ap, fmt);
		int len = vsnprintf(nullptr, 0, fmt, ap);
		va_end(ap);
		std::string buf(len + 1, '\0');
		va_start(ap, fmt);
		vsnprintf(&buf[0], buf.size(), fmt, ap);
		va_end(ap);
		buf.pop_back();
		return buf;
	}



}
