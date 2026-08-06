/**
* Some tools for string
* @verbatim
	History
	 -1.0 jdhuang	 2019-04-07 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gstring.cpp
* @brief		for string
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#include "gutils/gstring.h"
#include "ExportHeaders/ExportBase.h"
#include <iostream>
#include <sstream>
#include <stdarg.h>
#include <math.h>
namespace great
{
	std::string trim(std::string s)
	{
		if (s.empty())
		{
			return s;
		}

		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);
		return s;
	}

	string str2upper(const string& str)
	{
		string tmp = str;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		return tmp;
	}

	std::string str_repaceall(const std::string& str, const std::string& substr_old, const std::string& substr_new)
	{
		string::size_type pos = 0;
		string result = str;
		while ((pos = result.find(substr_old)) != string::npos) {
			result.replace(pos, substr_old.length(), substr_new);
		}
		return result;
	}

	string str2lower(const string & str)
	{
		string tmp = str;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
		return tmp;
	}


}