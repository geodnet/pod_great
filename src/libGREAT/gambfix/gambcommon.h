/**
*
* @verbatim
History
-1.0 glfeng BoWong xjhan 2019-02-26  creat the file.
-1.1 BoWong              2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambcommon.h
* @brief	declare some classes and some mathematical method.
*
* @author   glfeng BoWong xjhan, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/
#ifndef GAMBCOMMON_H
#define GAMBCOMMON_H 

#include "ExportHeaders/ExportAMB.h"
#include <string>
#include <map>
#include <memory>
#include <tuple>
#include <list>
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"

using namespace std;
using namespace gnut;

namespace great
{

	AMB_LIBRARY_EXPORT void getMeanWgtVec(vector<pair<double, double>>& wx, double& mean, double& sigma, double& mean_sig, int* EraseNum = nullptr);

}

#endif
