/**
*
* @verbatim
History
-1.0 glfeng BoWong xjhan 2019-02-26  creat the file.
-1.1 BoWong              2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambcommon.cpp
* @brief	declare some classes and some mathematical method.
*
* @author   glfeng BoWong xjhan, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/
#include "gambfix/gambcommon.h"
using namespace std;

namespace great
{


	void getMeanWgtVec(vector<pair<double, double>>& wx, double& mean, double& sigma, double& mean_sig,int* EraseNum)
	{
		if (EraseNum) {
			*EraseNum = 0;
		}
		if (wx.empty())
		{
			mean_sig = 0.0;
			return;
		}
		double wgt_sum;
		int erase_num = 0;
		bool loop_flag = true;
		while (loop_flag)
		{
			if (wx.size() - erase_num == 1)
			{
				sigma = 999.0;
				mean_sig = 999.0;
				return;
			}

			// mean
			mean = 0.0;
			wgt_sum = 0.0;
			for (auto itwx = wx.begin(); itwx != wx.end(); itwx++)
			{
				mean += itwx->first * itwx->second;
				wgt_sum += itwx->second;
			}
			mean /= wgt_sum;

			// sigma
			sigma = 0.0;
			wgt_sum = 0.0;
			for (auto itwx = wx.begin(); itwx != wx.end(); itwx++)
			{
				sigma += (itwx->first - mean) * (itwx->first - mean) * itwx->second;
				wgt_sum += itwx->second;
			}

			sigma = sqrt(sigma / (wx.size() - erase_num - 1.0));
				mean_sig = sigma / sqrt(wgt_sum);

			loop_flag = false;
			for (auto itwx = wx.begin(); itwx != wx.end();itwx++)
			{
				if (itwx->second >0 && fabs(itwx->first - mean) > 3 * (sigma))
				{
					itwx->second = 0;
					erase_num += 1;
					if (EraseNum) {
						*EraseNum += 1;
					}
					//itwx = wx.erase(itwx);
					loop_flag = true;
				}
			}
		}
		return;
	}

	

}