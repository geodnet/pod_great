/**
*
* @verbatim
History
-1.0         Liu Chengbo              2020-09-18 Creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University.All rights reserved.
*
* @file		      gfantennathrust.cpp
* @brief		    The class for antenna thrust.
*
* @author          Liu Chengbo, Wuhan University
* @version		    1.0.0
* @date		    2020-09-18
*
*/

#include "gforce/gfantennathrust.h"
#include "gutils/ginfolog.h"

using namespace std;
namespace great
{
	/** constructor. */
	t_gfantennathrust::t_gfantennathrust(t_init_arg* arg, t_glog* log) : t_gforce(log)
	{
		string function_name = "t_gfantennathrust::t_gfantennathrust";
		if (!log)
		{
			write_log(log, 0, function_name,
				"ERROR : the log point is empty");
		}

		t_init_antennathrust* ant = NULL;
		ant = dynamic_cast<t_init_antennathrust*>(arg);

		if (!ant)
		{

			write_log(log, 4, function_name, "ERROR : init failed ");
			write_log(log, 4, function_name, "ERROR : no init argument ");
		}
		else
		{
			_force_model = ant->force_model;
			model = ant->model;
		}
	}
	/**
	* defalut deconstructor.
	*/
	t_gfantennathrust::~t_gfantennathrust()
	{
	}
	/**value satellite transmit power.
	* @param[in] block satellite block name
	*/
	bool t_gfantennathrust::getTransmitPower(string block)
	{
		string function_name = "t_gfantennathrust::getTransmitPower";
		try
		{
			const char* blockc = block.c_str();
			if (strcmp(blockc, "BLOCKIIF") == 0) { antennapower = 240.0; }
			else if (strcmp(blockc, "BLOCKIIIA") == 0) { antennapower = 300.0; }
			else if (strcmp(blockc, "BLOCKIIR-B") == 0 || strcmp(blockc, "BLOCKIIR-A") == 0) { antennapower = 60.0; }
			else if (strcmp(blockc, "BLOCKIIR-M") == 0) { antennapower = 145.0; }
			else if (strcmp(blockc, "GLONASS-M") == 0) { antennapower = 50.0; }
			else if (strcmp(blockc, "GLONASS-K1") == 0) { antennapower = 120.0; }
			else if (strcmp(blockc, "BEIDOU-2M") == 0) { antennapower = 130.0; }
			else if (strcmp(blockc, "BEIDOU-2I") == 0) { antennapower = 185.0; }
			else if (strcmp(blockc, "GALILEO-2") == 0) { antennapower = 265.0; }
			else if (strcmp(blockc, "BEIDOU-3M-CAST") == 0) { antennapower = 310.0; }
			else if (strcmp(blockc, "BEIDOU-3M-SECM") == 0) { antennapower = 280.0; }
			else
			{
				write_log(_glog, 4, function_name, "ERROR : invalid satellite type");
				return false;
			}
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while getting satellite transmit power");
			return false;
		}
	}
	/**calculate the accleration caused by antenna thrust.
	* @param[in] argument Parameters to be used in calculation
	*/
	bool t_gfantennathrust::calculate(t_gforce_calarg* argument)
	{
		string function_name = "t_gfantennathrust::calculate";
		try 
		{
			double c = 299792458.0;
			t_gfcalarg_antennathrust* arg = NULL;
			arg = dynamic_cast<t_gfcalarg_antennathrust*>(argument);
			if (!arg)
			{
				write_log(_glog, 4, function_name, "ERROR : no calculate argument");
				return false;
			}
			ColumnVector ysat = *(arg->sat_pos);
			ColumnVector force(0.0, 3);
			ColumnVector radvec(0.0, 3);
			_mass = (arg->sat_mass);
			blockname = arg->sat_blocknum;
			getTransmitPower(blockname);
			double abspos = sqrt(ysat[0] * ysat[0] + ysat[1] * ysat[1] + ysat[2] * ysat[2]);
			for (int num = 0; num < 3; num++)
			{
				radvec[num] = ysat[num] / abspos;
			}
			double ATforce = antennapower / c;
			for (int i = 0; i < 3; i++)
			{
				force[i] = force[i] + ATforce * radvec[i];
			}
			_acc = force / _mass;
			return true;
		}
		catch (...)
		{
			write_log(_glog, 4, function_name, "ERROR : throw exception while calculating antenna thrust force");
			return false;
		}
	}
}