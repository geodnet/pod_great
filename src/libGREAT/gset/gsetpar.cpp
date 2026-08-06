/**
*
* The format :
* @verbatim
	History
	 -1.0 jdhuang	 2019-04-07 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetpar.cpp
* @brief		set pars from XML
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#include "gset/gsetpar.h"

//#define DEBUG
using namespace gnut;

great::t_gsetpar::t_gsetpar() : t_gsetbase()
{
	_set.insert(XMLKEY_PARS);
}

great::t_gsetpar::~t_gsetpar()
{
}

double great::t_gsetpar::sigCX()
{
	//get rkf model node
	string sigCX = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_GEO).attribute("sigCX").value();



	if (sigCX.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for the sigCX, use the defalut number : " + dbl2str(_default_sigCX));
		return _default_sigCX;
	}
	else
	{
		return str2dbl(sigCX);
	}
}

double great::t_gsetpar::sigCY()
{
	//get rkf model node
	string sigCY = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_GEO).attribute("sigCY").value();



	if (sigCY.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigCY, use the defalut number : " + dbl2str(_default_sigCY));
		return _default_sigCY;
	}
	else
	{
		return str2dbl(sigCY);
	}
}

double great::t_gsetpar::sigCZ()
{
	//get rkf model node
	string sigCZ = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_GEO).attribute("sigCZ").value();



	if (sigCZ.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigCZ, use the defalut number : " + dbl2str(_default_sigCZ));
		return _default_sigCZ;
	}
	else
	{
		return str2dbl(sigCZ);
	}
}

double great::t_gsetpar::sigXpole()
{
	//get rkf model node
	string sigXPOLE = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_ERP).attribute("sigXPOLE").value();


	if (sigXPOLE.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigXPOLE, use the defalut number : " + dbl2str(_default_sigXPOLE));
		return _default_sigXPOLE;
	}
	else
	{
		return str2dbl(sigXPOLE);
	}
}

double great::t_gsetpar::sigYpole()
{
	//get rkf model node
	string sigYPOLE = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_ERP).attribute("sigYPOLE").value();


	if (sigYPOLE.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigYPOLE, use the defalut number : " + dbl2str(_default_sigYPOLE));
		return _default_sigYPOLE;
	}
	else
	{
		return str2dbl(sigYPOLE);
	}
}

double great::t_gsetpar::sigDxpole()
{
	//get rkf model node
	string sigDXPOLE = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_ERP).attribute("sigDXPOLE").value();



	if (sigDXPOLE.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigDXPOLE, use the defalut number : " + dbl2str(_default_sigDXPOLE));
		return _default_sigDXPOLE;
	}
	else
	{
		return str2dbl(sigDXPOLE);
	}
}

double great::t_gsetpar::sigDypole()
{
	//get rkf model node
	string sigDYPOLE = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_ERP).attribute("sigDYPOLE").value();


	if (sigDYPOLE.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigDYPOLE, use the defalut number : " + dbl2str(_default_sigDXPOLE));
		return _default_sigDXPOLE;
	}
	else
	{
		return str2dbl(sigDYPOLE);
	}
}

double great::t_gsetpar::sigUt1()
{
	//get rkf model node
	string sigUT1 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_ERP).attribute("sigUT1").value();


	if (sigUT1.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigUT1, use the defalut number : " + dbl2str(_default_sigUT1));
		return _default_sigUT1;
	}
	else
	{
		return str2dbl(sigUT1);
	}
}

double great::t_gsetpar::sigDut1()
{
	//get rkf model node
	string sigDUT1 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).child(XMLKEY_PARS_ERP).attribute("sigDUT1").value();



	if (sigDUT1.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigDUT1, use the defalut number : " + dbl2str(_default_sigDUT1));
		return _default_sigDUT1;
	}
	else
	{
		return str2dbl(sigDUT1);
	}
}

double great::t_gsetpar::sigZtd(string sta_name)
{
	//get rkf model node
	xml_node station = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute("ID", sta_name.c_str());
	string   sigZTD = station.attribute("sigZTD").value();



	if (sigZTD.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigZTD, use the defalut number : " + dbl2str(_default_sigZTD));
		return _default_sigZTD;
	}
	else
	{
		return str2dbl(sigZTD);
	}
}


double great::t_gsetpar::sigTropPd(string sta_name)
{
	//get rkf model node
	xml_node station = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute("ID", sta_name.c_str());
	string   sigTropPd = station.attribute("sigTropPd").value();



	if (sigTropPd.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigTropPd, use the defalut number : " + dbl2str(_default_sigTropPd));
		return _default_sigTropPd;
	}
	else
	{
		return str2dbl(sigTropPd);
	}
}


double great::t_gsetpar::sigRecCLK(string sta_name)
{
	//get rkf model node
	xml_node station = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute("ID", sta_name.c_str());
	string sigCLK = station.attribute("sigCLK").value();


	if (sigCLK.empty())
	{
		if (_log) _log->comment(1, "WARNNING : no setting for sigCLK, use the defalut number : " + dbl2str(_default_sta_sigCLK));
		return _default_sta_sigCLK;
	}
	else
	{
		return str2dbl(sigCLK);
	}
}


double great::t_gsetpar::sigSatCLK(string sat_name)
{
	//get rkf model node
	xml_node satellite = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute("ID", sat_name.c_str());
	string   sigCLK = satellite.attribute("sigCLK").value();



	if (sigCLK.empty())
	{
		xml_node common_node = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute("ID", "XXX");
		string   common_sigCLK = common_node.attribute("sigCLK").value();

		if(!common_sigCLK.empty()) return str2dbl(common_sigCLK);

		if (_log) _log->comment(1, "WARNNING : no setting for sigCLK, use the defalut number : " + dbl2str(_default_sat_sigCLK));
		return _default_sat_sigCLK;
	}
	else
	{
		return str2dbl(sigCLK);
	}
}

map<string, double> great::t_gsetpar::sigSatPos(string sat_name)
{
	//get rkf model node
	xml_node satellite = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute(XMLKEY_PARS_SAT, "ID", sat_name.c_str());
	string   sigPOS = satellite.attribute("sigPOS").value();



	map<string, double> sigXYZ;
	if (!sigPOS.empty())
	{
		vector<string> sigmas;
		split(trim(sigPOS), "_", sigmas);
		if (sigmas.size() == 3)
		{
			sigXYZ.clear();
			sigXYZ.insert(pair<string, double>("sigPX", str2dbl(sigmas[0])));
			sigXYZ.insert(pair<string, double>("sigPY", str2dbl(sigmas[1])));
			sigXYZ.insert(pair<string, double>("sigPZ", str2dbl(sigmas[2])));
			return sigXYZ;
		}
	}

	if (_log) _log->comment(1, "WARNNING : no setting for sigPOS, use the defalut number : " + dbl2str(_default_sat_sigPOS));
	sigXYZ.clear();
	sigXYZ.insert(pair<string, double>("sigPX", _default_sat_sigPOS));
	sigXYZ.insert(pair<string, double>("sigPY", _default_sat_sigPOS));
	sigXYZ.insert(pair<string, double>("sigPZ", _default_sat_sigPOS));

	return sigXYZ;
}

map<string, double> great::t_gsetpar::sigSatVel(string sat_name)
{
	//get rkf model node
	xml_node satellite = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute(XMLKEY_PARS_SAT, "ID", sat_name.c_str());
	string   sigVEL = satellite.attribute("sigVEL").value();



	map<string, double> sigV;
	if (!sigVEL.empty())
	{
		vector<string> sigmas;
		split(trim(sigVEL), "_", sigmas);
		if (sigmas.size() == 3)
		{
			sigV.clear();
			sigV.insert(pair<string, double>("sigVX", str2dbl(sigmas[0])));
			sigV.insert(pair<string, double>("sigVY", str2dbl(sigmas[1])));
			sigV.insert(pair<string, double>("sigVZ", str2dbl(sigmas[2])));
			return sigV;
		}
	}

	if (_log) _log->comment(1, "WARNNING : no setting for sigVEL, use the defalut number : " + dbl2str(_default_sat_sigVEL));
	sigV.clear();
	sigV.insert(pair<string, double>("sigVX", _default_sat_sigVEL));
	sigV.insert(pair<string, double>("sigVY", _default_sat_sigVEL));
	sigV.insert(pair<string, double>("sigVZ", _default_sat_sigVEL));

	return sigV;
}

//
map<string, double> great::t_gsetpar::sigSatEcom(string sat_name)
{
	// todo : some problem
	//get rkf model node
	xml_node satellite = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute(XMLKEY_PARS_SAT, "ID", sat_name.c_str());
	
	string   sigECOM = satellite.attribute("sigECOM").value();

	map<string, double> sigPARS;
	if (!sigECOM.empty())
	{
		vector<string> sigmas;
		split(trim(sigECOM), "_", sigmas);
		if (sigmas.size() == 9)
		{
			sigPARS.clear();
			sigPARS.insert(pair<string, double>("sigECOM2_D0",  str2dbl(sigmas[0])));
			sigPARS.insert(pair<string, double>("sigECOM2_D2c", str2dbl(sigmas[1])));
			sigPARS.insert(pair<string, double>("sigECOM2_D2s", str2dbl(sigmas[2])));
			sigPARS.insert(pair<string, double>("sigECOM2_D4c", str2dbl(sigmas[3])));
			sigPARS.insert(pair<string, double>("sigECOM2_D4s", str2dbl(sigmas[4])));
			sigPARS.insert(pair<string, double>("sigECOM2_Y0",  str2dbl(sigmas[5])));
			sigPARS.insert(pair<string, double>("sigECOM2_X0",  str2dbl(sigmas[6])));
			sigPARS.insert(pair<string, double>("sigECOM2_X1c", str2dbl(sigmas[7])));
			sigPARS.insert(pair<string, double>("sigECOM2_X1s", str2dbl(sigmas[8])));
			return sigPARS;
		}
	}

	if (_log) _log->comment(1, "WARNNING : no setting for sigECOM, use the defalut number : " + dbl2str(_default_sat_sigECOM));
	sigPARS.clear();

	sigPARS.insert(pair<string, double>("sigECOM2_D0",  _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_D2c", _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_D2s", _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_D4c", _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_D4s", _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_Y0",  _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_X0",  _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_X1c", _default_sat_sigECOM));
	sigPARS.insert(pair<string, double>("sigECOM2_X1s", _default_sat_sigECOM));

	return sigPARS;
}

//
gnut::t_gtriple great::t_gsetpar::sigRecPosTriple(string sta_name)
{
	//get rkf model node
	xml_node satellite = _doc.child(XMLKEY_ROOT).child(XMLKEY_PARS).find_child_by_attribute(XMLKEY_PARS_STA, "ID", sta_name.c_str());
	string   sigPOS = satellite.attribute("sigPOS").value();



	gnut::t_gtriple sigXYZ(_default_sta_sigPOS, _default_sta_sigPOS, _default_sta_sigPOS);
	if (!sigPOS.empty())
	{
		vector<string> sigmas;
		split(trim(sigPOS), "_", sigmas);
		if (sigmas.size() == 3)
		{
			sigXYZ[0] = str2dbl(sigmas[0]);
			sigXYZ[1] = str2dbl(sigmas[1]);
			sigXYZ[2] = str2dbl(sigmas[2]);
			return sigXYZ;
		}
	}

	if (_log) _log->comment(1, "WARNNING : no setting for sigPOS, use the defalut number : " + dbl2str(_default_sta_sigPOS));
	return sigXYZ;
}

void great::t_gsetpar::help()
{
	cerr << "<!--> One example for this block : <!-->                       " << endl
		<< "<!--> The index of sigECOM is sigD0,sigDc,sigDs,sigY0,sigYC,sigYs,sigB0,sigBc,sigBs. <!-->" << endl
		<< "<parameters>" << endl
		<< "<GEO   sigCX    = \"0.001\"  sigCY    = \"0.001\" sigCZ  = \"0.001\" / >" << endl
		<< "<ERP   sigXPOLE = \"0.300\"  sigYPOLE = \"0.030\" sigUT1 = \"0.030\"  sigDXPOLE = \"0.0001\" sigDYPOLE = \"0.002\"  sigDUT1 = \"0.002\" />	 " << endl
		<< "<station   ID = \"AIRA\" sigCLK = \"9000\" sigZTD = \"0.201\"    sigION = \"9000\"     sigPOS = \"0.1_0.1_0.1\" / >" << endl
		<< "<satellite ID = \"G01\"  sigCLK = \"5000\" sigPOS = \"10_10_10\" sigVEL = \"10_10_10\"  sigECOM = \"0.1_0.1_0.1_0.1_0.1_0.1_0.1_0.1_0.1\"/>  " << endl
		<< "< / parameters>                          " << endl
		<< endl;
}

void  great::t_gsetpar::split(const std::string& s, std::string delim, std::vector< std::string >& ret)
{
	size_t last = 0;
	size_t index = s.find_first_of(delim, last);

	while (index != std::string::npos)
	{
		ret.push_back(s.substr(last, index - last));
		last = index + 1;
		index = s.find_first_of(delim, last);
	}
	if (index - last > 0)
	{
		ret.push_back(s.substr(last, index - last));
	}
}


