/**
*
* @verbatim
History
-1.0 glFeng   2019-02-27 creat the file.
-1.1 BoWong   2019-04-10 Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		   gsetamb.cpp
* @brief	   control set from XML
*
* @author      glFeng, Wuhan University
* @version	   1.0.0
* @date		   2019-04-10
*
*/

#include "gset/gsetamb.h"
#include <sstream>

using namespace std;
using namespace pugi;

namespace great
{
	/** @brief default constructor. */
	t_gsetamb::t_gsetamb()
	 : t_gsetbase()
	{
		_gmutex.lock();

		_set.insert(XMLKEY_AMBIGUITY);



		_gmutex.unlock();
	}

	/** @brief default destructor. */
	 t_gsetamb::~t_gsetamb()
	{

	}

	 /**
	 * @brief settings check.
	 * @return	  void
	 */
	 void t_gsetamb::check()
	 {
		 _gmutex.lock();
		 _gmutex.unlock();
		 return;
	 }

	 /**
	 * @brief settings help.
	 * @return	  void
	 */
	 void t_gsetamb::help()
	 {
		 _gmutex.lock();
		 cerr << "<ambiguity>\n"
			 << "<upd> EWL/WL/NL/NONE </upd>\n"
			 << "<fix_mode> ROUND/SEARCH/NO </fix_mode>\n"
			 << "<ratio> 2.0 </ratio>\n"
			 << "<all_baselines> NO </all_baselines>\n"
			 << "<min_common_time> 900 </min_common_time>\n"
			 << "<baseline_length_limit> 3500 </baseline_length_limit>\n"
			 << "<widelane_interval> 30 </widelane_interval>\n"
			 << "<widelane_decision     maxdev = \"0.25\" maxsig = \"0.15\" alpha = \"1000\"/>\n"
			 << "<narrowlane_decision   maxdev = \"0.15\" maxsig = \"0.33\" alpha = \"1000\"/>\n";

		 _gmutex.unlock();
		 return;
	 }
	 t_gsetamb::FIX_MODE t_gsetamb::fix_mode()
	 {
		 if (!_doc) str2fixmode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("fix_mode");
		 return str2fixmode(trim(tmp));
	 }

	 t_gsetamb::FIX_DD_MODE t_gsetamb::fix_dd_mode()
	 {
		if (!_doc){
			return FIX_DD_MODE::ALL;
		} 
		string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("fix_dd_mode");
		tmp=trim(tmp);
		if (tmp == "ALL"){
			return FIX_DD_MODE::ALL;
		}
		else if (tmp == "IND"){
			return FIX_DD_MODE::IND;
		}
		else{
			return FIX_DD_MODE::ALL;
		}
	 }
	 DD_MODEL t_gsetamb::dd_mode()
	 {
		 if (!_doc) str2fixmode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("dd_mode");

		 DD_MODEL mode;
		      if (tmp.find("IF_CB_WN")  != string::npos) mode = DD_MODEL::IF_CB_WN;
		 else if (tmp.find("RAW_CB_WN") != string::npos) mode = DD_MODEL::RAW_CB_WN;
		 else if (tmp.find("RAW_CB")    != string::npos) mode = DD_MODEL::RAW_CB;
		 else if (tmp.find("NONE" )     != string::npos) mode = DD_MODEL::NONE;
		 else
		 {
			 throw logic_error("wrong type of dd mode in node amb.");
		 }
			 
		 return mode;
	 }

	 double t_gsetamb::lambda_ratio()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("ratio");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }


	 bool t_gsetamb::independ_baseline()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("all_baselines"));
		 string tmp;
		 is >> tmp;
		 bool isindepend_baseline = (tmp == "NO" || tmp == "no");

		 _gmutex.unlock(); return isindepend_baseline;
	 }

	 double t_gsetamb::min_common_time()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("min_common_time");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 double t_gsetamb::max_baseline_length()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("baseline_length_limit");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }


	 double t_gsetamb::wl_interval() 
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("widelane_interval");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 map<string, double> t_gsetamb::get_amb_decision(string type)
	 {
		 _gmutex.lock();

		 xml_node tmp_set;
		 map<string, double> amb_decision;
		 if (type == "EWL") {
			 tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child("extra_widelane_decision");
			 string tmp = tmp_set.attribute("maxdev").value();
			 amb_decision["maxdev"] = tmp.empty() ? 0.07 : str2dbl(tmp);
			 tmp = tmp_set.attribute("maxsig").value();
			 amb_decision["maxsig"] = tmp.empty() ? 0.10 : str2dbl(tmp);
			 tmp = tmp_set.attribute("alpha").value();
			 amb_decision["alpha"] = tmp.empty() ? 1000 : str2dbl(tmp);
		 }
		 else if (type == "WL") {
			 tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child("widelane_decision");
			 string tmp = tmp_set.attribute("maxdev").value();
			 amb_decision["maxdev"] = tmp.empty() ? 0.25 : str2dbl(tmp);
			 tmp = tmp_set.attribute("maxsig").value();
			 amb_decision["maxsig"] = tmp.empty() ? 0.10 : str2dbl(tmp);
			 tmp = tmp_set.attribute("alpha").value();
			 amb_decision["alpha"] = tmp.empty() ? 1000 : str2dbl(tmp);
		 }
		 else if (type == "NL") {
			 tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child("narrowlane_decision");
			 string tmp = tmp_set.attribute("maxdev").value();
			 amb_decision["maxdev"] = tmp.empty() ? 0.25 : str2dbl(tmp);
			 tmp = tmp_set.attribute("maxsig").value();
			 amb_decision["maxsig"] = tmp.empty() ? 0.10 : str2dbl(tmp);
			 tmp = tmp_set.attribute("alpha").value();
			 amb_decision["alpha"] = tmp.empty() ? 1000 : str2dbl(tmp);
		 }

		 _gmutex.unlock(); return amb_decision;
	 }



	 /**
	  * @brief change string to FIXMODE.
	  * @param[in]  str   fix mode string
	  * @return	  FIXMODE
	  */
	 t_gsetamb::FIX_MODE t_gsetamb::str2fixmode(string str)
	 {
		 if (str == "NO") 
		 {
			 return t_gsetamb::FIX_MODE::NO;
		 }
		 else if (str == "ROUND") 
		 {
			 return t_gsetamb::FIX_MODE::ROUND;
		 }
		 else if (str == "SEARCH") 
		 {
			 return t_gsetamb::FIX_MODE::SEARCH;
		 }
		 else if (str == "HOLD")
		 {
			 return t_gsetamb::FIX_MODE::HOLD;
		 }
		 else
		 {
			 cout << "*** warning: not defined ambiguity fixing mode [" << str << "]\n"; cout.flush();
		 }
		 return t_gsetamb::FIX_MODE::NO;
	 }
}
