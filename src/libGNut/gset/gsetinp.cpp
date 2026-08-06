/**
*
* @verbatim
(c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
Ondrejov 244, 251 65, Czech Republic
@endverbatim
*
* @file		gsetinp.cpp
* @brief	implements input setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetinp.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut 
{

	// Convertor for INP formats
	// ----------
	IFMT t_gsetinp::str2ifmt(const string& s)
	{
		string tmp = s;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		if (tmp == "RINEXC") return RINEXC_INP;
		if (tmp == "CAS")    return RINEXC_INP; 
		if (tmp == "CAR")    return RINEXC_INP; 
		
		if (tmp == "RINEXC_PRD") return RINEXC_PRD;
		if (tmp == "RINEXC_REF") return RINEXC_REF;
		if (tmp == "RINEXO") return RINEXO_INP;
		if (tmp == "RNO") return RINEXO_INP;

		if (tmp == "RINEXN") return RINEXN_INP;
		if (tmp == "RNN") return RINEXN_INP;

		if (tmp == "RINEXN_PRD") return RINEXN_PRD;
		if (tmp == "RINEXN_REF") return RINEXN_REF;
		if (tmp == "RINEXM") return RINEXM_INP;
		if (tmp == "BNCOBS") return BNCOBS_INP;
		if (tmp == "BNCBRDC") return BNCBRDC_INP;
		if (tmp == "BNCCORR") return BNCCORR_INP;
		if (tmp == "BNCSEND") return BNCSEND_INP; 
		if (tmp == "BNCRTCM") return BNCRTCM_INP;
		if (tmp == "BNCRTCM_PRD") return BNCRTCM_PRD;
		if (tmp == "BNCRTCM_REF") return BNCRTCM_REF;
		if (tmp == "TROSINEX") return TROSINEX_INP;
		if (tmp == "TROSINEX0") return TROSINEX0_INP;
		if (tmp == "SINEX")     return SINEX_INP;
		if (tmp == "SNX")       return SINEX_INP; 
		if (tmp == "TROGRID")   return TROGRID_INP;
		if (tmp == "RTCM") return RTCM_INP;
		if (tmp == "TSMET") return TSMET_INP;
		if (tmp == "SP3") return SP3_INP;
		if (tmp == "SP3_PRD") return SP3_PRD;
		if (tmp == "SP3_REF") return SP3_REF;
		if (tmp == "ATX") return ATX_INP;
        if (tmp == "SRIF_VAR") return SRIF_VAR_INP;

		if (tmp == "OTL") return BLQ_INP; // BLQ == ocean load
		if (tmp == "BLQ") return BLQ_INP;
		if (tmp == "FCB") return FCB_INP;
		if (tmp == "RB") return RB_INP;
		if (tmp == "STA") return STA_INP;
		if (tmp == "KBR")    return KBR_INP;
		if (tmp == "GRIB_CHMI") return GRIB_CHMI_INP;
		if (tmp == "GRIB_HARM") return GRIB_HARM_INP;
		if (tmp == "GRIB_ERA") return GRIB_ERA_INP;
		if (tmp == "NCD3_ICS") return NCD3_ICS_INP;

		if (tmp == "RAO_TEXT") return RAO_TEXT_INP;
		if (tmp == "RAO_IGRA") return RAO_IGRA_INP;
		if (tmp == "RAO_BADC") return RAO_BADC_INP;

		if (tmp == "BIASINEX") return BIASINEX_INP;
		if (tmp == "BIABERN") return BIABERN_INP;
		if (tmp == "DCB")     return BIABERN_INP;

		if (tmp == "IONEX") return IONEX_INP;

		// special for  orb PPP-RTK
		if (tmp == "AUG")return AUG_INP;

		// special for  orb intergal
		if (tmp == "JPL") return DE_INP;
		if (tmp == "DE")  return DE_INP;

		if (tmp == "OCEANTIDE") return OCEANTIDE_INP;
		if (tmp == "DESAISCOPOLECOEF") return DESAISCOPOLECOEF_INP;
		if (tmp == "OTD")       return OCEANTIDE_INP;

		if (tmp == "EGM")       return EGM_INP;
		if (tmp == "SOLAR")return SOLAR_INP;
		if (tmp == "ICS")return ICS_INP;
		if (tmp == "ION")return ION_INP;

		// special for POD 
		if (tmp == "ORB") return ORB_INP;
		if (tmp == "erp")return ERP_INP;
		if (tmp == "ERP")    return ERP_INP; // optinal for xml node

		if (tmp == "LEP")    return LEAPSECOND_INP;
 		if (tmp == "LEAPSECOND")return LEAPSECOND_INP;

		if (tmp == "SATPARS")return SATPARS_INP;
		if (tmp == "SAT")    return SATPARS_INP; // optinal for xml node
		
		if (tmp == "ATTITUDE") return ATTITUDE_INP;
		if (tmp == "AMBCON") return AMBCON_INP;
		if (tmp == "CON")    return AMBCON_INP; // optinal for xml node
		if (tmp == "AMBUPD")  return AMBUPD_INP;
		
		if (tmp == "AMBFLAG")   return AMBFLAG12_INP;
		if (tmp == "AMBFLAG13") return AMBFLAG13_INP;
		if (tmp == "AMBFLAG14") return AMBFLAG14_INP;
		if (tmp == "AMBFLAG15") return AMBFLAG15_INP;

		if (tmp == "FLG")       return AMBFLAG12_INP;
		if (tmp == "FLG13")     return AMBFLAG13_INP;
		if (tmp == "FLG14")     return AMBFLAG14_INP;
		if (tmp == "FLG15")     return AMBFLAG15_INP;
		
		if (tmp == "UPD") return UPD_INP;
		if (tmp == "XTR") return XTR_INP;

		// special for recover pars
		if (tmp == "RECOVER") return RECOVER_INP;
		if (tmp == "RCV")     return RECOVER_INP; // optianl for xml node
		if (tmp == "PANDA_RECOVER") return PANDA_RECOVER_INP;
		if (tmp == "PREOBS") return PREOBS_INP;

		if (tmp == "IMU") return IMU_INP;  
		if (tmp == "CAM") return IMAGE_INP;
		if (tmp == "LIDAR") return LIDAR_INP;

		if (tmp == "POS") return IMU_POS_INP;
		
		if (tmp == "IFCB") return IFCB_INP;
		if (tmp == "LRA") return LRA_INP;
		if (tmp == "NPT") return NPT_INP;
		if (tmp == "SLRF") return SLRF_INP;
		return IFMT(-1);
	}


	// Convertor for INP formats
	// ----------
	string t_gsetinp::ifmt2str(const IFMT& f)
	{
		switch (f) {
		case RINEXO_INP:  return "RINEXO";
		case RINEXC_INP:  return "RINEXC";
		case RINEXC_PRD:  return "RINEXC_PRD";
		case RINEXC_REF:  return "RINEXC_REF";
		case RINEXN_INP:  return "RINEXN";
		case RINEXN_PRD:  return "RINEXN_PRD";
		case RINEXN_REF:  return "RINEXN_REF";
		case RINEXM_INP:  return "RINEXM";
		case BNCOBS_INP:  return "BNCOBS";
		case BNCBRDC_INP:  return "BNCBRDC";
		case BNCCORR_INP:  return "BNCCORR";
		case BNCSEND_INP:  return "BNCSEND"; 
		case BNCRTCM_INP: return "BNCRTCM";
		case BNCRTCM_PRD: return "BNCRTCM_PRD";
		case BNCRTCM_REF: return "BNCRTCM_REF";
		case TROSINEX_INP:return "TROSINEX";
		case TROSINEX0_INP:return "TROSINEX0";
		case SINEX_INP:   return "SINEX";
		case TROGRID_INP: return "TROGRID";
		case RTCM_INP:    return "RTCM";
		case TSMET_INP:   return "TSMET";
		case SP3_INP:     return "SP3";
		case PREOBS_INP:     return "PREOBS";
		case LRA_INP:     return "LRA";
		case RB_INP:      return "RB";
		case KBR_INP:     return "KBR";
		case STA_INP:      return "STA";
		case NPT_INP:     return "NPT";
		case SLRF_INP:     return "SLRF";
		case ION_INP:     return "ION";
		case ICS_INP:     return "ICS";
        case SRIF_VAR_INP:   return "SRIF_VAR";
		case SP3_PRD:     return "SP3_PRD";
		case SP3_REF:     return "SP3_REF";
		case ATX_INP:     return "ATX";
		case BLQ_INP:     return "BLQ";
		case ERP_INP: return "ERP";
		case FCB_INP:     return "FCB";
		case GRIB_CHMI_INP:  return "GRIB_CHMI";
		case GRIB_HARM_INP:  return "GRIB_HARM";
		case GRIB_ERA_INP:   return "GRIB_ERA";
		case NCD3_ICS_INP:   return "NCD3_ICS";
		case RAO_TEXT_INP:   return "RAO_TEXT";
		case RAO_IGRA_INP:   return "RAO_IGRA";
		case RAO_BADC_INP:   return "RAO_BADC";
		case BIASINEX_INP:   return "BIASINEX";
		case BIABERN_INP:    return "BIABERN";
		case IONEX_INP:      return "IONEX";
		case IMU_INP:		return "IMU";
		case IMAGE_INP:		return "CAM";
		case LIDAR_INP:  return "LIDAR";
		case IMU_POS_INP:	return "POS";
		case AMBUPD_INP:     return "AMBUPD";
		case AMBFLAG12_INP:    return "AMBFLAG";
		case AMBFLAG13_INP:  return "AMBFLAG13";
		case AMBFLAG14_INP:  return "AMBFLAG14";
		case AMBFLAG15_INP:  return "AMBFLAG15";
		case AMBCON_INP:		return "AMBCON";
		case UPD_INP:        return "UPD";
		case RECOVER_INP:	return "RECOVER";
		case PANDA_RECOVER_INP: return "PANDA_RECOVER";
		
		case IFCB_INP:        return "IFCB";
		case LEAPSECOND_INP:  return "LEAPSECOND";
		case SOLAR_INP:       return"SOLAR";
		case ATTITUDE_INP:    return"ATTITUDE";
		case DESAISCOPOLECOEF_INP:return"desaiscopolecoef";
		case DE_INP:     return "DE";
		default:             return "UNDEF";
		}
		return "UNDEF";
	}

	// Constructor
	// ----------
	t_gsetinp::t_gsetinp()
		: t_gsetbase()
	{
		_set.insert(XMLKEY_INP);
	}


	// Destructor
	// ----------
	t_gsetinp::~t_gsetinp()
	{}


	// Get formats input size
	// ----------
	int t_gsetinp::input_size(const string& fmt)
	{
		_gmutex.lock();

		int tmp = _inputs(str2ifmt(fmt)).size();

		_gmutex.unlock(); return tmp;
	}



	// Get formats inputs (all in multimap)
	// ----------
	multimap<IFMT, string> t_gsetinp::inputs_all()
	{
		_gmutex.lock();

		multimap<IFMT, string> map;

		set<string> ifmt = _iformats();
		set<string>::const_iterator itFMT = ifmt.begin();

		while (itFMT != ifmt.end())
		{
			string fmt = *itFMT;
			IFMT  ifmt = str2ifmt(fmt);
			vector<string> inputs = _inputs(fmt);//get file name in input node
			vector<string>::const_iterator itINP = inputs.begin();
			while (itINP != inputs.end()) 
			{
				map.insert(map.end(), pair<IFMT, string>(ifmt, *itINP));
				itINP++;
			}
			itFMT++;
		}
		_gmutex.unlock(); return map;
	}




	// Get formats inputs
	// ----------
	vector<string> t_gsetinp::inputs(const string& fmt)
	{

	
		IFMT ifmt = str2ifmt(fmt);
		return _inputs(ifmt);
	}

	vector<string> t_gsetinp::inputs(const IFMT & ifmt)
	{
		return _inputs(ifmt);
	}


	// Get input formats
	// ----------
	set<string> t_gsetinp::iformats()
	{
		return _iformats();
	}


	// Get formats inputs
	// ----------
	vector<string> t_gsetinp::_inputs(const string& fmt)
	{
		vector<string> tmp;
		set<string> list;
		string str;

		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
		{
			if (node.name() == fmt) 
			{
				istringstream is(node.child_value());
				while (is >> str && !is.fail())
				{
					if (str.find("://") == string::npos) str = GFILE_PREFIX + str;
					if (list.find(str) == list.end()) 
					{
						tmp.push_back(str);
						list.insert(str);
					}
					else 
					{
						if (_log) _log->comment(1, "gsetinp", "READ: " + str + " multiple request ignored");
					}
				}
			}
		}
		return tmp;
	}

	vector<string> t_gsetinp::_inputs(const IFMT &fmt)
	{
		vector<string> tmp;
		set<string> list;
		string str;

		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
		{
			// jdhuang
			IFMT ifmt = str2ifmt(node.name());
			if (ifmt == IFMT::UNDEF) continue;
			if (ifmt == fmt)
			{
				istringstream is(node.child_value());
				while (is >> str && !is.fail())
				{
					if (str.find("://") == string::npos)
						str = GFILE_PREFIX + str;
					if (list.find(str) == list.end())
					{
						tmp.push_back(str);
						list.insert(str);
					}
					else
					{
						if (_log)
							_log->comment(1, "gsetinp", "READ: " + str + " multiple request ignored");
					}
				}
			}
		}
		return tmp;
	}

	// Get input formats
	// ----------
	set<string> t_gsetinp::_iformats()
	{
		set<string> tmp;
		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
		{
			tmp.insert(node.name());
		}
		return tmp;
	}



	// settings check
	// ----------
	void t_gsetinp::check()
	{
		_gmutex.lock();

		// check existence of nodes/attributes
		xml_node parent = _doc.child(XMLKEY_ROOT);
		xml_node node = _default_node(parent, XMLKEY_INP);

		// check supported input formats (see IFMT enum !)
		set<string> ifmt = _iformats();
		set<string>::const_iterator itFMT = ifmt.begin();
		while (itFMT != ifmt.end()) {
			string fmt = *itFMT;
			IFMT  ifmt = str2ifmt(fmt);
			if (ifmt < 0) {
				_doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child(node.child(fmt.c_str()));
				if (_log) _log->comment(0, "Warning: " + fmt + " inp format not implemented [gsetinp::check()]!");
				itFMT++;
				continue;
			}
			// check application-specific output format
			if (_IFMT_supported.find(ifmt) == _IFMT_supported.end()) {
				_doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child(node.child(fmt.c_str()));
				if (_log) _log->comment(0, "Warning: " + fmt + " inp format not supported by this application!");
				else                cerr << "Warning: " + fmt + " inp format not supported by this application!\n";
			}
			itFMT++;
		}


		_gmutex.unlock(); return;
	}


	// settings help
	// ----------
	void t_gsetinp::help()
	{
		_gmutex.lock();


		cerr << " <inputs>\n"
			<< "   <rinexo> file://dir/name </rinexo> \t\t <!-- obs RINEX decoder -->\n"
			<< "   <rinexn> file://dir/name </rinexn> \t\t <!-- nav RINEX decoder -->\n"
			<< " </inputs>\n";

		cerr << "\t<!-- inputs description:\n"
			<< "\t <decoder> path1 path2 path3  </decoder>\n"
			<< "\t ... \n"
			<< "\t where path(i) contains [file,tcp,ntrip]:// depending on the application\n"
			<< "\t -->\n\n";

		_gmutex.unlock(); return;
	}

} // namespace
