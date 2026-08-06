/**
*
* @verbatim
	History
	 -1.0 GREAT  	 2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gcfg_podflt.h
* @brief		LSQ control set from XML
*
* @author		GREAT, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#include <iomanip>
#include <sstream>

#include "gcfg_podflt.h"

using namespace std;
using namespace pugi;

namespace great {

	// Constructor
	// ----------
	t_gcfg_podflt::t_gcfg_podflt()
		: t_gsetbase(),
		t_gsetgen(),
		t_gsetinp(),
		t_gsetout(),
		t_gsetgnss(),
		t_gsetproc(),
		t_gsetrec(),
		t_gsetpar(),
		t_gsetamb(),
		t_gsetforce()
	{
		_IFMT_supported.insert(RINEXO_INP);
		_IFMT_supported.insert(RINEXC_INP);
		_IFMT_supported.insert(RINEXN_INP);
		_IFMT_supported.insert(ATX_INP);
		_IFMT_supported.insert(BLQ_INP);
		_IFMT_supported.insert(SP3_INP);
		_IFMT_supported.insert(BIABERN_INP);
		_IFMT_supported.insert(IFCB_INP);
		_IFMT_supported.insert(SINEX_INP);
		_IFMT_supported.insert(ORB_INP);

		_IFMT_supported.insert(ICS_INP);
		_IFMT_supported.insert(ION_INP);
        _IFMT_supported.insert(SRIF_VAR_INP);

		_IFMT_supported.insert(AMBCON_INP);

		_IFMT_supported.insert(OCEANTIDE_INP);
		_IFMT_supported.insert(ERP_INP);
		_IFMT_supported.insert(SATPARS_INP);
		_IFMT_supported.insert(LEAPSECOND_INP);
		_IFMT_supported.insert(PREOBS_INP);
		_IFMT_supported.insert(SLRF_INP);
		_IFMT_supported.insert(LRA_INP);
		_IFMT_supported.insert(DE_INP);
		_IFMT_supported.insert(RB_INP);
		_IFMT_supported.insert(ATTITUDE_INP);
		_IFMT_supported.insert(EGM_INP);
		_IFMT_supported.insert(SOLAR_INP);
		_IFMT_supported.insert(DESAISCOPOLECOEF_INP);
        _IFMT_supported.insert(AMBFLAG12_INP);
        _IFMT_supported.insert(AMBFLAG13_INP);
        _IFMT_supported.insert(AMBFLAG14_INP);
        _IFMT_supported.insert(AMBFLAG15_INP);

		_OFMT_supported.insert(LOG_OUT);
		_OFMT_supported.insert(PPP_OUT);
		_OFMT_supported.insert(FLT_OUT);
        _OFMT_supported.insert(SRIF_VAR_OUT);
		_OFMT_supported.insert(SP3_OUT);
		_OFMT_supported.insert(RECOVER_OUT);
		_OFMT_supported.insert(RECCLK_OUT);
		_OFMT_supported.insert(SATCLK_OUT);
		_OFMT_supported.insert(ICS_OUT);

		_OFMT_supported.insert(RB_OUT);
		_OFMT_supported.insert(STA_OUT);
		_OFMT_supported.insert(ORB_OUT);
	}


	// Destructor
	// ----------
	t_gcfg_podflt::~t_gcfg_podflt()
	{}


	// settings check
	// ----------
	void t_gcfg_podflt::check()
	{
		t_gsetgen::check();
		t_gsetinp::check();
		t_gsetout::check();
		t_gsetrec::check();
		t_gsetproc::check();
		t_gsetgnss::check();
	
	}

	// settings help
	// ----------
	void t_gcfg_podflt::help()
	{
		t_gsetbase::help_header();
		t_gsetgen::help();
		t_gsetinp::help();
		t_gsetout::help();
		t_gsetrec::help();
		t_gsetproc::help();
		t_gsetgnss::help();
	
		t_gsetbase::help_footer();
	}

	bool t_gcfg_podflt::check_inp()
	{
		string out = "check your xml file with input node, there is no : ";


		// CHECK INPUTS
		bool valid = true;

		auto tmp = this->inputs_all();
		if (tmp.empty())
		{
			cerr << out + "file" << endl;
			valid = false;
		}

		if (this->input_size("orb") == 0 && this->input_size("sp3") == 0 && this->input_size("rnn")==0)
		{
			cerr << out + "orb/sp3/rinexn" << endl;
			valid = false;
		}

		if ((this->input_size("rinexc") == 0 && this->input_size("rnc") ==0)
			&& (this->input_size("rinexn") == 0 && this->input_size("rnn") ==0))
		{
			cerr << out + "rinexc or rinexn" << endl;
			valid = false;
		}

		if (this->input_size("rinexo") == 0 && this->input_size("rno")==0)
		{
			cerr << out + "rinexo" << endl;
			valid = false;
		}
		if (this->input_size("erp") == 0 && this->input_size("erp") ==0)
		{
			cerr << out + "erp" << endl;
			valid = false;
		}

		if (this->input_size("satpars") == 0 && this->input_size("sat")==0)
		{
			cerr << out + "satpars" << endl;
			valid = false;
		}
		if (this->input_size("atx") == 0)
		{
			cerr << out + "atx" << endl;
			valid = false;
		}
		return valid;
	}

	void t_gcfg_podflt::arg(int argc, char* argv[], bool add, bool thin)
	{
		string conf("");
		string save("");
		bool stdi = false; // ENFORCE STDINP
		bool thrd = false;

		//cerr << "command-line opt: " << argc << endl;  for( int i=1; i<argc; ++i ) cerr << argv[i] << " "; cerr << endl;

		if (argc == 1)
		{
			cerr << "This is a app used for pod of gns satellite.      " << endl
				<< "your can use  --help get how to use the cmd.      " << endl
				<< "         use  --X    get how to set the xml file. " << endl
				<< "         use  --x    to  use the xml file and cal." << endl
				<< "         use  --v    get the verb information.    " << endl;
			return;
		}

		for (int i = 1; i < argc; ++i)
		{
			if (!strcmp(argv[i], "--help") && i < argc) { usage(); exit(0); }
			else if (!strcmp(argv[i], "-h") && i < argc) { usage(); exit(0); }
			else if (!strcmp(argv[i], "-V") && i < argc) { cout << app() << endl; exit(0); }
			else if (!strcmp(argv[i], "-X") && i < argc) { help(); exit(0); }
			else if (!strcmp(argv[i], "-m") && i < argc) { thrd = true; }
			else if (!strcmp(argv[i], "-brdm") && i < argc) { isUseBrdm = true; }
			else if (!strcmp(argv[i], "-ambfix") && i < argc) { isAmbfix = true; }
			else if (!strcmp(argv[i], "-ion") && i < argc) { isUseIon = true; }
			else if (!strcmp(argv[i], "-x") && i + 1 < argc) { conf = argv[++i]; }
			else if (!strcmp(argv[i], "-Z") && i + 1 < argc) { save = argv[++i]; }
			else if (!strcmp(argv[i], "-v") && i + 1 < argc) { if (_log) { _log->verb(atoi(argv[++i])); } } // ZDE PRILIS BRZY!!
			else if (!strcmp(argv[i], "-l") && i + 1 < argc) { if (_log) { _log->mask(argv[++i]); } } // ZDE PRILIS BRZY!!
			//warning if no other commandline option is supported
			else if (!add) fprintf(stderr, "Unknown option: %s\n", argv[i]);
			else continue;
		}

		if (!conf.empty())
		{
			if (read(conf)) { exit(1); } // file XML
		}

		// add/replace command-line option in XML
		for (int i = 1; i < argc; ++i)
		{
			xml_node NODE = _default_node(_doc, "config", "", true);

			string s = argv[i]; size_t pos1 = 0, pos2 = 0;
			if (s[0] == '-') continue;

			
			bool isNode = s.find("node") != string::npos;
			bool isAttr = s.find("attr") != string::npos;
			
			if (isNode || isAttr) s = s.substr(5);

			while (((pos1 = s.find(':', pos2)) != string::npos ||    // i+1 < argc ... used later
				(pos1 = s.find('=', pos2)) != string::npos))     // ORDER IMPORTANT
			{
				if ((pos2 = s.find(':', pos1 + 1)) != string::npos)
				{
					
					// example : reset value
					// config:gen:int=30 reset int to 30
					//cout << s.substr(pos1 + 1, pos2 - pos1 - 1) << endl;
					NODE = _default_node(NODE, s.substr(pos1 + 1, pos2 - pos1 - 1).c_str(), "", true);
				}
				else if ((pos2 = s.find('=', pos1 + 1)) != string::npos)
				{
				
					//cout << s.substr(pos1 + 1, pos2 - pos1 - 1) << endl;
					//cout << s.substr(pos2 + 1, string::npos) << endl;
					if (isAttr) _default_attr(NODE, s.substr(pos1 + 1, pos2 - pos1 - 1).c_str(), s.substr(pos2 + 1, string::npos), true);
					else if (isNode) _default_node(NODE, s.substr(pos1 + 1, pos2 - pos1 - 1).c_str(), s.substr(pos2 + 1, string::npos).c_str(), true); //jdhuang
					else _default_node(NODE, s.substr(pos1 + 1, pos2 - pos1 - 1).c_str(), s.substr(pos2 + 1, string::npos).c_str(), true);

					pos2 = string::npos;

				}
				else
				{
					if (i + 1 < argc) _default_node(NODE, s.substr(pos1 + 1, string::npos).c_str(), argv[++i], true); // only if exists
					else           cerr << "Incomplete command-line argument\n";
					pos2 = string::npos;
				}
			}
		}

		if (thin)
		{
			xml_node node = _doc.child(XMLKEY_ROOT);
			_default_attr(node, "thin", "true", true);
		}

		this->check();

		if (!conf.empty()) { if (_log) _log->comment(0, "gsetbase", "xconfig: read from file " + _name); }
		else { if (_log) _log->comment(0, "gsetbase", "xconfig: read from istream"); }

		if (!save.empty()) {
			if (_log) _log->comment(0, "gsetbase", "xconfig: save to file " + _name);
			write(save);
		}

		// COMMENTS
		if (_log && !thin) {
			_log->comment(0, "gsetbase", "program: " + app());
			_log->comment(0, "gsetbase", "==========");
			_log->comment(1, "gsetbase", "setting: " + conf);
			_log->comment(1, "gsetbase", "logfile: " + _log->mask());
			_log->comment(1, "gsetbase", "verbose: " + int2str(_log->verb()));
			_log->comment(1, "gsetbase", "==========");
		}

		return;

	}

} // namespace
