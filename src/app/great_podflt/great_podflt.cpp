/**
*
* @verbatim
	History
	 -1.0 GREAT  	    2019-01-04 creat the file.
	 -1.1 Huang Jiande  2020-09-12 Add clock repair for gobs
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			great_podflt.cpp
* @brief		The main function of lsq.
*
* @author		GREAT, Wuhan University
* @version		1.0.0
* @date			2019-01-04
*
*/

#include "gproc/gpodfltIF.h"
#include "gcfg_podflt.h"
#include "gcdr_podflt.h"

#include <chrono>
#include <thread>
#include <stdio.h>
#include <signal.h>



using namespace std;
using namespace gnut;
using namespace great;
using namespace std::chrono;

void catch_signal(int) { cout << "Program interrupted by Ctrl-C [SIGINT,2]\n"; exit(1);}


// MAIN
// ----------
int main(int argc, char** argv)
{
	// ================================================================================================================================================================================================================
	signal(SIGINT, catch_signal);                                                                  // Only to cout the Reminder here



	// ================================================================================================================================================================================================================
	t_gcfg_podflt gset;                                                                            // Construct the gset class and init some values in the class
	gset.app("GREAT/POD_FLT", "0.9.0", "$Rev: 2448 $", "(WHU-SGG GREAT)", __DATE__, __TIME__);     // App information
	gset.arg(argc, argv, true, false);                                                             // Get the arguments from the command line
	
	t_glog glog;                                                                                   // Creat and set the log file : clk.log
	string log_filename = dynamic_cast<t_gsetout*>(&gset)->outputs("log");
	if (log_filename == "") log_filename = "LOGRT.app_log";
	substitute(log_filename, "LOGRT", "great_podlflt", false);
	glog.mask(log_filename);
	glog.cache_size(99);
	glog.tsys(t_gtime::GPS);
	glog.time_stamp(true);
	glog.verb(dynamic_cast<t_gsetout*>(&gset)->verb());
	if (glog.verb() > 0) {
		glog.compare(false);
	}
	
	// set ambfix by process node
	gset.isAmbfix = gset.isAmbfix || dynamic_cast<t_gsetproc*>(&gset)->ambfix();

	// ================================================================================================================================================================================================================
	t_gcdr_podflt gcdr; gcdr.log(&glog);   gcdr.set(&gset);                                       // coder class for reading file,and  add log file
	multimap<IFMT, string> inp = gset.inputs_all();                                               // Prepare input files list form gset
	if (!gset.check_inp()) return false;                                                          // then check the input file

	// DECLARATIONS/INITIALIZATIONS	
	map<t_gdata::ID_TYPE, t_gdata*> map_gdata;
//	t_gallrecover* grcv = new t_gallrecover(); grcv->glog(&glog);    map_gdata[grcv->id_type()] = grcv;
	t_gallpcv* gpcv = new t_gallpcv;  gpcv->glog(&glog);             map_gdata[gpcv->id_type()] = gpcv;   // gpcv for the atx data read from the atx file
	t_gallotl* gotl = new t_gallotl;  gotl->glog(&glog);             map_gdata[gotl->id_type()] = gotl;   // gotl for the blq data read from the blq file, which will be used for Ocean tidal corrections
	t_gallobj* gobj = new t_gallobj(gpcv, gotl); gobj->glog(&glog);  map_gdata[gobj->id_type()] = gobj;   // gobj for the gpcv and gotl, which means the model can be used by all satellites and stations
	
	// SET OBJECTS
	set<string> recs  = dynamic_cast<t_gsetgen*>(&gset)->rec_all();
	set<string> objs  = dynamic_cast<t_gsetrec*>(&gset)->objects();
	for (const auto& it_rec : recs) { t_gobj* ptrObj = dynamic_cast<t_gsetrec*>(&gset)->grec_all(it_rec); if (ptrObj) gobj->add(shared_ptr<t_gobj>(ptrObj)); }
	for (const auto& it_obj : objs) { shared_ptr<t_grec> ptrRec = dynamic_cast<t_gsetrec*>(&gset)->grec(it_obj, &glog); if (ptrRec) gobj->add(ptrRec); }

	t_gtime beg = dynamic_cast<t_gsetgen*>(&gset)->beg();                                          // Get beg and end time in xml.
	t_gtime end = dynamic_cast<t_gsetgen*>(&gset)->end();
	auto rom = dynamic_cast<t_gsetproc*>(&gset)->read_ofile_mode();
																								   // Read files
	int obs_true = 0; int obs_false = 0;
	for (const auto& it_inp : inp)
	{
		if (!gcdr.support(it_inp.first))
		{
			glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "can not support to read the file " + it_inp.second);
			if (it_inp.first == RINEXO_INP) { obs_false++; }
			continue;
		}
		gcdr.creat(it_inp.first, map_gdata);
		auto type = gcdr.ifmt2data(it_inp.first);

		if (rom == OFSTREAM_MODE::REAL_TIME) {
			if (it_inp.first == IFMT::RINEXO_INP || it_inp.first == IFMT::RINEXN_INP) {
				continue;
			}
			// read in process 
			//gcdr.coder(beg, beg + step, it_inp.first, it_inp.second, map_gdata[type], gobj); 
		}
		obs_true++; gcdr.coder(it_inp.first, it_inp.second, map_gdata[type], gobj);
	}
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "obs find is : " + int2str(obs_true) + " " + dbl2str((obs_true/(obs_false+obs_true)*100)) + " %");
	if (rom != OFSTREAM_MODE::REAL_TIME) {
		gobj->read_satinfo(beg);                                                                       	// set antennas for satllites (must be before PCV assigning)
		gobj->sync_pcvs();                                                                              // assigning PCV pointers to objects; PCV comes from atx file
	}
	// ADD DATA
	t_gallproc* data = new t_gallproc(); for (auto iter : map_gdata) { data->Add_Data(iter.second); }
	if ((*data)[t_gdata::ALLORB] == nullptr) {
		data->Add_Data(new t_gallorb());
	}
	
	// ================================================================================================================================================================================================================
	// orb PROCESSING
	t_gfltproc* vgorb = nullptr;
	set<string> mode = dynamic_cast<t_gsetgen*>(&gset)->mode();
	bool process_normal=true;

	t_gtime epo(t_gtime::GPS);
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "processing started [ POD Estimation ]");
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", beg.str_ymdhms("  beg: ") + end.str_ymdhms("  end : "));
	// ================================================================================================================================================================================================================

	// ================================================================================================================================================================================================================
	// GNS IF
	vgorb = new t_gpodfltIF(&gset, data, &glog);
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "NOTE   : IF Combination mode.");
	vgorb->setUseAmbcon(gset.isAmbfix);
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "Ambfix : " + bl2str(gset.isAmbfix));
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "Brdm   : " + bl2str(gset.isUseBrdm));

	// Start write log for the processing 
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", " processing started [ Satellite Clock Estimation ]");
	glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", beg.str_mjdsod("  beg: ") + end.str_mjdsod("  end: "));


	// =================================================================================
	auto sat_list = vgorb->sat_list();
	auto sys_list = vgorb->sys_list();
	auto rec_list = vgorb->rec_list();
	if (sat_list.empty() || sys_list.empty() || rec_list.empty())
	{
		glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", " skip for sat_list.empty() || sys_list.empty() || rec_list.empty()");
	}
	else
	{
		// =================================================================================
		// Mark for the begin time of processing
		// Main processing function
		// update clk file && ics file
		// Get the time of duration 
		t_gtime nav_runepoch = t_gtime::current_time(t_gtime::GPS);
		if (vgorb->ProcessBatch(data, beg, end))
		{
			t_gtime nav_lstepoch = t_gtime::current_time(t_gtime::GPS);

			// Write finished log
			glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", " processing finished : duration  " + dbl2str(nav_lstepoch.diff(nav_runepoch)) + " sec");
			glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "normal end");
		}
		else{
			process_normal=false;
			glog.comment(t_glog::LOG_LV::LOG_INFO, "great_podflt", "main", "Error!");
		}
	}
			
	// ==================================================================================
	// delete var.
	if (vgorb) { delete vgorb; vgorb = nullptr; }
	if (data)  { delete data;  data = nullptr;  }

	for (auto& iter : map_gdata) { if (iter.second) { delete iter.second; iter.second = nullptr;} }
	if (process_normal){
		return 0;
	}
	else{
		return 1;
	}
}