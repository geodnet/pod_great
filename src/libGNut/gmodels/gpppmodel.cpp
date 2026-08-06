
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *  
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include "gmodels/gpppmodel.h"

namespace gnut {

   // Constructors
   // -------------------
   t_gpppmodel::t_gpppmodel(string site, t_glog* glog, t_gsetbase* settings)
      : t_gsppmodel(site, glog, settings),
      _gallobj(0),
	  _attitude_last_matrix(Matrix(3,3))
   {
      gtrace("t_gpppmodel::t_gpppmodel");
      if (_trpModStr == EXTERN)  _tropoModel = make_shared<t_gtropo>();

      _phase = dynamic_cast<t_gsetproc*>(_settings)->phase();
      _grad_mf = dynamic_cast<t_gsetproc*>(_settings)->grad_mf();
	  _trop_corr = dynamic_cast<t_gsetproc*>(_settings)->trop_corr();
	  _attitudes = dynamic_cast<t_gsetproc*>(_settings)->attitudes();
   }

   t_gpppmodel::t_gpppmodel()
   {
   }

   // Destructor
   // --------------------------
   t_gpppmodel::~t_gpppmodel()
   {
      gtrace("t_gpppmodel::~t_gpppmodel");
   }


   // wind-up correction
   // ----------
   double t_gpppmodel::windUp(t_gsatdata& satdata, const ColumnVector& rRec)
   {
     gtrace("t_gpppmodel::windUp");

     t_gtime epoch = satdata.epoch();
     string prn = satdata.sat();
     ColumnVector rSat = satdata.satcrd().crd_cvect();
     //GSYS sys          = satdata.gsys();

     double Mjd = epoch.mjd(true) + epoch.sod(true) / 86400.0;

     // First time - initialize to zero
     // -------------------------------
     map<string, double>::iterator it = _windUpTime.find(prn);
     if (it == _windUpTime.end()){
       _windUpTime[prn] = Mjd;
       _windUpSum[prn] = 0.0;

     }
     
     // Compute the correction for new time
     // -----------------------------------
     else if (_windUpTime[prn] != Mjd) {
       
       _windUpTime[prn] = Mjd;
       
       ColumnVector rho = rRec - rSat;
       rho /= rho.norm_Frobenius();
       
       ColumnVector i(3);
       ColumnVector j(3);
       ColumnVector k(3);

       // attitude model
       string antype = "";
       if (_gallobj != 0){
         shared_ptr<t_gobj>  sat_obj = _gallobj->obj(satdata.sat());
         shared_ptr<t_gpcv>  sat_pcv;
         if (sat_obj != 0)  sat_pcv = sat_obj->pcv(satdata.epoch());
         if (sat_pcv != 0)  antype = sat_pcv->anten();
       }
	   if(_attitudes==YAW_NOMI)		   attitude(satdata, "", i, j, k);				//nominal modeling
       else if(_attitudes == YAW_RTCM) attitude(satdata, satdata.yaw(), i, j, k);	//value from RTCM used
	   else							   attitude(satdata, antype, i, j, k);			//default
	   //

       // Effective Dipole of the GPS Satellite Antenna
       // ---------------------------------------------
       ColumnVector dipSat = i - rho * DotProduct(rho, i) - crossproduct(rho, j);
       
       // Receiver unit Vectors rx, ry
       // ----------------------------
       ColumnVector rx(3);
       ColumnVector ry(3);
       
       double recEll[3]; xyz2ell(rRec.data(), recEll, false);
       double neu[3];
       
       neu[0] = 1.0;
       neu[1] = 0.0;
       neu[2] = 0.0;
       neu2xyz(recEll, neu, rx.data());
       
       neu[0] = 0.0;
       neu[1] = -1.0;
       neu[2] = 0.0;
       neu2xyz(recEll, neu, ry.data());

       // Effective Dipole of the Receiver Antenna
       // ----------------------------------------
       ColumnVector dipRec = rx - rho * DotProduct(rho, rx) + crossproduct(rho, ry);
       
       // Resulting Effect
       // ----------------
       double alpha = DotProduct(dipSat, dipRec) / (dipSat.norm_Frobenius() * dipRec.norm_Frobenius());
       
       if (alpha > 1.0) alpha = 1.0;
       if (alpha < -1.0) alpha = -1.0;
       
       double dphi = acos(alpha) / 2.0 / G_PI;  // in cycles

       if (DotProduct(rho, crossproduct(dipSat, dipRec)) < 0.0) dphi = -dphi;
   
       _windUpSum[prn] = floor(_windUpSum[prn] - dphi + 0.5) + dphi;
     }
     

     satdata.addwind(_windUpSum[prn]);
     return _windUpSum[prn];
   }
   

   // set gallobj pointer
   // --------------------
   void t_gpppmodel::setOBJ(t_gallobj* obj)
   {
     _gallobj = obj;
     if (obj){
       _grec = obj->obj(_site);         
     }
   }

  // model computed range value (phase/code)
  // ----------
  double t_gpppmodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& gsatdata, t_gobs& gobs, bool com)
  {
    gtrace("t_gpppmodel::cmpObs");

    double spp_model = t_gsppmodel::cmpObs(epoch, param, gsatdata, gobs);
    if (spp_model < 0) return -1;
	GSYS gs = gsatdata.gsys();
	GOBSBAND b1 = t_gsys::band_priority(gs, FREQ_1);
	GOBSBAND b2 = t_gsys::band_priority(gs, FREQ_2);
	// modified by glfeng consider tide correction, earth rotation
	// -------------------------------------------
    // Cartesian coordinates to ellipsodial coordinates
	t_gtriple xyz_r, xyz_s, ell;
	ColumnVector cRec(3);
	xyz_s = gsatdata.satcrd();
	if (param.getCrdParam(_site, xyz_r) < 0) 
	{
		xyz_r = gsatdata.reccrd();
		if (double_eq(xyz_r[0], 0.0) || double_eq(xyz_r[1], 0.0) || double_eq(xyz_r[2], 0.0)) {
			xyz_r = _grec->crd_arp(epoch);
		}
	}
	cRec = xyz_r.crd_cvect();
	xyz2ell(xyz_r, ell, false);

	// Earth rotation correction
	t_gtriple Rec_pco;
	double rho0 = sqrt(pow(xyz_r[0] - xyz_s[0], 2) + pow(xyz_r[1] - xyz_s[1], 2) + pow(xyz_r[2] - xyz_s[2], 2));
	double dPhi = OMEGA * rho0 / CLIGHT;
	Rec_pco[0] = xyz_r[0] * cos(dPhi) - xyz_r[1] * sin(dPhi);
	Rec_pco[1] = xyz_r[1] * cos(dPhi) + xyz_r[0] * sin(dPhi);
	Rec_pco[2] = xyz_r[2];

	// -------------------------------------------

    // Wind up correction 
    double wind = 0.0;
    if (gobs.is_phase()) 
	{        
      double wavelength = 0.0;
           
     wavelength = gsatdata.wavelength_L3(b1, b2);

	  _attitudes = DEF_YAWMODEL;   // add by glfeng
      if (fabs(gsatdata.wind()) > 0) wind = gsatdata.wind() * wavelength;
      else wind = windUp(gsatdata, cRec) * wavelength;
    }     
	_attitudes = YAW_NOMI;
    
    // Phase center variation correction
    double pcv_R = 0.0;  
    double pcv_S = 0.0;  
    double pco_R = 0.0;  
    double pco_S = 0.0;  
    if (_gallobj != 0)
	{
      shared_ptr<t_gobj>  sat_obj = _gallobj->obj(gsatdata.sat());
      shared_ptr<t_gobj> site_obj = _gallobj->obj(_site);
      
      shared_ptr<t_gpcv>  sat_pcv;
      shared_ptr<t_gpcv> site_pcv;

      if (sat_obj != 0)  sat_pcv  = sat_obj->pcv(epoch);
      if (site_obj != 0) site_pcv = site_obj->pcv(epoch);

      GOBS_LC lc = LC_IF;
  
      if(sat_pcv != 0 && com)
	  {
        // Satellite phase center variation
		if (gobs.is_phase()) sat_pcv->pcvS(pcv_S, gsatdata, xyz_r, lc,b1,b2);  // range don't consider PCV 
        // Satellite phase center offset
        t_gtriple pco(0,0,0);
        if(sat_pcv->pcoS(gsatdata, pco, lc,b1,b2) > 0) 
		{
          string antenna = sat_pcv->anten();
          t_gtriple dx(0,0,0);
          ColumnVector i(3);
          ColumnVector j(3);
          ColumnVector k(3);        
		  if (_attitudes == YAW_NOMI)	   attitude(gsatdata, "", i, j, k);				//nominal modeling
		  else if (_attitudes == YAW_RTCM) attitude(gsatdata, gsatdata.yaw(), i, j, k);	//value from RTCM used
		  else							   attitude(gsatdata, antenna, i, j, k);		//default
          dx[0] = pco[0]*i(1) + pco[1]*j(1) + pco[2]*k(1);
          dx[1] = pco[0]*i(2) + pco[1]*j(2) + pco[2]*k(2);
          dx[2] = pco[0]*i(3) + pco[1]*j(3) + pco[2]*k(3);
          sat_pcv->pco_proj(pco_S, gsatdata, Rec_pco, dx);// change xyz-> Rec_pco consider earth rotation 
		  gsatdata.addpco(dx);
        }
        
      }
      
      if(site_pcv != 0)
	  {
        // Receiver phase center variation
		if (gobs.is_phase()) site_pcv->pcvR(pcv_R, gsatdata, lc,b1,b2);// range don't consider PCV 
        // Receiver phase center offset
        t_gtriple dx(0.0, 0.0, 0.0);          
        if(site_pcv->pcoR(gsatdata, dx, xyz_r, lc) > 0) {
          site_pcv->pco_proj(pco_R, gsatdata, Rec_pco, dx);  // change xyz-> Rec_pco consider earth rotation 
        }
        pco_R *= -1; 
      }
    }

      // Return value
      return spp_model +
        wind  +
        pco_R +
        pco_S +        
        pcv_R +
        pcv_S;
	 
  }
  

   // Compute troposperic delay
   // -----------
   double t_gpppmodel::tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata)
   {
      gtrace("t_gpppmodel::tropoDelay");

      if (_tropoModel == 0){
         if (_log) _log->comment(0, "gppp", "Tropo Model setting is not correct. Default used! Check config.");
         else                cerr << "gppp - Tropo Model setting is not correct. Default used! Check config.\n";
         _tropoModel = make_shared<t_saast>();
      }

      double ele = satdata.ele();
      double azi = satdata.azi();

      double delay = 0.0;
      double zwd = 0.0;
      double zhd = 0.0;
      int i, j, k;
      i = param.getParam(_site, t_gpar::TRP, "");
      j = param.getParam(_site, t_gpar::GRD_N, "");
      k = param.getParam(_site, t_gpar::GRD_E, "");

      if (i >= 0){
         zwd = param[i].value();
         zhd = param[i].apriori();
		
      }
      else{
		  if (_tropoModel != 0) {
			  zwd = _tropoModel->getZWD(ell, epoch);
			  zhd = _tropoModel->getZHD(ell, epoch);
		  }
      }

      double mfh, mfw, dmfh, dmfw;

      t_gmf mf;
      mf.gmf(epoch.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
      mfh, mfw, dmfh, dmfw);


      satdata.addmfH(mfh);
      satdata.addmfW(mfw);

      delay = mfh * zhd + mfw * zwd;



      double grdN, grdE;
      grdN = grdE = 0.0;

      if (j >= 0 && k >= 0){
         double mfg = 1.0 / (sin(ele)*tan(ele) + 0.0032);
         grdN = param[j].value() * 1000.0 * mfg * cos(azi);  grdN /= 1000.0;
         grdE = param[k].value() * 1000.0 * mfg * sin(azi);  grdE /= 1000.0;
         satdata.addmfG(mfg);
         delay += grdN + grdE;
      }



      return delay;
   }

// satellite attitude model
// --------------------------------   
int t_gpppmodel::attitude_old(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    if(satdata.gsys() == GAL){      
         _ysm(satdata, i, j, k);
    }else{
      int irc = _yaw(satdata,antype, i, j, k);
      if (irc == 0)return 0;
   }
    
   return 1;
}

// satellite attitude model
// --------------------------------   
int t_gpppmodel::attitude(t_gsatdata& satdata, double yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
	_yaw2ijk(satdata, yaw, i, j, k);
	return 1;
}

// satellite attitude model
// --------------------------------   
int t_gpppmodel::attitude(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{ 

  ColumnVector satcrd = satdata.satcrd().crd_cvect();
  ColumnVector satvel = satdata.satvel().crd_cvect();       

  if(satcrd.NormFrobenius() == 0 || satvel.NormFrobenius() == 0){ return -1; }

  if(satdata.gsys() == GPS){
    if(antype == "BLOCK II"){
      _attitude_GPSIIA(satdata, antype, i, j, k);
    }else if(antype == "BLOCK IIA"){
      _attitude_GPSIIA(satdata, antype, i, j, k);
    }else if(antype.find("BLOCK IIR") != string::npos){
      _attitude_GPSIIR(satdata, antype, i, j, k);
    }else if(antype == "BLOCK IIF"){
      _attitude_GPSIIF(satdata, antype, i, j, k);
    }else if (antype == "BLOCK III"){  // added by yqyuan for GPSIII (BLOCK III/IIIA/IIIF???)
	  _attitude_GPSIII(satdata, antype, i, j, k);
	}
    else _ysm(satdata, i, j, k);  
  }
  else if(satdata.gsys() == GLO){    
    _attitude_GLO(satdata, antype, i, j, k);
         
    }
  else _ysm(satdata, i, j, k);
    
  
  if(i.NormFrobenius() == 0 || j.NormFrobenius() == 0 || k.NormFrobenius() == 0) return 0;
   

  return 1;
}   


// satellite attitude model ( with input Xsat, Vsat, Xsun) used in OI
// --------------------------------   

int t_gpppmodel::attitude(string antype, string prn, ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
	ColumnVector satcrd = xsat;
	ColumnVector satvel = vsat;

	double orb_angle = _orb_angle(xsat, vsat, xsun);
	double beta = _beta(xsat, vsat, xsun);

	if (satcrd.NormFrobenius() == 0 || satvel.NormFrobenius() == 0) { return -1; }

	if (prn.find("G") != string::npos) {
		if (antype == "BLOCK II") {
			//_attitude_GPSIIA(satdata, antype, i, j, k);
			_attitude_GPSIIA(antype, prn, xsat, vsat, xsun, i, j, k);
		}
		else if (antype == "BLOCK IIA") {
			//_attitude_GPSIIA(satdata, antype, i, j, k);
			_attitude_GPSIIA(antype, prn, xsat, vsat, xsun, i, j, k);
		}
		else if (antype.find("BLOCK IIR") != string::npos) {
			//_attitude_GPSIIR(satdata, antype, i, j, k);
			_attitude_GPSIIR(antype, prn, xsat, vsat, xsun, i, j, k);
		}
		else if (antype == "BLOCK IIF") {
			//_attitude_GPSIIF(satdata, antype, i, j, k);
			_attitude_GPSIIF(antype, prn, xsat, vsat, xsun, i, j, k);
		}
		else if (antype == "BLOCK III") {  
			_attitude_GPSIII(antype, prn, xsat, vsat, xsun, i, j, k);
		}
		else _ysm(prn, beta, orb_angle, xsat, vsat, xsun, i, j, k); 
	}
	else if (prn.find("R") != string::npos) {
		_attitude_GLO(antype, prn, xsat, vsat, xsun, i, j, k);
	}
	else _ysm(prn, beta, orb_angle, xsat, vsat, xsun, i, j, k);

	if (i.NormFrobenius() == 0 || j.NormFrobenius() == 0 || k.NormFrobenius() == 0) return 0;
	return 1;
}

double t_gpppmodel::_orb_angle(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun)
{
	double mi = 0.0;
	ColumnVector Sun = xsun;   //ICRF
	ColumnVector Satcrd = xsat;
	ColumnVector Satvel = vsat;

	ColumnVector n = crossproduct(Satcrd, Satvel);

	ColumnVector es = Satcrd / Satcrd.NormFrobenius();

	ColumnVector eSun = Sun / Sun.NormFrobenius();

	ColumnVector p = crossproduct(Sun, n);
	p /= p.NormFrobenius();

	double E = acos(DotProduct(es, p));

	double SunSat = acos(DotProduct(es, eSun));

	if (SunSat > G_PI / 2) {
		if (E <= G_PI / 2) mi = G_PI / 2 - E;
		else            mi = G_PI / 2 - E;
	}
	else {
		if (E <= G_PI / 2) mi = G_PI / 2 + E;
		else            mi = E - G_PI - G_PI / 2;
	}

	return mi;
}

double t_gpppmodel::_beta(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun)
{
	double beta = 0.0;

	ColumnVector Sun = xsun;

	ColumnVector Satcrd = xsat;
	ColumnVector Satvel = vsat;

	ColumnVector n = crossproduct(Satcrd, Satvel);

	n /= n.NormFrobenius();

	ColumnVector nSun = Sun / Sun.NormFrobenius();

	double cosa = DotProduct(nSun, n);


	beta = G_PI / 2.0 - acos(cosa);

	return beta;
}

   
// Yaw-steering mode attitude model
// --------------------------------   
void t_gpppmodel::_ysm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
   double MJD = satdata.epoch().dmjd();
   
   t_gtriple satcrd_t  = satdata.satcrd();   
   ColumnVector satcrd = satcrd_t.crd_cvect();

   // Satelite-Earth unit vector
   k = -satcrd;

   // Along solar panel unit vector  
   ColumnVector sun = _ephplan.sunPos(MJD).crd_cvect();
    j = crossproduct(k,sun);

   // complete to satelite fixed right-hand coord system
   i = crossproduct(j,k);
    
    i = i / i.NormFrobenius();
    j = j / j.NormFrobenius();
    k = k / k.NormFrobenius();    

    _last_beta[satdata.sat()] = satdata.beta();
    
    double yaw = atan2(-satdata.beta(), sin(satdata.orb_angle()));
    satdata.yaw(yaw);
}

void t_gpppmodel::_ysm(string prn, double bata, double mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	ColumnVector satcrd = xsat;

	// Satelite-Earth unit vector
	k = -satcrd;

	// Along solar panel unit vector  
	ColumnVector sun = xsun;
	j = crossproduct(k, sun);

	// complete to satelite fixed right-hand coord system
	i = crossproduct(j, k);

	i = i / i.NormFrobenius();
	j = j / j.NormFrobenius();
	k = k / k.NormFrobenius();

	_last_beta[prn] = bata;

	double yaw = atan2(-bata, sin(mi));
}


// yaw model of different satellite
// Nominal attitude (same as GPS BLOCK II/IIA) without yaw maneuver for MEO and IGSO satellites ();
// Yaw - fixed attitude mode used for GEO satellites
// --------------------------------   
int t_gpppmodel::_yaw(t_gsatdata& satdata,string antype, ColumnVector& xs, ColumnVector& ys, ColumnVector& zs)
{

   return 1;
}


// Attitude modelling for GPS Block IIA
void t_gpppmodel::_attitude_GPSIIA(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R_GPSIIA[] = {
         0.1046,0.1230,0.1255,0.1249,0.1003,0.1230,0.1136,0.1169,0.1253,0.0999,
         0.1230,0.1230,0.1230,0.1230,0.1092,0.1230,0.1230,0.1230,0.1230,0.1230,
         0.1230,0.1230,0.1230,0.0960,0.0838,0.1284,0.1183,0.1230,0.1024,0.1042,
         0.1230,0.1100,0.1230
    };
    int sat = str2int(satdata.sat().substr(1,2));
    double R = R_GPSIIA[sat - 1] * D2R;

    double beta0 = atan2(MUDOT_GPS, R);
    
    if(fabs(satdata.orb_angle()) > G_PI/2 && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn(satdata, R, i, j, k);
    }else if(fabs(satdata.orb_angle()) < G_PI/2 && fabs(satdata.beta()) < beta0){
         // midnight maneuver
         _midnight_turn_GPSIIA(satdata, R, i, j, k);
    }else{
         _ysm(satdata, i, j, k); 
    }  

}

void t_gpppmodel::_attitude_GPSIIA(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	const double R_GPSIIA[] = {
		0.1046,0.1230,0.1255,0.1249,0.1003,0.1230,0.1136,0.1169,0.1253,0.0999,
		0.1230,0.1230,0.1230,0.1230,0.1092,0.1230,0.1230,0.1230,0.1230,0.1230,
		0.1230,0.1230,0.1230,0.0960,0.0838,0.1284,0.1183,0.1230,0.1024,0.1042,
		0.1230,0.1100,0.1230
	};
	int sat = str2int(prn.substr(1, 2));
	double R = R_GPSIIA[sat - 1] * D2R;

	double beta0 = atan2(MUDOT_GPS, R);
	double orb_angle = _orb_angle(xsat, vsat, xsun);
	double beta = _beta(xsat, vsat, xsun);

	if (fabs(orb_angle) > G_PI / 2 && fabs(beta) < beta0) {
		// noon maneuver
		_noon_turn(prn, beta, orb_angle, xsat, vsat, xsun, R, i, j, k);
	}
	else if (fabs(orb_angle) < G_PI / 2 && fabs(beta) < beta0) {
		// midnight maneuver
		_midnight_turn_GPSIIA(prn, beta, orb_angle, xsat, vsat, xsun, R, i, j, k);
	}
	else {
		_ysm(prn, beta, orb_angle, xsat, vsat, xsun, i, j, k);
	}
}


// Attitude modelling for GPS Block IIR    
void t_gpppmodel::_attitude_GPSIIR(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R = 0.2*D2R;             // maximal yaw hardware rate
    double beta0 = atan2(MUDOT_GPS, R);
    
    if(fabs(satdata.orb_angle()) > G_PI/2 && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn(satdata, R, i, j, k);
    }else if(fabs(satdata.orb_angle()) < G_PI/2 && fabs(satdata.beta()) < beta0){
         // midnight maneuver
         _midnight_turn(satdata, R, i, j, k);
    }else{              
         _ysm(satdata, i, j, k);
    }  

    i = -1 * i;      // X away from the Sum
    j = -1 * j;
    satdata.yaw(satdata.yaw() + G_PI);  
    
}


void t_gpppmodel::_attitude_GPSIIR(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	const double R = 0.2*D2R;             // maximal yaw hardware rate
	double beta0 = atan2(MUDOT_GPS, R);
	double orb_angle = _orb_angle(xsat, vsat, xsun);
	double beta = _beta(xsat, vsat, xsun);

	if (fabs(orb_angle) > G_PI / 2 && fabs(beta) < beta0) {
		// noon maneuver
		_noon_turn(prn, beta, orb_angle, xsat, vsat, xsun, R, i, j, k);
	}
	else if (fabs(orb_angle) < G_PI / 2 && fabs(beta) < beta0) {
		// midnight maneuver
		_midnight_turn_GPSIIA(prn, beta, orb_angle, xsat, vsat, xsun, R, i, j, k);
	}
	else {
		_ysm(prn, beta, orb_angle, xsat, vsat, xsun, i, j, k);
	}

	i = -1 * i;      // X away from the Sum
	j = -1 * j;
}

    
// noon maneuver
// -------------------------------------------
void t_gpppmodel::_noon_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta0 = atan2(MUDOT_GPS, R);
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }

    if (beta >= 0) R *= -1;
    
    double mi_s = G_PI - sqrt( beta0 * fabs(beta) - beta*beta );
    double yaw;
    double yaw_nom = atan2(-tan(beta), sin(mi));

    if(mi >= mi_s || mi <= 0){
         if(mi < 0) mi += 2.0*G_PI;
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
         if ((beta >= 0 && yaw < yaw_nom) || (beta < 0 && yaw > yaw_nom)) yaw = yaw_nom;

    }else yaw = yaw_nom;
    
    _yaw2ijk(satdata, yaw, i, j, k);
}
    
// midnight maneuver
// -------------------------------------------
void t_gpppmodel::_midnight_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta0 = atan2(MUDOT_GPS, R);
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }
    
    if (beta < 0) R *= -1;
    
    double mi_s = - sqrt( beta0 * fabs(beta) - beta*beta );
    double yaw;
    double yaw_nom = atan2(-tan(beta), sin(mi));
    
    if(mi >= mi_s){
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
         if ((beta >= 0 && yaw > yaw_nom) || (beta < 0 && yaw < yaw_nom)) yaw = yaw_nom;
    }else yaw = yaw_nom;
    
    _yaw2ijk(satdata, yaw, i, j, k);       
}

void t_gpppmodel::_noon_turn(string _prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	double beta0 = atan2(MUDOT_GPS, R);
	double beta = _beta;
	double mi = _mi;

	// test of beta sign change
	auto itSAT = _last_beta.find(_prn);
	if (itSAT != _last_beta.end()) {
		if (itSAT->second * beta < 0) beta *= -1;
	}

	if (beta >= 0) R *= -1;

	double mi_s = G_PI - sqrt(beta0 * fabs(beta) - beta*beta);
	double yaw;
	double yaw_nom = atan2(-tan(beta), sin(mi));

	if (mi >= mi_s || mi <= 0) {
		if (mi < 0) mi += 2.0*G_PI;
		yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
		if ((beta >= 0 && yaw < yaw_nom) || (beta < 0 && yaw > yaw_nom)) yaw = yaw_nom;

	}
	else yaw = yaw_nom;

	_yaw2ijk(xsat, vsat, xsun, yaw, i, j, k);
}

void t_gpppmodel::_midnight_turn(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
}

// midnight maneuver
// -------------------------------------------
void t_gpppmodel::_midnight_turn_GPSIIA(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    double mi_s = - sqrt(EPS0_GPS*EPS0_GPS - beta*beta);
    double mi_e = - mi_s;         

    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }
    
    if (beta < 0) R *= -1;
    
    double yaw;
    
    if (mi_s <= mi && mi < mi_e) {
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
    }else if(mi_e <= mi && mi < mi_e + POST_SHADOW*MUDOT_GPS) satdata.setecl(true);
    else yaw = atan2(-tan(beta), sin(mi));
    
    _yaw2ijk(satdata, yaw, i, j, k);       
}   

void t_gpppmodel::_midnight_turn_GPSIIA(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	double beta = _beta;
	double mi = _mi;
	double mi_s = -sqrt(EPS0_GPS*EPS0_GPS - beta*beta);
	double mi_e = -mi_s;

	// test of beta sign change
	auto itSAT = _last_beta.find(prn);
	if (itSAT != _last_beta.end()) {
		if (itSAT->second * beta < 0) beta *= -1;
	}

	if (beta < 0) R *= -1;

	double yaw;

	if (mi_s <= mi && mi < mi_e) {
		yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
	}
	else if (mi_e <= mi && mi < mi_e + POST_SHADOW*MUDOT_GPS);//satdata.setecl(true);
	else yaw = atan2(-tan(beta), sin(mi));

	_yaw2ijk(xsat, vsat, xsun, yaw, i, j, k);
}
    

    
// midnight maneuver for GPS Block IIF
// -------------------------------------------
void t_gpppmodel::_midnight_turn_GPSIIF(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double  tan_beta = tan(beta);
    double mi   = satdata.orb_angle();
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }
    
    if (beta < 0) R *= -1;
    
    double mi_s = -acos(cos(EPS0_GPS) / cos(beta)); 
    double mi_e = - mi_s;
    double sin_mi_s = sin(mi_s);
    double mi_f = MUDOT_GPS*(atan2(-tan_beta, -sin_mi_s) - atan2(-tan_beta, sin_mi_s)) / R + mi_s;

    double yaw;
    
    if (mi_s <= mi && mi < mi_f) {
         yaw = atan2(-tan_beta, sin_mi_s) + R*(mi - mi_s) / MUDOT_GPS;
         _yaw2ijk(satdata, yaw, i, j, k);     
      }
      else if (mi_f <= mi && mi < mi_e) {
         yaw = atan2(-tan_beta, -sin_mi_s);
          _yaw2ijk(satdata, yaw, i, j, k);       
      }else _ysm(satdata, i, j, k);
    

}


void t_gpppmodel::_midnight_turn_GPSIIF(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	double beta = _beta;
	double  tan_beta = tan(beta);
	double mi = _mi;

	// test of beta sign change
	auto itSAT = _last_beta.find(prn);
	if (itSAT != _last_beta.end()) {
		if (itSAT->second * beta < 0) beta *= -1;
	}

	if (beta < 0) R *= -1;

	double mi_s = -acos(cos(EPS0_GPS) / cos(beta));
	double mi_e = -mi_s;
	double sin_mi_s = sin(mi_s);
	double mi_f = MUDOT_GPS*(atan2(-tan_beta, -sin_mi_s) - atan2(-tan_beta, sin_mi_s)) / R + mi_s;

	double yaw;

	if (mi_s <= mi && mi < mi_f) {
		yaw = atan2(-tan_beta, sin_mi_s) + R*(mi - mi_s) / MUDOT_GPS;
		_yaw2ijk(xsat, vsat, xsun, yaw, i, j, k);
	}
	else if (mi_f <= mi && mi < mi_e) {
		yaw = atan2(-tan_beta, -sin_mi_s);
		_yaw2ijk(xsat, vsat, xsun, yaw, i, j, k);
	}
	else _ysm(prn, beta, mi, xsat, vsat, xsun, i, j, k);

}

    
// midnight maneuver for GLONASS-M
// -------------------------------------------
void t_gpppmodel::_midnight_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double  tan_beta = tan(beta);
    double mi   = satdata.orb_angle();
    
    if (beta < 0) R *= -1;
    
    double mi_s = -acos(cos(EPS0_GLO) / cos(beta)); 
    double mi_e = - mi_s;
    double sin_mi_s = sin(mi_s);
    double mi_f = MUDOT_GLO*(atan2(-tan_beta, -sin_mi_s) - atan2(-tan_beta, sin_mi_s)) / R + mi_s;

	// debug by ZHJ initial
    double yaw=0.0;
    
    if (mi_s <= mi && mi < mi_f) {
         yaw = atan2(-tan_beta, sin_mi_s) + R*(mi - mi_s) / MUDOT_GLO;
      }
      else if (mi_f <= mi && mi < mi_e) {
         yaw = atan2(-tan_beta, -sin_mi_s);
      }      
    
    _yaw2ijk(satdata, yaw, i, j, k);
}   

void t_gpppmodel::_midnight_turn_GLOM(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	double beta = _beta;
	double  tan_beta = tan(beta);
	double mi = _mi;

	if (beta < 0) R *= -1;

	double mi_s = -acos(cos(EPS0_GLO) / cos(beta));
	double mi_e = -mi_s;
	double sin_mi_s = sin(mi_s);
	double mi_f = MUDOT_GLO*(atan2(-tan_beta, -sin_mi_s) - atan2(-tan_beta, sin_mi_s)) / R + mi_s;

	// debug by ZHJ initial
	double yaw = 0.0;

	if (mi_s <= mi && mi < mi_f) {
		yaw = atan2(-tan_beta, sin_mi_s) + R*(mi - mi_s) / MUDOT_GLO;
	}
	else if (mi_f <= mi && mi < mi_e) {
		yaw = atan2(-tan_beta, -sin_mi_s);
	}

	_yaw2ijk(xsat, vsat, xsun, yaw, i, j, k);
}


// noon maneuver
// -------------------------------------------
void t_gpppmodel::_noon_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    if (beta >= 0) R = -R;
    
    double mi_s, mi_e, sin_mi_s, B, yaw;
    int c = 0;

    for(mi_s = 178.6*D2R; c < 4; c++) {
         sin_mi_s = sin(mi_s);
         B = -beta*cos(mi_s) / (beta*beta + sin_mi_s*sin_mi_s);
         mi_s = (atan(-beta / sin_mi_s) + mi_s*B + G_PI*R / MUDOT_GLO - G_PI / 2.0) / (R / MUDOT_GLO + B);
    }
    if(beta >= 0) mi_s = 2.0*G_PI - mi_s;
      mi_e = 2.0*G_PI - mi_s;

      if (mi_s <= mi && mi < mi_e) {
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GLO;
      }
    
    _yaw2ijk(satdata, yaw, i, j, k);    
}

void t_gpppmodel::_noon_turn_GLOM(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	double beta = _beta;
	double mi = _mi;

	if (beta >= 0) R = -R;

	double mi_s, mi_e, sin_mi_s, B, yaw;
	int c = 0;

	for (mi_s = 178.6*D2R; c < 4; c++) {
		sin_mi_s = sin(mi_s);
		B = -beta*cos(mi_s) / (beta*beta + sin_mi_s*sin_mi_s);
		mi_s = (atan(-beta / sin_mi_s) + mi_s*B + G_PI*R / MUDOT_GLO - G_PI / 2.0) / (R / MUDOT_GLO + B);
	}
	if (beta >= 0) mi_s = 2.0*G_PI - mi_s;
	mi_e = 2.0*G_PI - mi_s;

	if (mi_s <= mi && mi < mi_e) {
		yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GLO;
	}

	_yaw2ijk(xsat, vsat, xsun, yaw, i, j, k);
}
    

    
// Attitude modelling for GPS Block IIR-M
void t_gpppmodel::_attitude_GPSIIRM(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    _attitude_GPSIIR(satdata, antype, i, j, k);
}
    
// Attitude modelling for GPS Block IIF
void t_gpppmodel::_attitude_GPSIIF(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R_noon = 0.11*D2R;             // maximal yaw hardware rate during noon turn
    const double R_midn = 0.06*D2R;             // maximal yaw hardware rate during midnight turn
    double beta0 = atan2(MUDOT_GPS, R_noon);

    if(fabs(satdata.orb_angle()) > G_PI/2 && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn(satdata, R_noon, i, j, k);
    }else if(fabs(satdata.orb_angle()) < EPS0_GPS && fabs(satdata.beta()) < EPS0_GPS){
         // midnight maneuver
         _midnight_turn_GPSIIF(satdata, R_midn, i, j, k);
    }else{
         _ysm(satdata, i, j, k); 
    }  
}

void t_gpppmodel::_attitude_GPSIIF(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	const double R_noon = 0.11*D2R;             // maximal yaw hardware rate during noon turn
	const double R_midn = 0.06*D2R;             // maximal yaw hardware rate during midnight turn
	double beta0 = atan2(MUDOT_GPS, R_noon);

	double orb_angle = _orb_angle(xsat, vsat, xsun);
	double beta = _beta(xsat, vsat, xsun);

	if (fabs(orb_angle) > G_PI / 2 && fabs(beta) < beta0) {
		// noon maneuver
		_noon_turn(prn, beta, orb_angle, xsat, vsat, xsun, R_noon, i, j, k);
	}
	else if (fabs(orb_angle) < EPS0_GPS && fabs(beta) < EPS0_GPS) {
		// midnight maneuver
		_midnight_turn_GPSIIF(prn, beta, orb_angle, xsat, vsat, xsun, R_midn, i, j, k);
	}
	else {
		_ysm(prn, beta, orb_angle, xsat, vsat, xsun, i, j, k);
	}
}
   
// Attitude modelling for GPS III, added by yqyuan but not realized.
void t_gpppmodel::_attitude_GPSIII(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k) {
   // to be added; 
   _ysm(satdata, i, j, k);
}

void t_gpppmodel::_attitude_GPSIII(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
}


    
// Attitude modelling for GLONASS
void t_gpppmodel::_attitude_GLO(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R = 0.25*D2R;             // maximal yaw hardware rate
    double beta0 = 2.0*D2R;
    double mi = satdata.orb_angle();
    double mi_s = 176.8*D2R;
    
    if(mi > mi_s && mi < 2.0*G_PI - mi_s && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn_GLOM(satdata, R, i, j, k);
    }else if(fabs(satdata.orb_angle()) < EPS0_GLO && fabs(satdata.beta()) < EPS0_GLO){
         // midnight maneuver
         _midnight_turn_GLOM(satdata, R, i, j, k);
    }else{
         _ysm(satdata, i, j, k); 
    }  
}   

void t_gpppmodel::_attitude_GLO(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k)
{
	const double R = 0.25*D2R;             // maximal yaw hardware rate
	double beta0 = 2.0*D2R;

	double beta = _beta(xsat, vsat, xsun);
	double mi = _orb_angle(xsat, vsat, xsun);
	double mi_s = 176.8*D2R;

	if (mi > mi_s && mi < 2.0*G_PI - mi_s && fabs(beta) < beta0) {
		// noon maneuver
		_noon_turn_GLOM(prn, beta, mi, xsat, vsat, xsun, R, i, j, k);
	}
	else if (fabs(mi) < EPS0_GLO && fabs(beta) < EPS0_GLO) {
		// midnight maneuver
		_midnight_turn_GLOM(prn, beta, mi, xsat, vsat, xsun, R, i, j, k);
	}
	else {
		_ysm(prn, beta, mi, xsat, vsat, xsun, i, j, k);
	}
}


// Calculate satellite-fixed vectors from yaw angle    
void t_gpppmodel::_yaw2ijk(t_gsatdata& satdata, double& yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
  satdata.yaw(yaw);   // store yaw angle

  ColumnVector satcrd = satdata.satcrd().crd_cvect();
  ColumnVector satvel = satdata.satvel().crd_cvect();       

  if(satcrd.NormFrobenius() == 0 || satvel.NormFrobenius() == 0){ return; }
    
  // ITRF -> ICRF velocity
  satvel(1) -= OMEGA*satcrd(2);
  satvel(2) += OMEGA*satcrd(1);
   
  ColumnVector n = crossproduct(satcrd, satvel);   
    
  // Satelite-Earth unit vector
  k = -satcrd;
  k /= k.NormFrobenius();    

  ColumnVector ex = crossproduct(n, satcrd);
  
  ex /= ex.NormFrobenius();
  n  /=  n.NormFrobenius();
    
       
  double cosy = cos(yaw);
  double siny = sin(yaw);
  for (int r = 1; r <= 3; r++) {
    i(r) = -siny*n(r) + cosy*ex(r);
    j(r) = -cosy*n(r) - siny*ex(r);
  }     
}

void t_gpppmodel::_yaw2ijk(ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, double& yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
	//satdata.yaw(yaw);   // store yaw angle
	ColumnVector satcrd = xsat;
	ColumnVector satvel = vsat;

	if (satcrd.NormFrobenius() == 0 || satvel.NormFrobenius() == 0) { return; }

	// ITRF -> ICRF velocity
	satvel(1) -= OMEGA*satcrd(2);
	satvel(2) += OMEGA*satcrd(1);

	ColumnVector n = crossproduct(satcrd, satvel);

	// Satelite-Earth unit vector
	k = -satcrd;
	k /= k.NormFrobenius();

	ColumnVector ex = crossproduct(n, satcrd);

	ex /= ex.NormFrobenius();
	n /= n.NormFrobenius();


	double cosy = cos(yaw);
	double siny = sin(yaw);
	for (int r = 1; r <= 3; r++) {
		i(r) = -siny*n(r) + cosy*ex(r);
		j(r) = -cosy*n(r) - siny*ex(r);
	}
}
   
// sign function: return the value of a with the sign of b; added by yqyuan
double t_gpppmodel::sign(double a, double b) {
	double value = fabs(a);
	int sign_of_b = 1;
	if (b >= 0) { ; }
	else { sign_of_b = -1; }
	return value*sign_of_b;
}
    
} // namespace
