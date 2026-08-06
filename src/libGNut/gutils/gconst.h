
#ifndef CONST_H
#define CONST_H 

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: defines constants
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include <map>
#include <string>
#include <vector>
#include "ExportHeaders/ExportBase.h"

using namespace std;

namespace gnut
{
#define G_PI        3.14159265358979311599796346854419e0    // pi
#define D2R         (G_PI/180.0)          // deg to rad
#define R2D         (180.0/G_PI)          // rad to deg
#define SEC2RAD     (D2R/3600.0)          // sec to rad
#define RAD2SEC     (R2D*3600.0)          // rad to sec
#define RAD2TSEC	(RAD2SEC/15.0)			
#define RHO_SEC     3600.0*D2R
#define CLIGHT      2.99792458e+8         // speed of light [m/s]
#define CLIGHT_2    8.9875517873681764e16 // pow(CLIGTH,2)
#define OMEGA       7292115.1467e-11
#define Aell        6378137.000           // semi major axes [m]
#define Finv        298.257223563         // inverse flattening [-]
#define MJD_J2000   51544.5               // J2000 MJD date
#define R_SPHERE    6371000.000           // [m] approx. sphere radius
#define POST_SHADOW 1800.0                // satellite post shadow period [s]

// GPS   
#define A_WGS       6378137.000           // [m] WGS84 semi-major axis
#define B_WGS       6356752.300           // [m] WGS84 semi-minor axis
#define E_WGS       0.081819              // [-] WGS84 eccentricity
#define F_WGS       0.003352811           // [-] WGS84 flatenning
#define MUDOT_GPS   (0.00836*D2R)         // avarage GPS satellite angular velocity [rad]
#define EPS0_GPS    (13.5*D2R)            // maximal GPS satellites crossing angle [deg]

// GLONASS
#define GM_PZ90     398.60044e12          // [] PZ90 earth's graviational constant
#define Aell_PZ90   6378136.000           // [m] PZ90 semi-major axes
#define C20_PZ90   -1082.62575e-6         // [] PZ90
#define MUDOT_GLO   (0.00888*D2R)         // avarage GLO satellite angular velocity [rad]   
#define EPS0_GLO    (14.2*D2R)            // maximal GLO satellites crossing angle [deg]	 
#define Finv_GLO    298.25784             // inverse flattening [-]
#define OMGE_DOT_GLO 7.292115e-5          // Mean angular velocity of the Earth [rad/sec]

// Galileo
#define GM_GAL       3.986004418e14      // Geocentric gravitational constant [m^3/s^2]
#define OMGE_DOT_GAL 7.2921151467e-5     // Mean angular velocity of the Earth [rad/sec]

// BeiDou
#define GM_CGCS      3.986004418e14      // []  CGCS2000 earth's graviational constant
#define OMGE_DOT_BDS 7.2921150e-5        // BDS value of the earth's rotation rate [rad/sec]// 
#define Aell_CGCS    6378137.000         // [m] CGCS2000 semi-major axes

// METEO
#define G_RATIO     0.003449787           // [-] gravity ratio
#define G_EQUA      9.7803253359          // [m/s2] equatorial gravity
#define G_POLE      9.8321849378          // [m/s2] polar gravity
#define G_WMO       9.80665               // [m/s2] gravity acceleration (WMO constant for lat=45)
#define KS_SOM      0.001931853           // [-] Somigliana's constant
#define Md          28.9644               // [g/mol] molar mass (mean molecular weight) of dry air
#define Mw          18.01528              // [g/mol] molar mass (mean molecular weight) of water
#define Ru          8.3144621             // [J/mol/K] universal gas constant
#define Rd          287.058               // [J/kg/K] spec. gas constant for Dry Air (Rd = 287.058)  Ru/(Md/1000)
#define Rv          461.495               // [J/kg/K] spec. gas constant for Wat Vap (Rv = 461.495)  Ru/(Mw/1000)
#define Eps         0.62198               // [-] molecular weight ratio (Mw/Md)     (Eps = 0.62198)
#define TC2TK       273.15                // conversion from deg of Celsius to Kelvins
#define TPOINT      273.16                // temperature of triple point [K], i.e. 0.01 Deg C

	typedef map<string, vector<double> > t_map_refr; // refractivity coefficients
	t_map_refr REFR_COEF();                          // static map of refractivity coefficients

#define K1_ESS      77.64                 // [K/hPa]   Essen and Froome (1951)
#define K2_ESS      64.68                 // [K/hPa]   Essen and Froome (1951)
#define K3_ESS      3.718e5               // [K^2/hPa] Essen and Froome (1951)

#define K1_SMW      77.607                // [K/hPa]   Smith and Weintraub (1953) // FULL RESOLUTION
#define K2_SMW      71.6                  // [K/hPa]   Smith and Weintraub (1953) // original paper 3-par formula!
#define K3_SMW      3.747e5               // [K^2/hPa] Smith and Weintraub (1953) // original paper 3-par formula!

#define K1_THA      77.604                // [K/hPa]   Thayer (1974)
#define K2_THA      64.79                 // [K/hPa]   Thayer (1974)
#define K3_THA      3.776e5               // [K^2/hPa] Thayer (1974)

#define K1_BEV      77.6                  // [K/hPa]   Bevis (1994)
#define K2_BEV      70.4                  // [K/hPa]   Bevis (1994)
#define K3_BEV      3.739e5               // [K^2/hPa] Bevis (1994)

#define K1_FOE      77.65                 // [K/hPa]   Foelsche (1999)
#define K2_FOE      65.99                 // [K/hPa]   Foelsche (1999)
#define K3_FOE      3.777e5               // [K^2/hPa] Foelsche (1999)

#define K1_RUE      77.6848               // [K/hPa]   Rueger (2002)
#define K2_RUE      71.2952               // [K/hPa]   Rueger (2002)
#define K3_RUE      3.75463e5             // [K^2/hPa] Rueger (2002)



// for calulate the ocean tide and else
#ifndef EARTH_GM
   #define EARTH_GM			398600.4415		///< GM of Earth
#endif // !EARTH_GM

#ifndef EARTH_R
   #define EARTH_R				6378.13646		///< radius of Earth(unit: km)
#endif

#ifndef SUN_GM
   #define SUN_GM				1.32712440e+11	///< GM of Sun
#endif // !SUN_GM

#ifndef SUN_R
	#define SUN_R				696000.0		///< radius of Sun(unit: km)
#endif // !SUN_R

#ifndef MOON_GM
	#define MOON_GM				4902.7989		///< GM of Moon
#endif // !MOON_GM

#ifndef MOON_R
	#define MOON_R              1738.090		///< radius of Moon(unit: km)
#endif // !MOON_R

//confidence 0.85
const double _norm[] ={
            0,2.072251e+00, 3.794240e+00, 5.317048e+00, 6.744883e+00, 8.115199e+00, 9.446103e+00, 1.074790e+01, 1.202707e+01, 1.328804e+01, 1.453394e+01,
            1.576710e+01, 1.698931e+01, 1.820198e+01, 1.940624e+01, 2.060301e+01, 2.179306e+01, 2.297703e+01, 2.415547e+01, 2.532885e+01, 2.649758e+01,
    };


//confidence 0.999
const double _chi2inv[] =
    {
        0,1.082757e+01, 1.381551e+01, 1.626624e+01, 1.846683e+01, 2.051501e+01, 2.245774e+01, 2.432189e+01, 2.612448e+01, 2.787716e+01, 2.958830e+01,
        3.126413e+01, 3.290949e+01, 3.452818e+01, 3.612327e+01, 3.769730e+01, 3.925235e+01, 4.079022e+01, 4.231240e+01, 4.382020e+01, 4.531475e+01,
        4.679704e+01, 4.826794e+01, 4.972823e+01, 5.117860e+01, 5.261966e+01, 5.405196e+01, 5.547602e+01, 5.689229e+01, 5.830117e+01, 5.970306e+01,
        6.109831e+01, 6.248722e+01, 6.387010e+01, 6.524722e+01, 6.661883e+01, 6.798517e+01, 6.934645e+01, 7.070289e+01, 7.205466e+01, 7.340196e+01,
        7.474494e+01, 7.608376e+01, 7.741858e+01, 7.874952e+01, 8.007673e+01, 8.140033e+01, 8.272042e+01, 8.403713e+01, 8.535056e+01, 8.666082e+01,
        8.796798e+01, 8.927215e+01, 9.057341e+01, 9.187185e+01, 9.316753e+01, 9.446054e+01, 9.575095e+01, 9.703883e+01, 9.832423e+01, 9.960723e+01,
        1.008879e+02, 1.021662e+02, 1.034424e+02, 1.047163e+02, 1.059881e+02, 1.072579e+02, 1.085256e+02, 1.097913e+02, 1.110551e+02, 1.123169e+02,
        1.135769e+02, 1.148351e+02, 1.160915e+02, 1.173462e+02, 1.185991e+02, 1.198503e+02, 1.211000e+02, 1.223480e+02, 1.235944e+02, 1.248392e+02,
        1.260826e+02, 1.273244e+02, 1.285648e+02, 1.298037e+02, 1.310412e+02, 1.322773e+02, 1.335121e+02, 1.347455e+02, 1.359776e+02, 1.372084e+02,
        1.384379e+02, 1.396661e+02, 1.408931e+02, 1.421189e+02, 1.433435e+02, 1.445670e+02, 1.457892e+02, 1.470104e+02, 1.482304e+02, 1.494493e+02,
        1.506671e+02, 1.518838e+02, 1.530995e+02, 1.543141e+02, 1.555277e+02, 1.567403e+02, 1.579518e+02, 1.591624e+02, 1.603721e+02, 1.615807e+02,
        1.627885e+02, 1.639953e+02, 1.652011e+02, 1.664061e+02, 1.676102e+02, 1.688133e+02, 1.700156e+02, 1.712171e+02, 1.724177e+02, 1.736174e+02,
        1.748164e+02, 1.760145e+02, 1.772118e+02, 1.784083e+02, 1.796040e+02, 1.807989e+02, 1.819930e+02, 1.831864e+02, 1.843791e+02, 1.855710e+02,
        1.867621e+02, 1.879526e+02, 1.891423e+02, 1.903313e+02, 1.915196e+02, 1.927072e+02, 1.938941e+02, 1.950803e+02, 1.962659e+02, 1.974508e+02,
        1.986350e+02, 1.998186e+02, 2.010015e+02, 2.021838e+02, 2.033655e+02, 2.045465e+02, 2.057270e+02, 2.069068e+02, 2.080860e+02, 2.092646e+02,
        2.104426e+02, 2.116200e+02, 2.127969e+02, 2.139732e+02, 2.151489e+02, 2.163240e+02, 2.174986e+02, 2.186726e+02, 2.198460e+02, 2.210190e+02,
        2.221914e+02, 2.233632e+02, 2.245345e+02, 2.257053e+02, 2.268756e+02, 2.280454e+02, 2.292146e+02, 2.303834e+02, 2.315516e+02, 2.327194e+02,
        2.338866e+02, 2.350534e+02, 2.362197e+02, 2.373855e+02, 2.385508e+02, 2.397157e+02, 2.408801e+02, 2.420440e+02, 2.432075e+02, 2.443705e+02,
        2.455330e+02, 2.466951e+02, 2.478568e+02, 2.490180e+02, 2.501788e+02, 2.513392e+02, 2.524991e+02, 2.536586e+02, 2.548177e+02, 2.559763e+02,
        2.571346e+02, 2.582924e+02, 2.594498e+02, 2.606068e+02, 2.617634e+02, 2.629197e+02, 2.640755e+02, 2.652309e+02, 2.663859e+02, 2.675405e+02,

   
        };


} // namespace

#endif

