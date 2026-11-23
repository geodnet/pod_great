# GREAT-PODFLT: GNSS Precise Orbit Determination Software by Wuhan University GREAT Group (Version 1.0)

## Overview

The GREAT (GNSS+ REsearch, Application and Teaching) software suite is designed and developed by the School of Geodesy and Geomatics, Wuhan University. It is a comprehensive platform for space geodesy data processing, precise positioning and orbit determination, as well as multi-source fusion navigation.

GREAT-PODFLT is a key module within the GREAT software, dedicated to GNSS Precise Orbit Determination (POD) using the Square Root Information Filter (SRIF). It processes ground-based GNSS observations together with navigation satellite dynamics models to perform batch/sequential estimation of precise satellite orbits. The software supports multi-constellation POD including GPS, GLONASS, Galileo, and BDS-2/3, with capabilities such as dual-frequency ionosphere-free combination processing, integer ambiguity resolution, Earth Orientation Parameter (EOP) estimation, and simultaneous receiver and satellite clock estimation. Core computation modules are implemented in C++ (C++14), using CMake for cross-platform builds on Linux, macOS, and Windows.

GREAT-PODFLT consists of two portable program libraries: LibGREAT and LibGNut. LibGREAT encapsulates the core precise orbit determination algorithms, including perturbation force models (non-spherical gravitation, N-body perturbation, solid Earth and ocean tides, ECOM2 solar radiation pressure, Earth albedo radiation, antenna thrust), Adams predictor-corrector numerical integrator, SRIF time/measurement update, and integer ambiguity resolution (LAMBDA method). LibGNut (derived from the G-Nut/Anubis library) provides fundamental GNSS data processing capabilities such as RINEX observation and navigation file reading, precise ephemeris and clock product parsing, antenna phase center correction, tropospheric and ionospheric models, and time/coordinate system transformations.

This open-sourced GREAT-PODFLT version 1.0 supports the following capabilities:

1. Supports multi-GNSS precise orbit determination for GPS, Galileo, BDS-2/3, and GLONASS based on Square Root Information Filtering (SRIF)

2. Supports integer ambiguity resolution

3. Supports dual-frequency ionosphere-free combination observations for both pseudorange and carrier phase measurements

4. Supports comprehensive perturbation force models including non-spherical gravitation, N-body perturbation, solid Earth and ocean tides, ECOM2 solar radiation pressure,and antenna thrust

5. Includes the orbit comparison tool (orbdif) for accuracy assessment, outputting orbit residuals in ACR (Along-track / Cross-track / Radial) coordinates against reference precise ephemerides

## Package Directory Structure

```
GREAT-PODFLT
  ├── .gitignore                 Ignores build artifacts and private run outputs
  ├── config/
  │   └── podflt_geodnet.xml.in  Template configuration for GEODNET NRT POD runs
  ├── scripts/
  │   ├── geodnet_prep.py        Prepare RTCM3/RINEX, SINEX, ERP, OSB, OTL BLQ, and podflt.xml
  │   ├── make_blq.py            Synthesize FES2004 CM ocean tide loading (.blq) for stations
  │   ├── make_satinfo.py        Update satinfo from IGS satellite metadata SINEX
  │   ├── great_product.py       Package filter orbits/clocks into IGS-named SP3 & CLK files
  │   └── hooks/pre-commit       Guard against committing IDs, coordinates, secrets, or local paths
  ├── src/                       Source code *
  │   ├── CMakeLists.txt         Top-level CMake build script
  │   ├── app/                   Main program —— GREAT-PODFLT POD filter *
  │   ├── libGNut/               G-Nut base library (includes built-in Newmat & pugixml) *
  │   └── libGREAT/              GREAT core library (perturbation models, SRIF, ambiguity fixing) *
  ├── third-party/               Third-party libraries *
  │   ├── Eigen/                 Linear algebra template library (header-only)
  │   └── zlib/                  Data compression library
  ├── sample/                    Example datasets *
  │   ├── podflt_ambfix_G_2021210/   GPS example
  │   ├── podflt_ambfix_E_2021210/   Galileo example
  │   ├── podflt_ambfix_C_2021210/   BDS example
  │   └── podflt_ambfix_R_2021210/   GLONASS example
  ├── utils/                     Utility tools
  │   ├── orbdif/                Orbit comparison tool (Windows & Linux versions)
  │   └── gnss_sys_file/         GNSS common system files (planetary ephemeris, ANTEX, gravity models, etc.)
  └── doc/                       Documentation related to GREAT-PODFLT *
      └── GREAT_PODFLT.pdf       GREAT-PODFLT Documentation
```

## Installation and Usage

### Build (Linux / macOS)

```sh
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

Pass `-DUSE_OPENMP=ON` on platforms with OpenMP available. See **doc/GREAT_PODFLT.pdf** for core algorithm details and sample datasets.

### GEODNET Near-Real-Time (NRT) SP3 & Clock Workflow

1. **Refresh `satinfo` (recommended for current constellations)** using the public IGS metadata SINEX (`igs_satellite_metadata.snx`):

   ```sh
   python3 scripts/make_satinfo.py igs_satellite_metadata.snx utils/gnss_sys_file/satinfo /path/to/work/satinfo
   ```

2. **Prepare an external working directory** (must stay outside this repository; private station IDs, ECEF coordinates, `.snx`, and `.blq` files are kept out of git):

   ```sh
   python3 scripts/geodnet_prep.py \
     --stations /path/to/stations.csv \
     --beg 2026-11-18T00 --end 2026-11-22T00 \
     --erp /path/to/bula.erp \
     --work /path/to/work \
     --rtcm-root /path/to/rtcm \
     --convbin /path/to/convbin \
     --atx /path/to/igs20.atx \
     --satinfo /path/to/work/satinfo \
     --bias /path/to/osb.bia \
     --otl-models /path/to/otl_models \
     --sys "GPS GAL" --int 300 --threads 8
   ```

   Key design conventions handled automatically by `geodnet_prep.py`:
   - **APC coordinates (no receiver PCO/PCV or H/E/N offset)**: GEODNET coordinates are referred to the GPS L1+L2 ionosphere-free antenna phase center (APC). Converted RINEX headers have `ANT # / TYPE` cleared (`ADVNULLANTENNA`) and `ANTENNA: DELTA H/E/N` zeroed (`0.0000 0.0000 0.0000`), overriding any RTCM-advertised receiver antenna metadata.
   - **Anonymized site aliases**: 12-character station IDs are mapped to 4-character aliases (`G000`...) in `<work>/alias_map.csv`.
   - **Ocean tide loading (`--otl-models`)**: Uses `scripts/make_blq.py` (FES2004 CM harmonic synthesis) to build `<work>/geodnet.blq`.
   - **Integrator/step consistency**: `--int` populates `@INT@`, `@STEP_OUT@`, inclusive `@END@` (`end - int`), and a dividing `@STEP_IN@` in `config/podflt_geodnet.xml.in`.

3. **Run `great_podflt`** from the generated run directory (broadcast-nav initialization enables automatically when no input SP3 is configured):

   ```sh
   cd /path/to/work/run_<tag>
   /path/to/build/Bin/great_podflt -x podflt.xml
   ```

4. **Package NRT SP3 & Clock products** (uses converged windows from `result/orbit.sp3` and `result/satclk.clk`):

   ```sh
   python3 scripts/great_product.py /path/to/work/run_<tag> --start 2026-11-20T00 --hours 24 --orb-step 300 --clk-step 300
   ```

   Outputs follow IGS filenames (`<AC>0<CMP><TYP>_<YYYYJJJHHMM>_<SPAN>_<SMP>_ORB.SP3` and `_CLK.CLK`, defaulting to `GRT0OPSNRT`). Note that real-time RTCM SSR streaming output (RT orbit/clock/phase-bias encoder) is not yet part of `great_podflt`.

### Public IGS Validation Benchmarks (2026-11-18 – 2026-11-21)

Evaluated with `IGS20` antenna phase centers and IERS Bulletin A EOPs against IGS Analysis Center public products (`IGS0OPSRAP`, `COD0OPSRAP`, `GFZ0MGXRAP`, `GFZ0OPSRAP`, `WUM0MGXNRT`):

- **100-station GPS (`300 s`, 5-day arc with FES2004 `.blq`, converged days 2–4)**:
  - **DOY 323**: `92.8 mm` 3D orbit RMS (`0.172 ns` clock STD, `43.3 mm` SISRE) vs `IGS0OPSRAP`; `92.9 mm` vs `GFZ0MGXRAP`; `94.5 mm` vs `WUM0MGXNRT`.
  - **DOY 324**: `74.1 mm` 3D (`0.164 ns` STD, `42.8 mm` SISRE) vs `IGS0OPSRAP`; `75.6 mm` vs `COD0OPSRAP`; `73.0 mm` vs `GFZ0MGXRAP`; `75.0 mm` vs `WUM0MGXNRT`.
  - **DOY 325**: `81.0 mm` 3D (`0.171 ns` STD, `48.3 mm` SISRE) vs `IGS0OPSRAP`; `84.0 mm` vs `COD0OPSRAP`; `82.1 mm` vs `GFZ0OPSRAP`.
- **100-station GPS (`30 s`, 3-day arc, without per-station `.blq`, DOY 324)**: `63.2 mm` 3D (`0.126 ns` STD, `29.8 mm` SISRE) vs `IGS0OPSRAP`; `65.1 mm` vs `COD0OPSRAP`; `62.9 mm` vs `GFZ0MGXRAP` (`29.2 mm` SISRE); `64.8 mm` vs `WUM0MGXNRT`.
- **40-station multi-GNSS `GPS+GAL` (`300 s`, DOY 324 18–24 h)**:
  - GPS (`32` sats): `198.1–199.9 mm` 3D (`0.206–0.212 ns` STD, `42.1–43.4 mm` SISRE) across `IGS0OPSRAP`, `COD0OPSRAP`, `GFZ0MGXRAP`, and `WUM0MGXNRT`.
  - GAL (`29–30` sats): `199.5–202.0 mm` 3D (`0.298–0.305 ns` STD, `52.1–52.5 mm` SISRE) across `COD0OPSRAP`, `GFZ0MGXRAP`, and `WUM0MGXNRT`.

### Public Repository Safety

Never store station IDs, coordinates, alias maps, SINEX/BLQ outputs, credentials, or local filesystem paths inside the repository. Install the included pre-commit hook before committing:

```sh
ln -s ../../scripts/hooks/pre-commit .git/hooks/pre-commit
```

## Contributing

### Developers

Wuhan University GREAT Team, School of Geodesy and Geomatics, Wuhan University.

### Third-party Libraries

GREAT-PODFLT uses the G-Nut library ([http://www.pecny.cz](http://www.pecny.cz)) Copyright (C) 2011–2016 GOP - Geodetic Observatory Pecny, RIGTC.

GREAT-PODFLT uses the pugixml library ([http://pugixml.org](http://pugixml.org)) Copyright (C) 2006–2014 Arseny Kapoulkine.

GREAT-PODFLT uses the Newmat library ([http://www.robertnz.net/nm_intro.htm](http://www.robertnz.net/nm_intro.htm)) Copyright (C) 2008: R B Davies.

GREAT-PODFLT uses the Eigen library ([https://eigen.tuxfamily.org](https://eigen.tuxfamily.org)) Copyright (C) 2008–2011 Gael Guennebaud & Eigen contributors.

GREAT-PODFLT uses the zlib library ([https://www.zlib.net](https://www.zlib.net)) Copyright (C) 1995–2017 Jean-loup Gailly & Mark Adler.

## Download

GitHub: [https://github.com/GREAT-WHU/GREAT_PODFLT](https://github.com/GREAT-WHU/GREAT_PODFLT)

## Others

You are welcome to join the QQ group (1009827379) for discussion and exchange.

WeChat Official Account: **GREAT智能导航实验室** — we will continue to share team updates.

bilibili Account: **GREAT智能导航实验室** — we will continue to publish software walkthrough videos.
