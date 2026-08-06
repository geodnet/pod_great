# GREAT-PODFLT: GNSS Precise Orbit Determination Software by Wuhan University GREAT Group (Version 1.0)

## Overview

The GREAT (GNSS+ REsearch, Application and Teaching) software suite is designed and developed by the School of Geodesy and Geomatics, Wuhan University. It is a comprehensive platform for space geodesy data processing, precise positioning and orbit determination, as well as multi-source fusion navigation.

GREAT-PODFLT is a key module within the GREAT software, dedicated to GNSS Precise Orbit Determination (POD) using the Square Root Information Filter (SRIF). It processes ground-based GNSS observations together with navigation satellite dynamics models to perform batch estimation of precise satellite orbits. The software supports multi-constellation POD including GPS, GLONASS, Galileo, and BDS-2/3, with capabilities such as dual-frequency ionosphere-free combination processing, integer ambiguity resolution, Earth Orientation Parameter (EOP) estimation, and simultaneous receiver and satellite clock estimation. Core computation modules are implemented in C++ (C++14/17), while auxiliary Python 3 scripts are provided for plotting results. GREAT-PODFLT uses CMake for build management, allowing users to flexibly choose mainstream C++ compilers such as GCC, Clang, and MSVC. It supports build and execution on both Windows and Linux.

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
  ├── util/                      Utility tools
  │   ├── orbdif/                Orbit comparison tool (Windows & Linux versions)
  │   └── gnss_sysfile/          GNSS common system files (planetary ephemeris, ANTEX, gravity models, etc.)
  └── doc/                       Documentation related to GREAT-PODFLT *
      └── GREAT_PODFLT.pdf       GREAT-PODFLT Documentation
```

## Installation and Usage

See **GREAT_PODFLT.pdf** included in the `./doc` directory, or follow our upcoming tutorial videos on bilibili.

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
