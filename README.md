[![GitHub (pre-)release](https://img.shields.io/github/release/vasilsaroka/TBpack/all.svg)](https://github.com/vasilsaroka/TBpack/releases)
[![Github All Releases](https://img.shields.io/github/downloads/vasilsaroka/TBpack/total.svg)](https://github.com/vasilsaroka/TBpack/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Support TBpack](https://img.shields.io/static/v1?label=support&message=5$&color=green&style=flat&logo=paypal)](https://paypal.me/vasilsaroka?locale.x=en_GB)

# TBpack
Tight-binding calculations in Mathematica

## Installation

 - [Download the latest release](https://github.com/vasilsaroka/TBpack/releases), distributed as a `.paclet` file, and install it using the `PacletInstall` function in Mathematica:

        Needs["PacletManager`"]
        PacletInstall["path2paclet"]
        
   Insert `path2paclet` via Mathematica's Insert → File Path... menu command.

 - Evaluate ``<<TBpack` `` to load the application into the Mathematica session.
 - Test TBpack using `AtomicStructure[Nanotube[10, 10]]`.
 - Test MaTeX using `MaTeX["x^2"]`.
 
   Packages [MaTeX](https://github.com/szhorvat/MaTeX/releases) by Szabolcs Horvát and [CustomTicks](https://library.wolfram.com/infocenter/Demos/5599/) by Mark A. Caprio  are integral parts of `TBpack`.

   Compilation to C code is used in TBpack to speed up optical absorption spectra calculations. See [how to make Mathematica working with C compiler on Windows](https://sites.google.com/site/sarokavasil/wolfram-mathematica). If compiler is not available Mathematica will run uncompiled function.

 - Open the documentation center and search for "Hamiltonian" to get started.

## Supporting the project
   TBpack is free to use and distribute and will stay so. However, since it is developed in my own time consider supporting the project to help me to explain my wife why I should spend our weekends on this :)
