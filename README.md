[![GitHub (pre-)release](https://img.shields.io/github/release/vasilsaroka/TBpack/all.svg)](https://github.com/vasilsaroka/TBpack/releases)
[![Github All Releases](https://img.shields.io/github/downloads/vasilsaroka/TBpack/total.svg)](https://github.com/vasilsaroka/TBpack/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Support TBpack](https://img.shields.io/static/v1?label=support&message=5$&color=green&style=flat&logo=paypal)](https://paypal.me/vasilsaroka?locale.x=en_GB)

# TBpack
Tight-binding calculations in Mathematica

## TBpack for fun
[<img align="right" src="https://github.com/vasilsaroka/TBpack/blob/master/Games/HOMOLUMOQuantumGame/HOMOLUMOQuantumGame.png" alt="HOMOLUMOQuantumGame" width="250"/>](https://t.me/HOMO_LUMO_Quantum_Game_bot)

Before putting your hands on **TBpack** application you can get its flavor and spirit by playing a couple of rounds with our Telegram bot. Put your knowledge and luck under ultimate test with Hostess Monica in HOMO-LUMO Quantum Game: [@HOMO_LUMO_Quantum_Game_bot](https://t.me/HOMO_LUMO_Quantum_Game_bot). Use `/rules` to see the rules of the game or have a look at our sample [game video](https://www.youtube.com/watch?v=DUy22OKLgvs) on YouTube.


## Hands-on Trial
Try our Cloud based demo -- **TBpackDemo**. **TBpackDemo** contains all fully functional core functions, excluding visualization and data currating functions. [<img align="right" src="https://github.com/vasilsaroka/TBpack/blob/master/TBpack/TBpackDemo_Logo.png" alt="TBpackDemo_Logo"/>](https://paypal.me/vasilsaroka?locale.x=en_GB)

   
   Get **TBpackDemo** immediately in your notebook from Wolfram Cloud by evaluating in Mathematica 10.0+:
   
       CloudGet["https://www.wolframcloud.com/obj/vasil.saroka/TBpack/Demo/TBpackDemo.wl"]
       
   **TBpackDemo** does not have documentation built-in into Wolfram Documentation Center, only user information messages can be invoked for the core fucntions. The available functions are `Hamiltonian`, `ElectronicStructure`, `ElectronicBands1D`. To get informaiton on the usage of `Hamiltonian` core fucntion, try
   
       ?Hamiltonian
   
   For a flawless start, you can use our Cloud based demo examples. Download [ExamplesDemo.nb](https://www.wolframcloud.com/obj/vasil.saroka/TBpack/Demo/ExamplesDemo.nb) with examples and follow the instructions inside.
       
       
       

## Installation guide
 - **The recommended option for TBpack paclet management in Mathematica 11.3+** is to use the [resource function](https://resources.wolframcloud.com/FunctionRepository/resources/InstallTBpack): [<img align="right" src="https://github.com/vasilsaroka/TBpack/blob/master/TBpack/TBpackPro_Logo.png" alt="TBpackPro_Logo"/>](https://paypal.me/vasilsaroka?locale.x=en_GB)
 
   To install the latest version of the TBpack from this GitHub repo, evaluate
   
       ResourceFunction["InstallTBpack"][]

   To find all installed versions on your PC and to uninstall chosen (all) version (s), evaluate
   
       ResourceFunction["InstallTBpack"][Method->"Uninstall"]
   
 - The automated management of TBpack paclet in Mathematica 10.0+:
 
   To install the latest version of the TBpack from this GitHub repo, evaluate
   
        CloudGet["https://www.wolframcloud.com/obj/vasil.saroka/TBpack/Services/InstallTBpack"];
        InstallTBpack[]
        
   To find all installed versions on your PC and to uninstall chosen (all) version (s), evaluate
   
        CloudGet["https://www.wolframcloud.com/obj/vasil.saroka/TBpack/Services/InstallTBpack"];
        InstallTBpack[Method->"Uninstall"]
        
 - The automated installation option for Mathematica 10.0+:
 
   Copy-paste the below function into a Mathematica notebook cell. Evaluate the cell to make the definition of this function known to Mathematica. In the next cell type and evaluate `InstallTBpack[]`.
   
        InstallTBpack[] := Block[{jsonreleases, info, url, message,tempfile,fpath},
        jsonreleases = Import["https://api.github.com/repos/vasilsaroka/TBpack/releases","JSON"];
        If[jsonreleases === $Failed, Return[$Failed]];
        info = "Downloading TBpack " <> First@Lookup[jsonreleases, "tag_name"];
        If[
           $Notebooks,
           PrintTemporary@Row[{info,ProgressIndicator[Appearance -> "Percolate"]},Frame -> True,RoundingRadius -> 9], 
           Print[info <> "..."]
        ];
        url = First@Lookup[First[Lookup[jsonreleases, "assets"]],"browser_download_url"];
        message = "TBpack is succefully installed.";
        If[
   	         $VersionNumber >= 12.1,
   	         Check[PacletInstall[url, ForceVersionInstall -> True],Return[$Failed]];
   	         Print[message],
   	         If[
    		   $VersionNumber >= 11.0,
                  tempfile = FileNameJoin[{$TemporaryDirectory,FileNameTake[url]}];
                  Check[fpath = URLDownload[url, tempfile],Return[$Failed]];
                  If[
                      FileExistsQ[fpath],
                      Check[PacletManager\`PacletInstall[fpath,"IgnoreVersion" -> True],Return[$Failed]];
                      DeleteFile[fpath];
                      Print[message],
                      $Failed
     		   ],
    		   If[
                      $VersionNumber >= 10.0,
                      tempfile = FileNameJoin[{$TemporaryDirectory,FileNameTake[url]}];
                      Check[fpath = URLSave[url,tempfile],Return[$Failed]];
                      If[
                           FileExistsQ[fpath],
                           Check[PacletManager\`PacletInstall[fpath,"IgnoreVersion" -> True],Return[$Failed]];
                           DeleteFile[fpath];
                           Print[message],
                           $Failed
                      ],
                      $Failed
                 ](* end If *)
              ](* end If *)
          ](* end If  *)
        ](* end Block *) 

 - An alternative manual intallation option for Mathematica 10.0+ is the following. [Download the latest release](https://github.com/vasilsaroka/TBpack/releases), distributed as a `.paclet` file, and install it using the `PacletInstall` function in Mathematica:

        Needs["PacletManager`"]
        PacletInstall["path2paclet"]
        
   Insert `path2paclet` via Mathematica's Insert → File Path... menu command.
   
 - In Mathematica 12.1+ ``PacletInstall["url"]`` can be used for the installation straight away:
        
        PacletInstall["https://github.com/vasilsaroka/TBpack/releases/download/v<version>/TBpack-<version>.paclet"]  
   where `<version>` stands for any existing version of the application. For example,
   
       PacletInstall["https://github.com/vasilsaroka/TBpack/releases/download/v0.2.0/TBpack-0.2.0.paclet"]

   
## Demo
 - Evaluate ``<<TBpack` `` to load the application into the Mathematica session.
 - Test TBpack using `AtomicStructure[Nanotube[10, 10]]`.
 - Test MaTeX using `MaTeX["x^2"]`.
 - Test CustomTicks using `Plot[Sin[x], {x, 0, 3}, Ticks -> {LinTicks[-1, 3, 1, 5], LinTicks[0, 1, 1, 5]}]`.
 - Test MaTeX and CustomTicks altogether using 
     
       Plot[Sin[x], {x, -1, 3}, 
            Axes -> {MaTeX["x"],MaTeX["\sin(x)"]},
            AxesStyle -> BlackFrame, 
            Ticks -> {
            LinTicks[-1, 3, 1, 5, MinorTickLength -> 0.015, MajorTickLength -> 0.025], 
            LinTicks[-1, 1, 1, 5, MinorTickLength -> 0.015, MajorTickLength -> 0.025]
            }
       ]
 
   Packages [MaTeX 1.7.8](https://github.com/szhorvat/MaTeX/releases) by Szabolcs Horvát and [CustomTicks 2.1.0](https://library.wolfram.com/infocenter/Demos/5599/) by Mark A. Caprio are integral parts of TBpack.
   
    **Note:** Make sure that [a TeX system](https://tug.org/begin.html) and Ghostscript 9.15 or later are installed so that MaTeX can work. The latest version of [Ghostscript](https://www.ghostscript.com/download/gsdnld.html) can be downloaded for Windows and Linux. On OS X, MacTeX 2015 and later already include a compatible version of Ghostscript. When MaTeX is not able to locate these installations automatically, you must provide paths to pdfLaTeX and Ghostscript executable files using `ConfigureMaTeX["pdfLaTeX"->"path2pdflatex","Ghostscript"->"path2gsfile"]`.
   
 - Test Sneg using

         snegfermionoperators[c];
         nc[c[AN, 1, UP], VACUUM] (* annihilates the vacuum producing zero result *)
         nc[c[CR, 1, UP], VACUUM] (* creates the state with fermion having spin-1/2 with +1/2 projection *)
     
   [Sneg 1.250](http://nrgljubljana.ijs.si/sneg/) library by Rok Zitko is an integral but independent part of TBpack starting from v0.4.0. Note that Sneg can be used as a stand-alone package.

 - Open the documentation center and search for `Hamiltonian` to get started.
 
 **Note:** Compilation to C code is used in TBpack to speed up optical absorption spectra calculations. See [how to make Mathematica working with C compiler on Windows](https://sites.google.com/view/vasilsaroka/wolfram-mathematica). If compiler is not available Mathematica will run uncompiled function.
 

## Supporting the project
   We believe everyone deserves access to knowledge that is grounded in science and integrity. That is why we keep our code open for all users, regardless of where they live or what they can afford to pay. This means more people can be better educated and inspired to make an impact on the global wellbeing. We have no shareholders or billionaire owner, meaning only your donations power our work and ensure it can remain open for all. Every contribution, however big or small, makes a real difference for TBpack future. If you find it useful, consider [supporting the project](https://paypal.me/vasilsaroka?locale.x=en_GB)
. 
