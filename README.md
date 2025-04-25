[![GitHub (pre-)release](https://img.shields.io/github/release/vasilsaroka/TBpack/all.svg)](https://github.com/vasilsaroka/TBpack/releases)
[![GitHub Latest Release](https://img.shields.io/github/downloads/vasilsaroka/TBpack/0.5.2.svg)](https://github.com/vasilsaroka/TBpack/releases/tag/v0.5.2)
[![Github All Releases](https://img.shields.io/github/downloads/vasilsaroka/TBpack/total.svg)](https://github.com/vasilsaroka/TBpack/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Support TBpack](https://img.shields.io/static/v1?label=support&message=5$&color=green&style=flat&logo=paypal)](https://paypal.me/vasilsaroka?locale.x=en_GB)

# TBpack
Tight-binding calculations in Mathematica

## Installation guide
 - **The recommended option for TBpack paclet management in Mathematica 11.3+** is to use the [resource function](https://resources.wolframcloud.com/FunctionRepository/resources/InstallTBpack): [<img align="right" src="https://github.com/vasilsaroka/TBpack/blob/master/TBpack/TBpackPro_Logo.png" alt="TBpackPro_Logo"/>](https://paypal.me/vasilsaroka?locale.x=en_GB)
 
   To install the latest version of the TBpack from this GitHub repo, evaluate
   
       ResourceFunction["InstallTBpack"][]

   To find all installed versions on your PC and to uninstall chosen (all) version (s), evaluate
   
       ResourceFunction["InstallTBpack"][Method->"Uninstall"]
   **Note: This and other options below may not work if the [latest release](https://github.com/vasilsaroka/TBpack/releases) is a special patch intended for example for a WLJS notebook. If this is the case, use the last option of this section with ``PacletInstall["url"]`` specifying explicitly a version of TBpack released for Mathematica.**

        
 - In Mathematica 12.1+ ``PacletInstall["url"]`` can be used for the installation straight away:
        
        PacletInstall["https://github.com/vasilsaroka/TBpack/releases/download/v<version>/TBpack-<version>.paclet"]  
   where `<version>` stands for any existing version of the application. For example,
   
       PacletInstall["https://github.com/vasilsaroka/TBpack/releases/download/v0.2.0/TBpack-0.2.0.paclet"]

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

## TBpack for fun
[<img align="right" src="https://github.com/vasilsaroka/TBpack/blob/master/Games/HOMOLUMOQuantumGame/HOMOLUMOQuantumGame.png" alt="HOMOLUMOQuantumGame" width="250"/>](https://t.me/HOMO_LUMO_Quantum_Game_bot)

**TBpack** application is fast enough to be used in games and Telegram bots. 

Put your knowledge and luck under ultimate test with Hostess Monica in HOMO-LUMO Quantum Game: [@HOMO_LUMO_Quantum_Game_bot](https://t.me/HOMO_LUMO_Quantum_Game_bot) **\***. Use `/rules` to see the rules of the game or have a look at our sample [game video](https://www.youtube.com/watch?v=DUy22OKLgvs) on YouTube.

**\*** *For now the game is not accessible due to being in the process of updating. Sorry for this inconvenience. You will be able to try it soon with new features.*


## TBpack for WLJS Notebooks - Electron application
   The [Electron application](https://github.com/JerryI/wolfram-js-frontend/releases) is an open-source javascript-based cross-platform graphical user interface for running Wolfram Language code via a [Wolfram Engine](https://www.wolfram.com/engine/). It is an actively developed [project](https://github.com/JerryI/wolfram-js-frontend) that is getting popularity due to its simplicity and ability to combine Wolfram Language with Javescript and Markdown. While WLJS notebooks are not superior to the original Mathematica notebooks, they do have some features such as gpt-based copilot, powerful command palette system, and enhanced animation capbilitites that may transform them into a useful companion if not a real contender of original notebooks. 
   
   We have release a special patch v0.5.1 of TBpack that can be installed into the Electron application.  
   ### Installation guide
   - Install [Wolfram Engine](https://www.wolfram.com/engine/) (14.0 or higher)
   - Install [Electron application](https://jerryi.github.io/wljs-docs/) (2.5.6 or higher)
   - Proceed to an automated installation of TBpack:

     Copy-paste the below function into a WLSJ notebook cell. Evaluate the cell, for example by pressing `Shift+Enter`, to make the definition of this function known to the Wolfram Engine. In the next cell, type `InstallTBpackIn2Electron[]` and evaluate it.

           Options[InstallTBpackIn2Electron]={Method->"Install" (* or "Uninstall" *) };
           InstallTBpackIn2Electron[OptionsPattern[]]:=Block[
           {
              method=OptionValue[Method],
              pacletsrepository, tbpackdirectory, message1, message2, jsonreleases, info, pos, url, message, tempfile, 
              version = "v0.5.1" (* special TBpack patch for Electron *),
              giturl = "https://api.github.com/repos/vasilsaroka/TBpack/releases",
              deletedirectory
           },
           (* auxiliary routine *)
           deletedirectory[dir_,successmessage_,failuremessage_]:=Block[{res},
           If[
                 DirectoryQ[dir],
                 res = DeleteDirectory[dir,DeleteContents->True];
                 If[res===Null, Print[successmessage], Print[failuremessage]]
            ];
            ](* end Block *);
            (* set up directory for the installation *)
            pacletsrepository = FileNameJoin[$InitialDirectory,"wl_packages"];
            tbpackdirectory = FileNameJoin[{pacletsrepository,"TBpack-"<>StringTake[version,{2,-1}]}];
            (* some general information messages *)
            (* success message *)
            message1 = "TBpack has been succefully uninstalled from Electron application.";
            (* failure message *)
            message2 = "Failed to uninstall the existing TBpack installation. Probaply due to its current usage by another application.";
            (* select the method *)
            Switch[
                    method,
                    "Install",
                    (* check that an installation folder does not contain TBpack version *)
                    (* clear the folder from the previously installed TBpack version of the same type as is to be installed *)
                    deletedirectory[tbpackdirectory,message1,message2];       
                    (* proceed to the installation/re-installation *)
                    jsonreleases = Import[giturl,"JSON"];
                    If[
                       jsonreleases === $Failed, 
                       message = "Cannot access GitHub releases: "<>giturl<>". Please, check your internet connection.";
                       Print[message]; Return[$Failed],
                       (* successful access to GitHub releases *)
                       info = "Looking for TBpack " <> version;
                       Print[info <> "..."];
                       pos=Position[Lookup[jsonreleases, "tag_name"],version][[1,1]];
                       url=First@Lookup[Lookup[jsonreleases, "assets"][[pos]],"browser_download_url"];
                       If[
                          $VersionNumber >= 14.0,
                          tempfile = FileNameJoin[{$TemporaryDirectory,FileNameTake[url]}];
                          Check[tempfile = URLDownload[url,tempfile],Return[$Failed]];
                          (* check if the downloaded paclet file still exist on disk and install it *)
                          If[
                             FileExistsQ[tempfile],
                             pacletsrepository=FileNameJoin[$InitialDirectory,"wl_packages"];          
                             Check[ExtractPacletArchive[tempfile,pacletsrepository],Return[$Failed]];
                             DeleteFile[tempfile];
                             message = "TBpack is succefully installed into Electron application. It is located at "<>pacletsrepository<>". Now close Electron application and open it again. After this use Needs[''TBpack`''], where '' is meant to be \" symbol (double quotation mark, entered without backslash), to activate TBpack for the session.";
                             Print[message],
                             message="Installation failed...";
                             Print[message]; Return[$Failed]
                          ],
                          message="It is recommended to use the latest Wolfram Engine that is version 14.0 or higher. Update Wolfram Engine and repeat installation of TBpack.";
                          Print[message]; Return[$Failed]
                       ](* end If version *)
                    ](* end If GitHub releases*),
                    "Uninstall",
                    (* clear the folder from the previously installed TBpack *)
                    deletedirectory[tbpackdirectory,message1,message2]; 
            ](* end Switch *)
          ](* end Block *)

      **Note:** After successful installation, you will have to restart Electron application to make TBpack visible for ``Needs["TBpack`"]``.

     To uninstall the installed TBpack version, evaluate

         InstallTBpackIn2Electron[Mehtod->"Uninstall"]

   - Manual installation of TBpack is possible (try only if automated option above has failed):
      - Download `.paclet` file of the [TBpack patch released](https://github.com/vasilsaroka/TBpack/releases/tag/v0.5.1) for the Electron application.
      - Extract the paclet file into the installation folder by evaluating in the Electron application cell the following line:

            ExtractPacletArchive["path2paclet",FileNameJoin[$InitialDirectory,"wl_packages"]]
        
      - Restart Electron application to make TBpack visible for ``Needs["TBpack`"]``.

   ### WLJS Notebook Demo
   - Evaluate ``Needs["TBpack`"]`` to load TBpack into the Wolfram Engine session run in Electron application.
   - Test TBpack using `AtomicStructure[Nanotube[5, 5]]`:
     ![WLJS Notebook Demo](https://github.com/user-attachments/assets/3a7b937f-36dd-42e8-85db-b3ec23bc5a72)

## Supporting the project

   We believe everyone deserves access to knowledge that is grounded in science and integrity. That is why we keep our code open for all users, regardless of where they live or what they can afford to pay. This means more people can be better educated and inspired to make an impact on the global wellbeing. We have no billionaire owner, meaning only your donations power our work and ensure it can remain open for all. Every contribution, however big or small, makes a real difference for TBpack future. If you find it useful, consider [supporting the project](https://paypal.me/vasilsaroka?locale.x=en_GB). 
