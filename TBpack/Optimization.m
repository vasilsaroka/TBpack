(* Wolfram Language Package *)

BeginPackage["TBpack`Optimization`"]

Unprotect[Evaluate[$Context<>"*"]]; (* taken from CustomTicks package *)

(* Exported symbols added here with SymbolName::usage *)  

(* for Functions *)
OptimizationByLAMMPS::usage = "OptimizationByLAMMPS[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\" , StyleBox[\"tr\",\"TI\"]}]\)] optimizes geometry of a 1D system presented by \!\(\*StyleBox[\"unitcell\",\"TI\"]\) and the translation vector \!\(\*StyleBox[\"tr\",\"TI\"]\) using LAMMPS package;
OptimizationByLAMMPS[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\" , StyleBox[\"tr\",\"TI\"], \",\" , StyleBox[\"options\",\"TI\"]}]\)] optimizes geometry using specified options.";
OptimizationByGULP::usage = "OptimizationByGULP[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\" , StyleBox[\"tr\",\"TI\"]}]\)] optimizes geometry of a 1D system presented by \!\(\*StyleBox[\"unitcell\",\"TI\"]\) and the translation vector \!\(\*StyleBox[\"tr\",\"TI\"]\) using GULP package;
OptimizationByGULP[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\" , StyleBox[\"tr\",\"TI\"], \",\", StyleBox[\"options\",\"TI\"]}]\)] optimizes geometry using specified options.";


(* for Options *)
AtomMass::usage = "Option setting the mass of an atom of the nanotructure in functions such as \*StyleBox[\"OptimizationByLAMMPS\",\"TI\"].";
EdgeAtomMass::usage = "Option setting the mass of an edge atom of the nanotructure in functions such as \*StyleBox[\"OptimizationByLAMMPS\",\"TI\"].";
EdgeBondLength::usage = "Option setting for a nanotructure the length of an edge bond in \[Angstrom].";
AtomLabel::usage = "Option setting a chemical element label for the base atom of a nanostructure.";
EdgeAtomLabel::usage = "Option setting a chemical element label for the edge atom of a nanostructure.";
BulkCoordinationNumber::usage = "Option setting the coordination number for atoms in the interior of a nanostructure.";
EdgePassivation::usage = "Option specifying whether dungling bonds at the edges of a nanostructure are passivated with edge atoms belonging to a different element compared to the atoms in the bulk. For example, hydrogen edge atoms for a carbon nanostructure.";

Begin["`Private`"] (* Begin Private Context *)


cleardir[path2dir_, basefiles_] := Block[
	{
		flist
	},
	Off[General::privv];
	flist = Complement[FileNames[All, path2dir], basefiles];
	While[flist =!= {},
		DeleteFile[flist];
		Pause[0.1];
		flist = Complement[FileNames[All, path2dir], basefiles];
		
	];
	On[General::privv];
](* end Block *);

EdgeAtomPositions[unitcell_, translationvector_, XXbondlength_, delta_, XYbondlength_, bulkcoordinationnumber_] := Block[
   {
    supercell, lenj, leni,
    atomposi, atomposj,
    nnatoms, counter,
    v, vec,
    edgeatomspos = {}
    },
   supercell = Flatten[Table[(# + i translationvector) & /@ unitcell, {i, -1, 1}], 1];
   lenj = Length@supercell;
   leni = Length@unitcell;
   Do[
    	nnatoms = {};
    	counter = 0;
    	atomposi = unitcell[[i]];
    	Do[
     		atomposj = supercell[[j]];
     		v = atomposi - atomposj;
     		If[
      			Abs[Sqrt[v.v] - XXbondlength] < delta,
      			AppendTo[nnatoms, atomposj];
      			counter++
      		];
     		If[
      			counter == bulkcoordinationnumber,
      			Break[]
      		],
     	{j, lenj}];
    	Switch[
     			Length@nnatoms,
     			2,
     			vec = -(nnatoms[[1]] + nnatoms[[2]] - 2 atomposi);
     			AppendTo[edgeatomspos, atomposi + vec/Sqrt[vec.vec] XYbondlength],
     			1,
     			vec = -(nnatoms[[1]] - atomposi);
     			v = (vec/Sqrt[vec.vec] XYbondlength);
     			AppendTo[edgeatomspos, atomposi + v.RotationMatrix[\[Pi]/3, {0, 0, 1}]];
     			AppendTo[edgeatomspos, atomposi + v.RotationMatrix[-\[Pi]/3, {0, 0, 1}]]
     	](* end Switch *),
    {i, leni}];
   edgeatomspos
](* end Block *);


(*Rotation to make x-direction the translation one*)
RotateTranslation2Direction[unitcell_, translationvector_, dir_] := 
  Block[
   {
    tv, rt,
    uc, rm
    },
   tv = Chop[translationvector];
   rt = Cross[tv, dir];
   If[
    	Chop@(rt.rt) == 0,
    	uc = unitcell,
    	rm = RotationMatrix[{tv, dir}];
    	uc = (rm.#) & /@ unitcell;
    	tv = rm.tv
    ](*end If*);
    {uc, tv}
](* end Block *);


(* 
   		This funcion creates read_data file to specify initial atom coordinates in the lammps input script
   		see https://docs.lammps.org/read_data.html for details
*)
(* It automatically passivates edge bonds of the structure with edge atoms *)
UnitCell2LammpsReadDataFile[unitcell_, translationvector_, XXbondlength_, delta_, XYbondlength_, atommass_, edgeatommass_, bulkcoordinationnumber_,
   			edgepassivation_, path_, fname_] := Module[
     {
      dir = {1, 0, 0},(* direction of translation *)
      uc,tv,
      edgeatomspos,
      vert, str, coo, min, max, text,
      firstline, headerkeywords, sectionkeywords, header, body,
      
      ps,nd
      },

     {uc,tv} = RotateTranslation2Direction[unitcell, translationvector, dir];
           
     edgeatomspos = If[edgepassivation, EdgeAtomPositions[uc, tv, XXbondlength, delta, XYbondlength, bulkcoordinationnumber], {}];
     
     (* atomic coordinates file creation *)
     vert = Flatten[{ConstantArray[1, Length@uc], ConstantArray[2, Length@edgeatomspos]}, 1];
     coo = Chop@Flatten[{uc, edgeatomspos}, 1];
     min = Min /@ Transpose@coo;
     max = Max /@ Transpose@coo;
     firstline = "Nanostructure";
     headerkeywords = {"atoms", "atom types", "xlo xhi", "ylo yhi", "zlo zhi"};
     sectionkeywords = {"Atoms \n", "Masses \n"};
     
     ps = Row[#, "\t"] &;(* parts separation function *)
     nd = PaddedForm[Chop@#, {13, 10}, ExponentFunction -> (If[-11 < # < 11, Null, #] &)] &; (* number display function *)
  

     header = {
     	ps@{Length@coo, headerkeywords[[1]]},
       	ps@{2, headerkeywords[[2]]},
       	ps@{nd@min[[1]], nd@(min[[1]] + Sqrt[tv.tv]), headerkeywords[[3]]},
       	ps@{nd@(min[[2]] - 2 (max[[2]] - min[[2]]) - 10), nd@(max[[2]] + 2 (max[[2]] - min[[2]]) + 10), headerkeywords[[4]]},
       	ps@{nd@(min[[3]] - 2 (max[[3]] - min[[3]]) - 10), nd@(max[[3]] + 2 (max[[2]] - min[[2]]) + 10), headerkeywords[[5]]}
       };
     
     body = {
     	sectionkeywords[[2]],
       	ps /@ {{1, atommass}, {2, edgeatommass}},
       	sectionkeywords[[1]],
       	ps /@ Table[{i, vert[[i]], nd@coo[[i, 1]], nd@coo[[i, 2]], nd@coo[[i, 3]]}, {i, 1, Length@coo}]
       };
     
     text = Flatten@{firstline, header, body};
     str = OpenWrite[FileNameJoin[{path, fname}], FormatType -> OutputForm, CharacterEncoding -> "ASCII"];
     Do[Write[str, text[[i]]], {i, 1, Length@text}];
     Close[str];
     tv
](* end Module *);
   

(* Error messages *)
OptimizationByLAMMPS::optfail = "Optimization failed";
OptimizationByLAMMPS::OS = "The optimization under this type of the operational system is not supported.";

Options[OptimizationByLAMMPS] = {
	TBpack`Path2File -> Automatic,
	TBpack`FileName -> Automatic,
	AtomMass -> 12.0107,
	EdgeAtomMass -> 1.0079,
	TBpack`LatticeConstant -> 1.42,
	EdgeBondLength -> 1.1,
	TBpack`BondLengthDelta -> 0.1,
	BulkCoordinationNumber -> 3,
	EdgePassivation -> True,
	TBpack`InputScriptTemplate -> "
units           metal
newton          on
dimension       3
boundary        p p p
atom_style      atomic
read_data       `readdatafname` 
pair_style      airebo 3.0
pair_coeff      * * CH.airebo C H
neigh_modify    every 1 delay 10 check yes
min_style       cg
min_modify      line quadratic
dump            1 all xyz 10000 `dumpfilename`
dump_modify     1 append yes
minimize        1.0e-4 1.0e-6 15000 30000
variable        e equal etotal
print           \"@ $e\""
};

OptimizationByLAMMPS[unitcell_List, translationvector_List, OptionsPattern[]] := 
 Catch[Module[
   {
   	path2lammps = OptionValue[TBpack`Path2File],
   	lammpsfile = OptionValue[TBpack`FileName],
   	atommass = OptionValue[AtomMass],
   	edgeatommass = OptionValue[EdgeAtomMass],
   	latticeconstant = OptionValue[TBpack`LatticeConstant],
   	edgebondlength = OptionValue[EdgeBondLength],
   	bondlengthdelta = OptionValue[TBpack`BondLengthDelta],
   	bulkcoordinationnumber = OptionValue[BulkCoordinationNumber],
   	edgepassivation = OptionValue[EdgePassivation],
   	inputscripttemplate = OptionValue[TBpack`InputScriptTemplate],
   	
    basefiles, rf, data, indata, outdata, 
    inunitcell, opunitcell,
    inedgeatoms, opedgeatoms,
    T, irm, dir = {1, 0, 0},
    lammpsinputfile = "lammps_input_relax_temp.txt",
    readdatafname = "mycoo.txt",
    dumpfilename = "out.xyz",
    stream
    },
    
    If[
   		path2lammps === Automatic,
   		path2lammps = FileNameJoin[{TBpack`$TBpackDirectory,"Optimization programs","LAMMPS"}]
   	];
   	
   	
   	If[
   		lammpsfile === Automatic,
   		lammpsfile = "lmp_win_no-mpi.exe"
   	];
   
    basefiles = FileNames[All, path2lammps];
   
   
   (* create the .txt file of atomic coordinates to read data from by lammps,
   the axis of translation symmetry is oriented along Ox-axis *)
   T = UnitCell2LammpsReadDataFile[
   			unitcell, 
   			translationvector, 
   			latticeconstant, 
   			bondlengthdelta, 
   			edgebondlength,
   			atommass,
   			edgeatommass,
   			bulkcoordinationnumber,
   			edgepassivation,
   			path2lammps,
   			readdatafname
   ];
   
   stream = OpenWrite[FileNameJoin[{path2lammps, lammpsinputfile}],FormatType -> OutputForm, CharacterEncoding -> "ASCII"];
   Write[stream,StringTemplate[inputscripttemplate][<|"readdatafname"-> readdatafname,"dumpfilename"-> dumpfilename|>]];
   Close[stream];
   
   (* Run optimization by lammps *)
   (* Check OS *)
   Switch[
   			$OperatingSystem,
   			"Windows",
   			rf = RunProcess[$SystemShell, "ExitCode", "cd /d "<>path2lammps<>" & "<>lammpsfile<>" -in "<>lammpsinputfile<>"
   			exit"],
     		_,
     		cleardir[path2lammps,basefiles];
     		Message[OptimizationByLAMMPS::OS];
     		Throw[$Failed]
   ](* end Switch *);
   
   If[ 
   		rf == 0,
    	(* Read optimized data from lammps dumpfile *)
    	data = ReadList[FileNameJoin[{path2lammps, dumpfilename}], Word, RecordLists -> True];
    	data = Map[Read[StringToStream[#], Real] &, Cases[data, {_, _, _, _}], {2}];
    
		indata = data[[1 ;; Length@data/2]];
    	outdata = data[[Length@data/2 + 1 ;; Length@data]];
    	opunitcell = Rest /@ Cases[outdata, {1., _, _, _}];
    	opedgeatoms = Rest /@ Cases[outdata,{2.,_,_,_}];
    	inunitcell = Rest /@ Cases[indata, {1., _, _, _}];
    	inedgeatoms = Rest /@ Cases[indata,{2.,_,_,_}];
    
	    (* Rotate structure back to its original orientation *)
		If[
			Chop[Cross[dir, translationvector].Cross[dir, translationvector]] != 0,
			irm = RotationMatrix[{dir, translationvector}];(* irm -- inverse rotation matrix *)
	    	opunitcell = irm.# & /@ opunitcell;
			T = irm.T;
			inedgeatoms = irm.# & /@ inedgeatoms;
			opedgeatoms = irm.# & /@ opedgeatoms
		];
		
		cleardir[path2lammps,basefiles];
    
	    {{unitcell, translationvector}, {opunitcell, T}, {inedgeatoms, opedgeatoms}}
	    ,
	    cleardir[path2lammps,basefiles];
		Message[OptimizationByLAMMPS::optfail];
		Throw[$Failed]
    ](* end If *)
 ](* end Module *)
](* end Catch *);
SyntaxInformation[OptimizationByLAMMPS] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

UnitCell2GULPinput[unitcell_, translationvector_, XXbondlength_, delta_, XYbondlength_, Xlabel_, Ylabel_, bulkcoordinationnumber_, edgepassivation_, inputscripttemplate_, path_, fname_] := Module[
    {
     dir = {1, 0, 0}(* do not change it! *),
     eta = 1.2,
     ldcutoff = 10 XXbondlength,
     uc, tv,
     blp, urp, 
     sx, sy, sz, ds, s,
     edgeatomspos,
     vert, str, wuc,
     optionkeywords,
     nanostructure,
     
     ps,nd
     },
    
    (* Rotation to make x-direction the translation one *)     
    {uc,tv} = RotateTranslation2Direction[unitcell, translationvector, dir];
           
    edgeatomspos = If[edgepassivation, EdgeAtomPositions[uc, tv, XXbondlength, delta, XYbondlength, bulkcoordinationnumber],{}];
    
    vert = Flatten[{ConstantArray[Xlabel, Length@uc], ConstantArray[Ylabel, Length@edgeatomspos]}, 1];
    
    (* Shift the whole unitcell to positive quater of the reference system *)
    wuc = Flatten[{uc, edgeatomspos}, 1]; (* whole unitcell including Y atoms on the edges *)
    blp = Min /@ Transpose@wuc; (* bottom left point of the rectangle enclosing the whole unitcell *)
    urp = Max /@ Transpose@wuc; (* up right point of the rectangle enclosing the whole unitcell *)
    {sx, sy, sz} = urp - blp; (* shift in positive x-, y-, and z-directions *)
    ds = (eta ldcutoff)/2;(* delta shift in positive y- and z-directions *)
    s = {(Sqrt[tv.tv] - sx)/2, ds, ds} - blp;
    wuc = (# + s) & /@ wuc;
    
    
    (* input file creation *)
    optionkeywords = {"cell", "cartesian"};
    ps = ToString@Row[#, " "] &;(* parts separation function *)
    nd = PaddedForm[1.0 Chop@#, {10, 7}, ExponentFunction -> (If[-11 < # < 11, Null, #] &)] &; (* number display function *)
    
    nanostructure = {
      optionkeywords[[1]],
      ps@({nd@Sqrt[tv.tv], nd@(2 ds + sy), nd@(2 ds + sz), 90, 90, 90, 1, 0, 0, 0, 0, 0}),
      optionkeywords[[2]],
      ps /@ Table[{vert[[i]], "core", nd@wuc[[i, 1]], nd@wuc[[i, 2]], nd@wuc[[i, 3]], 0, 1, 0, 1, 1, 1}, {i, 1, Length@wuc}]
      };
      
    nanostructure = StringRiffle[Flatten[nanostructure],"\n"];
    
    str =  OpenWrite[FileNameJoin[{path, fname}], FormatType -> OutputForm, CharacterEncoding -> "ASCII"];
    Write[str,StringTemplate[inputscripttemplate][<|"nanostructure"-> nanostructure,"AtomLabel"-> Xlabel,"EdgeAtomLabel"->Ylabel|>]];
    Close[str];
    s
](* end Module *);

(* Error messages *)
OptimizationByGULP::optfail = "Optimization failed";
OptimizationByGULP::OS = "The optimization under this type of the operational system is not supported.";

Options[OptimizationByGULP] = {
	TBpack`Path2File -> Automatic,
	TBpack`FileName -> Automatic,
	AtomLabel -> "C",
	EdgeAtomLabel -> "H",
	TBpack`LatticeConstant -> 1.42,
	EdgeBondLength -> 1.1,
	TBpack`BondLengthDelta -> 0.1,
	BulkCoordinationNumber -> 3,
	EdgePassivation -> True,
	TBpack`InputScriptTemplate -> "opti molq cart
`nanostructure`
brenner
lennard inter epsilon
`AtomLabel` `AtomLabel` 0.0024 3.370 1.9 15.0
maxcyc 30000
dump every 100 gulp.res"
};

OptimizationByGULP[unitcell_List, translationvector_List, OptionsPattern[]] := 
 Module[
  {
   path2gulp = OptionValue[TBpack`Path2File],
   gulpfile = OptionValue[TBpack`FileName],
   atomlabel = OptionValue[AtomLabel],
   edgeatomlabel = OptionValue[EdgeAtomLabel],
   latticeconstant = OptionValue[TBpack`LatticeConstant],
   edgebondlength = OptionValue[EdgeBondLength],
   bondlengthdelta = OptionValue[TBpack`BondLengthDelta],
   bulkcoordinationnumber = OptionValue[BulkCoordinationNumber],
   edgepassivation = OptionValue[EdgePassivation],
   inputscripttemplate = OptionValue[TBpack`InputScriptTemplate],
   
   basefiles, rf, data, T1, T2,
   begin, endin, begfi, endfi,
   indata, inunitcell,
   outdata, opunitcell,
   inedgeatoms, opedgeatoms,
   gulpinputfile = "input.txt",
   gulpoutputfile = "out.txt",
   
   s, rt, irm
   },
   
   If[
   		path2gulp === Automatic,
   		path2gulp = FileNameJoin[{TBpack`$TBpackDirectory,"Optimization programs","GULP"}]
   	];
   	
   	
   	If[
   		gulpfile === Automatic,
   		gulpfile = "gulp.exe"
   	];
   
    basefiles = FileNames[All, path2gulp];
    
    s = UnitCell2GULPinput[
  		unitcell, 
  		translationvector, 
  		latticeconstant,
  		bondlengthdelta,
  		edgebondlength, 
  		atomlabel,
  		edgeatomlabel,
  		bulkcoordinationnumber,
  		edgepassivation,
  		inputscripttemplate,
  		path2gulp,
  		gulpinputfile];
  
  (* Run optimization by gulp *)
  (* Check OS *)
   Switch[
   			$OperatingSystem,
   			"Windows",
   			rf = RunProcess[$SystemShell, "ExitCode", "cd /d " <> path2gulp <> " & "<>gulpfile<>"<"<>gulpinputfile<>">"<>gulpoutputfile<>"
   			exit"],
     		_,
     		cleardir[path2gulp,basefiles];
     		Message[OptimizationByGULP::OS];
     		Throw[$Failed]
   ](* end Switch *);
  
  If[
  		rf == 0,
   		(* Read data from gulp output .txt-file *)
   		data = ReadList[FileNameJoin[{path2gulp, gulpoutputfile}], Word, RecordLists -> True];
   		T1 = {ToExpression@data[[Position[data, {"Cell", "parameters", "(Angstroms/Degrees):"}][[1,1]] + 1, 3]], 0, 0};
   		T2 = {ToExpression@data[[Position[data, {"Final", "cell", "parameters", "and", "derivatives",":"}][[1, 1]] + 2, 2]], 0, 0};
   
   		begin = Position[data, {"Initial", "Cartesian", "coordinates", ":"}][[1, 1]] + 5;
   		endin = Position[data, {"Molecule", "list", "generated", "from", "bondlengths", ":"}][[1, 1]] - 2;
   		indata = data[[begin ;; endin]];
   		inunitcell = Map[Read[StringToStream[#], Real] &, #[[4;;6]] & /@ Cases[indata, {_, atomlabel, "c", _, _, _, _, _}], {2}];
   		inedgeatoms = Map[Read[StringToStream[#], Real] &, #[[4;;6]] & /@ Cases[indata, {_, edgeatomlabel, "c", _, _, _, _, _}], {2}];

   		
   		begfi = Position[data, {"Final", "Cartesian", "coordinates", ":"}][[1, 1]] + 5;
   		endfi = Position[data, {"Final", "Cartesian", "lattice", "vectors", "(Angstroms)", ":"}][[1, 1]] - 2;
   		outdata = data[[begfi ;; endfi]];
   		opunitcell = Map[Read[StringToStream[#], Real] &, #[[4;;6]] & /@ Cases[outdata, {_, atomlabel, "c", _, _, _, _, _}], {2}];
   		opedgeatoms = Map[Read[StringToStream[#], Real] &, #[[4;;6]] & /@ Cases[outdata, {_, edgeatomlabel, "c", _, _, _, _, _}], {2}];
   		
   		(* Shift back to the original position and rotate structure back to its original orientation *)
   		rt = Cross[T2, translationvector];
		If[
			Chop[rt.rt] != 0,
			irm = RotationMatrix[{T2, translationvector}];(* irm -- inverse rotation matrix *)
	    	opunitcell = irm.(#-s) & /@ opunitcell;
			T2 = irm.T2;
			inedgeatoms = irm.(#-s) & /@ inedgeatoms;
			opedgeatoms = irm.(#-s) & /@ opedgeatoms
		];
   		
   		cleardir[path2gulp,basefiles];
   		{{unitcell, translationvector}, {opunitcell, T2}, {inedgeatoms, opedgeatoms}},
   		cleardir[path2gulp,basefiles];
   		Message[OptimizationByGULP::optfail];
   		Throw[$Failed]
   ](* end If *)
](* end Module *)
SyntaxInformation[OptimizationByGULP] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};


End[] (* End Private Context *)

(Attributes[#] = {Protected, ReadProtected}) & /@ Names[Evaluate[$Context<>"*"]]

EndPackage[]