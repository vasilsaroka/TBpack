(* Wolfram Language Package *)

BeginPackage["TBpack`Optics`"]

Unprotect[Evaluate[$Context<>"*"]]; (* taken from CustomTicks package *)

(* Usage messages: use here String Representation of Boxes to get formatting similar to built-in functions *)

(* Functions *)
OpticalAbsorption1D::usage = "OpticalAbsorption1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"wavefunctions\",\"TI\"], \*StyleBox[\"velocityoperators\",\"TI\"], \*StyleBox[\"frequencyrange\",\"TI\"], \*StyleBox[\"numberofpoints\",\"TI\"], \*StyleBox[\"gaussianbroadening\",\"TI\"], \*StyleBox[\"polarization\",\"TI\"] \)] returns optical absorption spectrum for the 1D structure and linear polarization of the incident light.";


(* Options *)


Begin["`Private`"] (* Begin Private Context *) 


(* This function works for the zero temperature and the intrinsic position of the Fermi level *)

(* Note: this function uses input from the ElectronicBands1D function *)

OpticalAbsorption1D = Compile[
   {
    {bandstructure, _Real, 2},
    {wavefuncs, _Complex, 3},
    {velocityoperators, _Complex, 4},
    {frequencyrange, _Real, 1},
    {numberofpoints, _Integer},
    {gaussianbroadening, _Real},
    {polarizationindex, _Integer}
    },
   Module[
    {
     energybands,
     num, len,
     wavefunctions,
     velocities,
     numhalf,
     tol = 10^-6,
     tdata,
     energy, vme,
     x1, x2, xrange,
     tdatalen, r
     },
    energybands = Most[bandstructure];
    {num, len} = Dimensions[energybands];
    wavefunctions = wavefuncs;
    velocities = velocityoperators;
    
    numhalf = Floor[num/2];
    
    (* transition probability rates for all k points *)
    tdata = Flatten[
      Table[
       	energy = 
        energybands[[finalbandnumber, kindex]] - 
         energybands[[initialbandnumber, kindex]];
       	vme = (wavefunctions[[kindex, 
            finalbandnumber]])\[Conjugate].(velocities[[kindex, 
            1 ;; num, 1 ;; num, polarizationindex]].wavefunctions[[
            kindex, initialbandnumber]]);
       	{energy, If[energy < tol, 0, 1/energy Abs[vme]^2]},
       {kindex, len},
       {initialbandnumber, numhalf},
       {finalbandnumber, numhalf + 1, num}]
      , 2];
    
    (* reduction and smoothing of the optical absorption spectrum *)
    x1 = frequencyrange[[1]];(* min(left) value *)
    x2 = frequencyrange[[2]];(* max(right) value *)
    xrange = Range[x1, x2, (x2 - x1)/(numberofpoints - 1)];
    tdatalen = Length[tdata];
    Table[
     	r = 0.0 + I 0.0;
     	Do[
      		r = 
       r + tdata[[j, 
          
          2]] Exp[-((tdata[[j, 1]] - xrange[[i]])/
            gaussianbroadening)^2],
      	{j, tdatalen}];
     {xrange[[i]], Re@r},
     {i, numberofpoints}]
    ](* end Module *), 
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True
   ](* end Compile *);


End[] (* End Private Context *)

(Attributes[#] = {Protected, ReadProtected}) & /@ Names[Evaluate[$Context<>"*"]]

EndPackage[]