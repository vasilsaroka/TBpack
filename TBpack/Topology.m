(* Wolfram Language Package *)

BeginPackage["TBpack`Topology`"]

Unprotect[Evaluate[$Context<>"*"]]; (* taken from CustomTicks package *)

(* Usage messages: use here String Representation of Boxes to get formatting similar to built-in functions *)
 

(* Functions *)
BerryPhaseKudin::usage = "BerryPhaseKudin[\!\(\*StyleBox[\"data\",\"TI\"] \)] returns a Berry phase for \!\(\*StyleBox[\"data\",\"TI\"]\) that is an output of ElectronicStructure function run on a closed loop \!\(\*StyleBox[\"k\",\"TI\"]\)-path with an option EigenVectors\[Rule]True and HamiltonianGauge\[Rule]\"Periodic\".
See Sec. IV and Eq. (31) in \!\(\*TemplateBox[{\"K. N. Kudin, R. Car, and R. Resta, J. Chem.  Phys. 126, 23 (2007)\", \"https://doi.org/10.1063/1.2743018\"}, \"HyperlinkURL\"]\)";

SortBipartiteLattice::usage = "SortBipartiteLattice[\!\(\*StyleBox[\"unitcell\",\"TI\"], \*StyleBox[\"translations\",\"TI\"], \*StyleBox[\"nnvectors\",\"TI\"] \)] returns a sorted \!\(\*StyleBox[\"unitcell\",\"TI\"]\), wherein the sublattice with the nearest-neighbor \!\(\*StyleBox[\"nnvectors\",\"TI\"]\) comes first.
SortBipartiteLattice[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\", StyleBox[\"translations\",\"TI\"], \",\", StyleBox[\"nnvectors\",\"TI\"], \",\" , \"LatticeConstant\" , \"\[Rule]\" , StyleBox[\"value\",\"TI\"]}]\)] returns the sorted \!\(\*StyleBox[\"unitcell\",\"TI\"]\) with the characteristic bond length given by \!\(\*StyleBox[\"value\",\"TI\"]\) (defaults to 1.42 \[Angstrom]).
SortBipartiteLattice[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\", StyleBox[\"translations\",\"TI\"], \",\", StyleBox[\"nnvectors\",\"TI\"], \",\" , \"BondLengthDelta\" , \"\[Rule]\" , StyleBox[\"value\",\"TI\"]}]\)] returns the sorted \!\(\*StyleBox[\"unitcell\",\"TI\"]\), where bonds may differ from the characteristic bond length by \!\(\*StyleBox[\"value\",\"TI\"]\) (defaults to 0.05 \[Angstrom]).
SortBipartiteLattice[\!\(\*RowBox[{StyleBox[\"unitcell\",\"TI\"], \",\", StyleBox[\"translations\",\"TI\"], \",\", StyleBox[\"nnvectors\",\"TI\"], \",\" , \"SuperCellSize\" , \"\[Rule]\" , StyleBox[\"integer\",\"TI\"]}]\)] returns the sorted \!\(\*StyleBox[\"unitcell\",\"TI\"]\), where the nearest-neighbor bonds are traced throughout the supercell generated via \!\(\*StyleBox[\"translations\",\"TI\"]\) ranged from \!\(\*StyleBox[\"-integer\",\"TI\"]\) to \!\(\*StyleBox[\"+integer\",\"TI\"]\) (defaults to 1); useful for accounting for intercell bonds.
SortBipartiteLattice[\!\(\*StyleBox[\"unitcell\",\"TI\"]\), \!\(\*StyleBox[\"translations\",\"TI\"]\), \!\(\*StyleBox[\"nnvectors\",\"TI\"]\), \"TrackBondsConfiguration\" \[Rule] True] returns the sorted \!\(\*StyleBox[\"unitcell\",\"TI\"]\) acompanied by the bond configuration used in the sorting."

ZInvariantAIII1D::usage = "ZInvariantAIII1D[\!\(\*StyleBox[\"data\",\"TI\"] \)] returns a winding number for \!\(\*StyleBox[\"data\",\"TI\"]\) that is an output of ElectronicStructure function run on a closed loop \!\(\*StyleBox[\"k\",\"TI\"]\)-path with an option EigenVectors\[Rule]True and HamiltonianGauge\[Rule]\"Periodic\".
ZInvariantAIII1D[\!\(\*StyleBox[\"data\",\"TI\"]\), \"TrackDet\" \[Rule] True] returns the winding number accompanied by a plot of the determinant of Q-matrix as function of \!\(\*StyleBox[\"k\",\"TI\"]\)-points to track its smoothness and make sure quality of the winding number obtained.
See Sec. 2.2 and Eqs. (20,21) in \!\(\*TemplateBox[{\"Ryu, S., Schnyder, A. P., Furusaki, A. and Ludwig, A. W. W. New J. Phys. 12, 065010 (2010)\", \"https://doi.org/10.1088/1367-2630/12/6/065010\"}, \"HyperlinkURL\"]\)"

BerryCurvature2D::usage = "BerryCurvature2D[\!\(\*StyleBox[\"bandnum\",\"TI\"], \*StyleBox[\"data\",\"TI\"]\)] returns a list of Berry curvature data suitable for ListPlot3D for the given \!\(\*StyleBox[\"bandnum\",\"TI\"]\) band and \!\(\*StyleBox[\"data\",\"TI\"]\) that is an output of ElectronicStructure function run on a 2D Brillouin zone with an option EigenVectors\[Rule]True and HamiltonianGauge\[Rule]\"Canonical\".
See Eq. (10) in \!\(\*TemplateBox[{\"M. V. Berry, Proc. R. Soc. A 392, 45-57 (1984)\", \"https://doi.org/10.1098/rspa.1984.0023\"}, \"HyperlinkURL\"]\)"



(* Options *)


Begin["`Private`"] (* Begin Private Context *) 

(* Calculates unit cell invariant Berry-Zak phase [intercellular Zak phase] as prescribed in 

K. N. Kudin, R. Car, and R. Resta, J. Chem.  Phys. 126, 23 (2007)

*)
(* data is an output of ElectronicBands function for a closed-loop k-path, i.e. 
the first and the last points of this loop must be the same;
Periodic Hamiltonian gauge must be used 
*)
(* may be prone to overflows due to the logorithm of the product taken instead of sum of the logorithms *)
BerryPhaseKudin[data_] := Block[
   {
    ph = 1,
    bands, wfs, velops,
    klist, len, klen,
    wf1, wf2, S
    },
   
   {bands, wfs, velops} = data;
   klist = Last@bands;
   klen = Length[klist];
   len = Length[Most@bands]/2;
   
   Do[
    wf1 = wfs[[k]];
    wf2 = wfs[[k + 1]];
    S = Table[Conjugate[wf1[[n]]] . (wf2[[m]]), {n, len}, {m, len}];
    ph = ph Chop@Det[S],
    {k, klen - 1}];
   -Im[Log[ph]]
   ];
SyntaxInformation[BerryPhaseKudin] = {"ArgumentsPattern" -> {_}};


Options[SortBipartiteLattice] = {
   TBpack`LatticeConstant -> 1.42,
   TBpack`BondLengthDelta -> 0.05,
   TBpack`SuperCellSize -> 1,
   "TrackBondsConfiguration" -> False
   };
SortBipartiteLattice[unitcell_, translationvectors_, 
   sublatticeSiteBondVectors_, OptionsPattern[]] := Block[
   {
    a0 = OptionValue[TBpack`LatticeConstant],
    bldelta = OptionValue[TBpack`BondLengthDelta],
    supercellsize = OptionValue[TBpack`SuperCellSize],
    trackbonds = OptionValue["TrackBondsConfiguration"],
    
    codirQ, testQ,
    
    tlen, len, aplist, T1, T2, T3, alen,
    
    nnvectorlist, bondlength, bondvector,
    bondconfiglist, res,
    
    vfun
    },
   (* test for the alignment of two vectors *)
   codirQ = If[Sqrt[(#1 - #2).(#1 - #2)] < bldelta, True, False] &;
   (* A sublattice nnvectors: nearest-neighbor vectors for A sublattice *)
   testQ = (Or @@ 
       Flatten@Outer[codirQ, #[[2]], sublatticeSiteBondVectors, 1]) &;
   (* B sublattice, is defined as the one that is not an A sublattice *)
   
   (* supercell construction *)
   (* supercell is the initial unitcell + adjacent unitcells *)
   tlen = Length[translationvectors];
   len = Length[unitcell];
   
   (* atom position list *)
   aplist = Switch[
     tlen,
     0,
     unitcell,
     1,
     T1 = translationvectors[[1]];
     Flatten[
      Table[(#1 + i T1) & /@ unitcell, {i, -supercellsize, 
        supercellsize}], 1],
     2,
     T1 = translationvectors[[1]];
     T2 = translationvectors[[2]];
     Flatten[
      Table[(#1 + i T1 + j T2) & /@ unitcell, {i, -supercellsize, 
        supercellsize}, {j, -supercellsize, supercellsize}], 2],
     3,
     T1 = translationvectors[[1]];
     T2 = translationvectors[[2]];
     T3 = translationvectors[[3]];
     Flatten[
      Table[(#1 + i T1 + j T2 + k T3) & /@ 
        unitcell, {i, -supercellsize, 
        supercellsize}, {j, -supercellsize, 
        supercellsize}, {k, -supercellsize, supercellsize}], 3]
     ](* end Switch *);
   
   alen = Length[aplist];
   
   (* find bond configuration list *)
   bondconfiglist = Table[
     nnvectorlist = {};
     Do[
      bondvector = aplist[[j]] - unitcell[[i]];
      bondlength = Sqrt[bondvector . bondvector];
      If[
       Abs[bondlength - a0] < bldelta,
       AppendTo[nnvectorlist, bondvector]
       ],
      {j, alen}](*end Do*);
     {unitcell[[i]], nnvectorlist},
     {i, 1, len}](*end Table*);
   
   If[
    trackbonds,
    (* 2D projection to xy-plane *)
    vfun = {
       Red,
       Disk[Most[#[[1]]], 0.3],
       Gray, Arrowheads[0.02], Table[
        Arrow[{Most[#[[1]]], Most[#[[1]] + #[[2, i]]]}],
        {i, Length[#[[2]]]}]
       } &;
    Print[
     Graphics[{
       Opacity[0.3], Blue, Arrowheads[0.03],
       Table[
        Arrow[
         {
          Most[unitcell[[1]]], 
          Most[unitcell[[1]] + sublatticeSiteBondVectors[[i]]]
          }]
        , {i, Length[sublatticeSiteBondVectors]}],
       vfun /@ bondconfiglist
       }]
     ]
    ];
   
   res = Reap[If[testQ[#], Sow[#[[1]], "A"], Sow[#[[1]], "B"]] & /@ bondconfiglist][[2]];
   {Flatten[res, 1], Length /@ res}
   ](* end Block *);
SyntaxInformation[SortBipartiteLattice] = {"ArgumentsPattern" -> {_,_,_, OptionsPattern[]}};


(* Calculates winding number (Z invariant) as presented in Eqs.(20, 21) in 
Ryu, S., Schnyder, A. P., Furusaki, A. and Ludwig, A. W. W. New J. Phys. 12, 065010 (2010)
and in Eq. (4) in J. Jiang, and S. G. Louie, Nano Lett. 21, 197 (2021)

This function implies chiral structure of the Hamiltonian, i.e. we first enumerate atoms 
in sublattice 1 and then in the sublattice 2. 
k-path may be not exactly closed, i.e. the first and the last k-points can be different 
but very close to each other unlike in Kudin's approach
*)
Options[ZInvariantAIII1D] = {
   "TrackDet" -> False
   };
ZInvariantAIII1D[data_, OptionsPattern[]] := Block[
  {
   trackdet = OptionValue["TrackDet"],
   
   bands, wfs, velops,
   klen, len, s, wf,
   P, Q, q,
   
   nag, ag,
   
   dlist = {}
   },
  
  {bands, wfs, velops} = data;
  klen = Length[wfs];
  len = Length[Most@bands];
  s = 0;
  If[
   trackdet,
   Do[
    wf = Chop[wfs[[ind, 1 ;; len/2]]];
    P = Plus @@ (KroneckerProduct[#, #\[HermitianConjugate]] & /@ 
        wf);
    Q = Chop[IdentityMatrix[len] - 2 P];
    q = Q[[1 ;; len/2, len/2 + 1 ;; len]];
    nag = Arg[Det[q]];
    If[
     ind != 1,
     If[
      	nag - ag > \[Pi],
      	nag = nag - 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]);
      	s = s + nag - 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]) - ag,
      	If[
       		nag - ag < -\[Pi],
       		nag = nag + 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]);
       		s = s + nag + 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]) - ag,
       		s = s + nag - ag
       	]
      ]
     ];
    ag = nag;
    AppendTo[dlist, ag],
    {ind, klen}];
   Print[ListPlot[dlist]];
   s/(2 \[Pi]),
   Do[
    wf = Chop[wfs[[ind, 1 ;; len/2]]];
    P = Plus @@ (KroneckerProduct[#, #\[HermitianConjugate]] & /@ 
        wf);
    Q = Chop[IdentityMatrix[len] - 2 P];
    q = Q[[1 ;; len/2, len/2 + 1 ;; len]];
    nag = Arg[Det[q]];
    If[
     ind != 1,
     If[
      	nag - ag > \[Pi],
      	nag = nag - 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]);
      	s = s + nag - 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]) - ag,
      	If[
       		nag - ag < -\[Pi],
       		nag = nag + 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]);
       		s = s + nag + 2 \[Pi] (Round[Abs[nag - ag]/(2 \[Pi])]) - ag,
       		s = s + nag - ag
       	]
      ]
     ];
    ag = nag,
    {ind, klen}];
   s/(2 \[Pi])
   ](* end If track determinant *)
  ](* end Block *)
SyntaxInformation[ZInvariantAIII1D] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

BerryCurvature2D[i_, data_] := Block[
   {
    len, klen,
    evals, evecs, velops,
    vx, vy, vz,
    kx, ky,
    Bcurv,
    Eik, Emk,
    vki, vkm,
    vme
    },
   {evals, evecs, velops} = data;
   len = Length[evals] - 1;
   klen = Length[Last@evals];
   
   Table[
    vx = velops[[k, ;; , ;; , 1]];
    vy = velops[[k, ;; , ;; , 2]];
    vz = velops[[k, ;; , ;; , 3]];
    kx = evals[[-1, k, 1]];
    ky = evals[[-1, k, 2]];
    			
    Bcurv = {0.0, 0.0, 0.0};
    Do[
     	If[
      		m != i,
      		Eik = evals[[i, k]];
      		Emk = evals[[m, k]];
      		If[
       			Emk != Eik,
       			vki = evecs[[k, i]];
       			vkm = evecs[[k, m]];
       			
       vme = {Conjugate[vki] . (vx . vkm), 
         Conjugate[vki] . (vy . vkm), Conjugate[vki] . (vz . vkm)};
       			
       Bcurv = Bcurv + 
         Cross[vme, 
          Conjugate[vme]]/(Emk - Eik)^2
       		](* end If not degenerated energies *)
      ](* end If not the same bands *),
     {m, len}](* end Do *);
    {kx, ky, Chop[-Im[Bcurv . {0, 0, 1}]]},
    {k, klen}]
   ](* end Block *);
SyntaxInformation[BerryCurvature2D] = {"ArgumentsPattern" -> {_, _}};

End[] (* End Private Context *)

(Attributes[#] = {Protected, ReadProtected}) & /@ Names[Evaluate[$Context<>"*"]] 

EndPackage[]