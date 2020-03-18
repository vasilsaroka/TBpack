(* ::Package:: *)

(* Mathematica Package  *)

(* :Title: TBpack *)
(* :Author: Vasil A. Saroka <40.ovasil@gmail.com> *)
(* :Context: TBpack` *)
(* :Version: 0.0.1 *)
(* :Date: 2020-03-18 *)

(* :Mathematica Version: 10.0+ *)
(* :Copyright: (c) 2020 Vasil A. Saroka *)

(* Note that all cell must be initialization cells! Initialization Cell makes 
a cell auto-evaluate whenever the notebook is opened and the kernel is launched. *)
BeginPackage["TBpack`"];


Unprotect[Evaluate[$Context<>"*"]]; (* taken from CustomTicks package *)


(* Usage messages *)

(* first introduce all the functions from Nat. Commun. paper *)
(* Not@ValueQ is taken from FiniteFields package: C:\Program Files\Wolfram Research\Mathematica\12.0\AddOns\Packages\FiniteFields *)
(* for Functions *)
If[Not@ValueQ[Hamiltonian4TBpack::usage],Hamiltonian4TBpack::usage="Constructs tight-binding Hamiltonian from the atomic coordinates of the unit cell, nearest-neighbour distances and hopping integrals, etc."];

If[Not@ValueQ[ElectronicStructure4TBpack::usage],ElectronicStructure4TBpack::usage="Calculates electronic energy levels or band structure for the specified monoelement unit cell, translation vectors and tight-binding parameters, etc."];

(* for Options *)


Begin["`Private`"]


Hamiltonian4TBpack[unitcell_,translationvectors_,hoppingintegrals_,overlappingintegrals_,nnlengthlist_?(VectorQ[#,NumericQ]&),\[Beta]_(* strain *),k_,\[CapitalEpsilon]_(* electric field *),B_(* magnetic field *),edgecorrections_,gflag_(* gauge flag is not used any more *)]:=Catch[
Module[
{
tf1,tf2,
Tvec,len,
aunitcells,
cnlist,
aplist,counter,v,
bondvector,bondlength,\[CurlyPhi],
straining,hij,sij,\[CurlyEpsilon],
H,S,
\[Delta]=0.021(* or 0.05*),(* tolerance for the bond length deviation from ideal reference value: 
See for details `Method' section of V.A.Saroka,K.G.Batrakov,and L.A.Chernozatonskii,Phys.Solid State 56,2135 (2014).
free access link: https://drive.google.com/file/d/0B1PX-Ri2Z2h1SjM5NnZ0ZjhXNzA/view
*)
lim,r1,r2,r12,
t,\[Delta]t,tE,tSO,
s,

\[CapitalPhi],\[CurlyPhi]M,
\[CapitalPhi]0=2.067833758 10^5(* flux quantum T \[Angstrom]^2*),

(* subfunctions *)
U,
conjugate,

arg1,arg2,inconsarg
},

Hamiltonian4TBpack::arg1="The argument `1` should be a list of points in 3D space.";Hamiltonian4TBpack::arg2="The argument `1` should be a list numeric vectors of 3D space or a list of lists of such vectors up to 3 ones.";
Hamiltonian4TBpack::inconsarg="Arguments `1` , `2` must have the same lengths which is less than length of `3`.";

tf1=MatchQ[#1,{arg___?(VectorQ[#1,NumericQ]&&Length[#1]==3&)}]&;tf2=(MatchQ[#1,{arg___?(VectorQ[#1,NumericQ]&&Length[#1]==3&)}]&&Length[#1]<=3)&;

If[!And@@(tf1/@unitcell),Throw[Message[Hamiltonian4TBpack::arg1,unitcell]]];
If[!And@@(tf2/@translationvectors),Throw[Message[Hamiltonian4TBpack::arg2,translationvectors]]];
If[!(Length@hoppingintegrals==Length@overlappingintegrals<=Length@nnlengthlist),Throw[Message[Hamiltonian4TBpack::inconsarg,hoppingintegrals,overlappingintegrals,nnlengthlist]]];

(*-------------------- subfunctions --------------------*)
(* Electrostatic potential *)
U[r_,e_]:=-r.e;  (* Homogenious electric field *)

(* complex conjugation *)
conjugate[expr_]:=expr/.Complex[x_,y_]:> x-I y;


(*------------------------------------------------------*)

(* first argument structure: unitcell={idealunitcell,optimizedunitcell};
   second argument structure: translationvectors={idealtranslationvectors,optimizedtranslationvectors} *)
Tvec=translationvectors;
len=Dimensions[unitcell][[2]];

(* Adjacent unitcells *)
aunitcells=Table[
	Switch[
	(Tvec//Dimensions)[[2]],
	1,Table[(#1+i Tvec[[o,1]]&)/@(unitcell[[o]]),{i,-1,1}],
	2,Flatten[Table[(#1+i Tvec[[o,1]]+j Tvec[[o,2]]&)/@(unitcell[[o]]),{i,-1,1},{j,-1,1}],1],
	3,Flatten[Table[(#1+i Tvec[[o,1]]+j Tvec[[o,2]]+k Tvec[[o,3]]&)/@(unitcell[[o]]),{i,-1,1},{j,-1,1},{k,-1,1}],2]]
,{o,1,2}](* end Table *);


(* Coordination numbers list generation *)
aplist=Flatten[aunitcells[[1]],1];(* atom positions list *)

cnlist=Table[
counter=0;
Do[
	v=unitcell[[1,i]]-aplist[[j]];
	If[Abs[Sqrt[v.v]-nnlengthlist[[2]]]<\[Delta],counter++],
	{j,1,Length@aplist}
];(* end Do *)
counter,
{i,1,len}];(* end Table *)

(* Hamiltonian and Overlap matrices *)
H=ConstantArray[0,{len,len}];
S=ConstantArray[0,{len,len}];

(* Hamiltonian and Overlap matrix elements filling *)
Do[
	hij=0;
	sij=0;
	(* Summation over adjecent unitcells *)
	Do[
		r1=unitcell[[1,i]];
		r2=aunitcells[[1,c,j]];
		r12=r1-r2; (* bondvector ideal *)
		(* Summation over nearest neighbours of various orders *)
		Do[
			If[
				Abs[nnlengthlist[[l]]-Sqrt[r12.r12]]<\[Delta],
				r1=unitcell[[2,i]];
				r2=aunitcells[[2,c,j]];
				
				t=hoppingintegrals[[l]];
				s=overlappingintegrals[[l]];

				r12=r1-r2; (* bondvector real *)
				bondlength=Sqrt[r12.r12];
				\[CurlyEpsilon]=If[l==1,0,(bondlength-nnlengthlist[[l]])/nnlengthlist[[l]]];
(* straining is taken into account as in the paper: Ribeiro,R.M.,Pereira,V.M.,Peres,N.M.R.,Briddon,P.R.,and Castro Neto,A.H. "Strained graphene:tight-binding and density functional calculations", New Journal of Physics,11(11),115002 (2009). http://doi.org/10.1088/1367-2630/11/11/115002 *)


				(* Generalization of the Landau gauge to the arbitrary direction of the B :
				 eq.(24) from J.-C.Charlier and S.Roche,Rev.Mod.Phys.79,677 (2007).*)
				If[B//NumberQ,
				\[CurlyPhi]M=B,

				\[CapitalPhi]=(r2+r1)[[1]]/2 B[[3]] (r2-r1)[[2]]+((r2+r1)[[2]]/2  B[[1]]-(r2+r1)[[1]]/2  B[[2]]) (r2-r1)[[3]];(* magnetic flux *)
				
				\[CurlyPhi]M=\[Pi] \[CapitalPhi]/\[CapitalPhi]0;(* magnetic phase *)
				];(* end If *)
				\[CurlyPhi]=k.r12;(* space phase *)

				\[Delta]t=Switch[
							cnlist[[i]],
							1,edgecorrections[[l,1]],
							2,edgecorrections[[l,2]],
							3,0,
							_,0
				];(* edge correction *)

				(* electric field *)
				tE=If[l==1,U[r1,\[CapitalEpsilon]],0];(* electrostatic on-site energy *)
			
								
				hij=hij+Exp[-\[Beta] \[CurlyEpsilon]] (t+\[Delta]t+tE) Exp[I (\[CurlyPhi]-\[CurlyPhi]M)];
				sij=sij+Exp[-\[Beta] \[CurlyEpsilon]] s Exp[I (\[CurlyPhi]-\[CurlyPhi]M)];
	
				Break[];
			](* end If *),
		{l,1,Length[hoppingintegrals]}
		](* end Do *)
		,
		{c,1,Length[aunitcells[[1]]]}
	](* end Do *);
	H[[i,j]]=If[i==j,1/2 hij,hij];
	S[[i,j]]=If[i==j,1/2 sij,sij];
	,
{i,1,len},
{j,i,len}(* to save time we fill in only upper triangle of the matrix *)
](* end Do *);
{H+conjugate@Transpose@H,S+conjugate@Transpose@S}
](* end Module *)
](* end Catch *);


ElectronicStructure4TBpack[unitcell_,translationvectors_,hoppingintegrals_?(VectorQ[#,NumericQ]&),overlappingintegrals_?(VectorQ[#,NumericQ]&),nnlengthlist_,\[Beta]_?NumericQ,klist_,\[CapitalEpsilon]_,B_,edgecorrections_,gflag_,probflag_]:=Catch[
Module[
{(*--------------------------- Variables and constants ---------------------------*)
tf,len,klen,
bands,m,
k,rng,
vm,kx,ky,kz,

data,
epsol,epsolsorted,
eval,evec,
hsm,hmsorted,
v,vsorted
},

ElectronicBands4TBpack::klist="The argument must be a list of vectors three element long each.";

tf=VectorQ[#,(ListQ[#]&&Length[#]==3&)]&;

If[!tf@klist,Throw[Message[ElectronicBands4TBpack::klist,unitcell]]];

(*--------------------------- Body of the fucntion ---------------------------------*)
If[
	!probflag,
	len=Dimensions[unitcell][[2]];
	klen=Length[klist];
	m=Hamiltonian4TBpack[unitcell,translationvectors,hoppingintegrals,overlappingintegrals,nnlengthlist,\[Beta],k,\[CapitalEpsilon],B,edgecorrections,gflag];
	bands=ParallelTable[
				Join[Sort@Chop[Eigenvalues[m/.k-> klist[[i]]]],{klist[[i]]}],
				{i,1,klen}
	](* end ParallelTable *);
	Transpose[bands]
	,
	len=Dimensions[unitcell][[2]];
	klen=Length[klist];
	m=Hamiltonian4TBpack[unitcell,translationvectors,hoppingintegrals,overlappingintegrals,nnlengthlist,\[Beta],k,\[CapitalEpsilon],B,edgecorrections,gflag];
	(*=velocity operator in the gradient (effective mass) approximation=*)
	vm=D[m[[1]]/.k->{kx,ky,kz},{{kx,ky,kz}}];
	(*==================================================================*)
	data=Table[
				hsm=m/.k->klist[[i]];
				epsol=Chop[Eigensystem[hsm]];
				{eval,evec}=#[[Ordering[epsol[[1]]]]]&/@epsol;				
			
				v=vm/.{kx->klist[[i,1]],ky->klist[[i,2]],kz->klist[[i,3]]};
				(*========================================================*)
				
				{Append[eval,klist[[i]]],Map[Conjugate[#1]*#1&,evec],evec,v}
				,
				{i,1,klen}
	](* end ParallelTable *);
	data=Transpose[data];
	Join[{Transpose[First@data]},Rest@data]
](* end If probflag *)
](* end Module *)
](* end Catch *);


End[]


Protect[Evaluate[$Context<>"*"]]; (* this line is used following CustomTicks package, MaTeX package also protects its symbols *)
EndPackage[]
