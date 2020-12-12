(* Wolfram Language Package *)

BeginPackage["TBpack`UnitcellGenerators`"]
(* Exported symbols added here with SymbolName::usage *)  


(* Functions *)
Nanotube::usage = "Nanotube[\*RowBox[{StyleBox[\"n\",\"TI\"], \",\" , StyleBox[\"m\",\"TI\"]}]] retuns a list of atomic coordinates of the nanotube, its translation vector and the lattice constant used.";
Nanoribbon::usage = "Nanoribbon[\*StyleBox[\"n\",\"TI\"]] returns a list of atomic coordinate of a \*StyleBox[\"Zigzag\",\"TI\"] nanoribbon, its translation vector and the lattice constant used.
Nanoribbon[\*RowBox[{StyleBox[\"n\",\"TI\"], \",\" , StyleBox[\"options\",\"TI\"]}]] generates atomic coordinates with specified options.
Nanoribbon[\*StyleBox[\"n\",\"TI\"], EdgeType \[Rule] value] generates \*StyleBox[\"Zigzag\",\"TI\"], \*StyleBox[\"Armchair\",\"TI\"] and \*StyleBox[\"Bearded\",\"TI\"] nanoribbon for \*StyleBox[\"value\",\"TI\"] equal to 1, 2 and 3, respectively"; 
CNanoribbon::usage = "CNanoribbon[\*RowBox[{StyleBox[\"n\",\"TI\"], \",\", StyleBox[\"m\",\"TI\"]}]] returns a list of atomic coordinate of a nanoribbon by effectively unrolling a nanotube with the chirality vector \*RowBox[{\"(\", StyleBox[\"n\",\"TI\"], \",\", StyleBox[\"m\",\"TI\"] ,\")\"}].
CNanoribbon[\*RowBox[{StyleBox[\"n\",\"TI\"], \",\", StyleBox[\"m\",\"TI\"]}], RefinedEdge \[Rule] \*StyleBox[\"True\",\"TI\"]] refines the edge of the ribbon from the dangling atoms.";
ZigzagShapedNanoribbon::usage = "ZigzagShapedNanoribbon[\*RowBox[{StyleBox[SubscriptBox[\"l\",\"1\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"l\",\"2\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"w\",\"1\"],\"TI\"], \",\" , StyleBox[SubscriptBox[\"w\",\"2\"],\"TI\"] }]] returns a list of atomic coordinates of a \*StyleBox[\"Z60\",\"TI\"] zigzag-shaped nanoribbon, its translation vector and the lattice constant used.
ZigzagShapedNanoribbon[\*RowBox[{StyleBox[SubscriptBox[\"l\",\"1\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"l\",\"2\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"w\",\"1\"],\"TI\"], \",\" , StyleBox[SubscriptBox[\"w\",\"2\"],\"TI\"] }], \*StyleBox[\"options\",\"TI\"]] generates atomic coordinates with the specified \*StyleBox[\"options\",\"TI\"].
ZigzagShapedNanoribbon[\*RowBox[{StyleBox[SubscriptBox[\"l\",\"1\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"l\",\"2\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"w\",\"1\"],\"TI\"], \",\" , StyleBox[SubscriptBox[\"w\",\"2\"],\"TI\"] }], EdgeType \[Rule] \*StyleBox[\"value\",\"TI\"]] generates \*StyleBox[\"Z60\",\"TI\"], \*StyleBox[\"Z120\",\"TI\"], \*StyleBox[\"A60\",\"TI\"] and \*StyleBox[\"A120\",\"TI\"] nanoribbons for \*StyleBox[\"value\",\"TI\"] equal to 1, 2, 3, and 4, respectively.
ZigzagShapedNanoribbon[\*RowBox[{StyleBox[SubscriptBox[\"l\",\"1\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"l\",\"2\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"w\",\"1\"],\"TI\"], \",\" , StyleBox[SubscriptBox[\"w\",\"2\"],\"TI\"] }], TranslationAxis \[Rule] \*StyleBox[\"value\",\"TI\"]] sets the nanoribbon translation vector along \*StyleBox[\"Ox\",\"TI\"]- or \*StyleBox[\"Oy\",\"TI\"]-axis  for \*StyleBox[\"value\",\"TI\"] equal to 1 or 2, respectively.
ZigzagShapedNanoribbon[\*RowBox[{StyleBox[SubscriptBox[\"l\",\"1\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"l\",\"2\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"w\",\"1\"],\"TI\"], \",\" , StyleBox[SubscriptBox[\"w\",\"2\"],\"TI\"] }], ApexPoint \[Rule] \*StyleBox[\"integer\",\"TI\"]] sets the apex point to one of the four positions denoted by the \*StyleBox[\"integer\",\"TI\"].";
IVGroupQuantumDot::usage = "IVGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}]] returns a list of atomic coordinates for a graphene quantum dot, an effective translation vector (roughly double the size of the dot) and the lattice constant.
IVGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], \*StyleBox[\"options\",\"TI\"]] generates a list of atomic coordinates for a quantum dot with the specified options.
IVGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], ChemicalElement \[Rule] \*StyleBox[\"label\",\"TI\"]] returns a list of atomic coordinates for a quantum dot of the 2D material made of the element set by \*StyleBox[\"label\",\"TI\"].
IVGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], EdgeType \[Rule] \*StyleBox[\"value\",\"TI\"]] returns a list of atomic coordinates for a quantum dot with \*StyleBox[\"Zigzag\",\"TI\"] and \*StyleBox[\"Armchair\",\"TI\"] edges for \*StyleBox[\"value\",\"TI\"] equal to 1 and 2, respectively.
IVGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], Shape \[Rule] \*StyleBox[\"value\",\"TI\"]] returns a list of atomic coordinates for a quantum dot with \*StyleBox[\"triangular\",\"TI\"] and \*StyleBox[\"hexagonal\",\"TI\"] shapes for \*StyleBox[\"value\",\"TI\"] equal to \"TRI\" and \"HEX\", respectively.";
VGroupQuantumDot::usage = "VGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}]] returns a list of atomic coordinates for a phosphorene quantum dot, an effective translation vector (roughly double the size of the dot) and the lattice constant.
VGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], \*StyleBox[\"options\",\"TI\"]] generates a list of atomic coordinates for a quantum dot with the specified options.
VGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], ChemicalElement \[Rule] \*StyleBox[\"label\",\"TI\"]] returns a list of atomic coordinates for a quantum dot of the 2D material made of the element set by \*StyleBox[\"label\",\"TI\"].
VGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], EdgeType \[Rule] \*StyleBox[\"value\",\"TI\"]] returns a list of atomic coordinates for a quantum dot with \*StyleBox[\"Zigzag\",\"TI\"], \*StyleBox[\"Armchair\",\"TI\"] and \*StyleBox[\"Mixed\",\"TI\"] edges for \*StyleBox[\"value\",\"TI\"] equal to 1, 2 and 3, respectively.
VGroupQuantumDot[\*RowBox[{StyleBox[\"size\",\"TI\"], \",\" , StyleBox[\"numberoflayers\",\"TI\"]}], Shape \[Rule] \*StyleBox[\"value\",\"TI\"]] returns a list of atomic coordinates for a quantum dot with \*StyleBox[\"triangular\",\"TI\"], \*StyleBox[\"hexagonal\",\"TI\"] and \*StyleBox[\"circular\",\"TI\"] shapes for \*StyleBox[\"value\",\"TI\"] equal to \"TRI\", \"HEX\" and \"CIR\", respectively. \"(r)\" added to any shape string radomizes the edges with Koch curves.";



(* Options *)
LatticeConstant::usage = "Option specifying the scale of the 2D hexagonal lattice in \[Angstrom].";
EdgeType::usage = "Option specifying the type of the edge in some unit cell generators of nanostructures based on the 2D hexagonal lattice.";
RefinedEdge::usage = "Option taking values \*StyleBox[\"True\",\"TI\"] or \*StyleBox[\"False\",\"TI\"] specifying if the nanoribbon edge should be refined from dangling atoms.";
TranslationAxis::usage = "Option specifying the direction of translation in some unit cell generators of 1D structures.";
ApexPoint::usage = "Option specifying the arangement of a nanoribbon within the 2D hexagonal lattice, see V. A. Saroka and K. G. Batrakov, Russ. Phys. J. 59, 633 (2016).";
ChemicalElement::usage = "Option specifying the label of the element from the Periodic Table.";
Shape::usage = "Option specifying the shape of nanostructures in some unit cell generators.";


(* Constants *)
Zigzag::usage = "Constant equal to 1 setting the EdgeType option value.";
Armchair::usage = "Constant equal to 2 setting the EdgeType option value.";
Bearded::usage = "Constant equal to 3 setting the EdgeType option value.";
Mixed::usage = "Constant equal to 3 setting the EdgeType option value in \*StyleBox[\"VGroupQuantumDots\",\"TI\"]";
Ox::usage = "Constant equal to 1 setting the TranslationAxis option value.";
Oy::usage = "Constant equal to 2 setting the TranslationAxis option value.";
Oz::usage = "Constant equal to 3 setting the TranslationAxis option value.";
Z60::usage = "Constant equal to 1 and setting the zigzag-shaped nanoribbon superlattice type to \*StyleBox[\"Z60\",\"TI\"].";
Z120::usage = "Constant equal to 2 and setting the zigzag-shaped nanoribbon superlattice type to \*StyleBox[\"Z120\",\"TI\"].";
A60::usage = "Constant equal to 3 and setting the zigzag-shaped nanoribbon superlattice type to \*StyleBox[\"A60\",\"TI\"].";
A120::usage = "Constant equal to 4 and setting the zigzag-shaped nanoribbon superlattice type to \*StyleBox[\"A120\",\"TI\"].";


Begin["`Private`"] (* Begin Private Context *)
 
Options[Nanotube] = {
	LatticeConstant -> 1.42
	};
Nanotube[n_Integer, m_Integer, OptionsPattern[]]:=Catch[
Module[
{
(*-------------- functions ----------------*)
	RollingUp,
	CylindricalToCartesian,
	qfunc,lowqfunc,upqfunc,
(*-------------- variables ----------------*)
	a0 = OptionValue[LatticeConstant],
	i,j,
	a,a1,a2,
	L,r,Cc,theta,
	dR,T,t1,t2,
	NumberOfHexagons,
	R,p,q,
	lim,pqlist,
	initPoint,Shiftx,
	del1,del2,del3,
	Coord,LineCoord,V,
	max,u,cu,unew,dz

},
(* Initial parameters *)

RollingUp[list_]:={list[[1]],list[[2]]/r,list[[3]]};
CylindricalToCartesian[list_]:={list[[1]]*Cos[list[[2]]],list[[1]]*Sin[list[[2]]],list[[3]]};

a = a0 Sqrt[3];

(* primitive vectors of the 2D hexagonal lattice *)
a1 = {0,a Cos[Pi/6],a Sin[Pi/6]};
a2 = {0,a Cos[Pi/6],-a Sin[Pi/6]};

(* nanotube radius *)
L = a Sqrt[n^2 + m^2 + n m];
r = L/(2 Pi);

(* chirality vector of the nanotube *)
Cc = n a1 + m a2 + {r,0,0};

(* rotation of the reference system*)
theta = ArcTan[Cc[[3]]/Cc[[2]]];
a1 = RotationMatrix[-theta,{1,0,0}].a1;
a2 = RotationMatrix[-theta,{1,0,0}].a2;


(* translation vector of a nanotube *)
dR = GCD[2 n + m,2 m + n];
t1 = (2 m + n)/dR;
t2 = -(2 n + m)/dR;
T = t1 a1 + t2 a2;

NumberOfHexagons = (2 L^2)/(a^2 dR);
(* symmetry vector of the nanotube *)
If[
	t1 == 0,
	p = 0;
	q = -1;
	,
	lim = 100;
	qfunc[p_]:=(1 + p t2)/t1;
	upqfunc[p_]:=(m p)/n;
	lowqfunc[p_]:=(-NumberOfHexagons+m p)/n;
	pqlist = {};
	Do[If[IntegerQ[qfunc[p]],AppendTo[pqlist,{p,qfunc[p]}]],{p,-lim,lim}];
	Do[
		If[
			(lowqfunc[pqlist[[i,1]]]<pqlist[[i,2]]<upqfunc[pqlist[[i,1]]]),
			p = pqlist[[i,1]];	 
			q = pqlist[[i,2]];
		](* end If *),
		{i,Length[pqlist]}
	](* end Do *)
](* end If *);
R = p a1 + q a2;


(* shift of sublattices in 2D hexagonal lattice *)
initPoint = RotationMatrix[-theta,{1,0,0}].{0,a0,0};
(* shift from the origin of the coordinate system*)
Shiftx = {r,0,0};

(* vectors of dirrection to nearest neighbours *)
del1 = a0 RotationMatrix[-theta,{1,0,0}].{0,Cos[Pi/3],Sin[Pi/3]};
del2 = a0 RotationMatrix[-theta,{1,0,0}].{0,Cos[Pi/3],-Sin[Pi/3]};
del3 = a0 RotationMatrix[-theta,{1,0,0}].{0,-1,0};

(* list of coordinates *)
Coord = {};
(* list of bonds *)
LineCoord = {};
(* fill lists for the first sublattice *)
For[i = 1; j = 0; V = {0,0}, i <= NumberOfHexagons, i++,
	V = i R - j T + Shiftx;

	If[
		V[[3]] >= -(T[[2]]/T[[3]]) V[[2]]+(T[[2]]^2+T[[3]]^2)/T[[3]],
		j = j + 1;
		V = i R - j T + Shiftx
	];(* end If *)

	AppendTo[Coord,V];
	AppendTo[LineCoord,{V,V-del1}];
	AppendTo[LineCoord,{V,V-del2}];
	AppendTo[LineCoord,{V,V-del3}]
];(* end For i *)

(* Fill lists for the second sublattice *)
For[i = 1; j = 0; V = initPoint,i <= NumberOfHexagons, i++,
	V = i R - j T + initPoint + Shiftx;
	
	If[
		V[[3]] >= -(T[[2]]/T[[3]]) V[[2]]+(T[[2]]^2+T[[3]]^2)/T[[3]],
		j = j + 1; V = i R - j T + initPoint + Shiftx
	];(* end If *)
	
	AppendTo[Coord,V];
	AppendTo[LineCoord,{V,V+del1}];
	AppendTo[LineCoord,{V,V+del2}];
	AppendTo[LineCoord,{V,V+del3}]
];(* end For i *)

Coord = Map[RollingUp,Coord];
Coord = Map[CylindricalToCartesian,Coord];


If[
	m == 0,
	max = Max /@ Transpose[Coord];(* maxima of the atomic coordinates *)
	dz = 0.1 a0; (* delta step *)
	u = Select[Coord,  (#[[3]] >= max[[3]] - dz) &];(* selection of the edge atoms not further than dz from the "max" edge in z-direction *)
	cu = Complement[Coord,u]; (* atomic coordinates of the rest of the atoms *)
	unew = # - T&/@u; (* shift of the edge atoms onto opposite side by substracting translation vector T*)
	Coord = Join[cu,unew];
	Return[Chop@{Coord,T,a0}],
	Return[Chop@{Coord,T,a0}]
](* end If *)
](* end Module *);
](* Catch braket *)(* end Function Nanotube *)
SyntaxInformation[Nanotube] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};


Zigzag = 1;
Armchair = 2;
Bearded = 3;
Ox = 1;
Oy = 2;
Oz = 3;

Options[Nanoribbon] = {
	LatticeConstant -> 1.42,
	EdgeType -> Zigzag,
	TranslationAxis -> Oy
	};
Nanoribbon[numberofchains_Integer, OptionsPattern[]] := Module[
{
	nanoribbontype = OptionValue[EdgeType],
	a0 = OptionValue[LatticeConstant],
	dir = OptionValue[TranslationAxis],
	
	a, vec, p1, p2, p,
	tv, unitcell, T,
	pn,	rm,

	a1, a2
},
a = Sqrt[3] a0;
vec = {a0 {-(1/2),Sqrt[3]/2,0},a0 {1/2,Sqrt[3]/2,0}};
(* lattice primitive translations *)
a1 = a {1/2,Sqrt[3]/2,0};
a2 = a {1,0,0};

Switch[
		nanoribbontype,
		1,
		p1 = {{0,0,0},{0,0,0} + vec[[1]]};
		p2 = {{a0,0,0},{a0,0,0} + vec[[2]]};
		p = {p1,p2};
		tv = {3 a0,0,0};
		T = a {0,1,0},
		2,
		p1 = {{0,0,0},{a0, 0,0}};
		p2 = {{0,0,0} + vec[[1]],{a0,0,0} + vec[[2]]};
		p = {p1,p2};
		tv = a {0,1,0};
		T = {3 a0,0,0},
		3,
		p1 = {{0,0,0},{0,a0,0}};
		p2 = {{0,0,0} + a1,{0,a0,0} + a1};
		p = {p1,p2};
		tv = {0,3 a0,0};
		T = a2
](* end Switch *);
unitcell = Flatten[Table[(#+{(n-1)/2,(n-2)/2}[[Mod[n,2,1]]]tv)&/@p[[Mod[n,2,1]]] ,{n,numberofchains}],1];

(* a normal vector to the plane defined by T and Ox-, Oy- or Oz-axis *)
pn = Cross[T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]];
If[
	Sqrt[pn.pn] =!= 0,
	(* if the normal vector well defined then reorient the structure so that translation vector is directed along the chosen axis *)
	rm = If[
				Head[a0]===Symbol, 
				Simplify[RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}],a0>0],
				RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}]
	];
	unitcell = (rm.#) & /@ unitcell;
	T = rm.T
];

{unitcell,T,a0}

](* end Module *);
SyntaxInformation[Nanoribbon] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};


(* Chiral nanoribbon atomic coordinates by nanotube unrolling. Achiral nanoribbons can be also generated. *)
Options[CNanoribbon] = {
	LatticeConstant->1.42,
	RefinedEdge->True,
	TranslationAxis -> Oy
	};
CNanoribbon[n_Integer, m_Integer, OptionsPattern[]] := Module[
{
	(*-------------- functions ----------------*)
	qfunc,lowqfunc,upqfunc,
	(*-------------- variables ----------------*)
	a0 = OptionValue[LatticeConstant],
	edgeflag = OptionValue[RefinedEdge],
	dir = OptionValue[TranslationAxis],
	
	i,j,
	a,a1,a2,
	L,Cc,theta,
	dR,T,t1,t2,
	numberofhexagons,
	R,p,q,
	lim,pqlist,
	initPoint,
	coord,v,

	adjacentunitcells,
	counter,
	
	max,dy,u,cu,unew,
	pn, rm

},

(* parameters *)
a = a0 Sqrt[3];
 
(* Primitive vectors of the 2D hexagonal lattice *)
a1 = {a Cos[Pi/6],a Sin[Pi/6],0};
a2 = {a Cos[Pi/6],-a Sin[Pi/6],0};

(* nanotube radius *)
L = a Sqrt[n^2 + m^2 + n m];

(* chirality vector *)
Cc = n a1 + m a2;

(* rotation of the reference system *)
theta = ArcTan[Cc[[2]]/Cc[[1]]];
a1 = RotationMatrix[-theta,{0,0,1}].a1;
a2 = RotationMatrix[-theta,{0,0,1}].a2;


(* translation vector of a nanotube *)
dR = GCD[2 n + m,2 m + n];
t1 = (2 m + n)/dR;
t2 = -(2 n + m)/dR;
T = t1 a1 + t2 a2;

numberofhexagons = (2 L^2)/(a^2 dR);
(* Symmetry vector *)
If[
	t1 == 0,
	p = 0;
	q = -1;
	,
	lim = 100;
	qfunc[p_]:=(1 + p t2)/t1;
	upqfunc[p_]:=(m p)/n;
	lowqfunc[p_]:=(-numberofhexagons + m p)/n;
	pqlist = {};
	Do[If[IntegerQ[qfunc[p]],AppendTo[pqlist,{p,qfunc[p]}]],{p,-lim,lim}];
	Do[
		If[
			(lowqfunc[pqlist[[i,1]]]<pqlist[[i,2]]<upqfunc[pqlist[[i,1]]]),
			p = pqlist[[i,1]];	 
			q = pqlist[[i,2]];
		](* end If *),
		{i,Length[pqlist]}
	](* end Do *)
](* end If *);
R = p a1 + q a2;

(* shift of sublattices in graphene *)
initPoint = RotationMatrix[-theta,{0,0,1}].{a0,0,0};

(* list of coordinates *)
coord={};
(* generate a list for the first sublattice *)
For[i = 1; j = 0; v = {0,0,0}, i<= numberofhexagons, i++,
	v = i R - j T;
	If[
		v[[2]] >= -(T[[1]]/T[[2]]) v[[1]]+(T[[1]]^2+T[[2]]^2)/T[[2]],
		j = j + 1;
		v = i R - j T
	];(* end If *)
	AppendTo[coord,v]
];(* end For i *)

(* generate a list for the second sublattice *)
For[i = 1; j = 0; v = initPoint, i <= numberofhexagons, i++,
	v = i R - j T + initPoint;
	If[
		v[[2]] >= -(T[[1]]/T[[2]]) v[[1]]+(T[[1]]^2+T[[2]]^2)/T[[2]],
		j = j + 1;
		v = i R - j T + initPoint		
	];(* end If *)
	AppendTo[coord,v]
];(* end For i *)

If[
	edgeflag,
	adjacentunitcells = Flatten[Table[(# + i T)&/@coord,{i,-1,1}],1];
	coord = Reap[
		Do[
			counter=0;
			Do[
				v = coord[[i]] - adjacentunitcells[[j]];
				If[Abs[Sqrt[v.v] - a0] < 0.01 a0 , counter++];
				If[counter > 2,Break[]],
			{j,Length[adjacentunitcells]}];
			If[counter == 1,Sow[coord[[i]],"edge"]];
			If[counter > 1,Sow[coord[[i]]],"bulk"],
		{i,Length[coord]}]
	][[2]];
	
	(* correct the unit cell of armchair ribbons *)
	If[
		m == 0,
		max = Max /@ Transpose[coord[[2]]];(* maxima of the atomic coordinates *)
		dy = 0.2 a0; (* delta step *)
		u = Select[coord[[2]], (#[[2]] >= max[[2]] - dy) &];(* selection of the edge atoms not further than dz from the "max" edge in z-direction *)
		cu = Complement[coord[[2]],u]; (* atomic coordinates of the rest of the atoms *)
		unew = # - T&/@u; (* shift of the edge atoms onto opposite side by substracting translation vector T *)
		coord[[2]] = Join[cu,unew]
	];(* end If *)
	
	(* a normal vector to the plane defined by T and Ox-, Oy-, or Oz-axis *)
	pn = Cross[T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]];
	If[
		Sqrt[pn.pn] =!= 0,
		(* if the normal vector well defined then reorient the structure so that translation vector is directed along the chosen axis *)
		rm = If[
				Head[a0]===Symbol, 
				Simplify[RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}],a0>0],
				RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}]
		];
		coord = {(rm.#) & /@ coord[[1]], (rm.#) & /@ coord[[2]]};
		T = rm.T
	];
	
	Chop@{coord[[2]],T,a0,(* edge atoms *)coord[[1]]}
	,
	
	(* a normal vector to the plane defined by T and Ox-, Oy- or Oz-axis *)
	pn = Cross[T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]];
	If[
		Sqrt[pn.pn] =!= 0,
		(* if the normal vector well defined then reorient the structure so that translation vector is directed along the chosen axis *)
		rm = If[
				Head[a0]===Symbol, 
				Simplify[RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}],a0>0],
				RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}]
		];
		coord = (rm.#) & /@ coord;
		T = rm.T
	];
	
	Chop@{coord, T, a0}
]

](* end Module *);
SyntaxInformation[CNanoribbon] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};


(* Zigzag-shaped nanoribbon superlattices *)
Z60 = 1;
Z120 = 2;
A60 = 3;
A120 = 4;
Options[ZigzagShapedNanoribbon] = {
	LatticeConstant -> 1.42,
	EdgeType -> Z60,
	TranslationAxis -> Oy,
	ApexPoint -> 1
	};
ZigzagShapedNanoribbon[l1_Integer, l2_Integer, w1_Integer, w2_Integer, OptionsPattern[]] := Module[
{
  	rtype = OptionValue[EdgeType],
  	dir = OptionValue[TranslationAxis],
  	apexpoint = OptionValue[ApexPoint],
   	a0 = OptionValue[LatticeConstant], 
   	a, a1, a2,
   	l1c, l2c, L1c, L2c, Wc,
   	data, UnitCell, T, pn, rm
},

a = Sqrt[3] a0;
a1 = a {1/2, Sqrt[3]/2, 0};
a2 = a {1, 0, 0};
  
l1c = {{1, 0}, {1, -1}, {2, -1}, {2, -1}, {1, 0}, {1, 0}, {1, 
      0}, {1, 0}}[[rtype]];
l2c = {{0, 1}, {0, 1}, {1, 1}, {-1, 2}, {0, 1}, {0, 1}, {0, 1}, {0, 
      1}}[[rtype]];
(* sv = {{0, 0, 0}, (a1 + a2)/3, a0/2 {Sqrt[3]/2, -(1/2), 0}, 
     a0/2 {-(Sqrt[3]/2), -(1/2), 0}}[[smp]];(* shift vector setting the symmetry point *) *)

(*sv = {(a1 + a2)/3, -(a1/6) + (a2/3), -(a1 + a2)/6, {0, 0, 0}}[[apexpoint]];(* shift vector setting the apex point *)*)
  
L1c = l1 l1c;
L2c = l2 l2c;
  
Wc = If[rtype == 3 || rtype == 4, 1/3, 1] (w1 l1c + w2 l2c); (* width vector coordinates expressed by means of l1c, l2c basis *)
 
data = ZSNR[L1c[[1]], L1c[[2]], L2c[[1]], L2c[[2]], Wc[[1]], Wc[[2]], LatticeConstant -> a0, ApexPoint -> apexpoint];
UnitCell = data[[1]];
T = data[[2]];
  
(* a normal vector to the plane defined by T and Ox-, Oy- or Oz-axis *)
pn = Cross[T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]];
If[
	Sqrt[pn.pn] =!= 0,
	(* if the normal vector well defined then reorient the structure so that translation vector is directed along the chosen axis *)
	rm = RotationMatrix[{T, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}[[dir]]}];
	UnitCell = (rm.#) & /@ UnitCell;
	T = rm.T
];

Return[{UnitCell, T, a0}]
](* end Module *);
SyntaxInformation[ZigzagShapedNanoribbon] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};
  

Options[ZSNR] = {
	LatticeConstant -> 1.42,
	ApexPoint -> 1
};
ZSNR[l11_Integer, l12_Integer, l21_Integer, l22_Integer, w1_, w2_, OptionsPattern[]] := Catch[Module[
   {
    a0 = OptionValue[LatticeConstant],
    apexpoint = OptionValue[ApexPoint],
    
    \[Delta] = 0.05,
    a, a1, a2, L1, L2, W,
    lim, shift,
    A, B, graphenesheet,
    InCell, RPQ, l, r, NNQ,
    col, vec, d,
    t, lbl, dd, g
    },
   (* Errors *)
   ZSNR::veccoord = "The coordinates (`1`,`2`) are irrelevant";
   
   RPQ[p_, r0_, v_, \[Delta]_] := Module[
     {
      y
      },
     (* r=r0+nv *)
     If[
      v[[1]] != 0,
      y[x_] := v[[2]]/v[[1]] x + (r0[[2]] - v[[2]]/v[[1]] r0[[1]]);
      If[
       	v[[1]] > 0,
       	y[p[[1]]] - \[Delta] >= p[[2]],
       	y[p[[1]]] + \[Delta] <= p[[2]]
       ](* end If *),
      If[
       	v[[2]] > 0,
       	r0[[1]] + \[Delta] <= p[[1]],
       	r0[[1]] - \[Delta] >= p[[1]]
       ](* end If *)
      ](* end If *)
     ];(* end Module *)
   
   
   InCell[point_, L1_, L2_, W_] := Module[
     {
      ip, r01, v1, r02, v2, pol,
      list
      },
     ip = {0, 0, 0};
     r01 = {ip, L1 + W, W, L1};
     v1 = {L1, -L1, -W, W};
     r02 = {ip, W, L2 + W, L2};
     v2 = {W, L2, -W, -L2};
     pol = Transpose /@ {{r01, v1}, {r02, v2}}; (* polygons *)
     r = RPQ[point, #[[1]], #[[2]], \[Delta]] &;
     l = RPQ[point, #[[1]], #[[2]], -\[Delta]] &;
     list = {Flatten[{l /@ Most[pol[[1]]], r@Last[pol[[1]]]}], 
       l /@ pol[[2]]};
     Or @@ (And @@ # & /@ list)
     ];(* end Module *)
   
   NNQ[unitcell_, T_] := Module[
     {
      \[CapitalDelta] = 0.01,
      len, m, tfv
      },
     len = Length@unitcell;
     m = Table[
       unitcell[[i]] - (c T + unitcell[[j]]), {i, len}, {j, len}, {c, -1, 1}];
     tfv = 
      Or @@ # & /@ 
       Map[Or @@ (Abs[Sqrt[#.#] - a0] < \[CapitalDelta] & /@ #) &, 
        m, {2}];
     Reap[
       Scan[If[#[[1]], Sow[#[[2]]]] &, Transpose@{tfv, unitcell}]][[2,1]]
     ];(* end Module *)
   
   a = Sqrt[3] a0;
   a1 = a {1/2, Sqrt[3]/2, 0};
   a2 = a {1, 0, 0};
   
   L1 = l11 a1 + l12 a2;
   L2 = l21 a1 + l22 a2;
   W = w1 a1 + w2 a2;
   
   
   If[
    	Cross[L1, W].{0, 0, -1} < 0,
    	(*-------------------------*)
    	t[s_?StringQ] := Style[s, FontSize -> 20, FontFamily -> "Times"];
    	col = {Red, Green, Blue};
    	vec = {L1, W, L2};
    	lbl = 
     t /@ {"\!\(\*SubscriptBox[OverscriptBox[\(L\), \
\(\[RightVector]\)], \(1\)]\)", 
       "\!\(\*OverscriptBox[\(W\), \(\[RightVector]\)]\)", 
       "\!\(\*SubscriptBox[OverscriptBox[\(L\), \(\[RightVector]\)], \
\(2\)]\)"};
    	d = {col, vec}\[Transpose];
    	dd = {lbl, vec}\[Transpose];
    	g = Graphics[{{#[[1]], Arrow[{{0, 0}, Most@#[[2]]}]} & /@ d, 
       Text[#[[1]], Most@#[[2]], {-1, 0}] & /@ dd}];
    	(*-------------------------*)
    	Throw[Print@g; Message[ZSNR::veccoord, w1, w2]],
    	If[
     		Cross[L2, W].{0, 0, 1} < 0,
     		Throw[Print@g; Message[ZSNR::veccoord, w1, w2]],
     		If[
      			Cross[L1, L2].{0, 0, -1} < 0,
      			Throw[Print@g; Message[ZSNR::veccoord, l11, l12]]
      		](* end If *)
     	](* end If *)
    ](* end If *);
   
   (* shift = {{0, 0, 0}, (a1 + a2)/3, a0/2 {Sqrt[3]/2, -(1/2), 0}, 
      a0/2 {-(Sqrt[3]/2), -(1/2), 0}}[[smp]]; *)
   shift = {(a1 + a2)/3, -(a1/6) + (a2/3), -(a1 + a2)/6, {0, 0, 0}}[[apexpoint]];
   
   lim = Ceiling[Max[1.2 Sqrt[#.#] & /@ {L1 + W, L2 + W}]/a];
   A = Table[n a1 + m a2, {n, -lim, lim}, {m, -lim, lim}];
   B = Table[n a1 + m a2 + (a1 + a2)/3, {n, -lim, lim}, {m, -lim, lim}];
   graphenesheet = # + shift & /@ Flatten[{A, B}, 2];
   
   {NNQ[Select[graphenesheet, InCell[#, L1, L2, W] &], L1 - L2], 
    L1 - L2, a0}
   ](* end Module *)
  ](* end Catch *);
SyntaxInformation[ZSNR] = {"ArgumentsPattern" -> {_, _, _, _, _, _, OptionsPattern[]}};


(* Quantum dots group IV added on 26/07/2020 *)
(* original code from 24/04/2016 *)
(* Error messages *)
IVGroupQuantumDot::argval = "The argument `1` must be a positive non-zero integer.";

Options[IVGroupQuantumDot] = {
	ChemicalElement -> "C",
	EdgeType -> Zigzag,
	Shape -> "TRI"
};
IVGroupQuantumDot[size_Integer, numberoflayers_Integer, OptionsPattern[]] := Catch[
  Module[
   {
   	element = OptionValue[ChemicalElement],
   	edgetype = OptionValue[EdgeType],
   	shape = OptionValue[Shape],
   	
    a, a0, Deltalb,
    a1, a2, sz, 
    alist, blist, x0,
    unitcell,
    unitcellnew,
    n, rest, interlayerdistance(* AB-stacking *),
    
    l0, bha, bhb, bh, bpar, v1, v2,
    rsel,
    
    lfun(*layer producing function *)
    },
    
     (* ------------------------------Test of the arguments-----------------------------*)
    If[
     	size < 1,
     	Message[IVGroupQuantumDot::argval,size];
     	Throw[$Failed]
     ];
    
    If[
     	numberoflayers < 1,
     	Message[IVGroupQuantumDot::argval,numberoflayers];
     	Throw[$Failed]
     ];
    
    
   (* parameters were taken from S. Cahangirov et al., PRL, 102(23), 236804 (2009). http://doi.org/10.1103/PhysRevLett.102.236804 *)
   Switch[
    	element,
    	"C", a0 = 1.42; Deltalb = 0;(* \[Angstrom] *)
    	interlayerdistance = 3.35(* \[Angstrom] *),
    	"Si", a0 = 2.21; Deltalb = 0.44;(* \[Angstrom] *)
    	interlayerdistance = 3.35(* \[Angstrom] *),
    	"Ge", a0 = 2.29; Deltalb = 0.64;(* \[Angstrom] *)
    	interlayerdistance = 3.35(* \[Angstrom] *)
    ](* end Switch *);
   
   a = Sqrt[3] a0;
   (* Graphene lattice basic vectors *)
   a1 = a {1/2, Sqrt[3]/2, 0};
   a2 = a {1, 0, 0};
   
   Switch[
    	edgetype,
    	1,(* zigzag type QD *)
    
    	lfun = Flatten[Table[Function[{x}, x + ((i - 1) {0, 0, interlayerdistance} + If[EvenQ[i], (a1 + a2)/3, {0, 0, 0}])] /@ #,{i, numberoflayers}]
       , 1] &;(* end Flatten *)(* layer function *)
    
    	Switch[
    		 shape,
     		"TRI", 
     		sz = size + 1;
     		alist = Flatten[Table[a1 i + a2 j, {i, 0, sz}, {j, 0, sz}], 1];
     		blist = # + (a1 + a2)/3 + {0, 0, -(Deltalb/2)} & /@ alist;
     		alist = # + {0, 0, Deltalb/2} & /@ alist;
     		x0 = sz Sqrt[a2.a2];
     		unitcell = Select[Join[alist, blist], -Sqrt[3] (#[[1]] - x0) >= #[[2]] &];
     		n = Length@unitcell;
     		rest = unitcell[[{1, sz + 1, (sz^2 + 3 sz + 2)/2}]];
    
     		unitcellnew = 
      		lfun@Delete[unitcell, {{1}, {sz+1}, {(sz^2 + 3 sz + 2)/2}}];
      		{unitcellnew, {2 x0, 0, 0}, a0}
      		,
     		"HEX", 
     		sz = size - 1;
     		alist = Flatten[Table[a1 i + a2 j, {i, 0, sz}, {j, 0, sz}], 1];
     		blist = # + (a1 + a2)/3 + {0, 0, -(Deltalb/2)} & /@ alist;
     		alist = # + {0, 0, Deltalb/2} & /@ alist;
     		x0 = sz Sqrt[a2.a2];
     		unitcell = Select[Join[alist, blist], -Sqrt[3] (#[[1]] - x0) >= #[[2]] &];
     		n = Length@unitcell;
     		rest = unitcell[[{1, sz + 1, (sz^2 + 3 sz + 2)/2}]];
    
     		unitcellnew = lfun@Flatten[Table[RotationMatrix[\[Pi]/3 i, {0, 0, 1}].({1, 1, (-1)^(i - 1)} (# + (a1 + a2)/3)) & /@unitcell, {i, 6}], 1];
     		{unitcellnew, {4 x0, 0, 0}, a0}
     	](* end Switch shape *),
    	2,(* armchair type QD *)
    
	    lfun = Flatten[Table[# + ((i - 1) {0, 0, interlayerdistance} + If[EvenQ[i], ({-1, 1, 1} a1 + a2)/3, {0, 0, 0}]) & /@ #,{i, numberoflayers}]
       , 1] &;(* end Flatten *)(* layer function *)
    
	    Switch[
    	 	shape,
     		"TRI",
     		alist =  Flatten[Table[({-1, 1, 1} a1 + a2)/3 + a1 i + a2 j, {i, 0, size}, {j, 0, size}], 1];
     		blist = # + (a1 - a2)/3 + {0, 0, -(Deltalb/2)} & /@ alist;
     		alist = # + {0, 0, Deltalb/2} & /@ alist;
     		x0 = size Sqrt[a2.a2];
     	
     		unitcell = Flatten[Table[RotationMatrix[(2 \[Pi])/3 i, {0, 0, 1}].# & /@ Select[Join[alist, blist], -Sqrt[3]/3 (#[[1]] - x0) >= #[[2]] &], {i, 0, 2}], 
       		1]; (* triangle unitcell *)
     
	   		unitcellnew = lfun@unitcell;
       		{unitcellnew, {2 x0, 0, 0}, a0}
     		,
     		"HEX",
     		l0 = {0, 0, 0};(* origin *)
     		bha = # + {0, 0, Deltalb/2} & /@ {l0, l0 + a1, l0 + {-1, 1, 1} a1}; (* A atoms of basic hexagon *)
     		bhb = # + {0, 0, -(Deltalb/2)} & /@ {l0 + (a1 - a2)/3, l0 + (a1 - a2)/3 + a2, l0 + (a1 - a2)/3 + {-1, 1, 1} a1};(* B atoms of basic hexagon *)
     		bh = Join[bha, bhb];(* basic hexagon *)
     
     		v1 = a1 + {-1, 1, 1} a1;(* translation vector 1 *)
     		v2 = {-1, 1, 1} a1 + a2;(* translation vector 2 *)
     		bpar = Flatten[Table[(v1 (i - 1) + v2 (j - 1) + #) & /@ bh, {i, size}, {j, size}], 2];(* basic parallelogramm *)
     		(* size *)
     		x0 = size Sqrt[v1.v1];
     		rsel = Append[Select[bpar, ((-(Sqrt[3]/3) #[[1]] + x0 > #[[2]]) && #[[1]] > a0 Cos[\[Pi]/6]) &], {-1, 1, 1} a1 + {0, 0, Deltalb/2}];(* rotation symmetry element *)
     		unitcell = Flatten[Table[RotationMatrix[\[Pi]/3 i, {0, 0, 1}].({1, 1, (-1)^(i - 1)} (# - {0, a0, 0})) & /@ rsel, {i, 6}], 1];
     		unitcellnew = lfun@unitcell;
     		{unitcellnew, {4 x0, 0, 0}, a0}
     	](* end Switch shape *)
    
   ](* end Switch type *)
 ](* end Module *)
](* end Catch *);
SyntaxInformation[IVGroupQuantumDot] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

(* Quantum dots group V added 28/07/2020 *)
(* code from 01/10/2017 *)
Mixed = 3;

Options[VGroupQuantumDot] = {
	ChemicalElement -> "P",
	EdgeType -> Zigzag,
	Shape -> "TRI"
};

VGroupQuantumDot[size_Integer, numberoflayers_Integer, OptionsPattern[]] := Module[
   {
   	element = OptionValue[ChemicalElement],
   	edgetype = OptionValue[EdgeType],
   	shape = OptionValue[Shape],
   	
    a0, a, Deltalb, sft, Phi,
    Size,
    tol = 0.01(* tolerance *),
    a1, a2,
    alist, blist,
    x0, unitcell,
    phosphorenesheet,
    polygon,
    generation = 4,
    polygonfractalinfo,
    unitcellnew,
    interlayerdistance(* ABA-stacking *),
    shift,
    
    lfun(*layer producing function *),
    InPolygonQ,
    InPolygonFractalQ,
    PolygonFractal
    },
   (*---------------------------------------- subfunctions ----------------------------------------*)

  (*  plist_: points must be ordered in clockwise direction; the first point must not be duplicated *)
   InPolygonQ[point_, plist_, \[Delta]_] :=
    And @@ MapThread[If[
        #2[[1]] != 0,
        If[
         	#2[[1]] > 0,
         	#2[[2]]/#2[[1]] point[[
             1]] + (#1[[2]] - #2[[2]]/#2[[1]] #1[[1]]) - \[Delta] >= 
          point[[2]],
         	#2[[2]]/#2[[1]] point[[
             1]] + (#1[[2]] - #2[[2]]/#2[[1]] #1[[1]]) + \[Delta] <= 
          point[[2]]
         ](* end If *),
        If[
         	#2[[2]] > 0,(*#2\[Equal]0 && #3\[Equal]0 - skiped *)
         	#1[[1]] + \[Delta] <= point[[1]],
         	#1[[1]] - \[Delta] >= point[[1]]
         ](* end If *)
        ] &(* end If *), {plist, RotateLeft[plist, 1] - plist}];
   
   InPolygonFractalQ[point_, plist_, polfracinfo_] :=
     Module[
     {
      overlappingtriangles,
      lasttriangleinfo,
      \[Delta] = 0.01,
      
      InTriangleQ
      },
     InTriangleQ = 
      If[#[[2]] > 
         0, {InPolygonQ[point, #[[1]], \[Delta]], #[[2]], #[[
          3]]}, {InPolygonQ[point, Reverse[#[[1]]], \[Delta]], #[[
          2]], #[[3]]}] &;(* designed to be applied to PolygonFractal third argument (ilist) *) 
     
     (* selects those triangles (or initial polygon) that contains the given point, i.e. True elements in tinfo, 
     and chooses the one with the highes generation *)
     overlappingtriangles = Select[Join[{{InPolygonQ[point, plist, \[Delta]], 1, 0}}, Flatten[Map[InTriangleQ, polfracinfo, {2}], 1]], #[[1]] &];
     If[
      	overlappingtriangles == {},
      	False,
      	lasttriangleinfo = Sort[overlappingtriangles, #1[[3]] > #2[[3]] &][[1]];
      	(* if the highest generation triangle has the second parameter in the list equal to -1 then the result is converted to False,
      	 which means the point should be considered as being outside the triangle *)
      	And[lasttriangleinfo[[1]], Switch[lasttriangleinfo[[2]], -1, False, 1, True]]
      ](* end If *) 
     ](* end Module *);
     
     
     (* 
     	plist_: points must be ordered in clockwise direction; the first point must not be duplicated;
     	fractaltype_: 1 - ideal fractal; 2 - positive random; 3 - negative random; 4 - fully random;
     *)
     PolygonFractal[plist_, fractaltype_, gen_] := Module[
 	 {
   		list,
   		tlist,
   		ilist,
   		plen,
   
  		 (* subfunction *)
  		 EdgeFractal
   	},
  
  	EdgeFractal[point1_, point2_, type_, gn_] := Module[
    	{
     		vec,
     		p1, p2, p3, fac,
     		Phifac
     	},
    	vec = (point2 - point1)/3;
    	fac = Switch[type, 1, {1, 1}, 2, RandomReal[1, 2], 3, RandomReal[1, 2], 4, RandomReal[1, 2]];
    	Phifac = Switch[type, 1, 1, 2, 1, 3, -1, 4, (-1)^RandomInteger[{-1, 1}]^2];
    	p1 = point1 + fac[[1]] vec;
    	p2 = p1 + vec;
    	p3 = p1 + vec/2 + Sqrt[3]/2 fac[[2]] Phifac {-vec[[2]], vec[[1]]};
    
    	AppendTo[tlist, {p1, p3, p2}]; (* triangle list *)
    	AppendTo[ilist, {{p1, p3, p2}, Phifac, gn}]; (* information list *)
    
    	If[
    	 	gn == gen,
    	 	AppendTo[list, {point1, p1, p3, p2, point2}]
    	 ];
    
    	If[
     		gn < gen,
     		EdgeFractal[point1, p1, type, gn + 1];
     		EdgeFractal[p1, p3, type, gn + 1];
     		EdgeFractal[p3, p2, type, gn + 1];
     		EdgeFractal[p2, point2, type, gn + 1]
     	];
    
    ];(* end Module EdgeFractal *)
  
  	plen = Length[plist];
  	Transpose[Table[list = {}; tlist = {}; ilist = {}; EdgeFractal[plist[[i]], plist[[If[i == plen, 1, i + 1]]], fractaltype, 1]; {list, tlist, ilist}, {i, plen}]]
	](* end Module *);

   
   (*------------------------------------------------------------------------------------------------*)

   
   
   
   (* parameters were taken from 
   [1] M. Ezawa, New J. Phys. 16, 115004 (2014).  http://doi.org/10.1088/1367-2630/16/11/115004
   see also 
   [2] A. Castellanos-Gomez et al., 2D Mater. 1, 025001 (2014). http://doi.org/10.1088/2053-1583/1/2/025001
   *)
   Switch[
    	element,
    	"P",
    	a0 = 2.164; (* \[Angstrom] *)
    	a = 2.537;   (* \[Angstrom] *)
    	Deltalb = 2.1443;(* \[Angstrom]; vertical shift *)
    	sft = 0.5223;(* \[Angstrom]; horizontal shift *)
    	Phi = 0.7; (* radians; angle between a1 vector and Ox-axis *)
    	interlayerdistance = 5.445;(* \[Angstrom]; interlayer distance; taken from [2] *)
    	shift = {1.67857, 0, 0}; (* shift to combine the origin of reference system with one of the atoms of the lattice (in xOy plane) *)
    ](* end Switch *);
   
   
   (* Phosphorene lattice basic vectors *)
   a1 = a {Cos[Phi], Sin[Phi], 0};
   a2 = a {Cos[Phi], -Sin[Phi], 0};
   
   lfun = Flatten[
      Table[
       Function[{x}, x + ((i - 1) {0, 0, interlayerdistance} + If[EvenQ[i], {0, a Sin[Phi], 0}, {0, 0, 0}])] /@ #,
       {i, numberoflayers}]
      , 1] &;(* end Flatten *)(* layer function *)
   
   Switch[
    		edgetype,
    		1,(* zigzag type QD *)
    
    		Switch[
     				shape,
     				"TRI",
     				Size = size + 1;
     				x0 = (Size a2)[[1]];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, 0, Size}, {j, 0, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				
     				unitcell = Delete[Select[Join[alist, blist], x0 >= #[[1]] &], {{1}, {Size + 1}, {(Size^2 + 3 Size + 2)/2}}];
     				unitcellnew = lfun@unitcell,
     				"TRI(r)",
     				Size = Ceiling[4/3 size];
     				x0 = (Size a2)[[1]];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     				
     				(* cut initial phosphorene sheet to smaller one to facilitate fractal calculation *)
     				polygon = Most /@ (Size {-((a1 + a2)/2) , (3 a1 - a2)/4, a1 + a2, (3 a2 - a1)/4 });
     				phosphorenesheet = lfun@Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     				
     				(* bounding fractal generation *)
     				polygon = Most /@ (size {{0, 0, 0}, a1, a2});
     				polygonfractalinfo = PolygonFractal[polygon, 4, generation];
     				
     				(* selection of those atoms of the sheet which are within the bounding fractal *)
     				unitcell = Select[phosphorenesheet, InPolygonFractalQ[{#[[1]], #[[2]]}, polygon, polygonfractalinfo[[3]]] &];
     
     				(* number of layers does work with fractals at the moment! *)
     				unitcellnew = unitcell
     				,
     				"HEX",
     				Size = size;
     				x0 = Size Sqrt[a2.a2];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     
     				polygon = Most /@ (size {a1, a2, a2 - a1, -a1, -a2, a1 - a2});
     
     				unitcell = Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     				unitcellnew = lfun@unitcell
     				,
     				"HEX(r)",
     				Size = Ceiling[4.3/3 size];
     				x0 = Size Sqrt[a2.a2];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     				
     				(* cut initial phosphorene sheet to smaller one to facilitate fractal calculation *)
     				polygon = Most /@ (Size {a1, a2, a2 - a1, -a1, -a2, a1 - a2});
     				phosphorenesheet = lfun@Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     				
     				(* bounding fractal generation *)
     				polygon = Most /@ (size {a1, a2, a2 - a1, -a1, -a2, a1 - a2});
     				polygonfractalinfo = PolygonFractal[polygon, 4, generation];
     				
     				(* selection of those atoms of the sheet which are within the bounding fractal *)
     				unitcell = Select[phosphorenesheet, InPolygonFractalQ[{#[[1]], #[[2]]}, polygon, polygonfractalinfo[[3]]] &];
     				(* number of layers does work with fractals at the moment! *)
     				unitcellnew = unitcell
     		](* end Switch *);
    
    		Switch[
     				shape,
     				"TRI(r)",
     				{unitcellnew, {4 x0, 0, 0}, a0},
     				"HEX(r)",
     				{unitcellnew, {4 x0, 0, 0}, a0},
     				"TRI", 
     				{unitcellnew, {2 x0, 0, 0}, a0},
     				"HEX", 
     				{unitcellnew, {4 x0, 0, 0}, a0}
     		](* end Switch *),
    		2,(* armchair type QD *)
    
    		Switch[
     				shape,
     				"TRI",
     				Size = size + 1;
     				x0 = Size Sqrt[a2.a2];
     
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     
     				polygon = Most /@ ((size + 1/2) {a2, -a1, a1 - a2});
     				unitcell = Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     
     				unitcellnew = lfun@unitcell
     				,
     				"TRI(r)",
     				Size = Ceiling[4/3 size];
     				x0 = Size Sqrt[a2.a2];
     
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     				
     				(* cut initial phosphorene sheet to smaller one to facilitate fractal calculation *)
     				polygon = Most /@ (Size {a1, a2, a2 - a1, -a1, -a2, a1 - a2});
     				phosphorenesheet = lfun@Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     				
     				(* bounding fractal generation *)
     				polygon = Most /@ ((size + 1/2) {a2, -a1, a1 - a2});
     				polygonfractalinfo = PolygonFractal[polygon, 4, generation];
     				
     				(* selection of those atoms of the sheet which are within the bounding fractal *)
     				unitcell = Select[phosphorenesheet, InPolygonFractalQ[{#[[1]], #[[2]]}, polygon, polygonfractalinfo[[3]]] &];
     				(* number of layers does work with fractals at the moment! *)
     				unitcellnew = unitcell
     				,
     				"HEX",
     				Size = 2 size;
     				x0 = Size Sqrt[a2.a2];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     
     				polygon = Most /@ ((size - 1/2) {2 a1 - a2, (a1 + a2), 2 a2 - a1, a2 - 2 a1, -(a1 + a2), a1 - 2 a2});
     				unitcell = Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     
     				unitcellnew = lfun@unitcell
     				,
     				"HEX(r)",
     				Size = Ceiling[(2 + 2/3) size];
     				x0 = Size Sqrt[a2.a2];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     
     				(* cut initial phosphorene sheet to smaller one to facilitate fractal calculation *)
     				polygon = Most /@ (Size {a1, a2, a2 - a1, -a1, -a2, a1 - a2});
     				phosphorenesheet = lfun@Select[phosphorenesheet, InPolygonQ[#, polygon, tol] &];
     
     				(* bounding fractal generation *)
     				polygon = Most /@ ((size - 1/2) {2 a1 - a2, (a1 + a2), 2 a2 - a1, a2 - 2 a1, -(a1 + a2), a1 - 2 a2});
     				polygonfractalinfo = PolygonFractal[polygon, 4, generation];
     				
     				(* selection of those atoms of the sheet which are within the bounding fractal *)
     				unitcell = Select[phosphorenesheet, InPolygonFractalQ[{#[[1]], #[[2]]}, polygon, polygonfractalinfo[[3]]] &];
     				(* number of layers does work with fractals at the moment! *)
     				unitcellnew = unitcell
     
     		](* end Switch shape *);
    
    		Switch[
     				shape,
     				"TRI(r)",
     				{unitcellnew, {4 x0, 0, 0}, a0},
     				"HEX(r)",
     				{unitcellnew, {4 x0, 0, 0}, a0},
     				"TRI", 
     				{unitcellnew, {2 x0, 0, 0}, a0},
     				"HEX", 
     				{unitcellnew, {4 x0, 0, 0}, a0}
     		](* end Switch *),
    		3,(* Mixed edge type QQ *)
    		Switch[
     				shape,
     				"CIR",
     				Size = size + 1;
     				x0 = Size Sqrt[a2.a2];
     
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     
     				unitcell = Select[phosphorenesheet, #[[1]]^2 + #[[2]]^2 <= (size a)^2 &];
     				unitcellnew = lfun@unitcell
     				,
     				"CIR(r)",
     				Size = Ceiling[7/6 size];
     				x0 = Size Sqrt[a2.a2];
     				alist = Flatten[Table[a1 i + a2 j + {0, 0, (-1)^(i + j + 1) Deltalb/2}, {i, -Size, Size}, {j, -Size, Size}], 1];
     				blist = ({1, 1, -1} #) + {sft, 0, 0} & /@ alist;
     				phosphorenesheet = shift + # & /@Join[alist, blist];
     
     				(* cut initial phosphorene sheet to smaller one to facilitate fractal calculation *)
     				phosphorenesheet = Select[phosphorenesheet, #[[1]]^2 + #[[2]]^2 <= (Size a)^2 &];
     
     				(* bounding fractal generation, ring is approximated by 9 edges polygon *)
     				polygon = size a Reverse@CirclePoints[9];
     				polygonfractalinfo = PolygonFractal[polygon, 4, generation];
     
     				(* selection of those atoms of the sheet which are within the bounding fractal *)
     				unitcell = Select[phosphorenesheet, InPolygonFractalQ[{#[[1]], #[[2]]}, polygon, polygonfractalinfo[[3]]] &];
     				
     				(* number of layers does not work with fractals at the moment! *)
     				unitcellnew = unitcell
      		];(* end Switch shape *)
    
    		Switch[
     				shape,
     				"CIR",
     				{unitcellnew, {4 x0, 0, 0}, a0},
     				"CIR(r)",
     				{unitcellnew, {4 x0, 0, 0}, a0}
     		](* end Switch *)
    
    ](* end Switch type *)
   
   ](* end Module *);
SyntaxInformation[VGroupQuantumDot] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};


End[] (* End Private Context *)

EndPackage[]