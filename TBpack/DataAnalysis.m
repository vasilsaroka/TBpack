(* Wolfram Language Package *)

BeginPackage["TBpack`DataAnalysis`", {"TBpack`CustomTicks`", "TBpack`MaTeX`"}]
(* Exported symbols added here with SymbolName::usage *)  

(* Functions *)
ListOfBonds::usage = "ListOfBonds[\*StyleBox[\"unitcell\",\"TI\"], \*StyleBox[\"bondlength\",\"TI\"]] returns point pairs from \*StyleBox[\"unitcell\",\"TI\"] such that distance between the points is \*StyleBox[\"bondlength\",\"TI\"] \[PlusMinus] 0.05 \[Angstrom].
ListOfBonds[\*StyleBox[\"unitcell\",\"TI\"], \*StyleBox[\"bondlength\",\"TI\"], BondLengthDelta \[Rule] \*StyleBox[\"value\",\"TI\"]] returns point pairs from \*StyleBox[\"unitcell\",\"TI\"] such that distance between the points is \*StyleBox[\"bondlength\",\"TI\"] \[PlusMinus] \*StyleBox[\"value\",\"TI\"] \[Angstrom].
ListOfBonds[\*StyleBox[\"unitcell\",\"TI\"], \*StyleBox[\"bondlength\",\"TI\"], Dimentionality \[Rule] 2] returns a list of pairs \*RowBox[{\"{\", \"{\", StyleBox[SubscriptBox[\"x\",\"1\"],\"TI\"], \",\" , StyleBox[SubscriptBox[\"y\",\"1\"],\"TI\"] , \"}\" , \",\" , \"{\", StyleBox[SubscriptBox[\"x\",\"2\"],\"TI\"], \",\", StyleBox[SubscriptBox[\"y\",\"2\"],\"TI\"] ,\"}\" ,\"}\"}].";
AtomicStructure::usage = "AtomicStructure[{\*StyleBox[\"unitcell\",\"TI\"], \*StyleBox[\"tr\",\"TI\"], \*StyleBox[SubscriptBox[\"a\",\"0\"],\"TI\"]}] returns the ball-and-stick model of a \*StyleBox[\"system\",\"TI\"] specified by unit cell \*StyleBox[\"unitcell\",\"TI\"], translation vector \*StyleBox[\"tr\",\"TI\"] and the lattice constant \*StyleBox[SubscriptBox[\"a\",\"0\"],\"TI\"].
AtomicStructure[\*StyleBox[\"system\",\"TI\"], \*StyleBox[\"options\",\"TI\"]] uses option settings specified in \*StyleBox[\"options\",\"TI\"].
AtomicStructure[\*StyleBox[\"system\",\"TI\"], NumberOfUnitCells \[Rule] \*StyleBox[\"value\",\"TI\"]] returns the number of unit cells set by \*StyleBox[\"value\",\"TI\"]. 
AtomicStructure[\*StyleBox[\"system\",\"TI\"], Dimensionality \[Rule] 2] returns the \*StyleBox[\"system\",\"TI\"] projection to the \*StyleBox[\"xy\",\"TI\"]-plane.
AtomicStructure[\*StyleBox[\"system\",\"TI\"], AtomEnumeration \[Rule] \*StyleBox[\"True\",\"TI\"]] enumerates atom sites in the \*StyleBox[\"system\",\"TI\"].
AtomicStructure[\*StyleBox[\"system\",\"TI\"], FontSize \[Rule] \*StyleBox[\"value\",\"TI\"]] sets the font size for the atom enumeration labels to \*StyleBox[\"value\",\"TI\"].
AtomicStructure[\*StyleBox[\"system\",\"TI\"], ImageSize \[Rule] \*StyleBox[\"value\",\"TI\"]] sets the image size of the produced graphics to \*StyleBox[\"value\",\"TI\"].
AtomicStructure[\*StyleBox[\"system\",\"TI\"], AtomColor \[Rule] \*StyleBox[\"color\",\"TI\"]] sets the atom color to \*StyleBox[\"color\",\"TI\"] in displayed graphics.
AtomicStructure[\*StyleBox[\"system\",\"TI\"], BondColor \[Rule] \*StyleBox[\"color\",\"TI\"]] sets the bond color to \*StyleBox[\"color\",\"TI\"] in displayed graphics.
AtomicStructure[\*StyleBox[\"system\",\"TI\"], TextColor \[Rule] \*StyleBox[\"color\",\"TI\"]] sets the label text color to \*StyleBox[\"color\",\"TI\"] in displayed graphics.";
ReciprocalVectors::usage = "ReciprocalVectors[\*StyleBox[\"list\",\"TI\"]] returns a list of reciprical vectors for the given \*StyleBox[\"list\",\"TI\"] of the primitive translations.
ReciprocalVectors[{\*StyleBox[\"tr\",\"TI\"] }] for a 1D lattice with the primitive translation \*StyleBox[\"tr\",\"TI\"].
ReciprocalVectors[{\*StyleBox[SubscriptBox[\"tr\",\"1\"],\"TI\"], \*StyleBox[SubscriptBox[\"tr\",\"2\"],\"TI\"] }] for a 2D lattice with the primitive translations \*StyleBox[SubscriptBox[\"tr\",\"1\"],\"TI\"] and \*StyleBox[SubscriptBox[\"tr\",\"2\"],\"TI\"].
ReciprocalVectors[{\*StyleBox[SubscriptBox[\"tr\",\"1\"],\"TI\"], \*StyleBox[SubscriptBox[\"tr\",\"2\"],\"TI\"], \*StyleBox[SubscriptBox[\"tr\",\"3\"],\"TI\"] }] for a 3D lattice with the primitive translations \*StyleBox[SubscriptBox[\"tr\",\"1\"],\"TI\"], \*StyleBox[SubscriptBox[\"tr\",\"2\"],\"TI\"] and \*StyleBox[SubscriptBox[\"tr\",\"3\"],\"TI\"].";
FermiEnergy::usage = "FermiEnergy[\*StyleBox[\"bands\",\"TI\"]] estimates the intrinsic Fermi level for the \*StyleBox[\"bands\",\"TI\"] that are returned by ElectronicBands1D.";
EnergyGap::usage = "EnergyGap[\*StyleBox[\"bands\",\"TI\"]] estimates the energy band gap for the \*StyleBox[\"bands\",\"TI\"] that are returned by ElectronicBands1D.";
PlotElectronicBands1D::usage = "PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"]\)] plots the energy bands \*StyleBox[\"bands\",\"TI\"] shifting them by \*StyleBox[\"fermilevel\",\"TI\"] and normalizing \*StyleBox[\"k\",\"TI\"]-points by \*StyleBox[\"knorm\",\"TI\"] factor.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], \*StyleBox[\"options\",\"TI\"]\)] uses the specified options.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], Labels \[Rule] {\*StyleBox[\"title\",\"TI\"], {\*StyleBox[\"xlabel\",\"TI\"], \*StyleBox[\"ylabel\",\"TI\"]}, \*StyleBox[\"legend\",\"TI\"]}\)] uses the specified title, axes and legend labels.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], PlotAppearance \[Rule] \*StyleBox[\"integer\",\"TI\"]\)] uses the plot appearance wrapper specified by \*StyleBox[\"integer\",\"TI\"].
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], TitlePosition \[Rule] {\*StyleBox[\"x\",\"TI\"], \*StyleBox[\"y\",\"TI\"]}\)] positions the title inset within the plot at Scaled[{\*StyleBox[\"x\",\"TI\"], \*StyleBox[\"y\",\"TI\"]}].
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], TitleFrame \[Rule] \*StyleBox[\"True\",\"TI\"]\)] puts a black frame around the plot title.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], TitleBackground \[Rule] \*StyleBox[\"color\",\"TI\"]\)] sets the title background color to \*StyleBox[\"color\",\"TI\"].
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], NumberDisplay \[Rule] {\*StyleBox[SubscriptBox[\"integer\",\"1\"],\"TI\"], \*StyleBox[SubscriptBox[\"integer\",\"2\"],\"TI\"]}\)] diplays numbers in tick labels with \*StyleBox[SubscriptBox[\"integer\",\"1\"],\"TI\"] total digits and \*StyleBox[SubscriptBox[\"integer\",\"2\"],\"TI\"] digits after decimal point.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], YStep \[Rule] \*StyleBox[\"real\",\"TI\"]\)] sets major ticks along \*StyleBox[\"y\",\"TI\"]-axis at \*StyleBox[\"real\",\"TI\"] interval.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], XStep \[Rule] \*StyleBox[\"real\",\"TI\"]\)] sets major ticks along \*StyleBox[\"x\",\"TI\"]-axis at \*StyleBox[\"real\",\"TI\"] interval.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], YSubdivisions \[Rule] \*StyleBox[\"integer\",\"TI\"]\)] sets \*StyleBox[\"integer\",\"TI\"] number of minor ticks at \*StyleBox[\"y\",\"TI\"]-axis.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], XSubdivisions \[Rule] \*StyleBox[\"integer\",\"TI\"]\)] sets \*StyleBox[\"integer\",\"TI\"] number of minor ticks at \*StyleBox[\"x\",\"TI\"]-axis.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], YTickLen \[Rule] \*StyleBox[\"len\",\"TI\"]\)] sets the real number \*StyleBox[\"len\",\"TI\"] to be the major tick length at \*StyleBox[\"y\",\"TI\"]-axis.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], XTickLen \[Rule] \*StyleBox[\"len\",\"TI\"]\)] sets the real number \*StyleBox[\"len\",\"TI\"] to be the major tick length at \*StyleBox[\"x\",\"TI\"]-axis.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], YTickDirection \[Rule] \*StyleBox[\"value\",\"TI\"]\)] uses the inward, outward or double side ticks at \*StyleBox[\"y\",\"TI\"]-axis for \*StyleBox[\"Inward\",\"TI\"], \*StyleBox[\"Outward\",\"TI\"] and \*StyleBox[\"Both\",\"TI\"] values of the \*StyleBox[\"value\",\"TI\"].
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], XTickDirection \[Rule] \*StyleBox[\"value\",\"TI\"]\)] uses the inward, outward or double side ticks at \*StyleBox[\"x\",\"TI\"]-axis for \*StyleBox[\"Inward\",\"TI\"], \*StyleBox[\"Outward\",\"TI\"] and \*StyleBox[\"Both\",\"TI\"] values of the \*StyleBox[\"value\",\"TI\"].
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], TicksPattern \[Rule] {{\*StyleBox[\"left\",\"TI\"], \*StyleBox[\"right\",\"TI\"]}, {\*StyleBox[\"bottom\",\"TI\"], \*StyleBox[\"top\",\"TI\"]}}\)] uses Boolean values \*StyleBox[\"left\",\"TI\"], \*StyleBox[\"right\",\"TI\"], \*StyleBox[\"bottom\",\"TI\"] and \*StyleBox[\"top\",\"TI\"] to set up ticks at the plot frame sides.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], PlotRange \[Rule] {{\*StyleBox[\"xmin\",\"TI\"], \*StyleBox[\"xmax\",\"TI\"]}, {\*StyleBox[\"ymin\",\"TI\"], \*StyleBox[\"ymax\",\"TI\"]}}\)] uses the specified plot range.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], PlotStyle \[Rule] \*StyleBox[\"style\",\"TI\"]\)] uses the specified \*StyleBox[\"style\",\"TI\"] to plot \*StyleBox[\"bands\",\"TI\"], can use TBpackPallete for Automatic.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], FrameTicksStyle \[Rule] \*StyleBox[\"style\",\"TI\"]\)] uses the specified \*StyleBox[\"style\",\"TI\"] to format plot ticks.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], FrameStyle \[Rule] \*StyleBox[\"style\",\"TI\"]\)] uses the specified \*StyleBox[\"style\",\"TI\"] to format plot frame.
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], AspectRatio \[Rule] \*StyleBox[\"value\",\"TI\"]\)] sets plot aspect ratio to \*StyleBox[\"value\",\"TI\"].
PlotElectronicBands1D[\!\(\*StyleBox[\"bands\",\"TI\"], \*StyleBox[\"fermilevel\",\"TI\"], \*StyleBox[\"knorm\",\"TI\"], ImageSize \[Rule] \*StyleBox[\"size\",\"TI\"]\)] uses the specified \*StyleBox[\"size\",\"TI\"] to format ticks.";

PlotGrid::usage = "PlotGrid[\!\({{\*StyleBox[SubscriptBox[\"pl\",\"1\"],\"TI\"], \*StyleBox[SubscriptBox[\"pl\",\"2\"],\"TI\"], \[Ellipsis] }}\)] tiles the list of plots \*StyleBox[SubscriptBox[\"pl\",\"i\"],\"TI\"] into one image.
PlotGrid[\!\({{\*StyleBox[SubscriptBox[\"pl\",\"11\"],\"TI\"], \*StyleBox[SubscriptBox[\"pl\",\"12\"],\"TI\"], \[Ellipsis] },{\*StyleBox[SubscriptBox[\"pl\",\"21\"],\"TI\"], \*StyleBox[SubscriptBox[\"pl\",\"22\"],\"TI\"], \[Ellipsis] }}\)] tiles the list of plots \*StyleBox[SubscriptBox[\"pl\",\"ij\"],\"TI\"] into one image.
PlotGrid[\!\(\*StyleBox[\"plist\",\"TI\"], \*StyleBox[\"options\",\"TI\"] \)] uses the specified \*StyleBox[\"options\",\"TI\"].
PlotGrid[\!\(\*StyleBox[\"plist\",\"TI\"], ImageSize \[Rule] {\*StyleBox[\"width\",\"TI\"], \*StyleBox[\"height\",\"TI\"] }\)] produces image with the specified \*StyleBox[\"width\",\"TI\"] and \*StyleBox[\"height\",\"TI\"].
PlotGrid[\!\(\*StyleBox[\"plist\",\"TI\"], ImagePadding \[Rule] {{\*StyleBox[\"left\",\"TI\"], \*StyleBox[\"right\",\"TI\"]}, {\*StyleBox[\"bottom\",\"TI\"], \*StyleBox[\"top\",\"TI\"] }\)] produces image with the specified \*StyleBox[\"left\",\"TI\"], \*StyleBox[\"right\",\"TI\"], \*StyleBox[\"bottom\",\"TI\"] and \*StyleBox[\"top\",\"TI\"] image paddings around the tiled plots.";




(* Options *)
NumberOfUnitCells::usage = "Option setting the number of unit cell to be displayed.";
Dimensionality::usage = "Option setting the number of dimensions used in displaying an object.";
AtomEnumeration::usage = "Option taking values \*StyleBox[\"True\",\"TI\"] or \*StyleBox[\"False\",\"TI\"] specifying if the atom numbers should be displayed."
BondLengthDelta::usage = "Option setting in \[Angstrom] the maximum deviation of the bond length in the unit cell from the explicitly specified bond length value.";
AtomColor::usage = "Option setting the color of the atoms displayed in graphical objects.";
BondColor::usage = "Option setting the color of the bonds displayed in graphical objects."
TextColor::usage = "Option setting the color of the text displayed in graphical objects."
Labels::usage = "Option for plotting functions setting the labels such as plot title, axes label, frame label and legends.";
PlotAppearance::usage = "Option for plotting functions setting the way they are presented.";
TitlePosition::usage = "Option for plotting funcitons setting the position of the title within the plot.";
TitleFrame::usage = "Option for plotting funcitons specifying if the the plot title is framed.";
TitleBackground::usage = "Option for plotting funcitons setting the background style of the title label.";
NumberDisplay::usage = "Option setting the number of digits \*StyleBox[SubscriptBox[\"int\",\"1\"],\"TI\"] and the number of digits after decimal point \*StyleBox[SubscriptBox[\"int\",\"2\"],\"TI\"] in a string form of the number: {\*StyleBox[SubscriptBox[\"int\",\"1\"],\"TI\"], \*StyleBox[SubscriptBox[\"int\",\"2\"],\"TI\"]}.";
YStep::usage = "Option in some plotting functions setting the step for major ticks for y-axis.";
YSubdivisions::usage = "Option in some plotting functions setting the number of minor ticks for y-axis.";
XStep::usage = "Option in some plotting functions setting the step for major ticks for x-axis";
XSubdivisions::usage = "Option in some plotting functions setting the number of minor ticks for x-axis.";
YTickLen::usage = "Option in some plotting functions setting the length of ticks for y-axis scaled to the width of the plot.";
XTickLen::usage = "Option in some plotting functions setting the length of ticks for x-axis scaled to the height of the plot.";
YTickDirection::usage = "Option in some plotting functions setting \*StyleBox[\"Inward\",\"TI\"], \*StyleBox[\"Outward\",\"TI\"] or \*StyleBox[\"Both\",\"TI\"]  direction of ticks for y-axis."
XTickDirection::usage = "Option in some plotting functions setting \*StyleBox[\"Inward\",\"TI\"], \*StyleBox[\"Outward\",\"TI\"] or \*StyleBox[\"Both\",\"TI\"] direction of ticks for x-axis."
TicksPattern::usage = "Option is some plotting functions setting a pattern of the frame ticks displaying: {{\*StyleBox[\"left\",\"TI\"], \*StyleBox[\"right\",\"TI\"]}, {\*StyleBox[\"bottom\",\"TI\"], \*StyleBox[\"top\",\"TI\"]}}, where \*StyleBox[\"left\",\"TI\"], \*StyleBox[\"right\",\"TI\"], \*StyleBox[\"bottom\",\"TI\"] and \*StyleBox[\"top\",\"TI\"] are \*StyleBox[\"True\",\"TI\"] or \*StyleBox[\"False\",\"TI\"]."

(* Constants *)
TBpackPalette::usage = "The list of standard colors used by default in plotting functions.";


Begin["`Private`"] (* Begin Private Context *)

Options[ListOfBonds] = {
		Dimensionality -> 3,
		BondLengthDelta -> 0.05
		};
ListOfBonds[unitcell_List, bondlength_, OptionsPattern[]] := Module[
{
	dim = OptionValue[Dimensionality],
	delta = OptionValue[BondLengthDelta],
	len, vec,
	bondlist = {}
},

len = Length@unitcell;

Do[
	vec = unitcell[[i]] - unitcell[[j]];
	If[
		Abs[Sqrt[vec.vec] - bondlength] < delta,
		AppendTo[bondlist,Switch[
									dim,
									2,
									Most/@{unitcell[[i]],unitcell[[j]]},
									3,
									{unitcell[[i]],unitcell[[j]]}
							](* end Switch *)
		](* end AppendTo *)
	](* end If *),
{i,1,len-1},
{j,i+1,len}
](* end Do *);
Return[bondlist]
](* end Module *);
SyntaxInformation[ListOfBonds] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

Options[AtomicStructure] = {
	BondLengthDelta -> 0.05,
	NumberOfUnitCells -> 1,
	Dimensionality -> 3,
	AtomEnumeration -> False,
	AtomColor -> Gray,
	BondColor -> Gray,
	TextColor -> White,
	FontSize -> 10,
	ImageSize -> {300}
};

(* 2D structure visualization is not included *)
AtomicStructure[system_List, OptionsPattern[]]:=Module[
{
	bondlengthdelta = OptionValue[BondLengthDelta],
	lim = OptionValue[NumberOfUnitCells],
	dim = OptionValue[Dimensionality],
	acolor = OptionValue[AtomColor],
	bcolor = OptionValue[BondColor],
	tcolor = OptionValue[TextColor],
	imagesize = OptionValue[ImageSize],
	fontsize = OptionValue[FontSize],
	unitcell, T, a0,
	u,v,
	
	counter = 0,
	antiacolor,
	antibcolor,
	
	nd,
	AntiColor
},

AntiColor[color_] := Module[
  {
   list, gr,
   r, g, b
   },
  list = List @@ color;
  Switch[
   		Length@list,
   		1, {gr} = list; GrayLevel[1 - gr],
   		3, {r, g, b} = list; RGBColor @@ {1 - r, 1 - g, 1 - b}
   ](* end Switch *)](* end Module *);

nd = ToString@PaddedForm[1.0 #,{5,3},ExponentFunction->(If[-16<#<16,Null,#]&)]&; (* number display function *);

antiacolor = AntiColor[acolor];
antibcolor = AntiColor[bcolor];

(* system = {unitcell,translationvector,lattice constant} *)
unitcell = system[[1]];
T = system[[2]];
a0 = system[[3]];

(* TODO: generalize to 2 and 3 dimensions 12/08/2020 *)
u = Flatten[Table[(# + T (i-1))&/@unitcell,{i,lim}],1];

Switch[
		dim,
		3,
		Graphics3D[{
			Specularity[GrayLevel[1],100],
			Tooltip[
				counter++;
				Mouseover[
							{
								acolor,Sphere[#,0.3],
								If[
									OptionValue[AtomEnumeration],
									{tcolor,Text[Style[counter,fontsize],#]},
									{}
								](* end If AtomicEnumeration *)
							},
							{
								antiacolor,Sphere[#,0.3]
							}
						](* end Mouseover *)
				,"("<>nd[#[[1]]]<>","<>nd[#[[2]]]<>","<>nd[#[[3]]]<>") \[Angstrom]",
				TooltipStyle -> {Background -> White}
			](* end Tooltip *)&/@u,
			Tooltip[
				Mouseover[
					{bcolor,Tube[#,0.1]},
					{antibcolor,Tube[#,0.1]}
					](* end Mouseover *),
							v = #[[2]] - #[[1]];
							nd[Sqrt[v.v]]<>" \[Angstrom]",
							TooltipStyle -> {Background -> White }
			](* end Tooltip *)&/@ListOfBonds[u, a0, Dimensionality->3, BondLengthDelta->bondlengthdelta]
			},
			Boxed->False,
			Lighting-> "Neutral",
			ImageSize->imagesize],
		2,
		Graphics[{
			
			Tooltip[
				Mouseover[
					{bcolor,AbsoluteThickness[1],Line[#]},
					{antibcolor,AbsoluteThickness[3],Line[#],bcolor,AbsoluteThickness[1],Line[#]}
					](* end Mouseover *),
				v = #[[2]] - #[[1]];
				nd[Sqrt[v.v]]<>" \[Angstrom]",
				TooltipStyle -> {Background -> White}
			](*end Tooltip *)&/@ListOfBonds[u, a0, Dimensionality->2, BondLengthDelta->bondlengthdelta],
			Tooltip[
				counter++;
				Mouseover[
							{
								acolor,Disk[Most@#,0.3],
								If[
									OptionValue[AtomEnumeration],
									{tcolor,Text[Style[ToString[counter],fontsize],Most@#]},
									{}
								](* end If AtomicEnumeration *)
							},
							{
								EdgeForm[acolor],antiacolor,Disk[Most@#,0.3]
							}
						](* end Mouseover *),
				"("<>nd[#[[1]]]<>","<>nd[#[[2]]]<>") \[Angstrom]",
				TooltipStyle->{Background->White}
			](* end Tooltip *)&/@u
			},ImageSize->imagesize]
](* end Switch *)
](* end Module *);
SyntaxInformation[AtomicStructure] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};


(* Error messages *)
ReciprocalVectors::falseargumentstructure = 
     "The argument of the \!\(\*
StyleBox[\"ReciprocalVectors\",\nFontSlant->\"Italic\"]\) must be a \
list of 3-component vectors, numeric or symbolic.";
ReciprocalVectors::badargument = "The argument of the \!\(\*
StyleBox[\"ReciprocalVectors\",\nFontSlant->\"Italic\"]\) contains or \
was reduced to dependent vectors.";
ReciprocalVectors::somevaluesdroped = 
     "One or more last unzero vector components were droped during \
evaluation.";

ReciprocalVectors[listofvectors_] := Catch[
   Module[
    {
     list, len,
     k, b, eq,
     testfunction,
     f
     },
    
    testfunction = (MatchQ[#, 
         List[arg__?(VectorQ[#] && Length[#] == 3 &)]] && 
        Length[#] <= 3) &;
    If[
     	VectorQ[listofvectors] && Length[listofvectors] == 3 ,
     	list = {listofvectors},
     	list = listofvectors;
     ];
    
    If[
     	! testfunction[list],
     	Message[ReciprocalVectors::falseargumentstructure];
     	Throw[$Failed]
     ];
    
    len = Length[list];
    
    f = If[
       		! MatchQ[Take[#, len - 3], List[arg__?(# == 0 &)]],
       		Message[ReciprocalVectors::somevaluesdroped];
       		Drop[#, len - 3],
       		Drop[#, len - 3]
       	](* end If *)&;
    list = Map[f, list];
    
    If[
     	Det[list] == 0,
     	Message[ReciprocalVectors::badargument];
     	Throw[$Failed]
     ];
    
    k = Table[Symbol["k" <> ToString[i]], {i, len}];
    b = Permutations[Drop[{2 \[Pi], 0, 0}, len - 3]];
    Table[
     eq = 
      Table[With[{x = i}, Hold[k.list[[x]] == b[[j, x]]]], {i, 1, 
         len}] // ReleaseHold;
     Join[k /. Flatten@Solve[eq, k], ConstantArray[0, 3 - len]](* 
     end Join *)
     , {j, 1, len}]
    
  ](* end Module *)
](* end Catch *);




(* Error messages *)
FermiEnergy::argument = "The \!\(\*StyleBox[\"Most\",\"TI\"]\) part of the argument of \!\(\*StyleBox[\"FermiEnergy\",\nFontSlant->\"Italic\"]\) must be a numeric matrix.";

FermiEnergy[bands_List] := Catch[
   Module[
    {
     dim, x
     },
   
    dim = Dimensions[bands];
    
    (* ------------------------------Test of the arguments-----------------------------*)
    If[
     	(!MatrixQ[Most@bands, NumericQ]),
     	Message[FermiEnergy::argument];
     	Throw[$Failed]
     ];
    
    
    
    If[ 
     	(dim[[1]] - 1) // EvenQ,
     	(* an isolator or semiconductor because electrons fill all the bands (states) in pairs with spins up and down *)
     	x = (dim[[1]] - 1)/2;
     	(Max[bands[[x]]] + Min[bands[[x + 1]]])/2
     	,
     	(* metal because one band (state) should be half filled; the Fermi level is defined as the energy center of the band *)
     	x = bands[[Ceiling[(dim[[1]] - 1)/2]]];
     	(Max[x] + Min[x])/2
     ](* end If *)
  ](* end Module *)
](* end Catch *);



(* Error messages *)
EnergyGap::argument = "The \!\(\*StyleBox[\"Most\",\"TI\"]\) part of the argument of \!\(\*StyleBox[\"EnergyGap\",\nFontSlant->\"Italic\"]\) must be a numeric matrix.";

EnergyGap[bands_List] := Catch[
   Module[
    {
     dim, x
     },
    
    dim = Dimensions[bands];
    
    (* ------------------------------Test of the arguments-----------------------------*)
    If[
     	(!MatrixQ[Most@bands, NumericQ]),
     	Message[EnergyGap::argument];
     	Throw[$Failed]
     ];
    
    If[ 
     	(dim[[1]] - 1) // EvenQ,
     	(* an isolator or semiconductor because electrons fill all the bands (states) in pairs with spins up and down *)
     	x = (dim[[1]] - 1)/2;
     	Min[bands[[x + 1]]] - Max[bands[[x]]]
     	,
     	(* metal because one band (state) should be half filled *)
     	0
     ](* end If *)
  ](* end Module *)
](* end Catch *);



(* standard colors *)
TBpackPalette = {
	RGBColor[1, 0, 0],
	RGBColor[0, 1, 1], 
	RGBColor[0.8738063654339829, 0.5624850137259609, 0.03917956167261294], 
   	RGBColor[0.12619363456601707`, 0.43751498627403906`, 0.9608204383273871], 
   	RGBColor[0.8439079900024566, 0.3445418478490483, 0.8709555478667004], 
   	RGBColor[0.15609200999754336`, 0.6554581521509517, 0.12904445213329963`], 
   	RGBColor[0.022999318091469156`, 0.4213913494203896, 0.49057419000523383`], 
   	RGBColor[0.9770006819085308, 0.5786086505796104, 0.5094258099947662], 
   	RGBColor[0.21151508704203614`, 0.997254589585209, 0.8216256808736591], 
   	RGBColor[0.7884849129579639, 0.0027454104147910385`, 0.17837431912634094`]
};


(* plotting  electronic bands *)
(* Error messages *)
PlotElectronicBands1D::arguments = "The three arguments and PlotStyle option value must be same length lists.";

PlotElectronicBands1D::farg = "The first argument must be be a numeric matrix or a vector such matrices.";

PlotElectronicBands1D::notanumber = "The value `1` is not a number."

PlotElectronicBands1D::psoptnocolors = "The default PlotStyle pallete does not contain enough colors to present the given bands."

PlotElectronicBands1D::psoptnotstyle = "The PlotStyle option value is not RGBColor or GrayLevel color or Directive."

PlotElectronicBands1D::legendlbls = "The third part of Labels option does not contain enough labels to make a legend.";

PlotElectronicBands1D::plotrange = "The plot range cannot be set because either y tick length `1` or x tick length `2` devited by aspect ratio `3` is greater or equal 0.5.";



Options[PlotElectronicBands1D] = {
	Labels -> {"", {"k", "E-E_{\\mathrm{F}},\\text{ eV}"}, Automatic}(*{"", {"\*SubscriptBox[\"k\",\"\"]", "\*SubscriptBox[\"E\",\"\"]- \*SubscriptBox[\"E\",\"F\"], eV"}, Automatic}*),
	PlotAppearance -> 1,
	TitlePosition -> {0.5, 0.93},
	TitleFrame -> True,
	TitleBackground -> White,
	NumberDisplay -> {2,1},
	YStep -> 5,
	XStep -> 0.5,
	YSubdivisions -> 5,
	XSubdivisions -> 5,
	YTickLen -> 0.05,
	XTickLen -> 0.05,
	YTickDirection -> Outward,
	XTickDirection -> Outward,
	TicksPattern -> {{True, False}, {True, False}},
	PlotRange -> Automatic,
	PlotStyle -> Red,
	FrameTicksStyle -> (*Directive[AbsoluteThickness[0.5],Black]*)BlackFrame,
    FrameStyle -> (*Directive[AbsoluteThickness[0.5],Black]*)BlackFrame,
    AspectRatio -> 3,
	ImageSize -> {Automatic,380}
};
PlotElectronicBands1D[bands_List, fermilevel_, knormfactor_, OptionsPattern[]] := Catch[
   Module[
    {(*--------------------------- Variables and constants ---------------------------*)

	 labels = OptionValue[Labels],
	 paflag = OptionValue[PlotAppearance],
     plotrange = OptionValue[PlotRange],
     plotstyle = OptionValue[PlotStyle],
     titlepos = OptionValue[TitlePosition],
     aspectratio = OptionValue[AspectRatio],
     
     blen,
     
     linestyles = {Dashing[1], Dashing[0.025], Dashing[{0.02, 0.006, 0.006, 0.006}]},
     
     dy = OptionValue[YStep],
     dx = OptionValue[XStep],
     yticklen = OptionValue[YTickLen],
     xticklen = OptionValue[XTickLen],
     ytickdir = OptionValue[YTickDirection],
     xtickdir = OptionValue[XTickDirection],
     tickspattern = OptionValue[TicksPattern],
     ny = OptionValue[YSubdivisions],
     nx = OptionValue[XSubdivisions],
     
     bandslist,
     fermilevellist,
     knormfactorlist,
     
     xmin, xmax,
     ymin, ymax,
     yticks, xticks,
     fticks,
     
     
     showtittle,
     title, inset,
     axeslabels,
     
     plotrangepadding,
     len,p,
     
     
     bandsQ, plotstyleQ,
     ts, ls, nd,
     data2plotdata, lbfun, nolabel
     },
     
     (* test functions *)
     bandsQ = (MatrixQ[#, NumericQ]) &;
     plotstyleQ = MatchQ[#, RGBColor[_, _, _]] || MatchQ[#, GrayLevel[_]] || MatchQ[#, Directive[__]] &;
     
     (* ------------------------------Test of the arguments-----------------------------*)
	  
     If[
     		bandsQ[bands],
    		bandslist = {bands};
    		blen = 1,
    		If[
    			VectorQ[bands,bandsQ],
    			bandslist = bands;
    			blen = Length[bands],
    			Message[PlotElectronicBands1D::farg];
    			Throw[$Failed]
    		]
      ];
    	
      If[
     		Length[fermilevel] == 0,
     		If[
     			NumberQ[fermilevel],
     			fermilevellist = {fermilevel},
     			Message[PlotElectronicBands1D::notanumber,fermilevel];
     			Throw[$Failed]
     		],
     		fermilevellist = fermilevel
       ];
     	
     	If[
     		Length[knormfactor] == 0,
     		If[
     			NumberQ[knormfactor],
     			knormfactorlist = {knormfactor},
     			Message[PlotElectronicBands1D::notanumber,knormfactor];
     			Throw[$Failed]
     		],
     		knormfactorlist = knormfactor
     	];
     	
       
    	If[	
    		plotstyle === Automatic,
    		If[
    			blen <= 10,
    			plotstyle = Table[Directive[linestyles[[Mod[i,3,1]]],TBpackPalette[[i]]],{i,blen}],
    			Message[PlotElectronicBands1D::psoptnocolors];
    			Throw[$Failed]
    		]
    		,
       		If[
    			plotstyleQ[plotstyle],
    			plotstyle = {plotstyle},
     			If[
     				!VectorQ[plotstyle,plotstyleQ],
     				Message[PlotElectronicBands1D::psoptnotstyle];
     				Throw[$Failed]
     			]
    		]
    	  ];
     	
     	
		If[
    		Not[blen==Length[fermilevellist]==Length[knormfactorlist]==Length[plotstyle]],
    		Message[PlotElectronicBands1D::arguments];
    		Throw[$Failed]
    	];
    	    	
    	If[
     		labels[[3]] === Automatic,
     		labels[[3]] = Range[blen],
     		If[
     			Not[Length[labels[[3]]] == blen],
     			Message[PlotElectronicBands1D::legendlbls];
    			Throw[$Failed]
     		]
    	];
    	
    	If[
    		labels[[1]]===""||labels[[1]]===Null||labels[[1]]===None,
    		showtittle = False,
    		showtittle = True
    	];
    	
    	If[
    		plotrange === Automatic,
    	 	plotrange = {
    	 		{Min[MapThread[#1 Last[#2]&,{knormfactorlist,bandslist}]], Max[MapThread[#1 Last[#2]&,{knormfactorlist,bandslist}]]},
    	 		{Floor@Min[MapThread[First[#2]-#1 &,{fermilevellist,bandslist}]], Ceiling@Max[MapThread[#2[[-2]]-#1&,{fermilevellist,bandslist}]]}
    	 		}    		
    	];
    	

     (* ------------------------------ end Test of the arguments-----------------------------*)
    

    (* text styles *)
    ts[st_] := MaTeX(*Style*)[st, FontSize -> 21];(* plot title style *)
    ls[st_] := MaTeX(*Style*)[st, FontSize -> 20]; (* plot label style *)
    
    
    (* number display function *)
    nd = ls@ToString@PaddedForm[1.0 #,OptionValue[NumberDisplay],ExponentFunction->(If[-16<#<16,Null,#]&)]&;
    
    (* transforming bands to the ListLinePlot input format *)
    data2plotdata[d_, eff_, nfactor_] := Module[
      {
       ed, kd, dim
       },
      ed = Map[(# - eff) &, Most@d, {2}];(* shift the Fermi Energy to zero *)
      kd = Last@d;
      dim = Dimensions[ed];
      Table[
       					{nfactor kd[[j]], ed[[i, j]]},
       					{i, 1, dim[[1]]},
       					{j, 1, dim[[2]]}
       			](* end Table *)
      ](* end Module *);
      
    (* legend function *)
    lbfun = Labeled[#, LineLegend[plotstyle, ls /@ labels[[3]]], Right] &;
    
    (* removes labels from the ticks specification list *)
    nolabel = {#[[1]], "", #[[3]], #[[4]]} &;
    

    
    
    (* paflag -- plot appearance flag. 
    It returns reduced verstion of the plot. 1- full version, 2- 
    begin version, 3- middle version, 4- end version *)
    (* labels = {title, {Xaxeslabel, Yaxeslabel}, legendlabels} *)
    
    (* plot title as plot inset *)
    title = If[showtittle,ts@labels[[1]]];
    inset = If[
    			showtittle,
    			{Inset[Grid[{{title}}, Frame -> OptionValue[TitleFrame], FrameStyle -> Directive[AbsoluteThickness[0.5],Black]], 
    				Scaled[titlepos], 
    				Background -> OptionValue[TitleBackground]]},
    			{}
    ];
    
    
    axeslabels = Switch[
      				paflag,
      				1, {ls@labels[[2, 1]], ls@labels[[2, 2]]},
      				2, {ls@labels[[2, 1]], ls@labels[[2, 2]]},
      				3, {ls@labels[[2, 1]], None},
      				4, {ls@labels[[2, 1]], None}
      			](* end Switch *);
   
   
   {{xmin, xmax}, {ymin, ymax}} = plotrange;
    yticks = 
     LinTicks[ymin, ymax, dy, ny, 
      MajorTickLength -> yticklen, 
      MinorTickLength -> yticklen/2,
      TickDirection -> ytickdir,
      TickLabelFunction -> nd
      ];
      
    xticks = 
     LinTicks[xmin, xmax, dx, nx, 
      MajorTickLength -> xticklen, 
      MinorTickLength -> xticklen/2,
      TickDirection -> xtickdir,
      TickLabelFunction -> nd
      ];
      
    fticks = {Switch[
         				paflag,
         				1, {
         					If[tickspattern[[1,1]],yticks, None],
         					If[tickspattern[[1,2]],nolabel/@yticks, None]
         					},
         				2, {
         					If[tickspattern[[1,1]],yticks, None],
         					If[tickspattern[[1,2]],nolabel/@yticks, None]
         					},
         				3, {
         					If[tickspattern[[1,1]],nolabel /@ yticks, None], 
         					If[tickspattern[[1,2]],nolabel /@ yticks, None]
         					},
         				4, {
         					If[tickspattern[[1,1]],nolabel /@ yticks, None],
         					If[tickspattern[[1,2]],nolabel /@ yticks, None]
         					}
         		](* end Switch *),
    			{
    				If[tickspattern[[2,1]],xticks, None],
    				If[tickspattern[[2,2]],nolabel/@xticks, None]
    			}
    };
    
    
    (* making sure the ticks never overlap with the plotted data *)
    If[
    	(yticklen < 0.5) && (xticklen/aspectratio < 0.5),
    	plotrangepadding = {Scaled[yticklen],Scaled[xticklen/aspectratio]},
    	Message[PlotElectronicBands1D::plotrange, yticklen, xticklen, aspectratio];
    	Throw[$Failed]
    ];
   
    len = Length /@ bandslist;
    
    p = ListLinePlot[
      				Flatten[MapThread[data2plotdata,{bandslist, fermilevellist, knormfactorlist}],1],
      				PlotStyle -> Flatten[Table[plotstyle[[i]], {i, blen}, {j, (len[[i]] - 1)}], 1],
      				PlotRange -> plotrange,
      				PlotRangePadding -> plotrangepadding,
      				Axes -> False,
      				Frame -> True,
      				FrameStyle -> OptionValue[FrameStyle],
      				FrameTicksStyle -> OptionValue[FrameTicksStyle],
      				FrameLabel -> axeslabels,
      				FrameTicks -> fticks,
      				Epilog ->inset,
      				AspectRatio -> aspectratio,
      				ImageSize -> OptionValue[ImageSize]
      		](* end ListLinePlot *);
      		
    (* wrapping the plot with a legend *)  		
    Switch[paflag, 1, lbfun@p, 2, p, 3, p, 4, lbfun@p]
  ](* end Module *)
](* end Catch *);
SyntaxInformation[PlotElectronicBands1D] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

   
   
   
(* this function is modified code from stackexchange: 
https://mathematica.stackexchange.com/questions/6877/do-i-have-to-\code-each-case-of-this-grid-full-of-plots-separately/6882#6882
*)
  
(* 05/05/2019: this function uses explicitly provided values for the imagepadding parameters. 
It does not distort subplots for non-equal left and right or top and bottom image padding parameters.
*)

(* error messages *)
PlotGrid::argnotglist = "The argument is not a list of Graphics or Labeled graphics.";

Options[PlotGrid] = {
	ImageSize -> {300, 400},
	ImagePadding -> {{80,80},{80,80}}	
};
PlotGrid[l_?MatrixQ, OptionsPattern[]] := Catch[Module[
  {
   w = OptionValue[ImageSize][[1]],
   h = OptionValue[ImageSize][[2]],
   leftPadding = OptionValue[ImagePadding][[1,1]],
   rightPadding = OptionValue[ImagePadding][[1,2]], 
   bottomPadding = OptionValue[ImagePadding][[2,1]],
   topPadding = OptionValue[ImagePadding][[2,2]],
   
   g,fticks,aspectratio,
   plotrange,fac,
   
   nx, ny,
   widths, heights,
   positions,
   
   graphicsQ
},
   
  (* test functions *)
 (* graphicsQ = If[
   					Head[#] === Labeled, 
   					If[
   						Head[#[[1]]] === Graphics, 
   						#[[1]],
   						Message[PlotGrid::argnotglist];
   						Throw[$Failed]
   					],
   					If[
   						Head[#] === Graphics,
   						#,
   						Message[PlotGrid::argnotglist];
   						Throw[$Failed]
   					]
   			] &;
  (* ------------------------------Test of the arguments-----------------------------*)
  l = Map[graphicsQ, l, {2}];*)
     
  {ny, nx} = Dimensions[l];
  widths = (w - leftPadding - rightPadding)/nx Table[1, {nx}];
  widths[[1]] = widths[[1]] + leftPadding;
  widths[[-1]] = widths[[-1]] + rightPadding;
  heights = (h - topPadding - bottomPadding)/ny Table[1, {ny}];
  heights[[1]] = heights[[1]] + bottomPadding;
  heights[[-1]] = heights[[-1]] + topPadding;
  
  
  
  positions = Transpose@Partition[Tuples[Prepend[Accumulate[Most[#]], 0] & /@ {widths, heights}], ny];
  
  
  Graphics[
  			Table[
  					Inset[
     						Show[
     								(* adjust xticks of the plot *)
     								(* this does not work *)
     								g = l[[ny - j + 1, i]];
     								Options[g,AspectRatio][[1,2]];
     								fticks = Options[g, FrameTicks][[1, 2]];
     								plotrange = Options[g, PlotRange][[1, 2]];
     								aspectratio = Options[g, AspectRatio][[1, 2]];
     								fac = (plotrange[[2,2]] - plotrange[[2,1]])/ (2 aspectratio);
     								fticks[[2, 1]] = {#[[1]], #[[2]], fac #[[3]], #[[4]]} & /@ fticks[[2, 1]];
     								g /. (FrameTicks->_) -> (FrameTicks -> fticks), 
      								ImagePadding -> {
      											{If[i == 1, leftPadding, 0], If[i == nx, rightPadding, 0]},
      											 {If[j == 1, bottomPadding, 0], If[j == ny, topPadding, 0]}
      										},
      						 AspectRatio -> Full(*Full*)
      						](* end Show *), 
     						positions[[j, i]],
     						{Left, Bottom}, 
     						{widths[[i]], heights[[j]]}
     						](* end Inset *),
     			{i, 1, nx}, {j, 1, ny}](* end Table *), 
   				PlotRange -> {{0, w}, {0, h}}, 
   				ImageSize -> {w, h}
   	](* end Graphics *)
](* end Module *)]; (* end Catch *)
SyntaxInformation[PlotGrid] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

End[] (* End Private Context *)

EndPackage[]