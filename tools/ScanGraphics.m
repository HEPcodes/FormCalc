(*
	ScanGraphics.m
		special graphing routines for parameter scans
		last modified 20 Apr 01 th
*)


BeginPackage["ScanGraphics`"]

ScanPlot3D::usage = "ScanPlot3D[var1, var2, n, opts] makes a 3D plot of a
scan performed in the variables var1 and var2. The data to be plotted must
be stored in data sets (e.g. as by ReadData) and n specifies which data
sets are used (5 means 1...5, or {7, 10} means 7...10). The argument opts
may contain additional Graphics3D options."

ScanDensityPlot::usage = "ScanDensityPlot[var1, var2, n, opts] makes a
density plot of a scan performed in the variables var1 and var2. The data
to be plotted must be stored in data sets (e.g. as by ReadData) and n
specifies which data sets are used (5 means 1...5, or {7, 10} means
7...10). The argument opts may contain additional Graphics options."

MagnifyingFactor::usage = "MagnifyingFactor is an option of ScanPlot3D.
Specifying MagnifyingFactor -> m rescales all z-coordinates by m. This is
used to enhance very `flat' scans."

RenderMissing::usage = "RenderMissing is an option of ScanDensityPlot.
It determines the graphics primitive used for rendering missing values."

$MissingPoints::usage = "$MissingPoints contains a list of x-y-coordinates
which have been found to be missing in the scan. $MissingPoints is
updated by ScanPlot3D and ScanDensityPlot."

Assign::usage = "Assign[z, {var1, var2}, n] is called by ScanPlot3D and
ScanDensityPlot for every data set that is plotted to build a grid of plot
points. Assign is supposed to define a value for z[v1, v2], where v1 and
v2 are the numerical values of var1 and var2 as given by the nth parameter
set. By default, the one-loop result at the lowest energy in a data set is
taken, i.e. z[v1, v2] = Data[n][[1, 3]]. Assign can of course be redefined
to perform more complex tasks, like finding the minimum of all sets that
have the same values of var1 and var2."


Begin["`Private`"]

Complete[l_List] :=
Block[ {dif, c, delta},
  dif = Apply[#2 - #1 &, Partition[l, 2, 1], 1];
  c = Count[dif, #]&/@ dif;
  delta = dif[[ Position[c, Max[c], 1, 1][[1, 1]] ]];
  Range[Min[l], Max[l], delta]
]


Assign[z_, {var1_, var2_}, n_] :=
Block[ {v1, v2},
  {v1, v2} = {var1, var2} /. Global`Para[n];
  z[v1, v2] = Global`Data[n][[1, 3]]
]


ScanPlot3D::missing =
"Warning: there are missing points in your data set. See $MissingPoints."

Options[ScanPlot3D] = {MagnifyingFactor -> 1000.}

ScanPlot3D[var1_, var2_, n_, graphopts___Rule] :=
Block[ {mag, x, y, z, xyz, polys, graph, dummy},
  mag = MagnifyingFactor /. {graphopts} /. Options[ScanPlot3D];
  Do[Assign[z, {var1, var2}, x], Evaluate[Flatten[{x, n}]]];
  x = Complete[ Union[#[[1, 1, 1]]&/@ DownValues[z]] ];
  y = Complete[ Union[#[[1, 1, 2]]&/@ DownValues[z]] ];
  z[_, _] = dummy;
  xyz = Outer[{#1, #2, mag z[#1, #2]}&, x, y];

  $MissingPoints =
    Drop[#, -1]&/@ Select[Flatten[xyz, 1], !FreeQ[#, dummy]&];
  If[Length[$MissingPoints] =!= 0, Message[ScanPlot3D::missing]];

  polys = Table[
    Polygon[{xyz[[i, j]],
             xyz[[i + 1, j]],
             xyz[[i + 1, j + 1]],
             xyz[[i, j + 1]]}],
    {i, Length[x] - 1}, {j, Length[y] - 1} ]//Flatten;

  graph = Graphics3D[ Select[polys, FreeQ[#, dummy]&],
    Axes -> True,
    AxesLabel -> {ToString[var1], ToString[var2], "sigma"},
    PlotRange -> All,
    Sequence@@ DeleteCases[{graphopts}, MagnifyingFactor -> _]
  ];
  Show[graph];
  graph
]


ScanDensityPlot::missing = ScanPlot3D::missing

Options[ScanDensityPlot] = {RenderMissing -> RGBColor[0, 0, 1]}

ScanDensityPlot[var1_, var2_, n_, graphopts___Rule] :=
Block[ {miss, x, y, z, xyz, min, max, rects, gs, coorgs, graph, dummy},
  miss = RenderMissing /. {graphopts} /. Options[ScanDensityPlot];
  Do[Assign[z, {var1, var2}, x], Evaluate[Flatten[{x, n}]]];
  x = Complete[ Union[#[[1, 1, 1]]&/@ DownValues[z]] ];
  y = Complete[ Union[#[[1, 1, 2]]&/@ DownValues[z]] ];
  z[_, _] = dummy;
  xyz = Flatten[Outer[{#1, #2, z[#1, #2]}&, x, y], 1];

  $MissingPoints =
    Drop[#, -1]&/@ Select[xyz, !FreeQ[#, dummy]&];
  If[Length[$MissingPoints] =!= 0, Message[ScanDensityPlot::missing]];

  xyz = Cases[Last/@ xyz, _Real];
  min = Min[xyz];
  max = Max[xyz];

  rects = Flatten[Outer[
    ( gs = Flatten[Outer[z, ##] /. coord -> List];
      gs = (Plus@@ gs/Length[gs] - min)/(max - min);
      {GrayLevel[gs], Rectangle@@ Transpose[{##} /. coord -> List]} )&,
    Apply[coord, Partition[x, 2, 1], 1],
    Apply[coord, Partition[y, 2, 1], 1]
  ]] /. GrayLevel[x_] :> miss /; !FreeQ[x, dummy];

  graph = Graphics[rects,
    PlotRange -> All,
    Frame -> True,
    FrameLabel -> ToString/@ {var1, var2},
    Sequence@@ DeleteCases[{graphopts}, RenderMissing -> _]
  ];
  Show[graph];
  graph
]

End[]

EndPackage[]

