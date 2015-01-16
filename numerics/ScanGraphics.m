(*
	ScanGraphics.m
		special graphing routines for parameter scans
		last modified 1 Jun 00 th
*)


BeginPackage["ScanGraphics`"]

ScanPlot3D::usage = "ScanPlot3D[var1, var2, n, opts] makes a 3D plot of a
scan performed in the variables var1 and var2. It assumes that the data to
be plotted have been loaded with ReadData and are contained in the data
sets Data[1]...Data[n]. (n is the output of ReadData.) The argument opts
may contain additional Graphics3D options."

ScanDensityPlot::usage = "ScanDensityPlot[var1, var2, n, opts] makes a
density plot of a scan performed in the variables var1 and var2. It
assumes that the data to be plotted have been loaded with ReadData and are
contained in the data sets Data[1]...Data[n]. (n is the output of
ReadData.) The argument opts may contain additional Graphics options."

MagnifyingFactor::usage = "MagnifyingFactor is an option of ScanPlot3D.
Specifying MagnifyingFactor -> m rescales all z-coordinates by m. This is
used to enhance very `flat' scans."

RenderMissing::usage = "RenderMissing is an option of ScanDensityPlot.
It determines the graphics primitive used for rendering missing values."

$MissingPoints::usage = "$MissingPoints contains a list of x-y-coordinates
which have been found to be missing in the scan. $MissingPoints is
updated by ScanPlot3D and ScanDensityPlot."


Begin["`Private`"]

Complete[l_List] :=
Block[ {dif, c, delta},
  dif = Apply[#2 - #1 &, Partition[l, 2, 1], 1];
  c = Count[dif, #]&/@ dif;
  delta = dif[[ Position[c, Max[c], 1, 1][[1, 1]] ]];
  Range[Min[l], Max[l], delta]
]


ScanPlot3D::missing =
"Note: there are missing points in your data set. See $MissingPoints."

Options[ScanPlot3D] = {MagnifyingFactor -> 1000.}

ScanPlot3D[var1_, var2_, n_, graphopts___Rule] :=
Block[ {mag, x, y, z, xyz, polys, graph, dummy},
  mag = MagnifyingFactor /. {graphopts} /. Options[ScanPlot3D];
  Array[ Set@@ {z[var1, var2] /. Global`Para[#],
                mag Global`Data[#][[1, 3]]} &, n ];
  z[_, _] = dummy;
  x = Complete[ Array[var1 /. Global`Para[#] &, n]//Union ];
  y = Complete[ Array[var2 /. Global`Para[#] &, n]//Union ];
  xyz = Outer[{#1, #2, z[#1, #2]}&, x, y];

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
  Array[ Set@@ {z[var1, var2] /. Global`Para[#],
                Global`Data[#][[1, 3]]}&, n];
  z[_, _] = dummy;
  x = Complete[ Array[var1 /. Global`Para[#] &, n]//Union ];
  y = Complete[ Array[var2 /. Global`Para[#] &, n]//Union ];
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

