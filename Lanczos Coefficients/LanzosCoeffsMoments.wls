#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*"To compute the Lanczos coefficients we need a package called Combinatrica which contains the function " *)


(* ::Input:: *)
(*Needs["Combinatorica`"]*)


RetAmp4=ToExpression[Import["C:\\Users\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Return Amplitudes\\ReturnAmpNumL4TFIM.txt"]];


Plot[RetAmp4[[1]],{t,0,10}, PlotRange->{0,Automatic}]


Clear[\[Phi]n, NNND, NNNDTab, \[Phi]t, PPlot, NNNn];


Needs["Combinatorica`"]


nnn = 37;


NNNn[t_] := Simplify[Series[RetAmp4[[3]], {t, 0, nnn} ]]


\[Mu]Tab = Table[Factorial[j]SeriesCoefficient[ NNNn[t], j], {j, 0, nnn}];


MijR = Simplify[   Sum[Table[KroneckerDelta[i+j+1-mm]*\[Mu]Tab[[mm]], {i, 0, (nnn)/2}, {j, 0, (nnn)/2}], {mm, 1, Length[\[Mu]Tab]}]   ]; 
bProd = Table[Det[MijR[[1;;k, 1;;k]]      ]/Det[MijR[[1;;(k-1), 1;;(k-1)]]      ], {k, 2, nnn/2+1}];   


(* ::Input:: *)
(*bTab = Table[If[j==1, Sqrt[- bProd[[j]]   ],  Sqrt[- bProd[[j]]/bProd[[j-1]]    ]   ], {j, 1, Length[bProd]}];    aTab =Table[Simplify[-I(    Simplify[   Cofactor[MijR[[1;;k, 1;;k]], {k-1,k}]/Det[MijR[[1;;k-1, 1;;k-1]]]   ] -    If[k > 2, Simplify[   Cofactor[MijR[[1;;k-1, 1;;k-1]], {k-2,k-1}]/Det[MijR[[1;;k-2, 1;;k-2]]]   ]  , 0])   ]  , {k, 2, Length[bTab]}   ];  *)


aTab


bTab


ListLinePlot[bTab]


(* ::Section:: *)
(*L=5*)


RetAmp5=ToExpression[Import["C:\\Users\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Return Amplitudes\\ReturnAmpNumL5TFIM.txt"]];


Plot[RetAmp5[[3]],{t,0,10},PlotRange->All]


Clear[\[Phi]n, NNND, NNNDTab, \[Phi]t, PPlot, NNNn];


Needs["Combinatorica`"]


nnn = 33;


NNNn[t_] := Simplify[Series[RetAmp5[[3]], {t, 0, nnn} ]]


\[Mu]Tab = Table[Factorial[j]SeriesCoefficient[ NNNn[t], j], {j, 0, nnn}];


MijR = Simplify[   Sum[Table[KroneckerDelta[i+j+1-mm]*\[Mu]Tab[[mm]], {i, 0, (nnn)/2}, {j, 0, (nnn)/2}], {mm, 1, Length[\[Mu]Tab]}]   ]; 
bProd = Table[Det[MijR[[1;;k, 1;;k]]      ]/Det[MijR[[1;;(k-1), 1;;(k-1)]]      ], {k, 2, nnn/2+1}];   


(* ::Input:: *)
(*bTab = Table[If[j==1, Sqrt[- bProd[[j]]   ],  Sqrt[- bProd[[j]]/bProd[[j-1]]    ]   ], {j, 1, Length[bProd]}];    aTab =Table[Simplify[-I(    Simplify[   Cofactor[MijR[[1;;k, 1;;k]], {k-1,k}]/Det[MijR[[1;;k-1, 1;;k-1]]]   ] -    If[k > 2, Simplify[   Cofactor[MijR[[1;;k-1, 1;;k-1]], {k-2,k-1}]/Det[MijR[[1;;k-2, 1;;k-2]]]   ]  , 0])   ]  , {k, 2, Length[bTab]}   ];  *)


bTab
