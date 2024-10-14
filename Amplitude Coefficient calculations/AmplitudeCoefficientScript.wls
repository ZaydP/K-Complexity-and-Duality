#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]]


{as,bs}=ToExpression[Import["C:\\Users\\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Lanczos Coefficients\\bs_Lanczos_PC_L5TFIMAlt_periodic_5site_RanNormHerm12345.txt"]][[3;;4]];
bsLength=Length[bs]
(*bs=Chop[ToExpression[bs],10^-1000];
as=Chop[ToExpression[as],10^-1000];*)
bsPrec=Precision[bs[[12]]]
as=Table[as[[i]]/.x_/;Precision@x==MachinePrecision->N[Rationalize@as[[i]],bsPrec],{i,1,Length[as]}];
bs=Table[bs[[i]]/.x_/;Precision@x==MachinePrecision->N[Rationalize@bs[[i]],bsPrec],{i,1,Length[bs]}];


ListLinePlot[{as,bs}]


RetAmp=ToExpression[Import["C:\\Users\\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Return Amplitudes\\ReturnAmpJacoL5TFIMAlt_periodic_5site_RanNormHerm12345.txt"]];
RetAmpPrecision=Precision[RetAmp]


Plot[Re@RetAmp,{t,0,10}]


tMax=100;
tPts=500;
Prec=200;
tVals = Table[N[tMax*(i)/(tPts),Prec], {i, 1, tPts}];


InfoPA="ProbAmp_L5TFIMAlt_periodic_5site_RanNormHerm12345";
(*InfoK="K_L6_cdMode1";*)
Info2=StringJoin["bs length, precision: ", ToString[bsLength], ", ", ToString[bsPrec], ".  Return Amplitude precision: ", ToString[RetAmpPrecision]];
Info3=StringJoin["tMax: ",ToString[tMax], "  tPts: ",ToString[tPts],"  Prec: ", ToString[Prec]];


AmplitudeCoefficientScriptHPC[as_,Bs_,ReturnAmp_,tVals_,Prec_]:=Module[{bs,\[Phi]Pts,NNND,NNNDTab,Genf,FtoNSwap,\[Phi]t,cc,PTab},
bs=Bs[[1+LengthWhile[Bs,#==0&];;-1-LengthWhile[Reverse@Bs,#==0&]]]; (*removes tail of zeros and initial zeros*)
 \[Phi]Pts = Length[bs];
(*"Here the derivatives are computed at the time-points.  "*)
	 NNND[0] = ReturnAmp;
	Do[NNND[j+1] =  D[NNND[j], t]    , {j, 0, \[Phi]Pts}  ];//EchoTiming;
	NNNDTab = ParallelTable[NNND[j]/.t->tVals, {j, 0, \[Phi]Pts}, Method->"CoarsestGrained"];//EchoTiming;
(*"Here we compute the coefficeints cc "*)
Clear[cc];   cc[0,0]=1; cc[0,-1] = 0;
(*CHANGED TO HAVE N ONLY GO TO \[Phi]Pts-1*)
Do[  cc[m,n+1]  = N[(If[m>= 1, I*cc[m-1,n], 0]  -If[m<= n, as[[n+1]]*cc[m,n], 0]   - If[m<=n -1&& n>0, bs[[n]]*cc[m,n-1] , 0 ])/bs[[n+1]], Prec]  , {n, 0, \[Phi]Pts-1}, {m, 0, n+1}]//EchoTiming; 
PTab = Table[  Sum[cc[m,kk]*NNNDTab[[m+1]],{m,0,kk}] , {kk, 0, \[Phi]Pts}]//EchoTiming]


EchoTiming[AC=AmplitudeCoefficientScriptHPC[as,bs,RetAmp,tVals,Prec];,"Probability Amplitude time"]


Export[StringJoin[InfoPA,".m.gz"],ToString[{InfoPA,Info2,Info3,AC},InputForm]]


x=ToExpression[Import["ProbAmp_L7_periodic_sxMode1.txt"]][[4]]


Export[StringJoin[InfoPA,".m.gz"],{InfoPA,Info2,Info3,AC}]


EchoTiming[K=Table[{tVals[[j]],Sum[(k-1)*Abs[AC[[k,j]]^2],{k,1,21}]},{j,1,Length[tVals]}];,"Complexity sum time: "]
EchoTiming[PSum=Table[{tVals[[j]],Sum[Abs[AC[[k,j]]^2],{k,1,21}]},{j,1,Length[tVals]}];,"Probability sum time: "]
Export[StringJoin[InfoK,".m.gz"],{InfoK,Info2,Info3,K,PSum}]


Import["K_L7_cdMode3.m.gz"][[1;;3]]
Import["K_L7_cdMode4.m.gz"][[1;;3]]


cut=152
P=Table[{tVals[[j]],Total[Abs[AC[[1;;cut,j]]]^2]},{j,1,tPts}];
K=Table[{tVals[[j]],Range[cut] . Abs[AC[[1;;cut,j]]]^2},{j,1,tPts}];
ListLinePlot[P, PlotRange->All]
ListLinePlot[K,PlotRange->All]


Range[cut]


 


AC[[100,10]]


ListLinePlot[Import["K_L7_cdMode4.m.gz"][[4]], PlotRange->All]


xinfo=Import["ProbAmp_L7_cdMode4.m"][[1;;3]]


ACL7cd3=Import["ProbAmp_L7_cdMode3.m.gz"][[4]];


x=Import["ProbAmp_L7_cdMode4.m"][[4]]


pPlot=Table[{tVals[[j]],Sum[Abs[x[[k,j]]^2],{k,1,400}]},{j,1,Length[tVals]}];
ListLinePlot[pPlot]


Manipulate[ListLinePlot[Table[Abs[x[[k,tt]]]^2,{k,1,200}], PlotRange->{0,1}],{tt,0,300,1}]


pPlot=Table[{tVals[[j]],Abs[x[[1,j]]^2]},{j,1,Length[tVals]}];
ListLinePlot[pPlot, PlotRange->{0,1}]


Trunc=300;


K=Table[{tVals[[j]],Sum[k*Abs[x[[k,j]]^2],{k,1,Trunc}]},{j,1,Length[tVals]}];


ListLinePlot[K]


MaxKrylovBound=Table[{tVals[[j]],Length[bs]*(1-Sum[Abs[x[[k,j]]^2],{k,1,Trunc}])},{j,1,Length[tVals]}];
KUpper=Table[{tVals[[j]], K[[j,2]]+MaxKrylovBound[[j,2]]},{j,1,Length[tVals]}];
KLower=Table[{tVals[[j]], K[[j,2]]-MaxKrylovBound[[j,2]]},{j,1,Length[tVals]}];


ListLinePlot[{K,KUpper,KLower}]


KL7cd3=Table[{tVals[[j]],Sum[k*Abs[ACL7cd3[[k,j]]^2],{k,1,Trunc}]},{j,1,Length[tVals]}];


ListLinePlot[KL7cd3]


c=4.//Precision
d=Rationalize[c]//Precision
N[d]//Precision



