#!/usr/bin/env wolframscript
(* ::Package:: *)

{as,bs}=Import["bs_Lanczos_HPC_L7_cdMode1.txt","Table", HeaderLines->2];
as=Chop[as];
bsLength=Length[bs];
bsPrec=Precision[bs[[3]]];


RetAmp=ToExpression[Import["ReturnAmpJacoL7TFIMcd1234.txt"]][[1]];
RetAmpPrecision=Precision[RetAmp];


tMax=50;
tPts=500;
Prec=200;
tVals = Table[N[tMax*(i)/(tPts),Prec], {i, 1, tPts}];


InfoPA="ProbAmp_L7_cdMode1";
InfoK="K_L7_cdMode1";
Info2=StringJoin["bs length, precision: ", ToString[bsLength], ", ", ToString[bsPrec], ".  Return Amplitude precision: ", ToString[RetAmpPrecision]];
Info3=StringJoin["tMax: ",ToString[tMax], "  tPts: ",ToString[tPts],"  Prec: ", ToString[Prec]];


AmplitudeCoefficientScript[as_,Bs_,ReturnAmp_,tVals_,Prec_]:=Module[{bs,\[Phi]Pts,NNND,NNNDTab,Genf,FtoNSwap,\[Phi]t,cc,PTab},
bs=Bs[[1+LengthWhile[Bs,#==0&];;-1-LengthWhile[Reverse@Bs,#==0&]]]; (*removes tail of zeros and initial zeros*)
 \[Phi]Pts = Length[bs];
(*"Here the derivatives are computed at the time-points.  "*)
	 NNND[0] = ReturnAmp;
	Do[NNND[j+1] =  D[NNND[j], t]    , {j, 0, \[Phi]Pts}  ];//EchoTiming;
	NNNDTab = ParallelTable[NNND[j]/.t->tVals, {j, 0, \[Phi]Pts}, Method->"CoarsestGrained"];//EchoTiming;
(*"Here we compute the coefficeints cc "*)
Clear[cc];   cc[0,0]=1; cc[0,-1] = 0;
(*CHANGED TO HAVE N ONLY GO TO \[Phi]Pts-1*)
Do[  cc[m,n+1]  = N[(If[m>= 1, I*cc[m-1,n], 0]  +If[m<= n, as[[n+1]]*cc[m,n], 0]   - If[m<=n -1&& n>0, bs[[n]]*cc[m,n-1] , 0 ])/bs[[n+1]], Prec]  , {n, 0, \[Phi]Pts-1}, {m, 0, n+1}]//EchoTiming; 
PTab = Table[  Sum[ cc[m,kk]*NNNDTab[[m+1]], {m, 0, kk}] , {kk, 0, \[Phi]Pts}]//EchoTiming]


EchoTiming[AC=AmplitudeCoefficientScript[as,bs,RetAmp,tVals,Prec];,"Probability Amplitude time"]
Export[StringJoin[InfoPA,".m.gz"],{InfoPA,Info2,Info3,AC}]


EchoTiming[K=Table[{tVals[[j]],Sum[(k-1)*Abs[AC[[k,j]]^2],{k,1,21}]},{j,1,Length[tVals]}];,"Complexity sum time: "]
EchoTiming[PSum=Table[{tVals[[j]],Sum[Abs[AC[[k,j]]^2],{k,1,21}]},{j,1,Length[tVals]}];,"Probability sum time: "]
Export[StringJoin[InfoK,".m.gz"],{InfoK,Info2,Info3,K,PSum}]



