#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
MyInnerProd[A_,B_] := Module[{psi},psi = Tr[ConjugateTranspose[A] . B]]


(* ::Input::Initialization:: *)
comm[Op1_,Op2_]:= Module[{com},com= Op1 . Op2 - Op2 . Op1]


DocInnerProd[A_,B_] := Module[{psi},psi=Conjugate[Flatten[A]] . Flatten[B]]



(* ::Input::Initialization:: *)
LanczosDoc[H_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{b,Op1,Op2,Op3,B,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i},
(*In order to improve efficiency we initialize and pre-allocate arrays of length nmax to prevent dynamic resizing of arrays in the loop*)
(*Op and LOp are arrays of matrices which are the dimension of S*)
count=0;
b=ConstantArray[0.,nmax];


Op1=N[S,precision];
LOp=comm[H,Op1];
NormSq1=DocInnerProd[Op1,Op1];
Op2=LOp;

Do[(*A couple of dummy variables, to reduce repeat calculations. Worse space efficiency though*)
LOp=comm[H,Op2];
NormSq2=DocInnerProd[Op2,Op2];
B=DocInnerProd[Op1,LOp];
tolval=B/NormSq1;
Op3=LOp-tolval*Op1;

If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
b[[i]]=N[B/Sqrt[NormSq1*NormSq2],precision];

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;

count++;
,{i,2,nmax}];
Table[b[[i]],{i,1,count}]]


(* ::Input::Initialization:: *)
(*For Ising Matrix. Given a binary number occupation rep of spins, it gives the corresponding matrix*)
SpinStateMatrix[OccRep_, Spin_] := Module[{a}, 
a=IdentityMatrix[1];
Do[ If[OccRep[[i]]==0, a = KroneckerProduct[a,IdentityMatrix[Dimensions[Spin]]]];
If[OccRep[[i]]==1, a = KroneckerProduct[a, Spin]];
,{i, 1, Length[OccRep]}];
a]


(* ::Input::Initialization:: *)
(*For Ising Matrix. Given Length, Interaction Spins, Number of Interaction Spins (k-local), And the Spins for 2 magnetic fields, it Generates the Matrix representation.*)
IsingSum[L_,IntSpin_,IntSpinNo_,BSpin1_,BSpin2_]:=Module[{IntPerms, IntPieces, IntPart, B1Part, B2Part, SingleSpinPerms, B1Pieces, B2Pieces},
IntPerms = Table[Join[ConstantArray[0,i-1],ConstantArray[1,IntSpinNo],ConstantArray[0,L-(i-1)-IntSpinNo]],{i,1,L-IntSpinNo+1}];
SingleSpinPerms=Table[Join[ConstantArray[0,i-1],ConstantArray[1,1],ConstantArray[0,L-(i-1)-1]],{i,1,L}];
If[IntSpin == 0, IntPart = IntSpin, IntPart =IntSpin,
IntPieces = Table[SpinStateMatrix[IntPerms[[i]], IntSpin],{i, 1, Length[IntPerms]}];
IntPart = Total[IntPieces]];
If[BSpin1 == 0, B1Part=BSpin1,B1Part=BSpin1,
B1Pieces=Table[SpinStateMatrix[SingleSpinPerms[[i]],BSpin1],{i,1,Length[SingleSpinPerms]}];
B1Part = Total[B1Pieces]];
If[BSpin2 == 0, B2Part = BSpin2,B2Part=BSpin2,
B2Pieces =Table[SpinStateMatrix[SingleSpinPerms[[j]],BSpin2],{j,1,Length[SingleSpinPerms]}];
B2Part = Total[B2Pieces]];
IntPart+B1Part+B2Part]


(* ::Input::Initialization:: *)
NormKComplexityFromBs[b_,bmax_,tmax_]:=Module[{\[Phi],soln,bNew, bSlice,f,K,NormK},
(*Solve the particle hopping equation with Zero BC on n, and IC of Subscript[\[Phi], n](0)= Subscript[\[Delta], n](0)*)

bSlice=b[[;;bmax+1]];
bNew= Append[bSlice,0];
soln=NDSolve[Join[Table[D[\[Phi][n,t],t]==
 bNew[[n+1]]*\[Phi][n-1,t]- bNew[[n+2]]*\[Phi][n+1,t],{n,0,bmax}],
Table[\[Phi][n,0]==KroneckerDelta[n,0],{n,0,bmax}]],{Table[\[Phi][n,t],{n,0,bmax}]},{t,0,tmax}];
(*Extract the phi's from the solution*)
f=soln[[1]][[1]][[2]];
(*Calculate the Normalized K-Complexity from the phi's*)
NormK=Sum[n*Conjugate[f[[n+1]]]*f[[n+1]],{n,1,bmax}]/Sum[Conjugate[f[[n+1]]]*f[[n+1]],{n,1,bmax}]
]


(* ::Input::Initialization:: *)
L=4;
nmax=200;
prec=130;
tMax=50;
tStart=0.1;(*only used for KDiscrete*)
tInterval=1.0;(*interval for KDiscrete*)
btol=10^(-6);
bPrec=15;
H=IsingSum[L,PauliMatrix[3],2,PauliMatrix[3],PauliMatrix[1]];
(*Opx has is on x sites*)
Op1=IsingSum[L,0,2,PauliMatrix[3],0];
Op2=IsingSum[L,PauliMatrix[3],2,0,0];
Op3=IsingSum[L,PauliMatrix[3],3,0,0];

(*Calculate and Plot Lanczos coefficients*)
b1=LanczosDoc[H,Op1,nmax,prec,btol,bPrec];
(*Print["b1's complete: ", b1]*)
(*b2=LanczosDoc[H,Op2,nmax,prec,btol,bPrec];
Print["b2's complete: ", b2]
b3=LanczosDoc[H,Op3,nmax,prec,btol, bPrec];
Print["b3's complete: ", b3]
Plotb=ListLinePlot[{b1,b2,b3},AxesLabel->{"n","b_n"},PlotLabels->{"1-site","2-site","3-site"},PlotLabel->"L=7"]
Export["bs_KComplexity_test_output.jpg",Plotb]*)

(*Calculate and Plot K-Complexity*)
K1 = NormKComplexityFromBs[b1,Length[b1]-1,tMax];
K1Discrete=Table[K1,{t,tStart,tMax,tInterval}];
ListLinePlot[K1Discrete]
ts=Table[t,{t,tStart,tMax,tInterval}]
Points=Transpose[Join[{ts},{K1Discrete}]]

(*Plot[{K1,Interpolation[K1Discrete]},{t,0,tMax}]*)
(*Plot[K1-Interpolation[K1Discrete],{t,tStart,tMax}]*)
Print["K1 complete. K1 with Table params: ",{tStart,tMax,tInterval}]
Print[K1Discrete]
(*K2 = NormKComplexityFromBs[b2,Length[b2]-1,tmax];
Print["K2 complete: ", K2]
K3 = NormKComplexityFromBs[b3,Length[b3]-1,tmax];
Print["K3 complete: ", K3]
plotK=Plot[{K1,K2,K3},{t,0,tmax},PlotRange->{Full,{0,Automatic}},AxesLabel->{"t","K-Complexity"}, PlotLabels->{"1-site","2-site","3-site"},PlotLabel->"L=7"]
Export["Ks_KComplexity_test_output.jpg",plotK]*)






