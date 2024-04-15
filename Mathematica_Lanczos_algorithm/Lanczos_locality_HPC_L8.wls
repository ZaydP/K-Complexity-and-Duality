#!/usr/bin/env wolframscript
(* ::Package:: *)

comm[Op1_,Op2_]:= Module[{com},com= Op1 . Op2 - Op2 . Op1]


DocInnerProd[A_,B_] := Module[{psi},psi=Conjugate[Flatten[A]] . Flatten[B]]


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

If[Precision[tolval]<2,Print["Precision of tolval is less than 2"]];
If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
b[[i]]=N[B/Sqrt[NormSq1*NormSq2],bPrecision];

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;

count++;
,{i,2,nmax}];
Table[b[[i]],{i,1,count}]]


(*For Ising Matrix. Given a binary number occupation rep of spins, it gives the corresponding matrix*)
SpinStateMatrix[OccRep_, Spin_] := Module[{a}, 
a=IdentityMatrix[1];
Do[ If[OccRep[[i]]==0, a = KroneckerProduct[a,IdentityMatrix[Dimensions[Spin]]]];
If[OccRep[[i]]==1, a = KroneckerProduct[a, Spin]];
,{i, 1, Length[OccRep]}];
a]


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


SetDirectory[NotebookDirectory[]];
L=8;
nmax=32800;
prec=32800;
btol=10^(-6);
bPrec=20;
H=IsingSum[L,PauliMatrix[3],2,PauliMatrix[3],PauliMatrix[1]];
(*Opx has is on x sites*)
Op1=IsingSum[L,0,2,PauliMatrix[3],0];
Op2=IsingSum[L,PauliMatrix[3],2,0,0];
Op3=IsingSum[L,PauliMatrix[3],3,0,0];

(*Calculate and Plot Lanczos coefficients*)
b1timed=AbsoluteTiming[LanczosDoc[H,Op1,nmax,prec,btol,bPrec]];
b1=b1timed[[2]];
Print["b1's complete"]
Print["Time taken for operator b1s:", b1timed[[1]]]
(*b2=LanczosDoc[H,Op2,nmax,prec,btol,bPrec];
Print["b2's complete"]
b3=LanczosDoc[H,Op3,nmax,prec,btol, bPrec];
Print["b3's complete"]
Plotb=ListLinePlot[{b1,b2,b3},AxesLabel->{"n","b_n"},PlotLabels->{"1-site","2-site","3-site"},PlotLabel->Row[{"L=", L}]]*)

Info="L=8,nmax=32800,prec=32800,btol=10^(-6),bPrec=20";
Export["bs_Lanczos_locality_HPC_L8.jpg",Plotb];

b1info="H=IsingSum[L,PauliMatrix[3],2,PauliMatrix[3],PauliMatrix[1]];
Op1=IsingSum[L,0,2,PauliMatrix[3],0];";
Export["bs_Lanczos_locality_HPC_L8_Op1.txt",{Info,b1info,b1},"Table"];

(*b2info="H=IsingSum[L,PauliMatrix[3],2,PauliMatrix[3],PauliMatrix[1]];
Op2=IsingSum[L,PauliMatrix[3],2,0,0];";
Export["bs_Lanczos_locality_HPC_L8_Op2.txt",{Info,b2info,b2},"Table"];

b3info="H=IsingSum[L,PauliMatrix[3],2,PauliMatrix[3],PauliMatrix[1]];
Op3=IsingSum[L,PauliMatrix[3],3,0,0];";
Export["bs_Lanczos_locality_HPC_L8_Op3.txt",{Info,b3info,b3},"Table"];*)
