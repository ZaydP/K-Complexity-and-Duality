#!/usr/bin/env wolframscript
(* ::Package:: *)

L=5;
J=1;
nmax=130;
Prec=130;


(* ::Section::Closed:: *)
(*Functions to take tensor products of local operators*)


(* ::Input::Initialization:: *)
(*Given a binary number occupation rep of spins, it gives the corresponding matrix by taking a tensor product along the chain*)
ChainOccRepTensorProduct[OccRep_, Spin_] := Module[{a}, 
a=IdentityMatrix[1];
Do[ If[OccRep[[i]]==0, a = KroneckerProduct[a,IdentityMatrix[Dimensions[Spin]]]];
If[OccRep[[i]]==1, a = KroneckerProduct[a, Spin]];
,{i, 1, Length[OccRep]}];
a]


(* ::Input::Initialization:: *)
(*Takes Tensor product along the chain, with operator at position Position, and identities everywhere else*)
ChainTensorProduct[Length_,Operator_,Position_]:=Module[{list,chain},
list=ConstantArray[IdentityMatrix[2],Length];list[[Position]]=Operator;
chain=Apply[KroneckerProduct,list]]


(* ::Input::Initialization:: *)
(*Takes Tensor product along the chain, with operators at position Position and Position +1, and identities everywhere else*)
ChainNeighbourTensorProduct[Length_,Operator1_,Operator2_,Position_,PBC_:False]:=Module[{list,chain}, 
list=ConstantArray[IdentityMatrix[2],Length];
If[Position==Length, list[[Position]]=Operator1;list[[1]]=Operator2,list[[Position]]=Operator1;list[[Position+1]]=Operator2];
chain=Apply[KroneckerProduct,list]]


(* ::Input::Initialization:: *)
JordanWignerTensorProduct[Length_,StringOperator_,Operator_,Position_]:=Module[{list,chain},
list=Join[ConstantArray[StringOperator,Position-1],{Operator},ConstantArray[IdentityMatrix[2],Length-Position]];
chain=Apply[KroneckerProduct,list]]


(* ::Section::Closed:: *)
(*Functions to make Matrices for Ising Chain Hamiltonian and Operators*)


(* ::Input::Initialization:: *)
HIsingIntMode[n_]:=-J*ChainNeighbourTensorProduct[L,PauliMatrix[3],PauliMatrix[3],n];
HIsingBMode[n_]:=-J*ChainTensorProduct[L,PauliMatrix[1],n];


(* ::Input::Initialization:: *)
HIsing = Sum[HIsingIntMode[n],{n,1,L-1}]+Sum[HIsingBMode[n],{n,1,L}];
MatrixForm[HIsing];
HIsingNormed = HIsing/Sqrt[DocInnerProd[HIsing,HIsing]];


(* ::Subsection:: *)
(*Operators*)


(* ::Input::Initialization:: *)
\[Sigma]Plus0={{0,1},{0,0}};
\[Sigma]Minus0={{0,0},{1,0}};


\[Sigma]zMode[n_] := ChainTensorProduct[L, PauliMatrix[3], n];
\[Sigma]yMode[n_] := ChainTensorProduct[L, PauliMatrix[2], n];
\[Sigma]xMode[n_] := ChainTensorProduct[L, PauliMatrix[1], n];
\[Sigma]PlusMode[n_] := ChainTensorProduct[L, \[Sigma]Plus0, n];
\[Sigma]MinusMode[n_] := ChainTensorProduct[L, \[Sigma]Minus0, n];


(* ::Input::Initialization:: *)
\[Sigma]z = Sum[\[Sigma]zMode[n],{n,1,L}];
\[Sigma]y = Sum[\[Sigma]yMode[n],{n,1,L}];
\[Sigma]x = Sum[\[Sigma]xMode[n],{n,1,L}];
\[Sigma]Plus = Sum[\[Sigma]PlusMode[n],{n,1,L}]; 
\[Sigma]Minus = Sum[\[Sigma]MinusMode[n],{n,1,L}]; 


(* ::Subsection::Closed:: *)
(*Non-local Operators*)


(* ::Text:: *)
(*Subscript[c, i]= Subscript[\[CapitalPi], j<i] Subscript[(\[Sigma]^z), j] Subscript[(SuperMinus[\[Sigma]]), i]*)
(*Subscript[SuperPlus[c], i]= Subscript[\[CapitalPi], j<i] Subscript[(\[Sigma]^z), j] Subscript[(SuperPlus[\[Sigma]]), i]*)


(* ::Input::Initialization:: *)
ccMode[n_]:=JordanWignerTensorProduct[L,PauliMatrix[3],\[Sigma]Minus0,n];
cdMode[n_]:=JordanWignerTensorProduct[L,PauliMatrix[3],\[Sigma]Plus0,n];


ccc=Sum[ccMode[n],{n,1,L}];
cdc = Sum[cdMode[n], {n, 1, L}];
MatrixForm[ccc]
MatrixForm[cdc]


(* ::Section:: *)
(*Functions for Lanczos Coeffs, and Probability Amplitudes*)


DocInnerProd[A_,B_]:=Conjugate[Flatten[A]] . Flatten[B]


comm[A_,B_]:=A . B-B . A


(*Generate set with repeated actions of Liouvillian on state Op*)
InitialBasis[H_,Op_,n_,prec_]:=
Module[{B},B=N[NestList[comm[H,#]&,Op,n],prec]]


(*Generate Othogonal Basis*)
OrthoBasis[H_,Op_,n_,prec_]:=
Module[{B},B=Orthogonalize[InitialBasis[H,Op,n-1,prec],DocInnerProd, Method->"GramSchmidt",Tolerance->10^(-10)]]


(*Extract b_ns from the off-diagonals of the Louivillian in the Krylov Basis*)
LanczosOrthog[H_,Op_,nmax_,prec_]:=
Module[{B,a,b},
B=OrthoBasis[H,Op,nmax,prec];
a=Table[DocInnerProd[B[[n]],comm[H,B[[n]]]],{n,1,nmax-1}];
b=Table[DocInnerProd[B[[n-1]],comm[H,B[[n]]]],{n,2,nmax-1}];
{a,b}]


LanczosDoc[H_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{b,Op1,Op2,Op3,B,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i,Ops},
(*In order to improve efficiency we initialize and pre-allocate arrays of length nmax to prevent dynamic resizing of arrays in the loop*)
(*Op and LOp are arrays of matrices which are the dimension of S*)
count=1;
b=ConstantArray[0.,nmax];
Ops=ConstantArray[0.,{nmax,Length[H],Length[H]}];

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

(*If[Precision[tolval]<2,Print["Precision of tolval is less than 2"]];*)
If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
b[[i]]=N[B/Sqrt[NormSq1*NormSq2],bPrecision];
Ops[[i]]=Op1;

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;

count++;
(*If[Mod[count,100]==0,Print["Count is: ",count]];*)
,{i,2,nmax}];
{Table[b[[i]],{i,1,count}],Table[Ops[[i]],{i,1,count}]}]


orthout=LanczosOrthog[HIsing,\[Sigma]z,nmax,Prec]


docout=LanczosDoc[HIsing,\[Sigma]z,nmax,Prec,0.00001,20]


docout[[1]]


ListLinePlot[docout[[1]]]


OpOrthogTable=Table[DocInnerProd[docout[[2]][[i]],docout[[2]][[j]]],{i,1,42},{j,1,42}]//MatrixForm;


OpOrthogList=Table[DocInnerProd[docout[[2]][[4]],docout[[2]][[j]]],{j,1,Length[docout[[1]]]}]
ListLinePlot[OpOrthogList]



