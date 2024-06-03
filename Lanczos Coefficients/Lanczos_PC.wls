#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];


Get["C:/Users/daisy/OneDrive - University of Cape Town/0. Masters Project/Ising_matrices_package.wl"]


(* ::Section:: *)
(*Lanczos Algorithm*)


comm[Op1_,Op2_]:= Module[{com},com= Op1 . Op2 - Op2 . Op1]


DocInnerProd[A_,B_] := Module[{psi},psi=Conjugate[Flatten[A]] . Flatten[B]]


LanczosDoc[H0_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{b,Op1,Op2,Op3,B,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i},
(*In order to improve efficiency we initialize and pre-allocate arrays of length nmax to prevent dynamic resizing of arrays in the loop*)
(*Op and LOp are arrays of matrices which are the dimension of S*)
count=0;
b=ConstantArray[0.,nmax];
H=N[H0,precision];
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
Print[i];
count++;
If[Mod[count,10]==0,Print["Count is: ",count]];
,{i,2,nmax}];
Table[b[[i]],{i,1,count+1}]]


LanczosGeneralized[H0_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{a,b,Op1,Op2,Op3,A,B,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i,Ops},
(*In order to improve efficiency we initialize and pre-allocate arrays of length nmax to prevent dynamic resizing of arrays in the loop*)
(*Op and LOp are arrays of matrices which are the dimension of S*)
count=1;
a=ConstantArray[0.,nmax];
b=ConstantArray[0.,nmax];
Ops=ConstantArray[0.,{nmax,Length[H],Length[H]}];

H=N[H0,precision];
Op1=N[S,precision];
LOp=comm[H,Op1];
NormSq1=DocInnerProd[Op1,Op1];
Op2=LOp;

Do[(*A couple of dummy variables, to reduce repeat calculations. Worse space efficiency though*)
LOp=comm[H,Op2];
NormSq2=DocInnerProd[Op2,Op2];

A=DocInnerProd[Op2,LOp]/(NormSq2);
B=DocInnerProd[Op1,LOp];
tolval=B/NormSq1;

Op3=LOp-tolval*Op1-A*Op2;

(*If[Precision[tolval]<2,Print["Precision of tolval is less than 2"]];*)
If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
a[[i]]=N[A,bPrecision];
b[[i]]=N[B/Sqrt[NormSq1*NormSq2],bPrecision];
Ops[[i]]=Op1;

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;
Print[i];
count++;
(*If[Mod[count,100]==0,Print["Count is: ",count]];*)
,{i,2,nmax}];
{Table[a[[i]],{i,1,count}],Table[b[[i]],{i,1,count+1}],Table[Ops[[i]],{i,1,count}]}]


(* ::Section:: *)
(*Inputs*)


L=7;
nmax=100;
prec=300;
btol=10^(-6);
bPrec=40;
Info = "L=4,nmax=230,prec=130,btol=10^(-6),bPrec=30";


p1=PauliMatrix[1];p2=PauliMatrix[2];p3=PauliMatrix[3];id=IdentityMatrix[2];
sp={{0,1},{0,0}};
sm={{0,0},{1,0}};


Op1List={id,id,sp,sp,id,id,id};
Op1=Apply[KroneckerProduct,Op1List];


Op2List={id,sp,id,sp,id,id,id};
Op2=Apply[KroneckerProduct,Op2List];


H=TFIMHamiltonian[L,1,1];
p1=Op1;


EchoTiming[b1=LanczosGeneralized[H,Op1,nmax,prec,btol,bPrec]];


EchoTiming[b2=LanczosGeneralized[H,Op2,nmax,prec,btol,bPrec]];


(*b1info="H=IsingSum[L,PauliMatrix[3],2,PauliMatrix[3],PauliMatrix[1]];
p1=IsingSum[L,0,2,PauliMatrix[1],0];";
Export["bs_Lanczos_locality_HPC_L4_Op1.txt",{Info,b1info,b1},"Table"];*)



ListLinePlot[{b1[[2]],b2[[2]]}]


b2[[2]]


Table[Chop[N[DocInnerProd[b1[[3]][[i]],b1[[3]][[j]]],5]],{i,120,140},{j,100,130}]//MatrixForm


b1


Import["C:\\Users\\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Lanczos Coefficients\\bs_Lanczos_HPC_L8_cdMode3_Test.txt"]


{bs}=Import["C:\\Users\\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Lanczos Coefficients\\bs_Lanczos_HPC_L8_cdMode3_Test.txt","Table", HeaderLines->2];


bs


ListLinePlot[bs[[2]]]



