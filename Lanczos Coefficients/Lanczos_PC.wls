#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];


Get["C:/Users/daisy/OneDrive - University of Cape Town/0. Masters Project/Ising_matrices_package.wl"]


(* ::Section:: *)
(*Lanczos Algorithm*)


comm[Op1_,Op2_]:= Module[{com},com= Op1 . Op2 - Op2 . Op1]


DocInnerProd[A_,B_] := Module[{psi},psi=Conjugate[Flatten[A]] . Flatten[B]]


SubInnerProd[A_,B_] := Module[{subGroupFactor,psi},
subGroupFactor=Apply[KroneckerProduct,Table[PauliMatrix[3],Log[2,Length[A]]]];
psi=Conjugate[Flatten[A]] . Flatten[B . subGroupFactor]]+Conjugate[Flatten[A]] . Flatten[B]


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


LanczosGeneralized[H0_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{a,b,Op1,Op2,Op3,A,B,H,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i,Ops},
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
(*Chop[LOp,10^-(bPrec-1)];*)


A=DocInnerProd[Op2,LOp]/(NormSq2);
B=DocInnerProd[Op1,LOp];
tolval=B/NormSq1;


(*Chop[tolval,10^-(bPrecision-1)]*);
Op3=LOp-tolval*Op1-A*Op2;
(*If[Precision[tolval]<2,Print["Precision of tolval is less than 2"]];*)
If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
a[[i]]=N[A,bPrecision];
b[[i]]=N[B/Sqrt[NormSq1*NormSq2],bPrecision];

(*b[[i]]=Chop[B/Sqrt[NormSq1*NormSq2],10^-(bPrec-1)];*)

Ops[[i]]=Op1;

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;
count++;
If[Mod[count,50]==0,Print["Count is: ",count]];
,{i,2,nmax}];
{Table[a[[i]],{i,1,count}],Table[b[[i]],{i,1,count+1}],Table[Ops[[i]],{i,1,count}]}]


LanczosGeneralizedCorrected[H0_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{a,b,Op1,Op2,Op3,A,B,H,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i,Ops},
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
A=DocInnerProd[Op1,LOp]/(NormSq1);a[[1]]=A;

Op2=LOp - A*Op1;
Echo[Precision@Op1];
Do[(*A couple of dummy variables, to reduce repeat calculations. Worse space efficiency though*)
LOp=comm[H,Op2];

NormSq2=DocInnerProd[Op2,Op2];
(*Chop[LOp,10^-(bPrec-1)];*)


A=DocInnerProd[Op2,LOp]/(NormSq2);
B=DocInnerProd[Op1,LOp];
tolval=B/NormSq1;


(*Chop[tolval,10^-(bPrecision-1)]*);
Op3=LOp-tolval*Op1-A*Op2;
(*If[Precision[tolval]<2,Print["Precision of tolval is less than 2"]];*)
If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
a[[i]]=N[A,bPrecision];
b[[i]]=N[B/Sqrt[NormSq1*NormSq2],bPrecision];

(*b[[i]]=Chop[B/Sqrt[NormSq1*NormSq2],10^-(bPrec-1)];*)

(*Ops[[i]]=Op1;*)

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;

(*Echo[Precision@Op1];*)
count++;
If[Mod[count,50]==0,Print["Count is: ",count]];
,{i,2,nmax}];
{Table[a[[i]],{i,1,count}],Table[b[[i]],{i,1,count+1}],Table[Ops[[i]],{i,1,count}]}]


LanczosGeneralizedUnnormed[H0_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{a,b,Op1,Op2,Op3,A,B,H,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i,Ops},
(*In order to improve efficiency we initialize and pre-allocate arrays of length nmax to prevent dynamic resizing of arrays in the loop*)
(*Op and LOp are arrays of matrices which are the dimension of S*)
count=1;
a=ConstantArray[0.,nmax];
b=ConstantArray[0.,nmax];
Ops=ConstantArray[0.,{nmax,Length[H],Length[H]}];

H=N[H0,precision];
Op1=N[S,precision]/Sqrt[DocInnerProd[S,S]];
LOp=comm[H,Op1];
A=DocInnerProd[Op1,LOp];a[[1]]=A;

Op2=LOp - A*Op1;
B=Sqrt[DocInnerProd[Op2,Op2]];b[[2]]=B;
Op2=Op2/B;
LOp=comm[H,Op2];
A=DocInnerProd[Op2,LOp];a[[2]]=A;
Do[

Op3 = LOp - A*Op2 - B*Op1;
(*Op3=U_n*)

B=Sqrt[DocInnerProd[Op3,Op3]];
b[[i+1]]=N[B,bPrecision];


If[N[b[[i+1]]]<btol,Print["Cut-off reached at n= ",i,", since b and btol are ",b[[i+1]]," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
Ops[[i]]=Op1;

Op1=Op2;
Op2=Op3/B;

LOp=comm[H,Op2];

A=DocInnerProd[Op2,LOp];
a[[i+1]]=N[A,bPrecision];

count++;
If[Mod[count,50]==0,Print["Count is: ",count]];
,{i,2,nmax}];
{Table[a[[i]],{i,1,count+1}],Table[b[[i]],{i,1,count+1}], Table[Ops[[i]],{i,1,count}]}]


LanczosOthogonalize[H0_, O0_, nmax_, nPrec_]:=
Module[{Basis0, Basis,H,Op,as,bs,HKryl},H=N[H0,nPrec];Op=N[O0,nPrec];
Basis0=NestList[comm[H,#]&,Op,nmax];
Basis=Orthogonalize[Basis0,DocInnerProd];
bs=Table[DocInnerProd[Basis[[i]],comm[H,Basis[[i+1]]]],{i,1,nmax-1}];
as=Table[DocInnerProd[Basis[[i]],comm[H,Basis[[i]]]],{i,1,nmax}];
HKryl = Table[DocInnerProd[Basis[[i]],comm[H,Basis[[j]]]],{i,1,nmax},{j,1,nmax}];
{as,bs}]


(* ::Section::Closed:: *)
(*Nonsense*)


LanczosGeneralizedTesting[H0_,S_,nmax_,precision_,btol_,bPrecision_] := Module[{a,b,Op1,Op2,Op3,A,B,H,LOp ,NormSq1,NormSq2,tolval, count,OpDim, i,Ops},
(*In order to improve efficiency we initialize and pre-allocate arrays of length nmax to prevent dynamic resizing of arrays in the loop*)
(*Op and LOp are arrays of matrices which are the dimension of S*)
count=1;
a=ConstantArray[0.,nmax];
b=ConstantArray[0.,nmax];
Ops=ConstantArray[0.,{nmax,Length[H],Length[H]}];
time1=ConstantArray[0.,{nmax,2}];
timeDIP=ConstantArray[0.,{nmax,2}];
timeOp3=ConstantArray[0.,{nmax,2}];
mem=ConstantArray[0,nmax];
H=N[H0,precision];
Op1=N[S,precision];
LOp=comm[H,Op1];
NormSq1=DocInnerProd[Op1,Op1];
Op2=LOp;

Do[(*A couple of dummy variables, to reduce repeat calculations. Worse space efficiency though*)

mem[[i]]=MemoryInUse[];

time1[[i]]=AbsoluteTiming[
LOp=comm[H,Op2];];

timeDIP[[i]]=AbsoluteTiming[
NormSq2=DocInnerProd[Op2,Op2];
Chop[LOp,10^-(bPrecision-1)];
Chop[NormSq1,10^-(bPrecision-1)];

A=DocInnerProd[Op2,LOp]/(NormSq2);
B=DocInnerProd[Op1,LOp];];
tolval=B/NormSq1;

If[i==10,Echo[tolval, Chop[tolval,10^-(bPrecision-1)]] ];
Chop[tolval,10^-(bPrecision-1)];
Chop[A,10^-(bPrecision-1)];
Chop[B,10^-(bPrecision-1)];


timeOp3[[i]]=AbsoluteTiming[
Op3=LOp-tolval*Op1-A*Op2;];

Chop[Op3,10^(bPrecision-1)];

(*If[Precision[tolval]<2,Print["Precision of tolval is less than 2"]];*)
If[N[tolval]<btol,Print["Cut-off reached at n= ",i,", since tolval (\!\(\*FractionBox[\(b[i]\), \(NormSq[i - 1]\)]\)) and btol are ",tolval," and ",btol]; 
(* Break when the (approximate) Lanczos coefficients reach close to zero. We want to catch this and break before we hit a div by zero error in the b[[i]]*)
Break[]];
a[[i]]=N[A,bPrecision];
(*b[[i]]=N[B/Sqrt[NormSq1*NormSq2],bPrecision];*)

b[[i]]=Chop[B/Sqrt[NormSq1*NormSq2],10^-(bPrec-1)];

(*Ops[[i]]=Op1;*)

Op1=Op2;
Op2=Op3;
NormSq1=NormSq2;
(*Print[i];*)
count++;
(*If[Mod[count,100]==0,Print["Count is: ",count]];*)
,{i,2,nmax}];
Echo[Total[time1[[;;count,1]]],"1"];
Echo[Total[timeDIP[[;;count,1]]],"2"];
Echo[Total[timeOp3[[;;count,1]]],"3"];
Echo[mem[[count]]-mem[[2]]];

{Table[a[[i]],{i,1,count}],Table[b[[i]],{i,1,count+1}]}]


(* ::Section:: *)
(*Inputs*)


L=5;
nmax=1000;
prec=1000;
btol=10^(-6);
bPrec=400;
Info = "L=5,nmax=1000,prec=1000,btol=10^(-6),bPrec=400";


p1=PauliMatrix[1];p2=PauliMatrix[2];p3=PauliMatrix[3];id=IdentityMatrix[2];
sp={{0,1},{0,0}};
sm={{0,0},{1,0}};
spAlt=(-p3+I*p2)/2;
smAlt=(-p3-I*p2)/2;


Clear[RanNormHermOp]


RanNormHermOp[precision_]:=Module[{r,rHerm,rHermNorm},r=Table[RandomReal[WorkingPrecision->precision],{2},{2}];rHerm=r+ConjugateTranspose[r];rHermNorm=rHerm/Sqrt[DocInnerProd[rHerm,rHerm]]]


ran1=RanNormHermOp[prec];


ran2=RanNormHermOp[prec];


ran3=RanNormHermOp[prec];


ran4=RanNormHermOp[prec];


ran5=RanNormHermOp[prec];


H=TFIMHamiltonian[L,1,1];


HAlt=TFIMAltHamiltonian[L,1,1,True];


N@ran5


Op=ccMode[L,1];


Op1List={ran1,ran2,ran3,ran4,ran5};
Op1=Apply[KroneckerProduct,Op1List];


(*Calculate and Plot Lanczos coefficients*)
{a1,b1,ops1}=EchoTiming[LanczosGeneralizedCorrected[HAlt, Op1,nmax,prec,btol,bPrec], "Lanczos Algorithm Time: "];
Print["b1's complete: "]
Beep[]


binfo={"TFIMAlt periodic, different random normalized hermitian operators on sites 1, 2, 3, 4, 5"};


(*binfo=StringJoin["TFIMAlt periodic, different random normalized hermitian operators on sites 1, 2, 3";*)
Export["bs_Lanczos_PC_L5TFIMAlt_periodic_5site_RanNormHerm12345.txt",ToString[{Info,binfo,a1,b1,{ran1,ran2,ran3,ran4,ran5}}, InputForm]];


ListLinePlot[{a1,b1}]
(*ListLinePlot[Table[DocInnerProd[ops1[[i]],ops1[[2]]],{i,2,6}]]*)


x=ToExpression@Import["C:\\Users\\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Lanczos Coefficients\\bs_Lanczos_PC_L5TFIMAlt_periodic_5site_RanNormHerm12345.txt"];


(*Calculate and Plot Lanczos coefficients*)
{a2,b2,ops2}=EchoTiming[LanczosGeneralizedUnnormed[HAlt, Op,nmax,prec,btol,bPrec], "Lanczos Algorithm Time: "];
Print["b2's complete: "]
Beep[]


ListLinePlot[{a2,b2}]
(*ListLinePlot[Table[DocInnerProd[ops2[[i]],ops2[[2]]],{i,2,6}]]*)


{a3,b3}=LanczosOthogonalize[HAlt,Op1,60,100];b3=Prepend[b3,0.];//EchoTiming
ListLinePlot[{a3,b3}]


Precision[b1[[2;;27]]]


b1[[20;;27]]//N
b3[[20;;27]]//N


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
Do[  cc[m,n+1]  = N[( If[m>= 1, I*cc[m-1,n], 0]  -If[m<= n, as[[n+1]]*cc[m,n], 0]   - If[m<=n -1&& n>0, bs[[n]]*cc[m,n-1] , 0 ])/bs[[n+1]], Prec]  , {n, 0, \[Phi]Pts-1}, {m, 0, n+1}]//EchoTiming; 
PTab = Table[  Sum[ cc[m,kk]*NNNDTab[[m+1]], {m, 0, kk}] , {kk, 0, \[Phi]Pts}]//EchoTiming]


ReturnAmpJaco[H_, Operator_, NPrec_] := Module[{Nvals, Nvecs, VecInv},
EchoTiming[{Nvals, Nvecs} = Eigensystem[N[H, NPrec]]];
VecInv = Inverse[Nvecs];
Echo@MemoryInUse[];
EchoTiming[NNNexpr = Chop[Expand[DocInnerProd[Operator, ConjugateTranspose[Nvecs] . MatrixExp[-I*t*DiagonalMatrix[Nvals]] . ConjugateTranspose[VecInv] . Operator . VecInv . MatrixExp[I*t*DiagonalMatrix[Nvals]] . Nvecs]], 10^(-NPrec + 3)]];
Echo@MemoryInUse[];
NNNnorm = NNNexpr /. t -> 0;
NNNexpr = NNNexpr/NNNnorm]


retamp=ToExpression[Import["C:\\Users\\daisy\\OneDrive - University of Cape Town\\0. Masters Project\\Return Amplitudes\\ReturnAmpJacoL5TFIM_periodic_sp1.txt"]];


retamp=ReturnAmpJaco[HAlt,Op,50];


retamp2=ReturnAmpJaco[HAlt,Op,150];


a1=Table[0,Length[b1]];


tMax=50;
tPts=150;
Prec=120;
tVals = Table[N[tMax*(i)/(tPts),Prec], {i, 1, tPts}];


ac=AmplitudeCoefficientScript[a1,b1,retamp2,tVals,Prec];


Length[ac]


prob=Table[{tVals[[j]],Total[Abs[ac[[1;;152,j]]]^2]},{j,1,tPts}];
ListLinePlot[prob]


a1[[3]]//Precision


retamp//Precision


ListLinePlot[Table[{tVals[[j]],Total[Abs[ac[[1;;29,j]]]^2]},{j,1,tPts}]]


Total[Abs[ac[[1;;3,5]]]^2]


ListLinePlot[prob]


prob[[35]]


retamp2=ReturnAmpJaco[HAlt,ccMode[L,1],30];


retamp3=ReturnAmpJaco[HAlt,ccMode[L,2],30];


Plot[{Abs@retamp2,Abs@retamp3},{t,0,10}]
