(* ::Package:: *)

BeginPackage["IsingMatrixPackage`"];


TFIMHamiltonian::usage="TFIMHamiltonian[L,InteractionStrength,MagneticFieldStrength,PBC=False] generates matrix for Transverse Field Ising Model Hamiltonian with \[Sigma]z interaction and \[Sigma]x magnetic field. The strength constants have a negative coefficient. Open Boundary Conditions are default, set PBC=True for Periodic Boundary Conditions";
TFIMAltHamiltonian::usage="TFIMAltHamiltonian[L,InteractionStrength,MagneticFieldStrength,PBC=False] generates matrix for Transverse Field Ising Model Hamiltonian with \[Sigma]x interaction and -\[Sigma]z magnetic field. This is a (x,y,z)->(-z,y,x) rotation from the canonical form. The strength constants have a negative coefficient. Open Boundary Conditions are default, set PBC=True for Periodic Boundary Conditions";


(*Alt denotes that there has been a is a (x,y,z)->(-z,y,x) rotation from the canonical form *)


\[Sigma]x::usage="Sum of \[Sigma]x operators on each site over chain of length L";
\[Sigma]y::usage="Sum of \[Sigma]y operators on each site over chain of length L";
\[Sigma]z::usage="Sum of \[Sigma]z operators on each site over chain of length L";
\[Sigma]Plus::usage="Sum of \[Sigma]Plus {{0,1},{0,0}} operators on each site over chain of length L";
\[Sigma]Minus::usage="Sum of \[Sigma]Minus {{0,0},{1,0}} operators on each site over chain of length L";
cc::usage="Sum of annihilation(\[Sigma]Minus) Jordan-Wigner strings for each site over chain of length L";
cd::usage="Sum of creation (\[Sigma]Plus) Jordan-Wigner strings for each site over chain of length L";


(*\[Sigma]xMode::usage="\[Sigma]xMode[L_,n_] Produces operator representing a \[Sigma]x operator on the nth site of a chain of length L";
\[Sigma]yMode::usage="\[Sigma]yMode[L_,n_] Produces operator representing a \[Sigma]y operator on the nth site of a chain of length L";
\[Sigma]zMode::usage="\[Sigma]zMode[L_,n_] Produces operator representing a \[Sigma]z operator on the nth site of a chain of length L";
\[Sigma]PlusMode::usage="\[Sigma]PlusMode[L_,n_] Produces operator representing a \[Sigma]Plus operator on the nth site of a chain of length L";
\[Sigma]MinusMode::usage="\[Sigma]MinusMode[L_,n_] Produces operator representing a \[Sigma]Minus operator on the nth site of a chain of length L";
ccMode::usage="ccMode[L_,n_] Produces Jordan-Wigner string for an annihilation operator on the nth site of a chain of length L. It is a tensor product of \[Sigma]z operators for the first n-1 sites, \[Sigma]Minus on the nth site, and Identities for the last L-n sites";
cdMode::usage="ccMode[L_,n_] Produces Jordan-Wigner string for a creation operator on the nth site of a chain of length L. It is a tensor product of \[Sigma]z operators for the first n-1 sites, \[Sigma]Plus on the nth site, and Identities for the last L-n sites";*)


(* ::Section:: *)
(*PRIVATE CLASS*)


Begin["`Private`"];


(* ::Section:: *)
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
ChainNeighbourTensorProduct[Length_,Operator1_,Operator2_,Position_]:=Module[{list,chain}, 
list=ConstantArray[IdentityMatrix[2],Length];
If[Position==Length,
 list[[Position]]=Operator1;list[[1]]=Operator2,
list[[Position]]=Operator1;list[[Position+1]]=Operator2];
chain=Apply[KroneckerProduct,list]]


(* ::Input::Initialization:: *)
JordanWignerTensorProduct[Length_,StringOperator_,Operator_,Position_]:=Module[{list,chain},
list=Join[ConstantArray[StringOperator,Position-1],{Operator},ConstantArray[IdentityMatrix[2],Length-Position]];
chain=Apply[KroneckerProduct,list]]


End[];


(* ::Section:: *)
(*Functions to make Matrices for Ising Chain Hamiltonian and Operators*)


(* ::Input::Initialization:: *)
HIsingIntMode[L_,n_]:=ChainNeighbourTensorProduct[L,PauliMatrix[3],PauliMatrix[3],n];
HIsingBMode[L_,n_]:=ChainTensorProduct[L,PauliMatrix[1],n];


HIsingAltIntMode[L_,n_]:=ChainNeighbourTensorProduct[L,PauliMatrix[1],PauliMatrix[1],n];
HIsingAltBMode[L_,n_]:=ChainTensorProduct[L,-PauliMatrix[3],n];


HIsingAltIntMode[3,3]


KroneckerProduct[PauliMatrix[1],IdentityMatrix[2],PauliMatrix[1]]


(* ::Subsection:: *)
(*Local Operator Modes*)


(* ::Input::Initialization:: *)
\[Sigma]Plus0=(PauliMatrix[1]+I*PauliMatrix[2])/2;
\[Sigma]Minus0=(PauliMatrix[1]-I*PauliMatrix[2])/2;
\[Sigma]PlusAlt0=(-PauliMatrix[3]+I*PauliMatrix[2])/2;
\[Sigma]MinusAlt0=(-PauliMatrix[3]-I*PauliMatrix[2])/2;


\[Sigma]xMode[L_,n_] := ChainTensorProduct[L, PauliMatrix[1], n];
\[Sigma]yMode[L_,n_] := ChainTensorProduct[L, PauliMatrix[2], n];
\[Sigma]zMode[L_,n_] := ChainTensorProduct[L, PauliMatrix[3], n];
\[Sigma]PlusMode[L_,n_] := ChainTensorProduct[L, \[Sigma]Plus0, n];
\[Sigma]MinusMode[L_,n_] := ChainTensorProduct[L, \[Sigma]Minus0, n];
\[Sigma]PlusAltMode[L_,n_] := ChainTensorProduct[L, \[Sigma]PlusAlt0, n];
\[Sigma]MinusAltMode[L_,n_] := ChainTensorProduct[L, \[Sigma]MinusAlt0, n];


(* ::Subsection:: *)
(*Non-local Operator Modes*)


(* ::Input::Initialization:: *)
ccMode[L_,n_]:=JordanWignerTensorProduct[L,PauliMatrix[3],\[Sigma]Plus0,n];
cdMode[L_,n_]:=JordanWignerTensorProduct[L,PauliMatrix[3],\[Sigma]Minus0,n];


(* ::Input::Initialization:: *)
ccAltMode[L_,n_]:=JordanWignerTensorProduct[L,PauliMatrix[1],\[Sigma]PlusAlt0,n];
cdAltMode[L_,n_]:=JordanWignerTensorProduct[L,PauliMatrix[1],\[Sigma]MinusAlt0,n];


(* ::Subsection:: *)
(*Operators which are summed over the chain*)


TFIMHamiltonian[L_,IntStrength_,BStrength_,PBC_:False]:=If[PBC==True,-IntStrength*Sum[HIsingIntMode[L,n],{n,1,L}],-IntStrength*Sum[HIsingIntMode[L,n],{n,1,L-1}]] - BStrength*Sum[HIsingBMode[L,n],{n,1,L}];


TFIMAltHamiltonian[L_,IntStrength_,BStrength_,PBC_:False]:=If[PBC==True,-IntStrength*Sum[HIsingAltIntMode[L,n],{n,1,L}],-IntStrength*Sum[HIsingAltIntMode[L,n],{n,1,L-1}]] - BStrength*Sum[HIsingAltBMode[L,n],{n,1,L}];


\[Sigma]x[L_] := Sum[\[Sigma]xMode[L,n],{n,1,L}];


\[Sigma]y[L_] := Sum[\[Sigma]yMode[L,n],{n,1,L}];


\[Sigma]z[L_] := Sum[\[Sigma]zMode[L,n],{n,1,L}];


\[Sigma]Plus[L_] := Sum[\[Sigma]PlusMode[L,n],{n,1,L}]; 


\[Sigma]Minus[L_] := Sum[\[Sigma]MinusMode[L,n],{n,1,L}]; 


cc[L_] :=Sum[ccMode[L,n],{n,1,L}];


cd[L_] := Sum[cdMode[L,n], {n, 1, L}];


ccAlt[L_] :=Sum[ccAltMode[L,n],{n,1,L}];


cdAlt[L_] := Sum[cdAltMode[L,n], {n, 1, L}];


(* ::Section:: *)
(*END PRIVATE CLASS*)


EndPackage[];
