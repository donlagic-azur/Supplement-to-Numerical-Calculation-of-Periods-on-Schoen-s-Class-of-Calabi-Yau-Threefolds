read "ECHeader.mpl":
read "DEHeader.mpl":
read "PSHeader.mpl":
with(ECHeader);
with(DEHeader);
with(PSHeader);
"DOCUMENT START";
"Inputting data about the surfaces:";
Digits := 1000;
interface(rtablesize = 20);
Y := 1;
F := -X^3+X^2;
G := -X^2+X;
C_Tate := y^2+((1+G)*x+F)*y = F*x^2+x^3;
K := 1;
#
C1 := sort(sort(collect(expand(Wei(C_Tate)),x),X,descending),x,descending);
C2 := C1:
tempS := solve(1/J(C1) = 0):
S := {tempS}:
J(C1);
Delta(C1);
S := sort([seq(s,s in S)],(x, y) -> is(Re(x) < Re(y))):
print(Vector([A, B, C, D, E]) = evalf(Re(Vector(S)),10));
XA := S[1]:
XB := S[2]:
XC := S[3]:
XD := S[4]:
XE := S[5]:
txOFF := [`/`(1,16), `/`(3,64), `/`(5,64)];
"Plots of zeroes between singularities:";
with(plots):
p[1] := Re(S[1]):
p[2] := Re(S[2]):
p[3] := Re(S[3]):
p[4] := Re(S[4]):
p[5] := Re(S[5]):
p[0] := p[1]-`/`(1,2):
p[6] := p[5]+`/`(1,2):
M := 30:
q := ptmp:
R1 := [seq([seq([fsolve(0 = subs(X = q(i,k),C1),complex)],k = 1 .. M-1)],i = 1 .. 6)]:
display(Array(1 .. 6,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. M-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*M)][2]) = 0,black,red)),i = 1 .. 6)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := collect(GetDiffEq4(series(X^K*Period(C1)*Period(C2),X,100),100,11),D):
equ := sort(sort(collect(denom(equ)*equ,D),X,descending),D,descending);
sort(evalf(Re(Vector([solve(lcoeff(equ,D) = 0)])),12),(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around A for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,XA,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txA := Re(evalf(Vector(`+`~(XA,`/`(txOFF,2))))):
evalf(txA,10);
PerProdA := ExplEvalAndComb3(C1,C2,txA,Solutions,1,2,K):
IntA[1] := subs(Z = X-Revalf(XA),int(subs(X = Z+Revalf(XA),PerProdA[1]),Z)):
IntA[2] := subs(Z = X-Revalf(XA),int(subs(X = Z+Revalf(XA),PerProdA[2]),Z)):
IntA[3] := subs(Z = X-Revalf(XA),int(subs(X = Z+Revalf(XA),PerProdA[3]),Z)):
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,Revalf(XA)):
# Identifying the coincidences at point -0.158834:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResAB := CalcToPoint3(equ,IntA,[XA, -`/`(1,12), XB],350):
# The integral (between points -0.158834 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.158834 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.158834 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(collect(1/(X^K)*ResAB[3],X),Revalf(XB)):
# Identifying the coincidences at point 0:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB,txOFF);
MatB := Matrix(3,3,ExplEvalAndComb3(C1,C2,txB,ResAB[3],0,2,K)):
evalf(MatB,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatB := MTM:-round(MatB);
with(LinearAlgebra):
PerProdB := MatrixVectorMultiply(MatB,ResAB[3]):
IntB := MatrixVectorMultiply(MatB,ResAB[2]-ResAB[1]):
"=======================================";
"===== Main calculations for the segment BC: =====";
"=======================================";
OutputCoincidences(collect(1/(X^K)*PerProdB,X),Revalf(XB)):
# Identifying the coincidences at point 0:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResBC := CalcToPoint3(equ,IntB,[XB, `/`(1,12), `/`(1,6), `/`(3,12), `/`(1,2), `/`(9,14), `/`(5,7), `/`(11,14),
    XC],350):
# The integral (between points 0 and 0.862937) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 0.862937) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 0.862937) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],Revalf(XC)):
# Identifying the coincidences at point 0.862937:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Getting expansions around C for the four products of periods (by computing a transition matrix):";
txC := Re(evalf(Vector(`+`~(XC,`/`(txOFF,2))))):
evalf(txC,10);
MatC := Matrix(3,3,ExplEvalAndComb3(C1,C2,txC,ResBC[3],0,2,K)):
evalf(MatC,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatC := MTM:-round(MatC);
with(LinearAlgebra):
PerProdC := MatrixVectorMultiply(MatC,ResBC[3]):
IntC := MatrixVectorMultiply(MatC,ResBC[2]-ResBC[1]):
"=======================================";
"===== Main calculations for the segment CD: =====";
"=======================================";
OutputCoincidences(PerProdC,Revalf(XC)):
# Identifying the coincidences at point 0.862937:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResCD := CalcToPoint3(equ,IntC,[XC, `/`(13,14), XD],350):
# The integral (between points 0.862937 and 1) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0.862937 and 1) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0.862937 and 1) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],Revalf(XD)):
# Identifying the coincidences at point 1:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around D for the four products of periods (by computing a transition matrix):";
txD := Re(evalf(Vector(`+`~(XD,`/`(txOFF,2))))):
evalf(txD,10);
MatD := Matrix(3,3,ExplEvalAndComb3(C1,C2,txD,ResCD[3],0,2,K,[1, 1])):
evalf(MatD,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatD := MTM:-round(MatD);
with(LinearAlgebra):
PerProdD := MatrixVectorMultiply(MatD,ResCD[3]):
IntD := MatrixVectorMultiply(MatD,ResCD[2]-ResCD[1]):
"=======================================";
"===== Main calculations for the segment DE: =====";
"=======================================";
OutputCoincidences(PerProdD,Revalf(XD)):
# Identifying the coincidences at point 1:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResDE := CalcToPoint3(equ,IntD,[XD, `/`(15,14), `/`(8,7), `/`(17,14), `/`(10,7), `/`(23,14), `/`(16,7),
    `/`(22,7), XE],350):
# The integral (between points 1 and 7.2959) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 1 and 7.2959) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 1 and 7.2959) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDE[1]);
OutputCoincidences(ResDE[3],Revalf(XE)):
# Identifying the coincidences at point 7.2959:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Getting expansions around E for the four products of periods (by computing a transition matrix):";
txE := Re(evalf(Vector(`+`~(XE,`/`(txOFF,2))))):
evalf(txE,10);
MatE := Matrix(3,3,ExplEvalAndComb3(C1,C2,txE,ResDE[3],0,2,K)):
evalf(MatE,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatE := MTM:-round(MatE);
evalf(ResAB[1],40);
evalf(ResBC[1],40);
evalf(ResCD[1],40);
evalf(ResDE[1],40);
