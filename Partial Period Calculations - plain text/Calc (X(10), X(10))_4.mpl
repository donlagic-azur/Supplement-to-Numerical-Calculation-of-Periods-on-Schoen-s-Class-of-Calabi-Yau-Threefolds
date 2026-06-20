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
F := -X^3*(X^2-3*X+2)/(X^2-6*X+4)^2;
G := X*(X^2-3*X+2)/(X^2-6*X+4);
C_Tate := y^2+((1+G)*x+F)*y = F*x^2+x^3;
K := 4;
#
g2 := G2(Wei(C_Tate))*`^`(sqrt(6)*(X^2-6*X+4),4);
g3 := G3(Wei(C_Tate))*`^`(sqrt(6)*(X^2-6*X+4),6);
C1 := 4*x^3-g2*x-g3;
C2 := C1:
tempS := solve(1/J(C1) = 0):
S := {tempS}:
J(C1);
Delta(C1);
S := sort([seq(s,s in S)],(x, y) -> is(Re(x) < Re(y))):
print(Vector([A, B, C, 'D', E, 'F', 'G']) = evalf(Re(Vector(S)),10));
XA := S[1]:
XB := S[2]:
XC := S[3]:
XD := S[4]:
XE := S[5]:
XF := S[6]:
XG := S[7]:
txOFF := [`/`(1,16), `/`(3,64), `/`(5,64)];
"Plots of zeroes between singularities:";
with(plots):
p[1] := Re(S[1]):
p[2] := Re(S[2]):
p[3] := Re(S[3]):
p[4] := Re(S[4]):
p[5] := Re(S[5]):
p[6] := Re(S[6]):
p[7] := Re(S[7]):
p[0] := p[1]-`/`(1,2):
p[8] := p[7]+`/`(1,2):
N := 30:
q := ptmp:
R1 := [seq([seq([fsolve(0 = subs(X = q(i,k),C1),complex)],k = 1 .. N-1)],i = 1 .. 8)]:
display(Array(1 .. 8,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 8)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := collect(GetDiffEq4(series(X^K*Period(C1)*Period(C2),X,100),100,15),D):
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
# Identifying the coincidences at point -0.618034:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResAB := CalcToPoint3(equ,IntA,[XA, -`/`(1,3), XB],350):
# The integral (between points -0.618034 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.618034 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.618034 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(collect(1/(X^K)*ResAB[3],X),XB):
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
OutputCoincidences(collect(1/(X^K)*PerProdB,X),XB):
# Identifying the coincidences at point 0:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResBC := CalcToPoint3(equ,IntB,[XB, `/`(3,8), `/`(1,2), `/`(5,8), XC],350):
# The integral (between points 0 and 0.763932) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 0.763932) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 0.763932) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],Revalf(XC)):
# Identifying the coincidences at point 0.763932:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
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
# Identifying the coincidences at point 0.763932:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResCD := CalcToPoint3(equ,IntC,[XC, `/`(7,8), XD],350):
# The integral (between points 0.763932 and 1) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0.763932 and 1) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0.763932 and 1) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],Revalf(XD)):
# Identifying the coincidences at point 1:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around D for the four products of periods (by computing a transition matrix):";
txD := Re(evalf(Vector(`+`~(XD,`/`(txOFF,2))))):
evalf(txD,10);
MatD := Matrix(3,3,ExplEvalAndComb3(C1,C2,txD,ResCD[3],0,2,K)):
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
ResDE := CalcToPoint3(equ,IntD,[XD, `/`(17,15), `/`(9,7), `/`(7,5), XE],350):
# The integral (between points 1 and 1.61803) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 1 and 1.61803) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 1 and 1.61803) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDE[1]);
OutputCoincidences(ResDE[3],Revalf(XE)):
# Identifying the coincidences at point 1.61803:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Getting expansions around E for the four products of periods (by computing a transition matrix):";
txE := Re(evalf(Vector(`+`~(XE,`/`(txOFF,2))))):
evalf(txE,10);
MatE := Matrix(3,3,ExplEvalAndComb3(C1,C2,txE,ResDE[3],0,2,K)):
evalf(MatE,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatE := MTM:-round(MatE);
with(LinearAlgebra):
PerProdE := MatrixVectorMultiply(MatE,ResDE[3]):
IntE := MatrixVectorMultiply(MatE,ResDE[2]-ResDE[1]):
"=======================================";
"===== Main calculations for the segment EF: =====";
"=======================================";
OutputCoincidences(PerProdE,Revalf(XE)):
# Identifying the coincidences at point 1.61803:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResEF := CalcToPoint3(equ,IntE,[XE, `/`(9,5), XF],350):
# The integral (between points 1.61803 and 2) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 1.61803 and 2) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 1.61803 and 2) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResEF[1]);
OutputCoincidences(ResEF[3],Revalf(XF)):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around F for the four products of periods (by computing a transition matrix):";
txF := Re(evalf(Vector(`+`~(XF,`/`(txOFF,2))))):
evalf(txF,10);
MatF := Matrix(3,3,ExplEvalAndComb3(C1,C2,txF,ResEF[3],0,2,K)):
evalf(MatF,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatF := MTM:-round(MatF);
with(LinearAlgebra):
PerProdF := MatrixVectorMultiply(MatF,ResEF[3]):
IntF := MatrixVectorMultiply(MatF,ResEF[2]-ResEF[1]):
"=======================================";
"===== Main calculations for the segment FG: =====";
"=======================================";
OutputCoincidences(PerProdF,Revalf(XF)):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResFG := CalcToPoint3(equ,IntF,[XF, `/`(11,5), `/`(12,5), `/`(13,5), `/`(16,5), `/`(19,5), XG],350):
# The integral (between points 2 and 5.23607) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 2 and 5.23607) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 2 and 5.23607) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResFG[1]);
OutputCoincidences(ResFG[3],Revalf(XG)):
# Identifying the coincidences at point 5.23607:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Computing a transition matrix at G (for future reference -- checking solutions):";
txG := Re(evalf(Vector(`+`~(XG,`/`(txOFF,2))))):
evalf(txG,10);
MatG := Matrix(3,3,ExplEvalAndComb3(C1,C2,txG,ResFG[3],0,2,K)):
evalf(MatG,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatG := MTM:-round(MatG);
evalf(ResAB[1],40);
evalf(ResBC[1],40);
evalf(ResCD[1],40);
evalf(ResDE[1],40);
evalf(ResEF[1],40);
evalf(ResFG[1],40);
