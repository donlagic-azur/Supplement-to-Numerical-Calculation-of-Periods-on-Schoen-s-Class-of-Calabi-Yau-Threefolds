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
C_Tate := y^2+((1+G)*x+F)*y = F*x^2+x^3:
N := 3;
#<---- power of X in denominator
K := 4;
#
g2 := simplify(subs(X = -1/X,G2(Wei(C_Tate))*`^`(sqrt(6)*(X^2-6*X+4),4))*(X^N)^4);
g3 := simplify(subs(X = -1/X,G3(Wei(C_Tate))*`^`(sqrt(6)*(X^2-6*X+4),6))*(X^N)^6);
C1 := 4*x^3-g2*x-g3;
C2 := C1:
tempS := solve(1/J(C1) = 0):
S := {tempS}:
J(C1);
Delta(C1);
S := sort([seq(s,s in S)],(x, y) -> is(Re(x) < Re(y))):
print(Vector([C, 'D', E, 'F', 'G', 'H', 'A']) = evalf(Re(Vector(S)),10));
YC := S[1]:
YD := S[2]:
YE := S[3]:
YF := S[4]:
YG := S[5]:
YH := S[6]:
YA := S[7]:
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
M := 30:
q := ptmp:
R1 := [seq([seq([fsolve(0 = subs(X = q(i,k),C1),complex)],k = 1 .. M-1)],i = 1 .. 8)]:
display(Array(1 .. 8,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. M-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*M)][2]) = 0,black,red)),i = 1 .. 8)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := collect(GetDiffEq4(series(X^(2*N-2-K)*Period(C1)*Period(C2),X,100),100,14),D):
equ := sort(sort(collect(denom(equ)*equ,D),X,descending),D,descending);
sort(evalf(Re(Vector([solve(lcoeff(equ,D) = 0)])),12),(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around F for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,YF,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txF := Re(evalf(Vector(`+`~(YF,`/`(txOFF,2))))):
evalf(txF,10);
PerProdF := (-1)^K*ExplEvalAndComb3(C1,C2,txF,Solutions,1,2,2*N-2-K):
IntF[1] := subs(Z = X-Revalf(YF),int(subs(X = Z+Revalf(YF),PerProdF[1]),Z)):
IntF[2] := subs(Z = X-Revalf(YF),int(subs(X = Z+Revalf(YF),PerProdF[2]),Z)):
IntF[3] := subs(Z = X-Revalf(YF),int(subs(X = Z+Revalf(YF),PerProdF[3]),Z)):
"=======================================";
"===== Main calculations for the segment FG: =====";
"=======================================";
OutputCoincidences(PerProdF,Revalf(YF)):
# Identifying the coincidences at point -0.5:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResFG := CalcToPoint3(equ,IntF,[YF, -`/`(4,9), -`/`(4,11), -`/`(3,10), YG],350):
# The integral (between points -0.5 and -0.190983) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5 and -0.190983) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.5 and -0.190983) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResFG[1]);
OutputCoincidences(ResFG[3],YG):
# Identifying the coincidences at point -0.190983:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around G for the four products of periods (by computing a transition matrix):";
txG := `+`~(YG,txOFF);
MatG := (-1)^K*Matrix(3,3,ExplEvalAndComb3(C1,C2,txG,ResFG[3],0,2,2*N-2-K)):
evalf(MatG,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatG := MTM:-round(MatG);
with(LinearAlgebra):
PerProdG := MatrixVectorMultiply(MatG,ResFG[3]):
IntG := MatrixVectorMultiply(MatG,ResFG[2]-ResFG[1]):
"=======================================";
"===== Main calculations for the segment GH: =====";
"=======================================";
OutputCoincidences(PerProdG,YG):
# Identifying the coincidences at point -0.190983:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResGH := CalcToPoint3(equ,IntG,[YG, -`/`(1,10), YH],350):
# The integral (between points -0.190983 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.190983 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.190983 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResGH[1]);
OutputCoincidences(collect(1/(X^(2*N-2-K))*ResGH[3],X),Revalf(YH)):
# Identifying the coincidences at point 0:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around H for the four products of periods (by computing a transition matrix):";
txH := Re(evalf(Vector(`+`~(YH,`/`(txOFF,2))))):
evalf(txH,10);
MatH := (-1)^K*Matrix(3,3,ExplEvalAndComb3(C1,C2,txH,ResGH[3],0,2,2*N-2-K)):
evalf(MatH,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatH := MTM:-round(MatH);
with(LinearAlgebra):
PerProdH := MatrixVectorMultiply(MatH,ResGH[3]):
IntH := MatrixVectorMultiply(MatH,ResGH[2]-ResGH[1]):
"=======================================";
"===== Main calculations for the segment HA: =====";
"=======================================";
OutputCoincidences(collect(1/(X^(2*N-2-K))*PerProdH,X),Revalf(YH)):
# Identifying the coincidences at point 0:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResHA := CalcToPoint3(equ,IntH,[YH, `/`(1,10), `/`(1,5), `/`(3,10), `/`(3,5), `/`(9,10), YA],350):
# The integral (between points 0 and 1.61803) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 1.61803) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 1.61803) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResHA[1]);
OutputCoincidences(ResHA[3],Revalf(YA)):
# Identifying the coincidences at point 1.61803:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Computing a transition matrix at A (for future reference -- checking solutions):";
txA := Re(evalf(Vector(`+`~(YA,`/`(txOFF,2))))):
evalf(txA,10);
MatA := (-1)^K*Matrix(3,3,ExplEvalAndComb3(C1,C2,txA,ResHA[3],0,2,2*N-2-K)):
evalf(MatA,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatA := MTM:-round(MatA);
Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Matrix(3,3,[[1,
    6, 9], [-1, -5, -6], [1, 4, 4]]),Matrix(3,3,[[1, 0, 0], [-1, 1, 0], [1, -2, 1]]),true),Matrix(3,3,[[1, 10,
    25], [0, 1, 5], [0, 0, 1]]),true),Matrix(3,3,[[1, -2, 1], [0, 1, -1], [0, 0, 1]]),true),Matrix(3,3,[[1, 6, 9],
    [-1, -5, -6], [1, 4, 4]]),true),Matrix(3,3,[[1, 0, 0], [-1, 1, 0], [1, -2, 1]]),true),Matrix(3,3,[[1, 10, 25],
    [0, 1, 5], [0, 0, 1]]),true),Matrix(3,3,[[1, -2, 1], [0, 1, -1], [0, 0, 1]]),true);
evalf(ResFG[1],40);
evalf(ResGH[1],40);
evalf(ResHA[1],40);
