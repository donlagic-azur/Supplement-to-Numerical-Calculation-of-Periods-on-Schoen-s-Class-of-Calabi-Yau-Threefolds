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
F := -(X-1)*(2*X-1);
G := F/X;
C_Tate := y^2+((1+G)*x+F)*y = F*x^2+x^3:
C_Tate := subs(X = X+`/`(1,2),C_Tate):
N := 2;
#<---- power of X in denominator
K := 1;
#
g2 := simplify(subs(X = -1/X,G2(Wei(C_Tate))*`^`(sqrt(6)*(2*X+1),4))*(X^N)^4);
g3 := simplify(subs(X = -1/X,G3(Wei(C_Tate))*`^`(sqrt(6)*(2*X+1),6))*(X^N)^6);
C1 := 4*x^3-g2*x-g3;
C2 := C1:
tempS := solve(1/J(C1) = 0):
S := {tempS}:
J(C1);
Delta(C1);
S := sort([seq(s,s in S)],(x, y) -> is(Re(x) < Re(y))):
print(Vector(['D', E, 'F', 'A', B]) = evalf(Re(Vector(S)),10));
YD := S[1]:
YE := S[2]:
YF := S[3]:
YA := S[4]:
YB := S[5]:
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
equ := collect(GetDiffEq4(series(X^(2*N-2-K)*Period(C1)*Period(C2),X,100),100,11),D):
equ := sort(sort(collect(denom(equ)*equ,D),X,descending),D,descending);
sort(evalf(Re(Vector([solve(lcoeff(equ,D) = 0)])),12),(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around D for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,YD,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txD := Re(evalf(Vector(`+`~(YD,`/`(txOFF,2))))):
evalf(txD,10);
PerProdD := (-1)^K*ExplEvalAndComb3(C1,C2,txD,Solutions,1,2,2*N-2-K):
IntD[1] := subs(Z = X-Revalf(YD),int(subs(X = Z+Revalf(YD),PerProdD[1]),Z)):
IntD[2] := subs(Z = X-Revalf(YD),int(subs(X = Z+Revalf(YD),PerProdD[2]),Z)):
IntD[3] := subs(Z = X-Revalf(YD),int(subs(X = Z+Revalf(YD),PerProdD[3]),Z)):
"=======================================";
"===== Main calculations for the segment DE: =====";
"=======================================";
OutputCoincidences(PerProdD,Revalf(YD)):
# Identifying the coincidences at point -2.82843:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResDE := CalcToPoint3(equ,IntD,[YD, -`/`(12,5), YE],350):
# The integral (between points -2.82843 and -2) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -2.82843 and -2) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -2.82843 and -2) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDE[1]);
OutputCoincidences(ResDE[3],Revalf(YE)):
# Identifying the coincidences at point -2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around E for the four products of periods (by computing a transition matrix):";
txE := `+`~(YE,txOFF);
MatE := (-1)^K*Matrix(3,3,ExplEvalAndComb3(C1,C2,txE,ResDE[3],0,2,2*N-2-K)):
evalf(MatE,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatE := MTM:-round(MatE);
with(LinearAlgebra):
PerProdE := MatrixVectorMultiply(MatE,ResDE[3]):
IntE := MatrixVectorMultiply(MatE,ResDE[2]-ResDE[1]):
"=======================================";
"===== Main calculations for the segment EF: =====";
"=======================================";
OutputCoincidences(PerProdE,Revalf(YE)):
# Identifying the coincidences at point -2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResEF := CalcToPoint3(equ,IntE,[YE, -`/`(8,5), -`/`(6,5), -1, YF],350):
# The integral (between points -2 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -2 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -2 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResEF[1]);
OutputCoincidences(collect(1/(X^(2*N-2-K))*ResEF[3],X),Revalf(YF)):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around F for the four products of periods (by computing a transition matrix):";
txF := Re(evalf(Vector(`+`~(YF,`/`(txOFF,2))))):
evalf(txF,10);
MatF := (-1)^K*Matrix(3,3,ExplEvalAndComb3(C1,C2,txF,ResEF[3],0,2,2*N-2-K)):
evalf(MatF,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatF := MTM:-round(MatF);
with(LinearAlgebra):
PerProdF := MatrixVectorMultiply(MatF,ResEF[3]):
IntF := MatrixVectorMultiply(MatF,ResEF[2]-ResEF[1]):
"=======================================";
"===== Main calculations for the segment FA: =====";
"=======================================";
OutputCoincidences(collect(1/(X^(2*N-2-K))*PerProdF,X),Revalf(YF)):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResFA := CalcToPoint3(equ,IntF,[YF, 1, `/`(6,5), `/`(8,5), YA],350):
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResFA[1]);
OutputCoincidences(ResFA[3],Revalf(YA)):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Computing a transition matrix at A (for future reference -- checking solutions):";
txA := Re(evalf(Vector(`+`~(YA,`/`(txOFF,2))))):
evalf(txA,10);
MatA := (-1)^K*Matrix(3,3,ExplEvalAndComb3(C1,C2,txA,ResFA[3],0,2,2*N-2-K)):
evalf(MatA,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatA := MTM:-round(MatA);
Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Typesetting[delayDotProduct](Matrix(3,3,[[1,
    0, 0], [-1, 1, 0], [1, -2, 1]]),Matrix(3,3,[[1, 8, 16], [0, 1, 4], [0, 0, 1]]),true),Matrix(3,3,[[1, -2, 1],
    [0, 1, -1], [0, 0, 1]]),true),Matrix(3,3,[[9, 12, 4], [-6, -7, -2], [4, 4, 1]]),true),Matrix(3,3,[[0, 0, 1],
    [0, -1, 1], [1, -2, 1]]),true),Matrix(3,3,[[1, 8, 16], [0, 1, 4], [0, 0, 1]]),true);
evalf(ResDE[1],40);
evalf(ResEF[1],40);
evalf(ResFA[1],40);
