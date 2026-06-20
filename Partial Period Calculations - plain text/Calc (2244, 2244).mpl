read "ECHeader.mpl":
read "DEHeader.mpl":
read "PSHeader.mpl":
with(ECHeader);
with(DEHeader);
with(PSHeader);
"DOCUMENT START";
"Inputting data about the surfaces:";
Digits := 1000;
Y := 1;
g2_1 := `*`(12,X^4-X^2*Y^2+Y^4);
g3_1 := `*`(4,2*X^6-3*X^4*Y^2-3*X^2*Y^4+2*Y^6);
C1 := 4*x^3-g2_1*x-g3_1;
tempS := solve(1/J(C1) = 0):
tempS;
S := {tempS}:
g2_2 := g2_1;
g3_2 := g3_1;
C2 := C1;
tempS := solve(1/J(C2) = 0):
tempS;
S := S union {tempS}:
S := sort([seq(s,s in S)],(a, b) -> is(a < b)):
print('[A, B, C]' = S);
XA := S[1];
XB := S[2];
XC := S[3];
txOFF := [`/`(3,64), `/`(1,16), `/`(5,64)];
"Plots of zeroes between singularities:";
with(plots):
p[1] := S[1]:
p[2] := S[2]:
p[3] := S[3]:
p[0] := p[1]-`/`(1,2):
p[4] := p[3]+`/`(1,2):
N := 30:
q := ptmp:
R1 := [seq([seq([solve(0 = subs(X = q(i,k),C1))],k = 1 .. N-1)],i = 1 .. 4)]:
display(Array(1 .. 4,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 4)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := collect(GetDiffEq4(series(Period(C1)*Period(C2),X,100),100,6),D):
equ := sort(sort(collect(denom(equ)*equ,D),X,descending),D,descending);
sort([solve(lcoeff(equ,D) = 0)],(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around A for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,XA,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txA := `+`~(XA,`/`(txOFF,2));
PerProdA := ExplEvalAndComb3(C1,C2,txA,Solutions,1,2,0,[1, 1]):
IntA[1] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[1]),Z)):
IntA[2] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[2]),Z)):
IntA[3] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[3]),Z)):
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,XA):
# Identifying the coincidences at point -1:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResAB := CalcToPoint3(equ,IntA,[XA, -`/`(1,2), XB],350):
# The integral (between points -1 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -1 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -1 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],XB):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB,`/`(txOFF,4));
MatB := Matrix(3,3,ExplEvalAndComb3(C1,C2,txB,ResAB[3],0,1)):
evalf(MatB,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatB := MTM:-round(MatB);
with(LinearAlgebra):
PerProdB := MatrixVectorMultiply(MatB,ResAB[3]):
IntB := MatrixVectorMultiply(MatB,ResAB[2]-ResAB[1]):
"=======================================";
"===== Main calculations for the segment BC: =====";
"=======================================";
OutputCoincidences(PerProdB,XB):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResBC := CalcToPoint3(equ,IntB,[XB, `/`(1,2), XC],350):
# The integral (between points 0 and 1) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 1) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 1) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],XC):
# Identifying the coincidences at point 1:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Computing a transition matrix at C (for future reference -- checking solutions):";
txC := `+`~(XC,txOFF);
MatC := Matrix(3,3,ExplEvalAndComb3(C1,C2,txC,ResBC[3],0,1)):
evalf(MatC,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatC := MTM:-round(MatC);
