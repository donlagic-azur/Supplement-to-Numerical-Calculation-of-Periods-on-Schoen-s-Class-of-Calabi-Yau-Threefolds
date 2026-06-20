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
g2_1 := `*`(12,X^4-4*X^3*Y+2*X*Y^3+Y^4);
g3_1 := `*`(4,2*X^6-12*X^5*Y+12*X^4*Y^2+14*X^3*Y^3+3*X^2*Y^4+6*X*Y^5+2*Y^6);
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
S := sort([seq(s,s in S)]):
print('[A, B, C]' = S);
XA := S[1];
XB := S[2];
XC := S[3];
txOFF := [`/`(1,16), `/`(3,64), `/`(5,64)];
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
equ := (4*X^6-28*X^5+33*X^4+56*X^3+16*X^2)*D^3+(36*X^5-210*X^4+198*X^3+252*X^2+48*X)*D^2+(76*X^4-352*X^3+222*X^2+200*X+16)*D+32*X^3-108*X^2+24*X+16;
sort(identify(evalf([solve(lcoeff(equ,D) = 0)],12)),(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around A for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,XA,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txA := `+`~(XA,txOFF);
PerProdA := ExplEvalAndComb3(C1,C2,txA,Solutions,1,1):
IntA[1] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[1]),Z)):
IntA[2] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[2]),Z)):
IntA[3] := subs(Z = X-XA,int(subs(X = Z+XA,PerProdA[3]),Z)):
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,XA):
# Identifying the coincidences at point -0.5:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResAB := CalcToPoint3(equ,IntA,[XA, -`/`(1,4), XB],350):
# The integral (between points -0.5 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.5 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],XB):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around B for the four products of periods (by computing a transition matrix):";
txB := `+`~(XB,txOFF);
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
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResBC := CalcToPoint3(equ,IntB,[XB, `/`(1,4), `/`(1,2), `/`(3,4), `/`(3,2), 2, `/`(5,2), `/`(13,4), `/`(7,2),
    `/`(15,4), XC],350):
# The integral (between points 0 and 4) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 4) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 4) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResBC[1]);
OutputCoincidences(ResBC[3],XC):
# Identifying the coincidences at point 4:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Computing a transition matrix at C (for future reference -- checking solutions):";
txC := `+`~(XC,txOFF);
MatC := Matrix(3,3,ExplEvalAndComb3(C1,C2,txC,ResBC[3],0,1)):
evalf(MatC,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatC := MTM:-round(MatC);
