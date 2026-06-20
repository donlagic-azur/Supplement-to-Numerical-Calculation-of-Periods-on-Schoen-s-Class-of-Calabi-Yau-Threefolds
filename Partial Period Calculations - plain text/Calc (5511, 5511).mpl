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
g2 := `*`(3,X^4-12*X^3*Y+14*X^2*Y^2+12*X*Y^3+Y^4);
g3 := X^6-18*X^5*Y+75*X^4*Y^2+75*X^2*Y^4+18*X*Y^5+Y^6;
g2_1 := expand(X^4*subs(X = -1/X,g2));
g3_1 := expand(X^6*subs(X = -1/X,g3));
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
print('[C, D, A]' = S);
YC := S[1];
YD := S[2];
YA := S[3];
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
equ := (X^6-22*X^5+119*X^4+22*X^3+X^2)*D^3+(9*X^5-165*X^4+714*X^3+99*X^2+3*X)*D^2+(19*X^4-276*X^3+842*X^2+78*X+1)*D+8*X^3-84*X^2+128*X+6;
sort([solve(lcoeff(equ,D) = 0)],(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around C for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,YC,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txC := `+`~(YC,`/`(txOFF,2));
PerProdC := ExplEvalAndComb3(C1,C2,txC,Solutions,1,1):
IntC[1] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[1]),Z)):
IntC[2] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[2]),Z)):
IntC[3] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[3]),Z)):
"=======================================";
"===== Main calculations for the segment CD: =====";
"=======================================";
OutputCoincidences(PerProdC,YC):
# Identifying the coincidences at point -0.0901699:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResCD := CalcToPoint3(equ,IntC,[YC, -`/`(1,20), YD],350):
# The integral (between points -0.0901699 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.0901699 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.0901699 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],YD):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around D for the four products of periods (by computing a transition matrix):";
txD := `+`~(YD,`/`(txOFF,4));
MatD := Matrix(3,3,ExplEvalAndComb3(C1,C2,txD,ResCD[3],0,1)):
evalf(MatD,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatD := MTM:-round(MatD);
with(LinearAlgebra):
PerProdD := MatrixVectorMultiply(MatD,ResCD[3]):
IntD := MatrixVectorMultiply(MatD,ResCD[2]-ResCD[1]):
"=======================================";
"===== Main calculations for the segment DA: =====";
"=======================================";
OutputCoincidences(PerProdD,YD):
# Identifying the coincidences at point 0:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
ResDA := CalcToPoint3(equ,IntD,[YD, `/`(1,20), `/`(3,20), `/`(3,20), `/`(7,20), `/`(7,20), `/`(3,4), `/`(3,4),
    `/`(10,7), `/`(10,7), `/`(5,2), `/`(5,2), 4, 4, 7, 7, YA],350):
# The integral (between points 0 and 11.0902) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 11.0902) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 11.0902) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDA[1]);
OutputCoincidences(ResDA[3],YA):
# Identifying the coincidences at point 11.0902:
# Zeroes 3 and 1 of the 1st curve coincide at this point.
# Zeroes 3 and 1 of the 2nd curve coincide at this point.
"Computing a transition matrix at A (for future reference -- checking solutions):";
txA := `+`~(YA,txOFF);
MatA := Matrix(3,3,ExplEvalAndComb3(C1,C2,txA,ResDA[3],0,1)):
evalf(MatA,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatA := MTM:-round(MatA);
#checking if the product (from D to A) of the transition matrices we got is Identity &#8722; this reflects the correctness of the solutions, as the matrices were not only instrumental to the calculation of the results, but were also calculated using the information availableatthatmoment
((Matrix(3,3,[[4, -4, 1], [2, -3, 1], [1, -2, 1]]) . Matrix(3,3,[[1, 0, 0], [3, -1, 0], [9, -6,
    1]])) . Matrix(3,3,[[4, -4, 1], [2, -3, 1], [1, -2, 1]])) . Matrix(3,3,[[1, 0, 0], [3, -1, 0], [9, -6, 1]]);
