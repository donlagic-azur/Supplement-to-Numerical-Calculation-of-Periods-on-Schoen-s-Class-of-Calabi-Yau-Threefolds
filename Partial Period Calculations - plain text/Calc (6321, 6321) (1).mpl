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
g2 := `*`(12,X^4-4*X^3*Y+2*X*Y^3+Y^4);
g3 := `*`(4,2*X^6-12*X^5*Y+12*X^4*Y^2+14*X^3*Y^3+3*X^2*Y^4+6*X*Y^5+2*Y^6);
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
S := sort([seq(s,s in S)]):
print('[C, D, A]' = S);
YC := S[1];
YD := S[2];
YA := S[3];
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
equ := (16*X^6-56*X^5+33*X^4+28*X^3+4*X^2)*D^3+(144*X^5-420*X^4+198*X^3+126*X^2+12*X)*D^2+(304*X^4-704*X^3+222*X^2+100*X+4)*D+128*X^3-216*X^2+24*X+8;
sort([solve(lcoeff(equ,D) = 0)],(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around C for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,YC,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
txC := `+`~(YC,txOFF);
PerProdC := ExplEvalAndComb3(C1,C2,txC,Solutions,1,1):
IntC[1] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[1]),Z)):
IntC[2] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[2]),Z)):
IntC[3] := subs(Z = X-YC,int(subs(X = Z+YC,PerProdC[3]),Z)):
"=======================================";
"===== Main calculations for the segment CD: =====";
"=======================================";
CircC := Matrix(3,3,GetMonodromy3(`-`~(2*YC,txC),PerProdC)):
evalf(CircC,5);
#monodromy matrix at C &#8722;&#8722; should be integral
CircC := MTM:-round(CircC);
OutputCoincidences(PerProdC,YC):
# Identifying the coincidences at point -0.25:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
ResCD := CalcToPoint3(equ,IntC,[YC, -`/`(1,8), YD],350):
# The integral (between points -0.25 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResCD[1]);
OutputCoincidences(ResCD[3],YD):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around D for the four products of periods (by computing a transition matrix):";
txD := `+`~(YD,txOFF);
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
CircD := Matrix(3,3,GetMonodromy3(`-`~(2*YD,txD),PerProdD)):
evalf(CircD,5);
#monodromy matrix at D &#8722;&#8722; should be integral
CircD := MTM:-round(CircD);
OutputCoincidences(PerProdD,YD):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResDA := CalcToPoint3(equ,IntD,[YD, `/`(1,8), `/`(1,4), `/`(3,8), `/`(3,4), 1, `/`(5,4), `/`(13,8), `/`(7,4),
    `/`(15,8), YA],350):
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDA[1]);
OutputCoincidences(ResDA[3],YA):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Computing a transition matrix at A (for future reference -- checking solutions):";
txA := `+`~(YA,txOFF);
MatA := Matrix(3,3,ExplEvalAndComb3(C1,C2,txA,ResDA[3],0,1)):
evalf(MatA,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatA := MTM:-round(MatA);
#checking if the product (from D to A) of the transition matrices we got is Identity &#8722; this reflects the correctness of the solutions, as the matrices were not only instrumental to the calculation of the results, but were also calculated using the information availableatthatmoment
((Matrix(3,3,[[1, 0, 0], [-3, 1, 0], [9, -6, 1]]) . Matrix(3,3,[[1, 0, 0], [1, -1, 0], [1, -2,
    1]])) . Matrix(3,3,[[1, -2, 1], [-2, 3, -1], [4, -4, 1]])) . Matrix(3,3,[[1, 2, 1], [0, 1, 1], [0, 0, 1]]);
`-`~(2*XC,txC);
