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
g2_2 := expand(X^4*subs(X = -1/X,simplify((2*X+1)^4*subs(X = -X/(2*X+1),g2))));
g3_2 := expand(X^6*subs(X = -1/X,simplify((2*X+1)^6*subs(X = -X/(2*X+1),g3))));
C2 := 4*x^3-g2_2*x-g3_2;
tempS := solve(1/J(C2) = 0):
tempS;
S := S union {tempS}:
S := sort([seq(s,s in S)]):
print('[D, E, A, B]' = S);
YD := S[1];
YE := S[2];
YA := S[3];
YB := S[4];
txOFF := [`/`(1,16), `/`(3,64), `/`(5,64), `/`(1,32)];
"Plots of zeroes between singularities:";
with(plots):
p[1] := S[1]:
p[2] := S[2]:
p[3] := S[3]:
p[4] := S[4]:
p[0] := p[1]-`/`(1,2):
p[5] := p[4]+`/`(1,2):
N := 30:
q := ptmp:
R1 := [seq([seq([solve(0 = subs(X = q(i,k),C1))],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R1[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R1[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
R2 := [seq([seq([solve(0 = subs(X = q(i,k),C2))],k = 1 .. N-1)],i = 1 .. 5)]:
display(Array(1 .. 5,[seq(pointplot([seq(seq([q(i,k),
    Re(R2[i][k][j])],k = 1 .. N-1),j = 1 .. 3)],color = ifelse(Im(R2[i][floor(.5*N)][2]) = 0,black,red)),i = 1 .. 5)]),size = [1000,
    1200]);
"Setting up the differential equation and checking its points of degeneracy:";
with(DETools):
equ := (1024*X^11-11264*X^10+52352*X^9-133248*X^8+198436*X^7-163324*X^6+47969*X^5+23283*X^4-12312*X^3-2916*X^2)*D^4+(15360*X^10-153600*X^9+657280*X^8-1571840*X^7+2254140*X^6-1841000*X^5+613585*X^4+138780*X^3-101250*X^2-14580*X)*D^3+(65536*X^9-589824*X^8+2309248*X^7-5154688*X^6+6989304*X^5-5419736*X^4+1809754*X^3+178002*X^2-175932*X-11664)*D^2+(86016*X^8-688128*X^7+2452864*X^6-5083392*X^5+6362120*X^4-4404256*X^3+1253558*X^2+85716*X-49248)*D+24576*X^7-172032*X^6+559104*X^5-1075200*X^4+1201920*X^3-681216*X^2+122112*X+20736;
sort(identify(evalf([solve(lcoeff(equ,D) = 0)],9)),(x, y) -> is(Re(x) < Re(y)));
"Getting expansions around D for the four products of periods (by solving the differential equation):";
SerSol := nformal_sol(equ,YD,500):
Solutions[1] := convert(SerSol[1],polynom):
Solutions[2] := convert(SerSol[2],polynom):
Solutions[3] := convert(SerSol[3],polynom):
Solutions[4] := convert(SerSol[4],polynom):
txD := `+`~(YD,txOFF);
PerProdD := ExplEvalAndComb(C1,C2,txD,Solutions,1,1):
IntD[1] := subs(Z = X-YD,int(subs(X = Z+YD,PerProdD[1]),Z)):
IntD[2] := subs(Z = X-YD,int(subs(X = Z+YD,PerProdD[2]),Z)):
IntD[3] := subs(Z = X-YD,int(subs(X = Z+YD,PerProdD[3]),Z)):
IntD[4] := subs(Z = X-YD,int(subs(X = Z+YD,PerProdD[4]),Z)):
"=======================================";
"===== Main calculations for the segment DE: =====";
"=======================================";
OutputCoincidences(PerProdD,YD):
# Identifying the coincidences at point -0.25:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# The 2nd curve doesn't have a singularity at this point.
ResDE := CalcToPoint4(equ,IntD,[YD, -`/`(1,8), YE],350):
# The integral (between points -0.25 and 0) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points -0.25 and 0) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResDE[1]);
OutputCoincidences(ResDE[3],YE):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Getting expansions around E for the four products of periods (by computing a transition matrix):";
txE := `+`~(YE,txOFF);
MatE := Matrix(4,4,ExplEvalAndComb(C1,C2,txE,ResDE[3],0,1)):
evalf(MatE,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatE := MTM:-round(MatE);
with(LinearAlgebra):
PerProdE := MatrixVectorMultiply(MatE,ResDE[3]):
IntE := MatrixVectorMultiply(MatE,ResDE[2]-ResDE[1]):
"=======================================";
"===== Main calculations for the segment EA: =====";
"=======================================";
OutputCoincidences(PerProdE,YE):
# Identifying the coincidences at point 0:
# Zeroes 1 and 2 of the 1st curve coincide at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
#ResEA&#8788; CalcToPoint4:
ResEA := CalcToPoint4(equ,IntE,[YE, `/`(1,8), `/`(1,4), `/`(3,8), `/`(3,4), 1+`/`(I,8), `/`(5,4), `/`(13,8),
    `/`(7,4), `/`(15,8), YA],350):
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 0 and 2) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResEA[1]);
OutputCoincidences(ResEA[3],YA):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
"Getting expansions around A for the four products of periods (by computing a transition matrix):";
txA := `+`~(YA,txOFF);
MatA := Matrix(4,4,ExplEvalAndComb(C1,C2,txA,ResEA[3],0,1)):
evalf(MatA,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatA := MTM:-round(MatA);
with(LinearAlgebra):
PerProdA := MatrixVectorMultiply(MatA,ResEA[3]):
IntA := MatrixVectorMultiply(MatA,ResEA[2]-ResEA[1]):
"=======================================";
"===== Main calculations for the segment AB: =====";
"=======================================";
OutputCoincidences(PerProdA,YA):
# Identifying the coincidences at point 2:
# Zeroes 2 and 3 of the 1st curve coincide at this point.
# Zeroes 1 and 2 of the 2nd curve coincide at this point.
ResAB := CalcToPoint4(equ,IntA,[YA, `/`(17,8), YB],350):
# The integral (between points 2 and 2.25) of the product of period 1, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 2 and 2.25) of the product of period 1, 1st curve, and period 2, 2nd curve, is:
# The integral (between points 2 and 2.25) of the product of period 2, 1st curve, and period 1, 2nd curve, is:
# The integral (between points 2 and 2.25) of the product of period 2, 1st curve, and period 2, 2nd curve, is:
ZRankRes(ResAB[1]);
OutputCoincidences(ResAB[3],YB):
# Identifying the coincidences at point 2.25:
# The 1st curve doesn't have a singularity at this point.
# Zeroes 2 and 3 of the 2nd curve coincide at this point.
"Computing a transition matrix at B (for future reference -- checking solutions):";
txB := `+`~(YB,txOFF);
MatB := Matrix(4,4,ExplEvalAndComb(C1,C2,txB,ResAB[3],0,1,0,`+`~(-sqrt(I),[0, 0, 0, 0]),[0, 1])):
evalf(MatB,5);
#(the above matrix needs to be integral, otherwise there has been an error in the explicit numerical calculation of integrals over some complex branches)
MatB := MTM:-round(MatB);
#checking if the product (from E to A) of the transition matrices we got is Identity - this reflects the correctness of the solutions, as the matrices were not only instrumental to the calculation of the results, but were also calculated using the information (analytically continued solutions) available at that moment
(((Matrix(4,4,[[1, 1, 0, 0], [0, 1, 0, 0], [-3, -3, 1, 1], [0, -3, 0, 1]]) . Matrix(4,4,[[1, 0, 0, 0], [0, 1,
    0, 0], [1, 0, -1, 0], [0, 1, 0, -1]])) . Matrix(4,4,[[-1, -1, 1, 1], [1, 2, -1, -2], [2, 2, -1, -1], [-2, -4,
    1, 2]])) . Matrix(4,4,[[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 0]])) . Matrix(4,4,[[1, 0, 1, 0],
    [-3, 1, -3, 1], [0, 0, 1, 0], [0, 0, -3, 1]]);
