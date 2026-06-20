read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
ResAB := EnDim(Vector(Vector[column](3,[1.03860014114220952265814983848121766837556143,
    -.305566310563530204248717900406070003960883489*I, -.116934128606181712674666514152865215891382444]))):
ResBC := EnDim(Vector(Vector[column](3,[.338697409736872703982343519367920847455676400e-1-.178996600523806853500416650664944373628019977*I,
    -.201071826106491977430534488119543067515933446,
    .338697409736872703982343519367920847455676400e-1+.178996600523806853500416650664944373628019977*I]))):
ResCD := EnDim(Vector(Vector[column](3,[.338697409736872703982343519367920847455676420e-1-.178996600523806853500416650664944373628019979*I,
    -.201071826106491977430534488119543067515933447+.1101677772e-44*I,
    .338697409736872703982343519367920847455676407e-1+.178996600523806853500416650664944373628019976*I]))):
ResDE := EnDim(Vector(Vector[column](3,[1.03860014114220952265814983848121766837556143,
    -.305566310563530204248717900406070003960883489*I, -.116934128606181712674666514152865215891382444]))):
ResEF := EnDim(Vector(Vector[column](3,[1.70007946345877704073737779731396334352648195,
    -.737702331166783759245737050553265638254630491*I, -.484134729977662574658065326704554545698945498]))):
ResFA := EnDim(Vector(Vector[column](3,[1.70007946345877704073737779731396334352648195,
    -.737702331166783759245737050553265638254630491*I, -.484134729977662574658065326704554545698945498]))):
MatA := EnDim(Matrix(3,3,[[1, 8, 16], [0, 1, 4], [0, 0, 1]]));
MatB := EnDim(Matrix(3,3,[[0, 0, 1], [0, -1, 1], [1, -2, 1]]));
MatC := EnDim(Matrix(3,3,[[9, 12, 4], [-6, -7, -2], [4, 4, 1]]));
MatD := EnDim(Matrix(3,3,[[1, -2, 1], [0, 1, -1], [0, 0, 1]]));
MatE := EnDim(Matrix(3,3,[[1, 8, 16], [0, 1, 4], [0, 0, 1]]));
MatF := EnDim(Matrix(3,3,[[1, 0, 0], [-1, 1, 0], [1, -2, 1]]));
ResFA+MatF . ResEF+(MatF . MatE) . ResDE+((MatF . MatE) . MatD) . ResCD+(((MatF . MatE) . MatD) . MatC) . ResBC+((((MatF . MatE) . MatD) . MatC) . MatB) . ResAB;
MonodA := EnDim(Matrix(3,3,[[1, 16, 64], [0, 1, 8], [0, 0, 1]]));
MonodB := EnDim(Matrix(3,3,[[1, 0, 0], [-1, 1, 0], [1, -2, 1]]));
MonodC := EnDim(Matrix(3,3,[[49, -224, 256], [28, -127, 144], [16, -72, 81]]));
MonodD := EnDim(Matrix(3,3,[[529, -2944, 4096], [207, -1151, 1600], [81, -450, 625]]));
MonodE := EnDim(Matrix(3,3,[[529, -3312, 5184], [184, -1151, 1800], [64, -400, 625]]));
MonodF := EnDim(Matrix(3,3,[[49, -448, 1024], [14, -127, 288], [4, -36, 81]]));
qAB := ResAB:
qBC := 1/MatB . ResBC:
qCD := (1/MatB . (1/MatC)) . ResCD:
qDE := ((1/MatB . (1/MatC)) . (1/MatD)) . ResDE:
qEF := (((1/MatB . (1/MatC)) . (1/MatD)) . (1/MatE)) . ResEF:
qFA := ((((1/MatB . (1/MatC)) . (1/MatD)) . (1/MatE)) . (1/MatF)) . ResFA:
wA := Vector([wtf_1, wtf_2, wtf_3, wtf_4]):
wB := wA+qAB:
wC := wB+qBC:
wD := wC+qCD:
wE := wD+qDE:
wF := wE+qEF:
((((`^`(MonodA,0) . `^`(MonodB,0)) . `^`(MonodC,0)) . `^`(MonodD,0)) . `^`(MonodE,0)) . `^`(MonodF,0);
qAB+qBC+qCD+qDE+qEF+qFA;
((((MonodA . MonodB) . MonodC) . MonodD) . MonodE) . MonodF;
qAB+MonodB . qBC+(MonodB . MonodC) . qCD+((MonodB . MonodC) . MonodD) . qDE+(((MonodB . MonodC) . MonodD) . MonodE) . qEF+((((MonodB . MonodC) . MonodD) . MonodE) . MonodF) . qFA;
NV := Array([]):
NV0 := Array([]):
NVnumb := Array([]):
Monod := [MonodA, MonodB, MonodC, MonodD, MonodE, MonodF]:
for M in Monod do
    # tM ~ NullSpace:
    tM := (M-1)^2:
    tM := [(M-1)[1], (M-1)[2], (M-1)[3], (M-1)[4], RedVec(tM[1]), RedVec(tM[2]), RedVec(tM[3]), RedVec(tM[4])]:
    tL := LLL(Matrix([seq([seq(IdentityMatrix(8)[tx,ty],ty = 1 .. 8),
        `^`(10,20)*seq(tM[tx][tz],tz = 1 .. 4)],tx = 1 .. 8)]));
    nv := Array([]):
    for s from 1 to 8 do
        if not (tL[s,8+1]=0 and tL[s,8+2]=0 and tL[s,8+3]=0 and tL[s,8+4]=0) then
            Append(nv,Vector([seq(`^`(10,-20)*tL[s,tx+8],tx = 1 .. 4)]));
        end if:
    end do:
    nv0 := nv;
    Append(NVnumb,numelems(nv));
    Append(NV0,nv);
    nv := add(d(k)[temps]*nv[temps],temps = 1 .. numelems(nv));
    Append(NV,nv);
    print(nv0,"--->",nv);
end do:
NW := [wA, wB, wC, wD, wE, wF]:
A := Array([]):
nv := [seq(subs(k = j,NV[j]),j = 1 .. 6)]:
v := nv[6]+nv[5]+nv[4]+nv[3]+nv[2]+nv[1]:
tD := [seq(seq(d(k)[j] = 0,j = 1 .. NVnumb[k]),k = 1 .. 6)]:
tn := numelems(tD):
tM := Matrix([seq(seq([seq(IdentityMatrix(tn)[add(NVnumb[k],k = 1 .. tx-1)+ty,tz],tz = 1 .. tn),
    `^`(10,20)*seq(subs(tD,subs(d(tx)[ty] = 1,v[tz])),tz = 1 .. 4)],ty = 1 .. NVnumb[tx]),tx = 1 .. 6)]):
tL := LLL(tM,integer):
print(tL);
for s from 1 to tn do
    if tL[s,tn+1]=0 and tL[s,tn+2]=0 and tL[s,tn+3]=0 and tL[s,tn+4]=0 then
        v := [seq(seq(d(tx)[ty] = tL[s,ty+add(NVnumb[j],j = 1 .. tx-1)],ty = 1 .. NVnumb[tx]),tx = 1 .. 6)];
        # print;
        ta := subs(v,nv[1] . NW[1]+nv[2] . NW[2]+nv[3] . NW[3]+nv[4] . NW[4]+nv[5] . NW[5]+nv[6] . NW[6]):
        if ta<>0 then
            print(ta):
            Append(A,ta):
        end if:
    else
        break;
    end if:
end do:
numelems(A);
VV := QBasisFast(A,35,1000000000000):
Eres1 := BasisLLL(convert(A,list),100,VV,35,1000000000000,`^`(10,20)):
for s from 1 to numelems(Eres1) do
    if abs(Re(Eres1[s]))<10^(-30) then
        Eres1[s] := I*Im(Eres1[s]);
    end if:
    if abs(Im(Eres1[s]))<10^(-30) then
        Eres1[s] := Re(Eres1[s]);
    end if:
end do:
evalf(Vector([Eres1]),30);
NV := Array([]):
NV0 := Array([]):
NVnumb := Array([]):
Monod := [MonodA, MonodB, MonodC, MonodD, MonodE, MonodF]:
for M in Monod do
    tM := (M-1)^(1+Rank((M-1)^2)):
    tL := LLL(Matrix([seq([seq(IdentityMatrix(4)[tx,ty],ty = 1 .. 4),
        `^`(10,20)*seq(tM[tx,tz],tz = 1 .. 4)],tx = 1 .. 4)])):
    nv := Array([]):
    for s from 1 to 4 do
        if tL[s,4+1]=0 and tL[s,4+2]=0 and tL[s,4+3]=0 and tL[s,4+4]=0 then
            Append(nv,Vector([seq(tL[s,tx],tx = 1 .. 4)]));
        end if:
    end do:
    nv0 := nv;
    Append(NVnumb,numelems(nv));
    Append(NV0,nv);
    nv := add(d(k)[temps]*nv[temps],temps = 1 .. numelems(nv));
    Append(NV,nv);
    print(nv0,"--->",nv);
end do:
NW := [wA, wB, wC, wD, wE, wF]:
A := Array([]):
nv := [seq(subs(k = j,NV[j]),j = 1 .. 6)]:
v := nv[6]+nv[5]+nv[4]+nv[3]+nv[2]+nv[1]:
tD := [seq(seq(d(k)[j] = 0,j = 1 .. NVnumb[k]),k = 1 .. 6)]:
tn := numelems(tD):
tM := Matrix([seq(seq([seq(IdentityMatrix(tn)[add(NVnumb[k],k = 1 .. tx-1)+ty,tz],tz = 1 .. tn),
    `^`(10,20)*seq(subs(tD,subs(d(tx)[ty] = 1,v[tz])),tz = 1 .. 4)],ty = 1 .. NVnumb[tx]),tx = 1 .. 6)]):
tL := LLL(tM,integer):
print(tL);
for s from 1 to tn do
    if tL[s,tn+1]=0 and tL[s,tn+2]=0 and tL[s,tn+3]=0 and tL[s,tn+4]=0 then
        v := [seq(seq(d(tx)[ty] = tL[s,ty+add(NVnumb[j],j = 1 .. tx-1)],ty = 1 .. NVnumb[tx]),tx = 1 .. 6)];
        # print;
        ta := subs(v,nv[1] . NW[1]+nv[2] . NW[2]+nv[3] . NW[3]+nv[4] . NW[4]+nv[5] . NW[5]+nv[6] . NW[6]):
        if ta<>0 then
            print(ta):
            Append(A,ta):
        end if:
    else
        break;
    end if:
end do:
numelems(A);
VV := QBasisFast(A,35,1000000000000):
Eres2 := BasisLLL(convert(A,list),100,VV,35,1000000000000,`^`(10,20)):
for s from 1 to numelems(Eres2) do
    if abs(Re(Eres2[s]))<10^(-30) then
        Eres2[s] := I*Im(Eres2[s]);
    end if:
    if abs(Im(Eres2[s]))<10^(-30) then
        Eres2[s] := Re(Eres2[s]);
    end if:
end do:
evalf(Vector([Eres2]),30);
coeff1 := 1.187241682451252893525893275023958089437:
tn1 := numelems(Eres1):
tn2 := numelems(Eres2):
rM := Vector(tn1):
for tt from 1 to tn1 do
    tM := Matrix([seq([seq(0,tz = 1 .. tn2), `^`(10,20)*(Re(Eres1[tt])+coeff1*Im(Eres1[tt]))],tx = 1 .. 1),
        seq([seq(IdentityMatrix(tn2)[tx,tz],tz = 1 .. tn2),
        seq(`^`(10,20)*(Re(Eres2[tx])+coeff1*Im(Eres2[tx])),tz = 1 .. 1)],tx = 1 .. tn2)]):
    # print;
    tM := LLL(tM);
    # print;
    rM[tt] := [seq(tM[1,tz],tz = 1 .. tn2)];
end do:
rM := Matrix(convert(rM,list));
# Vector, Vector
evalf(`*`~(rM . Vector(Eres2)-Vector(Eres1),rM . Vector(Eres2)+Vector(Eres1)),45);
round(Determinant(rM));
Vector([1, 2, 3]) *~ Vector([1, 2, 3]);
