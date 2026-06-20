read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
ResAB := EnDim(Vector(Vector[column](3,[-3.89485612473563672469880861636022750152128988,
    1.82133932962127051276404395566239914905074455*I, .994980355698474189519572883670708104639781313]))):
ResBC := EnDim(Vector(Vector[column](3,[20.1380946081475774749039948306567479222277000-35.5234108355035296015557697651552252063379483*I,
    -46.4177387781349775096153856149643639397484692+.2479945950e-44*I,
    20.1380946081475774749039948306567479222277000+35.5234108355035296015557697651552252063379483*I]))):
ResCD := EnDim(Vector(Vector[column](3,[78.4181744582046772249705207557830070934627783+.4350266782e-44*I,
    -.8535267925e-44-26.1756189776264677755574343630350759519811234*I,
    -11.1329363262649204719892964575080220711912960+.1262152606e-44*I]))):
ResDE := EnDim(Vector(Vector[column](3,[93.9121288870246399683617367289291265185923625-43.1981240407631357976211210600951834966699176*I,
    -108.331136841020932313822207007231157921191997,
    93.9121288870246399683617367289291265185923625+43.1981240407631357976211210600951834966699176*I]))):
ResEF := EnDim(Vector(Vector[column](3,[4.25040152016606085740823753852265295536676889,
    -.1307175381e-44-2.01601727300853258526863169180757999611524014*I,
    -1.11543257995285391148455495158346211163537917]))):
ResFA := EnDim(Vector(Vector[column](3,[-12.7697265404985156976772248307988897371380173+12.7957925737323407790218355772997667053938184*I,
    19.7082130129243074959862503171462089226977394-.1796245976e-43*I,
    -12.7697265404985156976772248307988897371380173-12.7957925737323407790218355772997667053938183*I]))):
MatA := EnDim(Matrix(3,3,[[1, -2, 1], [0, 1, -1], [0, 0, 1]]));
MatB := EnDim(Matrix(3,3,[[1, 8, 16], [-1, -7, -12], [1, 6, 9]]));
MatC := EnDim(Matrix(3,3,[[1, -2, 1], [0, 1, -1], [0, 0, 1]]));
MatD := EnDim(Matrix(3,3,[[1, 8, 16], [-1, -7, -12], [1, 6, 9]]));
MatE := EnDim(Matrix(3,3,[[1, -2, 1], [0, 1, -1], [0, 0, 1]]));
MatF := EnDim(Matrix(3,3,[[1, 8, 16], [-1, -7, -12], [1, 6, 9]]));
ResFA+MatF . ResEF+(MatF . MatE) . ResDE+((MatF . MatE) . MatD) . ResCD+(((MatF . MatE) . MatD) . MatC) . ResBC+((((MatF . MatE) . MatD) . MatC) . MatB) . ResAB;
MonodA := EnDim(Matrix(3,3,[[1, 0, 0], [-1, 1, 0], [1, -2, 1]]));
MonodB := EnDim(Matrix(3,3,[[1, 14, 49], [0, 1, 7], [0, 0, 1]]));
MonodC := EnDim(Matrix(3,3,[[225, 1470, 2401], [-60, -391, -637], [16, 104, 169]]));
MonodD := EnDim(Matrix(3,3,[[484, 2772, 3969], [-154, -881, -1260], [49, 280, 400]]));
MonodE := EnDim(Matrix(3,3,[[484, 2156, 2401], [-198, -881, -980], [81, 360, 400]]));
MonodF := EnDim(Matrix(3,3,[[225, 840, 784], [-105, -391, -364], [49, 182, 169]]));
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
