read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
qAB := Vector([7.021010683705331299486227855323363057191+`*`(0,I),
    `+`(0,-`*`(9.167731910493613273137103164339106378094,I)),
    `+`(0,-`*`(8.409793987192955595486626337043473527717,I)),
    -14.40613737063808927966155725216831517849+`*`(0,I)]):
qBC := Vector([1.045223406789223777398452158191662475146+`*`(0,I),
    3.135670220367671332195356474574987425437-`*`(2.550618687964099305945383169915946892447,I),
    1.567835110183835666097678237287493712719-`*`(2.515172948817846991424615817358262517707,I),
    -1.953389349926397620504091601939292854923-`*`(11.37144687839968993319192220694870789179,I)]):
qCD := Vector([4.930563870126883744689323538940038106899+`*`(11.27187065606882523280180628170808142417,I),
    14.79169161038065123406797061682011432070+`*`(11.27187065606882523280180628170808142417,I),
    17.92736183074832256626332709139510174614+`*`(11.27187065606882523280180628170808142417,I),
    32.71905344112897380033129770821521606683+`*`(0,I)]):
qDE := Vector([-5.611671273688680841398674155610111517933+`*`(6.305655241617743635821923219674498481638,I),
    -3.521224460110233286601769839226786567642+`*`(8.856273929581842941767306389590445374085,I),
    -5.089059570294068952699448076514280280361+`*`(8.820828190435590627246539037032760999345,I),
    -1.953389349926397620504091601939292854923+`*`(11.37144687839968993319192220694870789179,I)]):
qEA := Vector([-7.385126686932757980175329396844952121302-`*`(17.57752589768656886862372950138257990581,I),
    -14.40613737063808927966155725216831517849-`*`(8.409793987192955595486626337043473527717,I),
    -14.40613737063808927966155725216831517849-`*`(9.167731910493613273137103164339106378094,I),
    -14.40613737063808927966155725216831517849+`*`(0,I)]):
MonodA := Matrix(4,4,[[1, 2, 2, 4], [0, 1, 0, 2], [0, 0, 1, 2], [0, 0, 0, 1]]);
MonodB := Matrix(4,4,[[1, 0, 0, 0], [-6, 1, 0, 0], [-3, 0, 1, 0], [18, -3, -6, 1]]);
MonodC := Matrix(4,4,[[-2, 1, 0, 0], [-9, 4, 0, 0], [0, 0, -2, 1], [0, 0, -9, 4]]);
MonodD := Matrix(4,4,[[-5, 0, 4, 0], [0, -5, 0, 4], [-9, 0, 7, 0], [0, -9, 0, 7]]);
MonodE := Matrix(4,4,[[10, -15, -12, 18], [15, -20, -18, 24], [12, -18, -14, 21], [18, -24, -21, 28]]);
wA := Vector([wtf_1, wtf_2, wtf_3, wtf_4]):
wB := wA+qAB:
wC := wB+qBC:
wD := wC+qCD:
wE := wD+qDE:
(((`^`(MonodA,0) . `^`(MonodB,0)) . `^`(MonodC,0)) . `^`(MonodD,0)) . `^`(MonodE,0);
qAB+qBC+qCD+qDE+qEA;
(((MonodA . MonodB) . MonodC) . MonodD) . MonodE;
qAB+MonodB . qBC+(MonodB . MonodC) . qCD+((MonodB . MonodC) . MonodD) . qDE+(((MonodB . MonodC) . MonodD) . MonodE) . qEA;
NV := Array([]):
NV0 := Array([]):
NVnumb := Array([]):
Monod := [MonodA, MonodB, MonodC, MonodD, MonodE]:
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
NW := [wA, wB, wC, wD, wE]:
A := Array([]):
nv := [seq(subs(k = j,NV[j]),j = 1 .. 5)]:
v := nv[5]+nv[4]+nv[3]+nv[2]+nv[1]:
tD := [seq(seq(d(k)[j] = 0,j = 1 .. NVnumb[k]),k = 1 .. 5)]:
tn := numelems(tD):
tM := Matrix([seq(seq([seq(IdentityMatrix(tn)[add(NVnumb[k],k = 1 .. tx-1)+ty,tz],tz = 1 .. tn),
    `^`(10,20)*seq(subs(tD,subs(d(tx)[ty] = 1,v[tz])),tz = 1 .. 4)],ty = 1 .. NVnumb[tx]),tx = 1 .. 5)]):
tL := LLL(tM,integer):
print(tL);
for s from 1 to tn do
    if tL[s,tn+1]=0 and tL[s,tn+2]=0 and tL[s,tn+3]=0 and tL[s,tn+4]=0 then
        v := [seq(seq(d(tx)[ty] = tL[s,ty+add(NVnumb[j],j = 1 .. tx-1)],ty = 1 .. NVnumb[tx]),tx = 1 .. 5)];
        # print;
        ta := subs(v,nv[1] . NW[1]+nv[2] . NW[2]+nv[3] . NW[3]+nv[4] . NW[4]+nv[5] . NW[5]):
        if ta<>0 then
            print(ta):
            Append(A,ta):
        end if:
    else
        break;
    end if:
end do:
numelems(A);
Eres1 := Basis2MaxSearchCount(convert(A,list)):
for s from 1 to numelems(Eres1) do
    if abs(Re(Eres1[s]))<10^(-30) then
        Eres1[s] := I*Im(Eres1[s]);
    end if:
    if abs(Im(Eres1[s]))<10^(-30) then
        Eres1[s] := Re(Eres1[s]);
    end if:
end do:
evalf(Vector([Eres1]),30);
# order of magnitude of error is about 10^-34.699
NV := Array([]):
NV0 := Array([]):
NVnumb := Array([]):
Monod := [MonodA, MonodB, MonodC, MonodD, MonodE]:
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
NW := [wA, wB, wC, wD, wE]:
A := Array([]):
nv := [seq(subs(k = j,NV[j]),j = 1 .. 5)]:
v := nv[5]+nv[4]+nv[3]+nv[2]+nv[1]:
tD := [seq(seq(d(k)[j] = 0,j = 1 .. NVnumb[k]),k = 1 .. 5)]:
tn := numelems(tD):
tM := Matrix([seq(seq([seq(IdentityMatrix(tn)[add(NVnumb[k],k = 1 .. tx-1)+ty,tz],tz = 1 .. tn),
    `^`(10,20)*seq(subs(tD,subs(d(tx)[ty] = 1,v[tz])),tz = 1 .. 4)],ty = 1 .. NVnumb[tx]),tx = 1 .. 5)]):
tL := LLL(tM,integer):
print(tL);
for s from 1 to tn do
    if tL[s,tn+1]=0 and tL[s,tn+2]=0 and tL[s,tn+3]=0 and tL[s,tn+4]=0 then
        v := [seq(seq(d(tx)[ty] = tL[s,ty+add(NVnumb[j],j = 1 .. tx-1)],ty = 1 .. NVnumb[tx]),tx = 1 .. 5)];
        # print;
        ta := subs(v,nv[1] . NW[1]+nv[2] . NW[2]+nv[3] . NW[3]+nv[4] . NW[4]+nv[5] . NW[5]):
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
