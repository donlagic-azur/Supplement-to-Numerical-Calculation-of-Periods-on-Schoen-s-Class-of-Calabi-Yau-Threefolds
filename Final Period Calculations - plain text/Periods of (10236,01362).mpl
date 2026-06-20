read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
qAB := Vector([3.627504369507369709643261636914544939650-`*`(8.321124798398730259912583416922026981361,I),
    3.627504369507369709643261636914544939650+`*`(8.321124798398730259912583416922026981361,I),
    -8.907977187423762623758158364640315243364-`*`(4.745569190638751838446696313134790419722,I),
    8.907977187423762623758158364640315243364-`*`(4.745569190638751838446696313134790419722,I)]):
qBC := Vector([1.011064598800979536257835504776719334573+`*`(0,I),
    1.011064598800979536257835504776719334573+`*`(2.485957425617978130013102288818581002831,I),
    `+`(0,-`*`(2.447862660430063128391385611680462075986,I)),
    6.524474829567187467416569153323383490256-`*`(2.447862660430063128391385611680462075986,I)]):
qCD := Vector([9.315429976166523001733606139025089057071+`*`(0,I),
    -18.63085995233304600346721227805017811414+`*`(11.98691846509117440061248747147718427860,I),
    13.97314496424978450260040920853763358561-`*`(9.679389276686793096851184213633833504955,I),
    -11.64708066571794611115593778068454257009+`*`(37.33915625101034779462109963448344342781,I)]):
qDE := Vector([3.945604474573995832278421144924697994659+`*`(5.218005089806316545661470529109392710108,I),
    -10.82334136538602377500025115363882350311-`*`(8.036064942566177558259264680065464694567,I),
    5.918406711860993748417631717387046991988+`*`(5.218005089806316545661470529109392710108,I),
    -14.76894583996001960727867229856352149777-`*`(6.836092324042949791727426490988804331744,I)]):
qEA := Vector([-17.89960341904886807991312442564105133053+`*`(3.103119708592413714251112887812634279575,I),
    24.81563234941072053256636628999773733845-`*`(14.75793574654170523227890849715232757655,I),
    -10.98357448868701562725988256128436532260+`*`(11.65481603794929151802779560933969329697,I),
    10.98357448868701562725988256128436532260-`*`(23.30963207589858303605559121867938659394,I)]):
MonodA := Matrix(4,4,[[-2, -3, -4, -6], [3, 4, 6, 8], [0, 0, -2, -3], [0, 0, 3, 4]]);
MonodB := Matrix(4,4,[[2, -1, 0, 0], [1, 0, 0, 0], [0, 0, 2, -1], [0, 0, 1, 0]]);
MonodC := Matrix(4,4,[[1, 0, 0, 0], [6, 1, 0, 0], [-3, 0, 1, 0], [-18, -3, 6, 1]]);
MonodD := Matrix(4,4,[[-5, 0, 4, 0], [0, -5, 0, 4], [-9, 0, 7, 0], [0, -9, 0, 7]]);
MonodE := Matrix(4,4,[[15, 10, -18, -12], [-40, -25, 48, 30], [18, 12, -21, -14], [-48, -30, 56, 35]]);
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
# order of magnitude of error is about 10^-34.301
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
