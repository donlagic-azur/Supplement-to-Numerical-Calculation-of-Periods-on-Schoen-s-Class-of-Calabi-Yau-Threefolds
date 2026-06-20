read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
qAB := Vector([.4932005593217494790348026431155872493323+`*`(0,I),
    `+`(0,-`*`(1.304501272451579136415367632277348177527,I)),
    `+`(0,-`*`(.2999931546308069416329595472691650907060,I)),
    -.7330331040595080276108522159473568784473+`*`(0,I)]):
qBC := Vector([1.729007232590463113163310466089171501982+`*`(2.137924082339219450944170680381764754349,I),
    7.932915319943143246304902038588605836596+`*`(3.500068237530335473325563138810370938805,I),
    `+`(0,-`*`(2.137924082339219450944170680381764754349,I)),
    -2.745893622171753906814970640321091330649-`*`(6.413772247017658352832512041145294263048,I)]):
qCD := Vector([-3.133870389232783083350355000388715049805+`*`(.8938889019399946053664717759468091408860,I),
    -7.174616870842408594111525410006066335665+`*`(3.868059003479671775711089406124125029192,I),
    4.947622573986467938171985818845987521916+`*`(1.478895693379381313856876380620586072183,I),
    12.61587342510346315857641786537788375200+`*`(.8775101871590800627356069070106653969459,I)]):
qDE := Vector([-.2527661497002448840644588761941798336434-`*`(.6119656651075157820978464029201155189966,I),
    -.7582984491007346521933766285825395009301-`*`(1.223931330215031564195692805840231037993,I),
    -.8728202582910622146607656597483063716339+`*`(2.457386351727041878796814780964991807697,I),
    -.9873420674813897771281546909140732423378+`*`(5.536262059858578290096905134134628866103,I)]):
qEA := Vector([1.164428747020815375216700767378136132134-`*`(2.419847319171698274212796053408458376239,I),
    `+`(0,-`*`(4.839694638343396548425592106816916752478,I)),
    -4.074802315695405723511220159097681150282-`*`(1.498364808136396800076560933934648034825,I),
    -8.149604631390811447022440318195362300564+`*`(0,I)]):
MonodA := Matrix(4,4,[[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]]);
MonodB := Matrix(4,4,[[1, 0, 2, 0], [-6, 1, -12, 2], [0, 0, 1, 0], [0, 0, -6, 1]]);
MonodC := Matrix(4,4,[[-20, 8, -15, 6], [-72, 28, -54, 21], [15, -6, 10, -4], [54, -21, 36, -14]]);
MonodD := Matrix(4,4,[[4, 0, 1, 0], [0, 4, 0, 1], [-9, 0, -2, 0], [0, -9, 0, -2]]);
MonodE := Matrix(4,4,[[-5, 3, 0, 0], [-12, 7, 0, 0], [30, -18, -5, 3], [72, -42, -12, 7]]);
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
