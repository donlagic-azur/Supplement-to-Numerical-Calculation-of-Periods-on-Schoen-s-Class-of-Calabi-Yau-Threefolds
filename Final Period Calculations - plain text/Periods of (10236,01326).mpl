read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
ResAB := Vector([-2.355366068381452168163886096942358108328208959057987181629030165167534483366372686983055914313535321+`*`(4.960896220901313523851118929173559000011942113975462985751071530786581688790278060524743364270837589,I),
    2.355366068381452168163886096942358108328208959057987181629030165167534483366372686983055914313535321+`*`(4.960896220901313523851118929173559000011942113975462985751071530786581688790278060524743364270837589,I),
    5.908217577212413828147015285043517624150281335825048269420009340217576846789042610819807254857964429+`*`(3.489233495349459617160392475124949108692413333897809108876641860803914385061136333999751306782638234,I),
    5.908217577212413828147015285043517624150281335825048269420009340217576846789042610819807254857964429-`*`(3.489233495349459617160392475124949108692413333897809108876641860803914385061136333999751306782638234,I)]):
ResBC := Vector([2.649586944747155220636984411582236289568004821398566297853484188647729445022220005386429712746883021,
    -`*`(3.570375633469668177638593084553728512736306770276689946616871311580346700885895587285759250912205513,I),
    -`*`(3.570375633469668177638593084553728512736306770276689946616871311580346936922647139288409154662517912,I),
    -3.867674320183361994383077335340326379596548207454397645279566114491965140291173469103433124614938621]):
ResCD := Vector([2.355366068381452168163886096942358108328208959057987181629030165167533547746992536172320966008642510-`*`(4.960896220901313523851118929173559000011942113975462985751071530786581959575252504409914711695757906,I),
    -5.908217577212413828147015285043517624150281335825048269420009340217578169132630369453752491426977546-`*`(3.489233495349459617160392475124949108692413333897809108876641860803914122368108028010988206574903652,I),
    2.355366068381452168163886096942358108328208959057987181629030165167536638273042103810856966153008689+`*`(4.960896220901313523851118929173559000011942113975462985751071530786581272758489907587778253165174670,I),
    5.908217577212413828147015285043517624150281335825048269420009340217579237828911733811044635457256325-`*`(3.489233495349459617160392475124949108692413333897809108876641860803916259185990782511265524046742278,I)]):
ResDE := Vector([.9032645640838064393461447765189232262540675553684947899374949864023128701450578495490806566810647971,
    -`*`(2.768002132936499137809505655250057212339410801651366992551881670238450001150266976715622538399888923,I),
    -`*`(2.274036226085022946519813883352281191716424748702976830939419311324564214730341961612913265294250990,I),
    -7.842054737304094825338553952713680813948555439552247092059792397463559452937293717766209394570017356]):
ResEA := Vector([.9032645640838064393461447765189232262540675553684947899374949864023128701450578495490806566810647971,
    -`*`(2.274036226085022946519813883352281191716424748702976830939419311324564214730341961612913265294250990,I),
    -`*`(2.768002132936499137809505655250057212339410801651366992551881670238450001150266976715622538399888923,I),
    -7.842054737304094825338553952713680813948555439552247092059792397463559452937293717766209394570017356]):
MatA := Matrix(4,4,[[-1, -1, 0, 0], [1, 0, 0, 0], [0, 0, -1, -1], [0, 0, 1, 0]]):
MatB := Matrix(4,4,[[-1, 1, -1, 1], [1, -2, 1, -2], [0, 0, -1, 1], [0, 0, 1, -2]]):
MatC := Matrix(4,4,[[-1, -1, 1, 1], [0, -1, 0, 1], [2, 2, -1, -1], [0, 2, 0, -1]]):
MatD := Matrix(4,4,[[1, 0, 0, 0], [0, 1, 0, 0], [1, 0, -1, 0], [0, 1, 0, -1]]):
MatE := Matrix(4,4,[[1, 0, 0, 0], [-3, 1, 0, 0], [-3, 0, 1, 0], [9, -3, -3, 1]]):
ResEA+MatE . ResDE+(MatE . MatD) . ResCD+((MatE . MatD) . MatC) . ResBC+(((MatE . MatD) . MatC) . MatB) . ResAB;
MonodA := Matrix(4,4,[[2, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 2, 1], [0, 0, -1, 0]]);
MonodB := Matrix(4,4,[[-2, 3, -4, 6], [-3, 4, -6, 8], [0, 0, -2, 3], [0, 0, -3, 4]]);
MonodC := Matrix(4,4,[[-12, 32, -9, 24], [-8, 20, -6, 15], [9, -24, 6, -16], [6, -15, 4, -10]]);
MonodD := Matrix(4,4,[[4, 0, 1, 0], [0, 4, 0, 1], [-9, 0, -2, 0], [0, -9, 0, -2]]);
MonodE := Matrix(4,4,[[1, 6, 0, 0], [0, 1, 0, 0], [-6, -36, 1, 6], [0, -6, 0, 1]]);
qAB := ResAB:
qBC := 1/MatB . ResBC:
qCD := (1/MatB . (1/MatC)) . ResCD:
qDE := ((1/MatB . (1/MatC)) . (1/MatD)) . ResDE:
qEA := (((1/MatB . (1/MatC)) . (1/MatD)) . (1/MatE)) . ResEA:
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
# order of magnitude of error is about 10^-34.222
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
