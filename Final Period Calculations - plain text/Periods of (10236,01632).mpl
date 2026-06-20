read "PeriodHeader.mpl":
with(PeriodHeader);
with(ArrayTools):
with(ListTools):
with(IntegerRelations):
with(LinearAlgebra):
interface(rtablesize = 70):
Digits := 100;
qAB := Vector([.9679511179195162876369199870857504430095,
    `+`(0,-`*`(2.572789479004615250672507360957092966287,I)),
    `+`(0,-`*`(.5868488182113089051184060385694197826852,I)), -1.441193947493577250079607713185223078578]):
qBC := Vector([1.944151666397092721764622877051537279735+`*`(3.696377193925902896465690468059260698818,I),
    `+`(0,-`*`(7.392754387851805792931380936118521397636,I)),
    5.088751150221702139142118152798543123405-`*`(2.185434077233118267230861812475025741428,I),
    -10.17750230044340427828423630559708624681+`*`(0,I)]):
qCD := Vector([-3.888303332794185443529245754103074559469+`*`(4.370868154466236534461723624950051482857,I),
    11.66490999838255633058773726230922367841-`*`(8.741736308932473068923447249900102965714,I),
    -.7437038489695760261517504783560687157986+`*`(10.25267942562525769815827590548433792310,I),
    7.319862697130430217597369587866749270801-`*`(24.20173604517641829278224227902793654503,I)]):
qDE := Vector([-3.377096183332609825353447687356723964597-`*`(1.399091842581997440435695283818253400917,I),
    8.690094602504252225980735348884948815213+`*`(5.370973164168610131543897928593599768121,I),
    -4.345047301252126112990367674442474407607-`*`(.8122430243706885353172892452488336182316,I),
    11.59394795626280108889149531014220014424+`*`(4.197275527745992321307085851454760202750,I)]):
qEA := Vector([4.353296731810186259481150577322510801322-`*`(6.668153505810141990491718809191058780758,I),
    -20.35500460088680855656847261119417249362+`*`(13.33630701162028398098343761838211756152,I),
    `+`(0,-`*`(6.668153505810141990491718809191058780758,I)),
    -7.295114405456249778125020879226640089654+`*`(20.00446051743042597147515642757317634227,I)]):
MonodA := Matrix(4,4,[[1, 0, 2, 0], [-6, 1, -12, 2], [0, 0, 1, 0], [0, 0, -6, 1]]);
MonodB := Matrix(4,4,[[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]]);
MonodC := Matrix(4,4,[[7, 3, 0, 0], [-12, -5, 0, 0], [-21, -9, 7, 3], [36, 15, -12, -5]]);
MonodD := Matrix(4,4,[[-5, 0, 4, 0], [0, -5, 0, 4], [-9, 0, 7, 0], [0, -9, 0, 7]]);
MonodE := Matrix(4,4,[[-35, -10, 42, 12], [90, 25, -108, -30], [-42, -12, 49, 14], [108, 30, -126, -35]]);
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
# order of magnitude of error is about 10^-34.398
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
