

AttachSpec("../ArtinAlgebras/ArtinAlgebras.spec");

load "../ProjAction.m";
load "../WeightAction.m";
load "../Space.m";
load "../Hecke.m";
load "../PeriodPols.m";
load "../H2.m";


d:=11;

level:=[1,0];
weight:=[10,10, 0, 0 ];
char:=0;
HB:=30;
chi:=1;
type:="GL";


time W:=PolSpace(d,level,weight,char,type,chi);

K:=W`field;
ZK:=W`ord;
level:=W`level;


if W`dim eq 0 then
  print "No EV systems found";
else
  if char eq 0 then
    HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level) eq 1*ZK];
  else
    HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level*char) eq 1*ZK];
  end if;
  //HP:=PrimesUpTo(HB,K);
  //HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level) eq 1*ZK];

  HNF:=[HNF_basis(J): J in HP];
  ParallelSort(~HNF,~HP);

  print "finding Hecke matrices...";
  HH,HHB:=GetHeckeMatrices(W,HP);
  print "finding eigenvalue systems...";
  EV_systems,pol_vals:=GetPolVals(W,HH,HP);

  EV_systems;
end if;


// now we do computations with the period polynomials
F:=NumberField(EV_systems[Index([Degree(u[1]) : u in EV_systems],4)][1]);
ZF:=MaximalOrder(F);

// the period polynomial of delta
r_Delta:=VecToPol(W,pol_vals[Index([u[3][1] : u in EV_systems],252)][1]);

irr_pols:=[ pol_vals[i] : i in [1..W`dim] | IsIsomorphic(CoefficientRing(pol_vals[i][1]),F) ];
// the period polynomials of F1 and F2
r_F1:=VecToPol(W,irr_pols[Index([u[3][1] : u in EV_systems | IsIsomorphic(Parent(u[3][1]),F)],2*F.1^2 + 3375)][1]);
r_F2:=VecToPol(W,irr_pols[Index([u[3][1] : u in EV_systems | IsIsomorphic(Parent(u[3][1]),F)],-2*F.1^2 - 4075)][1]);

// creates the H2 space
H2,m:=H2quo(W`spec);
// we only need one hecke operator to split up the space
HH3:=H2Hecke(W`spec,H2,m,HP[1]);

// now we collect the various corresponding polynomials from H2
v_Delta:=VecToPol(W,Inverse(m)(Kernel(HH3-252).1));

// for compatibility reasons, Inverse(m) won't accept a vector defined over F
// so we can cheat and do its job for it
v1:=Kernel(HH3-(2*F.1^2 + 3375)).1;
v_F1:=VecToPol(W,&+[v1[i]*ChangeRing(Inverse(m)(H2.i),F) : i in [1..4]]);

v2:=Kernel(HH3-(-2*F.1^2 - 4075)).1;
v_F2:=VecToPol(W,&+[v2[i]*ChangeRing(Inverse(m)(H2.i),F) : i in [1..4]]);


// first we do Delta
tt,gg:=IsPrincipal(ScalePol(weight,r_Delta)*ScalePol(weight,v_Delta));
D:=Rationals()!gg;

DeltaPair:=PairPol(weight,r_Delta,v_Delta) / D;
print "Pairing of r_Delta and v_Delta:", DeltaPair;
print "Factorization of numerator:", Factorization(Integers()!Numerator(DeltaPair));



// next we do F1
F1Pair:=PairPol(weight,r_F1,v_F1)*ZF / (ScalePol(weight,r_F1) * ScalePol(weight,v_F1));
print "(Norm of numerator of) Pairing of r_F1 and v_F1:", Factorization(Integers()!Numerator(Norm(F1Pair)));


// and F2, which is essentially the same
F2Pair:=PairPol(weight,r_F2,v_F2)*ZF / (ScalePol(weight,r_F2) * ScalePol(weight,v_F2));
print "(Norm of numerator of) Pairing of r_F2 and v_F2:", Factorization(Integers()!Numerator(Norm(F2Pair)));


// for the eiseinstein-geuine-cusp congruences, we note the following:
print "(Norm of) Denominator of leading coeff of Delta:", Factorization(Numerator(Norm(ScalePol(weight,r_Delta))));
print "(Norm of) Denominator of leading coeff of F1:", Factorization(Numerator(Norm(ScalePol(weight,r_F1))));


//
