

<<<<<<< HEAD
AttachSpec("../../ArtinAlgebras/ArtinAlgebras.spec");
=======
AttachSpec("ArtinAlgebras/ArtinAlgebras.spec");
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4

load "ProjAction.m";
load "WeightAction.m";
load "Space.m";
load "Hecke.m";
load "PeriodPols.m";


<<<<<<< HEAD
<<<<<<< HEAD
d:=2;

level:=[-7,4];
weight:=[0,0,0,0];
=======
d:=11;

level:=[1,0];
weight:=[10, 10, 0, 0];
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
char:=0;
HB:=30;
chi:=1;
type:="GL";
<<<<<<< HEAD
optimizeHF:=true;
=======
d:=1;

level:=[4,-7];
weight:=[0,0, 0, 0];
char:=0;
HB:=30;
chi:=1;
type:="SL";
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8


print "computing space of period polynomials...";
=======


<<<<<<< HEAD
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
=======
print "computing space of period polynomials...";
>>>>>>> f6cc33ad9b6bf601e9db3ff5120b8e427f4ef87a
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

<<<<<<< HEAD
<<<<<<< HEAD
  print "computing Hecke matrices...";
  HH,HHB:=GetHeckeMatrices(W,HP);
  HH,F:=MakeHeckeFieldSmall(HH : Optimize:=optimizeHF);
  print "computing eigenvalue systems and algebraic eigen-polynomials...";
  EV_systems,pol_vals:=GetPolVals(W,HH,HP);

  print "found the following eigenvalue systems:";
=======
=======
  print "computing Hecke matrices...";
>>>>>>> f6cc33ad9b6bf601e9db3ff5120b8e427f4ef87a
  HH,HHB:=GetHeckeMatrices(W,HP);
  print "computing eigenvalue systems and algebraic eigen-polynomials...";
  EV_systems,pol_vals:=GetPolVals(W,HH,HP);

<<<<<<< HEAD
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
=======
  print "found the following eigenvalue systems:";
>>>>>>> f6cc33ad9b6bf601e9db3ff5120b8e427f4ef87a
  EV_systems;
end if;



//
