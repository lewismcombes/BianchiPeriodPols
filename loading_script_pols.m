

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


print "computing space of period polynomials...";
=======


>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
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
  print "computing Hecke matrices...";
  HH,HHB:=GetHeckeMatrices(W,HP);
  HH,F:=MakeHeckeFieldSmall(HH : Optimize:=optimizeHF);
  print "computing eigenvalue systems and algebraic eigen-polynomials...";
  EV_systems,pol_vals:=GetPolVals(W,HH,HP);

  print "found the following eigenvalue systems:";
=======
  HH,HHB:=GetHeckeMatrices(W,HP);
  EV_systems,pol_vals:=GetPolVals(W,HH,HP);

>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
  EV_systems;
end if;



//
