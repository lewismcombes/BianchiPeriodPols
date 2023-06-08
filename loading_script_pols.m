

AttachSpec("ArtinAlgebras/ArtinAlgebras.spec");

load "ProjAction.m";
load "WeightAction.m";
load "Space.m";
load "Hecke.m";
load "PeriodPols.m";


d:=11;

level:=[1,0];
weight:=[10, 10, 0, 0];
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

  HH,HHB:=GetHeckeMatrices(W,HP);
  EV_systems,pol_vals:=GetPolVals(W,HH,HP);

  EV_systems;
end if;



//
