

/*

AttachSpec("/mnt/data/smp19lmc_code/ArtinAlgebras/ArtinAlgebras.spec");
ChangeDirectory("/mnt/data/smp19lmc_code/period_pols/");

*/

AttachSpec("ArtinAlgebras/ArtinAlgebras.spec");

load "ProjAction.m";
load "WeightAction.m";
load "Space.m";
load "Hecke.m";
load "PeriodPols.m";


d:=11;

for a,b,c,e in [0..4] do

    level:=[2,3];
    weight:=[a,b,c,e];
    char:=5;
    chi:=1;
    HB:=30;
    type:="GL";
    optimizeHF:=true;


    W:=PolSpace(d,level,weight,char,type,chi);

    K:=W`field;
    ZK:=W`ord;
    level:=W`level;


    if W`dim eq 0 then

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

      print "computing Hecke matrices...";
      HH,HHB:=GetHeckeMatrices(W,HP);
      HH,F:=MakeHeckeFieldSmall(HH : Optimize:=optimizeHF);
      print "computing eigenvalue systems and algebraic eigen-polynomials...";
      EV_systems,pol_vals:=GetPolVals(W,HH,HP);


      print "found the following eigenvalue systems:";
      EV_systems;
      chi;
    end if;
    weight;
    "";
    //chi;
end for;

//
