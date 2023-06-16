

AttachSpec("ArtinAlgebras/ArtinAlgebras.spec");

load "ProjAction.m";
load "WeightAction.m";
load "Space.m";
load "Hecke.m";
load "H2.m";


d:=11;

level:=[1,0];
weight:=[10, 10, 0, 0];
HB:=30;
chi:=1;

K:=QuadFld(d);
ZK:=MaximalOrder(K);
level:=(K!level)*ZK;
PL,r:=ProjectiveLine(quo<ZK|level>);
chi:=Elements(DirichletGroup(level))[chi];


SpecData:= recformat <
  d                 : RngIntElt,
  level             : RngOrdIdl,
  weight            : SeqEnum,
  chi               : GrpDrchNFElt,
  PL                : SetIndx,
  r                 : UserProgram,
  field             : FldNum >;


spec:=rec<SpecData | d:=d,
                     level:=level,
                     weight:=weight,
                     chi:=chi,
                     PL:=PL,
                     r:=r,
                     field:=Compositum(Codomain(chi),K) >;



print "computing H2...";
time H2,m:=H2quo(spec);




if Dimension(H2) eq 0 then
  print "No EV systems found";
else
  HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level) eq 1*ZK];

  HNF:=[HNF_basis(J): J in HP];
  ParallelSort(~HNF,~HP);

  HH:=[];
  print "finding Hecke matrices...";
  for P in HP do
    time Append(~HH,H2Hecke(spec,H2,m,P));
  end for;
  print "finding eigenvalue systems...";
  EV_systems:=GET_EV(HH);
  print "found the following eigenvalue systems:";
  print EV_systems;
end if;



//
