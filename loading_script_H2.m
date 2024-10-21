

AttachSpec("ArtinAlgebras/ArtinAlgebras.spec");

load "ProjAction.m";
load "WeightAction.m";
load "Space.m";
load "Hecke.m";
load "H2.m";


<<<<<<< HEAD
<<<<<<< HEAD
d:=1;

level:=[1,-7];
weight:=[0, 0, 0, 0];
HB:=30;
chi:=11;
=======
d:=11;
=======
d:=1;
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8

level:=[1,-7];
weight:=[0, 0, 0, 0];
HB:=30;
<<<<<<< HEAD
chi:=1;
>>>>>>> f6cc33ad9b6bf601e9db3ff5120b8e427f4ef87a
=======
chi:=11;
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8

K:=QuadFld(d);
ZK:=MaximalOrder(K);
level:=(K!level)*ZK;
PL,r:=ProjectiveLine(quo< ZK|level >);
chi:=Elements(DirichletGroup(level))[chi];


SpecData:= recformat <
  d                 : RngIntElt,
  level             : RngOrdIdl,
  weight            : SeqEnum,
  chi               : GrpDrchNFElt,
  PL                : SetIndx,
  r                 : UserProgram,
<<<<<<< HEAD
<<<<<<< HEAD
  char_field        : FldNum,
  field        : FldNum >;
=======
  field             : FldNum >;
>>>>>>> f6cc33ad9b6bf601e9db3ff5120b8e427f4ef87a
=======
  char_field        : FldNum,
  field        : FldNum >;
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8


spec:=rec<SpecData | d:=d,
                     level:=level,
                     weight:=weight,
                     chi:=chi,
                     PL:=PL,
                     r:=r,
<<<<<<< HEAD
<<<<<<< HEAD
                     char_field:=Compositum(Codomain(chi),K),
                     field:=K >;
=======
                     field:=Compositum(Codomain(chi),K) >;
>>>>>>> f6cc33ad9b6bf601e9db3ff5120b8e427f4ef87a
=======
                     char_field:=Compositum(Codomain(chi),K),
                     field:=K >;
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8



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
