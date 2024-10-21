
//
// Combines the action on the projective line and the action on the polynomial
// space V_{k,l}(C). Covers both trivial cases.
//
Act:=function(spec,mat)

  PM:=ProjMatChi(spec,mat);
  WM:=WeightMat(spec,mat);

  return TensorProduct(PM,WM);

end function;



//
// Handles setting up the various spaces for d in {1,2,3,7,11}
//
StandardMats:=function(spec)

  ZK:=CoefficientRing(spec`PL[1]);
  K<z>:=NumberField(ZK);

  d:=spec`d;

  if d eq 1 then

    T:=Matrix(ZK,2,2,[1,1,0,1]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    Tz:=Matrix(ZK,2,2,[1,z,0,1]);
    L:=Matrix(ZK,2,2,[z,0,0,-z]);
    J:=Matrix(ZK,2,2,[z,0,0,1]);

    TT:=Act(spec,T);
    TS:=Act(spec,S);
    TTz:=Act(spec,Tz);
    TL:=Act(spec,L);
    TJ:=Act(spec,J);
    ID:=TT^0;

    TE:=Act(spec,Tz*S*L);

    return [
    ID+TS, ID-TL, ID+TT*TS+(TT*TS)^2, ID+TE+TE^2, ID-TJ
    ];

  elif d eq 2 then

    T:=Matrix(ZK,2,2,[1,1,0,1]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    Tz:=Matrix(ZK,2,2,[1,z,0,1]);
    Tzi:=Matrix(ZK,2,2,[1,-z,0,1]);
    J:=Matrix(ZK,2,2,[-1,0,0,1]);

    TT:=Act(spec,T);
    TS:=Act(spec,S);
    TTz:=Act(spec,Tz);
    TTzi:=Act(spec,Tzi);
    TJ:=Act(spec,J);
    ID:=TT^0;

    TTzS:=TTz*TS; //a little bit of time saving

    return [
    ID+TS, ID+TT*TS+(TT*TS)^2, ID + TS*TTz + TTzS + TTzi*TS*TTzS, ID-TJ
    ];

  elif d eq 3 then

    T:=Matrix(ZK,2,2,[1,1,0,1]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    Tz:=Matrix(ZK,2,2,[1,z,0,1]);
    L:=Matrix(ZK,2,2,[z,0,0,-z^2]);
    J:=Matrix(ZK,2,2,[z,0,0,1]);

    TT:=Act(spec,T);
    TS:=Act(spec,S);
    TTz:=Act(spec,Tz);
    TL:=Act(spec,L);
    TJ:=Act(spec,J);
    ID:=TT^0;

    TE:=Act(spec,T^-1*Tz*S*L);

    return [
    ID+TS, ID-TL, ID+TT*TS+(TT*TS)^2, ID+TE+TE^2, ID-TJ
    ];

  elif d eq 7 then

    T:=Matrix(ZK,2,2,[1,1,0,1]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    Tz:=Matrix(ZK,2,2,[1,z,0,1]);
    J:=Matrix(ZK,2,2,[-1,0,0,1]);

    TT:=Act(spec,T);
    TS:=Act(spec,S);
    TTz:=Act(spec,Tz);
    TJ:=Act(spec,J);
    ID:=TT^0;

    return [
    ID+TS, ID+TT*TS+(TT*TS)^2, TT+TS*TTz+TTz*TS*TT+TS*TTz^-1*TS*TTz, ID-TJ
    ];

  elif d eq 11 then

    T:=Matrix(ZK,2,2,[1,1,0,1]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    Tz:=Matrix(ZK,2,2,[1,z,0,1]);
    J:=Matrix(ZK,2,2,[-1,0,0,1]);
    E:=Tz^-1*S*Tz*S*T;
    F:=Matrix(ZK,2,2,[-1,0,0,-1]);

    TT:=Act(spec,T);
    TS:=Act(spec,S);
    TTz:=Act(spec,Tz);
    TE:=Act(spec,E);
    TEi:=Act(spec,E^-1);
    TJ:=Act(spec,J);
    ID:=TT^0;

    return [
    ID+TS, ID+TT*TS+(TT*TS)^2, TT + TS*TTz + TT*TE + TS*TTz*TEi + TTz*TS*TT+TS*TTz^-1*TS*TTz, ID-TJ
    ];


  else
    print "bad choice for Q(sqrt(-d))";
    return -1;
  end if;
end function;



//
// Defines the space of polynomials invariant under the action of Gamma_0(level).
// Uses the matrices in Sengun ExpMath.
//
PolSpace:=function(d,level,weight,char,type,chi)

  // this will be the final format we output the modular form space in
  PolData:= recformat <
    space             : ModTupFld,
    level             : RngOrdIdl,
    chi               : GrpDrchNFElt,
    weight            : SeqEnum,
    dim               : RngIntElt,
    id_index          : RngIntElt,
    d                 : RngIntElt,
    field             : FldNum,
    ord               : RngOrd,
    char              : RngIntElt,
    down              : Map,
    spec              : Rec >;

  // this format type allows us to package up all the data specifying which space
  // we care about: d for the field, level, weight, character, projective line,
  // and the field all of the matrices are going to live in
  SpecData:= recformat <
    d                 : RngIntElt,
    level             : RngOrdIdl,
    weight            : SeqEnum,
    chi               : GrpDrchNFElt,
    PL                : SetIndx,
    r                 : UserProgram,
    field             : FldNum >;




  K<z>:=QuadFld(d);
  ZK:=MaximalOrder(K);
  level:=ZK!level;

  PL,r:=ProjectiveLine(quo<ZK|level>);
  // for compatibility reasons these lines needs to be here
  // otherwise ProjActionChi can't evaluate chi at the
<<<<<<< HEAD
<<<<<<< HEAD
  // scalar given by r.
=======
  // scalar given by r. 
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
=======
  // scalar given by r.
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8
  level:=Parent(1*CoefficientRing(PL[1]))!level;
  ZK:=Order(level);
  K:=NumberField(ZK);


  DG:=Elements(DirichletGroup(level));
  chi:=DG[chi];

  spec:=rec<SpecData | d:=d,
                       level:=level,
                       weight:=weight,
                       chi:=chi,
                       PL:=PL,
                       r:=r,
<<<<<<< HEAD
                       field:=Compositum(ImageField(level,chi),K) >;
=======
                       field:=Compositum(Codomain(chi),K) >;
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4


  M:=StandardMats(spec);

  // in characteristic p we want to send all of our matrices that give the action
  // to the residue class field of a prime over p*ZK. it doesn't really matter which
  if char ne 0 then
<<<<<<< HEAD
<<<<<<< HEAD
    KK:=Compositum(BaseRing(M[1]),Codomain(chi));
    ZKK:=MaximalOrder(KK);
    //M:=[ChangeRing(m,ZKK) : m in M];
=======
    KK:=Compositum(K,Codomain(chi));
    ZKK:=MaximalOrder(KK);
    M:=[ChangeRing(m,ZKK) : m in M];
>>>>>>> 349a8eb58628495365844c9a2a6782c577e390d4
=======
    KK:=Compositum(BaseRing(M[1]),Codomain(chi));
    ZKK:=MaximalOrder(KK);
    //M:=[ChangeRing(m,ZKK) : m in M];
>>>>>>> 9726bf3cd58d932c5d58afb9fa987fbd9bec1eb8
    F,down:=ResidueClassField(Factorization(char*ZKK)[1,1]);
  else
    down:=IdentityHomomorphism(Parent(M[1]));
  end if;

  // we always put the [e,0,0,1] relation at the end of the list
  // so we can do SL calculations if we want by just omitting the last relation
  // here e generates the unit group of ZK
  if type eq "GL" then
    W:=&meet [Kernel(down(u)) : u in M];
  elif type eq "SL" then
    W:=&meet [Kernel(down(u)) : u in M[1..#M-1]];
  end if;




  W:=rec<PolData | space:=W,
                   level:=level,
                   chi:=chi,
                   weight:=weight,
                   dim:=Dimension(W),
                   id_index:=IdIndex(level),
                   d:=d,
                   field:=K,
                   ord:=ZK,
                   char:=char,
                   down:=down,
                   spec:=spec
                   >;

  return W;

end function;










//
