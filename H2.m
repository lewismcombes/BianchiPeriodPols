


// computes H2
H2quo:=function(spec)

  ZK:=MaximalOrder(spec`field);
  z:=spec`field.1;
  d:=spec`d;

  if d eq 1 then

    E:=Matrix(ZK,2,2,[-1,z,z,0]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    SL:=Matrix(ZK,2,2,[0,z,z,0]);
    U:=Matrix(ZK,2,2,[1,-1,1,0]);
    J:=Matrix(ZK,2,2,[z,0,0,1]);

    TE:=Act(spec,E);
    TS:=Act(spec,S);
    TSL:=Act(spec,SL);
    TU:=Act(spec,U);
    TJ:=Act(spec,J);

    KER:=Kernel(1-TE) + Kernel(1-TS) + Kernel(1-TSL) + Kernel(1-TU) + Kernel(1+TJ);

  elif d eq 2 then

    A:=Matrix(ZK,2,2,[1,z,z,-1]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    U:=Matrix(ZK,2,2,[1,-1,1,0]);
    Tz:=Matrix(ZK,2,2,[1,z,0,1]);
    J:=Matrix(ZK,2,2,[-1,0,0,1]);

    TA:=Act(spec,A);
    TS:=Act(spec,S);
    TU:=Act(spec,U);
    TTzi:=Act(spec,Tz^(-1));
    TJ:=Act(spec,J);

    KER:=Kernel(1-TA)*(1-TTzi) + Kernel(1-TS) + Kernel(1-TU) + Kernel(1+TJ);

  elif d eq 3 then

    LS:=Matrix(ZK,2,2,[0,-z,1-z,0]);
    U:=Matrix(ZK,2,2,[1,-1,1,0]);
    SL:=Matrix(ZK,2,2,[0,-1+z,z,0]);
    J:=Matrix(ZK,2,2,[z,0,0,1]);

    TLS:=Act(spec,LS);
    TU:=Act(spec,U);
    TSL:=Act(spec,SL);
    TJ:=Act(spec,J);

    KER:=Kernel(1-TLS) + Kernel(1-TU) + Kernel(1-TSL) + Kernel(1+TJ);

  elif d eq 7 or d eq 11 then

    if d eq 7 then
      A:=Matrix(ZK,2,2,[1,-1+z,z,-1]);
    else
      A:=Matrix(ZK,2,2,[1,-1+z,z,-2]);
    end if;
    g:=Matrix(ZK,2,2,[0,-1,1,-z]);
    S:=Matrix(ZK,2,2,[0,-1,1,0]);
    U:=Matrix(ZK,2,2,[1,-1,1,0]);
    J:=Matrix(ZK,2,2,[-1,0,0,1]);

    TA:=Act(spec,A);
    TS:=Act(spec,S);
    TU:=Act(spec,U);
    Tg:=Act(spec,g);
    TJ:=Act(spec,J);

    KER:=Kernel(1-TS) + Kernel(1-TU) + Kernel(1-TA)*(1+Tg) + Kernel(1+TJ);

  end if;

  V:=RSpace(CoefficientRing(KER),Degree(KER));
  H2,m:=quo<V|KER>;

  return H2,m;

end function;

// computes hecke operators on H2
H2Hecke:=function(spec,H2,m,P)

  T:=Heilbronn(P,spec`level);
  H:=[Act(spec,Matrix(spec`field,2,2,t)) : t in T];
  M:= Matrix([m(Inverse(m)(H2.i) * &+H) : i in [1..Dimension(H2)]]);

  return M;
end function;




//
