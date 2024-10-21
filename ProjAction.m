// does the same thing as BaseRing which doesn't work for dirichlet chars over 
// number fields 
ImageField:=function(level,chi)

  // there HAS to be a better way to do this -_-
  ZK:=Parent(Generators(level)[1]);
  QQ,down:=quo<ZK|level>;

  images:=[chi(Inverse(down)(u)) : u in Exclude(Elements(QQ),0)];
  pols:=[MinimalPolynomial(u) : u in images];
  if #pols eq 0 then 
    return Rationals();
  else 
    return SplittingField(&*pols);
  end if;
end function;


//
// Gives the action of M on the projective line P^1(ZK/J).
//
ProjMat:=function(J,M)

  ZK:=CoefficientRing(Parent(M));

  if Type(ZK) ne RngOrd then
    ZK:=MaximalOrder(ZK);
  end if;

  SS:=MatrixAlgebra(ZK,2);
  M:=SS!M;

  // this line covers weird magma compatibility problems that can sometimes arise
  t:=IsIsomorphic(NumberField(CoefficientRing(Parent(M))),NumberField(Order(J)));

  PL,r:=ProjectiveLine(quo<ZK|J>);

  proj_action:=function(M,i)
    t,im :=r(PL[i]*M,true,false);
    if t then
      return Index(PL,im);
    else
      return -1; //this picks out the bad elements
    end if;
  end function;

  perm:=[];
  for i in [1..#PL] do
    Append(~perm,proj_action(M,i));
  end for;

  l:=#PL;
  Mat:=[];
  for i in [1..l] do
    new_row:=[];
    for j in [1..l] do
      if j eq perm[i] then
        Append(~new_row,1);
      else
        //we ignore the elements with problems at p, per Wang 94
        // this allows for computing Hecke operators at primes dividing J
        Append(~new_row,0);
      end if;
    end for;
    Append(~Mat,new_row);
  end for;

  return ChangeRing(Matrix(Mat),NumberField(ZK));

end function;




//
// This is the usual ProjMat we know and love. Gives the action of M on the
// projective line P^1(ZK/J).
// We supply the projective line and its decision function r because
// Magma sometimes generates different orderings of elements in the PL
//
ProjMatChi:=function(spec,M)

  ZK:=CoefficientRing(Parent(M));

  if Type(ZK) ne RngOrd then
    ZK:=MaximalOrder(ZK);
  end if;

  SS:=MatrixAlgebra(ZK,2);
  M:=SS!M;

  PL:=spec`PL;
  r:=spec`r;

  proj_action:=function(M,i)
    t,im,scal :=r(PL[i]*M,true,true);
    if t then
      return Index(PL,im),scal;
    else
      return -1,0; //this picks out the bad elements
    end if;
  end function;

  perm:=[];
  scalars:=[];
  for i in [1..#PL] do
    pa,sc:=proj_action(M,i);
    Append(~perm,pa);
    Append(~scalars,sc);
  end for;

  l:=#PL;
  Mat:=[];
  for i in [1..l] do
    new_row:=[];
    for j in [1..l] do
      if j eq perm[i] then
        Append(~new_row,spec`field!spec`chi(scalars[i]));
      else
        Append(~new_row,spec`field!0); //we ignore the elements with problems at p, per Wang 94
      end if;
    end for;
    Append(~Mat,new_row);
  end for;

  //return ChangeRing(Matrix(Mat),NumberField(ZK));
  return ChangeRing(Matrix(Mat),spec`field);

end function;




//
// A simple function that returns the index of (0 : 1) so we know where to look for our period polynomial
//
IdIndex:=function(level)

  ZK:=Ring(Parent(level));
  PL,r:=ProjectiveLine(quo<ZK|level>);
  return Index(PL,Vector([ZK!0,ZK!1]));

end function;


//
// for compatibility reasons, we use the same fields that the LMFDB uses
//
QuadFld:=function(d)

    _<x>:=PolynomialRing(Integers());
  if d in [1,2] then
    return NumberField(x^2+d);
  elif d in [3,7,11] then
    return NumberField(x^2-x-Integers()!((-d-1)/4));
  else
    return "d = " cat Sprint(d) cat " currently not supported";
  end if;
end function;









//
