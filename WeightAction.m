
//
// Computes the action of mat on the weight module. Variable weight is a 4 entry list
// [k,l,a,b], where k and l are the degrees of the polynomials and a and b are the
// twists on the two spaces V_k, V_l respectively.
//
WeightMat:=function(spec,mat);

  weight:=spec`weight;

  R:=CoefficientRing(Parent(mat));
  K:=NumberField(R);

  P<x,y>:=PolynomialRing(spec`field,2);

  Symm:=function(k,d,T)
    ST:=ZeroMatrix(spec`field,k+1,k+1);
    for i in [0..k] do
      Q:=(T[1,1]*x+T[1,2]*y)^(k-i)*(T[2,1]*x+T[2,2]*y)^(i);
      for j in [0..k] do
        ST[i+1,j+1]:=MonomialCoefficient(Q,x^(k-j)*y^(j));
      end for;
    end for;
    return Determinant(T)^d*ST;
  end function;


  k:=weight[1];
  l:=weight[2];
  a:=weight[3];
  b:=weight[4];

  con:=Automorphisms(K)[2];


  matc:=Matrix(R,2,2,[R!con(mat[1][1]),R!con(mat[1][2]),R!con(mat[2][1]),R!con(mat[2][2])]);
  TM:=TensorProduct(Symm(k,a,mat),Symm(l,b,matc));

  return TM;

end function;
