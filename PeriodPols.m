



PeriodPol:=function(W,vec : ind:=W`id_index)

  // if the user wants a different component polynomial, they can ask for it, otherwise
  // we give them the component associated to the identity element
  // this is the same as the polynomial attached to {0,\infty}
  if not assigned ind then
    ind:=W`id_index;
  end if;
  weight:=W`weight;
  k:=weight[1];
  l:=weight[2];
  R:=CoefficientRing(W`space);

  D:=(k+1)*(l+1);

  coeffs:=[vec[D*(ind-1)+i] : i in [1..D]];


  M:=Matrix(R,k+1,l+1,coeffs); //turns the coeffs into a k+1-by-l+1 matrix, which is more pleasant to look at
  return M;



end function;


//
// Returns the central coefficient(s) of the period polynomial
//
CentralCoeffs:=function(W,vec)
  weight:=W`weight;
  k:=weight[1];
  l:=weight[2];

  pol,M:=PeriodPol(W,vec);

  k2:=k mod 2;
  l2:=l mod 2;

  i:=(k-k2) div 2 + 1;
  j:=(l-l2) div 2 + 1;

  return Submatrix(M,i,j,k2+1,l2+1);

end function;


//
// pairs two monomials
//
PairMon:=function(weight,mon1,mon2)

  PP1:=Parent(mon1);
  PP2:=Parent(mon2);
  degx1:=Degree(mon1,PP1.1);
  degxb1:=Degree(mon1,PP1.3);

  degx2:=Degree(mon2,PP2.1);
  degxb2:=Degree(mon2,PP2.3);

  k:=weight[1];
  l:=weight[2];

  if degx1+degx2 eq k and degxb1+degxb2 eq l then
    return (-1)^(degx1+degxb1) * Binomial(k,degx1)^(-1) * Binomial(k,degxb1)^(-1);
  else
    return 0;
  end if;
end function;

//
// breaks two polynomials up into monomials, pairs them all, then adds it all up
//
PairPol:=function(weight,pol1,pol2)

  mons1:=Monomials(pol1);
  mons2:=Monomials(pol2);

  coeffs1:=Coefficients(pol1);
  coeffs2:=Coefficients(pol2);

  pair:=0;

  for i in [1..#mons1] do
    for j in [1..#mons2] do
      pair+:=coeffs1[i]*coeffs2[j]*PairMon(weight,mons1[i],mons2[j]);
    end for;
  end for;

  return pair;
end function;


//
// input a vector (from pol_vals, or elsewhere) and get a polynomial back
// to do: make the index variable do anything. right now this only works for
// level 1 forms.
//
VecToPol:=function(W,vec : ind:=W`id_index)

  P<X,Y,Xb,Yb>:=PolynomialRing(CoefficientRing(vec),4);
  c:=ElementToSequence(vec);
  pol:=0;
  k:=W`spec`weight[1];
  l:=W`spec`weight[2];
  for i in [1..k+1] do
    for j in [1..l+1] do
      pol+:=c[(i-1)*k+j+(i-1)]*X^(k-i+1)*Y^(i-1)*Xb^(l-j+1)*Yb^(j-1);
    end for;
  end for;
  return pol;
end function;


ScalePol:=function(weight,pol)
  F:=CoefficientRing(Parent(pol));
  ZF:=MaximalOrder(F);

  // here we make sure to not include the first and last terms in the scaling.
  // this needs to work for any polynomial, not just those coming from pol_vals,
  // hence we define the two monomials, then take them away
  PP:=Parent(pol);
  first_mon:=(PP.1*PP.3)^weight[1];
  last_mon:=(PP.2*PP.4)^weight[2];

  pol2:=pol - MonomialCoefficient(pol,first_mon)*first_mon - MonomialCoefficient(pol,last_mon)*last_mon;

  coeffs:=Coefficients(pol2);
  ideals :=[ideal<ZF|coeffs[i]> : i in [1..#coeffs]];
  D1:=GCD(ideals[1],ideals[2]);
  for i in [3..#ideals] do
    D1:=GCD(D1,ideals[i]);
  end for;
  return D1;
end function;





//
