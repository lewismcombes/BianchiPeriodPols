
// finds the smallest field over which a given matrix has all its eigenvalues 
// used when the matrix is defined over some number field 
SmallestFieldNF:=function(mat)
  F:=BaseRing(mat);
  Q,m:=sub< F|[F!1] >;
  S:=Subfields(F) cat [<QNF(),m>];
  pol:=CharacteristicPolynomial(mat);

  // the fields do not come to us sorted by degree, so we must sort them 
  degs:=[Degree(F)/Degree(u[1]) : u in S];
  ParallelSort(~degs,~S);

  // we run over the subfields and find the one with smallest degree such that 
  // the polynomial splits completely (meaning deg_sum = 0)
  best:=F;
  for s in S do 
    // this bit is for weird compatibility reasons 
    if IsIsomorphic(s[1],QNF()) then 
      matIsIn:=&and[u in Q : u in Eltseq(mat)];
    else
      matIsIn:=&and[u in s[1] : u in Eltseq(mat)];
    end if;

    if matIsIn then 
      pol_sub:=ChangeRing(pol,s[1]);
      fac:=Factorization(pol_sub);
      deg_sum:=&+[Degree(u[1])-1 : u in fac];
      if deg_sum eq 0 and Degree(s[1]) lt Degree(best) then 
        best:=s[1];
      end if;
    end if;

  end for;

  return best;
end function;


// same as above, but for when the matrix is defined over a finite field 
SmallestFieldFF:=function(mat)
  p:=Characteristic(BaseRing(mat));
  F1:=Parent(mat[1,1]);
  if IsPrime(#BaseRing(mat)) then
    deg:=1;
  else 
    deg:=2;
  end if;
  // I don't love how big this number could get but I don't see how to not use it
  // polynomials can split into factors of any length smaller than the dimension of the matrix
  d:=deg*LCM([1..Ncols(mat)]);
  F:=ext<F1|d>;

  pol:=CharacteristicPolynomial(mat);

  best:=F;
  for s in Reverse(Exclude(Divisors(d),d)) do 
    FF:=ext<F1|s>;
    if &and[u in FF : u in Eltseq(mat)] then 
      
      pol_sub:=ChangeRing(pol,FF);
      fac:=Factorization(pol_sub);
      deg_sum:=&+[Degree(u[1])-1 : u in fac];
      if deg_sum eq 0 then 
        best:=FF;
      end if;

    end if;
  end for;

  return best;
end function;



// given a list of matrices, we make the smallest field they all split over 
MakeHeckeFieldSmall:=function(list : Optimize:=false)
  if Characteristic(BaseRing(list[1])) eq 0 then 
    F:=QNF();
    for u in list do 
      FF:=SmallestFieldNF(u);
      F:=Compositum(F,FF);
    end for;

    if Optimize then 
      ZF:=MaximalOrder(F);
      F:=OptimizedRepresentation(F);
    end if;

  else
    F:=BaseRing(list[1]);
    for u in list do 
      FF:=SmallestFieldFF(u);
      if Degree(FF) gt Degree(F) then 
        F:=FF;
      end if;
    end for;
  end if;

  return [ChangeRing(u,F) : u in list],F;
end function;


//
// Haluk's code for computing Heilbronn matrices
// previously used magma's QuadraticField functionality
// now featuring some extra stuff to make it work with FldNum
//
Heilbronn:=function(J,level)

        ZK:=Order(level);
        K:=NumberField(ZK);

        // this section translates between number field and quadratic field
        DD:=Discriminant(K);
        Z:=Integers();
        if DD in [-4,-8] then
          K1:=QuadraticField(Z!(DD/4));
        else
          K1:=QuadraticField(Z!DD);
        end if;
        O:=MaximalOrder(K1);

        t,m:=IsIsomorphic(K,K1);

        t,g:=IsPrincipal(J);
        w:=O.2;
        J:=O*m(K!g);



        t,pi:=IsPrincipal(J);
        pi:=O!pi;
        a:=Z!pi[1];
        b:=Z!pi[2];
        q:=Norm(J);

        List:= [**];



          if IsPrime(q) then
	       for  s in [0..q-1]  do
		 x1 := pi; x2 := -s;  y1 := 0 ; y2 := 1;
		 Append(~List, [x1,x2, y1, y2]);

                 a:=-pi; b:=s;

		 while b ne 0 do
                        r:= a mod b;
                        q:= O! ((a-r) / b);
                        x3:= -x1+q*x2;
                        x1 := x2;  x2 := x3;
                        y3 := -y1+q*y2;
                        y1 := y2;  y2 := y3;

			Append(~List,[ x1, x2, y1, y2]);

                        a:=-b; b:=r;
	          end while;
               end for;
               Append(~List, [1, 0, 0,  pi]);
          else
               p:=PrimeDivisors(q)[1];
               for  s,t in [0..p-1]  do
		  x1 := pi; x2 := -(s+t*w);  y1 := 0 ; y2 := 1;
		  Append(~List, [x1,x2, y1, y2]);

                  a:=(-pi); b:=(s+t*w);

		  while b ne 0 do
                        r:= a mod b;
                        q:= O! ((a-r) / b);
                        x3:= -x1+q*x2;
                        x1 := x2;  x2 := x3;
                        y3 := -y1+q*y2;
                        y1 := y2;  y2 := y3;

			Append(~List,[ x1, x2, y1, y2]);

                        a:=-b; b:=r;
	          end while;
                end for;
                Append(~List, [1, 0, 0,  pi]);
          end if;



// this translates back
List2:=[];
mm:=Inverse(m);
for l in List do
  Append(~List2,[mm(l[1]),mm(l[2]),mm(l[3]),mm(l[4])]);
end for;


return List2;
end function;


//
// Computes the action of the Hecke operator at the prime ideal P. returns the
// operator as it acts on the basis of W given by (1 0 ... 0), (0 1 0 ... 0), ...
// (0 0 ... 1), and the operator as it acts on the actual vectors of W as well. Both
// are useful depending on what you want to do with them, so we keep both.
//
Hecke:=function(W,P)

  T:=Heilbronn(P,W`level);
  H:=[Act(W`spec,Matrix(MaximalOrder(W`field),2,2,[t[4],-t[2],-t[3],t[1]])) : t in T];

  HPB:=&+H;
  // forces compatibility, should always come back true
  tt:=IsIsomorphic(CoefficientRing(Parent(HPB)),NumberField(CoefficientRing(Domain(W`down))));

  if W`char ne 0 then
    HPB:=W`down(HPB);
  end if;

  BW:=Basis(W`space);

  HP:=Matrix([Solution(Matrix(BW),b*HPB) : b in BW ]);

  // in case one wants to do things mod p
  if W`char eq 0 then
    F:=BaseRing(W`space);
  else
    F:=CoefficientRing(W`space);
  end if;
  return ChangeRing(HP,F),HPB;
end function;



//
// Simple iterating function that returns the Hecke matrices attahed to W associated
// to the prime ideals in list
//
GetHeckeMatrices:=function(W,list)

  HH:=[];
  HHB:=[];

  for P in list do
    time T,TB:=Hecke(W,P);
    Append(~HH,T);
    Append(~HHB,TB);
  end for;

  return HH,HHB;

end function;


//
// Given a list of commuting matrices (coming in our case from Hecke operators),
// returns the simultaneous eigenvalues of the matrices, as well as their field of
// definition and their multiplicities. Uses Wiese's ArtinAlgebras.
// can sometimes take a long time to run, if there are eigenvalues in a large degree field.
// this is partly because it is not written with efficiency in mind.
//
GET_EV:=function(list)

  char:=Characteristic(CoefficientRing(list[1]));

  A:=MatrixAlgebra(list);


  if Dimension(A) eq 0 then // return the 0 eigenvalue system with maximum multiplicity
    F:=CoefficientRing(list[1]);
    P<x>:=PolynomialRing(F);
    EV_list:=[<x-1,Ncols(list[1]),[F!0 : i in [1..#list]]>];
    basis:=Basis(Kernel(list[1]));
  else

    D:=Decomposition(A);
    EV_list:=<>;

    if #D ne 1 then
      pieces:=[**];
      for j in [1..#D] do
        TT:=[BaseChange(T,D[j]): T in list];
        if char eq 0 then
          F:=OptimizedRepresentation(SplittingField(&*[MinimalPolynomial(u) : u in TT]));
        else
          F<w>:=SplittingField(&*[MinimalPolynomial(u) : u in TT]);
        end if;
        TT:=[ChangeRing(T,F) : T in TT];
        // we re-call this function but with one fewer "piece"
        // this recursiveness is probably not good for the efficiency of the program
        // but this was the least-peverse way I could think of to get it working
        EVs:=$$(TT);
        for e in EVs do
          Append(~pieces,e);
        end for;
      end for;

      return pieces;

    else

      // in the case there is only one eigenvalue system, this part is used
      for j in [1..#D] do
        TT:=[BaseChange(T,D[j]): T in list];
        m:=NumberOfRows(TT[1]);
        ev_size:=[];
        for M in TT do
          total:=0;
          for item in Eigenvalues(M) do
            total:=total+item[2];
          end for;
          Append(~ev_size,total);
        end for;

          if m eq Min(ev_size) then
            Zx<x>:=PolynomialRing(Integers());
              if char eq 0 then
                fx:=Zx!DefiningPolynomial(BaseRing(TT[1]));
              else
                fx:=DefiningPolynomial(BaseRing(TT[1]));
              end if;
            mult:=SetToSequence(Eigenvalues(TT[1]))[1][2];
            ev:=[SetToSequence(Eigenvalues(M))[1][1] : M in TT];
            Append(~EV_list,<fx,mult,ev>);
          end if;
      end for;
      basis:=&cat[Rows(D[i,1]) : i in [1..#D]];
    end if;
  end if;

  return EV_list;
end function;



OptimizeEigenvalues:=function(EVs)
  if Characteristic(Parent(EVs[1])) eq 0 then 
    F:=Parent(EVs[3][1]);
    Q,m:=sub<F|[F!1]>;
    S:=Subfields(F) cat [<QNF(),m>];
    degs:=[Degree(F)/Degree(u[1]) : u in S];
    ParallelSort(~degs,~S);

    best:=F;
    for u in S do 
      if IsIsomorphic(QNF(),u[1]) then
        testField:=Q;
      else
        testField:=u[1];
      end if;
      if &and[v in testField : v in EVs[3]] then 
        best:=testField;
      end if;
    end for;

    return <DefiningPolynomial(best),EVs[2],[best!u : u in EVs[3]]>;
  else 
    deg:=Degree(EVs[1]);
    p:=Characteristic(Parent(EVs[1]));

    best:=GF(p^deg);
    for u in Reverse(Exclude(Divisors(deg),deg)) do 
      FF:=GF(p^u);
      if &and[u in FF : u in EVs[3]] then 
        best:=FF;
      end if;
    end for;

    if Degree(best) eq Degree(EVs[1]) then 
      return EVs;
    else 
      return <DefiningPolynomial(best),EVs[2],[best!u : u in EVs[3]]>;
    end if;
  end if;

end function;





//
// Given a list of commuting matrices and vals in the form <pol, multiplicity, eigs>, finds
// a basis of simultaneous generalised eigenvectors.
// this function has gone awry many times. if there is a problem it is probably hiding here
//
GenEigVecs:=function(W,HH,vals)

  n:=Ncols(HH[1]); //matrices in list are always sqaure and the same size
  d:=vals[2]; // dimension of the space we want
  e:=vals[3];

  FF:=SplittingField(ChangeRing(vals[1],CoefficientRing(W`space)));
  HHe:=[ChangeRing(HH[i],FF) : i in [1..#HH]];
  SP:= [ChangeRing(HH[i],Parent(e[i]))-e[i]*ChangeRing(HH[i]^0,Parent(e[i])) : i in [1..#HHe]];

  MM:=[];
  for M in SP do
    M1:=M;
    i:=1;
    while Rank(M1) gt (n-d) or i lt 100 do
      M1*:=M;
      i+:=1; // just in case
    end while;
    Append(~MM,M1);
  end for;

  return Basis(&meet [Kernel(m) : m in MM]);
end function;






//
// Given Hecke matrices HH and their associated primes HP, returns vectors of W
// along with their Hecke eigenvalues at the primes in HP, and the generators of
// the primes in HP. This uses Wiese's ArtinAlgebras
//
GetPolVals:=function(W,HH,HP)


  KK:=Kernel(ChangeRing(HH[1]-HH[1],CoefficientRing(W`space))); //this is kinda dumb but it gets the job done
  h:=hom<W`space -> KK | [KK.i : i in [1..Dimension(KK)]] >;
  g:=Inverse(h); // map from <(1 0 ... 0) ... (0 0 ... 1)> to W

  // first we get the eigenvalues themselves, with no talk of multiplicity
  EVs:=GET_EV(HH);
  EVs_nomult:=[**];
  for e in EVs do
    for i in [1..e[2]] do
      Append(~EVs_nomult,e[3]);
    end for;
  end for;

  list:=[**];

  // now we gather the (generalised) eigenvectors corresponding to each EV system
  for e in EVs do
    LL:=GenEigVecs(W,HH,e);
    for l in LL do
      Append(~list,l);
    end for;
  end for;

  //gathers the generators of the hecke operators
  gens:=[];
  for i in HP do
    t,gen:=IsPrincipal(i);
    Append(~gens,gen);
  end for;

  // puts everything together
  pol_vals:=[];
  for i in [1..#list] do
    l:=list[i];
    WW:=ChangeRing(W`space,Parent(EVs_nomult[i,1]));
    ll:=ElementToSequence(l);
    Append(~pol_vals,[*&+[ll[j]*WW.j : j in [1..W`dim]],EVs_nomult[i],gens*]);
  end for;

  EV_systems:=[**];
  for e in EVs do
    Append(~EV_systems,OptimizeEigenvalues(<e[1],e[2],e[3]>));
  end for;



  return EV_systems,pol_vals;
end function;



// some utility functions for putting ideals into normal form
detect_hnf:=function(J,M)
      N:=Norm(J);
      a:=M[1,1];  d:=M[1,2];
	  b:=M[2,1];  c:=M[2,2];
	  if (d eq 0) and (N eq a*c) and (b in [0..a-1]) then
	     return 1;
	  else
         return 0;
      end if;
end function;

HNF_basis:=function(J)
      N:=Norm(J);
      M:=BasisMatrix(J);
	  Mt:=Matrix(Integers(),2,2,[M[1,2],M[1,1],M[2,2],M[2,1]]);
      HN:=HermiteForm(Mt);
	  H:=Matrix(Integers(),2,2,[HN[2,2],HN[2,1],HN[1,2],HN[1,1]]);
      c:=H[2,2];  b:=H[2,1];   a:=H[1,1];
	  if c lt 0 then
	     c:=-c;   b:=-b;
	  end if;

	  b:=(b mod a);
	  assert 1 eq detect_hnf(J,Matrix(Integers(),2,2,[a,0,b,c]));
	  return [Norm(J), b,c];
end function;


//
