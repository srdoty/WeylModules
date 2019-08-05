#############################################################################
InstallMethod(SchurAlgebraWeylModule, "for a prime and a partition", true,
[IsInt, IsList], 0, 
function(p,mu)
   # returns a characteristic p Weyl module for a Schur algebra S(n,r) 
   # where r = |mu| and n = length(mu)
   local V, W, n, r;
   n:= Length(mu); 
   if n=1 then n:=2; Add(mu, 0); fi; #n must be at least 2!
   r:= Sum(mu);
   V:= WeylModule(p, CompositionToWeight(mu), "A", n-1);
   W:= Objectify(NewType( FamilyObj(V), IsSchurAlgebraWeylModule ), 
              rec(prime:=p,partition:=mu,degree:=r,size:=n,module:=V) 
              );
  return(W);
end );

#############################################################################
InstallMethod(PrintObj, "for Schur algebra Weyl modules", true, 
[IsSchurAlgebraWeylModule], 0, 
function(W)
  Print("<Schur algebra Weyl module of highest weight ", 
         W!.partition, " at prime p = ", W!.prime, ">");
end );

#############################################################################
InstallMethod(DecompositionNumbers, "for a Schur algebra Weyl module", true, 
[IsSchurAlgebraWeylModule], 0, 
function(W)
   # returns the decomposition numbers of W using partition notation
   local V, d, k, dec;
   d:= W!.degree;
   V:= W!.module;
   dec:= DecompositionNumbers(V);
   for k in [2,4..Length(dec)] do
      dec[k-1]:= WeightToComposition(d, dec[k-1]);
   od;
   return dec;
end );

#############################################################################
InstallMethod(SchurAlgebraDecompositionMatrix, 
"for a Schur algebra S(n,r) in characteristic p", true, 
[IsInt, IsInt, IsInt], 0, 
function(p,n,r)
   # returns the decomposition matrix of a Schur algebra S(n,r) 
   # in characteristic p
   local ptns, size, V, dn, mu, row, k, pos, mat, time;

   ptns:= BoundedPartitions(n,r,r);
   size:= Length(ptns);
   mat:= [];
   for mu in ptns do
      #Print("Creating Weyl module of highest weight ", mu, "\n");
      #time:= Runtimes().user_time;
      V:= SchurAlgebraWeylModule(p, mu);
      #time:= Runtimes().user_time - time; 
      #Print("Last step took ", StringTime(time), " seconds\n");
      #Print("Computing decomposition numbers...\n");
      #time:= Runtimes().user_time;
      dn:= DecompositionNumbers(V);
      #time:= Runtimes().user_time - time; 
      #Print("Last step took ", StringTime(time), " seconds\n");
      row:= List([1..size], i -> 0);
         for k in [2,4..Length(dn)] do
            pos:= Position(ptns, dn[k-1]);
            if pos = fail then 
               return fail; 
            else 
               row[pos]:= dn[k];
            fi;
         od;
      Add(mat, row);
   od;
   return mat;
end );

#############################################################################
InstallMethod(SymmetricGroupDecompositionNumbers,
"for a prime and a partition", true, [IsInt, IsList], 0, 
function(p, mu)
   # returns the decomposition numbers of the dual Specht
   # module labeled by a partition mu, in characteristic p
   local deg, dnums, k, ptn, wt, ans, strip, V;

   strip:= function(x)
      local c;
      c:= ShallowCopy(x);
      while c[Length(c)] = 0 do Remove(c); od;
      return(c);
   end;

   deg:= Sum(mu);
   while Length(mu) < deg do Add(mu,0); od;
   V:= SchurAlgebraWeylModule(p, mu);
   dnums:= DecompositionNumbers(V);
   ans:= [];
   for k in [2,4..Length(dnums)] do
      ptn:= dnums[k-1];
      wt:= CompositionToWeight(ptn);
      if IsRestrictedWeight(p, wt) then 
         Add(ans, strip(ptn));
         Add(ans, dnums[k]);
      fi;
   od;
   return ans;
end );

#############################################################################
InstallMethod(SymmetricGroupDecompositionMatrix,
"for a prime and a positive integer", true, [IsInt, IsInt], 0, 
function(p, n)
   # returns the decomposition matrix for the symmetric group
   # on n letters, in characteristic p
   local mat, ptns, out, t, k, wt;
   ptns:= BoundedPartitions(n,n);
   mat:= SchurAlgebraDecompositionMatrix(p,n,n);
   out:= [];
   t:= TransposedMat(mat);
   for k in [1..Length(mat)] do
      wt:= CompositionToWeight(ptns[k]);
      if IsRestrictedWeight(p, wt) then 
         Add(out, t[k]);
      fi;
   od;
   return TransposedMat(out);
end );

#############################################################################
