#############################################################################
##
#W  weylmod.gi                   GAP package              S.R. Doty
##
##
#Y  Copyright (C)  2009,  S.R. Doty
##
##  This file contains the implementation of methods for
##  Weyl modules, simple characters, etc.
##
#############################################################################
InstallMethod(WeylModule, 
"for a prime <p>, a list <wt>, and a simple Lie algebra of type <t>, rank <r>",
true, [IsPosInt, IsList, IsString, IsPosInt], 0, 
function(p, wt, t, r) 
  local V, W, L;
  L:= SimpleLieAlgebra(t,r,Rationals);
  V:= HighestWeightModule(L, wt);
  W:= Objectify(NewType( FamilyObj(V), IsWeylModule ), 
              rec(prime:=p,highestWeight:=wt,type:=t,rank:=r,LieAlgebra:=L,
                  module:=V,maximalVecs:=[],maximalVecsAmbiguous:=[],
		  simpleQuotient:=[]) 
              );
  return(W);
end );

#############################################################################
InstallMethod(WeylModule, 
"for a prime <p>, a list <wt>, and a simple Lie algebra of type <t>, rank <r>",
true, [IsWeylModule,IsList], 0, 
function(M, wt)
# returns a Weyl module of highest weight <wt> of the same type as <M>
# in particular, with the SAME underlying Lie algebra  
  local V, W, L,p,t,r;
  L:= M!.LieAlgebra;
  p:= M!.prime;
  t:= M!.type;
  r:= M!.rank;
  V:= HighestWeightModule(L, wt);
  W:= Objectify(NewType( FamilyObj(V), IsWeylModule ), 
              rec(prime:=p,highestWeight:=wt,type:=t,rank:=r,LieAlgebra:=L,
                  module:=V,maximalVecs:=[],maximalVecsAmbiguous:=[],
		  simpleQuotient:=[])  
      );
  return(W);
end );

#############################################################################
InstallMethod(PrintObj, "for Weyl modules", true, 
[IsWeylModule], 0, 
function(W)
  Print("<Type ",  W!.type, W!.rank, " Weyl module of highest weight ", 
         W!.highestWeight, " at prime p = ", W!.prime, ">");
end );

#############################################################################
InstallMethod(IsAmbiguous,  "for a Weyl module", true, 
[IsWeylModule], 0, 
function(W)
  if Length(W!.maximalVecsAmbiguous) > 0 then return true; fi;
end );

#############################################################################
InstallMethod(AmbiguousMaxVecs,  "for a Weyl module", true, 
[IsWeylModule], 0, 
function(W)
  return(W!.maximalVecsAmbiguous);
end );

#############################################################################
InstallMethod(TheLieAlgebra,  "for a Weyl module", true, 
[IsWeylModule], 0, 
function(W)
  return(W!.LieAlgebra);
end );

#############################################################################
InstallMethod(TheCharacteristic,  "for a Weyl module", true, 
[IsWeylModule], 0, 
function(W)
  return(W!.prime);
end );

#############################################################################
InstallMethod(BasisVecs, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(W)
 # returns the standard basis of <W>
 return BasisVectors(Basis(W!.module));
end );

#############################################################################
InstallMethod(Generator, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(W)
 # returns the generator of <W>
 return BasisVecs(W)[1];
end );

#############################################################################
InstallMethod(Dim, "for Weyl modules", true, 
[IsWeylModule], 0, 
function(W)
 # returns the dimension of <W>
 return Length(BasisVecs(W));
end );

#############################################################################
InstallMethod(Weights, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns a list of the weight space labels of <V>
 return( DuplicateFreeList(List(BasisVecs(V), Weight)) );
end );

#############################################################################
InstallMethod(DominantWeights, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns a list of the dominant weight space labels of <V>
 local out,item;
 out:= [];
 for item in Weights(V) do
   if IsDominant(item) then Add(out, item); fi;
 od;
 return(out);
end );

#############################################################################
InstallMethod(WeightSpace, "for a Weyl module", true, 
[IsWeylModule,IsList], 0, 
function(V,wt)
 # returns a basis for the given weight space in <V> of weight <wt>
 local bb,out,b;
 bb:=BasisVecs(V); out:=[ ];
 for b in bb do
   if Weight(b) = wt then Add(out,b); fi;
 od;
 return(out);
end );

#############################################################################
InstallMethod(WeightSpaces, "for a Weyl module", true,
[IsWeylModule], 0,
function(V)
 # returns a list of the weights and weight spaces of <V>
 local wts,out,w;

 wts:=Weights(V);
 out:= [ ];
 for w in wts do
   Add(out, w);
   Add(out, WeightSpace(V,w));
 od;
 return(out);
end );

#############################################################################
InstallMethod(DominantWeightSpaces, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns a list of the dominant weight spaces of V
 local out, k, ws;
 out:= [ ];
 ws:= WeightSpaces(V);
 for k in [2,4..Length(ws)] do
   if IsDominant( ws[k-1] ) then 
     Add(out, ws[k-1]);
     Add(out, ws[k]);
   fi;
 od;
 return(out);
end );

#############################################################################
InstallMethod(ActOn, 
"for a Weyl module, an element of the hyperalgebra, and a vector",
true, [IsWeylModule, IsUEALatticeElement, IsLeftAlgebraModuleElement], 0, 
function(V, h, v) 
 # returns the result of acting by <h> on <v>
 local p;
 p:= TheCharacteristic(V);
 return ((h^v) mod p);
end );

#############################################################################
InstallMethod(Character, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns the character (a list of the weights and their multiplicities)
 # for the given WeyModule V
 local ws, k, out;
 out:= [ ];
 ws:= WeightSpaces(V);
 for k in [2,4..Length(ws)] do
   Add(out, ws[k-1]);
   Add(out, Length(ws[k]));
 od;
 return(out);
end );
   
#############################################################################
InstallMethod(RowVec,
"for a Weyl module and a given weight vector", true,
[IsWeylModule,IsLeftAlgebraModuleElement], 0, 
function(V,vec)
  # converts <vec> to a row vector with mod p coefficients
  local e,k,out,p;
  p:= V!.prime;
  out := [ ];
  for k in [1..Dim(V)] do Add(out, 0*Z(p)^0); od;
  e:= ExtRepOfObj(ExtRepOfObj(vec));
  for k in [2,4..Length(e)] do
     out[ e[k-1][1] ] := e[k]*Z(p)^0;
  od;
  return( out );
  end );

#############################################################################
InstallMethod(SimpleLieAlgGens,
"for a WeylModule", true,
[IsWeylModule], 0,
function(V)
   # returns xsimple & ysimple
   local L, rank, noPosRoots, k, g, xsimple, ysimple;
   L:=V!.LieAlgebra;
   g:= LatticeGeneratorsInUEA( L );
   rank:= Length(CanonicalGenerators(RootSystem(L))[1]);
   noPosRoots:= Length(ChevalleyBasis(L)[1]);
   xsimple:= [ ]; ysimple:= [ ];
   for k in [1..rank] do 
     Add(ysimple, g[k]);
     Add(xsimple, g[k+noPosRoots]);
   od;
   return( [xsimple, ysimple] );
 end );

#############################################################################
InstallMethod(IsMaximalVector, "for a Weyl module and weight vector", true, 
[IsWeylModule,IsLeftAlgebraModuleElement], 0, 
function(V,vec)
 # Tests <vec> to see if it is maximal in <V>. See below for a relative
 # version of this function.

 local rank,j,k,p,zerovec,height,xy,xsimple,ysimple;

 p:= V!.prime;
 xy:= SimpleLieAlgGens(V); xsimple:=xy[1]; ysimple:=xy[2];
 zerovec:= 0*vec;
 rank:= Length(xsimple);
 height:= HighestPrimePower(p, Sum(V!.highestWeight));
 for j in [1..rank] do 
     for k in [0..height] do 
       if not (xsimple[j]^(p^k)/Factorial(p^k))^vec mod p = zerovec then 
         return(false); fi;
     od;
 od;
 return(true);
end );

#############################################################################
# At the moment, operations depending on the following are not likely
# to work if we obtain more than one independent maximal vector in the
# given weight space. In such a case, the multiple maximal vectors in
# question are stored in the Weyl module, and a warning is printed.
#############################################################################
InstallMethod(MaximalVectors, "for a Weyl module and weight", true,
[IsWeylModule,IsList], 0, function(V,wt)
# Returns maximal vectors of a given weight space in a module, if possible

 local rank,i,j,k,p,height,vec,rowvec,finalmatrix,wtspace,
       item,outlist,result,join,S,matrix,xy,xsimple,ysimple,z;

 join:= function(A,B)
   # Given matrices A,B compute an augmented matrix M such that the
   # left nullspace of M equals the intersection of the left nullspaces
   # of A and B. 
   local tA,tB,CS,vec,p; 
   p:= V!.prime;
   tA:= TransposedMatMutable(A); tB:= TransposedMatMutable(B);
   if Length(tA) = 0 then tA:= tB; fi;
   CS:= VectorSpace(GF(p),tA);
   for vec in tB do
     if not vec in CS then 
       Add(tA, vec);
       CS:= VectorSpace(GF(p),tA);
     fi;
   od;
   return( TransposedMatMutable(tA) );
   end; 

 p:= V!.prime;
 wtspace:= WeightSpace(V,wt);
 xy:= SimpleLieAlgGens(V); xsimple:=xy[1]; ysimple:=xy[2];
 rank:= Length(xsimple);
 height:= HighestPrimePower(p, Sum(V!.highestWeight));
 finalmatrix:= [];
 for j in [1..rank] do 
   for k in [0..height] do
     matrix:= []; 
     for vec in wtspace do
       rowvec:= RowVec(V, (xsimple[j]^(p^k)/Factorial(p^k))^vec mod p );
       Add(matrix, rowvec);
     od;
     finalmatrix:= join(finalmatrix, matrix); 
   od;
 od;
 
 z:= NullspaceMatDestructive(finalmatrix);
 outlist:= [];
 for item in z do
   result:= 0*wtspace[1];
   for i in [1..Length(item)] do
     result:= result + IntFFESymm(item[i])*wtspace[i];
   od;
   Add(outlist, result);
 od;
 if Length(outlist) > 1 then
    Add(V!.maximalVecsAmbiguous, outlist);
    if Length(V!.maximalVecsAmbiguous) = 1 then # first time
       Print("***** WARNING: Ambiguous module detected *****\n");
    fi;
 fi; 
 return(outlist);
end );  


#############################################################################
InstallMethod(MaximalVectors, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # Finds a maximal vector in each weight space of given Weyl module <V>, 
 # if possible 
 local k,max,wts,out;
 if V!.maximalVecs = [] then   
   wts:= DominantWeights(V);
   out:= [ ];
   for k in [1..Length(wts)] do
     max:= MaximalVectors(V,wts[k]);
     if Length(max) > 0 then
       Append(out, max);
     fi;
   od;
   V!.maximalVecs:= out; #remember for next time
   return( out );
 else
   return( V!.maximalVecs );
 fi;
end );

#############################################################################
InstallMethod(SocleSeries, "for a Weyl module", true, [IsWeylModule], 0, 
function(V)
 # Returns the socle series of V

 local S, ans;

 S:= SocleWeyl(V);
 ans:= [S];
 
 while Dim(S) < Dim(V) do
   S:= NextSocle(S);
   Add(ans,S);
 od;
 return ans;
end );

#############################################################################
InstallMethod(ExtWeyl, "for a Weyl module and a list of vectors", true, 
[IsWeylModule,IsList], 0, 
function(V,gens)
 # Returns a list of relative maximal (i.e. primitive or maximal) 
 # vectors that generate the socle of <V>/<S>, where <S> is the
 # submodule generated by the given <gens>. 

 local outlist, j, v, w, mvecs, s, Q, subv, vlist, p;
 p:= TheCharacteristic(V);
 s:= SubWeylModule(V,gens);
 Q:= QuotientWeylModule(V,s);
 outlist:= SocleWeyl(Q);
 # Sometimes (e.g. V(2,0) at p=2 for Type G2) the above obtains a 
 # bad choice for a primitive vector, so we may need to replace it 
 # by a better choice, when available.
 for v in gens do
    subv:= SubWeylModule(V,v);
    Q:= QuotientWeylModule(V,subv);
    vlist:= SocleWeyl(Q);
    for w in vlist do
        for j in [1..Length(outlist)] do
           if w <> outlist[j] and Weight(w) = Weight(outlist[j]) 
              and IsWithin(V,s,w-outlist[j] mod p) then
              outlist[j]:= w;
           fi;
        od;
    od;
 od;
 return outlist;
end );

#############################################################################
InstallMethod(PrimitiveVectors, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns a list of the primitive vectors of the given Weyl module

 local prim, s, i,j,sub, vec;

 prim:=[];
 s:= SocleWeyl(V);
 Append(prim, s);
 
 while s[1] <> Generator(V) do
   s:= ExtWeyl(V,prim);
   Append(prim, s);
 od;

 return(prim);
end );

#############################################################################
InstallMethod(SimpleQuotient, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns the quotient of the given Weyl module by the maximal submodule
 local mvecs,S,Q;

 Q:= V!.simpleQuotient;
 if Q <> [ ] then return(Q[1]); fi;
 # else compute and save it as follows
 S:=SubWeylModule(V,[ ]); # zero module
 mvecs:= MaximalVectors(V);
 if Length(mvecs) = 1 then
    return(QuotientWeylModule(S));
 fi;
    
 while Length(mvecs) > 1 do
    S:=SubWeylModule(S,mvecs[2]);
    Q:=QuotientWeylModule(S);
    mvecs:=MaximalVectors(Q);
 od;
 Add(V!.simpleQuotient, Q); # save it
 return(Q);
end );

#############################################################################
InstallMethod(MaximalSubmodule, "for a Weyl module", true,
[IsWeylModule], 0,
function(V)
 # returns the (unique) maximal submodule of <V>
 local Q;
 Q:=SimpleQuotient(V);
 return(DefiningKernel(Q));
end );

#############################################################################
InstallMethod(SimpleTopFactorCharacter, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns the character of the simple top factor of the given 
 # Weyl module
 return( Character(SimpleQuotient(V)) );
end );


#############################################################################
InstallMethod(SimpleTopFactorDim, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns the dimension of the simple top factor of the given 
 # Weyl module

 return( CharacterDim(SimpleTopFactorCharacter(V)) );
end );
 
#############################################################################
InstallMethod(SubmoduleStructure, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # deduces the structure of the given Weyl module

 local prim, s, i,j,sub, level, index, vec;

 index:=0; level:=1; 
 prim:=[]; 
 s:= SocleWeyl(V);
 Append(prim, s);
 Print("Level ", level,"\n");
 for vec in s do
   index:=index+1;
   if IsMaximalVector(V,vec) then
     Print("-maximal vector v", index, 
        " = ", vec, " of weight ", Weight(vec),"\n");
   else
     Print("-primitive vector v", index, 
        " = ", vec, " of weight ", Weight(vec),"\n");
   fi;
 od;

 while s[1] <> Generator(V) do
   level:=level+1;
   s:= ExtWeyl(V,prim);  
   Append(prim, s);
   Print("Level ", level, "\n");
   for vec in s do
     index:=index+1;
     if IsMaximalVector(V,vec) then
       Print("-maximal vector v", index, 
          " = ", vec, " of weight ", Weight(vec),"\n");
     else
       Print("-primitive vector v", index, 
          " = ", vec, " of weight ", Weight(vec),"\n");
     fi;
   od;
 od;

 for i in [1..Length(prim)] do
   sub:= SubWeylModule(V, prim[i]);
   Print("The submodule <v", i, "> contains ");
   for j in [1..Length(prim)] do
     if IsWithin(V,sub,prim[j]) then 
       Print("v", j, " ");
     fi;
   od;
   Print("\n");
 od;
 return prim;
end );

#############################################################################
InstallMethod(DecompositionNumbers, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns the decomposition numbers of the given Weyl module

 local ch, m, wt, mult, ms, decnums,
       maximalWeight, multiple, dominates;

 multiple:= function(scalar, char)
    # returns a <scalar> multiple of the given <char>
    local k, out;
    out:= [];
    for k in [2,4..Length(char)] do
       Add(out, char[k-1]);
       Add(out, scalar*char[k]);
    od;
    return out;
 end;

 dominates:= function(lambda, mu)
    # returns true if <lambda> dominates <mu>, false otherwise
    local cf, R, L, B, space, bas, zero, c;
    L:= TheLieAlgebra(V);
    R:= RootSystem(L);
    bas:= SimpleSystem(R); #the simple roots
    space:= Rationals^Length(bas);
    B:= Basis(space, bas);
    cf:= Coefficients(B, lambda - mu);
    zero:= [];
    for c in cf do
       Add(zero, 0);
    if cf = zero then return false; fi;
    od;
    for c in cf do
       if not c in NonnegativeIntegers then 
          return false; fi;
    od;
    return true;
 end;

 maximalWeight:= function(char)
    # returns a maximal [weight, mult] pair of the given <char>
    local winner, winnermult, k;
    winner:= char[1]; winnermult:= char[2];

    for k in [4,6..Length(char)] do
        if dominates(char[k-1], winner) then
           winner:= char[k-1];
           winnermult:= char[k];
        fi;
    od;
    return [winner, winnermult];
 end;   
 
 decnums:= [];
 ch:= Character(V);
 while ch <> [] do 
     m:= maximalWeight(ch);
     wt:= m[1]; mult:= m[2];
     Add(decnums, wt);
     Add(decnums, mult);
     ms:= multiple(mult, SimpleCharacter(V,wt));
     ch:= DifferenceCharacter(ch, ms); 
 od;
 return decnums;
end );     

#############################################################################