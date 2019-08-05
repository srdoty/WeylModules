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


# The following function is not documented at this time
#############################################################################
InstallMethod(SortedCharacter, "for a list of weights and multiplicities",
true, [IsList], 0,
function(char)
 local k, wts, mults, ans;
 wts:= []; mults:=[];
 for k in [2,4..Length(char)] do
    Add(wts, char[k-1]);
    Add(mults, char[k]);
 od;
 SortParallel(wts,mults);
 ans:=[];
 for k in [1..Length(mults)] do
    Add(ans, wts[k]);
    Add(ans, mults[k]);
 od;
 return ans;
end );
    

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
InstallMethod(Weight, "for weight vectors", true, 
[IsLeftAlgebraModuleElement], 0, 
function(wtvec)
 # returns the weight of <wtvec>, if possible
 local u,rep,i,wt,u1,u2,v1,v2;
 rep:= ExtRepOfObj(wtvec);
 u:= ExtRepOfObj( rep );
 if IsWeightRepElement(rep) then
    wt:= u[1][3];    
    if Length(u) = 2 then  return( wt ); fi;
    for i in [2,4..Length(u)] do 
      if u[i-1][3] <> wt then return fail; fi; 
    od;
 elif IsTensorElement(rep) then
    v1:= u[1][1]; v2:= u[1][2];
    u1:= ExtRepOfObj( ExtRepOfObj(v1) );
    u2:= ExtRepOfObj( ExtRepOfObj(v2) );
    wt:= u1[1][3] + u2[1][3];
    if Length(u) = 2 then  return( wt ); fi;
    for i in [2,4..Length(u)] do 
      v1:= u[i-1][1]; v2:= u[i-1][2];
      u1:= ExtRepOfObj( ExtRepOfObj(v1) );
      u2:= ExtRepOfObj( ExtRepOfObj(v2) );
      if u1[1][3] + u2[1][3] <> wt then return fail; fi; 
    od;
 else return fail;
 fi;

 return( wt );
end );

#############################################################################
InstallMethod(Weights, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # returns a list of the weight space labels of <V>
 return( DuplicateFreeList(List(BasisVecs(V), Weight)) );
end );

#############################################################################
InstallMethod(IsDominant, "for a weight (i.e. a list)", true, 
[IsList], 0, 
function(wt)
 # returns true iff the given weight is dominant
 local k;
 for k in [1..Length(wt)] do
   if wt[k] < 0 then return false; fi;
 od;
 return true;
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
InstallMethod(IsRestrictedWeight, 
"for a prime <p> and a dominant weight <wt>",
true, [IsPosInt, IsList], 0, 
function(p, wt)
 local k;
 for k in wt do
   if k >= p then return(false); fi;
 od;
 return(true);
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
InstallMethod(DifferenceCharacter, "for two characters", true, 
[IsList,IsList], 0, 
function(ch1,ch2)
 # returns the difference of two characters or "fail" if the 
 # difference is not a character
 local j,new, out;
 out:= [ ];
 for j in [2,4..Length(ch2)] do
   if not ch2[j-1] in ch1 then return fail; fi;
 od;
 for j in [2,4..Length(ch1)] do
   if ch1[j-1] in ch2 then
     new:= ch1[j] - ch2[Position(ch2,ch1[j-1])+1];
     if new < 0 then return fail; fi;
     if new > 0 then 
       Add(out,ch1[j-1]);
       Add(out,new);
     fi;
   else
     Add(out,ch1[j-1]);
     Add(out,ch1[j]);
   fi;
 od;
 return(out);
end );

#############################################################################
InstallMethod(CharacterDim, "for a character", true, 
[IsList], 0, 
function(ch)
 # returns the dimension of the given character, i.e. the sum
 # of the dimensions of the weight spaces
 local out, k;
 out:= 0;
 for k in [2,4..Length(ch)] do
    out:= out + ch[k];
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
InstallMethod(IsWithin, 
"for a sub Weyl module and a given weight vector", true, 
[IsSubWeylModule,IsLeftAlgebraModuleElement], 0, 
function(S,vec)
 # Returns True if the given <vec> lies within the given submodule S
 local V,rowbasis, v, rowvec, SS; 
 V:=S!.weylModule;
 rowbasis:=S!.repbasis;
 SS:= VectorSpace(GF(V!.prime), rowbasis);
 rowvec:= RowVec(V,vec);
 if rowvec in SS then
    return true;
 else
    return false;
 fi;
end );  

#############################################################################
InstallMethod(HighestPrimePower,
"for a prime and a positive integer", true,
[IsPosInt,IsInt], 0,
function(p,n)
   #returns 1 plus the highest power of <p> that is less than or equal to <n>
   local k;
   k:= 0;
   while p^k <= n do k:= k+1; od;
   return(k);
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
InstallMethod(IsMaximalVector, 
        "for a sub Weyl module and a weight vector", true, 
        [IsSubWeylModule,IsLeftAlgebraModuleElement], 0, 
function(S,vec)
 # Tests <vec> to see if it is maximal in <V>/<S>, for the 
 # given submodule <S>. In other words,
 # if this returns "true" then <vec> is primitive in <V>.
 
 local V,rank,j,k,p,zerovec,height,xy,xsimple,ysimple;
 V:=S!.weylModule;
 p:= V!.prime;
 xy:= SimpleLieAlgGens(V); xsimple:=xy[1]; ysimple:=xy[2];
 zerovec:= 0*vec;
 rank:= Length(xsimple);
 height:= HighestPrimePower(p, Sum(V!.highestWeight));
 for j in [1..rank] do 
     for k in [0..height] do 
       if not IsWithin(S,(xsimple[j]^(p^k)/Factorial(p^k))^vec mod p) then
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
InstallMethod(SubWeylModule, "for a Weyl module and vector", true, 
[IsWeylModule,IsLeftAlgebraModuleElement], 0, 
function(W,vec)
 # constructs the submodule of <W> generated by <vec>
 local i,j,k,rowbasis,lbasis,newvec,rowvec,S,L,
       p,extendbasis,g,gens,noPosRoots,height,comps,wt,SM,v;
  
 if vec = Generator(W) then
    lbasis:= BasisVecs(W);
    rowbasis:=[];
    for v in lbasis do
       Add(rowbasis, RowVec(W,v));
    od;   
    SM:=Objectify(NewType(FamilyObj(W), IsSubWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=[vec],weylModule:=W) );
    return(SM);
 fi;

 p:= W!.prime;
 wt:= W!.highestWeight;
 height:= HighestPrimePower(p, Sum(wt));
 L:= W!.LieAlgebra;
 g:= LatticeGeneratorsInUEA( L );
 noPosRoots:= Length(ChevalleyBasis(L)[1]);
 gens:= g{ [1..2*noPosRoots] };
 for i in [1..height] do
    for j in [1..2*noPosRoots] do
       Add(gens, gens[j]^(p^i)/Factorial(p^i));
    od;
 od;

 rowbasis:= [ RowVec(W,vec) ]; lbasis:= [ vec ];
 S:= VectorSpace(GF(p), rowbasis);   

 extendbasis:= function(v)  # helper function
   local i;
   for i in gens do
       newvec:= (i^v) mod p;
       if newvec <> 0*v then  
         rowvec:= RowVec(W,newvec);
         if not (rowvec in S) then
            Add( rowbasis, rowvec);
            S:= VectorSpace(GF(p), rowbasis);
            Add( lbasis, newvec );
         fi;
       fi;
   od;
 end;

 extendbasis(vec); 
 j:=2;
 while j <= Length(lbasis) do
    extendbasis(lbasis[j]); 
    j:= j+1;
 od;
 
 SM:=Objectify(NewType(FamilyObj(W), IsSubWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=[vec],weylModule:=W) );
 return(SM);
end );

#############################################################################
InstallMethod(SubWeylModule, "for a sub Weyl module and a vector", true, 
[IsSubWeylModule, IsLeftAlgebraModuleElement], 0, 
function(W,vec)
 # Returns the submodule of the ambient Weyl module generated by the
 # given submodule <W> and given vector <vec>.

 local V, p, U, rowbasis, lbasis, S, generators, ulbasis, urowbasis, i, SM;

 V:=W!.weylModule; p:=V!.prime;
 lbasis:=ShallowCopy(W!.eltbasis);
 if Length(lbasis) = 0 then  # W is the zero module
    return( SubWeylModule(V,vec) );
 fi;
 rowbasis:=ShallowCopy(W!.repbasis);
 generators:=ShallowCopy(W!.gens);
 Add(generators, vec);

 U:=SubWeylModule(V,vec);
 ulbasis:= U!.eltbasis;
 urowbasis:=U!.repbasis;
 for i in [1..Length(ulbasis)] do
   S:=VectorSpace(GF(p), rowbasis);
   if not (urowbasis[i] in S) then
     Add(lbasis, ulbasis[i]);
     Add(rowbasis, urowbasis[i]);
   fi;
 od;

 SM:=Objectify(NewType(FamilyObj(V), IsSubWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=generators,
     weylModule:=V) );
 return(SM);
end );

#############################################################################
InstallMethod(SubWeylModule, "for a Weyl module and list of vectors", true, 
[IsWeylModule, IsList], 0, 
function(V,vecs)
 # Returns the sub Weyl module generated by the given vectors

 local S, v, SM;
 if Length(vecs) = 0 then
   SM:=Objectify(NewType(FamilyObj(V), IsSubWeylModule),
     rec(eltbasis:=[],repbasis:=[],gens:=[],
     weylModule:=V) );
   return(SM);
 fi;
 
 S:= SubWeylModule(V, vecs[1]); # to start
 for v in vecs{ [2..Length(vecs)] } do
   S:= SubWeylModule(S,v);
 od;
 return(S);
end );

#############################################################################
InstallMethod(SubWeylModuleDirectSum, "for a list of sub Weyl modules", true,
[IsList], 0,
function(inlist)
 # returns the direct sum of the inputs (which are assumed independent)
 # Note: No checking of the assumption is done!

 local V, S, lbasis, rowbasis, generators, SM;

 lbasis:= []; rowbasis:= []; generators:= [];
 V:=inlist[1]!.weylModule;
 for S in inlist do
    Append(lbasis, S!.eltbasis);
    Append(rowbasis, S!.repbasis);
    Append(generators, S!.gens);
    Assert(0,V = S!.weylModule);
 od;

 SM:= Objectify(NewType(FamilyObj(V), IsSubWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=generators,
     weylModule:=V) );
 return(SM);
end );

#############################################################################
InstallMethod(PrintObj, "for a sub Weyl module", true,
[IsSubWeylModule], 0,
function(S)
 Print("<SubWeylModule of dimension ", Dim(S), ", generated by elements ",
 S!.gens, " of weights ", List(S!.gens, Weight), ">, in\n", S!.weylModule);
end );

InstallMethod(Generators, "for a sub Weyl module", true,
[IsSubWeylModule], 0,
function(S)
 return(S!.gens);
end );

InstallMethod(BasisVecs, "for a sub Weyl module", true,
[IsSubWeylModule], 0,
function(S)
 # returns the basis vectors of <S>
 return(S!.eltbasis);
end );

InstallMethod(Dim, "for a sub Weyl module", true,
[IsSubWeylModule], 0,
function(S)
 # returns the dimension of <S>
 return(Length(S!.eltbasis));
end );

InstallMethod(AmbientWeylModule, "for a sub Weyl module", true,
[IsSubWeylModule], 0,
function(S)
 # returns the ambient Weyl module
 return(S!.weylModule);
end );

InstallMethod(Weights, "for a sub Weyl module", true, 
[IsSubWeylModule], 0, 
function(S)
 # returns a list of the weight space labels of <S>
 return( DuplicateFreeList(List(BasisVecs(S), Weight)) );
end );

InstallMethod(WeightSpaces, "for a sub Weyl module", true,
[IsSubWeylModule], 0,
function(S)
 # returns a list of the weight spaces of <S>
 local wts,out,w;

 wts:=Weights(S);
 out:= [ ];
 for w in wts do
   Add(out, w);
   Add(out, WeightSpace(S,w));
 od;
 return(out);
end );

InstallMethod(Character, "for a sub Weyl module", true, 
[IsSubWeylModule], 0, 
function(S)
 # returns the character (a list of the weights and their multiplicities)
 # for the given sub Wey module S
 local ws, k, out;
 out:= [ ];
 ws:= WeightSpaces(S);
 for k in [2,4..Length(ws)] do
   Add(out, ws[k-1]);
   Add(out, Length(ws[k]));
 od;
 return(out);
end );

InstallMethod(DominantWeightSpaces, "for a sub Weyl module", true, 
[IsSubWeylModule], 0, 
function(S)
 # returns a list of the dominant weight spaces of S
 local out, k, ws;
 out:= [ ];
 ws:= WeightSpaces(S);
 for k in [2,4..Length(ws)] do
   if IsDominant( ws[k-1] ) then 
     Add(out, ws[k-1]);
     Add(out, ws[k]);
   fi;
 od;
 return(out);
end );

InstallMethod(WeightSpace, "for a sub Weyl module", true, 
[IsSubWeylModule,IsList], 0, 
function(S,wt)
 # returns a basis for the given weight space in <S> of weight <wt>
 local bb,out,b;
 bb:=BasisVecs(S); out:=[ ];
 for b in bb do
   if Weight(b) = wt then Add(out,b); fi;
 od;
 return(out);
end );




# The following is for testing purposes only, and is not documented.
#############################################################################
InstallMethod(SubWeylModule, 
"for a Weyl module, target dimensions, and a list of weight vectors", true, 
[IsWeylModule,IsPosInt,IsPosInt,IsList], 0, 
function(V,low,high,wtspace)
 # Returns a submodule of dimension not more than <dim>, 
 # generated by a linear combination of the given list of weight vectors.
 # Returns 'fail' if no such submodule is found.
 local zerovec, testcases,coeff,vec,p,sub;
 p:= V!.prime;
 zerovec:= 0*wtspace[1];
 testcases:= Tuples( [0..p-1], Length(wtspace) );
 for coeff in testcases do
   vec:= LinearCombination(wtspace, coeff); 
   if not vec = zerovec then 
      sub:= SubWeylModule(V,vec);
      if Length(sub) <= high and low <= Length(sub) then return(sub); fi; 
   fi;
 od;
 return( fail );
end );  


########################################################################
InstallOtherMethod(\^,
"for a UEA element and a highest wt module element",
true, 
[ IsUEALatticeElement, IsLeftAlgebraModuleElement ], 0, 
function(x,vec)
# apply a UEA element <x> to a vector <vec> in a highest weight module,
# or in a tensor product of highest weight modules
 local list,sum,i,n,a,coeff,applymon;

 applymon:= function(a,vec)
 # apply a single monomial <a> to <vec>
  local rep,rep2,L,UU,k,i,e,cb,N,n,posrtvec,negrtvec,cartanvec,zerov;
  zerov:= 0*vec;
  if vec = zerov then return vec; fi;
  rep:= ExtRepOfObj(vec);
  if IsWeightRepElement(rep) 
  then L:= FamilyObj(rep)!.algebra; 
  elif IsTensorElement(rep)
    then rep2:= ExtRepOfObj( ExtRepOfObj(rep)[1][1] );
         L:= FamilyObj(rep2)!.algebra; 
  fi;  
  cb:= ChevalleyBasis(L); 
  posrtvec:=cb[1]; negrtvec:=cb[2]; cartanvec:=cb[3];
  N:=Length(posrtvec); 
  n:=Length(a);
  while n>0 do 
   e:=a[n]; i:=a[n-1];
   if i>2*N then
      for k in [1..e] do vec:=cartanvec[i-2*N]^vec - (k-1)*vec; od;
      vec:= (1/Factorial(e))*vec;
    elif i>N then
      for k in [1..e] do vec:=posrtvec[i-N]^vec; od;
      vec:= (1/Factorial(e))*vec;
    else
      for k in [1..e] do vec:=negrtvec[i]^vec; od;
      vec:= (1/Factorial(e))*vec;
    fi;
   n:=n-2; 
   od;
  return vec;
  end;


 list:=ExtRepOfObj(x);
 i:=1; n:=Length(list); sum:=0*vec;
 while i<n do
  a:=list[i]; coeff:=list[i+1];
  sum:=sum+coeff*applymon(a,vec);
  i:=i+2;
  od;
 return sum;
 end );

################################################################# 
InstallOtherMethod(\mod,
"for a highest wt module element and a positive integer",
true, 
[ IsLeftAlgebraModuleElement, IsPosInt ], 0, 
function( elt, p )
# reduce coefficents of a highest weight or tensor product module element mod p
 local u,lu,lelt,k,res,out;

 u:= ExtRepOfObj( elt );
 lu:= ShallowCopy( u![1] );
 out:= [ ];

 for k in [2,4..Length(lu)] do
     res:= lu[k] mod p;
     if res <> 0 then 
        Add( out, lu[k-1] );
        Add( out, res );
     fi;
 od;


 if out = [ ] and IsTensorElement(u) then out:= [ [ ], 0]; fi;  
 lelt:= ObjByExtRep( FamilyObj( u ), out );
 return ObjByExtRep( FamilyObj( elt ), lelt );
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
InstallMethod(QuotientWeylModule, 
"for a sub Weyl module", true, [IsSubWeylModule], 0, 
function(S) 
 # creates the quotient V/S, where V is the ambient Weyl module
 local V,p,submodBasis, rowbasis,v,rowvec,
      VV,SS,Q,h,record,b,bb,cosetreps,vv,i,result;

 V:= S!.weylModule;
 if Dim(S) = 0 then
     record:= Objectify(NewType( FamilyObj(V), IsQuotientWeylModule ), 
             rec(WeylModule:=V,kernel:=S,quotient:=V,homomorphism:=[], 
             cosetReps:=BasisVecs(V),maximalVecs:=[],maximalVecsAmbiguous:=[]) 
      );
      return(record);
 fi;

 p:= V!.prime;
 submodBasis:= BasisVecs(S);
 rowbasis := S!.repbasis;
 SS:= VectorSpace(GF(p), rowbasis);
 VV:= GF(p)^Dim(V);
 Q:= VV/SS;
 h:= NaturalHomomorphismBySubspace(VV, SS);
 
 cosetreps:=[];
 b:= BasisVectors( Basis(Q) );
 bb:= BasisVecs(V);
 for v in b do
   vv:= PreImagesRepresentative(h, v);
   result:= 0*bb[1];
   for i in [1..Length(vv)] do
     result:= result + IntFFESymm(vv[i])*bb[i];
   od;
   Add(cosetreps, result);
 od;

 record:= Objectify(NewType( FamilyObj(V), IsQuotientWeylModule ), 
             rec(WeylModule:=V,kernel:=S,quotient:=Q,homomorphism:=h, 
             cosetReps:=cosetreps,maximalVecs:=[],maximalVecsAmbiguous:=[]) 
      );
 return(record);
end );

#############################################################################
InstallMethod(PrintObj, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(Q)
  local W;
  W:= Q!.WeylModule;
  Print("<Quotient of Type ",  W!.type, W!.rank, 
  " Weyl module of highest weight ", 
  W!.highestWeight, " at prime p = ", W!.prime, ">");
end );

#############################################################################
InstallMethod(AmbientWeylModule, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(Q)
 # returns the ambient Weyl module of <Q>
 return( Q!.WeylModule );
end );

#############################################################################
InstallMethod(DefiningKernel, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(Q)
 # returns the defining kernel of <Q>
 return( Q!.kernel );
end );

#############################################################################
InstallMethod(IsAmbiguous,  "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
function(W)
  if Length(W!.maximalVecsAmbiguous) > 0 then return true; fi;
end );

#############################################################################
InstallMethod(AmbiguousMaxVecs,  "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
function(W)
  return(W!.maximalVecsAmbiguous);
end );

#############################################################################
InstallMethod(TheLieAlgebra, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(W)
 # returns the underlying Lie algebra of <W>
 return( TheLieAlgebra(W!.WeylModule) );
end );

#############################################################################
InstallMethod(TheCharacteristic, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(W)
 # returns the characteristic of the base field
 return( TheCharacteristic(W!.WeylModule) );
end );

#############################################################################
InstallMethod(BasisVecs, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(W)
 # returns the standard basis of <W>
 return(W!.cosetReps);
end );

#############################################################################
InstallMethod(Generator, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(W)
 # returns the generator of <W>
 return BasisVecs(W)[1];
end );

#############################################################################
InstallMethod(Dim, "for quotient Weyl modules", true, 
[IsQuotientWeylModule], 0, 
function(W)
 # returns the dimension of <W>
 return Length(BasisVecs(W));
end );

#############################################################################
InstallMethod(Weights, "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
function(V)
 # returns a list of the weight space labels of <V>
 return( DuplicateFreeList(List(BasisVecs(V), Weight)) );
end );

#############################################################################
InstallMethod(DominantWeights, "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
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
InstallMethod(WeightSpaces, "for a quotient Weyl module", true,
[IsQuotientWeylModule], 0,
function(Q)
 # returns a list of the weight spaces of <Q>
 local wts,out,w;

 wts:=Weights(Q);
 out:= [ ];
 for w in wts do
   Add(out, w);
   Add(out, WeightSpace(Q,w));
 od;
 return(out);
end );

#############################################################################
InstallMethod(Character, "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
function(V)
 # returns a list of the weights and multiplicities for V
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
InstallMethod(DominantWeightSpaces, "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
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
InstallMethod(WeightSpace, "for a quotient Weyl module", true, 
[IsQuotientWeylModule,IsList], 0, 
function(Q,wt)
 # returns a basis for the given weight space in <Q> of weight <wt>
 local bb,out,b;
 bb:=BasisVecs(Q); out:=[ ];
 for b in bb do
   if Weight(b) = wt then Add(out,b); fi;
 od;
 return(out);
end );

#############################################################################
InstallMethod(ActOn, 
"for a quotient Weyl module, UEA element, and coset rep vector", true, 
[IsQuotientWeylModule,IsUEALatticeElement,IsLeftAlgebraModuleElement], 0, 
function(Q,x,v)
 # returns the result of acting by x on v
 local h,p,V,result,row,bb,vv,ans,i;

 V:= AmbientWeylModule(Q);
 h:= Q!.homomorphism;
 p:= V!.prime;
 result:= Image(h, RowVec(V, x^v mod p));
 
 bb:= BasisVecs(V);
 vv:= PreImagesRepresentative(h, result);
 ans:= 0*bb[1];
 for i in [1..Length(vv)] do
   ans:= ans + IntFFESymm(vv[i])*bb[i];
 od;
 return(ans);
end );

#############################################################################
# At the moment, operations depending on the following are not likely
# to work if we obtain more than one independent maximal vector in the
# given weight space. Such ambiguous vectors are stored in the module,
# and a warning message is issued in such a case.
#############################################################################
InstallMethod(MaximalVectors, "for a quotient Weyl module and a weight", true, 
[IsQuotientWeylModule,IsList], 0, 
function(Q,wt)
 # Returns a relative maximal vector of a given weight space in a 
 # quotient Weyl module, if possible - relative to the <submodule> spanned
 # by the given <basis>. In other words, it produces a primitive vector, 
 # if one exists relative to <submodule>, in the given weight space. 
 
 local rank,i,j,k,p,height,highestPrimePower,vec,rowvec,finalmatrix,wtspace,
       item,outlist,result,join,row,S,matrix,simpleGens,xy,xsimple,ysimple,
       z,V;

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

 V:= AmbientWeylModule(Q);
 wtspace:= WeightSpace(Q,wt);
 p:= V!.prime;
 xy:= SimpleLieAlgGens(V); xsimple:=xy[1]; ysimple:=xy[2];
 rank:= Length(xsimple);
 height:= HighestPrimePower(p, Sum(V!.highestWeight));
 finalmatrix:= [];
 for j in [1..rank] do 
   for k in [0..height] do
     matrix:= []; 
     for vec in wtspace do
       rowvec:= RowVec(V, ActOn(Q, (xsimple[j]^(p^k)/Factorial(p^k)), vec) );
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
       Print("***** WARNING: Ambiguous quotient module detected *****\n");
    fi;
 fi; 
 return(outlist);
end );  


#############################################################################
InstallMethod(MaximalVectors, "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
function(Q)
 # Finds a maximal vector in each weight space of given quotient Weyl 
 # module <Q>, if possible.
 
 local k,max,wts,out;
 if Q!.maximalVecs = [] then   
   wts:= DominantWeights(Q);
   out:= [ ];
   for k in [1..Length(wts)] do
     max:= MaximalVectors(Q,wts[k]);
     if Length(max) > 0 then
       Append(out, max);
     fi;
   od;
   Q!.maximalVecs:= out; #remember for next time
   return( out );
 else
   return Q!.maximalVecs; 
 fi;   
end );

#############################################################################
InstallMethod(NextSocle, "for a sub Weyl module", true,
[IsSubWeylModule], 0, 
function(S)
 # Returns the next socle of V (the submodule that maps onto soc(V/S))

 local V, Q,mvecs, m, sub, lam, p, t, r, ans;
 V:= AmbientWeylModule(S);
 ans:= S; # initially
 p:=V!.prime; t:=V!.type; r:=V!.rank;

 Q:= QuotientWeylModule(S);
 mvecs:= MaximalVectors(Q);
 for m in mvecs do
    sub:= SubWeylModule(Q, m);
    lam:= Weight(m);
    if SortedCharacter(Character(sub))
       = SortedCharacter(SimpleCharacter(p,lam,t,r)) then
       ans:= SubWeylModule(ans, m);
    fi;
 od;
 return ans;
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
InstallMethod(IsWithin, 
"for a sub quotient Weyl module and a given weight vector", true, 
[IsSubQuotientWeylModule,IsLeftAlgebraModuleElement], 0, 
function(S,vec)
 # Returns true if the given <vec> lies within the given submodule S
 local Q,h,V,rowbasis, v, rowvec, SS; 
 Q:=S!.quotientWeylModule;
 V:=AmbientWeylModule(Q);
 h:= Q!.homomorphism;

 rowbasis:=S!.repbasis;
 SS:= VectorSpace(GF(V!.prime), rowbasis);
 rowvec:= Image(h, RowVec(V,vec));
 if rowvec in SS then
    return true;
 else
    return false;
 fi;
end );  



#############################################################################
InstallMethod(SubWeylModule, "for a quotient Weyl module and vector", true, 
[IsQuotientWeylModule,IsLeftAlgebraModuleElement], 0, 
function(Q,vec)
 # finds a basis for the submodule of <Q> generated by <vec>
 local i,j,k,rowbasis,lbasis,newvec,rowvec,S,L,
       p,zerotensor,extendbasis,g,gens,noPosRoots,v,
       height,comps,wt,W,h,SM;

 W:= AmbientWeylModule(Q);
 if vec = Generator(W) then
    lbasis:= BasisVecs(Q);
    rowbasis:=[];
    for v in lbasis do
       Add(rowbasis, RowVec(W,v));
    od;   
    SM:=Objectify(NewType(FamilyObj(Q), IsSubQuotientWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=[vec],
     quotientWeylModule:=Q) );
    return(SM);
 fi;

 h:= Q!.homomorphism;
 p:= W!.prime;
 wt:= W!.highestWeight;
 height:= HighestPrimePower(p, Sum(wt));

 L:= W!.LieAlgebra;
 g:= LatticeGeneratorsInUEA( L );
 noPosRoots:= Length(ChevalleyBasis(L)[1]);
 gens:= g{ [1..2*noPosRoots] };
 for i in [1..height] do
    for j in [1..2*noPosRoots] do
       Add(gens, gens[j]^(p^i)/Factorial(p^i));
    od;
 od;

 rowbasis:= [ Image(h, RowVec(W,vec)) ]; lbasis:= [ vec ];
 S:= VectorSpace(GF(p), rowbasis);   

 extendbasis:= function(v)
   local i;
   for i in gens do
       newvec:= ActOn(Q, i, v);
       if newvec <> 0*v then  
         rowvec:= Image(h, RowVec(W,newvec));
         if not (rowvec in S) then
            Add( rowbasis, rowvec);
            S:= VectorSpace(GF(p), rowbasis);
            Add( lbasis, newvec );
         fi;
       fi;
   od;
 end;

 extendbasis(vec); 
 j:=2;
 while j <= Length(lbasis) do
    extendbasis(lbasis[j]); 
    j:= j+1;
 od;

 SM:=Objectify(NewType(FamilyObj(Q), IsSubQuotientWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=[vec],
     quotientWeylModule:=Q) );
 return(SM);
end );


#############################################################################
InstallMethod(SubWeylModule, "for a sub quotient and a vector", true, 
[IsSubQuotientWeylModule, IsLeftAlgebraModuleElement], 0, 
function(W,vec)
 # Returns the submodule of the ambient quotient Weyl module generated by the
 # given submodule <W> and given vector <vec>.

 local Q, p, U, rowbasis, lbasis, S, generators, ulbasis, urowbasis, i, SM;

 Q:=W!.quotientWeylModule; p:=AmbientWeylModule(Q)!.prime;
 lbasis:=ShallowCopy(W!.eltbasis);
 if Length(lbasis) = 0 then  # W is the zero module
    return( SubWeylModule(Q,vec) );
 fi;
 rowbasis:=ShallowCopy(W!.repbasis);
 generators:=ShallowCopy(W!.gens);
 Add(generators, vec);

 U:=SubWeylModule(Q,vec);
 ulbasis:= U!.eltbasis;
 urowbasis:=U!.repbasis;
 for i in [1..Length(ulbasis)] do
   S:=VectorSpace(GF(p), rowbasis);
   if not (urowbasis[i] in S) then
     Add(lbasis, ulbasis[i]);
     Add(rowbasis, urowbasis[i]);
   fi;
 od;

 SM:=Objectify(NewType(FamilyObj(Q), IsSubQuotientWeylModule),
     rec(eltbasis:=lbasis,repbasis:=rowbasis,gens:=generators,
     quotientWeylModule:=Q) );
 return(SM);
end );

#############################################################################
InstallMethod(SubWeylModule, 
"for a quotient Weyl module and list of vectors", true, 
[IsQuotientWeylModule, IsList], 0, 
function(Q,vecs)
 # Returns the sub quotient module generated by the given vectors

 local S,v,SM;
 if Length(vecs) = 0 then
    SM:=Objectify(NewType(FamilyObj(Q), IsSubQuotientWeylModule),
     rec(eltbasis:=[],repbasis:=[],gens:=[],
     quotientWeylModule:=Q) );
 return(SM);
 fi;
 
 S:= SubWeylModule(Q, vecs[1]);
 for v in vecs{ [2..Length(vecs)] } do
   S:= SubWeylModule(S,v);
 od;
 return( S );
end );

#############################################################################
InstallMethod(PrintObj, "for a sub quotient Weyl module", true,
[IsSubQuotientWeylModule], 0,
function(S)
 Print("<SubQuotient WeylModule of dimension ", Dim(S),
 ", generated by elements ", S!.gens, " of weights ",
 List(S!.gens, Weight), ">, in\n", S!.quotientWeylModule);
end );

InstallMethod(Generators, "for a sub quotient Weyl module", true,
[IsSubQuotientWeylModule], 0,
function(S)
 return(S!.gens);
end );

InstallMethod(BasisVecs, "for a sub quotient Weyl module", true,
[IsSubQuotientWeylModule], 0,
function(S)
 # returns the basis vectors of <S>
 return(S!.eltbasis);
end );

InstallMethod(Dim, "for a sub quotient Weyl module", true,
[IsSubQuotientWeylModule], 0,
function(S)
 # returns the dimension of <S>
 return(Length(S!.eltbasis));
end );

InstallMethod(AmbientQuotient, "for a sub quotient Weyl module", true,
[IsSubQuotientWeylModule], 0,
function(S)
 # returns the ambient quotient Weyl module
 return(S!.quotientWeylModule);
end );

InstallMethod(Weights, "for a sub quotient Weyl module", true, 
[IsSubQuotientWeylModule], 0, 
function(S)
 # returns a list of the weight space labels of <S>
 return( DuplicateFreeList(List(BasisVecs(S), Weight)) );
end );

InstallMethod(WeightSpaces, "for a sub quotient Weyl module", true,
[IsSubQuotientWeylModule], 0,
function(S)
 # returns a list of the weight spaces of <S>
 local wts,out,w;

 wts:=Weights(S);
 out:= [ ];
 for w in wts do
   Add(out, w);
   Add(out, WeightSpace(S,w));
 od;
 return(out);
end );

InstallMethod(Character, "for a sub quotient Weyl module", true, 
[IsSubQuotientWeylModule], 0, 
function(S)
 # returns the character (a list of the weights and their multiplicities)
 # for the given sub quotient Wey module S
 local ws, k, out;
 out:= [ ];
 ws:= WeightSpaces(S);
 for k in [2,4..Length(ws)] do
   Add(out, ws[k-1]);
   Add(out, Length(ws[k]));
 od;
 return(out);
end );

InstallMethod(DominantWeightSpaces, "for a sub quotient Weyl module", true, 
[IsSubQuotientWeylModule], 0, 
function(S)
 # returns a list of the dominant weight spaces of S
 local out, k, ws;
 out:= [ ];
 ws:= WeightSpaces(S);
 for k in [2,4..Length(ws)] do
   if IsDominant( ws[k-1] ) then 
     Add(out, ws[k-1]);
     Add(out, ws[k]);
   fi;
 od;
 return(out);
end );

InstallMethod(WeightSpace, "for a sub quotient Weyl module", true, 
[IsSubQuotientWeylModule,IsList], 0, 
function(S,wt)
 # returns a basis for the given weight space in <S> of weight <wt>
 local bb,out,b;
 bb:=BasisVecs(S); out:=[ ];
 for b in bb do
   if Weight(b) = wt then Add(out,b); fi;
 od;
 return(out);
end );


#############################################################################
InstallMethod(SimpleCharacter, 
"for a prime <p>, a highest weight <wt>, a type <t>, and a rank <r>",
true, [IsPosInt, IsList, IsString, IsPosInt], 0, 
function(p, wt, t, r) 
 local frob, product, rmdr, quot, char;
  
 frob:= function(lambda)
   # Frobenius twist of SimpleCharacter of highest weight lambda
   local ch, j;
   ch:= SimpleCharacter(p,lambda,t,r);
   for j in [2,4..Length(ch)] do ch[j-1]:= p*ch[j-1]; od;
   return(ch);
 end;

 product:= function(a,b)
   # multiplies two characters
   local i,j,k,w,mult,out,found;
   out:=[ ];
   for j in [2,4..Length(a)] do
     for k in [2,4..Length(b)] do
       w:= a[j-1]+b[k-1]; mult:= a[j]*b[k];
       # is w already in output list?
       found:= 0;
       for i in [1..Length(out)] do
         if out[i] = w then found:= i; break; fi;
       od;
       if found = 0 then
         Add(out,w); Add(out,mult); 
       else
         out[found+1]:= out[found+1] + mult;
       fi;
     od;
   od;
   return(out);
 end;

 if IsRestrictedWeight(p,wt) then 
   return( SimpleTopFactorCharacter( WeylModule(p,wt,t,r) ) );
 else
   # the weight is not restricted - apply Steinberg TP thm
   rmdr:= wt mod p; quot:= (wt - rmdr)/p; 
   char:= SimpleTopFactorCharacter( WeylModule(p,rmdr,t,r) );
   return( product(char, frob(quot)) ); 
 fi;
end);
 
#############################################################################
InstallMethod(SimpleCharacter, 
"for a WeylModule <V> and a highest weight <wt>",
true, [IsWeylModule, IsList], 0, 
function(V, wt)
 # Returns the simple char of h.w. <wt> of the same type as V  
 return SimpleCharacter(V!.prime, wt, V!.type, V!.rank);
end );

#############################################################################
InstallMethod(SocleWeyl, "for a quotient Weyl module", true, 
[IsQuotientWeylModule], 0, 
function(Q)
 # Returns a list of maximal vectors that generate the socle of <Q> 
 
 local outlist, v, mvecs, s, V, p, b, dima, dimb;
 
 V:= Q!.WeylModule; p:=V!.prime;
 outlist:= []; 
 mvecs:= MaximalVectors(Q);
 for v in mvecs do
     s:= SubWeylModule(Q,v); 
     dima:= Dim(s);
     b:= SimpleCharacter(p,Weight(v),V!.type,V!.rank); dimb:=CharacterDim(b);
     if dima = dimb then
       Add(outlist, s);
     fi;
 od;
 return SubWeylModuleDirectSum(outlist);
end );

#############################################################################
InstallMethod(SocleWeyl, "for a Weyl module", true, 
[IsWeylModule], 0, 
function(V)
 # Returns a list of maximal vectors that generate the socle of <V> 
 local mvecs,outlist,p,s,b,dima,dimb,v;
 p:= V!.prime;
 outlist:= []; 
 mvecs:= MaximalVectors(V);
 for v in mvecs do
     s:= SubWeylModule(V,v); 
     dima:= Dim(s);
     b:= SimpleCharacter(p,Weight(v),V!.type,V!.rank); dimb:=CharacterDim(b);
     if dima = dimb then
         Add(outlist, s);
     fi;
 od;
 return SubWeylModuleDirectSum(outlist);
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
InstallMethod(CompositionToWeight, "for a list", true, 
[IsList], 0, 
function(mu)
 # returns the weight corresponding to the given composition
 local diffs, k;
 diffs:= [];
 for k in [1..Length(mu)-1] do
    Add(diffs, mu[k] - mu[k+1]);
 od;
 return diffs;
end );

#############################################################################
InstallMethod(WeightToComposition, "for an integer and a list", true, 
[IsInt, IsList], 0, 
function(deg, mu)
 # returns the composition corresponding to the given weight, 
 # in the given degree <deg>, if possible.
 local ans, k, s, n;
 ans:= [];
 for k in [1..Length(mu)] do
    Add(ans, Sum( mu{[k..Length(mu)]} ));
 od;
 Add(ans, 0);  #extend the length
 s:= Sum(ans);
 n:= Length(ans);
 if (deg - s) mod n = 0 then
    ans:= (deg-s)/n + ans;
    if Minimum(ans) < 0 then
       return(fail);
    else
       return(ans);
    fi;
 else
    return(fail);
 fi;
end );

#############################################################################
InstallMethod(BoundedPartitions, "for integers n,r,s", true,
[IsInt, IsInt, IsInt], 0, 
function(n,r,s)
   # returns a list of n-part partitions of r where each
   # part lies in the interval [0,s]
   local k, m, ans, a, b;

   #base cases
   if n*s < r then return []; fi;
   if n*s = r then return [ List( [1,2..n], i -> s ) ]; fi;
   if r = 0 then return [ List( [1,2..n], i -> 0 ) ]; fi;
   if n = 1 then return [ [r] ]; fi;

   #general case
   ans:= [];
   m:= Maximum(0,r-s);
   for k in [m,m+1..r-1] do
      for a in BoundedPartitions(n-1,k,Minimum(s,r-k)) do
         b:= a;
         Add(b, r-k, 1);  #put r-k in front of a
         Add(ans, b);
         od;
      od;
   return ans;
end );

#############################################################################
InstallMethod(BoundedPartitions, "for integers n,r", true,
[IsInt, IsInt], 0, 
function(n,r)
   return BoundedPartitions(n, r, r);
end );

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
InstallMethod(pRestrictedPartitions,
"for a prime p and a positive integer n", true, [IsInt, IsInt], 0, 
function(p, n)
   # returns a list of the p-restricted partitions of n
   local ptns, out, mu, strip, wt;

   strip:= function(x)
      local c;
      c:= ShallowCopy(x);
      while c[Length(c)] = 0 do Remove(c); od;
      return(c);
   end;

   ptns:= BoundedPartitions(n,n);
   out:= [];
   for mu in ptns do
      wt:= CompositionToWeight(mu);
      if IsRestrictedWeight(p, wt) then 
         Add(out, strip(mu));
      fi;
   od;
   return out;
end );

#############################################################################
InstallMethod(AllPartitions,
"for a positive integer n", true, [IsInt], 0, 
function(n)
   # returns a list of the partitions of n
   local strip, ptns, mu, out;

   strip:= function(x)
      local c;
      c:= ShallowCopy(x);
      while c[Length(c)] = 0 do Remove(c); od;
      return(c);
   end;

   ptns:= BoundedPartitions(n,n);
   out:= [];
   for mu in ptns do
      Add(out, strip(mu));
   od;
   return out;
end );

#############################################################################
InstallMethod(ProductCharacter,
"for two characters", true, [IsList,IsList], 0, 
function(a,b)
   # returns the product of two characters
   local i,j,k,w,mult,out,found;
   out:=[ ];
   for j in [2,4..Length(a)] do
     for k in [2,4..Length(b)] do
       w:= a[j-1]+b[k-1]; mult:= a[j]*b[k];
       # is w already in output list?
       found:= 0;
       for i in [1..Length(out)] do
         if out[i] = w then found:= i; break; fi;
       od;
       if found = 0 then
         Add(out,w); Add(out,mult); 
       else
         out[found+1]:= out[found+1] + mult;
       fi;
     od;
   od;
   return(out);
end );


#############################################################################
InstallMethod(DecomposeCharacter,
"for a given character, characteristic, and root system", true, 
[IsList,IsPosInt,IsString,IsPosInt], 0, 
function(ch,p,type,rank)
 # returns the decomposition numbers of a character
 local m, wt, mult, ms, decnums,
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
    L:= SimpleLieAlgebra(type,rank,Rationals);
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
 while ch <> [] do 
     m:= maximalWeight(ch);
     wt:= m[1]; mult:= m[2];
     Add(decnums, wt);
     Add(decnums, mult);
     ms:= multiple(mult, SimpleCharacter(p,wt,type,rank));
     ch:= DifferenceCharacter(ch, ms); 
 od;
 return decnums;
end );     


#############################################################################
InstallMethod(Conjugate, "for a given partition", true, [IsList], 0, 
function(mu)
  # returns the conjugate of a partition <mu>
  local j,ans,len;
  
  len:= function(mu, k)
    #returns the largest index n such that mu[i]>=k, all i<=n
    local n;
    if mu[Length(mu)]>=k then return Length(mu); fi;
    n:= 1;
    while mu[n] >= k do n:= n+1; od;
    return n-1;
  end;

  if Length(mu) = 0 then return []; fi;
  ans:= [];
  for j in [1..mu[1]] do Add(ans, len(mu,j)); od;
  return ans;
end );

#############################################################################
InstallMethod(pRestricted, "for a positive integer and a given partition", 
true, [IsPosInt,IsList], 0, 
function(e,a)
  # returns true iff <a> is e-restricted
  local i,ere;
  ere:=true;
  if Length(a)>0 then
      if a[Length(a)]>=e then
            ere:=false;
      elif Length(a)>1 then
            for i in [1..Length(a)-1] do
                  if a[i]-a[i+1]>=e then ere:=false;fi;
            od;
      fi;
  fi;
  return ere;
end );

#############################################################################
InstallMethod(pRegular, "for a positive integer and a given partition", 
true, [IsPosInt,IsList], 0, 
function(e,a)
  # returns true iff <a> is e-regular
  return pRestricted(e,Conjugate(a));
end );


#############################################################################
InstallMethod(Mullineux, "for a positive integer and a given partition", 
true, [IsPosInt,IsList], 0, 
function(e,pa)
  local r,newr,prim,newpa,row,i;
  #applies the e-Mullineux map to pa
  if pRegular(e,pa) =false then
       Print("e-singular!\n");
  elif pa=[] or e=2 then
       return pa;
  else
       r:=Length(pa);
       prim:=0;
       newpa:=ShallowCopy(pa);
       row:=1;
       while row<=Length(pa) and newpa[row]>0 do
               newpa[row]:=newpa[row]-1;
               prim:=prim+1;
               if prim mod e=0 or (row<Length(newpa) and
                  newpa[row]<newpa[row+1]) then row:=row+1;fi;
       od;
       for i in [1..Length(pa)] do
               if newpa[i]=0 then Unbind(newpa[i]);fi;
       od;
       newpa:= Mullineux(e,newpa);
       if prim mod e=0 then newr:=prim-r; else newr:=prim+1-r;fi;
       row:=newr;
       if Length(newpa)<newr then
               for i in [Length(newpa)+1..newr] do
                       newpa[i]:=0;
               od;
       fi;
       while prim>0 do
               newpa[row]:=newpa[row]+1;
               prim:=prim-1;
               if prim mod e=0 or (row>1 and newpa[row]>newpa[row-1]) 
                  then row:=row-1;fi;
       od;
       return newpa;
  fi;
end );


#############################################################################
InstallMethod(pRegularPartitions, 
"for a positive integer and a given partition", true, [IsPosInt,IsPosInt], 0, 
function(p,n)
    local ptns, f;
    f:= function(mu)
        return Mullineux(p,Conjugate(mu));
    end;
    
    ptns:= pRestrictedPartitions(p,n);
    Apply(ptns, f);
    return ptns;
end );

#############################################################################
