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
InstallMethod(IsRestrictedWeight, 
"for a prime <p> and a dominant weight <wt>",
true, [IsPosInt, IsList], 0, 
function(p, wt)
 local k;
 for k in wt do
   if k >= p or k < 0 then return(false); fi;
 od;
 return(true);
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

