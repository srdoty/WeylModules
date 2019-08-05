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
