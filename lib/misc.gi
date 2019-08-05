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

