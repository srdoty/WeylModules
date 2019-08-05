DeclareCategory( "IsSubWeylModule", 
    CategoryCollections(IsLeftAlgebraModuleElement) );
    
DeclareOperation( "IsWithin", [IsSubWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation("IsMaximalVector", 
        [IsSubWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule", [IsWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule",
     [IsSubWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation("SubWeylModule", [IsWeylModule,IsList]);

DeclareOperation("SubWeylModuleDirectSum", [IsList]);

DeclareOperation("Generators", [IsSubWeylModule]);

DeclareOperation("BasisVecs", [IsSubWeylModule]);

DeclareOperation("Dim", [IsSubWeylModule]);

DeclareOperation("AmbientWeylModule", [IsSubWeylModule]);

DeclareOperation("Weights", [IsSubWeylModule]);

DeclareOperation("WeightSpaces", [IsSubWeylModule]);

DeclareOperation("Character", [IsSubWeylModule]);

DeclareOperation("DominantWeightSpaces", [IsSubWeylModule]);

DeclareOperation("WeightSpace", [IsSubWeylModule,IsList]);

# The following command is not documented at this time
DeclareOperation("SubWeylModule", [IsWeylModule,IsPosInt,IsPosInt,IsList]);

DeclareOperation("NextSocle", [IsSubWeylModule]);



