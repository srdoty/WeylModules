DeclareCategory( "IsSubWeylModule", 
    CategoryCollections(IsLeftAlgebraModuleElement) );
    
DeclareOperation( "SubWeylModule", [IsWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule",
     [IsSubWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation("SubWeylModule", [IsWeylModule,IsList]);

DeclareOperation("SubWeylModuleDirectSum", [IsWeylModule,IsList]);

DeclareOperation( "IsAmbiguous",  [IsSubWeylModule]);

DeclareOperation( "AmbiguousMaxVecs",  [IsSubWeylModule]);

DeclareOperation("Generators", [IsSubWeylModule]);

DeclareOperation("BasisVecs", [IsSubWeylModule]);

DeclareOperation("Dim", [IsSubWeylModule]);

DeclareOperation("AmbientWeylModule", [IsSubWeylModule]);

DeclareOperation("Weights", [IsSubWeylModule]);

DeclareOperation("DominantWeights", [IsSubWeylModule]);

DeclareOperation("WeightSpaces", [IsSubWeylModule]);

DeclareOperation("Character", [IsSubWeylModule]);

DeclareOperation("DominantWeightSpaces", [IsSubWeylModule]);

DeclareOperation("WeightSpace", [IsSubWeylModule,IsList]);

# The following command is not documented at this time
DeclareOperation("SubWeylModule", [IsWeylModule,IsPosInt,IsPosInt,IsList]);

DeclareOperation("IsMaximalVector", 
        [IsSubWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "MaximalVectors", [IsSubWeylModule,IsList]);

DeclareOperation( "MaximalVectors", [IsSubWeylModule]);

DeclareOperation( "IsWithin", [IsSubWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation("SocleWeyl", [IsWeylModule]);

DeclareOperation("NextSocle", [IsSubWeylModule]);

DeclareOperation("GensNextSocle", [IsSubWeylModule]);

