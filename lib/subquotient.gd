DeclareCategory( "IsSubQuotientWeylModule", 
  CategoryCollections(IsLeftAlgebraModuleElement) );

DeclareOperation( "IsWithin", 
            [IsSubQuotientWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule", 
            [IsQuotientWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule", 
            [IsSubQuotientWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule", 
            [IsQuotientWeylModule,IsList]);

DeclareOperation("SubWeylModuleDirectSum", [IsQuotientWeylModule,IsList]);

DeclareOperation("Generators", [IsSubQuotientWeylModule]);

DeclareOperation("BasisVecs", [IsSubQuotientWeylModule]);

DeclareOperation("Dim", [IsSubQuotientWeylModule]);

DeclareOperation("AmbientQuotient", [IsSubQuotientWeylModule]);

DeclareOperation("Weights", [IsSubQuotientWeylModule]);

DeclareOperation("WeightSpaces", [IsSubQuotientWeylModule]);

DeclareOperation("Character", [IsSubQuotientWeylModule]);

DeclareOperation("DominantWeightSpaces", [IsSubQuotientWeylModule]);

DeclareOperation("WeightSpace", [IsSubQuotientWeylModule,IsList]);

DeclareOperation("SocleWeyl", [IsQuotientWeylModule]);

DeclareOperation("NextSocle", [IsSubQuotientWeylModule]);
