DeclareCategory( "IsQuotientWeylModule", 
  CategoryCollections(IsLeftAlgebraModuleElement) );

DeclareOperation( "QuotientWeylModule", [IsSubWeylModule] );

DeclareOperation( "AmbientWeylModule", [IsQuotientWeylModule]);

DeclareOperation( "DefiningKernel", [IsQuotientWeylModule]);

DeclareOperation( "IsAmbiguous",  [IsQuotientWeylModule]);

DeclareOperation( "AmbiguousMaxVecs",  [IsQuotientWeylModule]);

DeclareOperation( "TheLieAlgebra", [IsQuotientWeylModule]);

DeclareOperation( "TheCharacteristic", [IsQuotientWeylModule]);

DeclareOperation( "BasisVecs", [IsQuotientWeylModule]);

DeclareOperation( "Generator", [IsQuotientWeylModule]);

DeclareOperation( "Dim", [IsQuotientWeylModule]);

DeclareOperation( "Weights", [IsQuotientWeylModule]);

DeclareOperation( "DominantWeights", [IsQuotientWeylModule]);

DeclareOperation( "WeightSpaces", [IsQuotientWeylModule]);

DeclareOperation( "Character", [IsQuotientWeylModule]);

DeclareOperation( "DominantWeightSpaces", [IsQuotientWeylModule]);

DeclareOperation( "WeightSpace", [IsQuotientWeylModule,IsList]);

DeclareOperation( "ActOn", [IsQuotientWeylModule,IsUEALatticeElement,
           IsLeftAlgebraModuleElement]);

DeclareOperation( "MaximalVectors", [IsQuotientWeylModule,IsList]);

DeclareOperation( "MaximalVectors", [IsQuotientWeylModule]);

DeclareOperation( "SocleSeries", [IsQuotientWeylModule]);