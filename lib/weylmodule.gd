#############################################################################
##
#W  weylmod.gd                   GAP package              S.R. Doty
##
##
#Y  Copyright (C)  2009,  S.R. Doty
##
##  This file contains the declaration of attributes, properties, and
##  operations for Weyl modules, simple characters, etc.
##
#############################################################################

DeclareCategory( "IsWeylModule", 
  CategoryCollections(IsLeftAlgebraModuleElement) );

DeclareOperation( "WeylModule", [IsPosInt, IsList, IsString, IsPosInt] );

DeclareOperation( "WeylModule", [IsWeylModule,IsList] );

DeclareOperation( "IsAmbiguous",  [IsWeylModule]);

DeclareOperation( "AmbiguousMaxVecs",  [IsWeylModule]);

DeclareOperation( "TheLieAlgebra",  [IsWeylModule]);

DeclareOperation( "TheCharacteristic", [IsWeylModule]); 

DeclareOperation( "BasisVecs", [IsWeylModule]);

DeclareOperation( "Generator", [IsWeylModule]);

DeclareOperation( "Dim", [IsWeylModule]);

DeclareOperation( "Weights", [IsWeylModule]);

DeclareOperation( "DominantWeights", [IsWeylModule]);

DeclareOperation( "WeightSpace", [IsWeylModule,IsList]);

DeclareOperation( "WeightSpaces", [IsWeylModule]);

DeclareOperation("ActOn", 
     [IsWeylModule, IsUEALatticeElement, IsLeftAlgebraModuleElement]);

DeclareOperation( "DominantWeightSpaces", [IsWeylModule]);

DeclareOperation( "Character", [IsWeylModule]);

DeclareOperation( "RowVec", [IsWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SimpleLieAlgGens", [IsWeylModule]);

DeclareOperation( "IsMaximalVector", [IsWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "MaximalVectors", [IsWeylModule,IsList]);

DeclareOperation( "MaximalVectors", [IsWeylModule]);

DeclareOperation("SocleSeries", [IsWeylModule]);

DeclareOperation("SimpleQuotient", [IsWeylModule]);

DeclareOperation("MaximalSubmodule", [IsWeylModule]);

# The following command is not documented at this time
DeclareOperation("SimpleTopFactorCharacter", [IsWeylModule]);

# The following command is not documented at this time
DeclareOperation("SimpleTopFactorDim", [IsWeylModule]);

DeclareOperation("DecompositionNumbers", [IsWeylModule]);



