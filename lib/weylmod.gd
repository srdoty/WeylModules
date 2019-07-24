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

DeclareOperation( "TheLieAlgebra",  [IsWeylModule]);

DeclareOperation( "TheCharacteristic", [IsWeylModule]); 

DeclareOperation( "BasisVecs", [IsWeylModule]);

DeclareOperation( "Generator", [IsWeylModule]);

DeclareOperation( "Dim", [IsWeylModule]);

DeclareOperation( "Weight", [IsLeftAlgebraModuleElement]);

DeclareOperation( "Weights", [IsWeylModule]);

# The following command is not documented at this time
DeclareOperation( "IsDominant", [IsList]);

DeclareOperation( "DominantWeights", [IsWeylModule]);

# The following command is not documented at this time
DeclareOperation( "IsRestrictedWeight", [IsPosInt, IsList] );

DeclareOperation( "WeightSpaces", [IsWeylModule]);

DeclareOperation( "Character", [IsWeylModule]);

DeclareOperation( "Character", [IsList]);

DeclareOperation( "DifferenceCharacter", [IsList,IsList]);

# The following command is not documented at this time
DeclareOperation( "CharacterDim", [IsList]);

DeclareOperation( "DominantWeightSpaces", [IsWeylModule]);

DeclareOperation( "WeightSpace", [IsWeylModule,IsList]);

DeclareOperation( "IsWithin", [IsWeylModule,IsList,IsLeftAlgebraModuleElement]);

DeclareOperation( "IsMaximalVector", [IsWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation("IsMaximalVector", 
        [IsWeylModule,IsList,IsLeftAlgebraModuleElement]);

DeclareOperation( "MaximalVectors", [IsWeylModule,IsList]);

DeclareOperation( "MaximalVectors", [IsWeylModule]);

DeclareOperation( "SubWeylModule", [IsWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation("SubWeylModule", [IsWeylModule,IsList]);

# The following command is not documented at this time
DeclareOperation("SubWeylModule", [IsWeylModule,IsPosInt,IsPosInt,IsList]);

DeclareOperation("ActOn", 
     [IsWeylModule, IsUEALatticeElement, IsLeftAlgebraModuleElement]);

############################################################################

DeclareCategory( "IsQuotientWeylModule", 
        CategoryCollections(IsLeftAlgebraModuleElement) );

DeclareOperation( "QuotientWeylModule", [IsWeylModule, IsList] );

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

############################################################################

DeclareOperation("PrimitiveVectors", [IsWeylModule]);

DeclareOperation("MaximalSubmodule", [IsWeylModule]);

# The following command is not documented at this time
DeclareOperation("SimpleTopFactorCharacter", [IsWeylModule]);

# The following command is not documented at this time
DeclareOperation("SimpleTopFactorDim", [IsWeylModule]);

DeclareOperation( "IsWithin", 
            [IsQuotientWeylModule,IsList,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule", 
            [IsQuotientWeylModule,IsLeftAlgebraModuleElement]);

DeclareOperation( "SubWeylModule", 
            [IsQuotientWeylModule,IsList]);


DeclareOperation( "SimpleCharacter", [IsPosInt, IsList, IsString, IsPosInt] );

DeclareOperation( "SimpleCharacter", [IsWeylModule, IsList] );

DeclareOperation("SocleWeyl", [IsQuotientWeylModule]);

DeclareOperation("SocleWeyl", [IsWeylModule]);

DeclareOperation("ExtWeyl", [IsWeylModule,IsList]);

DeclareOperation("SubmoduleStructure", [IsWeylModule]);

DeclareOperation("DecompositionNumbers", [IsWeylModule]);



DeclareOperation("CompositionToWeight", [IsList]);

DeclareOperation("WeightToComposition", [IsInt, IsList]);

DeclareOperation("BoundedPartitions", [IsInt, IsInt, IsInt]);

DeclareOperation("BoundedPartitions", [IsInt, IsInt]);

DeclareCategory("IsSchurAlgebraWeylModule", IsWeylModule );

DeclareOperation("SchurAlgebraWeylModule", [IsInt, IsList]);


#
#note - the following operation has ALREADY been declared, so no need 
#to declare it again...

#DeclareOperation("DecompositionNumbers", [IsSchurAlgebraWeylModule]);

#doing so leads to a warning error, but everything still works.
#

DeclareOperation("SchurAlgebraDecompositionMatrix", [IsInt, IsInt, IsInt]);

DeclareOperation("SymmetricGroupDecompositionNumbers", [IsInt, IsList]);

DeclareOperation("SymmetricGroupDecompositionMatrix", [IsInt, IsInt]);

DeclareOperation("pRestrictedPartitions", [IsInt, IsInt]);

DeclareOperation("AllPartitions", [IsInt]);

DeclareOperation("ProductCharacter", [IsList,IsList]);

DeclareOperation("DecomposeCharacter", [IsList, IsPosInt,
    IsString, IsPosInt]);


DeclareOperation("Conjugate", [IsList]);

DeclareOperation("pRestricted", [IsPosInt, IsList]);

DeclareOperation("pRegular", [IsPosInt, IsList]);

DeclareOperation("Mullineux", [IsPosInt, IsList]);

DeclareOperation("pRegularPartitions", [IsPosInt, IsPosInt]);

