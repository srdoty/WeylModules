# The following command is not documented at this time
DeclareOperation( "SortedCharacter", [IsList] );

DeclareOperation( "Weight", [IsLeftAlgebraModuleElement]);

# The following command is not documented at this time
DeclareOperation( "IsDominant", [IsList]);

DeclareOperation( "DifferenceCharacter", [IsList,IsList]);

# The following command is not documented at this time
DeclareOperation( "CharacterDim", [IsList]);

# The following command is not documented at this time
DeclareOperation( "IsRestrictedWeight", [IsPosInt, IsList] );

DeclareOperation("ProductCharacter", [IsList,IsList]);


DeclareOperation("DecomposeCharacter", [IsList, IsPosInt,
    IsString, IsPosInt]);


DeclareOperation( "SimpleCharacter", [IsPosInt, IsList, IsString, IsPosInt] );

DeclareOperation( "SimpleCharacter", [IsWeylModule, IsList] );
