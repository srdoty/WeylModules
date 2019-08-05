DeclareCategory("IsSchurAlgebraWeylModule", IsWeylModule );

DeclareOperation("SchurAlgebraWeylModule", [IsInt, IsList]);


###########################################################################
#note - the following operation has ALREADY been declared, so no need 
#to declare it again...

#DeclareOperation("DecompositionNumbers", [IsSchurAlgebraWeylModule]);

#doing so leads to a warning error, but everything still works.
###########################################################################

DeclareOperation("SchurAlgebraDecompositionMatrix", [IsInt, IsInt, IsInt]);

DeclareOperation("SymmetricGroupDecompositionNumbers", [IsInt, IsList]);

DeclareOperation("SymmetricGroupDecompositionMatrix", [IsInt, IsInt]);
