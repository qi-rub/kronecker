TestInternals := module()
option package;
export Run, ModuleApply;
local Tests;

Tests := table();

##############################################################################
# List and vector utilities
##############################################################################
Tests[ListDot] := proc()
  return Kronecker:-ListDot([1,2,3],[a,b,c]) = a + 2*b + 3*c;
end proc;

Tests[ListPermute] := proc()
  return Kronecker:-ListPermute([2, 1], [a, b]) = [b, a]
     and Kronecker:-ListPermute([2, 3, 1], [a, b, c]) = [b, c, a];
end proc;

Tests[ListMakePrimitiveAndUnique] := proc()
  return Kronecker:-ListMakePrimitiveAndUnique([-4, 2, 8]) = [2, -1, -4]
     and Kronecker:-ListMakePrimitiveAndUnique([-4/7, 2/7, 8/7]) = [2, -1, -4];
end proc;

Tests[RemoveListFromList] := proc()
  return Kronecker:-RemoveListFromList([2, 4, 2, 1], [1, 2, 2, 4, 4, 4, 1]) = [4, 4, 1];
end proc;

Tests[IsZeroVector] := proc()
  return not Kronecker:-IsZeroVector([0,0,2])
     and not Kronecker:-IsZeroVector(Vector([0,0,2]))
     and Kronecker:-IsZeroVector([0,0])
     and Kronecker:-IsZeroVector(Vector([0,0]));
end proc;

Tests[IsSameVector] := proc()
  return not Kronecker:-IsSameVector([2,1], [3,1])
     and not Kronecker:-IsSameVector(Vector([2,1]), Vector([3,1]))
     and not Kronecker:-IsSameVector([2,1], Vector([3,1]))
     and Kronecker:-IsSameVector([2,1], [2,1])
     and Kronecker:-IsSameVector(Vector([2,1]), [2,1])
     and Kronecker:-IsSameVector(Vector([2,1]), Vector([2,1]));
end proc;

Tests[UniqueVectorList] := proc()
  return numelems(Kronecker:-UniqueVectorList([Vector([2, 1]), Vector([2, 1]), [2, 1]])) = 1;
end proc;

Tests[RemoveVectorsFromList] := proc()
  return map(convert, Kronecker:-RemoveVectorsFromList([Vector([1, 0, -1])], map(Vector, [[1, -1, 0], [1, 0, -1], [0, 1, -1]])), list) = [[1, -1, 0], [0, 1, -1]];
end proc;

Tests[IsVectorSubset] := proc()
  return Kronecker:-IsVectorSubset([Vector([1, 0, -1])], map(Vector, [[1, -1, 0], [1, 0, -1], [0, 1, -1]]))
     and not Kronecker:-IsVectorSubset([Vector([-1, 0, 1])], map(Vector, [[1, -1, 0], [1, 0, -1], [0, 1, -1]]));
end proc;

Tests[BasisVector] := proc()
  return Kronecker:-IsSameVector(Kronecker:-BasisVector(2, 5), [0,1,0,0,0]);
end proc;

Tests[PolarizeVectors] := proc()
  local B;
  B := map(Vector, [[1, -1, 0], [1, 0, -1], [0, 1, -1]]);
  return map(convert, Kronecker:-PolarizeVectors([Vector([1, 0, -1])], B), list) = [[1, 0, -1]]
     and map(convert, Kronecker:-PolarizeVectors([Vector([-1, 0, 1])], B), list) = [[1, 0, -1]];
end proc;


##############################################################################
# Partitions, Highest Weights, Roots
##############################################################################

Tests[IsPartition] := proc()
  return Kronecker:-IsPartition([3, 2, 2, 1])
     and Kronecker:-IsPartition([2, 1, 0])
     and not Kronecker:-IsPartition([2, 1, 0, -1])
     and not Kronecker:-IsPartition([1, 2, 0]);
end proc;

Tests[TrimPartition] := proc()
  return Kronecker:-TrimPartition([3, 2, 1, 0, 0, 0, 0, 0, 0]) = [3, 2, 1]
     and Kronecker:-TrimPartition([3, 2, 1, 0]) = [3, 2, 1]
     and Kronecker:-TrimPartition([3, 2, 1]) = [3, 2, 1];
end proc;

Tests[ParseParts] := proc()
  return [Kronecker:-ParseParts([[3, 2, 1], [3, 1, 1, 1], [4, 2]])] = [8, 3, [4, 2], [3, 2, 1, 0, 0, 0, 0, 0], [3, 1, 1, 1, 4, 2]];
end proc;

Tests[PositiveRoots] := proc()
  return map(convert, Kronecker:-PositiveRoots(3), list) = [[1, -1, 0], [1, 0, -1], [0, 1, -1]];
end proc;

Tests[PositiveRootsProduct] := proc()
  return map(convert, Kronecker:-PositiveRootsProduct([3, 2]), list) = [[1, -1, 0, 0, 0], [1, 0, -1, 0, 0], [0, 1, -1, 0, 0], [0, 0, 0, 1, -1]];
end proc;

Tests[RestrictionMatrix] := proc()
  local cols;
  cols := [LinearAlgebra:-Column(Kronecker:-RestrictionMatrix([3, 2]), [1..-1])];
  return map(convert, cols, list) = [[1, 0, 0, 1, 0], [0, 1, 0, 1, 0], [0, 0, 1, 1, 0], [1, 0, 0, 0, 1], [0, 1, 0, 0, 1], [0, 0, 1, 0, 1]];
end proc;

Tests[GLToSLMatrix] := proc()
  local M;
  M := Kronecker:-GLToSLMatrix([3]);
  return Kronecker:-IsSameVector(M . Vector([1,-1,0]), [1,0])
     and Kronecker:-IsSameVector(M . Vector([0,1,-1]), [0,1])
     and Kronecker:-IsSameVector(M . Vector([1,1,1]), [0,0]);
end proc;


##############################################################################
# Orlik-Solomon bases
##############################################################################
Tests[OS] := proc()
  local A4, CH4, OS4;
  A4 := [[1, 0, 0, 0], [1, 0, -1, 0], [0, 1, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, 0], [1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1], [0, 0, 0, 1]];
  CH4 := {[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, 0], [0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, 0], [1, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 1, 0], [1, 1, 1, 1]};
  OS4 := Kronecker:-OS(A4, CH4);
  return numelems(OS4) = factorial(4);
end proc;

#Tests[GoodParts] := proc()
#end proc;


##############################################################################
# Hyperplanes and Indices
##############################################################################
Tests[HyperplaneNormals] := proc()
  local ns, Psi_nodup;
  ns := [3, 2];
  Psi_nodup := Kronecker:-UniqueVectorList([seq(Kronecker:-GLToSLMatrix(ns) . Kronecker:-RestrictionMatrix(ns) . alpha, alpha=Kronecker:-PositiveRoots(mul(ns)))]);
  return map(convert, Kronecker:-HyperplaneNormals(ns, Psi_nodup, use_database = true), list)
       = map(convert, Kronecker:-HyperplaneNormals(ns, Psi_nodup, use_database = false), list);
end proc;

Tests[Indices] := proc()
  local ns, GLToSL, R, DeltaPlusK_SL, AlphasPolarized_nodup, RootLatticeMatrix;
  ns := [3, 2];
  GLToSL := Kronecker:-GLToSLMatrix(ns);
  R := Kronecker:-RestrictionMatrix(ns);
  DeltaPlusK_SL := [seq(GLToSL . alpha, alpha=Kronecker:-PositiveRootsProduct(ns))];
  AlphasPolarized_nodup := Kronecker:-UniqueVectorList(Kronecker:-RemoveVectorsFromList(DeltaPlusK_SL, [seq(GLToSL . R . Vector(alpha), alpha=Kronecker:-PositiveRoots(mul(ns)))], strict=false));
  RootLatticeMatrix := LinearAlgebra:-IdentityMatrix(add(ns) - numelems(ns));
  return [Kronecker:-Indices(mul(ns), false, ns, AlphasPolarized_nodup, RootLatticeMatrix, use_database = true)]
       = [Kronecker:-Indices(mul(ns), false, ns, AlphasPolarized_nodup, RootLatticeMatrix, use_database = false)];
end proc;

##############################################################################
# Weyl groups
##############################################################################
Tests[WeylGroups] := proc()
  return numelems(Kronecker:-WeylGroups(9, 3)[3]) = 504
     and numelems(Kronecker:-WeylGroups(9, 3, rect=true)[3]) = 84;
end proc;


##############################################################################
# Iterated residues
##############################################################################
Tests[VanishingOrderAtZero] := proc()
  return Kronecker:-VanishingOrderAtZero(12, x) = 0
     and Kronecker:-VanishingOrderAtZero(3*x, x) = 1
     and Kronecker:-VanishingOrderAtZero(3*x^7, x) = 7;
end proc;

Tests[IteratedResidueAtOne] := proc()
  return Kronecker:-IteratedResidueAtOne((1 + 3*(x[1] - 1)) / (x[1]-1)^2, 1, x) = 3;
end proc;


##############################################################################
# Kronecker coefficients
##############################################################################
Tests[Kronecker] := proc()
  return Kronecker([[3, 2, 1], [3, 2, 1], [6, 0]]) = 1;
#     and simplify(Kronecker([[5, 3, 2, 1], [6, 5], [6, 5]], k) - ((1/4)*k^2+5/8+(3/8)*(-1)^k+(1/2)*k)) = 0   # stretched
#     and Kronecker([[1, 1], [1, 1], [1, 1]], coset=0, k) = 1   # stretched on cosets
#     and Kronecker([[1, 1], [1, 1], [1, 1]], coset=1, k) = 0
#     and simplify(Kronecker([[288, 192, 174, 120, 30, 6], [343, 270, 197], [654, 156]], [lambda, mu, nu], forceperturb = true) - 1/5040 * mul(lambda[1] + lambda[2] + lambda[3] - nu[1] + i, i=1..7) * (lambda[1] + lambda[2] + lambda[4] + lambda[5] - mu[1] - mu[2] + 1)) = 0;  # symbolic
end proc;


##############################################################################
# Run all internal tests
##############################################################################
Run := proc()
  local tname;
  for tname in indices(Tests, nolist, indexorder) do
    printf("Testing %s...", tname);
    if Tests[tname]() then
      printf("OK\n");
    else
      error("TEST FAILED");
    end if;
  end do;
  printf("ALL %d TESTS SUCCEEDED.\n", numelems(Tests));
end proc;

ModuleApply := Run;

end module;


TestInternals();
