Kronecker := module()
option package;
export VERSION;
export ListDot, ListPermute, RemoveListFromList, ListMakePrimitiveAndUnique, IsZeroVector, IsSameVector, UniqueVectorList, IsVectorSubset, RemoveVectorsFromList, BasisVector, PolarizeVectors;
export IsPartition, TrimPartition, ParseParts, PositiveRoots, PositiveRootsProduct, RestrictionMatrix, GLToSLMatrix;
export OS, GoodParts, AllowRect;
export HyperplaneNormals, Indices;
export WeylGroups;
export VanishingOrderAtZero, IteratedResidueAtOne;
export Kronecker, ModuleApply;

VERSION := "v1.0";

##############################################################################
# Initialization
##############################################################################
if not interface(quiet) then
  printf("Kronecker %s by Velleda Baldoni, Michèle Vergne, Michael Walter ", VERSION);
  printf("(see https://github.com/catch22/kronecker for more information).\n");
end if;


##############################################################################
# List and vector utilities
##############################################################################

# dot product of lists
ListDot := proc(X, Y, $)
  return add(X[i] * Y[i], i=1..numelems(X));
end proc;

# apply permutation to list
ListPermute := proc(w, X, $)
  return [seq(X[w[i]], i=1..numelems(X))];
end proc;

# make list primitive in standard lattice and force first nonzero entry to be positive to make it unique
ListMakePrimitiveAndUnique := proc(X, $)
  local Y;

  # make primitive
  Y := convert(X, list);
  Y := Y * ilcm(seq(map(denom, Y)));
  Y := Y / igcd(seq(Y));

  # force first nonzero entry to be positive
  if ListTools:-SelectFirst(x -> x <> 0, Y) < 0 then
    Y := -Y;
  end if;
  return Y;
end proc;

# remove A from B
RemoveListFromList := proc(A, B, { [strict, strict_] := true }, $)
  local B_remaining, a, i;
  B_remaining := B;
  for a in A do
    if member(a, B_remaining, 'i') then
      B_remaining := subsop(i=NULL, B_remaining);
    elif strict then
      error "B not a subset of A";
    end if;
  end do;
  return B_remaining;
end proc;

# zero vector?
IsZeroVector := proc(X, $)
  return {op(convert(X, list))} = {0};
end proc;

# same vectors?
IsSameVector := proc(X, Y, $)
  return {seq(X[i] - Y[i], i=1..numelems(X))} = {0};
end proc;

# make list of vectors unique
UniqueVectorList := proc(L, $)
  return [seq(Vector(psi), psi={seq(convert(psi, list), psi=L)})];
end proc;

# remove A from B
RemoveVectorsFromList := proc(A, B, { strict := true }, $)
  return map(Vector, RemoveListFromList(map(convert, A, list), map(convert, B, list), strict_ = strict));
end proc;

# does A contain a subset of the vectors in B?
IsVectorSubset := proc(A, B, $)
  return {seq(convert(a, list), a=A)} subset {seq(convert(b, list), b=B)};
end proc;

# i-th basis vector of R^n
BasisVector := proc(i, n, $)
  return Vector(subsop(i=1, [seq(0, j=1..n)]));
end proc;

# polarize list of vectors in A with respect to B
PolarizeVectors := proc(A, B)
  local L, Bs, X;
  L := [];
  Bs := map(convert, B, list);
  for X in A do
    if convert(X, list) in Bs then
      L := [op(L), X];
    elif convert(-X, list) in Bs then
      L := [op(L), -X];
    else
      error "Internal error: cannot polarization";
    end if;
  end do;
  return L;
end proc;


##############################################################################
# Partitions, Highest Weights, Roots
##############################################################################

# return whether given list is a partition
IsPartition := proc(part, $)
  local diffs;
  diffs := zip((x,y) -> x-y, part, [op(2..,part), 0]);
  return evalb(numelems(select(x -> x < 0, diffs)) = 0);
end proc;

# remove all nonzero parts from partition
TrimPartition := proc(part, $)
  return select(x -> x <> 0, part);
end proc;

ParseParts := proc(parts, $)
  local trimmed_first, degrees, M, nz, ns, lambda, mu;

  # verify that parts are partitions, and trim first
  if numelems(select(p -> not IsPartition(p), parts)) > 0 then
    error "Expected partitions.";
  end if;
  trimmed_first := TrimPartition(parts[1]);

  # degrees should match
  degrees := map(add, parts);
  if numelems({degrees[]}) > 1 then
    error "Number of boxes should be the same in each partition.";
   end if;

  # first partition too large?
  ns := map(numelems, parts[2..]);
  nz := numelems(trimmed_first);
  M := mul(ns);
  if nz > M then
    error "First partition is too large";
  end if;

  # set up highest weights for GL(M) and GL(n[1]) x ... x GL(n[s])
  lambda := [seq(trimmed_first), seq(0, i=1..M-nz)];
  mu := [seq(op(p), p = parts[2..])];

  return M, nz, ns, lambda, mu;
end proc;

# postitive roots of GL(n)
PositiveRoots := proc(n, $)
  return [seq(seq(BasisVector(i, n) - BasisVector(j, n), j=i+1..n), i=1..n-1)];
end proc;

# positive roots of product of GL(n)'s
PositiveRootsProduct := proc(ns, $)
  local roots, left, right, k, alpha;
  roots := [];
  for k from 1 to numelems(ns) do
    for alpha in PositiveRoots(ns[k]) do
      left := LinearAlgebra:-ZeroVector(add(ns[..k-1]));
      right := LinearAlgebra:-ZeroVector(add(ns[k+1..]));
      alpha := Vector([left, alpha, right]);
      roots := [op(roots), alpha];
    end do;
  end do;
  return roots;
end proc;

# matrix for restricting weights from GL(prod(ns)) to product of GL(n)'s
RestrictionMatrix := proc(ns, $)
  local M, LL, T, i, j, k, indices;
  M := Matrix(add(ns), mul(ns), storage=sparse);
  #LL := [seq([seq(1..n)], n = ns)];
  LL := ListTools:-Reverse([seq([seq(1..n)], n = ns)]);  # TODO: without the reverse, the code is much slower - why?
  T := combinat:-cartprod(LL);
  i := 1;
  while not T:-finished do
    #indices := T:-nextvalue();
    indices := ListTools:-Reverse(T:-nextvalue());  # TODO: without the reverse, the code is much slower - why?
    for k from 1 to numelems(ns) do
      j := add(ns[1..k-1]) + indices[k];
      M[j, i] := 1;
    end do;
    i := i + 1;
  end do;
  return M;
end proc;

# convert from GL(n) weights to basis of simple roots for SL(n)
GLToSLMatrix := proc(ns, $)
  local Ms, M, n, i, j;
  Ms := [];
  for n in ns do
    M := Matrix(n-1, n);
    for i from 1 to n-1 do
      for j from 1 to i do
        M[i, j] := 1 - i/n;
      end do;
      for j from i+1 to n do
        M[i, j] := -i/n;
      end do;
    end do;
    Ms := [op(Ms), M];
  end do;
  return LinearAlgebra:-DiagonalMatrix(Ms);
end proc;


##############################################################################
# Orlik-Solomon bases
##############################################################################

# Returns all Orlik-Solomon bases for Psi.
#   The list Psi_nodup is an ordered set of vectors and should be without redundancies.
#   The set A must contain the normal vectors of all linear hyperplanes spanned by vectors in Psi.
#   If the vector opt_v is specified then only bases adapted to (the tope defined by) v are returned.
#   If the integer opt_dim is specified then it should contain the dimension of the vector space spanned by the Psi.
OS := proc(Psi_nodup_, A_, { opt_v_ := NULL, opt_dim_ := NULL }, $)
  # determine dimension of vector space spanned by Psi
  local Psi_nodup, A, opt_v, dim, todos, i, bases, GrowOS;

  # convert to list
  Psi_nodup := map(convert, Psi_nodup_, list);
  A := map(convert, A_, list);

  # need to compute rank?
  dim := opt_dim_;
  if dim = NULL then
    dim := LinearAlgebra:-Rank(Matrix(Psi_nodup));
  else
    ASSERT(dim = LinearAlgebra:-Rank(Matrix(Psi_nodup)));
  end if;

  # if vector passed...
  opt_v := NULL;
  if opt_v_ <> NULL then
    opt_v := Vector(opt_v_);

    # check that vector defines a unique tope
    if numelems(select(X -> ListDot(opt_v, X) = 0, A)) > 0 then
      error "Vector does not determine a unique tope.";
    end if;

    # TODO: possible optimization would be to check that opt_v is in the convex cone spanned by the Psi (probably not an optimization in practice)
  end if;

  # grow partial OS basis by one
  GrowOS := proc(pbasis, Psi, A, $)
    local s, psi, new_pbasis, new_A, other_A, new_todos, new_Psi, X, x;
    s := numelems(pbasis);
    psi := Psi[1];
    new_pbasis := [op(pbasis), psi];
    # SLOW:
    #new_A := select(X -> X . psi = 0, A);
    #if opt_v <> NULL then
    #  other_A := select(X -> X . psi * X . opt_v > 0, A);
    #else
    #  other_A := select(X -> X . psi <> 0, A);
    #end if;
    new_A := [];
    other_A := [];
    for X in A do
      x := ListDot(X, psi);
      if x = 0 then
        new_A := [op(new_A), X];
      elif opt_v <> NULL then
        if x * ListDot(X, opt_v) > 0 then
          other_A := [op(other_A), X];
        end if;
      else
        other_A := [op(other_A), X];
      end if;
    end do;

    new_todos := [];
    for X in other_A do
      # SLOW:
      #new_Psi := select(phi -> ListDot(X, phi) = 0, Psi);
      new_Psi := [];
      for psi in Psi do
        if ListDot(X, psi) = 0 then
          new_Psi := [op(new_Psi), psi];
        end if;
      end do;

      if numelems(new_Psi) >= dim - (s + 1) then
        new_todos := [op(new_todos), [new_pbasis, new_Psi, new_A]];
      end if;
    end do;
    return new_todos;
  end proc;

  # each todo item is a triple (partial basis, candidate vectors, normal vectors)
  todos := [[[], Psi_nodup, A]];
  for i from 1 to dim do
    todos := [seq(seq(GrowOS(todo[])), todo in todos)];
  end do;
  bases := [seq(todo[1], todo in todos)];

  # check that we actually created bases
  ASSERT(numelems(select(b -> LinearAlgebra:-Rank(Matrix(b)) <> dim, bases)) = 0, "Internal error: OS did create a non-basis.");
  return bases;
end proc;

# allow rectangular optimiation for given (nz, ns)? the requirement is that the corresponding slice of the moment polytope still needs to be maximal dimensional.
AllowRect := proc(nz, ns)
  return (nz = 2 and ns = [2,2])
      or (nz = 3 and ns = [3,3])
      or (nz = 4 and ns = [3,3])
      or (nz = 4 and ns = [4,2])
      or (nz = 4 and ns = [4,3])
      or (nz = 4 and ns = [4,4])
      or (nz = 2 and ns = [2,2,2])
      or (nz = 2 and ns = [2,2,2,2]);
end proc;

# return partition tuples that determine a unique tope for given (nz, rect, ns)
GoodParts := proc(nz, rect, ns, $)
  # rectangular scenarios?
  if rect then
    if nz = 2 and ns = [2,2] then
      return [[1/2, 1/2], [1/2+1/21,1/2-1/21], [1/2+4/21,1/2-4/21]];
    elif nz = 3 and ns = [3,3] then
      return [[1/3, 1/3, 1/3], [1/3+11/13650, 1/3-11/13650+1/1300, 1/3-1/1300], [1/3+97/6825, 1/3-97/6825+27/2275, 1/3-27/2275]];
    elif nz = 4 and ns = [3,3] then
      return [[1/4, 1/4, 1/4, 1/4], [1/3+911/13650, 1/3-911/13650+34/975, 1/3-34/975], [1/3+487/6825, 1/3-487/6825+779/13650, 1/3-779/13650]];
    elif nz = 2 and ns = [2,2,2] then
      return [[1/2, 1/2], [1/2+7/9600,1/2-7/9600], [1/2+5/1920,1/2-5/1920], [1/2+7/1600, 1/2-7/1600]];
    elif nz = 2 and ns = [2,2,2,2] then
      return [[1/2, 1/2], [1/2+127/245760, 1/2-127/245760], [1/2+1927/245760,1/2-1927/245760], [1/2+6553/245760,1/2-6553/245760], [1/2+5461/122880,1/2-5461/122880]];
    elif nz = 4 and ns = [4,2] then
      return [[1/4, 1/4, 1/4, 1/4], [263/1024, 65/256, 251/1024, 125/512], [7623/10240, 2617/10240]];
    elif nz = 4 and ns = [4,3] then
      return [[1/4, 1/4, 1/4, 1/4], [263/1024, 65/256, 251/1024, 125/512], [6601/10240, 2417/10240, 611/5120]];
    elif nz = 4 and ns = [4,4] then
      return [[1/4, 1/4, 1/4, 1/4], [263/1024, 65/256, 251/1024, 125/512], [18803/30720, 6751/30720, 3607/30720, 1559/30720]];
    end if;
  end if;

  #elif nz <= 2 and ns = [2,2] then
  #  return [[200/729, 80/729], [455/1458, 35/486], [217/729, 7/81]];
  #elif nz <= 3 and ns = [2,2] then
  #  return [[9/31, 6/31, 3/62], [23/62, 5/31], [23/62, 5/31]];
  if nz <= 4 and ns = [2,2] then
    return [[45/142, 15/71, 15/142, 3/142], [61/142, 16/71], [61/142, 16/71]];
  #elif nz <= 3 and ns = [3,2] then
  #  return [[133/850, 43/850, 1/850], [709/4250, 87/2125, 1/2125], [753/4250, 66/2125]];
  elif nz <= 6 and ns = [3,2] then
    return [[23/753, 16/753, 14/753, 5/753, 2/753, 1/753], [472/11295, 298/11295, 29/2259], [79/1506, 43/1506]];
  elif nz <= 3 and ns = [3,3] then
    return [[5333/308990, 153/14045, 1129/308990], [5281/308990, 3147/308990, 140/30899], [542/30899, 3277/308990, 1131/308990]];
  elif nz <= 2 and ns = [2,2,2] then
    return [[349/2908, 17/727], [377/2908, 10/727], [381/2908, 9/727], [191/1454, 35/2908]];
  elif nz <= 2 and ns = [2,2,2,2] then
    return [[39060/879107, 12915/879107], [181335/3516428, 26565/3516428], [180873/3516428, 27027/3516428], [180675/3516428, 27225/3516428], [2695/52484, 27335/3516428]];
  else
    error "Unknown scenario: %1, %2, %3", nz, rect, ns;
  end if;
end proc;


##############################################################################
# Hyperplanes and Indices
##############################################################################

# compute hyperplane normals for given system of restricted roots
HyperplaneNormals := proc(ns, Psi_nodup, { use_database := true }, $)
  local dim, normals, S, M, ker, X;

  # check dimension of ambient space
  dim := add(ns) - numelems(ns);
  if dim <> numelems(Psi_nodup[1]) then
    error "Dimension mismatch.";
  end if;

  # one of the known cases?
  if use_database then
    if ns = [2,2] then
      return map(Vector, {[0, 1], [1, -1], [1, 0], [1, 1]});
    elif ns = [3,2] then
      return map(Vector, {[0, 0, 1], [0, 1, -1], [0, 1, 0], [0, 1, 1], [1, -2, -1], [1, -2, 1], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, 0, -1], [1, 0, 0], [1, 0, 1], [1, 1, -1], [1, 1, 1], [2, -1, -1], [2, -1, 1]});
    elif ns = [2,2,2] then
      return map(Vector, {[0, 0, 1], [0, 1, -1], [0, 1, 0], [0, 1, 1], [1, -2, -1], [1, -2, 1], [1, -1, -2], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, -1, 2], [1, 0, -1], [1, 0, 0], [1, 0, 1], [1, 1, -2], [1, 1, -1], [1, 1, 0], [1, 1, 1], [1, 1, 2], [1, 2, -1], [1, 2, 1], [2, -1, -1], [2, -1, 1], [2, 1, -1], [2, 1, 1]});
    elif ns = [3,3] then
      return map(Vector, {[0, 0, 0, 1], [0, 0, 1, -1], [0, 0, 1, 0], [0, 1, -2, 1], [0, 1, -1, -1], [0, 1, -1, 0], [0, 1, -1, 1], [0, 1, -1, 2], [0, 1, 0, -1], [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, -2], [0, 1, 1, -1], [0, 1, 1, 0], [0, 1, 1, 1], [0, 1, 2, -1], [1, -3, -2, 1], [1, -3, -1, -1], [1, -3, -1, 2], [1, -3, 1, -2], [1, -3, 1, 1], [1, -3, 2, -1], [1, -2, -3, 1], [1, -2, -3, 2], [1, -2, -2, -1], [1, -2, -2, 1], [1, -2, -2, 3], [1, -2, -1, -2], [1, -2, -1, -1], [1, -2, -1, 0], [1, -2, -1, 1], [1, -2, -1, 2], [1, -2, -1, 3], [1, -2, 0, -1], [1, -2, 0, 1], [1, -2, 1, -3], [1, -2, 1, -2], [1, -2, 1, -1], [1, -2, 1, 0], [1, -2, 1, 1], [1, -2, 1, 2], [1, -2, 2, -3], [1, -2, 2, -1], [1, -2, 2, 1], [1, -2, 3, -2], [1, -2, 3, -1], [1, -1, -2, 1], [1, -1, -1, -1], [1, -1, -1, 0], [1, -1, -1, 1], [1, -1, -1, 2], [1, -1, 0, -1], [1, -1, 0, 0], [1, -1, 0, 1], [1, -1, 1, -2], [1, -1, 1, -1], [1, -1, 1, 0], [1, -1, 1, 1], [1, -1, 2, -1], [1, 0, -2, 1], [1, 0, -1, -1], [1, 0, -1, 0], [1, 0, -1, 1], [1, 0, -1, 2], [1, 0, 0, -1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, -2], [1, 0, 1, -1], [1, 0, 1, 0], [1, 0, 1, 1], [1, 0, 2, -1], [1, 1, -3, 1], [1, 1, -3, 2], [1, 1, -2, -1], [1, 1, -2, 1], [1, 1, -2, 3], [1, 1, -1, -2], [1, 1, -1, -1], [1, 1, -1, 0], [1, 1, -1, 1], [1, 1, -1, 2], [1, 1, -1, 3], [1, 1, 0, -1], [1, 1, 0, 1], [1, 1, 1, -3], [1, 1, 1, -2], [1, 1, 1, -1], [1, 1, 1, 0], [1, 1, 1, 1], [1, 1, 1, 2], [1, 1, 2, -3], [1, 1, 2, -1], [1, 1, 2, 1], [1, 1, 3, -2], [1, 1, 3, -1], [1, 2, -2, 1], [1, 2, -1, -1], [1, 2, -1, 2], [1, 2, 1, -2], [1, 2, 1, 1], [1, 2, 2, -1], [2, -3, -2, 1], [2, -3, -1, -1], [2, -3, -1, 2], [2, -3, 1, -2], [2, -3, 1, 1], [2, -3, 2, -1], [2, -1, -3, 1], [2, -1, -3, 2], [2, -1, -2, -1], [2, -1, -2, 1], [2, -1, -2, 3], [2, -1, -1, -2], [2, -1, -1, -1], [2, -1, -1, 0], [2, -1, -1, 1], [2, -1, -1, 2], [2, -1, -1, 3], [2, -1, 0, -1], [2, -1, 0, 1], [2, -1, 1, -3], [2, -1, 1, -2], [2, -1, 1, -1], [2, -1, 1, 0], [2, -1, 1, 1], [2, -1, 1, 2], [2, -1, 2, -3], [2, -1, 2, -1], [2, -1, 2, 1], [2, -1, 3, -2], [2, -1, 3, -1], [2, 1, -2, 1], [2, 1, -1, -1], [2, 1, -1, 2], [2, 1, 1, -2], [2, 1, 1, 1], [2, 1, 2, -1], [3, -2, -2, 1], [3, -2, -1, -1], [3, -2, -1, 2], [3, -2, 1, -2], [3, -2, 1, 1], [3, -2, 2, -1], [3, -1, -2, 1], [3, -1, -1, -1], [3, -1, -1, 2], [3, -1, 1, -2], [3, -1, 1, 1], [3, -1, 2, -1]});
    elif ns = [3,2,2] then
      return map(Vector, {[0, 0, 0, 1], [0, 0, 1, -1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, -2, -1], [0, 1, -2, 1], [0, 1, -1, -2], [0, 1, -1, -1], [0, 1, -1, 0], [0, 1, -1, 1], [0, 1, -1, 2], [0, 1, 0, -1], [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, -2], [0, 1, 1, -1], [0, 1, 1, 0], [0, 1, 1, 1], [0, 1, 1, 2], [0, 1, 2, -1], [0, 1, 2, 1], [0, 2, -1, -1], [0, 2, -1, 1], [0, 2, 1, -1], [0, 2, 1, 1], [1, -4, -2, -1], [1, -4, -2, 1], [1, -4, -1, -2], [1, -4, -1, 2], [1, -4, 1, -2], [1, -4, 1, 2], [1, -4, 2, -1], [1, -4, 2, 1], [1, -3, -2, -1], [1, -3, -2, 1], [1, -3, -1, -2], [1, -3, -1, -1], [1, -3, -1, 1], [1, -3, -1, 2], [1, -3, 1, -2], [1, -3, 1, -1], [1, -3, 1, 1], [1, -3, 1, 2], [1, -3, 2, -1], [1, -3, 2, 1], [1, -2, -3, -2], [1, -2, -3, -1], [1, -2, -3, 1], [1, -2, -3, 2], [1, -2, -2, -3], [1, -2, -2, -1], [1, -2, -2, 1], [1, -2, -2, 3], [1, -2, -1, -3], [1, -2, -1, -2], [1, -2, -1, -1], [1, -2, -1, 0], [1, -2, -1, 1], [1, -2, -1, 2], [1, -2, -1, 3], [1, -2, 0, -1], [1, -2, 0, 1], [1, -2, 1, -3], [1, -2, 1, -2], [1, -2, 1, -1], [1, -2, 1, 0], [1, -2, 1, 1], [1, -2, 1, 2], [1, -2, 1, 3], [1, -2, 2, -3], [1, -2, 2, -1], [1, -2, 2, 1], [1, -2, 2, 3], [1, -2, 3, -2], [1, -2, 3, -1], [1, -2, 3, 1], [1, -2, 3, 2], [1, -1, -2, -1], [1, -1, -2, 1], [1, -1, -1, -2], [1, -1, -1, -1], [1, -1, -1, 0], [1, -1, -1, 1], [1, -1, -1, 2], [1, -1, 0, -1], [1, -1, 0, 0], [1, -1, 0, 1], [1, -1, 1, -2], [1, -1, 1, -1], [1, -1, 1, 0], [1, -1, 1, 1], [1, -1, 1, 2], [1, -1, 2, -1], [1, -1, 2, 1], [1, 0, -2, -1], [1, 0, -2, 1], [1, 0, -1, -2], [1, 0, -1, -1], [1, 0, -1, 0], [1, 0, -1, 1], [1, 0, -1, 2], [1, 0, 0, -1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, -2], [1, 0, 1, -1], [1, 0, 1, 0], [1, 0, 1, 1], [1, 0, 1, 2], [1, 0, 2, -1], [1, 0, 2, 1], [1, 1, -3, -2], [1, 1, -3, -1], [1, 1, -3, 1], [1, 1, -3, 2], [1, 1, -2, -3], [1, 1, -2, -1], [1, 1, -2, 1], [1, 1, -2, 3], [1, 1, -1, -3], [1, 1, -1, -2], [1, 1, -1, -1], [1, 1, -1, 0], [1, 1, -1, 1], [1, 1, -1, 2], [1, 1, -1, 3], [1, 1, 0, -1], [1, 1, 0, 1], [1, 1, 1, -3], [1, 1, 1, -2], [1, 1, 1, -1], [1, 1, 1, 0], [1, 1, 1, 1], [1, 1, 1, 2], [1, 1, 1, 3], [1, 1, 2, -3], [1, 1, 2, -1], [1, 1, 2, 1], [1, 1, 2, 3], [1, 1, 3, -2], [1, 1, 3, -1], [1, 1, 3, 1], [1, 1, 3, 2], [1, 2, -2, -1], [1, 2, -2, 1], [1, 2, -1, -2], [1, 2, -1, -1], [1, 2, -1, 1], [1, 2, -1, 2], [1, 2, 1, -2], [1, 2, 1, -1], [1, 2, 1, 1], [1, 2, 1, 2], [1, 2, 2, -1], [1, 2, 2, 1], [1, 3, -2, -1], [1, 3, -2, 1], [1, 3, -1, -2], [1, 3, -1, 2], [1, 3, 1, -2], [1, 3, 1, 2], [1, 3, 2, -1], [1, 3, 2, 1], [2, -4, -3, -1], [2, -4, -3, 1], [2, -4, -1, -3], [2, -4, -1, -1], [2, -4, -1, 1], [2, -4, -1, 3], [2, -4, 1, -3], [2, -4, 1, -1], [2, -4, 1, 1], [2, -4, 1, 3], [2, -4, 3, -1], [2, -4, 3, 1], [2, -3, -2, -1], [2, -3, -2, 1], [2, -3, -1, -2], [2, -3, -1, -1], [2, -3, -1, 1], [2, -3, -1, 2], [2, -3, 1, -2], [2, -3, 1, -1], [2, -3, 1, 1], [2, -3, 1, 2], [2, -3, 2, -1], [2, -3, 2, 1], [2, -2, -1, -1], [2, -2, -1, 1], [2, -2, 1, -1], [2, -2, 1, 1], [2, -1, -3, -2], [2, -1, -3, -1], [2, -1, -3, 1], [2, -1, -3, 2], [2, -1, -2, -3], [2, -1, -2, -1], [2, -1, -2, 1], [2, -1, -2, 3], [2, -1, -1, -3], [2, -1, -1, -2], [2, -1, -1, -1], [2, -1, -1, 0], [2, -1, -1, 1], [2, -1, -1, 2], [2, -1, -1, 3], [2, -1, 0, -1], [2, -1, 0, 1], [2, -1, 1, -3], [2, -1, 1, -2], [2, -1, 1, -1], [2, -1, 1, 0], [2, -1, 1, 1], [2, -1, 1, 2], [2, -1, 1, 3], [2, -1, 2, -3], [2, -1, 2, -1], [2, -1, 2, 1], [2, -1, 2, 3], [2, -1, 3, -2], [2, -1, 3, -1], [2, -1, 3, 1], [2, -1, 3, 2], [2, 0, -1, -1], [2, 0, -1, 1], [2, 0, 1, -1], [2, 0, 1, 1], [2, 1, -2, -1], [2, 1, -2, 1], [2, 1, -1, -2], [2, 1, -1, -1], [2, 1, -1, 1], [2, 1, -1, 2], [2, 1, 1, -2], [2, 1, 1, -1], [2, 1, 1, 1], [2, 1, 1, 2], [2, 1, 2, -1], [2, 1, 2, 1], [2, 2, -3, -1], [2, 2, -3, 1], [2, 2, -1, -3], [2, 2, -1, -1], [2, 2, -1, 1], [2, 2, -1, 3], [2, 2, 1, -3], [2, 2, 1, -1], [2, 2, 1, 1], [2, 2, 1, 3], [2, 2, 3, -1], [2, 2, 3, 1], [3, -4, -2, -1], [3, -4, -2, 1], [3, -4, -1, -2], [3, -4, -1, 2], [3, -4, 1, -2], [3, -4, 1, 2], [3, -4, 2, -1], [3, -4, 2, 1], [3, -2, -2, -1], [3, -2, -2, 1], [3, -2, -1, -2], [3, -2, -1, -1], [3, -2, -1, 1], [3, -2, -1, 2], [3, -2, 1, -2], [3, -2, 1, -1], [3, -2, 1, 1], [3, -2, 1, 2], [3, -2, 2, -1], [3, -2, 2, 1], [3, -1, -2, -1], [3, -1, -2, 1], [3, -1, -1, -2], [3, -1, -1, -1], [3, -1, -1, 1], [3, -1, -1, 2], [3, -1, 1, -2], [3, -1, 1, -1], [3, -1, 1, 1], [3, -1, 1, 2], [3, -1, 2, -1], [3, -1, 2, 1], [3, 1, -2, -1], [3, 1, -2, 1], [3, 1, -1, -2], [3, 1, -1, 2], [3, 1, 1, -2], [3, 1, 1, 2], [3, 1, 2, -1], [3, 1, 2, 1], [4, -3, -2, -1], [4, -3, -2, 1], [4, -3, -1, -2], [4, -3, -1, 2], [4, -3, 1, -2], [4, -3, 1, 2], [4, -3, 2, -1], [4, -3, 2, 1], [4, -2, -3, -1], [4, -2, -3, 1], [4, -2, -1, -3], [4, -2, -1, -1], [4, -2, -1, 1], [4, -2, -1, 3], [4, -2, 1, -3], [4, -2, 1, -1], [4, -2, 1, 1], [4, -2, 1, 3], [4, -2, 3, -1], [4, -2, 3, 1], [4, -1, -2, -1], [4, -1, -2, 1], [4, -1, -1, -2], [4, -1, -1, 2], [4, -1, 1, -2], [4, -1, 1, 2], [4, -1, 2, -1], [4, -1, 2, 1]});
    elif ns = [2,2,2,2] then
      return map(Vector, {[0, 0, 0, 1], [0, 0, 1, -1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, -2, -1], [0, 1, -2, 1], [0, 1, -1, -2], [0, 1, -1, -1], [0, 1, -1, 0], [0, 1, -1, 1], [0, 1, -1, 2], [0, 1, 0, -1], [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, -2], [0, 1, 1, -1], [0, 1, 1, 0], [0, 1, 1, 1], [0, 1, 1, 2], [0, 1, 2, -1], [0, 1, 2, 1], [0, 2, -1, -1], [0, 2, -1, 1], [0, 2, 1, -1], [0, 2, 1, 1], [1, -4, -3, -2], [1, -4, -3, 2], [1, -4, -2, -3], [1, -4, -2, -1], [1, -4, -2, 1], [1, -4, -2, 3], [1, -4, -1, -2], [1, -4, -1, 2], [1, -4, 1, -2], [1, -4, 1, 2], [1, -4, 2, -3], [1, -4, 2, -1], [1, -4, 2, 1], [1, -4, 2, 3], [1, -4, 3, -2], [1, -4, 3, 2], [1, -3, -4, -2], [1, -3, -4, 2], [1, -3, -2, -4], [1, -3, -2, -2], [1, -3, -2, -1], [1, -3, -2, 1], [1, -3, -2, 2], [1, -3, -2, 4], [1, -3, -1, -2], [1, -3, -1, -1], [1, -3, -1, 1], [1, -3, -1, 2], [1, -3, 1, -2], [1, -3, 1, -1], [1, -3, 1, 1], [1, -3, 1, 2], [1, -3, 2, -4], [1, -3, 2, -2], [1, -3, 2, -1], [1, -3, 2, 1], [1, -3, 2, 2], [1, -3, 2, 4], [1, -3, 4, -2], [1, -3, 4, 2], [1, -2, -4, -3], [1, -2, -4, -1], [1, -2, -4, 1], [1, -2, -4, 3], [1, -2, -3, -4], [1, -2, -3, -2], [1, -2, -3, -1], [1, -2, -3, 1], [1, -2, -3, 2], [1, -2, -3, 4], [1, -2, -2, -3], [1, -2, -2, -1], [1, -2, -2, 1], [1, -2, -2, 3], [1, -2, -1, -4], [1, -2, -1, -3], [1, -2, -1, -2], [1, -2, -1, -1], [1, -2, -1, 0], [1, -2, -1, 1], [1, -2, -1, 2], [1, -2, -1, 3], [1, -2, -1, 4], [1, -2, 0, -1], [1, -2, 0, 1], [1, -2, 1, -4], [1, -2, 1, -3], [1, -2, 1, -2], [1, -2, 1, -1], [1, -2, 1, 0], [1, -2, 1, 1], [1, -2, 1, 2], [1, -2, 1, 3], [1, -2, 1, 4], [1, -2, 2, -3], [1, -2, 2, -1], [1, -2, 2, 1], [1, -2, 2, 3], [1, -2, 3, -4], [1, -2, 3, -2], [1, -2, 3, -1], [1, -2, 3, 1], [1, -2, 3, 2], [1, -2, 3, 4], [1, -2, 4, -3], [1, -2, 4, -1], [1, -2, 4, 1], [1, -2, 4, 3], [1, -1, -4, -2], [1, -1, -4, 2], [1, -1, -3, -2], [1, -1, -3, -1], [1, -1, -3, 1], [1, -1, -3, 2], [1, -1, -2, -4], [1, -1, -2, -3], [1, -1, -2, -2], [1, -1, -2, -1], [1, -1, -2, 0], [1, -1, -2, 1], [1, -1, -2, 2], [1, -1, -2, 3], [1, -1, -2, 4], [1, -1, -1, -3], [1, -1, -1, -2], [1, -1, -1, -1], [1, -1, -1, 0], [1, -1, -1, 1], [1, -1, -1, 2], [1, -1, -1, 3], [1, -1, 0, -2], [1, -1, 0, -1], [1, -1, 0, 0], [1, -1, 0, 1], [1, -1, 0, 2], [1, -1, 1, -3], [1, -1, 1, -2], [1, -1, 1, -1], [1, -1, 1, 0], [1, -1, 1, 1], [1, -1, 1, 2], [1, -1, 1, 3], [1, -1, 2, -4], [1, -1, 2, -3], [1, -1, 2, -2], [1, -1, 2, -1], [1, -1, 2, 0], [1, -1, 2, 1], [1, -1, 2, 2], [1, -1, 2, 3], [1, -1, 2, 4], [1, -1, 3, -2], [1, -1, 3, -1], [1, -1, 3, 1], [1, -1, 3, 2], [1, -1, 4, -2], [1, -1, 4, 2], [1, 0, -2, -1], [1, 0, -2, 1], [1, 0, -1, -2], [1, 0, -1, -1], [1, 0, -1, 0], [1, 0, -1, 1], [1, 0, -1, 2], [1, 0, 0, -1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, -2], [1, 0, 1, -1], [1, 0, 1, 0], [1, 0, 1, 1], [1, 0, 1, 2], [1, 0, 2, -1], [1, 0, 2, 1], [1, 1, -4, -2], [1, 1, -4, 2], [1, 1, -3, -2], [1, 1, -3, -1], [1, 1, -3, 1], [1, 1, -3, 2], [1, 1, -2, -4], [1, 1, -2, -3], [1, 1, -2, -2], [1, 1, -2, -1], [1, 1, -2, 0], [1, 1, -2, 1], [1, 1, -2, 2], [1, 1, -2, 3], [1, 1, -2, 4], [1, 1, -1, -3], [1, 1, -1, -2], [1, 1, -1, -1], [1, 1, -1, 0], [1, 1, -1, 1], [1, 1, -1, 2], [1, 1, -1, 3], [1, 1, 0, -2], [1, 1, 0, -1], [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 0, 2], [1, 1, 1, -3], [1, 1, 1, -2], [1, 1, 1, -1], [1, 1, 1, 0], [1, 1, 1, 1], [1, 1, 1, 2], [1, 1, 1, 3], [1, 1, 2, -4], [1, 1, 2, -3], [1, 1, 2, -2], [1, 1, 2, -1], [1, 1, 2, 0], [1, 1, 2, 1], [1, 1, 2, 2], [1, 1, 2, 3], [1, 1, 2, 4], [1, 1, 3, -2], [1, 1, 3, -1], [1, 1, 3, 1], [1, 1, 3, 2], [1, 1, 4, -2], [1, 1, 4, 2], [1, 2, -4, -3], [1, 2, -4, -1], [1, 2, -4, 1], [1, 2, -4, 3], [1, 2, -3, -4], [1, 2, -3, -2], [1, 2, -3, -1], [1, 2, -3, 1], [1, 2, -3, 2], [1, 2, -3, 4], [1, 2, -2, -3], [1, 2, -2, -1], [1, 2, -2, 1], [1, 2, -2, 3], [1, 2, -1, -4], [1, 2, -1, -3], [1, 2, -1, -2], [1, 2, -1, -1], [1, 2, -1, 0], [1, 2, -1, 1], [1, 2, -1, 2], [1, 2, -1, 3], [1, 2, -1, 4], [1, 2, 0, -1], [1, 2, 0, 1], [1, 2, 1, -4], [1, 2, 1, -3], [1, 2, 1, -2], [1, 2, 1, -1], [1, 2, 1, 0], [1, 2, 1, 1], [1, 2, 1, 2], [1, 2, 1, 3], [1, 2, 1, 4], [1, 2, 2, -3], [1, 2, 2, -1], [1, 2, 2, 1], [1, 2, 2, 3], [1, 2, 3, -4], [1, 2, 3, -2], [1, 2, 3, -1], [1, 2, 3, 1], [1, 2, 3, 2], [1, 2, 3, 4], [1, 2, 4, -3], [1, 2, 4, -1], [1, 2, 4, 1], [1, 2, 4, 3], [1, 3, -4, -2], [1, 3, -4, 2], [1, 3, -2, -4], [1, 3, -2, -2], [1, 3, -2, -1], [1, 3, -2, 1], [1, 3, -2, 2], [1, 3, -2, 4], [1, 3, -1, -2], [1, 3, -1, -1], [1, 3, -1, 1], [1, 3, -1, 2], [1, 3, 1, -2], [1, 3, 1, -1], [1, 3, 1, 1], [1, 3, 1, 2], [1, 3, 2, -4], [1, 3, 2, -2], [1, 3, 2, -1], [1, 3, 2, 1], [1, 3, 2, 2], [1, 3, 2, 4], [1, 3, 4, -2], [1, 3, 4, 2], [1, 4, -3, -2], [1, 4, -3, 2], [1, 4, -2, -3], [1, 4, -2, -1], [1, 4, -2, 1], [1, 4, -2, 3], [1, 4, -1, -2], [1, 4, -1, 2], [1, 4, 1, -2], [1, 4, 1, 2], [1, 4, 2, -3], [1, 4, 2, -1], [1, 4, 2, 1], [1, 4, 2, 3], [1, 4, 3, -2], [1, 4, 3, 2], [2, -4, -3, -1], [2, -4, -3, 1], [2, -4, -1, -3], [2, -4, -1, -1], [2, -4, -1, 1], [2, -4, -1, 3], [2, -4, 1, -3], [2, -4, 1, -1], [2, -4, 1, 1], [2, -4, 1, 3], [2, -4, 3, -1], [2, -4, 3, 1], [2, -3, -4, -1], [2, -3, -4, 1], [2, -3, -2, -1], [2, -3, -2, 1], [2, -3, -1, -4], [2, -3, -1, -2], [2, -3, -1, -1], [2, -3, -1, 1], [2, -3, -1, 2], [2, -3, -1, 4], [2, -3, 1, -4], [2, -3, 1, -2], [2, -3, 1, -1], [2, -3, 1, 1], [2, -3, 1, 2], [2, -3, 1, 4], [2, -3, 2, -1], [2, -3, 2, 1], [2, -3, 4, -1], [2, -3, 4, 1], [2, -2, -3, -1], [2, -2, -3, 1], [2, -2, -1, -3], [2, -2, -1, -1], [2, -2, -1, 1], [2, -2, -1, 3], [2, -2, 1, -3], [2, -2, 1, -1], [2, -2, 1, 1], [2, -2, 1, 3], [2, -2, 3, -1], [2, -2, 3, 1], [2, -1, -4, -3], [2, -1, -4, -1], [2, -1, -4, 1], [2, -1, -4, 3], [2, -1, -3, -4], [2, -1, -3, -2], [2, -1, -3, -1], [2, -1, -3, 1], [2, -1, -3, 2], [2, -1, -3, 4], [2, -1, -2, -3], [2, -1, -2, -1], [2, -1, -2, 1], [2, -1, -2, 3], [2, -1, -1, -4], [2, -1, -1, -3], [2, -1, -1, -2], [2, -1, -1, -1], [2, -1, -1, 0], [2, -1, -1, 1], [2, -1, -1, 2], [2, -1, -1, 3], [2, -1, -1, 4], [2, -1, 0, -1], [2, -1, 0, 1], [2, -1, 1, -4], [2, -1, 1, -3], [2, -1, 1, -2], [2, -1, 1, -1], [2, -1, 1, 0], [2, -1, 1, 1], [2, -1, 1, 2], [2, -1, 1, 3], [2, -1, 1, 4], [2, -1, 2, -3], [2, -1, 2, -1], [2, -1, 2, 1], [2, -1, 2, 3], [2, -1, 3, -4], [2, -1, 3, -2], [2, -1, 3, -1], [2, -1, 3, 1], [2, -1, 3, 2], [2, -1, 3, 4], [2, -1, 4, -3], [2, -1, 4, -1], [2, -1, 4, 1], [2, -1, 4, 3], [2, 0, -1, -1], [2, 0, -1, 1], [2, 0, 1, -1], [2, 0, 1, 1], [2, 1, -4, -3], [2, 1, -4, -1], [2, 1, -4, 1], [2, 1, -4, 3], [2, 1, -3, -4], [2, 1, -3, -2], [2, 1, -3, -1], [2, 1, -3, 1], [2, 1, -3, 2], [2, 1, -3, 4], [2, 1, -2, -3], [2, 1, -2, -1], [2, 1, -2, 1], [2, 1, -2, 3], [2, 1, -1, -4], [2, 1, -1, -3], [2, 1, -1, -2], [2, 1, -1, -1], [2, 1, -1, 0], [2, 1, -1, 1], [2, 1, -1, 2], [2, 1, -1, 3], [2, 1, -1, 4], [2, 1, 0, -1], [2, 1, 0, 1], [2, 1, 1, -4], [2, 1, 1, -3], [2, 1, 1, -2], [2, 1, 1, -1], [2, 1, 1, 0], [2, 1, 1, 1], [2, 1, 1, 2], [2, 1, 1, 3], [2, 1, 1, 4], [2, 1, 2, -3], [2, 1, 2, -1], [2, 1, 2, 1], [2, 1, 2, 3], [2, 1, 3, -4], [2, 1, 3, -2], [2, 1, 3, -1], [2, 1, 3, 1], [2, 1, 3, 2], [2, 1, 3, 4], [2, 1, 4, -3], [2, 1, 4, -1], [2, 1, 4, 1], [2, 1, 4, 3], [2, 2, -3, -1], [2, 2, -3, 1], [2, 2, -1, -3], [2, 2, -1, -1], [2, 2, -1, 1], [2, 2, -1, 3], [2, 2, 1, -3], [2, 2, 1, -1], [2, 2, 1, 1], [2, 2, 1, 3], [2, 2, 3, -1], [2, 2, 3, 1], [2, 3, -4, -1], [2, 3, -4, 1], [2, 3, -2, -1], [2, 3, -2, 1], [2, 3, -1, -4], [2, 3, -1, -2], [2, 3, -1, -1], [2, 3, -1, 1], [2, 3, -1, 2], [2, 3, -1, 4], [2, 3, 1, -4], [2, 3, 1, -2], [2, 3, 1, -1], [2, 3, 1, 1], [2, 3, 1, 2], [2, 3, 1, 4], [2, 3, 2, -1], [2, 3, 2, 1], [2, 3, 4, -1], [2, 3, 4, 1], [2, 4, -3, -1], [2, 4, -3, 1], [2, 4, -1, -3], [2, 4, -1, -1], [2, 4, -1, 1], [2, 4, -1, 3], [2, 4, 1, -3], [2, 4, 1, -1], [2, 4, 1, 1], [2, 4, 1, 3], [2, 4, 3, -1], [2, 4, 3, 1], [3, -4, -2, -1], [3, -4, -2, 1], [3, -4, -1, -2], [3, -4, -1, 2], [3, -4, 1, -2], [3, -4, 1, 2], [3, -4, 2, -1], [3, -4, 2, 1], [3, -2, -4, -1], [3, -2, -4, 1], [3, -2, -2, -1], [3, -2, -2, 1], [3, -2, -1, -4], [3, -2, -1, -2], [3, -2, -1, -1], [3, -2, -1, 1], [3, -2, -1, 2], [3, -2, -1, 4], [3, -2, 1, -4], [3, -2, 1, -2], [3, -2, 1, -1], [3, -2, 1, 1], [3, -2, 1, 2], [3, -2, 1, 4], [3, -2, 2, -1], [3, -2, 2, 1], [3, -2, 4, -1], [3, -2, 4, 1], [3, -1, -4, -2], [3, -1, -4, 2], [3, -1, -2, -4], [3, -1, -2, -2], [3, -1, -2, -1], [3, -1, -2, 1], [3, -1, -2, 2], [3, -1, -2, 4], [3, -1, -1, -2], [3, -1, -1, -1], [3, -1, -1, 1], [3, -1, -1, 2], [3, -1, 1, -2], [3, -1, 1, -1], [3, -1, 1, 1], [3, -1, 1, 2], [3, -1, 2, -4], [3, -1, 2, -2], [3, -1, 2, -1], [3, -1, 2, 1], [3, -1, 2, 2], [3, -1, 2, 4], [3, -1, 4, -2], [3, -1, 4, 2], [3, 1, -4, -2], [3, 1, -4, 2], [3, 1, -2, -4], [3, 1, -2, -2], [3, 1, -2, -1], [3, 1, -2, 1], [3, 1, -2, 2], [3, 1, -2, 4], [3, 1, -1, -2], [3, 1, -1, -1], [3, 1, -1, 1], [3, 1, -1, 2], [3, 1, 1, -2], [3, 1, 1, -1], [3, 1, 1, 1], [3, 1, 1, 2], [3, 1, 2, -4], [3, 1, 2, -2], [3, 1, 2, -1], [3, 1, 2, 1], [3, 1, 2, 2], [3, 1, 2, 4], [3, 1, 4, -2], [3, 1, 4, 2], [3, 2, -4, -1], [3, 2, -4, 1], [3, 2, -2, -1], [3, 2, -2, 1], [3, 2, -1, -4], [3, 2, -1, -2], [3, 2, -1, -1], [3, 2, -1, 1], [3, 2, -1, 2], [3, 2, -1, 4], [3, 2, 1, -4], [3, 2, 1, -2], [3, 2, 1, -1], [3, 2, 1, 1], [3, 2, 1, 2], [3, 2, 1, 4], [3, 2, 2, -1], [3, 2, 2, 1], [3, 2, 4, -1], [3, 2, 4, 1], [3, 4, -2, -1], [3, 4, -2, 1], [3, 4, -1, -2], [3, 4, -1, 2], [3, 4, 1, -2], [3, 4, 1, 2], [3, 4, 2, -1], [3, 4, 2, 1], [4, -3, -2, -1], [4, -3, -2, 1], [4, -3, -1, -2], [4, -3, -1, 2], [4, -3, 1, -2], [4, -3, 1, 2], [4, -3, 2, -1], [4, -3, 2, 1], [4, -2, -3, -1], [4, -2, -3, 1], [4, -2, -1, -3], [4, -2, -1, -1], [4, -2, -1, 1], [4, -2, -1, 3], [4, -2, 1, -3], [4, -2, 1, -1], [4, -2, 1, 1], [4, -2, 1, 3], [4, -2, 3, -1], [4, -2, 3, 1], [4, -1, -3, -2], [4, -1, -3, 2], [4, -1, -2, -3], [4, -1, -2, -1], [4, -1, -2, 1], [4, -1, -2, 3], [4, -1, -1, -2], [4, -1, -1, 2], [4, -1, 1, -2], [4, -1, 1, 2], [4, -1, 2, -3], [4, -1, 2, -1], [4, -1, 2, 1], [4, -1, 2, 3], [4, -1, 3, -2], [4, -1, 3, 2], [4, 1, -3, -2], [4, 1, -3, 2], [4, 1, -2, -3], [4, 1, -2, -1], [4, 1, -2, 1], [4, 1, -2, 3], [4, 1, -1, -2], [4, 1, -1, 2], [4, 1, 1, -2], [4, 1, 1, 2], [4, 1, 2, -3], [4, 1, 2, -1], [4, 1, 2, 1], [4, 1, 2, 3], [4, 1, 3, -2], [4, 1, 3, 2], [4, 2, -3, -1], [4, 2, -3, 1], [4, 2, -1, -3], [4, 2, -1, -1], [4, 2, -1, 1], [4, 2, -1, 3], [4, 2, 1, -3], [4, 2, 1, -1], [4, 2, 1, 1], [4, 2, 1, 3], [4, 2, 3, -1], [4, 2, 3, 1], [4, 3, -2, -1], [4, 3, -2, 1], [4, 3, -1, -2], [4, 3, -1, 2], [4, 3, 1, -2], [4, 3, 1, 2], [4, 3, 2, -1], [4, 3, 2, 1]});
    elif ns = [4,2] then
      return map(Vector, {[0, 0, 0, 1], [0, 0, 1, -1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, -2, -1], [0, 1, -2, 1], [0, 1, -1, -1], [0, 1, -1, 0], [0, 1, -1, 1], [0, 1, 0, -1], [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, -1], [0, 1, 1, 1], [0, 2, -1, -1], [0, 2, -1, 1], [1, -3, 1, -1], [1, -3, 1, 1], [1, -2, -1, -1], [1, -2, -1, 1], [1, -2, 0, -1], [1, -2, 0, 1], [1, -2, 1, -1], [1, -2, 1, 1], [1, -2, 2, -1], [1, -2, 2, 1], [1, -2, 3, -1], [1, -2, 3, 1], [1, -1, -1, -1], [1, -1, -1, 1], [1, -1, 0, -1], [1, -1, 0, 0], [1, -1, 0, 1], [1, -1, 1, -1], [1, -1, 1, 0], [1, -1, 1, 1], [1, -1, 2, -1], [1, -1, 2, 1], [1, 0, -2, -1], [1, 0, -2, 1], [1, 0, -1, -1], [1, 0, -1, 0], [1, 0, -1, 1], [1, 0, 0, -1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, -1], [1, 0, 1, 1], [1, 1, -3, -1], [1, 1, -3, 1], [1, 1, -2, -1], [1, 1, -2, 1], [1, 1, -1, -1], [1, 1, -1, 1], [1, 1, 0, -1], [1, 1, 0, 1], [1, 1, 1, -1], [1, 1, 1, 1], [1, 2, -1, -1], [1, 2, -1, 1], [2, -3, 2, -1], [2, -3, 2, 1], [2, -2, 1, -1], [2, -2, 1, 1], [2, -1, -2, -1], [2, -1, -2, 1], [2, -1, -1, -1], [2, -1, -1, 1], [2, -1, 0, -1], [2, -1, 0, 1], [2, -1, 1, -1], [2, -1, 1, 1], [2, -1, 2, -1], [2, -1, 2, 1], [2, 0, -1, -1], [2, 0, -1, 1], [2, 1, -2, -1], [2, 1, -2, 1], [3, -2, 1, -1], [3, -2, 1, 1], [3, -1, -1, -1], [3, -1, -1, 1]});
    else
      printf("Warning: Should precompute hyperplanes for %a and add to Kronecker:-HyperplaneNormals.\n", ns);
    end if;
  end if;

  # ...otherwise brute-force over all subsets of "codimension-one" cardinality
  normals := {};
  for S in combinat:-choose(Psi_nodup, dim - 1) do
    # check that they span a hyperplane
    M := LinearAlgebra:-Transpose(Matrix(S));
    ker := LinearAlgebra:-NullSpace(M);
    if numelems(ker) > 1 then
      next;
    end if;
    if numelems(ker) = 0 then
      error "Matrix cannot be injective.";
    end if;

    # add normal vector
    normals := normals union {ListMakePrimitiveAndUnique(ker[1])};
  end do;
  return map(Vector, normals);
end proc;


# compute indices of bases in Psi_nodup with respect to the lattice with basis vectors L, and their lcm
Indices := proc(nz, rect, ns, Psi_nodup, L, { use_database := true }, $)
  local dim, ds, S, M, R, d;

  # check dimension of ambient space
  dim := add(ns) - numelems(ns);
  if dim <> numelems(Psi_nodup[1]) then
    error "Dimension mismatch.";
  end if;

  # implement shortcuts for known cases
  if use_database then
    if ns = [2,2] then
      return {1,2}, 2;
    elif ns = [3,2] and nz <= 6 then
      return {1,2,3}, 6;
    elif ns = [2,2,2,2] and nz <= 2 then
      return {1,2,3,4}, 12;
    elif ns = [3,3] and nz <= 3 then
      return {1,2,3,4}, 12;
    #ns = [2,2,2,2] and nz=2 then out:=5;
    #ns = [2,2,2] and nz=3 then out:=4;
    #ns = [3,2,2] and nz=3 then out:=5;
    #ns = [3,3,2] and nz=3 then out:=8; #to check again
    #ns = [3,3,3] and nz=3 then out:=5; # I am not sure I computed
    #ns = [3,2] and nz=6 then out:=3;#I added this row
    elif ns = [2,2,2] and nz = 2 and rect then
      return {1,2,3}, 6;
    elif ns = [4,2] and nz = 4 and rect then
      return {1,2,3,4}, 12;
    else
      printf("Warning: Should precompute indices for %a, %a, %a and add to Kronecker:-Indices.\n", ns, rect, nz);
    end if;
  end if;

  # ...otherwise compute by brute force
  ds := {};
  for S in combinat:-choose(Psi_nodup, dim) do
    # check that the vectors form a basis
    M := Matrix(S);
    if LinearAlgebra:-Rank(M) < dim then
      next;
    end if;

    # compute minimal stretching factor d such that d*Psi is contained in Z[sigma]
    R := LinearAlgebra:-LinearSolve(M, L); #Matrix(RemoveVectorsFromList(S, Psi_nodup, strict=true)));
    d := ilcm(seq(map(denom, R)));
    ds := ds union {d};
  end do;
  return ds, ilcm(seq(ds));
end proc;


##############################################################################
# Weyl groups
##############################################################################

# return Weyl groups of GL(M), of the stabilizer for given (nz, rect), and a set of representatives of the left cosets for the latter
WeylGroups := proc(M, nz, { rect := false }, $)
  local WeylG, PermsL, WeylL;
  WeylG := GroupTheory:-SymmetricGroup(M);
  PermsL := [seq(Perm([[i, i+1]]), i=nz+1..M-1)];
  if rect then
    PermsL := [seq(PermsL), seq(Perm([[i, i+1]]), i=1..nz-1)];
  end if;
  WeylL := GroupTheory:-Subgroup(PermsL, WeylG);
  return WeylG, WeylL, map(Representative, GroupTheory:-LeftCosets(WeylL, WeylG));
end proc;


##############################################################################
# Iterated residues
##############################################################################

# determine vanishing order of function f at x=0
VanishingOrderAtZero := proc(f_, x)
  local f, k;
  f := f_;
  k := 0;
  while subs(x=0, f) = 0 do
    f := simplify(f / x);
    k := k + 1;
  od;
  return k;
end:

# compute iterated residue of function f at x=1
IteratedResidueAtOne := proc(f, dim, x, $)
  local res, i, order;

  # SLOW
  #res := f;
  #for i from 1 to dim do
  #  res := residue(res, x[dim+1-i] = 1);
  #end do;
  #return res;

  res := f;
  for i from 1 to dim do
    res := subs(x[dim+1-i] = 1+t, res);
    res := simplify(res);
    order := VanishingOrderAtZero(denom(res), t);
    if order = 0 then
      return 0;
    end if;
    res := series(res, t=0, order+1);  # XXX: ONE MORE THAN NEEDED.
    res := coeff(res, t, -1);
    res := simplify(res);
  end do;
  return res;
end proc;


##############################################################################
# Kronecker coefficients
##############################################################################

# Compute Kronecker coefficient of given partitions.
Kronecker := proc(parts, sym := 1, { forceperturb := false, coset := NULL, keeptheta := false }, $)
  local M, nz, ns, lambda_numeric, mu_numeric, onerect;
  local R, GLToSL, dim;
  local DeltaPlusG, DeltaPlusL, DeltaU, DeltaPlusK_SL, RootLatticeMatrix;
  local WeylG, WeylL, WeylCosets;
  local Psi_nodup, zero;
  local A, q, qCosets, T, Alphas, AlphasPolarized_nodup, d, ds;
  local lambda_good, mu_good, lambda_tope, mu_tope, lambda, mu, lambda_coset, mu_coset, vars, assumptions;
  local V_tope, V_numeric;
  local w, w_index, V, V_coset, AlphasPolarized_trimmed, qCosets_trimmed, sigmas, sigma_abs_dets, sigma_invs, sigma_invs_times_abs_det_transposed, sigma, Z_from_U, g, res, gamma, gamma_inv, gamma_vec, gamma_res, i, j, PI, thetaq;
  local st;

  # startup
  userinfo(2, Kronecker, "parts =", parts, ", sym =", sym, ", forceperturb =", forceperturb, ", coset =", coset, ", keeptheta =", keeptheta);
  st := time():

  # compute highest weights
  M, nz, ns, lambda_numeric, mu_numeric := ParseParts(parts);
  userinfo(3, Kronecker, "M =", M, ", nz =", nz, ", ns =", ns, ", lambda_numeric =", lambda_numeric, ", mu_numeric =", mu_numeric);

  # use rectangular optimization?
  if lambda_numeric[1] = lambda_numeric[nz] and (is(sym, numeric) or is(sym, symbol)) then
    onerect := AllowRect(nz, ns);
    if onerect then
      userinfo(3, Kronecker, "Rectangular mode is enabled.");
    else
      userinfo(3, Kronecker, "First partition is rectangle but rectangular mode is *not* enabled.");
    end if;
  else
    onerect := false;
  end if;
  
  # compute restriction and GL->SL matrices
  userinfo(5, Kronecker, "Computing restriction and GL->SL matrices...");
  R := RestrictionMatrix(ns);
  GLToSL := GLToSLMatrix(ns);
  dim := add(ns) - numelems(ns);

  # compute positive roots & those that are orthogonal to highest weight
  userinfo(5, Kronecker, "Computing restricted root systems...");
  DeltaPlusG := PositiveRoots(M);
  DeltaPlusL := [seq(seq(BasisVector(i, M) - BasisVector(j, M), j=i+1..M), i=nz+1..M)];
  if onerect then
    DeltaPlusL := [seq(DeltaPlusL), seq(seq(BasisVector(i, M) - BasisVector(j, M), j=i+1..nz), i=1..nz)];
  end if;
  # TODO: allow for additional Sigma (with the user promising that C_{K,K}^Sigma remains solid)
  DeltaU := RemoveVectorsFromList(DeltaPlusL, DeltaPlusG, strict=true);
  DeltaPlusK_SL := [seq(GLToSL . alpha, alpha=PositiveRootsProduct(ns))];
  RootLatticeMatrix := LinearAlgebra:-IdentityMatrix(dim);

  # compute Weyl groups and cosets
  userinfo(5, Kronecker, "Computing Weyl groups and coset...");
  WeylG, WeylL, WeylCosets := WeylGroups(M, nz, rect=onerect);
  if numelems(select(w -> ListPermute(w, lambda_numeric) <> lambda_numeric, GroupTheory:-Generators(WeylL))) > 0 then
    error "Internal error in computation of WeylL.";
  end if;

  # compute ordered set of nonzero restricted roots
  Psi_nodup := [seq(GLToSL . R . alpha, alpha=DeltaPlusG)];
  Psi_nodup := UniqueVectorList(Psi_nodup);
  if numelems(select(IsZeroVector, Psi_nodup)) > 0 then
    error "TODO: Zero vector after restriction. Should be fine, but I want to debug this once.";
  end if;

  # verify that the restricted roots span the root lattice (which, in our conventions, is Z^dim)
  if not IsVectorSubset([LinearAlgebra:-Column(RootLatticeMatrix, 1..-1)], Psi_nodup) then
    error "Internal error: Restricted roots should contain DeltaPlusK";
  end if;
  if numelems(select(x -> not is(x, integer), [seq(op(convert(psi, list)), psi=Psi_nodup)])) > 0 then
    error "Internal error: Restricted roots should span root lattice";
  end if;

  # compute hyperplane normals
  userinfo(5, Kronecker, "Computing hyperplane normals...");
  A := HyperplaneNormals(ns, Psi_nodup);

  # compute index (we do not use Psi_nodup but the following smaller system; note that it is a subset of Psi_nodup and hence automatically polarized with respect to the latter)
  Alphas := [seq(GLToSL . R . Vector(alpha), alpha=DeltaU)];
  AlphasPolarized_nodup := UniqueVectorList(RemoveVectorsFromList(DeltaPlusK_SL, Alphas, strict=false));
  ds, q := Indices(nz, onerect, ns, AlphasPolarized_nodup, RootLatticeMatrix);
  qCosets := {};
  for d in ds do
    T := combinat:-cartprod([seq([q/d*seq(0..d-1)], i = 1 .. dim)]);
    while not T:-finished do
      gamma := T:-nextvalue();

      # optimization: for each gamma we either compute gamma or its inverse
      gamma_inv := map(x -> (q - x) mod q, gamma);
      if not gamma_inv in qCosets then
        qCosets := qCosets union {gamma};
      end if;
    end do;
  end do;

  # determine tope by shifting vector (if appropriate) & setup symbolic partition (if appropriate)
  if is(sym, numeric) or is(sym, symbol) then
    # deform towards a tope
    lambda_good, mu_good := ParseParts(GoodParts(nz, onerect, ns))[4..5];
    lambda_tope, mu_tope := lambda_numeric + lambda_good, mu_numeric + mu_good;

    # stretch (this makes the partition symbolic, if necessary)
    lambda, mu := Vector(lambda_numeric) * sym, Vector(mu_numeric) * sym;
    if is(sym, symbol) then
      assumptions := [sym :: nonnegint];
    end if;
  elif is(sym, list) then
    if numelems(sym) <> numelems(parts) then
      error "Expect as many symbolic variables as there are partitions: %1", sym;
    end if;

    # in this case we only deform if explicitly asked
    if forceperturb then
      lambda_good, mu_good := ParseParts(GoodParts(nz, onerect, ns))[4..5];
      lambda_tope, mu_tope := lambda_numeric + lambda_good, mu_numeric + mu_good;
    else
      lambda_tope, mu_tope := lambda_numeric, mu_numeric;
    end if;

    # introduce symbolic variables for the first partition
    lambda := Vector([seq(sym[1][i], i=1..nz), seq(0, i=nz+1..M)]);
    mu := Vector([seq(op([seq(sym[1+k][i], i=1..ns[k]-1), add(lambda) - add(seq(sym[1+k][i], i=1..ns[k]-1))]), k=1..numelems(ns))]);

    # assume that the variables are all nonnegative integers
    vars := [seq(sym[1][i], i=1..nz), seq(seq(sym[1+k][i], i=1..ns[k]-1), k=1..numelems(ns))];
    assumptions := [seq(v :: nonnegint, v=vars)];
  else
    error "Unexpected parameter %1 of type %2.", sym, whattype(sym);
  end if;

  # force coset?
  if coset <> NULL then
    if is(sym, symbol) then
      lambda_coset, mu_coset := Vector(lambda_numeric) * coset, Vector(mu_numeric) * coset;
    else
      error "Unexpected coset parameter; can only use when computing stretching (quasi)polynomial.";
    end if;
  else
    lambda_coset, mu_coset := lambda, mu;
  end if;

  # verify that, for all permutations w, (w, lambda_tope, mu_tope) determine a tope which contains (w, lambda_numeric, mu_numeric) in its closure
  for w in WeylCosets do
    # compute the two projections
    V_tope := GLToSL . (R . Vector(ListPermute(w, lambda_tope)) - Vector(mu_tope));
    V_numeric := GLToSL . (R . Vector(ListPermute(w, lambda_numeric)) - Vector(mu_numeric));

    # V_tope should have strict inequalities only
    if numelems(select(X -> V_tope . X = 0, A)) > 0 then
      if is(sym, list) then
        error "Unsure whether partitions determine a unique tope (use forceperturb=true to perturb anyways)";
      else
        error "Internal error: V_tope does not determine a unique tope";
      end if;
    end if;

    # V_numeric should be in the closure of the tope determined by V_tope
    if numelems(select(X -> (V_tope . X) * (V_numeric . X) < 0, A)) > 0 then
      error "Internal error: V_numeric is not in closure of tope determined by V_tope";
    end;
  end do;

  # TODO: skip if V is not in the root lattice: we can only check this in the non-symbolic case, and it should suffice to check this for w = id

  # compute the residue formula
  userinfo(4, Kronecker, "Entering main loop,", numelems(WeylCosets), "Weyl cosets...");
  g := 0;
  w_index := 0;
  PI := proc(V, $) return mul(x[i]^V[i], i=1..dim); end proc;
  for w in WeylCosets do
    w_index := w_index + 1;
    userinfo(4, Kronecker, "Iteration", w_index, "of", numelems(WeylCosets), ", w =", convert(w, list));

    # determine projection and tope
    V := GLToSL . (R . Vector(ListPermute(w, lambda)) - Vector(mu));
    V_coset := GLToSL . (R . Vector(ListPermute(w, lambda_coset)) - Vector(mu_coset));
    V_tope := GLToSL . (R . Vector(ListPermute(w, lambda_tope)) - Vector(mu_tope));
    Alphas := [seq(GLToSL . R . Vector(ListPermute(w, alpha)), alpha = DeltaU)];

    # prune the cosets using the basis condition for AlphasPolarized_nodup (TODO: possible optimization: prune by checking that V_tope is contained the cone spanned by the AlphasPolarized; break once a basis has been found)
    AlphasPolarized_nodup := UniqueVectorList(RemoveVectorsFromList(DeltaPlusK_SL, PolarizeVectors(Alphas, Psi_nodup), strict=false));
    if numelems(OS(AlphasPolarized_nodup, A, opt_v_=V_tope, opt_dim_=dim)) = 0 then
      next;
    end if;

    # iterate over cosets...
    for gamma in qCosets do
      # setup residue formula
      userinfo(5, Kronecker, "  gamma =", gamma);
      gamma_inv := map(x -> (q - x) mod q, gamma);
      gamma_res := 0;

      # collect polarized alphas with alpha . gamma = 0 and check that they can possibly form a basis
      gamma_vec := Vector(gamma);
      AlphasPolarized_trimmed := select(alpha -> modp(alpha . gamma_vec, q) = 0, AlphasPolarized_nodup);
      if numelems(AlphasPolarized_trimmed) = 0 then
        next;
      end if;
      # SLOWER?
      #if LinearAlgebra:-Rank(Matrix(AlphasPolarized_trimmed)) < dim then
      #  next;
      #end if;

      # compute OS bases
      sigmas := map(Matrix, OS(AlphasPolarized_trimmed, A, opt_v_=V_tope, opt_dim_=dim));  # Transpose because sigma contains column vectors!
      sigma_abs_dets := map(abs, map(LinearAlgebra:-Determinant, sigmas));
      sigma_invs := map(LinearAlgebra:-MatrixInverse, sigmas);
      sigma_invs_times_abs_det_transposed := map(LinearAlgebra:-Transpose, zip(`*`, sigma_abs_dets, sigma_invs));
      if numelems(sigmas) = 0 then
        next;
      end if;

      # iterate over OS bases...
      userinfo(5, Kronecker, "    Computing residues for ", numelems(sigmas), "OS bases");
      for j from 1 to numelems(sigmas) do
        # setup function of which to calculate the residue and base change Z = det(Sigma) Sigma^{-1} log(X) (don't forget Jacobian - this is a volume form!)
        res := mul(seq(1 - PI(-sigma_invs_times_abs_det_transposed[j] . beta) * qRoot^(-modp(beta . gamma_vec, q)), beta = DeltaPlusK_SL))
           * PI(sigma_invs_times_abs_det_transposed[j] . V) * qRoot^(modp(V_coset . gamma_vec, q))
           / mul(seq(1 - PI(-sigma_invs_times_abs_det_transposed[j] . alpha) * qRoot^(-modp(alpha . gamma_vec, q)), alpha = Alphas));
        res := res * sigma_abs_dets[j]^(dim-1) / mul(x[i], i=1..dim);

        # SLOWER:
        #Z_from_U := sigma_invs[j] . U;
        #Z_from_X := sigma_abs_dets[j] * sigma_invs[j] . map(log, X);
        #res := mul(seq(1 - exp(-beta . Z) * evala(qRoot^(-beta . gamma_vec)), beta = DeltaPlusK_SL))
        #   * exp(V . Z) * evala(qRoot^(V . gamma_vec))
        #   / mul(seq(1 - exp(-alpha . Z) * evala(qRoot^(-alpha . gamma_vec)), alpha = Alphas));
        #res := subs(seq(Z[i] = Z_from_X[i], i=1..dim), res) * sigma_abs_dets[j]^(dim-1) / mul(X[i], i=1..dim);
        #res := simplify(res);

        # compute iterated residue
        res := IteratedResidueAtOne(res, dim, x);

        if res <> 0 then
          gamma_res := gamma_res + res;
          userinfo(5, Kronecker, "    Nonzero residue for w =", convert(w, list), ", gamma =", gamma, "=> res =", res);
        end if;
      end do;

      # add residues, and also the term that we would have obtained from the inverse
      if gamma <> gamma_inv then
        ASSERT(not gamma_inv in qCosets);
        gamma_res := gamma_res + subs(qRoot=1/qRoot, gamma_res);
      end if;

      # simplify roots and exponents as much as possible
      gamma_res := simplify(gamma_res);
      thetaq := q;
      if q = 2 then
        gamma_res := subs(qRoot = -1, gamma_res);
      elif q = 4 then
        gamma_res := subs(qRoot = I, gamma_res);
      elif q = 6 then
        gamma_res := applyrule(qRoot^a::algebraic = (-1)^a * theta[3]^'modp(a, 3)', gamma_res);
        thetaq := 3;
      elif q = 12 then
        gamma_res := applyrule(qRoot^a::algebraic = I^a * theta[3]^'modp(a, 3)', gamma_res);
        thetaq := 3;
      else
        gamma_res := applyrule(qRoot^a::algebraic = theta[q]^'modp(a, q)', gamma_res);
      end if;
      if thetaq = 3 then
        gamma_res := algsubs(theta[3]^2 = -theta[3]-1, gamma_res);
      end if;

      # plug in root of cyclotomic polynomial?
      if not keeptheta then
        gamma_res := subs(theta[thetaq] = RootOf(numtheory:-cyclotomic(thetaq, t)), gamma_res);
      end if;
      gamma_res := simplify(evala(simplify(gamma_res))) assuming seq(assumptions);

      # add to coeff
      g := g + gamma_res;
    end do;
  end do;
  userinfo(5, Kronecker, "Simplifying result...");
  userinfo(5, Kronecker, "Finished in", time() - st, "seconds.");
  return simplify(g);
end proc;

ModuleApply := Kronecker;

end module;
