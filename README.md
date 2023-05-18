# A Maple package for Kronecker coefficients [![arXiv](https://img.shields.io/badge/arxiv-1601.04325-blue)](http://arxiv.org/abs/1601.04325) [![download](https://img.shields.io/badge/download-zip-red)](https://github.com/qi-rub/kronecker/zipball/master)

`Kronecker` is a Maple package for computing Kronecker coefficients g(λ,μ,ν). It implements a highly optimized version of the algorithm proposed by [Baldoni-Vergne-Walter (2017)](http://arxiv.org/abs/1601.04325). `Kronecker` works just as well symbolically, computing *dilated* or *stretched* Kronecker coefficients g(kλ,kμ,kν) and even the entire *quasipolynomial* on a polyhedral chamber of (λ,μ,ν)'s.


# Getting started

Download the zip file [**here**](https://github.com/qi-rub/kronecker/zipball/master), unpack it, and run the `GettingStarted.mw` worksheet.
You should see the following code and output at the top of the file:
```maple
> restart:
  read "Kronecker.mpl":

Kronecker v1.0 by Velleda Baldoni, Michèle Vergne, Michael Walter (see https://github.com/qi-rub/kronecker for more information).
```
If you instead see an error message please let us now by [opening a new bug report](https://github.com/qi-rub/kronecker/issues/new).
The worksheet also contains the following illustrative examples demonstrating the functionality of the Kronecker package:

## Calculate individual coefficients

Pass a list of partitions [λ,μ,ν] to `Kronecker` to compute the corresponding Kronecker coefficient g(λ,μ,ν):

```maple
> Kronecker([[131, 39, 19, 11], [110, 90], [120, 80]]);

70

> Kronecker([[15, 10, 9, 4, 3, 2], [21, 14, 8], [27, 16]]);

148

> Kronecker([[400000, 200000, 100000], [500000, 100000, 100000], [300000, 200000, 200000]]);

1596403568763802677125206373630515086
```

This last example is a branching problem for SL(9) and hence might take a little while to compute (around 6 minutes).

## Calculate dilated coefficients

By providing an additional symbolic parameter to `Kronecker`, we can instead compute the dilated Kronecker coefficient, i.e., the quasipolynomial function p(k) = g(kλ,kμ,kν):

```maple
> Kronecker([[10, 6, 2], [10, 8], [11, 7]], k);

3/8 + (5/8)*(-1)^k + (3/2)*k + (7/4)*k^2

> Kronecker([[1, 1, 1, 1, 1, 1], [2, 2, 2], [3, 3]], k);

1
```

We can also probe individual cosets by using the `coset` option:

```maple
> Kronecker([[1, 1, 1, 1], [1, 1, 1, 1], [2, 2]], k, coset = 0);

(1/6)*k + 1

> Kronecker([[1, 1, 1, 1], [1, 1, 1, 1], [2, 2]], k, coset = 1);

-1/6 + (1/6)*k
```

```maple
> Kronecker([[1, 1], [1, 1], [1, 1], [1, 1]], k, coset = 0);

(2/3)*k + (1/6)*k^2 + (1/72)*k^3 + 1

> Kronecker([[1, 1], [1, 1], [1, 1], [1, 1]], k, coset = 1);

5/18 + (13/24)*k + (1/6)*k^2 + (1/72)*k^3
```

The last example shows that the `Kronecker` package also supports iterated Kronecked coefficients indexed by four or more partitions.

## Calculate quasipolynomials on a chamber

Lastly, by providing three symbolic variables we can compute the entire quasipolynomial on the corresponding chamber:

```maple
> Kronecker([[132, 38, 19, 11], [110, 90], [120, 80]], [lambda, mu, nu]);

3/4 + (1/8)*(-1)^(lambda[4]+mu[1]+lambda[2]+nu[1]) - (1/2)*lambda[4]*mu[1] - (1/2)*lambda[4]*lambda[2] + (1/4)*lambda[4]^2 + (1/2)*lambda[4]*nu[1] + (1/2)*lambda[2] + (1/2)*lambda[3] - lambda[4] + (1/2)*mu[1] - (1/2)*nu[1] - (1/4)*lambda[3]^2 - (1/2)*lambda[3]*nu[1] + (1/2)*lambda[3]*mu[1] + (1/2)*lambda[3]*lambda[2] + (1/8)*(-1)^(lambda[3]+mu[1]+lambda[2]+nu[1])
```

Sometimes, the chamber and hence the quasipolynomial are not uniquely determined, because the given partition triple sits on a wall between two chambers.
In other cases, such as in the following, the chamber is uniquely determined, but `Kronecker` is not aware of it since the partition triple sits on the boundary of the Kirwan cone.
In either case the `forceperturb` option can be used to go ahead and select a chamber:

```maple
> factor(Kronecker([[288, 192, 174, 120, 30, 6], [343, 270, 197], [654, 156]], [lambda, mu, nu], forceperturb = true));

1/5040 * (-nu[1]+7+lambda[1]+lambda[2]+lambda[3]) * (-nu[1]+6+lambda[3]+lambda[2]+lambda[1]) * (-nu[1]+5+lambda[3]+lambda[2]+lambda[1]) * (-nu[1]+4+lambda[3]+lambda[2]+lambda[1]) * (-nu[1]+3+lambda[3]+lambda[2]+lambda[1]) * (-nu[1]+2+lambda[3]+lambda[2]+lambda[1]) * (-nu[1]+1+lambda[3]+lambda[2]+lambda[1]) * (lambda[1]+lambda[2]+lambda[4]+lambda[5]-mu[1]-mu[2]+1)
```


# Attribution

If you find this software useful in your research please consider citing our paper:

```
@article{kronecker,
  title={Computation of dilated Kronecker coefficients},
  author={Baldoni, Velleda and Vergne, Mich{\`e}le},
  journal={Journal of Symbolic Computation},
  year={2017},
  doi={10.1016/j.jsc.2017.03.005},
  eprint={1601.04325},
  note={In press. With an appendix by Michael Walter. Software available at \url{https://qi-rub.github.io/kronecker/}.},
}
```


# See also

- [barvikron](https://github.com/qi-rub/barvikron): a Python package for efficiently computing Kronecker coefficients (using Barvinok's algorithm to evaluate characters rather than the iterated residues used here)
- [LiE](http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/): a computer algebra system for reductive Lie group computations
- [SageMath](http://sagemath.org/): a computer algebra system which includes support for [symmetric functions](http://doc.sagemath.org/html/en/reference/combinat/sage/combinat/sf/sfa.html#sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.itensor)
