# Gram Matrix — Mathematica / Wolfram Language

This document implements the Gram matrix effective dimension computation in Wolfram Language.
Run as a script: `math -script gram_matrix.wl`, or paste sections into a Mathematica notebook.

---

## Setup

```mathematica
Needs["Developer`"]
```

---

## Helpers

Convert any list of points to a packed numerical `(N × d)` matrix.

```mathematica
toMatrix[pts_] :=
  Developer`ToPackedArray[
    If[VectorQ[pts, NumericQ], List /@ N[pts], N[pts]]
  ]
```

Squared Euclidean distance matrix via the BLAS dot-product identity
`‖x−y‖² = ‖x‖² + ‖y‖² − 2 x·y` — no `DistanceMatrix` dependency.

```mathematica
squaredDistances[p_] :=
  Module[{norms = Total[p^2, {2}]},
    Developer`ToPackedArray[
      Outer[Plus, norms, norms] - 2. (p . Transpose[p])
    ]
  ]
```

L1 distance matrix.

```mathematica
l1Distances[p_] :=
  Developer`ToPackedArray[
    Total[Abs[Outer[Subtract, p, p, 1]], {3}]
  ]
```

---

## Kernel Functions

Each builder takes a list of points and returns the full `N × N` Gram matrix using bulk numerical operations — no element-wise kernel calls.

### Linear  `K(x,y) = x · y`

```mathematica
gramLinear[pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[p . Transpose[p]]
  ]
```

### Polynomial  `K(x,y) = (x · y + c)^d`

```mathematica
gramPolynomial[c_, d_][pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[(p . Transpose[p] + c)^d]
  ]
```

### Gaussian  `K(x,y) = exp(−‖x−y‖² / 2σ²)`

```mathematica
gramGaussian[sigma_][pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[Exp[-squaredDistances[p] / (2. sigma^2)]]
  ]
```

### Laplacian  `K(x,y) = exp(−‖x−y‖₁ / σ)`

```mathematica
gramLaplacian[sigma_][pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[Exp[-l1Distances[p] / sigma]]
  ]
```

### Custom kernel (fallback)

```mathematica
gramCustom[kernel_][pts_] :=
  Developer`ToPackedArray[N[Outer[kernel, pts, pts, 1]]]
```

---

## Discretize Domain

```mathematica
discretize[domain : {_?NumericQ, _?NumericQ}, nPoints_Integer,
           method_ : "Grid"] :=
  discretize[{domain}, nPoints, method]

discretize[domain : {{_?NumericQ, _?NumericQ} ..}, nPoints_Integer,
           method_ : "Grid"] :=
  Module[{axes, tuples},
    Switch[method,
      "Grid",
        axes = Table[N @ Subdivide[d[[1]], d[[2]], nPoints - 1], {d, domain}];
        tuples = Tuples[axes];
        If[Length[domain] == 1, Flatten[tuples], tuples],
      "Random",
        Transpose[Table[RandomReal[{d[[1]], d[[2]]}, nPoints], {d, domain}]],
      _, Message[discretize::badmethod, method]; $Failed
    ]
  ]

discretize::badmethod = "Unknown method `1`. Use \"Grid\" or \"Random\".";
```

---

## Effective Dimension

Counts eigenvalues above a threshold. Uses `Method -> "LAPACK"` to force the real-symmetric solver.

- `relative = True` (default): cutoff = `threshold × λ_max`
- `relative = False`: absolute cutoff

Returns `{effectiveDim, eigenvalues}` in descending order.

```mathematica
effectiveDimension[G_?MatrixQ, threshold_ : 1*^-6, relative_ : True] :=
  Module[{eigs, cutoff},
    eigs   = Eigenvalues[N[G], Method -> "LAPACK"];
    cutoff = If[relative, threshold * First[eigs], threshold];
    {Count[eigs, e_ /; e > cutoff], eigs}
  ]
```

---

## Example

Compare six kernels on a 50-point uniform grid over `[0, 1]`.

```mathematica
n         = 50;
domain    = {0.0, 1.0};
threshold = 1*^-6;

points = discretize[domain, n, "Grid"];

kernelBuilders = {
  {"Linear",          gramLinear},
  {"Polynomial d=2",  gramPolynomial[1, 2]},
  {"Polynomial d=5",  gramPolynomial[1, 5]},
  {"Gaussian s=1.0",  gramGaussian[1.0]},
  {"Gaussian s=0.1",  gramGaussian[0.1]},
  {"Gaussian s=0.01", gramGaussian[0.01]}
};

results = Table[
  Module[{name, builder, G, dim, eigs, nonzero, lambdaMin},
    {name, builder} = row;
    G       = builder[points];
    {dim, eigs} = effectiveDimension[G, threshold, True];
    nonzero = Select[eigs, # > threshold * First[eigs] &];
    lambdaMin = If[nonzero =!= {}, Last[nonzero], 0.0];
    {name, dim, First[eigs], lambdaMin}
  ],
  {row, kernelBuilders}
];

Print[""];
Print[StringJoin[
  StringPadRight["Kernel", 22], "  ",
  StringPadLeft["Eff. dim", 10], "  ",
  StringPadLeft["lambda_max", 14], "  ",
  StringPadLeft["lambda_min", 14]]];
Print[StringRepeat["-", 66]];
Scan[
  Function[row,
    Print[StringJoin[
      StringPadRight[row[[1]], 22], "  ",
      StringPadLeft[ToString[row[[2]]], 10], "  ",
      StringPadLeft[ToString[NumberForm[row[[3]], {8, 4}]], 14], "  ",
      StringPadLeft[ToString[ScientificForm[row[[4]], 3]], 14]]]
  ],
  results
];
```
