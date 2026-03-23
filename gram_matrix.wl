(* gram_matrix.wl
   Gram matrix effective dimension — Wolfram Language implementation.
   Run with: math -script gram_matrix.wl
   Or open in Mathematica and evaluate all cells.

   Speed notes:
   - Distances computed via ||x-y||^2 = ||x||^2 + ||y||^2 - 2 x.y  (pure BLAS)
   - Arrays forced to packed form with Developer`ToPackedArray
   - Eigenvalues uses the real-symmetric solver (La, not general complex path)
   - No DistanceMatrix dependency (unavailable in some script-mode kernels)
*)

Needs["Developer`"]

(* ── Helpers ──────────────────────────────────────────────────────────── *)

(* Convert points to a packed numerical 2D matrix (N x d) *)
toMatrix[pts_] :=
  Developer`ToPackedArray[
    If[VectorQ[pts, NumericQ], List /@ N[pts], N[pts]]
  ]

(* Squared Euclidean distance matrix via BLAS dot-product trick:
   D2[i,j] = ||x_i||^2 + ||x_j||^2 - 2 x_i.x_j                *)
squaredDistances[p_] :=
  Module[{norms = Total[p^2, {2}]},
    Developer`ToPackedArray[
      Outer[Plus, norms, norms] - 2. (p . Transpose[p])
    ]
  ]

(* L1 distance matrix: sum of |x_i_k - x_j_k| over dimensions *)
l1Distances[p_] :=
  Developer`ToPackedArray[
    Total[Abs[Outer[Subtract, p, p, 1]], {3}]
  ]


(* ── Vectorized Gram matrix builders ─────────────────────────────────── *)

gramLinear[pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[p . Transpose[p]]
  ]

gramPolynomial[c_, d_][pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[(p . Transpose[p] + c)^d]
  ]

gramGaussian[sigma_][pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[Exp[-squaredDistances[p] / (2. sigma^2)]]
  ]

gramLaplacian[sigma_][pts_] :=
  Module[{p = toMatrix[pts]},
    Developer`ToPackedArray[Exp[-l1Distances[p] / sigma]]
  ]

(* Custom scalar kernel fallback *)
gramCustom[kernel_][pts_] :=
  Developer`ToPackedArray[N[Outer[kernel, pts, pts, 1]]]


(* ── Discretize domain ────────────────────────────────────────────────── *)

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


(* ── Effective dimension ──────────────────────────────────────────────── *)

effectiveDimension[G_?MatrixQ, threshold_ : 1*^-6, relative_ : True] :=
  Module[{eigs, cutoff},
    (* Symmetric real matrix: force the real-symmetric LAPACK path *)
    eigs = Eigenvalues[N[G], Method -> "LAPACK"];
    (* Eigenvalues returns descending order — no Sort needed *)
    cutoff = If[relative, threshold * First[eigs], threshold];
    {Count[eigs, e_ /; e > cutoff], eigs}
  ]


(* ── Example ──────────────────────────────────────────────────────────── *)

n         = 50;
domain    = {0.0, 1.0};
threshold = 1*^-6;

points = discretize[domain, n, "Grid"];

kernelBuilders = {
  {"Linear",              gramLinear},
  {"Polynomial d=2",      gramPolynomial[1, 2]},
  {"Polynomial d=5",      gramPolynomial[1, 5]},
  {"Gaussian s=1.0",      gramGaussian[1.0]},
  {"Gaussian s=0.1",      gramGaussian[0.1]},
  {"Gaussian s=0.01",     gramGaussian[0.01]}
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

(* Print results table *)
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
