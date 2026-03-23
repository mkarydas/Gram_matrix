(* gram_matrix.wl
   Gram matrix effective dimension — Wolfram Language implementation.
   Run with: math -script gram_matrix.wl
   Or open in Mathematica and evaluate all cells.

   Optimized: kernel matrices are built via vectorized numerical operations
   (DistanceMatrix, tensor dot products) instead of element-wise Outer calls.
*)

(* ── Helpers ──────────────────────────────────────────────────────────── *)

(* Ensure points are a numerical 2D matrix: (N, d) *)
toMatrix[pts_] := If[VectorQ[pts, NumericQ], List /@ N[pts], N[pts]]


(* ── Vectorized Gram matrix builders ─────────────────────────────────── *)
(* Each function takes a list of points and returns the full N x N matrix
   using bulk numerical operations — no element-wise kernel calls.         *)

gramLinear[pts_] :=
  Module[{p = toMatrix[pts]}, p . Transpose[p]]

gramPolynomial[c_, d_][pts_] :=
  Module[{p = toMatrix[pts]}, (p . Transpose[p] + c)^d]

gramGaussian[sigma_][pts_] :=
  Exp[-DistanceMatrix[toMatrix[pts]]^2 / (2. sigma^2)]

gramLaplacian[sigma_][pts_] :=
  Exp[-DistanceMatrix[toMatrix[pts],
        DistanceFunction -> ManhattanDistance] / sigma]

(* For a custom scalar kernel f[x,y]: falls back to Outer, still numerical *)
gramCustom[kernel_][pts_] :=
  N[Outer[kernel, pts, pts, 1]]


(* ── Discretize domain ────────────────────────────────────────────────── *)

(* 1D: domain = {a, b}
   nD: domain = {{a1,b1},{a2,b2},...}
   method: "Grid" or "Random" *)

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

(* threshold: cutoff value
   relative -> True:  cutoff = threshold * lambda_max
   relative -> False: cutoff = threshold (absolute)
   Returns {effectiveDim, eigenvalues (descending)}                        *)

effectiveDimension[G_?MatrixQ, threshold_ : 1*^-6, relative_ : True] :=
  Module[{eigs, cutoff},
    (* eigvalsh equivalent: symmetric matrix, real eigenvalues *)
    eigs = Sort[Eigenvalues[N[G]], Greater];
    cutoff = If[relative, threshold * First[eigs], threshold];
    {Count[eigs, e_ /; e > cutoff], eigs}
  ]


(* ── Example ──────────────────────────────────────────────────────────── *)

n         = 50;
domain    = {0.0, 1.0};
threshold = 1*^-6;

points = discretize[domain, n, "Grid"];

(* Each entry: {display name, vectorized gram-matrix builder} *)
kernelBuilders = {
  {"Linear",              gramLinear},
  {"Polynomial d=2",      gramPolynomial[1, 2]},
  {"Polynomial d=5",      gramPolynomial[1, 5]},
  {"Gaussian \[Sigma]=1.0",   gramGaussian[1.0]},
  {"Gaussian \[Sigma]=0.1",   gramGaussian[0.1]},
  {"Gaussian \[Sigma]=0.01",  gramGaussian[0.01]}
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
Print[StringForm["`` ``  ``  ``",
  StringPadRight["Kernel", 22],
  StringPadLeft["Effective dim", 14],
  StringPadLeft["\[Lambda]_max", 12],
  StringPadLeft["\[Lambda]_min (non-zero)", 20]]];
Print[StringRepeat["-", 72]];
Scan[
  Function[row,
    Print[StringForm["`` ``  ``  ``",
      StringPadRight[row[[1]], 22],
      StringPadLeft[ToString[row[[2]]], 14],
      StringPadLeft[ToString[NumberForm[row[[3]], {6, 4}]], 12],
      StringPadLeft[ToString[ScientificForm[row[[4]], 3]], 20]]]
  ],
  results
];
