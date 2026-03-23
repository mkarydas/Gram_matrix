(* gram_matrix.wl
   Gram matrix effective dimension — Wolfram Language implementation.
   Run with: math -script gram_matrix.wl
   Or open in Mathematica and evaluate all cells.
*)

(* ── Kernels ──────────────────────────────────────────────────────────── *)

linearKernel[x_, y_] := x . y

polynomialKernel[c_, d_][x_, y_] := (x . y + c)^d

gaussianKernel[sigma_][x_, y_] := Exp[-(x - y) . (x - y) / (2 sigma^2)]

laplacianKernel[sigma_][x_, y_] := Exp[-Total[Abs[x - y]] / sigma]


(* ── Discretize domain ────────────────────────────────────────────────── *)

(* 1D: domain = {a, b}
   nD: domain = {{a1,b1},{a2,b2},...}
   method: "Grid" or "Random" *)

discretize[domain:{_?NumericQ, _?NumericQ}, nPoints_Integer, method_:"Grid"] :=
  discretize[{domain}, nPoints, method]

discretize[domain:{{_?NumericQ, _?NumericQ}..}, nPoints_Integer, method_:"Grid"] :=
  Module[{axes, tuples},
    Switch[method,
      "Grid",
        axes = Table[Subdivide[d[[1]], d[[2]], nPoints - 1], {d, domain}];
        tuples = Tuples[axes];
        If[Length[domain] == 1, Flatten[tuples], tuples],
      "Random",
        Transpose[Table[
          RandomReal[{d[[1]], d[[2]]}, nPoints], {d, domain}]],
      _, Message[discretize::badmethod, method]; $Failed
    ]
  ]

discretize::badmethod = "Unknown method `1`. Use \"Grid\" or \"Random\".";


(* ── Build Gram matrix ────────────────────────────────────────────────── *)

buildGramMatrix[kernel_, points_List] :=
  Outer[kernel, points, points, 1]


(* ── Effective dimension ──────────────────────────────────────────────── *)

(* threshold: cutoff value
   relative -> True:  cutoff = threshold * lambda_max
   relative -> False: cutoff = threshold (absolute) *)

effectiveDimension[G_?MatrixQ, threshold_:1*^-6, relative_:True] :=
  Module[{eigs, cutoff},
    eigs = Sort[Eigenvalues[N[G]], Greater];   (* descending *)
    cutoff = If[relative, threshold * First[eigs], threshold];
    {Count[eigs, e_ /; e > cutoff], eigs}
  ]


(* ── Example ──────────────────────────────────────────────────────────── *)

n = 50;
domain = {0.0, 1.0};
threshold = 1*^-6;

points = discretize[domain, n, "Grid"];

kernels = {
  {"Linear",           linearKernel},
  {"Polynomial d=2",   polynomialKernel[1, 2]},
  {"Polynomial d=5",   polynomialKernel[1, 5]},
  {"Gaussian \[Sigma]=1.0",  gaussianKernel[1.0]},
  {"Gaussian \[Sigma]=0.1",  gaussianKernel[0.1]},
  {"Gaussian \[Sigma]=0.01", gaussianKernel[0.01]}
};

results = Table[
  Module[{name, kernel, G, dim, eigs, nonzero, lambdaMin},
    {name, kernel} = row;
    G = buildGramMatrix[kernel, points];
    {dim, eigs} = effectiveDimension[G, threshold, True];
    nonzero = Select[eigs, # > threshold * First[eigs] &];
    lambdaMin = If[nonzero =!= {}, Last[nonzero], 0.0];
    {name, dim, First[eigs], lambdaMin}
  ],
  {row, kernels}
];

(* Print table *)
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
      StringPadLeft[ScientificForm[row[[4]], 3], 20]]]
  ],
  results
];
