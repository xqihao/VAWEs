(* ::Package:: *)

Linspace[x0_, x1_, n_] := Range[x0, x1, (x1 - x0)/(n - 1)];


Options[plotGrid] = {ImagePadding -> 20};
plotGrid[l_List, w_, h_, opts : OptionsPattern[]] := 
  Module[{nx, ny, sidePadding = OptionValue[plotGrid, ImagePadding], 
    topPadding = 0, widths, heights, dimensions, positions, 
    frameOptions = 
     FilterRules[{opts}, 
      FilterRules[Options[Graphics], 
       Except[{ImagePadding, Frame, FrameTicks}]]]}, {ny, nx} = 
    Dimensions[l];
   widths = (w - 2 sidePadding)/nx Table[1, {nx}];
   widths[[1]] = widths[[1]] + sidePadding;
   widths[[-1]] = widths[[-1]] + sidePadding;
   heights = (h - 2 sidePadding)/ny Table[1, {ny}];
   heights[[1]] = heights[[1]] + sidePadding;
   heights[[-1]] = heights[[-1]] + sidePadding;
   positions = 
    Transpose@
     Partition[
      Tuples[Prepend[Accumulate[Most[#]], 0] & /@ {widths, heights}], 
      ny];
   Graphics[
    Table[Inset[
      Show[l[[ny - j + 1, i]], 
       ImagePadding -> {{If[i == 1, sidePadding, 0], 
          If[i == nx, sidePadding, 0]}, {If[j == 1, sidePadding, 0], 
          If[j == ny, sidePadding, topPadding]}}, 
       AspectRatio -> Full], 
      positions[[j, i]], {Left, Bottom}, {widths[[i]], 
       heights[[j]]}], {i, 1, nx}, {j, 1, ny}], 
    PlotRange -> {{0, w}, {0, h}}, ImageSize -> {w, h}, 
    Evaluate@Apply[Sequence, frameOptions]]];


VAthomsen[c11_, c12_, c13_, c22_, c23_, c33_] := 
  Module[{vp0, vnmo1, 
    vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, \[Epsilon]1, \[Delta]1, \
\[Epsilon]2, \[Delta]2, \[Delta]3, 
    Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q, c11R, c12R, c13R, c22R, c23R, c33R, iQ11, iQ12, iQ13, 
    iQ22, iQ23, iQ33},
   c11R = Re[c11];
   c12R = Re[c12];
   c13R = Re[c13];
   c22R = Re[c22];
   c23R = Re[c23];
   c33R = Re[c33];
   
   iQ11 = -(Im[c11]/Re[c11]);
   iQ12 = -(Im[c12]/Re[c12]);
   iQ13 = -(Im[c13]/Re[c13]);
   iQ22 = -(Im[c22]/Re[c22]);
   iQ23 = -(Im[c23]/Re[c23]);
   iQ33 = -(Im[c33]/Re[c33]);
   vp0 = Sqrt[c33R];
   \[Epsilon]1 = (c22R - c33R)/(2 c33R);
   \[Delta]1 = (c23R^2 - c33R^2)/(2 c33R^2);
   \[Epsilon]2 = (c11R - c33R)/(2 c33R);
   \[Delta]2 = (c13R^2 - c33R^2)/(2 c33R^2);
   \[Delta]3 = (-c11R^2 + c12R^2)/(2 c11R^2);
   Ap0 = iQ33/(1 + Sqrt[1 + iQ33^2]);
   \[Epsilon]1Q = -1 + iQ22/iQ33;
   \[Delta]1Q = (2 c23R^2 (iQ23 - iQ33))/(c33R^2 iQ33);
   \[Epsilon]2Q = -1 + iQ11/iQ33;
   \[Delta]2Q = (2 c13R^2 (iQ13 - iQ33))/(c33R^2 iQ33);
   \[Delta]3Q = -((2 c12R^2 (iQ11 - iQ12))/(c11R^2 iQ11));
   vnmo1 = vp0 Sqrt[1 + 2 \[Delta]1];
   vnmo2 = vp0 Sqrt[1 + 2 \[Delta]2];
   \[Eta]1 = (\[Epsilon]1 - \[Delta]1)/(1 + 2 \[Delta]1);
   \[Eta]2 = (\[Epsilon]2 - \[Delta]2)/(1 + 2 \[Delta]2);
   \[Eta]3 = -((\[Delta]3 - \[Epsilon]1 + \[Epsilon]2 + 
      2 \[Delta]3 \[Epsilon]2)/(
     1 + 2 \[Delta]3 + 2 \[Epsilon]2 + 4 \[Delta]3 \[Epsilon]2));
   
   Return[{vp0, vnmo1, vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
     Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q}];
   ];


Ap2stfORT[vp0_, 
   vs0_, \[Epsilon]1_, \[Delta]1_, \[Gamma]1_, \[Epsilon]2_, \
\[Delta]2_, \[Gamma]2_, \[Delta]3_] := 
  Module[{c11R, c12R, c13R, c22R, c23R, c33R, c44R, c55R, c66R},
   (* convert Thomsen pars of P- and SV-waves to the density-
   normalized stiffness *)
   c11R = vp0^2 (1 + 2 \[Epsilon]2);
   c12R = -vs0^2 (1 + 
        2 \[Gamma]1) + \[Sqrt]((-vs0^2 (1 + 2 \[Gamma]1) + 
          vp0^2 (1 + 2 \[Epsilon]2)) (-vs0^2 (1 + 2 \[Gamma]1) + 
          vp0^2 (1 + 2 \[Delta]3) (1 + 2 \[Epsilon]2)));
   c13R = -vs0^2 + 
     Sqrt[(vp0^2 - vs0^2) (-vs0^2 + vp0^2 (1 + 2 \[Delta]2))];
   c22R = vp0^2 (1 + 2 \[Epsilon]1);
   c23R = -((vs0^2 (1 + 2 \[Gamma]1))/(
      1 + 2 \[Gamma]2)) + \[Sqrt](1/(1 + 
          2 \[Gamma]2)^2 (-vs0^2 (1 + 2 \[Gamma]1) + 
           vp0^2 (1 + 2 \[Gamma]2)) (-vs0^2 (1 + 2 \[Gamma]1) + 
           vp0^2 (1 + 2 \[Gamma]2) (1 + 2 \[Delta]1)));
   c33R = vp0^2;
   c44R = (vs0^2 (1 + 2 \[Gamma]1))/(1 + 2 \[Gamma]2);
   c55R = vs0^2;
   c66R = vs0^2 (1 + 2 \[Gamma]1);
   
   Return[{c11R, c12R, c13R, c22R, c23R, c33R, c44R, c55R, c66R}];
   ];


Ap2stfAORT[vp0_, 
   vs0_, \[Epsilon]1_, \[Delta]1_, \[Gamma]1_, \[Epsilon]2_, \
\[Delta]2_, \[Gamma]2_, \[Delta]3_, Ap0_, 
   As0_, \[Epsilon]1Q_, \[Gamma]1Q_, \[Delta]1Q_, \[Epsilon]2Q_, \
\[Delta]2Q_, \[Gamma]2Q_, \[Delta]3Q_] := 
  Module[{c11R, c12R, c13R, c22R, c23R, c33R, c44R, c55R, c66R, iQ11, 
    iQ12, iQ13, iQ22, iQ23, iQ33, iQ44, iQ55, iQ66, c11, c12, c13, 
    c22, c23, c33, c44, c55, c66},
   c11R = vp0^2 (1 + 2 \[Epsilon]2);
   c12R = -vs0^2 (1 + 
        2 \[Gamma]1) + \[Sqrt]((-vs0^2 (1 + 2 \[Gamma]1) + 
          vp0^2 (1 + 2 \[Epsilon]2)) (-vs0^2 (1 + 2 \[Gamma]1) + 
          vp0^2 (1 + 2 \[Delta]3) (1 + 2 \[Epsilon]2)));
   c13R = -vs0^2 + 
     Sqrt[(vp0^2 - vs0^2) (-vs0^2 + vp0^2 (1 + 2 \[Delta]2))];
   c22R = vp0^2 (1 + 2 \[Epsilon]1);
   c23R = -((vs0^2 (1 + 2 \[Gamma]1))/(
      1 + 2 \[Gamma]2)) + \[Sqrt](1/(1 + 
          2 \[Gamma]2)^2 (-vs0^2 (1 + 2 \[Gamma]1) + 
           vp0^2 (1 + 2 \[Gamma]2)) (-vs0^2 (1 + 2 \[Gamma]1) + 
           vp0^2 (1 + 2 \[Gamma]2) (1 + 2 \[Delta]1)));
   c33R = vp0^2;
   c44R = (vs0^2 (1 + 2 \[Gamma]1))/(1 + 2 \[Gamma]2);
   c55R = vs0^2;
   c66R = vs0^2 (1 + 2 \[Gamma]1);
   
   iQ11 = -((2 Ap0 (1 + \[Epsilon]2Q))/(-1 + Ap0^2));
   iQ12 = (-As0 (c11R + c12R)^2 c66R (1 + \[Gamma]1Q) + 
       Ap0^2 As0 (c11R + c12R)^2 c66R (1 + \[Gamma]1Q) - 
       Ap0 (-1 + As0^2) (-c12R c66R (c12R + 2 c66R) + 
          c11R^3 \[Delta]3Q + c11R^2 (c66R - 2 c66R \[Delta]3Q) + 
          c11R (2 c12R^2 + 4 c12R c66R + 
             c66R^2 \[Delta]3Q)) (1 + \[Epsilon]2Q))/((-1 + 
         Ap0^2) (-1 + As0^2) c12R (c11R - c66R) (c12R + c66R));
   iQ13 = (-As0 (c13R + c33R)^2 c55R + 
       Ap0^2 As0 (c13R + c33R)^2 c55R - 
       Ap0 (-1 + As0^2) (c13R^2 (2 c33R - c55R) + 
          2 c13R (2 c33R - c55R) c55R + 
          c33R (c33R^2 \[Delta]2Q + c55R^2 \[Delta]2Q + 
             c33R (c55R - 2 c55R \[Delta]2Q))))/((-1 + Ap0^2) (-1 + 
         As0^2) c13R (c33R - c55R) (c13R + c55R));
   iQ22 = -((2 Ap0 (1 + \[Epsilon]1Q))/(-1 + Ap0^2));
   iQ23 = (-As0 (c23R + c33R)^2 c44R (1 + \[Gamma]1Q) + 
       Ap0^2 As0 (c23R + c33R)^2 c44R (1 + \[Gamma]1Q) - 
       Ap0 (-1 + As0^2) (1 + \[Gamma]2Q) (c23R^2 (2 c33R - c44R) + 
          2 c23R (2 c33R - c44R) c44R + 
          c33R (c33R^2 \[Delta]1Q + c44R^2 \[Delta]1Q + 
             c33R (c44R - 2 c44R \[Delta]1Q))))/((-1 + Ap0^2) (-1 + 
         As0^2) c23R (c33R - c44R) (c23R + c44R) (1 + \[Gamma]2Q));
   iQ33 = -((2 Ap0)/(-1 + Ap0^2));
   iQ44 = -((
     2 As0 (1 + \[Gamma]1Q))/((-1 + As0^2) (1 + \[Gamma]2Q)));
   iQ55 = -((2 As0)/(-1 + As0^2));
   iQ66 = -((2 As0 (1 + \[Gamma]1Q))/(-1 + As0^2));
   
   c11 = c11R (1 - I iQ11);
   c12 = c12R (1 - I iQ12);
   c13 = c13R (1 - I iQ13);
   c22 = c22R (1 - I iQ22);
   c23 = c23R (1 - I iQ23);
   c33 = c33R (1 - I iQ33);
   c44 = c44R (1 - I iQ44);
   c55 = c55R (1 - I iQ55);
   c66 = c66R (1 - I iQ66);
   
   Return[{c11, c12, c13, c22, c23, c33, c44, c55, c66}];
   ];


Acousticp2stfAORT[vp0_, vn1_, vn2_, \[Eta]1_, \[Eta]2_, \[Eta]3_, 
   Ap0_, \[Epsilon]1Q_, \[Delta]1Q_, \[Epsilon]2Q_, \[Delta]2Q_, \
\[Delta]3Q_] := 
  Module[{c11R, c12R, c13R, c22R, c23R, c33R, iQ11, iQ12, iQ13, iQ22, 
    iQ23, iQ33, c11, c12, c13, c22, c23, c33, \[Xi], kQ},
   \[Xi] = Sqrt[((1 + 2 \[Eta]1) (1 + 2 \[Eta]2))/(1 + 2 \[Eta]3)];
   kQ = Ap0/(1 - Ap0^2);
   (* real part of the complex stiffness coefficients *)
   
   c11R = vn2^2 (1 + 2 \[Eta]2);
   c12R = vn1 vn2 \[Xi];
   c13R = vn2 vp0;
   c22R = vn1^2 (1 + 2 \[Eta]1);
   c23R = vn1 vp0;
   c33R = vp0^2;
   
   iQ11 = 2 kQ (1 + \[Epsilon]2Q);
   iQ12 = 
    kQ (1 + \[Epsilon]2Q) (2 + (\[Delta]3Q (vn2 + 2 vn2 \[Eta]2)^2)/(
       vn1^2 \[Xi]^2));
   iQ13 = kQ (2 + (vp0^2 \[Delta]2Q)/vn2^2);
   iQ22 = 2 kQ (1 + \[Epsilon]1Q);
   iQ23 = kQ (2 + (vp0^2 \[Delta]1Q)/vn1^2);
   iQ33 = 2 kQ;
   
   c11 = c11R (1 - I iQ11);
   c12 = c12R (1 - I iQ12);
   c13 = c13R (1 - I iQ13);
   c22 = c22R (1 - I iQ22);
   c23 = c23R (1 - I iQ23);
   c33 = c33R (1 - I iQ33);
   
   Return[{c11, c12, c13, c22, c23, c33}];
   ];


KelvinVoigt[f_, f0_, vp0_, vnmo1_, 
   vnmo2_, \[Eta]1_, \[Eta]2_, \[Eta]3_, 
   Ap0_, \[Epsilon]1Q_, \[Delta]1Q_, \[Epsilon]2Q_, \[Delta]2Q_, \
\[Delta]3Q_] := 
  Module[{\[Xi], kQ, c11R, c12R, c13R, c22R, c23R, c33R, iQ11, iQ12, 
    iQ13, iQ22, iQ23, iQ33,
    m11, m12, m13, m22, m23, 
    m33, \[Eta]11, \[Eta]12, \[Eta]13, \[Eta]22, \[Eta]23, \[Eta]33, \
\[Omega], c11, c12, c13, c22, c23, c33,
    vp0f, vnmo1f, vnmo2f, \[Eta]1f, \[Eta]2f, \[Eta]3f, 
    Ap0f, \[Epsilon]1Qf, \[Delta]1Qf, \[Epsilon]2Qf, \[Delta]2Qf, \
\[Delta]3Qf},
   \[Xi] = Sqrt[((1 + 2 \[Eta]1) (1 + 2 \[Eta]2))/(1 + 2 \[Eta]3)];
   kQ = Ap0/(1 - Ap0^2);
   (* real part of the complex stiffness coefficients *)
   
   c11R = vnmo2^2 (1 + 2 \[Eta]2);
   c12R = vnmo1 vnmo2 \[Xi];
   c13R = vnmo2 vp0;
   c22R = vnmo1^2 (1 + 2 \[Eta]1);
   c23R = vnmo1 vp0;
   c33R = vp0^2;
   iQ11 = 2 kQ (1 + \[Epsilon]2Q);
   iQ12 = 
    kQ (1 + \[Epsilon]2Q) (2 + (\[Delta]3Q (vnmo2 + 
          2 vnmo2 \[Eta]2)^2)/(vnmo1^2 \[Xi]^2));
   iQ13 = kQ (2 + (vp0^2 \[Delta]2Q)/vnmo2^2);
   iQ22 = 2 kQ (1 + \[Epsilon]1Q);
   iQ23 = kQ (2 + (vp0^2 \[Delta]1Q)/vnmo1^2);
   iQ33 = 2 kQ;
   
   m11 = c11R;
   m12 = c12R;
   m13 = c13R;
   m22 = c22R;
   m23 = c23R;
   m33 = c33R;
   \[Eta]11 = 1/(2 Pi f0) c11R iQ11;
   \[Eta]12 = 1/(2 Pi f0) c12R iQ12;
   \[Eta]13 = 1/(2 Pi f0) c13R iQ13;
   \[Eta]22 = 1/(2 Pi f0) c22R iQ22;
   \[Eta]23 = 1/(2 Pi f0) c23R iQ23;
   \[Eta]33 = 1/(2 Pi f0) c33R iQ33;
   
   \[Omega] = 2 Pi f;
   
   c11 = m11 - I \[Omega] \[Eta]11;
   c12 = m12 - I \[Omega] \[Eta]12;
   c13 = m13 - I \[Omega] \[Eta]13;
   c22 = m22 - I \[Omega] \[Eta]22;
   c23 = m23 - I \[Omega] \[Eta]23;
   c33 = m33 - I \[Omega] \[Eta]33;
   
   Return[{c11, c12, c13, c22, c23, c33}];
   ];


Maxwell[f_, f0_, vp0_, vnmo1_, vnmo2_, \[Eta]1_, \[Eta]2_, \[Eta]3_, 
   Ap0_, \[Epsilon]1Q_, \[Delta]1Q_, \[Epsilon]2Q_, \[Delta]2Q_, \
\[Delta]3Q_] := 
  Module[{\[Xi], kQ, c11R, c12R, c13R, c22R, c23R, c33R, iQ11, iQ12, 
    iQ13, iQ22, iQ23, iQ33, m11, m12, m13, m22, m23, 
    m33, \[Eta]11, \[Eta]12, \[Eta]13, \[Eta]22, \[Eta]23, \[Eta]33, \
\[Omega], c11, c12, c13, c22, c23, c33, vp0f, vnmo1f, 
    vnmo2f, \[Eta]1f, \[Eta]2f, \[Eta]3f, 
    Ap0f, \[Epsilon]1Qf, \[Delta]1Qf, \[Epsilon]2Qf, \[Delta]2Qf, \
\[Delta]3Qf},
   \[Xi] = Sqrt[((1 + 2 \[Eta]1) (1 + 2 \[Eta]2))/(1 + 2 \[Eta]3)];
   kQ = Ap0/(1 - Ap0^2);
   (* real part of the complex stiffness coefficients *)
   
   c11R = vnmo2^2 (1 + 2 \[Eta]2);
   c12R = vnmo1 vnmo2 \[Xi];
   c13R = vnmo2 vp0;
   c22R = vnmo1^2 (1 + 2 \[Eta]1);
   c23R = vnmo1 vp0;
   c33R = vp0^2;
   
   iQ11 = 2 kQ (1 + \[Epsilon]2Q);
   iQ12 = 
    kQ (1 + \[Epsilon]2Q) (2 + (\[Delta]3Q (vnmo2 + 
          2 vnmo2 \[Eta]2)^2)/(vnmo1^2 \[Xi]^2));
   iQ13 = kQ (2 + (vp0^2 \[Delta]2Q)/vnmo2^2);
   iQ22 = 2 kQ (1 + \[Epsilon]1Q);
   iQ23 = kQ (2 + (vp0^2 \[Delta]1Q)/vnmo1^2);
   iQ33 = 2 kQ;
   
   m11 = c11R (1 + iQ11^2);
   m12 = c12R (1 + iQ12^2);
   m13 = c13R (1 + iQ13^2);
   m22 = c22R (1 + iQ22^2);
   m23 = c23R (1 + iQ23^2);
   m33 = c33R (1 + iQ33^2);
   \[Eta]11 = 1/(2 Pi f0) c11R (1/iQ11 + iQ11);
   \[Eta]12 = 1/(2 Pi f0) c12R (1/iQ12 + iQ12);
   \[Eta]13 = 1/(2 Pi f0) c13R (1/iQ13 + iQ13);
   \[Eta]22 = 1/(2 Pi f0) c22R (1/iQ22 + iQ22);
   \[Eta]23 = 1/(2 Pi f0) c23R (1/iQ23 + iQ23);
   \[Eta]33 = 1/(2 Pi f0) c33R (1/iQ33 + iQ33);
   
   \[Omega] = 2 Pi f;
   
   c11 = (1/m11 +I/(\[Omega] \[Eta]11))^-1;
   c12 = (1/m12 +I/(\[Omega] \[Eta]12))^-1;
   c13 = (1/m13 +I/(\[Omega] \[Eta]13))^-1;
   c22 = (1/m22 +I/(\[Omega] \[Eta]22))^-1;
   c23 = (1/m23 +I/(\[Omega] \[Eta]23))^-1;
   c33 = (1/m33 +I/(\[Omega] \[Eta]33))^-1;
   
   Return[{c11, c12, c13, c22, c23, c33}];
   ];


StandardLinear[f_, f0_, vp0_, vnmo1_, 
   vnmo2_, \[Eta]1_, \[Eta]2_, \[Eta]3_, 
   Ap0_, \[Epsilon]1Q_, \[Delta]1Q_, \[Epsilon]2Q_, \[Delta]2Q_, \
\[Delta]3Q_] := 
  Module[{\[Xi], kQ, c11R, c12R, c13R, c22R, c23R, c33R, iQ11, iQ12, 
    iQ13, iQ22, iQ23, iQ33, n11, n12, n13, n22, n23, 
    n33, \[Tau]11, \[Tau]12, \[Tau]13, \[Tau]22, \[Tau]23, \[Tau]33, \
\[Tau]\[Sigma], \[Omega], c11, c12, c13, c22, c23, c33, vp0f, vnmo1f, 
    vnmo2f, \[Eta]1f, \[Eta]2f, \[Eta]3f, 
    Ap0f, \[Epsilon]1Qf, \[Delta]1Qf, \[Epsilon]2Qf, \[Delta]2Qf, \
\[Delta]3Qf},
   \[Xi] = Sqrt[((1 + 2 \[Eta]1) (1 + 2 \[Eta]2))/(1 + 2 \[Eta]3)];
   kQ = Ap0/(1 - Ap0^2);
   (* real part of the complex stiffness coefficients *)
   
   c11R = vnmo2^2 (1 + 2 \[Eta]2);
   c12R = vnmo1 vnmo2 \[Xi];
   c13R = vnmo2 vp0;
   c22R = vnmo1^2 (1 + 2 \[Eta]1);
   c23R = vnmo1 vp0;
   c33R = vp0^2;
   
   iQ11 = 2 kQ (1 + \[Epsilon]2Q);
   iQ12 = 
    kQ (1 + \[Epsilon]2Q) (2 + (\[Delta]3Q (vnmo2 + 
          2 vnmo2 \[Eta]2)^2)/(vnmo1^2 \[Xi]^2));
   iQ13 = kQ (2 + (vp0^2 \[Delta]2Q)/vnmo2^2);
   iQ22 = 2 kQ (1 + \[Epsilon]1Q);
   iQ23 = kQ (2 + (vp0^2 \[Delta]1Q)/vnmo1^2);
   iQ33 = 2 kQ;
   
   n11 = c11R (1 - iQ11);
   n12 = c12R (1 - iQ12);
   n13 = c13R (1 - iQ13);
   n22 = c22R (1 - iQ22);
   n23 = c23R (1 - iQ23);
   n33 = c33R (1 - iQ33);
   \[Tau]\[Sigma] = 1/(2 Pi f0);
   \[Tau]11 = \[Tau]\[Sigma] ((1 + iQ11)/(1 - iQ11));
   \[Tau]12 = \[Tau]\[Sigma]  ((1 + iQ12)/(1 - iQ12));
   \[Tau]13 = \[Tau]\[Sigma] ((1 + iQ13)/(1 - iQ13));
   \[Tau]22 = \[Tau]\[Sigma]  ((1 + iQ22)/(1 - iQ22));
   \[Tau]23 = \[Tau]\[Sigma] ((1 + iQ23)/(1 - iQ23));
   \[Tau]33 = \[Tau]\[Sigma]  ((1 + iQ33)/(1 - iQ33));
   
   \[Omega] = 2 Pi f;
   
   c11 = n11 ((1 - I \[Omega] \[Tau]11)/(
      1 - I \[Omega] \[Tau]\[Sigma]));
   c12 = n12 ((1 - I \[Omega] \[Tau]12)/(
      1 - I \[Omega] \[Tau]\[Sigma]));
   c13 = n13 ((1 - I \[Omega] \[Tau]13)/(
      1 - I \[Omega] \[Tau]\[Sigma]));
   c22 = n22 ((1 - I \[Omega] \[Tau]22)/(
      1 - I \[Omega] \[Tau]\[Sigma]));
   c23 = n23 ((1 - I \[Omega] \[Tau]23)/(
      1 - I \[Omega] \[Tau]\[Sigma]));
   c33 = n33 ((1 - I \[Omega] \[Tau]33)/(
      1 - I \[Omega] \[Tau]\[Sigma]));
   
   Return[{c11, c12, c13, c22, c23, c33}];
   ];


Kjartansson[f_, f0_, vp0_, vnmo1_, 
   vnmo2_, \[Eta]1_, \[Eta]2_, \[Eta]3_, 
   Ap0_, \[Epsilon]1Q_, \[Delta]1Q_, \[Epsilon]2Q_, \[Delta]2Q_, \
\[Delta]3Q_] := 
  Module[{\[Xi], kQ, c11R, c12R, c13R, c22R, c23R, c33R, iQ11, iQ12, 
    iQ13, iQ22, iQ23, iQ33, m11, m12, m13, m22, m23, 
    m33, \[Gamma]11, \[Gamma]12, \[Gamma]13, \[Gamma]22, \[Gamma]23, \
\[Gamma]33, \[Tau]\[Sigma], \[Omega], c11, c12, c13, c22, c23, c33, 
    vp0f, vnmo1f, vnmo2f, \[Eta]1f, \[Eta]2f, \[Eta]3f, 
    Ap0f, \[Epsilon]1Qf, \[Delta]1Qf, \[Epsilon]2Qf, \[Delta]2Qf, \
\[Delta]3Qf},
   \[Xi] = Sqrt[((1 + 2 \[Eta]1) (1 + 2 \[Eta]2))/(1 + 2 \[Eta]3)];
   kQ = Ap0/(1 - Ap0^2);
   (* real part of the complex stiffness coefficients *)
   
   c11R = vnmo2^2 (1 + 2 \[Eta]2);
   c12R = vnmo1 vnmo2 \[Xi];
   c13R = vnmo2 vp0;
   c22R = vnmo1^2 (1 + 2 \[Eta]1);
   c23R = vnmo1 vp0;
   c33R = vp0^2;
   
   iQ11 = 2 kQ (1 + \[Epsilon]2Q);
   iQ12 = 
    kQ (1 + \[Epsilon]2Q) (2 + (\[Delta]3Q (vnmo2 + 
          2 vnmo2 \[Eta]2)^2)/(vnmo1^2 \[Xi]^2));
   iQ13 = kQ (2 + (vp0^2 \[Delta]2Q)/vnmo2^2);
   iQ22 = 2 kQ (1 + \[Epsilon]1Q);
   iQ23 = kQ (2 + (vp0^2 \[Delta]1Q)/vnmo1^2);
   iQ33 = 2 kQ;
   
   \[Gamma]11 = 1/Pi ArcTan[iQ11];
   \[Gamma]12 = 1/Pi ArcTan[iQ12];
   \[Gamma]13 = 1/Pi ArcTan[iQ13];
   \[Gamma]22 = 1/Pi ArcTan[iQ22];
   \[Gamma]23 = 1/Pi ArcTan[iQ23];
   \[Gamma]33 = 1/Pi ArcTan[iQ33];
   
   m11 = c11R/Cos[Pi \[Gamma]11];
   m12 = c12R/Cos[Pi \[Gamma]12];
   m13 = c13R/Cos[Pi \[Gamma]13];
   m22 = c22R/Cos[Pi \[Gamma]22];
   m23 = c23R/Cos[Pi \[Gamma]23];
   m33 = c33R/Cos[Pi \[Gamma]33];
   \[Tau]\[Sigma] = 1/(2 Pi f0);
   
   \[Omega] = 2 Pi f;
   
   c11 = m11 (((-I \[Omega] )/(2 Pi f0)))^(2 \[Gamma]11);
   c12 = m12 ((-I \[Omega] )/(2 Pi f0))^(2 \[Gamma]12);
   c13 = m13 ((-I \[Omega] )/(2 Pi f0))^(2 \[Gamma]13);
   c22 = m22 ((-I \[Omega] )/(2 Pi f0))^(2 \[Gamma]22);
   c23 = m23 ((-I \[Omega] )/(2 Pi f0))^(2 \[Gamma]23);
   c33 = m33 ((-I \[Omega] )/(2 Pi f0))^(2 \[Gamma]33);
   
   Return[{c11, c12, c13, c22, c23, c33}];
   ];


VrayORT[N1_, N2_, N3_, Stiffpars_] := 
  Module[{c11, c12, c13, c22, c23, c33, c44, c55, c66, slowness, p1, 
    p2, p3, p10, p20, p30, eq1, eq2, eq3, sol, px, py, 
    pz, \[CapitalGamma]11 , \[CapitalGamma]12, \[CapitalGamma]13, \
\[CapitalGamma]22, \[CapitalGamma]23, \[CapitalGamma]33, D11, D22, 
    D33, D12, D13, D21, D23, D31, D32, DD, pN, C1ijklDjknl, 
    C2ijklDjknl, C3ijklDjknl, g1g1, g1g2, g1g3, g2g2, g2g3, g3g3, Vx, 
    Vy, Vz, Vray},
   (* complex ray velocity as a function of ray direction (N1,N2,
   N3) *)
   
   c11 = Stiffpars[[1]];
   c12 = Stiffpars[[2]];
   c13 = Stiffpars[[3]];
   c22 = Stiffpars[[4]];
   c23 = Stiffpars[[5]];
   c33 = Stiffpars[[6]];
   c44 = Stiffpars[[7]];
   c55 = Stiffpars[[8]];
   c66 = Stiffpars[[9]];
   
   (* slowness of P-waves in the elastic, 
   elliptically VTI background *)
   
   slowness = Sqrt[N1^2/Re[c11] + N2^2/Re[c22] + N3^2/Re[c33]];
   
   p10 = slowness N1;
   p20 = slowness N2;
          p30 = slowness N3;
   
   \[CapitalGamma]11 = c11 p1^2 + c66 p2^2 + c55 p3^2;
   \[CapitalGamma]12 = (c12 + c66) p1 p2;
   \[CapitalGamma]22 = c66 p1^2 + c22 p2^2 + c44 p3^2;
   \[CapitalGamma]13 = (c13 + c55) p1 p3;
   \[CapitalGamma]33 = c55 p1^2 + c44 p2^2 + c33 p3^2;
   \[CapitalGamma]23 = (c23 + c44) p2 p3;
   
   D11 = (\[CapitalGamma]22 - 1) (\[CapitalGamma]33 - 
        1) - \[CapitalGamma]23^2;
   D22 = (\[CapitalGamma]11 - 1) (\[CapitalGamma]33 - 
        1) - \[CapitalGamma]13^2;
   D33 = (\[CapitalGamma]11 - 1) (\[CapitalGamma]22 - 
        1) - \[CapitalGamma]12^2;
   D12 = D21 = \[CapitalGamma]13 \[CapitalGamma]23 - \
\[CapitalGamma]12 (\[CapitalGamma]33 - 1);
   D13 = D31 = \[CapitalGamma]12 \[CapitalGamma]23 - \
\[CapitalGamma]13 (\[CapitalGamma]22 - 1);
   D23 = D32 = \[CapitalGamma]12 \[CapitalGamma]13 - \
\[CapitalGamma]23 (\[CapitalGamma]11 - 1);
   DD = D11 + D22 + D33;
   
   C1ijklDjknl = 
    c11 p1 D11 + c12 p2 D12 + c13 p3 D13 + c66 p2 D21 + c66 p1 D22 + 
     c55 p3 D31 + c55 p1 D33;
   C2ijklDjknl = 
    c66 p2 D11 + c66 p1 D12 + c12 p1 D21 + c22 p2 D22 + c23 p3 D23 + 
     c44 p3 D32 + c44 p2 D33;
   C3ijklDjknl = 
    c55 p3 D11 + c55 p1 D13 + c44 p3 D22 + c44 p2 D23 + c13 p1 D31 + 
     c23 p2 D32 + c33 p3 D33;
   pN = p1 N1 + p2 N2 + p3 N3;
   eq1 = DD N1 - C1ijklDjknl pN;
   eq2 = DD N2 - C2ijklDjknl pN;
   eq3 = DD N3 - C3ijklDjknl pN;
   
   sol = FindRoot[{eq1 == 0, eq2 == 0, 
      eq3 == 0}, {{p1, p10}, {p2, p20}, {p3, p30}}];
   
   px = p1 /. sol;
   py = p2 /. sol;
   pz = p3 /. sol;
   
   g1g1 = (D11/DD) /. sol;
   g2g2 = (D22/DD) /. sol;
   g3g3 = (D33/DD) /. sol;
   g1g2 = (D12/DD) /. sol;
   g1g3 = (D13/DD) /. sol;
   g2g3 = (D23/DD) /. sol;
   
   (* Ray velocity *)
   
   Vx = px c11 g1g1 + px c66 g2g2 + px c55 g3g3 + 
     py (c12 + c66) g1g2 + pz (c13 + c55) g1g3;
   Vy = py c66 g1g1 + py c22 g2g2 + py c44 g3g3 + 
     px (c12 + c66) g1g2 + pz (c23 + c44) g2g3;
   Vz = pz c55 g1g1 + pz c44 g2g2 + pz c33 g3g3 + 
     px (c13 + c55) g1g3 + py (c23 + c44) g2g3;
   
   Vray = Sqrt[Vx^2 + Vy^2 + Vz^2];
   
   Return[Vray];
   ];


VrayAray[N1_, N2_, N3_, Stiffpars_] := 
  Module[{Vgroup, ReV, ImV, Vray, Aray},
   (* Ray velocity and ray attenuation as a function of ray direction \
*)
   
   Vgroup = VrayORT[N1, N2, N3, Stiffpars];
   
   ReV = Re[Vgroup];
   ImV = Im[Vgroup];
   
   Vray = (ReV^2 + ImV^2)/ReV;
   Aray = - (ImV/(ReV^2 + ImV^2));
   
   Return[{Vray, Aray}];
   
   ];


viscoGreenQuantities[\[Alpha]_, \[Beta]_, Stiffpars_, freq_] := 
  Module[{N1, N2, N3, c11, c12, c13, c22, c23, c33, c44, c55, c66, 
    slowness, p1, p2, p3, p10, p20, p30, eq1, eq2, eq3, sol, px, py, 
    pz, \[CapitalGamma]11 , \[CapitalGamma]12, \[CapitalGamma]13, \
\[CapitalGamma]22, \[CapitalGamma]23, \[CapitalGamma]33, D11, D22, 
    D33, D12, D13, D21, D23, D31, D32, DD, pN, C1ijklDjknl, 
    C2ijklDjknl, C3ijklDjknl, vp1, vp2, vp3, M11, M12, M21, M13, M31, 
    M22, M23, M32, M33, dVEdvp1, dVEdvp2, dVEdvp3, dEdp3, 
    d2VEdvp1dvp1, d2VEdvp1dvp2, d2VEdvp2dvp1, d2VEdvp1dvp3, 
    d2VEdvp3dvp1, d2VEdvp2dvp2, d2VEdvp2dvp3, d2VEdvp3dvp2, 
    d2VEdvp3dvp3, d2Edp1dp1, d2Edp1dp2, d2Edp2dp1, d2Edp2dp2, 
    d2p3dp1dp1, d2p3dp1dp2, d2p3dp2dp1, d2p3dp2dp2, K, absK, lamda, 
    vr, vi, Vray, Aray, absdEdp3, absSqrtKdEdp3},
   (* complex ray velocity as a function of ray direction (N1,N2,
   N3) *)
   (* 
   Gaussian Curvature as a function of the ray direction, \[Alpha]:
   polar angle, \[Beta]:azimuthal angle *)
   
   N1 = Cos[\[Beta]] Sin[\[Alpha]];
   N2 = Sin[\[Beta]] Sin[\[Alpha]];
   N3 = Cos[\[Alpha]];
   
   c11 = Stiffpars[[1]];
   c12 = Stiffpars[[2]];
   c13 = Stiffpars[[3]];
   c22 = Stiffpars[[4]];
   c23 = Stiffpars[[5]];
   c33 = Stiffpars[[6]];
   c44 = 0;
   c55 = 0;
   c66 = 0;
   
   (* slowness of P-waves in the elastic, 
   elliptically VTI background *)
   
   slowness = Sqrt[N1^2/Re[c11] + N2^2/Re[c22] + N3^2/Re[c33]];
   
   p10 = slowness N1;
   p20 = slowness N2;
   p30 = slowness N3;
   
   \[CapitalGamma]11 = c11 p1^2 + c66 p2^2 + c55 p3^2;
   \[CapitalGamma]12 = (c12 + c66) p1 p2;
   \[CapitalGamma]22 = c66 p1^2 + c22 p2^2 + c44 p3^2;
   \[CapitalGamma]13 = (c13 + c55) p1 p3;
   \[CapitalGamma]33 = c55 p1^2 + c44 p2^2 + c33 p3^2;
   \[CapitalGamma]23 = (c23 + c44) p2 p3;
   
   D11 = (\[CapitalGamma]22 - 1) (\[CapitalGamma]33 - 
        1) - \[CapitalGamma]23^2;
   D22 = (\[CapitalGamma]11 - 1) (\[CapitalGamma]33 - 
        1) - \[CapitalGamma]13^2;
   D33 = (\[CapitalGamma]11 - 1) (\[CapitalGamma]22 - 
        1) - \[CapitalGamma]12^2;
   D12 = D21 = \[CapitalGamma]13 \[CapitalGamma]23 - \
\[CapitalGamma]12 (\[CapitalGamma]33 - 1);
   D13 = D31 = \[CapitalGamma]12 \[CapitalGamma]23 - \
\[CapitalGamma]13 (\[CapitalGamma]22 - 1);
   D23 = D32 = \[CapitalGamma]12 \[CapitalGamma]13 - \
\[CapitalGamma]23 (\[CapitalGamma]11 - 1);
   DD = D11 + D22 + D33;
   
   C1ijklDjknl = 
    c11 p1 D11 + c12 p2 D12 + c13 p3 D13 + c66 p2 D21 + c66 p1 D22 + 
     c55 p3 D31 + c55 p1 D33;
   C2ijklDjknl = 
    c66 p2 D11 + c66 p1 D12 + c12 p1 D21 + c22 p2 D22 + c23 p3 D23 + 
     c44 p3 D32 + c44 p2 D33;
   C3ijklDjknl = 
    c55 p3 D11 + c55 p1 D13 + c44 p3 D22 + c44 p2 D23 + c13 p1 D31 + 
     c23 p2 D32 + c33 p3 D33;
   pN = p1 N1 + p2 N2 + p3 N3;
   eq1 = DD N1 - C1ijklDjknl pN;
   eq2 = DD N2 - C2ijklDjknl pN;
   eq3 = DD N3 - C3ijklDjknl pN;
   
   sol = FindRoot[{eq1 == 0, eq2 == 0, 
      eq3 == 0}, {{p1, p10}, {p2, p20}, {p3, p30}}];
   
   vp1 = p1 /. sol;
   vp2 = p2 /. sol;
   vp3 = p3 /. sol;
   
   (* - slowness components in the rotated coordiante system ---- *)
 \
  ({
      {p1},
      {p2},
      {p3}
     }) = \!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{
RowBox[{"Cos", "[", "\[Alpha]", "]"}], " ", 
RowBox[{"Cos", "[", "\[Beta]", "]"}]}], 
RowBox[{
RowBox[{"Cos", "[", "\[Alpha]", "]"}], " ", 
RowBox[{"Sin", "[", "\[Beta]", "]"}]}], 
RowBox[{"-", 
RowBox[{"Sin", "[", "\[Alpha]", "]"}]}]},
{
RowBox[{"-", 
RowBox[{"Sin", "[", "\[Beta]", "]"}]}], 
RowBox[{"Cos", "[", "\[Beta]", "]"}], "0"},
{
RowBox[{
RowBox[{"Cos", "[", "\[Beta]", "]"}], " ", 
RowBox[{"Sin", "[", "\[Alpha]", "]"}]}], 
RowBox[{
RowBox[{"Sin", "[", "\[Alpha]", "]"}], " ", 
RowBox[{"Sin", "[", "\[Beta]", "]"}]}], 
RowBox[{"Cos", "[", "\[Alpha]", "]"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\).({
       {vp1},
       {vp2},
       {vp3}
      });
   
   (* ======= output pars: ray velocity and attenuation ======= *)
   
   vr = Re[1/p3];
   vi = Im[1/p3];
   Vray = (vr^2 + vi^2)/vr;
   Aray = -(vi/(vr^2 + vi^2));
   
   (* ------------- dEdp3 ------------------------ *)
   ({
      {M11, M12, M13},
      {M21, M22, M23},
      {M31, M32, M33}
     }) = ({
      {Cos[\[Beta]] Cos[\[Alpha]], -Sin[\[Beta]], 
       Cos[\[Beta]] Sin[\[Alpha]]},
      {Sin[\[Beta]] Cos[\[Alpha]], Cos[\[Beta]], 
       Sin[\[Beta]] Sin[\[Alpha]]},
      {-Sin[\[Alpha]], 0, Cos[\[Alpha]]}
     });
   
   dVEdvp1 = 
    2 vp1 (2 c12 c13 c23 vp2^2 vp3^2 + c13^2 (1 - c22 vp2^2) vp3^2 + 
       c12^2 vp2^2 (1 - c33 vp3^2) + 
       c11 (1 - c33 vp3^2 - c23^2 vp2^2 vp3^2 + 
          c22 vp2^2 (-1 + c33 vp3^2)));
   dVEdvp2 = 
    2 vp2 (2 c12 c13 c23 vp1^2 vp3^2 + c23^2 (1 - c11 vp1^2) vp3^2 + 
       c12^2 vp1^2 (1 - c33 vp3^2) + 
       c22 (1 - c33 vp3^2 - c13^2 vp1^2 vp3^2 + 
          c11 vp1^2 (-1 + c33 vp3^2)));
   dVEdvp3 = 
    2 (2 c12 c13 c23 vp1^2 vp2^2 + c23^2 (1 - c11 vp1^2) vp2^2 + 
       c13^2 vp1^2 (1 - c22 vp2^2) + 
       c33 (1 - c22 vp2^2 - c12^2 vp1^2 vp2^2 + 
          c11 vp1^2 (-1 + c22 vp2^2))) vp3;
   (* ======= output pars: dI'dp3' *)
   dEdp3 = (({
          {dVEdvp1, dVEdvp2, dVEdvp3}
         }).({
          {M13},
          {M23},
          {M33}
         }))[[1]][[1]];
   (* ---------------d2p3dpidpj-------------------------------------*)

      d2VEdvp1dvp1 = -2 c23 (-c12 c13 + c11 c23) vp2^2 vp3^2 + 
     2 c13 (c13 - c13 c22 vp2^2 + c12 c23 vp2^2) vp3^2 - 
     2 (c11 + c12^2 vp2^2 - c11 c22 vp2^2) (-1 + c33 vp3^2);
   d2VEdvp1dvp2 = 
    d2VEdvp2dvp1 = -4 vp1 vp2 (c13^2 c22 vp3^2 - 2 c12 c13 c23 vp3^2 +
         c12^2 (-1 + c33 vp3^2) + 
        c11 (c22 + c23^2 vp3^2 - c22 c33 vp3^2));
   d2VEdvp1dvp3 = 
    d2VEdvp3dvp1 = -4 vp1 (-2 c12 c13 c23 vp2^2 + c12^2 c33 vp2^2 + 
        c13^2 (-1 + c22 vp2^2) + 
        c11 (c33 + c23^2 vp2^2 - c22 c33 vp2^2)) vp3;
   d2VEdvp2dvp2 = 
    2 c13 (-c13 c22 + c12 c23) vp1^2 vp3^2 + 
     2 c23 (c23 + c12 c13 vp1^2 - c11 c23 vp1^2) vp3^2 - 
     2 (c22 + c12^2 vp1^2 - c11 c22 vp1^2) (-1 + c33 vp3^2);
   d2VEdvp2dvp3 = 
    d2VEdvp3dvp2 = -4 (-2 c12 c13 c23 vp1^2 + c12^2 c33 vp1^2 + 
        c23^2 (-1 + c11 vp1^2) + 
        c22 (c33 + c13^2 vp1^2 - c11 c33 vp1^2)) vp2 vp3;
   d2VEdvp3dvp3 = 
    2 (c23 (c23 + c12 c13 vp1^2 - c11 c23 vp1^2) vp2^2 + 
       c13 vp1^2 (c13 - c13 c22 vp2^2 + c12 c23 vp2^2) + 
       c33 (1 - c22 vp2^2 - c12^2 vp1^2 vp2^2 + 
          c11 vp1^2 (-1 + c22 vp2^2)));
   d2Edp1dp1 = (({
          {M11, M21, M31}
         }).({
          {d2VEdvp1dvp1, d2VEdvp1dvp2, d2VEdvp1dvp3},
          {d2VEdvp2dvp1, d2VEdvp2dvp2, d2VEdvp2dvp3},
          {d2VEdvp3dvp1, d2VEdvp3dvp2, d2VEdvp3dvp3}
         }).({
          {M11},
          {M21},
          {M31}
         }))[[1]][[1]];
   d2Edp1dp2 = d2Edp2dp1 = (({
           {M11, M21, M31}
          }).({
           {d2VEdvp1dvp1, d2VEdvp1dvp2, d2VEdvp1dvp3},
           {d2VEdvp2dvp1, d2VEdvp2dvp2, d2VEdvp2dvp3},
           {d2VEdvp3dvp1, d2VEdvp3dvp2, d2VEdvp3dvp3}
          }).({
           {M12},
           {M22},
           {M32}
          }))[[1]][[1]];
   d2Edp2dp2 = (({
          {M12, M22, M32}
         }).({
          {d2VEdvp1dvp1, d2VEdvp1dvp2, d2VEdvp1dvp3},
          {d2VEdvp2dvp1, d2VEdvp2dvp2, d2VEdvp2dvp3},
          {d2VEdvp3dvp1, d2VEdvp3dvp2, d2VEdvp3dvp3}
         }).({
          {M12},
          {M22},
          {M32}
         }))[[1]][[1]];
   
   d2p3dp1dp1 = d2Edp1dp1/dEdp3;
   d2p3dp1dp2 = d2p3dp2dp1 = d2Edp1dp2/dEdp3;
   d2p3dp2dp2 = d2Edp2dp2/dEdp3;
   
   (* ======\[Equal] Gaussian Curvature K=====================*)
   
   K = (d2p3dp1dp1 d2p3dp2dp2 - d2p3dp1dp2 d2p3dp2dp1);
   
   absK = Abs[K];
   
   lamda = -(1/2) Arg[K] - Arg[dEdp3];
   absdEdp3 = Abs[dEdp3];
   
   absSqrtKdEdp3 = Sqrt[absK] absdEdp3;
   
   Return[{absK, lamda, Vray, Aray, absdEdp3, absSqrtKdEdp3}];
   ];


viscoKVdEdp3[\[Alpha]_, \[Beta]_, f_, f0_, AcoustPars0_, opt_] := 
  Module[{vp0, vnmo1, vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
    Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q, c110, c120, c130, c220, c230, c330, c44B, c55B, c66B, 
    Stiffpars, K, absK, lamda, vr, vi, Vray, Aray, absdEdp3, p3, 
    absSqrtKdEdp3},
   
   vp0 = AcoustPars0[[1]];
   vnmo1 = AcoustPars0[[2]];
   vnmo2 = AcoustPars0[[3]];
   \[Eta]1 = AcoustPars0[[4]];
   \[Eta]2 = AcoustPars0[[5]];
   \[Eta]3 = AcoustPars0[[6]];
   Ap0 = AcoustPars0[[7]];
   \[Epsilon]1Q = AcoustPars0[[8]];
   \[Delta]1Q = AcoustPars0[[9]];
   \[Epsilon]2Q = AcoustPars0[[10]];
   \[Delta]2Q = AcoustPars0[[11]];
   \[Delta]3Q = AcoustPars0[[12]];
   
   If[opt == 1, {
     {c110, c120, c130, c220, c230, c330} = 
       KelvinVoigt[f, f0, vp0, vnmo1, 
        vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     }];
   If[opt == 2, {
     {c110, c120, c130, c220, c230, c330} = 
       Maxwell[f, f0, vp0, vnmo1, vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     }];
   If[opt == 3, {
     {c110, c120, c130, c220, c230, c330} = 
       StandardLinear[f, f0, vp0, vnmo1, 
        vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     }];
   If[opt == 4, {
     {c110, c120, c130, c220, c230, c330} = 
       Kjartansson[f, f0, vp0, vnmo1, 
        vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     }];
   
   Which[opt == 1, {
     {c110, c120, c130, c220, c230, c330} = 
       KelvinVoigt[f, f0, vp0, vnmo1, 
        vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     },
    opt == 2, {
     {c110, c120, c130, c220, c230, c330} = 
       Maxwell[f, f0, vp0, vnmo1, vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     },
    opt == 3, {
     {c110, c120, c130, c220, c230, c330} = 
       StandardLinear[f, f0, vp0, vnmo1, 
        vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     },
    opt == 4, {
     {c110, c120, c130, c220, c230, c330} = 
       Kjartansson[f, f0, vp0, vnmo1, 
        vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
        Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
     },
    True, {
     If[f > 0,
      {c110, c120, c130, c220, c230, c330} = 
        Acousticp2stfAORT[vp0, vnmo1, 
         vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, 
         Ap0, \[Epsilon]1Q, \[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \
\[Delta]3Q];
      ];
     If[f < 0,
      {c110, c120, c130, c220, c230, c330} = 
        Acousticp2stfAORT[vp0, vnmo1, 
         vnmo2, \[Eta]1, \[Eta]2, \[Eta]3, -Ap0, \[Epsilon]1Q, \
\[Delta]1Q, \[Epsilon]2Q, \[Delta]2Q, \[Delta]3Q];
      ];
     }
    ];
   
   Stiffpars = {c110, c120, c130, c220, c230, c330};
   
   {absK, lamda, Vray, Aray, absdEdp3, absSqrtKdEdp3} = 
    viscoGreenQuantities[\[Alpha], \[Beta], Stiffpars, f];
   Return[{absK, lamda, Vray, Aray, absdEdp3, absSqrtKdEdp3}]
   ];


WaveletSource[fm_, t_] := Module[{ft, t0, ft0},
   t0 = 1/fm;
   
   ft = (1 - 2 Pi^2 fm^2 (t - t0)^2) Exp[-Pi^2 fm^2 (t - t0)^2];
   
   Return[ft];
   ];


Synthesis[\[Alpha]_,\[Beta]_,R_,fc_,np_,dt_,AcoustPars_,opt_]:=Module[{f,df,tdata,fdata,ip,absK,lamda,Vray,Aray,absdIdp3,p3,omega,tdata0},
(* ploar angle: \[Alpha], azimuth: \[Beta] *)
df=1/(np dt);

tdata=Table[WaveletSource[fc,(it-1)*dt],{it,1,np}];
fdata=Fourier[tdata,FourierParameters->{-1,1}];

(* --------- nonzero frequency --------- *)
For[ip=2,ip<=np,ip++,{
If[ip<=np/2+1,{
f=(ip-1)*df;
{absK,lamda,Vray,Aray,absdIdp3,p3}=viscoKVdEdp3[\[Alpha],\[Beta],f,fc,AcoustPars,opt];
},{
(* 3-ip *)
f=(ip-np-1)*df;
{absK,lamda,Vray,Aray,absdIdp3,p3}=viscoKVdEdp3[\[Alpha],\[Beta],f,fc,AcoustPars,opt];
}
];
omega=2Pi f;
fdata[[ip]]=1/(2Pi*R*Sqrt[absK]) fdata[[ip]]/absdIdp3 Exp[-omega Aray] Exp[I (omega R /Vray+lamda)];

}
];

(* --------- zero frequency --------- *)
ip=1;
fdata[[ip]]=0.;

tdata0=InverseFourier[fdata,FourierParameters->{-1,1}];
tdata=Re[tdata0];

Return[tdata];
];
