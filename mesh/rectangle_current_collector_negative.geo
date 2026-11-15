// -----------------------------------------------------------------------------
// Negative current collector: structured quadrilateral mesh with top protrusion (tab)
// Generates rectangle_current_collector_negative.msh:
//   gmsh rectangle_current_collector_negative.geo -2 -o rectangle_current_collector_negative.msh
//
// The domain (surfaces 1-6) is split exactly as in rectangle_current_collector_positive.geo, so the
// domain elements are numbered in the same order in both meshes, as batree requires (the coupling
// pairs element e of the current space with element e of each collector). The tab sits at
// x = 0.07..0.09 m; its base points split the top edge of surface 6.
// -----------------------------------------------------------------------------
SetFactory("OpenCASCADE");
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;
Mesh.MshFileVersion = 2.2;

// -----------------------------------------------------------------------------
// Mesh size
// -----------------------------------------------------------------------------
h = 0.01;

// -----------------------------------------------------------------------------
// Geometry
// -----------------------------------------------------------------------------
Lx = 0.1;      // Width of main domain
Ly = 0.12;     // Height of main domain
w  = 0.02;     // Tab width
hp = 0.02;     // Tab height
xL = 0.01;     // Domain split as in the positive mesh (left edge of the positive tab)
xR = xL + w;   // Domain split as in the positive mesh (right edge of the positive tab)
xT = 0.07;     // Left edge of the negative tab
yM = Ly/2;     // Horizontal split

// -----------------------------------------------------------------------------
// Points
// -----------------------------------------------------------------------------
// Bottom
Point(1) = {0, 0, 0, h};
Point(2) = {xL, 0, 0, h};
Point(3) = {xR, 0, 0, h};
Point(4) = {Lx, 0, 0, h};

// Middle
Point(5) = {0, yM, 0, h};
Point(6) = {xL, yM, 0, h};
Point(7) = {xR, yM, 0, h};
Point(8) = {Lx, yM, 0, h};

// Top of main rectangle
Point(9)  = {0, Ly, 0, h};
Point(10) = {xL, Ly, 0, h};
Point(11) = {xR, Ly, 0, h};
Point(12) = {Lx, Ly, 0, h};

// Top of the tab
Point(13) = {xT, Ly + hp, 0, h};
Point(14) = {xT + w, Ly + hp, 0, h};

// Base of the tab, on the top edge of the main rectangle
Point(15) = {xT, Ly, 0, h};
Point(16) = {xT + w, Ly, 0, h};

// -----------------------------------------------------------------------------
// Horizontal lines
// -----------------------------------------------------------------------------
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {7,8};

Line(7) = {9,10};
Line(8) = {10,11};
Line(9) = {11,15};   // top edge, from x = xR to the tab
Line(21) = {15,16};  // top edge under the tab
Line(22) = {16,12};  // top edge, from the tab to x = Lx

Line(10) = {14,13};  // tab top, right to left (the tab outline below runs counter-clockwise)

// -----------------------------------------------------------------------------
// Vertical lines
// -----------------------------------------------------------------------------
Line(11) = {1,5};
Line(12) = {5,9};

Line(13) = {2,6};
Line(14) = {6,10};

Line(16) = {3,7};
Line(17) = {7,11};

Line(19) = {4,8};
Line(20) = {8,12};

Line(15) = {13,15};  // tab left side
Line(18) = {16,14};  // tab right side

// -----------------------------------------------------------------------------
// Surfaces
// -----------------------------------------------------------------------------
Curve Loop(1) = {1,13,-4,-11};
Plane Surface(1) = {1};

Curve Loop(2) = {2,16,-5,-13};
Plane Surface(2) = {2};

Curve Loop(3) = {3,19,-6,-16};
Plane Surface(3) = {3};

Curve Loop(4) = {4,14,-7,-12};
Plane Surface(4) = {4};

Curve Loop(5) = {5,17,-8,-14};
Plane Surface(5) = {5};

Curve Loop(6) = {6,20,-22,-21,-9,-17};
Plane Surface(6) = {6};

// Surface 7 (tab), counter-clockwise, so the tab elements are too
Curve Loop(7) = {21,18,10,15};
Plane Surface(7) = {7};

// -----------------------------------------------------------------------------
// Structured mesh divisions
// -----------------------------------------------------------------------------
nxL = 1;
nxC = 2;
nxR = 7;

nxA = 4;   // top edge, x = xR .. xT
nxT = 2;   // under the tab
nxB = 1;   // top edge, x = xT + w .. Lx  (nxA + nxT + nxB = nxR)

nyB = 6;
nyT = 6;
nyP = 2;

// Horizontal divisions
Transfinite Curve{1,4,7} = nxL + 1;
Transfinite Curve{2,5,8} = nxC + 1;
Transfinite Curve{3,6} = nxR + 1;
Transfinite Curve{9} = nxA + 1;
Transfinite Curve{21,10} = nxT + 1;
Transfinite Curve{22} = nxB + 1;

// Vertical divisions
Transfinite Curve{11,13,16,19} = nyB + 1;
Transfinite Curve{12,14,17,20} = nyT + 1;
Transfinite Curve{15,18} = nyP + 1;

// -----------------------------------------------------------------------------
// Transfinite surfaces
// -----------------------------------------------------------------------------
Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Surface{4};
Transfinite Surface{5};
Transfinite Surface{6} = {7, 8, 12, 11};   // corners: its top edge is in three pieces
Transfinite Surface{7};

Recombine Surface{1:7};

// -----------------------------------------------------------------------------
// Physical groups
// -----------------------------------------------------------------------------
Physical Surface("domain") = {1:6};
Physical Surface("tab") = {7};
Physical Curve("protrusion_top") = {10};

Mesh.MshFileVersion = 2.2;
