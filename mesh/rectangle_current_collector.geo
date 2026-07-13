// -----------------------------------------------------------------------------
// Structured 7-patch quadrilateral mesh with top protrusion
// Main domain: 100 x 120
// Protrusion: 10 x 10
// Protrusion starts 10 units from the left edge
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
Ly = 0.12;      // Height of main domain

w  = 0.02;       // Protrusion width
hp = 0.02;       // Protrusion height

xL = 0.02;       // Left edge of protrusion
xR = xL + w;   // Right edge of protrusion

yM = Ly/2;     // Horizontal split (60)

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

// Top of protrusion
Point(13) = {xL, Ly + hp, 0, h};
Point(14) = {xR, Ly + hp, 0, h};

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
Line(9) = {11,12};

Line(10) = {13,14};

// -----------------------------------------------------------------------------
// Vertical lines
// -----------------------------------------------------------------------------

Line(11) = {1,5};
Line(12) = {5,9};

Line(13) = {2,6};
Line(14) = {6,10};
Line(15) = {10,13};

Line(16) = {3,7};
Line(17) = {7,11};
Line(18) = {11,14};

Line(19) = {4,8};
Line(20) = {8,12};

// -----------------------------------------------------------------------------
// Surfaces
// -----------------------------------------------------------------------------

// Surface 1
Curve Loop(1) = {1,13,-4,-11};
Plane Surface(1) = {1};

// Surface 2
Curve Loop(2) = {2,16,-5,-13};
Plane Surface(2) = {2};

// Surface 3
Curve Loop(3) = {3,19,-6,-16};
Plane Surface(3) = {3};

// Surface 4
Curve Loop(4) = {4,14,-7,-12};
Plane Surface(4) = {4};

// Surface 5
Curve Loop(5) = {5,17,-8,-14};
Plane Surface(5) = {5};

// Surface 6
Curve Loop(6) = {6,20,-9,-17};
Plane Surface(6) = {6};

// Surface 7 (protrusion)
Curve Loop(7) = {10,-18,-8,15};
Plane Surface(7) = {7};

// -----------------------------------------------------------------------------
// Structured mesh divisions
// -----------------------------------------------------------------------------

nxL = 2;
nxC = 2;
nxR = 6;

nyB = 6;
nyT = 6;
nyP = 2;

// Horizontal divisions
Transfinite Curve{1,4,7} = nxL + 1;
Transfinite Curve{2,5,8,10} = nxC + 1;
Transfinite Curve{3,6,9} = nxR + 1;

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
Transfinite Surface{6};
Transfinite Surface{7};

Recombine Surface{1:7};

// -----------------------------------------------------------------------------
// Physical groups
// -----------------------------------------------------------------------------

Physical Surface("domain") = {1:6};
Physical Surface("tab") = {7};
Physical Curve("protrusion_top") = {10}

Mesh.MshFileVersion = 2.2;