// ------- user-defined refinement factor -----------------
If(!Exists(ref))
  p = 0 ;            // ← 1 = original size, 0.5 = half, 0.25 = quarter ...
EndIf
// --------------------------------------------------------

ref = 2^(-p);

// reference sizes (keep your old numbers here)
lcCore_ref  = 0.10 ;
hFirst_ref  = 0.004 ;

// actual sizes used by the mesh
lcCore  = lcCore_ref  * ref ;
hFirst  = hFirst_ref  * ref ;

growth  = 1.25 ;
nLayers = 8 ;
Thickness = hFirst * (1 - growth^nLayers) / (1 - growth);

Mesh.CharacteristicLengthMin = hFirst;
Mesh.CharacteristicLengthMax = lcCore;

// geometry --------------------------------------------------------
R1 = 2./3.;  R2 = 5./3.;
Point(1)    = {0,0,0,lcCore};
Point(1000) = { R1, 0, 0, lcCore};
Point(1001) = {-R1, 0, 0, lcCore};
Point(1100) = { R2, 0, 0, lcCore};
Point(1101) = {-R2, 0, 0, lcCore};

Circle(2000) = {1000,1,1001};   // inner
Circle(2001) = {1001,1,1000};
Circle(2100) = {1100,1,1101};   // outer
Circle(2101) = {1101,1,1100};

Line Loop(30) = {2100,2101};    // outer loop  (positive)
Line Loop(31) = {2001,2000};    // inner loop  (negative → hole)
Plane Surface(40) = {30,31};    // << ONE surface with a hole

// boundary-layer field -------------------------------------------
Field[1] = BoundaryLayer;
Field[1].EdgesList = {2000,2001};  // only the circle
Field[1].hwall_n   = hFirst;
Field[1].ratio     = growth;
Field[1].thickness = Thickness;
Field[1].hfar      = lcCore;
// Field[1].Quads     = 1;

BoundaryLayer Field = 1;
