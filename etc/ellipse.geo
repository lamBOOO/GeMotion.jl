SetFactory("OpenCASCADE");

// middle point
Point(9) = {0.5, 0.5, 0.5};

// major axis
Point(10) = {0.7, 0.5, 0.5};
Point(11) = {0.5, 0.5, 0.7};
Point(12) = {0.3, 0.5, 0.5};
Point(13) = {0.5, 0.5, 0.3};
Point(14) = {0.5, 1.5, 0.5};
Point(15) = {0.5, -0.5, 0.5};

// ellipsoid, lines
Ellipse(13) = {10, 9, 11, 11};
Ellipse(14) = {12, 9, 11, 11};
Ellipse(15) = {12, 9, 13, 13};
Ellipse(16) = {13, 9, 10, 10};

Ellipse(17) = {14, 9, 10, 10};
Ellipse(18) = {14, 9, 11, 11};
Ellipse(19) = {14, 9, 12, 12};
Ellipse(20) = {14, 9, 13, 13};

Ellipse(21) = {15, 9, 10, 10};
Ellipse(22) = {15, 9, 11, 11};
Ellipse(23) = {15, 9, 12, 12};
Ellipse(24) = {15, 9, 13, 13};


// line loops
Line Loop(1) = {16,17,-20};
Surface(1) = {1};
