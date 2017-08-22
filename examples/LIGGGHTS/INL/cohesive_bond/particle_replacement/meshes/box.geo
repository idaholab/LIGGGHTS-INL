cl_1 = 0.0025;
l = 0.05;  // half length
w = 0.005;  // half width
h = 0.005;  // half height

// front
Point(1) = {-w, -l, -h, cl_1}; // counter cockwise
Point(2) = { w, -l, -h, cl_1};
Point(3) = { w, -l,  h, cl_1};
Point(4) = {-w, -l,  h, cl_1};
Point(5) = {-w,  l, -h, cl_1}; // clockwise
Point(6) = {-w,  l,  h, cl_1};
Point(7) = { w,  l,  h, cl_1};
Point(8) = { w,  l, -h, cl_1};
// front
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
// left and right
Line(9) = {5, 1};
Line(10) = {4, 6};
Line(11) = {2, 8};
Line(12) = {7, 3};
// planes
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};
Line Loop(17) = {9, -4, 10, -5};
Plane Surface(18) = {17};
Line Loop(19) = {11, -7, 12, -2};
Plane Surface(20) = {19};
Line Loop(21) = {1, 11, 8, 9};
Plane Surface(22) = {21};
Line Loop(23) = {3, 10, 6, 12};
Plane Surface(24) = {23};
