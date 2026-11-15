// Scaled positive current collector mesh: rectangle_current_collector_positive.geo with every
// coordinate divided by L0 = sqrt(A_area) = sqrt(0.012 m^2), batree's collector length scale.
// Generates rectangle_current_collector_positive_scaled.msh:
//   gmsh rectangle_current_collector_positive_scaled.geo -2 -o rectangle_current_collector_positive_scaled.msh
Include "rectangle_current_collector_positive.geo";
Mesh.ScalingFactor = 1 / Sqrt(0.012);
