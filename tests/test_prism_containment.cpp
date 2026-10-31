// stl includes
#include <array>

// testing includes
#include <catch2/catch_test_macros.hpp>

#include "xdg/prism_contain.h"
#include "xdg/vec3da.h"

using namespace xdg;

TEST_CASE("Prism/Wedge Containment Unit Test")
{
  // Define prism/wedge vertices in MBCN Prism6 ordering.
  const std::array<Vertex, 6> verts = {
    Vertex(0.0, 0.0, 0.0),
    Vertex(1.0, 1.0, 0.0),
    Vertex(0.0, 1.0, 0.0),
    Vertex(0.0, 0.0, 1.0),
    Vertex(1.0, 1.0, 1.0),
    Vertex(0.0, 1.0, 1.0)
  };

  // Test points
  Position inside_point(0.4, 0.6, 0.5); // Inside the prism
  Position outside_point(2.0, 2.0, 2.0); // Clearly outside
  Position boundary_point(0.5, 0.5, 0.5); // On the diagonal face

  // Check containment
  CHECK(prism_containment_test(inside_point, verts) == true);
  CHECK(prism_containment_test(outside_point, verts) == false);
  CHECK(prism_containment_test(boundary_point, verts) == true);

  // Check point that is co-planar with one of the faces, but outside the hex.
  Position coplanar_exterior_point(1.5, 1.0, 0.5);
  CHECK(prism_containment_test(coplanar_exterior_point, verts) == false);

  // Test points near the faces but not coplanar.
  Position near_bottom_face(0.25, 0.25, 0.1);
  Position near_top_face(0.25, 0.25, 0.9);
  Position near_left_face(0.05, 0.25, 0.5);
  Position near_back_face(0.25, 0.95, 0.5);
  Position near_diag_face(0.5, 0.55, 0.5);

  CHECK(prism_containment_test(near_bottom_face, verts) == true);
  CHECK(prism_containment_test(near_top_face, verts) == true);
  CHECK(prism_containment_test(near_left_face, verts) == true);
  CHECK(prism_containment_test(near_diag_face, verts) == true);
  CHECK(prism_containment_test(near_back_face, verts) == true);

  // Test points just outside the prism.
  Position just_outside_bottom_face(0.25, 0.25, -0.1);
  Position just_outside_top_face(0.25, 0.25, 1.1);
  Position just_outside_left_face(-0.05, 0.25, 0.5);
  Position just_outside_back_face(0.25, 1.05, 0.5);
  Position just_outside_diag_face(0.5, 0.45, 0.5);

  CHECK(prism_containment_test(just_outside_bottom_face, verts) == false);
  CHECK(prism_containment_test(just_outside_top_face, verts) == false);
  CHECK(prism_containment_test(just_outside_left_face, verts) == false);
  CHECK(prism_containment_test(just_outside_diag_face, verts) == false);
  CHECK(prism_containment_test(just_outside_back_face, verts) == false);

  // Check lower-dimensional boundaries.
  Position edge_point(0.5, 0.5, 0.0); // on a diagonal edge
  CHECK(prism_containment_test(edge_point, verts) == true);
  Position vertex_point(1.0, 1.0, 1.0);
  CHECK(prism_containment_test(vertex_point, verts) == true);

  // Create a second prism with different positions and dimensions as a sanity check.
  const std::array<Vertex, 6> verts_2 = {
    Vertex(1.0, 1.0, 1.0 ),
    Vertex(3.0, 1.0, 1.0 ),
    Vertex(3.0, 3.0, 1.0 ),
    Vertex(1.0, 1.0, 10.0),
    Vertex(3.0, 1.0, 10.0),
    Vertex(3.0, 3.0, 10.0)
  };

  // Check containment for the second prism.
  Position inside_point_2(2.0, 1.5, 4.5);
  Position outside_point_2(5.0, 5.0, 11.0);
  Position boundary_point_2(2.0, 1.0, 4.5);

  CHECK(prism_containment_test(inside_point_2, verts_2) == true);
  CHECK(prism_containment_test(outside_point_2, verts_2) == false);
  CHECK(prism_containment_test(boundary_point_2, verts_2) == true);

  // Test points near faces of second prism.
  Position near_face_2(2.9, 1.5, 4.5);
  Position just_outside_face_2(3.1, 1.5, 4.5);

  CHECK(prism_containment_test(near_face_2, verts_2) == true);
  CHECK(prism_containment_test(just_outside_face_2, verts_2) == false);

  // Test points that are in one prism but not the other.
  CHECK(prism_containment_test(inside_point, verts_2) == false);
  CHECK(prism_containment_test(inside_point_2, verts) == false);
}
