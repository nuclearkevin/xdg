#include "xdg/prism_contain.h"

#include "xdg/constants.h"
#include "xdg/error.h"
#include "xdg/geometry/face_common.h"
#include "xdg/geometry/plucker.h"
#include "xdg/ray_tracing_interface.h"
#include "xdg/ray.h"
#include "xdg/vec3da.h"

namespace xdg
{

bool prism_containment_test(const Position& point,
                            const std::array<Vertex, 6>& verts)
{
  using VertArray = std::array<Vertex, 6>;
  static constexpr std::array<std::array<int, 3>, 2> k_prism_tri_faces = {
    std::array<int, 3>{0, 2, 1},
    std::array<int, 3>{3, 4, 5}
  };

  static constexpr std::array<std::array<int, 4>, 3> k_prism_quad_faces = {
    std::array<int, 4>{0, 1, 4, 3},
    std::array<int, 4>{1, 2, 5, 4},
    std::array<int, 4>{2, 0, 3, 5},
  };

  Position centroid {0.0, 0.0, 0.0};
  for (const auto& v : verts) {
    centroid += v;
  }
  centroid = centroid / 6.0;

  Direction direction = centroid - point;
  const double distance_to_centroid = direction.length();
  if (distance_to_centroid <= TINY_BIT) {
    return true;
  }
  direction /= distance_to_centroid;

  // Check the triangle faces.
  for (const auto& face : k_prism_tri_faces) {
    if (crosses_boundary_triangle<VertArray>(point, centroid, direction, verts, face[0], face[1], face[2])) {
      return false;
    }
  }

  // Check the quad faces.
  for (const auto& face : k_prism_quad_faces) {
    // diagonal is selected based on a lexicographic comparison
    // of vertex coordinates so diagonals are chosen consistently
    // regardless of face connectivity ordering
    if (select_diagonal<VertArray>(verts, face)) {
      if (crosses_boundary_triangle<VertArray>(point, centroid, direction, verts, face[0], face[1], face[2]) ||
          crosses_boundary_triangle<VertArray>(point, centroid, direction, verts, face[0], face[2], face[3])) {
        return false;
      }
    } else {
      if (crosses_boundary_triangle<VertArray>(point, centroid, direction, verts, face[1], face[2], face[3]) ||
          crosses_boundary_triangle<VertArray>(point, centroid, direction, verts, face[1], face[3], face[0])) {
        return false;
      }
    }
  }

  return true;
}

// Embree call back functions for element search
void PrismIntersectionFunc(RTCIntersectFunctionNArguments* args)
{
  const VolumeElementsUserData* user_data = (const VolumeElementsUserData*)args->geometryUserPtr;
  const MeshManager* mesh_manager = user_data->mesh_manager;

  const PrimitiveRef primitive_ref = user_data->prim_ref_buffer[args->primID];
  auto vertices = mesh_manager->element_vertices(primitive_ref.primitive_id);
  if (vertices.size() != 6) {
    fatal_error("PrismIntersectionFunc expected 6 vertices, got {}", vertices.size());
  }

  std::array<Vertex, 6> verts;
  for (size_t i = 0; i < verts.size(); ++i) {
    verts[i] = vertices[i];
  }

  RTCDualRayHit* rayhit = (RTCDualRayHit*)args->rayhit;
  RTCSurfaceDualRay& ray = rayhit->ray;

  Position ray_origin = {ray.dorg[0], ray.dorg[1], ray.dorg[2]};

  bool inside = prism_containment_test(ray_origin, verts);
  if (!inside) return;

  rayhit->hit.u = 0.0;
  rayhit->hit.v = 0.0;
  rayhit->hit.Ng_x = 0.0;
  rayhit->hit.Ng_y = 0.0;
  rayhit->hit.Ng_z = 0.0;
  rayhit->hit.geomID = args->geomID;
  rayhit->hit.primID = args->primID;
}

void PrismOcclusionFunc(RTCOccludedFunctionNArguments* args)
{
  const VolumeElementsUserData* user_data = (const VolumeElementsUserData*)args->geometryUserPtr;
  const MeshManager* mesh_manager = user_data->mesh_manager;

  const PrimitiveRef primitive_ref = user_data->prim_ref_buffer[args->primID];
  auto vertices = mesh_manager->element_vertices(primitive_ref.primitive_id);
  if (vertices.size() != 6) {
    fatal_error("PrismOcclusionFunc expected 6 vertices, got {}", vertices.size());
  }

  std::array<Vertex, 6> verts;
  for (size_t i = 0; i < verts.size(); ++i) {
    verts[i] = vertices[i];
  }

  RTCElementDualRay* ray = (RTCElementDualRay*)args->ray;
  Position ray_origin = {ray->dorg[0], ray->dorg[1], ray->dorg[2]};

  bool inside = prism_containment_test(ray_origin, verts);
  if (!inside) return;

  // indicate that the query is complete due to a successful containment result
  ray->element = primitive_ref.primitive_id;
  ray->set_tfar(-INFTY);
}

} // namespace xdg
