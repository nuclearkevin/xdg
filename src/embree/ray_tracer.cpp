#include "xdg/embree/ray_tracer.h"
#include "xdg/error.h"
#include "xdg/geometry_data.h"
#include "xdg/ray.h"

namespace xdg {

void error(void* dum, RTCError code, const char* str) {
  if (code != RTC_ERROR_NONE)
    fatal_error("Embree error: {}", str);
}

EmbreeRayTracer::EmbreeRayTracer()
{
  device_ = rtcNewDevice(nullptr);
  rtcSetDeviceErrorFunction(device_, (RTCErrorFunction)error, nullptr);
}

EmbreeRayTracer::~EmbreeRayTracer()
{
  rtcReleaseDevice(device_);
}

void EmbreeRayTracer::init()
{

}

RTCScene EmbreeRayTracer::create_embree_scene() {
  RTCScene rtcscene = rtcNewScene(device_);
  rtcSetSceneFlags(rtcscene, RTC_SCENE_FLAG_ROBUST);
  rtcSetSceneBuildQuality(rtcscene, RTC_BUILD_QUALITY_HIGH);
  return rtcscene;
}

TreeID
EmbreeRayTracer::register_volume(const std::shared_ptr<MeshManager> mesh_manager,
                           MeshID volume_id)
{
  TreeID tree = next_tree_id();
  trees_.push_back(tree);
  auto volume_scene = this->create_embree_scene();

  // allocate storage for this volume
  auto volume_faces = mesh_manager->get_volume_faces(volume_id);
  this->primitive_ref_storage_[volume_scene].resize(volume_faces.size());
  auto& triangle_storage = this->primitive_ref_storage_[volume_scene];

  auto volume_surfaces = mesh_manager->get_volume_surfaces(volume_id);
  int storage_offset {0};
  for (auto& surface_id : volume_surfaces) {
    // get the sense of this surface with respect to the volume
    Sense triangle_sense {Sense::UNSET};
    auto surf_to_vol_senses = mesh_manager->get_parent_volumes(surface_id);
    if (volume_id == surf_to_vol_senses.first) triangle_sense = Sense::FORWARD;
    else if (volume_id == surf_to_vol_senses.second) triangle_sense = Sense::REVERSE;
    else fatal_error("Volume {} is not a parent of surface {}", volume_id, surface_id);

    auto surface_faces = mesh_manager->get_surface_faces(surface_id);
    for (int i = 0; i < surface_faces.size(); ++i) {
      auto& primitive_ref = triangle_storage[i + storage_offset];
      primitive_ref.primitive_id = surface_faces[i];
      primitive_ref.sense = triangle_sense;
    }
    storage_offset += surface_faces.size();
 }

  PrimitiveRef* tri_ref_ptr = triangle_storage.data();

  auto bump = bounding_box_bump(mesh_manager, volume_id);
/*
  TODO: none of the above is ray tracer specific. This can be a part of the common register_volume
  implementation. Then another virtual function can be called register_volume_RT_specific.
  in this scenario register_volume isn't a virtual function but instead calls a virtual function.
  That virtual function will be overrided to do the RT specific things in registering the volume.
  Primitive_ref_storage is the only non-local variable used in this context that is a member of the derived class.
  However when implementing GPRT the above may need to be written differently with GPU buffers in mind.
  So we will leave it as is for now.
*/

  // create a new geometry for each surface
  int buffer_start = 0;
  for (auto surface : volume_surfaces) {
    auto surface_triangles = mesh_manager->get_surface_faces(surface);
    RTCGeometry surface_geometry = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_USER);
    rtcSetGeometryUserPrimitiveCount(surface_geometry, surface_triangles.size());
    unsigned int embree_surface = rtcAttachGeometry(volume_scene, surface_geometry);
    this->surface_to_geometry_map_[surface] = surface_geometry;

    std::shared_ptr<GeometryUserData> surface_data = std::make_shared<GeometryUserData>();
    surface_data->box_bump = bump;
    surface_data->surface_id = surface;
    surface_data->mesh_manager = mesh_manager.get();
    surface_data->prim_ref_buffer = tri_ref_ptr + buffer_start;
    user_data_map_[surface_geometry] = surface_data;

    // TODO: This could be a problem if user_data_map_ is reallocated?
    rtcSetGeometryUserData(surface_geometry, user_data_map_[surface_geometry].get());

    for (int i = 0; i < surface_triangles.size(); ++i) {
      auto& triangle_ref = surface_data->prim_ref_buffer[i];
    }
    buffer_start += surface_triangles.size();

    // set the bounds function
    rtcSetGeometryBoundsFunction(surface_geometry, (RTCBoundsFunction)&TriangleBoundsFunc, nullptr);
    rtcSetGeometryIntersectFunction(surface_geometry, (RTCIntersectFunctionN)&TriangleIntersectionFunc);
    rtcSetGeometryOccludedFunction(surface_geometry, (RTCOccludedFunctionN)&TriangleOcclusionFunc);

    rtcCommitGeometry(surface_geometry);
  }
  rtcCommitScene(volume_scene);

  tree_to_scene_map_[tree] = volume_scene;
  return tree;
}

bool EmbreeRayTracer::point_in_volume(TreeID tree,
                                const Position& point,
                                const Direction* direction,
                                const std::vector<MeshID>* exclude_primitives) const
{
  RTCScene scene = tree_to_scene_map_.at(tree);
  RTCDRayHit rayhit; // embree specfic rayhit struct (payload?)
  rayhit.ray.set_org(point);
  if (direction != nullptr) rayhit.ray.set_dir(*direction);
  else rayhit.ray.set_dir({1. / std::sqrt(2.0), 1 / std::sqrt(2.0), 0.0});
  rayhit.ray.rf_type = RayFireType::VOLUME;
  rayhit.ray.orientation = HitOrientation::ANY;
  rayhit.ray.set_tfar(INFTY);
  rayhit.ray.set_tnear(0.0);

  if (exclude_primitives != nullptr) rayhit.ray.exclude_primitives = exclude_primitives;

  {
    rtcIntersect1(scene, (RTCRayHit*)&rayhit);
  }

  // if the ray hit nothing, the point is outside of the volume
  if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;

  // use the hit triangle normal to determine if the intersection is
  // exiting or entering
  return rayhit.ray.ddir.dot(rayhit.hit.dNg) > 0.0;
}

std::pair<double, MeshID>
EmbreeRayTracer::ray_fire(TreeID tree,
                    const Position& origin,
                    const Direction& direction,
                    const double dist_limit,
                    HitOrientation orientation,
                    std::vector<MeshID>* const exclude_primitves)
{
  RTCScene scene = tree_to_scene_map_.at(tree);
  RTCDRayHit rayhit;
  // set ray data
  rayhit.ray.set_org(origin);
  rayhit.ray.set_dir(direction);
  rayhit.ray.set_tfar(dist_limit);
  rayhit.ray.set_tnear(0.0);
  rayhit.ray.rf_type = RayFireType::VOLUME;
  rayhit.ray.orientation = orientation;
  rayhit.ray.mask = -1; // no mask
  if (exclude_primitves != nullptr) rayhit.ray.exclude_primitives = exclude_primitves;

  // fire the ray
  {
    rtcIntersect1(scene, (RTCRayHit*)&rayhit);
    // TODO: I don't quite understand this...
    rayhit.hit.Ng_x *= -1.0;
    rayhit.hit.Ng_y *= -1.0;
    rayhit.hit.Ng_z *= -1.0;
  }

  if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID)
    return {INFTY, ID_NONE};
  else
    if (exclude_primitves) exclude_primitves->push_back(rayhit.hit.primitive_ref->primitive_id);
    return {rayhit.ray.dtfar, rayhit.hit.surface};
}

void EmbreeRayTracer::closest(TreeID tree,
                        const Position& point,
                        double& distance,
                        MeshID& triangle)
{
  RTCScene scene = tree_to_scene_map_.at(tree);
  RTCDPointQuery query;
  query.set_point(point);

  RTCPointQueryContext context;
  rtcInitPointQueryContext(&context);

  rtcPointQuery(scene, &query, &context, (RTCPointQueryFunction)&TriangleClosestFunc, &scene);

  if (query.geomID == RTC_INVALID_GEOMETRY_ID) {
    distance = INFTY;
    return;
  }

  distance = query.dradius;
  triangle = query.primitive_ref->primitive_id;
}

void EmbreeRayTracer::closest(TreeID scene,
                        const Position& point,
                        double& distance)
{
  MeshID triangle;
  closest(scene, point, distance, triangle);
}

bool EmbreeRayTracer::occluded(TreeID tree,
                         const Position& origin,
                         const Direction& direction,
                         double& distance) const
{
  RTCScene scene = tree_to_scene_map_.at(tree);
  RTCDRay ray;
  ray.set_org(origin);
  ray.set_dir(direction);
  ray.set_tfar(INFTY);
  ray.set_tnear(0.0);
  ray.rf_type = RayFireType::FIND_VOLUME;
  ray.orientation = HitOrientation::ANY;
  ray.flags = 0;
  ray.mask = -1; // no mask

  // fire the ray
  {
    rtcOccluded1(scene, (RTCRay*)&ray);
  }

  distance = ray.dtfar;
  return distance != INFTY;
}

} // namespace xdg