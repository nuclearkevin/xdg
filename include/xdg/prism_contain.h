#ifndef XDG_PRISM_CONTAIN_H
#define XDG_PRISM_CONTAIN_H

#include <array>

#include "xdg/vec3da.h"

// Forward declarations for Embree callback argument types
struct RTCIntersectFunctionNArguments;
struct RTCOccludedFunctionNArguments;

namespace xdg {

/**
 * @brief Determines if a point is inside or on the boundary of a prism / wedge
 * using MBCN canonical Wedge6 ordering (Tautges 2010).
 *
 * The prism is defined by six vertices ordered per the MBCN
 * canonical numbering. The containment test casts rays toward the prism
 * centroid and uses Plucker edge-coordinate sign tests against triangulated
 * faces; points on the boundary are considered inside.
 */
bool prism_containment_test(const Position& point,
                            const std::array<Vertex, 6>& verts);

// Embree call back functions for element search
void PrismIntersectionFunc(RTCIntersectFunctionNArguments* args);
void PrismOcclusionFunc(RTCOccludedFunctionNArguments* args);

} // namespace xdg

#endif // include guard
