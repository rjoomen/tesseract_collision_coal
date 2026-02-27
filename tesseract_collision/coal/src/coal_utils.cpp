/**
 * @file coal_utils.cpp
 * @brief Tesseract Coal Utility Functions.
 *
 * @author Roelof Oomen, Levi Armstrong
 * @date Dec 18, 2017
 *
 * @copyright Copyright (c) 2017, Southwest Research Institute
 *
 * @par License
 * Software License Agreement (BSD)
 * @par
 * All rights reserved.
 * @par
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * @par
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 * @par
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <coal/collision_data.h>
#include <coal/collision.h>
#include <coal/distance.h>
#include <coal/BVH/BVH_model.h>
#include <coal/shape/geometric_shapes.h>
#include <coal/shape/geometric_shapes_utility.h>
#include <coal/narrowphase/support_functions.h>
#include <coal/narrowphase/minkowski_difference.h>
#include <coal/narrowphase/gjk.h>
#include <coal/shape/convex.h>
#include <coal/data_types.h>
#include <coal/octree.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_collision/coal/coal_utils.h>
#include <tesseract_collision/coal/coal_collision_geometry_cache.h>
#include <tesseract_collision/coal/coal_casthullshape.h>
#include <tesseract_collision/coal/coal_casthullshape_utility.h>
#include <tesseract_geometry/geometries.h>

namespace tesseract_collision::tesseract_collision_coal
{
namespace
{
CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Plane::ConstPtr& geom)
{
  return std::make_shared<coal::Plane>(geom->getA(), geom->getB(), geom->getC(), geom->getD());
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Box::ConstPtr& geom)
{
  return std::make_shared<coal::Box>(geom->getX(), geom->getY(), geom->getZ());
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Sphere::ConstPtr& geom)
{
  return std::make_shared<coal::Sphere>(geom->getRadius());
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Cylinder::ConstPtr& geom)
{
  return std::make_shared<coal::Cylinder>(geom->getRadius(), geom->getLength());
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Cone::ConstPtr& geom)
{
  return std::make_shared<coal::Cone>(geom->getRadius(), geom->getLength());
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Capsule::ConstPtr& geom)
{
  return std::make_shared<coal::Capsule>(geom->getRadius(), geom->getLength());
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Mesh::ConstPtr& geom)
{
  const int vertex_count = geom->getVertexCount();
  const int triangle_count = geom->getFaceCount();
  const tesseract_common::VectorVector3d& vertices = *(geom->getVertices());
  const Eigen::VectorXi& triangles = *(geom->getFaces());

  auto g = std::make_shared<coal::BVHModel<coal::OBBRSS>>();
  if (vertex_count > 0 && triangle_count > 0)
  {
    std::vector<coal::Triangle> tri_indices(static_cast<size_t>(triangle_count));
    for (int i = 0; i < triangle_count; ++i)
    {
      assert(triangles[4L * i] == 3);
      tri_indices[static_cast<size_t>(i)] = coal::Triangle(static_cast<size_t>(triangles[(4 * i) + 1]),
                                                           static_cast<size_t>(triangles[(4 * i) + 2]),
                                                           static_cast<size_t>(triangles[(4 * i) + 3]));
    }

    g->beginModel();
    g->addSubModel(vertices, tri_indices);
    g->endModel();

    return g;
  }

  CONSOLE_BRIDGE_logError("The mesh is empty!");
  return nullptr;
}

// Coal polygon type (modelled after TriangleTpl)
template <typename IndexType_>
struct PolygonTpl : Eigen::Matrix<IndexType_, -1, 1>
{
  using IndexType = IndexType_;
  using size_type = int;

  // template <typename OtherIndexType>
  // friend class Polygon;

  /// @brief Default constructor
  PolygonTpl() = default;

  /// @brief Copy constructor
  PolygonTpl(const PolygonTpl& other) : Eigen::Matrix<IndexType_, -1, 1>(other) {}

  /// @brief Copy constructor from another vertex index type.
  template <typename OtherIndexType>
  PolygonTpl(const PolygonTpl<OtherIndexType>& other)
  {
    *this = other;
  }

  /// @brief Copy operator
  PolygonTpl& operator=(const PolygonTpl& other)
  {
    this->_set(other);
    return *this;
  }

  /// @brief Copy operator from another index type.
  template <typename OtherIndexType>
  PolygonTpl& operator=(const PolygonTpl<OtherIndexType>& other)
  {
    *this = other.template cast<OtherIndexType>();
    return *this;
  }

  template <typename OtherIndexType>
  PolygonTpl<OtherIndexType> cast() const
  {
    PolygonTpl<OtherIndexType> res;
    res._set(*this);
    return res;
  }
};

using Polygon = PolygonTpl<std::uint32_t>;

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::ConvexMesh::ConstPtr& geom)
{
  const auto vertex_count = geom->getVertexCount();
  const auto face_count = geom->getFaceCount();
  const auto& faces = *geom->getFaces();

  if (vertex_count > 0 && face_count > 0)
  {
    auto vertices = std::const_pointer_cast<tesseract_common::VectorVector3d>(geom->getVertices());

    auto new_faces = std::make_shared<std::vector<Polygon>>();
    new_faces->reserve(face_count);
    for (int i = 0; i < faces.size(); ++i)
    {
      Polygon new_face;
      // First value of each face is the number of vertices
      new_face.resize(faces[i]);
      for (std::uint32_t& j : new_face)
      {
        ++i;
        j = faces[i];
      }
      new_faces->emplace_back(new_face);
    }
    assert(new_faces->size() == face_count);

    return std::make_shared<coal::Convex<Polygon>>(vertices, vertex_count, new_faces, face_count);
  }

  CONSOLE_BRIDGE_logError("The mesh is empty!");
  return nullptr;
}

CollisionGeometryPtr createShapePrimitive(const tesseract_geometry::Octree::ConstPtr& geom)
{
  switch (geom->getSubType())
  {
    case tesseract_geometry::OctreeSubType::BOX:
    {
      return std::make_shared<coal::OcTree>(geom->getOctree());
    }
    default:
    {
      CONSOLE_BRIDGE_logError("This Coal octree sub shape type (%d) is not supported for geometry octree",
                              static_cast<int>(geom->getSubType()));
      return nullptr;
    }
  }
}
}  // namespace

CollisionGeometryPtr createShapePrimitiveHelper(const CollisionShapeConstPtr& geom)
{
  switch (geom->getType())
  {
    case tesseract_geometry::GeometryType::PLANE:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Plane>(geom));
    }
    case tesseract_geometry::GeometryType::BOX:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Box>(geom));
    }
    case tesseract_geometry::GeometryType::SPHERE:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Sphere>(geom));
    }
    case tesseract_geometry::GeometryType::CYLINDER:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Cylinder>(geom));
    }
    case tesseract_geometry::GeometryType::CONE:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Cone>(geom));
    }
    case tesseract_geometry::GeometryType::CAPSULE:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Capsule>(geom));
    }
    case tesseract_geometry::GeometryType::MESH:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Mesh>(geom));
    }
    case tesseract_geometry::GeometryType::CONVEX_MESH:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::ConvexMesh>(geom));
    }
    case tesseract_geometry::GeometryType::OCTREE:
    {
      return createShapePrimitive(std::static_pointer_cast<const tesseract_geometry::Octree>(geom));
    }
    case tesseract_geometry::GeometryType::COMPOUND_MESH:
    {
      throw std::runtime_error("CompundMesh type should not be passed to this function!");
    }
    default:
    {
      CONSOLE_BRIDGE_logError("This geometric shape type (%d) is not supported using Coal yet",
                              static_cast<int>(geom->getType()));
      return nullptr;
    }
  }
}

CollisionGeometryPtr createShapePrimitive(const CollisionShapeConstPtr& geom)
{
  CollisionGeometryPtr shape = CoalCollisionGeometryCache::get(geom);
  if (shape != nullptr)
    return shape;

  shape = createShapePrimitiveHelper(geom);
  CoalCollisionGeometryCache::insert(geom, shape);
  return shape;
}

constexpr double COAL_SUPPORT_FUNC_TOLERANCE = 0.01;
constexpr double COAL_LENGTH_TOLERANCE = 0.001;

/**
 * @brief Compute the support point for a COAL shape along a direction.
 *
 * Delegates to coal::details::getSupport which handles all shape types
 * including ConvexBase32 (returns a single extreme vertex).
 */
void GetSupport(const coal::ShapeBase* shape,
                const coal::Vec3s& normal,
                double& out_support,
                coal::Vec3s& out_point)
{
  int hint = 0;
  out_point = coal::details::getSupport(shape, normal, hint);
  out_support = normal.dot(out_point);
}

inline bool needsCollisionCheck(const CollisionObjectWrapper* cd1,
                                const CollisionObjectWrapper* cd2,
                                const std::shared_ptr<const tesseract_common::ContactAllowedValidator>& validator,
                                bool verbose)
{
  return cd1->m_enabled && cd2->m_enabled && (cd2->m_collisionFilterGroup & cd1->m_collisionFilterMask) &&  // NOLINT
         (cd1->m_collisionFilterGroup & cd2->m_collisionFilterMask) &&                                      // NOLINT
         !isContactAllowed(cd1->getName(), cd2->getName(), validator, verbose);
}

/**
 * @brief Populate continuous collision fields (cc_time, cc_type, cc_transform) on a ContactResult.
 *
 * Uses support-function-based approach matching Bullet's calculateContinuousData:
 * finds the shape's extreme points along the contact normal at t=0 and t=1, then
 * classifies the collision time based on which pose has greater support.
 *
 * For objects that are not CastHullShapes (static objects), the fields are left at their defaults.
 */
void populateContinuousCollisionFields(ContactResult& contact,
                                       const coal::CollisionObject* o1,
                                       const coal::CollisionObject* o2)
{
  const std::array<const coal::CollisionObject*, 2> objects = { o1, o2 };
  const Eigen::Vector3d pt_world = (contact.nearest_points[0] + contact.nearest_points[1]) / 2.0;

  for (std::size_t i = 0; i < 2; ++i)
  {
    const auto* cast_shape = dynamic_cast<const CastHullShape*>(objects[i]->collisionGeometry().get());
    if (cast_shape == nullptr)
      continue;

    // cc_transform = pose2 = pose1 * relative_motion
    // castTransform stores tf1.inverseTimes(tf2), so pose1 * castTransform = pose2
    const auto& ct = cast_shape->getCastTransform();
    Eigen::Isometry3d cast_eigen;
    cast_eigen.linear() = ct.getRotation();
    cast_eigen.translation() = ct.getTranslation();
    contact.cc_transform[i] = contact.transform[i] * cast_eigen;

    // Shape world transforms at t=0 and t=1
    const coal::Transform3s& tf_world0 = objects[i]->getTransform();
    coal::Transform3s tf_world1 = tf_world0 * ct;

    // Normal pointing from current object toward the other (matching Bullet convention)
    // COAL contact.normal points from o1 to o2; Bullet's m_normalWorldOnB points from o2 to o1
    const coal::Vec3s normal_world = (i == 0) ? coal::Vec3s(contact.normal) : coal::Vec3s(-contact.normal);

    // Transform normal into local frames at t=0 and t=1
    coal::Vec3s normal_local0 = tf_world0.getRotation().transpose() * normal_world;
    coal::Vec3s normal_local1 = tf_world1.getRotation().transpose() * normal_world;

    // Get support points on the underlying shape at both local normals
    const coal::ShapeBase* underlying = cast_shape->getUnderlyingShape().get();
    coal::Vec3s pt_local0;
    double sup_local0 = 0;
    GetSupport(underlying, normal_local0, sup_local0, pt_local0);
    coal::Vec3s pt_world0 = tf_world0.transform(pt_local0);

    coal::Vec3s pt_local1;
    double sup_local1 = 0;
    GetSupport(underlying, normal_local1, sup_local1, pt_local1);
    coal::Vec3s pt_world1 = tf_world1.transform(pt_local1);

    // Compare support projections along the contact normal
    double shape_sup0 = normal_world.dot(pt_world0);
    double shape_sup1 = normal_world.dot(pt_world1);

    if (shape_sup0 - shape_sup1 > COAL_SUPPORT_FUNC_TOLERANCE)
    {
      contact.cc_time[i] = 0;
      contact.cc_type[i] = ContinuousCollisionType::CCType_Time0;
    }
    else if (shape_sup1 - shape_sup0 > COAL_SUPPORT_FUNC_TOLERANCE)
    {
      contact.cc_time[i] = 1;
      contact.cc_type[i] = ContinuousCollisionType::CCType_Time1;
    }
    else
    {
      // Between: interpolate based on distances from average contact point to support points
      double l0c = (pt_world - Eigen::Vector3d(pt_world0)).norm();
      double l1c = (pt_world - Eigen::Vector3d(pt_world1)).norm();

      // Update nearest_points_local to averaged support point (matching Bullet)
      Eigen::Isometry3d link_tf_inv = contact.transform[i].inverse();
      Eigen::Isometry3d shape_tf0;
      shape_tf0.linear() = tf_world0.getRotation();
      shape_tf0.translation() = tf_world0.getTranslation();
      contact.nearest_points_local[i] = link_tf_inv * (shape_tf0 * (((pt_local0 + pt_local1) / 2.0).eval()));

      contact.cc_type[i] = ContinuousCollisionType::CCType_Between;

      if (l0c + l1c < COAL_LENGTH_TOLERANCE)
        contact.cc_time[i] = 0.5;
      else
        contact.cc_time[i] = l0c / (l0c + l1c);
    }
  }
}

/**
 * @brief Custom GetSupportFunction for CastHullShape collisions.
 *
 * Uses the Schulman et al. (2013) approach: the support of a swept shape is
 * max(support_start(d), support_end(d)) using the underlying shape's exact
 * support function. Shape0 is always the CastHullShape; shape1 may also be one.
 */
void castHullGetSupportFunc(const coal::details::MinkowskiDiff& md,
                            const coal::Vec3s& dir,
                            coal::Vec3s& support0,
                            coal::Vec3s& support1,
                            coal::support_func_guess_t& hint,
                            coal::details::ShapeSupportData data[2])
{
  // Shape0 is the CastHullShape — use the Schulman support function
  const auto* cast_hull0 = static_cast<const CastHullShape*>(md.shapes[0]);
  coal::details::getShapeSupport<coal::details::SupportOptions::NoSweptSphere>(
      cast_hull0, dir, support0, hint[0], data[0]);

  // Negate direction for shape1 per Minkowski difference convention
  // (coal's getSupportTpl uses -dir for shape1: support1 = s_S1(-dir))
  coal::Vec3s neg_dir1 = -(md.oR1.transpose() * dir);

  // Check if shape1 is also a CastHullShape — if so, use Schulman support
  // directly to avoid incorrect dispatch through coal's generic getSupport
  // (CastHullShape is a ConvexBase32 which coal's dispatcher may mishandle)
  const auto* cast_hull1 = dynamic_cast<const CastHullShape*>(md.shapes[1]);
  if (cast_hull1 != nullptr)
  {
    coal::details::getShapeSupport<coal::details::SupportOptions::NoSweptSphere>(
        cast_hull1, neg_dir1, support1, hint[1], data[1]);
  }
  else
  {
    support1 = coal::details::getSupport<coal::details::SupportOptions::NoSweptSphere>(
        md.shapes[1], neg_dir1, hint[1]);
  }
  support1 = md.oR1 * support1 + md.ot1;
}

/**
 * @brief Perform narrowphase collision for CastHullShape using custom support functions.
 *
 * Bypasses coal::ComputeCollision to inject the Schulman support function into
 * the MinkowskiDiff, then runs GJK/EPA directly.
 *
 * Uses DefaultGJK (no acceleration) because PolyakAcceleration's momentum
 * interferes with the discontinuous Schulman support function, which switches
 * between start and end configurations.
 *
 * @param o1 First collision object
 * @param o2 Second collision object
 * @param request The collision request parameters
 * @param result Output collision result
 * @return true if collision was detected
 */
bool castHullCollide(coal::CollisionObject* o1,
                     coal::CollisionObject* o2,
                     const coal::CollisionRequest& request,
                     coal::CollisionResult& result)
{
  // Determine which object is the CastHullShape; ensure it's shape0
  const auto* geom0 = o1->collisionGeometry().get();
  const auto* geom1 = o2->collisionGeometry().get();
  const auto* cast0 = dynamic_cast<const CastHullShape*>(geom0);
  const auto* cast1 = dynamic_cast<const CastHullShape*>(geom1);

  bool swapped = false;
  const coal::ShapeBase* shape0 = nullptr;
  const coal::ShapeBase* shape1 = nullptr;
  coal::Transform3s tf0, tf1;

  if (cast0 != nullptr)
  {
    shape0 = static_cast<const coal::ShapeBase*>(cast0);
    shape1 = static_cast<const coal::ShapeBase*>(geom1);
    tf0 = o1->getTransform();
    tf1 = o2->getTransform();
  }
  else
  {
    // cast1 must be CastHullShape; swap so CastHullShape is shape0
    swapped = true;
    shape0 = static_cast<const coal::ShapeBase*>(cast1);
    shape1 = static_cast<const coal::ShapeBase*>(geom0);
    tf0 = o2->getTransform();
    tf1 = o1->getTransform();
  }

  // Set up MinkowskiDiff with standard transform initialization
  coal::details::MinkowskiDiff md;
  md.set<coal::details::SupportOptions::NoSweptSphere>(shape0, shape1, tf0, tf1);

  // Replace the support function with our custom CastHullShape support
  md.getSupportFunc = castHullGetSupportFunc;

  // Configure GJK — use DefaultGJK (no acceleration).
  // PolyakAcceleration's momentum interferes with the discontinuous
  // Schulman support, causing convergence failures for shapes with
  // continuous support functions (e.g. Sphere, Cylinder, Capsule).
  coal::details::GJK gjk(request.gjk_max_iterations, request.gjk_tolerance);
  gjk.setDistanceEarlyBreak(request.distance_upper_bound);

  // Compute initial guess — use the center difference between shapes as a
  // better starting direction than an arbitrary axis.
  coal::Vec3s guess = tf0.getTranslation() - tf1.getTranslation();
  if (guess.squaredNorm() < 1e-12)
    guess = coal::Vec3s(1, 0, 0);
  coal::support_func_guess_t support_hint = coal::support_func_guess_t::Zero();
  if (request.gjk_initial_guess == coal::CachedGuess)
  {
    guess = request.cached_gjk_guess;
    support_hint = request.cached_support_func_guess;
  }

  // Run GJK
  gjk.evaluate(md, guess.cast<coal::SolverScalar>(), support_hint);

  // Cache GJK guess for subsequent calls
  result.cached_gjk_guess = gjk.ray.cast<coal::Scalar>();
  result.cached_support_func_guess = gjk.support_hint;

  // If no collision, return
  if (gjk.status != coal::details::GJK::Collision &&
      gjk.status != coal::details::GJK::CollisionWithPenetrationInformation)
  {
    return false;
  }

  // Extract witness points from GJK if it already has penetration info
  coal::Scalar distance;
  coal::Vec3s p1, p2, normal;

  if (gjk.status == coal::details::GJK::CollisionWithPenetrationInformation)
  {
    coal::Vec3ps p1_, p2_, normal_;
    gjk.getWitnessPointsAndNormal(md, p1_, p2_, normal_);
    p1 = p1_.cast<coal::Scalar>();
    p2 = p2_.cast<coal::Scalar>();
    normal = normal_.cast<coal::Scalar>();
    distance = gjk.distance;
  }
  else if (request.enable_contact)
  {
    // Run EPA for penetration depth
    coal::details::EPA epa(request.epa_max_iterations, request.epa_tolerance);
    epa.evaluate(gjk, (-guess).cast<coal::SolverScalar>());

    result.cached_gjk_guess = -(epa.depth * epa.normal).cast<coal::Scalar>();
    result.cached_support_func_guess = epa.support_hint;

    distance = std::min(coal::Scalar(0), -coal::Scalar(epa.depth));
    coal::Vec3ps p1_, p2_, normal_;
    epa.getWitnessPointsAndNormal(md, p1_, p2_, normal_);
    p1 = p1_.cast<coal::Scalar>();
    p2 = p2_.cast<coal::Scalar>();
    normal = normal_.cast<coal::Scalar>();
  }
  else
  {
    // Collision detected but no penetration info requested
    distance = -(std::numeric_limits<coal::Scalar>::max)();
    p1 = p2 = normal = coal::Vec3s::Constant(std::numeric_limits<coal::Scalar>::quiet_NaN());
  }

  // Transform witness points to world frame (same pattern as coal's EPAExtractWitnessPointsAndNormal)
  coal::Vec3s p = tf0.transform(0.5 * (p1 + p2));
  normal = tf0.getRotation() * normal;
  p1.noalias() = p - 0.5 * distance * normal;
  p2.noalias() = p + 0.5 * distance * normal;

  // Create contact
  coal::Contact contact;
  contact.o1 = o1->collisionGeometry().get();
  contact.o2 = o2->collisionGeometry().get();
  contact.penetration_depth = distance;

  if (swapped)
  {
    // Swap witness points and negate normal to match original ordering
    contact.nearest_points[0] = p2;
    contact.nearest_points[1] = p1;
    contact.normal = -normal;
  }
  else
  {
    contact.nearest_points[0] = p1;
    contact.nearest_points[1] = p2;
    contact.normal = normal;
  }

  result.addContact(contact);
  return true;
}

bool CollisionCallback::collide(coal::CollisionObject* o1, coal::CollisionObject* o2)
{
  if (cdata->done)
    return true;

  const auto* cd1 = static_cast<const CollisionObjectWrapper*>(o1->getUserData());
  const auto* cd2 = static_cast<const CollisionObjectWrapper*>(o2->getUserData());

  if (!needsCollisionCheck(cd1, cd2, cdata->validator, false))
    return false;

  std::size_t num_contacts = (cdata->req.contact_limit > 0) ? static_cast<std::size_t>(cdata->req.contact_limit) :
                                                              std::numeric_limits<std::size_t>::max();
  if (cdata->req.type == ContactTestType::FIRST)
    num_contacts = 1;
  const auto security_margin = cdata->collision_margin_data.getCollisionMargin(cd1->getName(), cd2->getName());

  // Detect CastHullShape and use custom narrowphase with Schulman support function
  const bool is_cast_hull = (dynamic_cast<const CastHullShape*>(o1->collisionGeometry().get()) != nullptr) ||
                            (dynamic_cast<const CastHullShape*>(o2->collisionGeometry().get()) != nullptr);

  coal::CollisionResult col_result;

  if (is_cast_hull)
  {
    // CastHullShape path: use custom GJK with Schulman support function.
    // Uses DefaultGJK (no acceleration) because PolyakAcceleration's momentum
    // interferes with the discontinuous Schulman support, causing convergence
    // failures for shapes with continuous support functions (Sphere, etc.).
    CollisionObjectPair object_pair = std::make_pair(o1, o2);
    auto col_request_it = cdata->collision_cache->find(object_pair);
    if (col_request_it == cdata->collision_cache->end())
    {
      coal::CollisionRequest col_request;
      // Use DefaultGJK — do NOT set gjk_variant, convergence_criterion, or
      // convergence_criterion_type, leaving them at coal defaults.
      col_request.gjk_initial_guess = coal::BoundingVolumeGuess;
      col_request.enable_contact = cdata->req.calculate_penetration;
      col_request.num_max_contacts = num_contacts;
      col_request.security_margin = security_margin;
      // CastHullShape swept volumes can be much larger than security_margin.
      // Using security_margin as the GJK early-break threshold would cause GJK to
      // return NoCollisionEarlyStopped before the simplex encloses the origin.
      col_request.distance_upper_bound = (std::numeric_limits<coal::Scalar>::max)();
      // Dummy ComputeCollision functor (not used for CastHullShape, but needed for cache type)
      auto col_functor = coal::ComputeCollision(o1->collisionGeometryPtr(), o2->collisionGeometryPtr());
      col_request_it =
          cdata->collision_cache->try_emplace(object_pair, std::move(col_functor), std::move(col_request)).first;
    }
    else
    {
      auto& cached_request = col_request_it->second.second;
      cached_request.enable_contact = cdata->req.calculate_penetration;
      cached_request.num_max_contacts = num_contacts;
      cached_request.security_margin = security_margin;
      cached_request.distance_upper_bound = (std::numeric_limits<coal::Scalar>::max)();
    }

    auto& cached_request = col_request_it->second.second;
    castHullCollide(o1, o2, cached_request, col_result);

    if (cached_request.gjk_initial_guess != coal::CachedGuess)
    {
      cached_request.gjk_initial_guess = coal::CachedGuess;
      cached_request.cached_gjk_guess = col_result.cached_gjk_guess;
      cached_request.cached_support_func_guess = col_result.cached_support_func_guess;
    }
  }
  else
  {
    // Standard path for non-CastHullShape: use coal's ComputeCollision functor
    CollisionObjectPair object_pair = std::make_pair(o1, o2);
    auto col_request_it = cdata->collision_cache->find(object_pair);
    if (col_request_it == cdata->collision_cache->end())
    {
      coal::CollisionRequest col_request;
      col_request.gjk_variant = coal::GJKVariant::PolyakAcceleration;
      col_request.gjk_convergence_criterion = coal::GJKConvergenceCriterion::DualityGap;
      col_request.gjk_convergence_criterion_type = coal::GJKConvergenceCriterionType::Absolute;
      col_request.gjk_initial_guess = coal::BoundingVolumeGuess;
      col_request.enable_contact = cdata->req.calculate_penetration;
      col_request.num_max_contacts = num_contacts;
      col_request.security_margin = security_margin;
      col_request.distance_upper_bound = security_margin + col_request.gjk_tolerance;
      auto col_functor = coal::ComputeCollision(o1->collisionGeometryPtr(), o2->collisionGeometryPtr());
      col_request_it =
          cdata->collision_cache->try_emplace(object_pair, std::move(col_functor), std::move(col_request)).first;
    }
    else
    {
      auto& cached_request = col_request_it->second.second;
      cached_request.enable_contact = cdata->req.calculate_penetration;
      cached_request.num_max_contacts = num_contacts;
      cached_request.security_margin = security_margin;
      cached_request.distance_upper_bound = security_margin + cached_request.gjk_tolerance;
    }

    auto& [functor, cached_request] = col_request_it->second;
    functor(o1->getTransform(), o2->getTransform(), cached_request, col_result);

    if (cached_request.gjk_initial_guess != coal::CachedGuess)
    {
      cached_request.gjk_initial_guess = coal::CachedGuess;
      cached_request.cached_gjk_guess = col_result.cached_gjk_guess;
      cached_request.cached_support_func_guess = col_result.cached_support_func_guess;
    }
  }

  if (!col_result.isCollision())
    return false;

  TESSERACT_THREAD_LOCAL tesseract_common::LinkNamesPair link_pair;
  tesseract_common::makeOrderedLinkPair(link_pair, cd1->getName(), cd2->getName());

  const Eigen::Isometry3d& tf1 = cd1->getCollisionObjectsTransform();
  const Eigen::Isometry3d& tf2 = cd2->getCollisionObjectsTransform();
  Eigen::Isometry3d tf1_inv = tf1.inverse();
  Eigen::Isometry3d tf2_inv = tf2.inverse();

  for (size_t i = 0; i < col_result.numContacts(); ++i)
  {
    const coal::Contact& coal_contact = col_result.getContact(i);
    ContactResult contact;
    contact.link_names[0] = cd1->getName();
    contact.link_names[1] = cd2->getName();
    contact.shape_id[0] = CollisionObjectWrapper::getShapeIndex(o1);
    contact.shape_id[1] = CollisionObjectWrapper::getShapeIndex(o2);
    contact.subshape_id[0] = static_cast<int>(coal_contact.b1);
    contact.subshape_id[1] = static_cast<int>(coal_contact.b2);
    contact.nearest_points[0] = coal_contact.nearest_points[0];
    contact.nearest_points[1] = coal_contact.nearest_points[1];
    contact.nearest_points_local[0] = tf1_inv * contact.nearest_points[0];
    contact.nearest_points_local[1] = tf2_inv * contact.nearest_points[1];
    contact.transform[0] = tf1;
    contact.transform[1] = tf2;
    contact.type_id[0] = cd1->getTypeID();
    contact.type_id[1] = cd2->getTypeID();
    contact.distance = coal_contact.penetration_depth;
    contact.normal = coal_contact.normal;

    populateContinuousCollisionFields(contact, o1, o2);
    const auto it = cdata->res->find(link_pair);
    const bool found = (it != cdata->res->end() && !it->second.empty());

    processResult(*cdata, contact, link_pair, found);
  }

  return cdata->done;
}

CollisionObjectWrapper::CollisionObjectWrapper(std::string name,
                                               const int& type_id,
                                               CollisionShapesConst shapes,
                                               tesseract_common::VectorIsometry3d shape_poses)
  : name_(std::move(name)), type_id_(type_id), shapes_(std::move(shapes)), shape_poses_(std::move(shape_poses))
{
  assert(!shapes_.empty());                       // NOLINT
  assert(!shape_poses_.empty());                  // NOLINT
  assert(!name_.empty());                         // NOLINT
  assert(shapes_.size() == shape_poses_.size());  // NOLINT

  m_collisionFilterGroup = CollisionFilterGroups::KinematicFilter;
  m_collisionFilterMask = CollisionFilterGroups::StaticFilter | CollisionFilterGroups::KinematicFilter;

  collision_geometries_.reserve(shapes_.size());
  collision_objects_.reserve(shapes_.size());
  collision_objects_raw_.reserve(shapes_.size());
  for (std::size_t i = 0; i < shapes_.size(); ++i)  // NOLINT
  {
    if (shapes_[i]->getType() == tesseract_geometry::GeometryType::COMPOUND_MESH)
    {
      const auto& meshes = std::static_pointer_cast<const tesseract_geometry::CompoundMesh>(shapes_[i])->getMeshes();
      for (const auto& mesh : meshes)
      {
        const CollisionGeometryPtr subshape = createShapePrimitive(mesh);
        if (subshape != nullptr)
        {
          collision_geometries_.push_back(subshape);
          auto co = std::make_shared<CoalCollisionObjectWrapper>(subshape);
          co->setUserData(this);
          co->setShapeIndex(static_cast<int>(i));
          co->setTransform(coal::Transform3s(shape_poses_[i].rotation(), shape_poses_[i].translation()));
          co->updateAABB();
          collision_objects_.push_back(co);
          collision_objects_raw_.push_back(co.get());
        }
      }
    }
    else
    {
      const CollisionGeometryPtr subshape = createShapePrimitive(shapes_[i]);
      if (subshape != nullptr)
      {
        collision_geometries_.push_back(subshape);
        auto co = std::make_shared<CoalCollisionObjectWrapper>(subshape);
        co->setUserData(this);
        co->setShapeIndex(static_cast<int>(i));
        co->setTransform(coal::Transform3s(shape_poses_[i].rotation(), shape_poses_[i].translation()));
        co->updateAABB();
        collision_objects_.push_back(co);
        collision_objects_raw_.push_back(co.get());
      }
    }
  }
}

int CollisionObjectWrapper::getShapeIndex(const coal::CollisionObject* co)
{
  return static_cast<const CoalCollisionObjectWrapper*>(co)->getShapeIndex();
}

}  // namespace tesseract_collision::tesseract_collision_coal
