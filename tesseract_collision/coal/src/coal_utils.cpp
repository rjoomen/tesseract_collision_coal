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
constexpr double COAL_EPSILON = 1e-6;

/**
 * @brief Compute the average support point for a shape along a direction.
 *
 * Matches Bullet's GetAverageSupport algorithm: for polyhedral shapes
 * (ConvexBase32), iterates all vertices, finds the maximum support value,
 * and averages all vertices within COAL_EPSILON of that maximum. This
 * handles degenerate cases where multiple vertices have equal support.
 *
 * For non-polyhedral shapes (Sphere, Capsule, etc.), falls back to coal's
 * getSupport<WithSweptSphere> which includes the shape's radius.
 */
void GetAverageSupport(const coal::ShapeBase* shape,
                       const coal::Vec3s& localNormal,
                       double& outsupport,
                       coal::Vec3s& outpt)
{
  const auto* convex = dynamic_cast<const coal::ConvexBase32*>(shape);
  if (convex != nullptr && convex->points && !convex->points->empty())
  {
    coal::Vec3s ptSum = coal::Vec3s::Zero();
    double ptCount = 0;
    double maxSupport = -1e10;

    for (const auto& pt : *convex->points)
    {
      double sup = pt.dot(localNormal);
      if (sup > maxSupport + COAL_EPSILON)
      {
        ptCount = 1;
        ptSum = pt;
        maxSupport = sup;
      }
      else if (sup >= maxSupport - COAL_EPSILON)
      {
        ptCount += 1;
        ptSum += pt;
      }
    }
    outsupport = maxSupport;
    outpt = ptSum / ptCount;
  }
  else
  {
    // For primitive shapes (Sphere, etc.), use coal's standard support function.
    // WithSweptSphere mode is required so that Sphere returns radius*normalize(dir)
    // instead of zero (NoSweptSphere treats Sphere as a point + inflation).
    int hint = 0;
    outpt = coal::details::getSupport<coal::details::SupportOptions::WithSweptSphere>(shape, localNormal, hint);
    outsupport = localNormal.dot(outpt);
  }
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
 * Matches Bullet's calculateContinuousData algorithm:
 * finds the shape's extreme points along the contact normal at t=0 and t=1, then
 * classifies the collision time based on which pose has greater support.
 *
 * For the CCType_Between case:
 * - cc_time is computed by projecting the GJK witness point onto the shape's
 *   center trajectory (center0 → center1).
 * - nearest_points_local uses the average of the support points at t=0 and t=1,
 *   matching Bullet's approach.
 *
 * For objects that are not CastHullShapes (static objects), the fields are left at their defaults.
 */
void populateContinuousCollisionFields(ContactResult& contact,
                                       const coal::CollisionObject* o1,
                                       const coal::CollisionObject* o2)
{
  const std::array<const coal::CollisionObject*, 2> objects = { o1, o2 };

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

    // Normal pointing from current object toward the other (matching Bullet convention:
    // for link 0, Bullet uses -m_normalWorldOnB = A→B direction;
    // for link 1, Bullet uses m_normalWorldOnB = A→B = toward link 1's "other")
    // COAL contact.normal points from o1 to o2, same as Bullet's m_normalWorldOnB
    const coal::Vec3s normal_world = (i == 0) ? coal::Vec3s(contact.normal) : coal::Vec3s(-contact.normal);

    // Transform normal into local frames at t=0 and t=1
    coal::Vec3s normal_local0 = tf_world0.getRotation().transpose() * normal_world;
    coal::Vec3s normal_local1 = tf_world1.getRotation().transpose() * normal_world;

    // Get support points on the underlying shape at both local normals
    const coal::ShapeBase* underlying = cast_shape->getUnderlyingShape().get();
    coal::Vec3s pt_local0;
    double sup_local0 = 0;
    GetAverageSupport(underlying, normal_local0, sup_local0, pt_local0);
    coal::Vec3s pt_world0 = tf_world0.transform(pt_local0);

    coal::Vec3s pt_local1;
    double sup_local1 = 0;
    GetAverageSupport(underlying, normal_local1, sup_local1, pt_local1);
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
      contact.cc_type[i] = ContinuousCollisionType::CCType_Between;

      // Use the GJK witness point on THIS shape (matching Bullet's pt_world convention)
      Eigen::Vector3d pt_on_shape = Eigen::Vector3d(contact.nearest_points[i]);

      // Interpolate cc_time by projecting the witness point onto the shape's
      // center trajectory (center0 → center1). The witness point's position
      // along the sweep encodes the collision time.
      Eigen::Vector3d center0 = Eigen::Vector3d(tf_world0.getTranslation());
      Eigen::Vector3d center1 = Eigen::Vector3d(tf_world1.getTranslation());
      Eigen::Vector3d sweep = center1 - center0;
      double sweep_sq = sweep.squaredNorm();

      if (sweep_sq < COAL_LENGTH_TOLERANCE * COAL_LENGTH_TOLERANCE)
        contact.cc_time[i] = 0.5;
      else
        contact.cc_time[i] = sweep.dot(pt_on_shape - center0) / sweep_sq;

      // nearest_points_local: average of the two local support points,
      // transformed to world via shape_tf0, then to link-local coordinates.
      // Matches Bullet's calculateContinuousData: (shape_ptLocal0 + shape_ptLocal1) / 2.0
      coal::Vec3s avg_pt_local = (pt_local0 + pt_local1) / 2.0;
      Eigen::Isometry3d link_tf_inv = contact.transform[i].inverse();
      Eigen::Isometry3d shape_tf0;
      shape_tf0.linear() = tf_world0.getRotation();
      shape_tf0.translation() = tf_world0.getTranslation();
      contact.nearest_points_local[i] = link_tf_inv * (shape_tf0 * Eigen::Vector3d(avg_pt_local));
    }
  }
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

  // Detect CastHullShape for distance_upper_bound and post-processing
  const bool is_cast_hull = (dynamic_cast<const CastHullShape*>(o1->collisionGeometry().get()) != nullptr) ||
                            (dynamic_cast<const CastHullShape*>(o2->collisionGeometry().get()) != nullptr);

  coal::CollisionResult col_result;

  // Use coal's standard ComputeCollision for all shapes, including CastHullShape.
  // CastHullShape is a ConvexBase32 subclass with proper convex hull vertices and
  // neighbor data, so coal's standard GJK operates on the hull vertices directly.
  // This matches Bullet's approach (btConvexHullShape with tessellated vertices)
  // and produces consistent witness points for cc_time computation.
  {
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
      if (is_cast_hull)
      {
        // CastHullShape swept volumes can be much larger than security_margin.
        // Using security_margin as the GJK early-break threshold would cause GJK to
        // return NoCollisionEarlyStopped before the simplex encloses the origin.
        col_request.distance_upper_bound = (std::numeric_limits<coal::Scalar>::max)();
      }
      else
      {
        col_request.distance_upper_bound = security_margin + col_request.gjk_tolerance;
      }
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
      if (is_cast_hull)
        cached_request.distance_upper_bound = (std::numeric_limits<coal::Scalar>::max)();
      else
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
