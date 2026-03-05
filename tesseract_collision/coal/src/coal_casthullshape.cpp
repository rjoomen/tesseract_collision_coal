/**
 * @file coal_casthullshape.cpp
 * @brief Tesseract Coal Utility Functions.
 *
 * @author Roelof Oomen
 * @date Aug 04, 2025
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
#include <coal/narrowphase/support_functions.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_collision/coal/coal_casthullshape.h>

namespace tesseract_collision::tesseract_collision_coal
{
CastHullShape::CastHullShape(std::shared_ptr<coal::ShapeBase> shape, const coal::Transform3s& castTransform)
  : shape_(std::move(shape))
  , castTransform_(castTransform)
  , castTransformInv_(coal::Transform3s(castTransform).inverse())
{
}

void CastHullShape::computeLocalAABB()
{
  const coal::AABB aabb = computeSweptAABB();
  aabb_local = aabb;
  aabb_center = aabb_local.center();
  aabb_radius = (aabb_local.max_ - aabb_local.min_).norm() * coal::Scalar(0.5);
}

double CastHullShape::computeVolume() const
{
  // Get the base volume
  double baseVolume = shape_->computeVolume();

  // Check if transform is identity (no sweeping needed)
  coal::Vec3s translation = castTransform_.getTranslation();
  double translation_length = translation.norm();

  // Check for rotation component
  coal::Matrix3s rotation = castTransform_.getRotation();
  bool has_rotation = !rotation.isIdentity(1e-6);

  // If no transformation, return base volume
  if (translation_length < 1e-6 && !has_rotation)
  {
    return baseVolume;
  }

  // Approximate using the swept AABB volume
  const coal::AABB swept_aabb = computeSweptAABB();
  return std::max(baseVolume, swept_aabb.volume());
}

bool CastHullShape::isEqual(const coal::CollisionGeometry& _other) const
{
  const auto* other_ptr = dynamic_cast<const CastHullShape*>(&_other);
  if (other_ptr == nullptr)
    return false;
  const CastHullShape& other = *other_ptr;

  return shape_ == other.shape_ && castTransform_ == other.castTransform_;
}

void CastHullShape::updateCastTransform(const coal::Transform3s& castTransform)
{
  castTransform_ = castTransform;
  castTransformInv_ = coal::Transform3s(castTransform).inverse();
  // Recompute aabb_local so that subsequent calls to computeAABB() on the
  // enclosing CollisionObject produce a correct world-space AABB.
  computeLocalAABB();
}

void CastHullShape::computeShapeSupport(const coal::Vec3s& dir,
                                        coal::Vec3s& support,
                                        int& hint,
                                        coal::details::ShapeSupportData& /*data*/) const
{
  const coal::ShapeBase* shape = shape_.get();

  // Support at pose 0 (identity — shape in its own local frame)
  const coal::Vec3s s0 =
      coal::details::getSupport<coal::details::SupportOptions::NoSweptSphere>(shape, dir, hint);

  // Support at pose 1 (cast transform applied).
  // Transform the query direction into the local frame of the end pose,
  // compute the support there, then transform the result back.
  const coal::Vec3s dir_cast = castTransformInv_.getRotation() * dir;
  const coal::Vec3s s1_local =
      coal::details::getSupport<coal::details::SupportOptions::NoSweptSphere>(shape, dir_cast, hint);
  const coal::Vec3s s1 = castTransform_.transform(s1_local);

  // Return whichever endpoint is furthest along dir.
  support = (dir.dot(s0) >= dir.dot(s1)) ? s0 : s1;
}

coal::AABB CastHullShape::computeSweptAABB() const
{
  int hint = 0;
  coal::details::ShapeSupportData data;
  coal::Vec3s support;
  coal::Vec3s min_pt;
  coal::Vec3s max_pt;

  for (int i = 0; i < 3; ++i)
  {
    coal::Vec3s dir = coal::Vec3s::Zero();
    dir[i] = coal::Scalar(1);
    computeShapeSupport(dir, support, hint, data);
    max_pt[i] = support[i];
    computeShapeSupport(-dir, support, hint, data);
    min_pt[i] = support[i];
  }

  coal::AABB aabb;
  aabb.min_ = min_pt;
  aabb.max_ = max_pt;
  return aabb;
}

}  // namespace tesseract_collision::tesseract_collision_coal
