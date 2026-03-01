/**
 * @file coal_casthullshape.cpp
 * @brief CastHullShape implementation.
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
#include <coal/broadphase/broadphase_collision_manager.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_collision/coal/coal_casthullshape.h>
#include <tesseract_collision/coal/coal_casthullshape_utility.h>

namespace tesseract_collision::tesseract_collision_coal
{
CastHullShape::CastHullShape(std::shared_ptr<coal::ShapeBase> shape, const coal::Transform3s& castTransform)
  : shape_(std::move(shape)), castTransform_(castTransform), castTransformInv_(coal::Transform3s(castTransform).inverse())
{
  computeLocalAABB();
}

void CastHullShape::computeLocalAABB()
{
  coal::Transform3s tf = coal::Transform3s::Identity();
  coal::AABB aabb;
  coal::computeBV<coal::AABB, CastHullShape>(*this, tf, aabb);

  aabb_local = aabb;
  aabb_center = aabb_local.center();
  aabb_radius = (aabb_local.min_ - aabb_center).norm();
}

// NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
CastHullShape* CastHullShape::clone() const { return new CastHullShape(shape_, castTransform_); }

double CastHullShape::computeVolume() const
{
  double baseVolume = shape_->computeVolume();

  coal::Vec3s translation = castTransform_.getTranslation();
  double translation_length = translation.norm();

  coal::Matrix3s rotation = castTransform_.getRotation();
  bool has_rotation = !rotation.isIdentity(1e-6);

  if (translation_length < 1e-6 && !has_rotation)
    return baseVolume;

  coal::Transform3s identity = coal::Transform3s::Identity();
  coal::AABB swept_aabb;
  coal::computeBV<coal::AABB, CastHullShape>(*this, identity, swept_aabb);

  return std::max(baseVolume, swept_aabb.volume());
}

bool CastHullShape::isEqual(const coal::CollisionGeometry& _other) const
{
  const auto* other_ptr = dynamic_cast<const CastHullShape*>(&_other);
  if (other_ptr == nullptr)
    return false;

  return shape_ == other_ptr->shape_ && castTransform_ == other_ptr->castTransform_;
}

void CastHullShape::updateCastTransform(const coal::Transform3s& castTransform)
{
  castTransform_ = castTransform;
  castTransformInv_ = coal::Transform3s(castTransform).inverse();
  computeLocalAABB();
}

}  // namespace tesseract_collision::tesseract_collision_coal
