/**
 * @file coal_casthullshape.h
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

#ifndef TESSERACT_COLLISION_COAL_CASTHULLSHAPE_H
#define TESSERACT_COLLISION_COAL_CASTHULLSHAPE_H

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <memory>
#include <coal/shape/geometric_shapes.h>
#include <coal/narrowphase/support_data.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

namespace tesseract_collision::tesseract_collision_coal
{
class CastHullShape : public coal::ShapeBase
{
public:
  CastHullShape(std::shared_ptr<coal::ShapeBase> shape, const coal::Transform3s& castTransform);

  void computeLocalAABB() override;

  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  coal::ShapeBase* clone() const override { return new CastHullShape(*this); }

  double computeVolume() const override;

  bool isEqual(const coal::CollisionGeometry& _other) const override;

  void updateCastTransform(const coal::Transform3s& castTransform);

  const std::shared_ptr<coal::ShapeBase>& getUnderlyingShape() const { return shape_; }

  const coal::Transform3s& getCastTransform() const { return castTransform_; }

  const coal::Transform3s& getCastTransformInverse() const { return castTransformInv_; }

  /// @brief GJK/EPA support function (required for Coal native collision/distance).
  /// @param dir support direction, always unit-length when called by Coal.
  /// @param support output: the point in the swept shape maximising dot(p, dir).
  /// @param hint warm-start hint for convex shapes.
  /// @param data temporary data for support computation.
  /// @note Do NOT add the swept-sphere radius here; Coal applies it externally.
  void computeShapeSupport(const coal::Vec3s& dir,
                           coal::Vec3s& support,
                           int& hint,
                           coal::details::ShapeSupportData& data) const override;

private:
  std::shared_ptr<coal::ShapeBase> shape_;
  coal::Transform3s castTransform_;
  coal::Transform3s castTransformInv_;

  /// @brief Helper: compute the AABB of the swept volume using 6 support queries.
  coal::AABB computeSweptAABB() const;
};

}  // namespace tesseract_collision::tesseract_collision_coal
#endif  // TESSERACT_COLLISION_COAL_CASTHULLSHAPE_H
