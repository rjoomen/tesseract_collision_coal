/**
 * @file get_collision_objects_transform_unit.cpp
 * @brief Verifies getCollisionObjectsTransform(LinkId) round-trips the setter for the Coal backends.
 *
 * @copyright Copyright (c) 2026, Southwest Research Institute
 *
 * @par License
 * Software License Agreement (Apache License)
 * @par
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 */
#include <tesseract/common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <gtest/gtest.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract/collision/test_suite/get_collision_objects_transform_unit.hpp>
#include <tesseract/collision/coal/coal_discrete_managers.h>
#include <tesseract/collision/coal/coal_cast_managers.h>

using namespace tesseract::collision;
using namespace tesseract::collision::tesseract_collision_coal;

TEST(TesseractCollisionUnit, CoalDiscreteBVHGetCollisionObjectsTransform)  // NOLINT
{
  CoalDiscreteBVHManager checker;
  test_suite::runDiscreteGetCollisionObjectsTransformUnit(checker);
}

TEST(TesseractCollisionUnit, CoalCastBVHGetCollisionObjectsTransform)  // NOLINT
{
  CoalCastBVHManager checker;
  test_suite::runContinuousGetCollisionObjectsTransformUnit(checker);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
