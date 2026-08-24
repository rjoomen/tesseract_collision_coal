#include <tesseract/common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <gtest/gtest.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract/collision/test_suite/collision_clone_unit.hpp>
#include <tesseract/collision/bullet/bullet_discrete_simple_manager.h>
#include <tesseract/collision/bullet/bullet_discrete_bvh_manager.h>
#include <tesseract/collision/bullet/bullet_cast_simple_manager.h>
#include <tesseract/collision/bullet/bullet_cast_bvh_manager.h>
#include <tesseract/collision/fcl/fcl_discrete_managers.h>
#include <tesseract/collision/coal/coal_discrete_managers.h>
#include <tesseract/collision/coal/coal_cast_managers.h>

using namespace tesseract::collision;

TEST(TesseractCollisionUnit, BulletDiscreteSimpleCollisionCloneUnit)  // NOLINT
{
  tesseract::collision::BulletDiscreteSimpleManager checker;
  test_suite::runTest(checker, 0.001, 0.001, 0.001);
}

TEST(TesseractCollisionUnit, BulletDiscreteBVHCollisionCloneUnit)  // NOLINT
{
  tesseract::collision::BulletDiscreteBVHManager checker;
  test_suite::runTest(checker, 0.001, 0.001, 0.001);
}

TEST(TesseractCollisionUnit, FCLDiscreteBVHCollisionCloneUnit)  // NOLINT
{
  tesseract::collision::FCLDiscreteBVHManager checker;
  test_suite::runTest(checker, 0.001, 0.001, 0.001);
}

TEST(TesseractCollisionUnit, CoalDiscreteBVHCollisionCloneUnit)  // NOLINT
{
  tesseract_collision_coal::CoalDiscreteBVHManager checker;
  test_suite::runTest(checker, 0.001, 0.001, 0.001);
}

TEST(TesseractCollisionUnit, BulletDiscreteSimpleCollisionCloneOrderUnit)  // NOLINT
{
  tesseract::collision::BulletDiscreteSimpleManager checker;
  test_suite::runCloneOrderTest(checker);
}

TEST(TesseractCollisionUnit, BulletDiscreteBVHCollisionCloneOrderUnit)  // NOLINT
{
  tesseract::collision::BulletDiscreteBVHManager checker;
  test_suite::runCloneOrderTest(checker);
}

TEST(TesseractCollisionUnit, FCLDiscreteBVHCollisionCloneOrderUnit)  // NOLINT
{
  tesseract::collision::FCLDiscreteBVHManager checker;
  test_suite::runCloneOrderTest(checker);
}

TEST(TesseractCollisionUnit, CoalDiscreteBVHCollisionCloneOrderUnit)  // NOLINT
{
  tesseract_collision_coal::CoalDiscreteBVHManager checker;
  test_suite::runCloneOrderTest(checker);
}

TEST(TesseractCollisionUnit, BulletCastSimpleCollisionCloneOrderUnit)  // NOLINT
{
  tesseract::collision::BulletCastSimpleManager checker;
  test_suite::runCloneOrderTest(checker);
}

TEST(TesseractCollisionUnit, BulletCastBVHCollisionCloneOrderUnit)  // NOLINT
{
  tesseract::collision::BulletCastBVHManager checker;
  test_suite::runCloneOrderTest(checker);
}

TEST(TesseractCollisionUnit, CoalCastBVHCollisionCloneOrderUnit)  // NOLINT
{
  tesseract_collision_coal::CoalCastBVHManager checker;
  test_suite::runCloneOrderTest(checker);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
