#include "gtest/gtest.h"

#include <Eigen/Core>
#include "MUQ/Utilities/VectorSlice.h"

using namespace std;
using namespace muq::Utilities;

TEST(Utilities, VectorSlice_Eigen)
{
  Eigen::VectorXd v(7);
  for(int i=0; i<v.size(); ++i)
    v(i) = i;

  auto vSlice = GetSlice(v,0,7,2);
  EXPECT_EQ(4,vSlice.size());
  EXPECT_EQ(0,vSlice(0));
  EXPECT_EQ(0,vSlice[0]);
  EXPECT_EQ(2,vSlice(1));
  EXPECT_EQ(2,vSlice[1]);
  EXPECT_EQ(4,vSlice(2));
  EXPECT_EQ(4,vSlice[2]);
  EXPECT_EQ(6,vSlice(3));
  EXPECT_EQ(6,vSlice[3]);

  auto vSlice2 = GetSlice(vSlice,3,-1, -1);
  EXPECT_EQ(4,vSlice2.size());
  EXPECT_EQ(6,vSlice2(0));
  EXPECT_EQ(6,vSlice2[0]);
  EXPECT_EQ(4,vSlice2(1));
  EXPECT_EQ(4,vSlice2[1]);
  EXPECT_EQ(2,vSlice2(2));
  EXPECT_EQ(2,vSlice2[2]);
  EXPECT_EQ(0,vSlice2(3));
  EXPECT_EQ(0,vSlice2[3]);

  vSlice = GetSlice(v,2,-1,-1);
  EXPECT_EQ(3,vSlice.size());
  EXPECT_EQ(2,vSlice(0));
  EXPECT_EQ(2,vSlice[0]);
  EXPECT_EQ(1,vSlice(1));
  EXPECT_EQ(1,vSlice[1]);
  EXPECT_EQ(0,vSlice(2));
  EXPECT_EQ(0,vSlice[2]);


}



TEST(Utilities, VectorSlice_StdVec)
{
  std::vector<double> v(7);
  for(int i=0; i<v.size(); ++i)
    v[i] = i;

  auto vSlice = GetSlice(v,0,7,2);
  EXPECT_EQ(4,vSlice.size());
  EXPECT_EQ(0,vSlice(0));
  EXPECT_EQ(0,vSlice[0]);
  EXPECT_EQ(2,vSlice(1));
  EXPECT_EQ(2,vSlice[1]);
  EXPECT_EQ(4,vSlice(2));
  EXPECT_EQ(4,vSlice[2]);
  EXPECT_EQ(6,vSlice(3));
  EXPECT_EQ(6,vSlice[3]);

  auto vSlice2 = GetSlice(vSlice, 3,-1, -1);
  EXPECT_EQ(4,vSlice2.size());
  EXPECT_EQ(6,vSlice2(0));
  EXPECT_EQ(6,vSlice2[0]);
  EXPECT_EQ(4,vSlice2(1));
  EXPECT_EQ(4,vSlice2[1]);
  EXPECT_EQ(2,vSlice2(2));
  EXPECT_EQ(2,vSlice2[2]);
  EXPECT_EQ(0,vSlice2(3));
  EXPECT_EQ(0,vSlice2[3]);

  vSlice = GetSlice(v,2,-1,-1);
  EXPECT_EQ(3,vSlice.size());
  EXPECT_EQ(2,vSlice(0));
  EXPECT_EQ(2,vSlice[0]);
  EXPECT_EQ(1,vSlice(1));
  EXPECT_EQ(1,vSlice[1]);
  EXPECT_EQ(0,vSlice(2));
  EXPECT_EQ(0,vSlice[2]);
}
