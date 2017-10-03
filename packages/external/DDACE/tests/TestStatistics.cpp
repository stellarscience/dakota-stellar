#include "TestStatistics.h"

TestStatistics::TestStatistics()
{
}

TestStatistics::~TestStatistics()
{
}

void TestStatistics::run()
{
  testAll();
}

void TestStatistics::testAll()
{
  std::vector<double> data;
  data.push_back(20.48);
  data.push_back(20.32);
  data.push_back(20.57);
  data.push_back(19.73);

  _test( closeEnough(Statistics::sum(data), 81.1, 0.0001));

  _test( closeEnough(Statistics::average(data), 20.275, 0.0001));

  _test( closeEnough(Statistics::sumOfSquares(data,20.275), 0.4281, 0.0001));

  std::cout << Statistics::sumOfSquares(data, 20.275) << std::endl;
  _test( closeEnough(Statistics::variance(data), 0.1427, 0.0001));

}
