## AllTests is dependent on these files being available (in the event of a
## VPATH build, they should be copied to the build directory).
#check_DATA = TestDDaceUserInputSamplerData

set(AllTests_srcs
  AllTests.cpp 
  arrcmp.cpp 
  suite.cpp 
  TestDDaceArraySampler.cpp 
  TestDDaceBoxBehnkenSampler.cpp 
  TestDDaceCentralCompositeSampler.cpp 
  TestDDaceFactorialSampler.cpp 
  TestDDaceLHSampler.cpp 
  TestDDaceOALHSampler.cpp 
  TestDDaceOASampler.cpp 
  TestDDaceRandomSampler.cpp 
  TestDDaceSamplePoint.cpp  
  TestDDaceSampler.cpp 
  TestDDaceUserInputSampler.cpp 
  TestDistribution.cpp 
  TestFactor.cpp 
  TestMainEffectsAnalyzer.cpp 
  TestMainEffectsConverter.cpp 
  TestMainEffectsExcelOutput.cpp
  TestMean.cpp 
  TestNormalDistribution.cpp
  TestPseudoRandom.cpp 
  test.cpp 
  TestResponse.cpp 
  TestStdDeviation.cpp 
  TestUniformDistribution.cpp
)

if(DDACE_ENABLE_STANDALONE)
  list(APPEND AllTests_srcs TestMarsAnalyzer.cpp TestOneWayANOVA.cpp)
endif()

add_executable(AllTests ${AllTests_srcs})
target_link_libraries(AllTests ddace)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/TestDDaceUserInputSamplerData
  ${CMAKE_CURRENT_BINARY_DIR}/TestDDaceUserInputSamplerData
  COPYONLY)
