set(CMAKE_REQUIRED_LIBRARIES ${SPDLOG_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${SPDLOG_INCLUDE_DIRS})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include \"spdlog/spdlog.h\"

  int main()
  {
      spdlog::info(\"Welcome to spdlog!\");
      spdlog::error(\"Some error message with arg: {}\", 1);

      spdlog::set_level(spdlog::level::debug); // Set global log level to debug
      spdlog::debug(\"This message should be displayed..\");

      // change log pattern
      spdlog::set_pattern(\"[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v\");
  }
  "
  SPDLOG_COMPILES)

if (NOT SPDLOG_COMPILES)
  set(SPDLOG_TEST_FAIL 1)
else ()
  set(SPDLOG_TEST_FAIL 0)
endif()
