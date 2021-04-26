add_library(wake SHARED wake.cpp ring.cpp)
target_compile_features(wake PUBLIC cxx_constexpr)
