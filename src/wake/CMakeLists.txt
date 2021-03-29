add_library(wake SHARED wake.cpp)
target_compile_features(wake PUBLIC cxx_constexpr)
