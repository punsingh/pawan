add_library(io SHARED io.cpp)
target_compile_features(io PUBLIC cxx_constexpr)
