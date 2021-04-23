add_library(integration SHARED integration.cpp)
target_compile_features(integration PUBLIC cxx_constexpr)
