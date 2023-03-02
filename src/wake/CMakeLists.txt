add_library(wake SHARED wake.cpp ring.cpp square.cpp vring.cpp rotor_ll.cpp)
target_compile_features(wake PUBLIC cxx_constexpr)
