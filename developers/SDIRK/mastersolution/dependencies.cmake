#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/sdirk_main.cc
  ${DIR}/sdirk.h
  ${DIR}/sdirk.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
