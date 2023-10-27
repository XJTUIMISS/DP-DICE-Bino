cmake_minimum_required(VERSION 3.15)

# set the project name
project(FMS VERSION 1.0.1)

# add the executable
set(SRC_LIST core.cpp FMSsketch.cpp MurmurHash3/MurmurHash3.cpp FMsketch/FMsketch.cpp batch_est.cpp LLSketch/LLSketch.cpp noise/noise.cpp)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(FMS pthread)

configure_file(FMSConfig.h.in FMSConfig.h)
string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)

target_include_directories(${PROJECT_NAME} PUBLIC
                           ${PROJECT_BINARY_DIR}
                           )