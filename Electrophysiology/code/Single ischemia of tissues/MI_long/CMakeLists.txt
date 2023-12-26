cmake_minimum_required(VERSION 3.25)
project(MI_long CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(MI_long model.cu model.cuh)

set_target_properties(MI_long PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
