cmake_minimum_required(VERSION 3.25)
project(Ischemia_1b CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(Ischemia_1b model.cu model.cuh)

set_target_properties(Ischemia_1b PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
