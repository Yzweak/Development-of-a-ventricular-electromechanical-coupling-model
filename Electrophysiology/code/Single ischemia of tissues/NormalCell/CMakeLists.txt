cmake_minimum_required(VERSION 3.25)
project(NormalCell CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(NormalCell model.cu model.cuh)

set_target_properties(NormalCell PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
