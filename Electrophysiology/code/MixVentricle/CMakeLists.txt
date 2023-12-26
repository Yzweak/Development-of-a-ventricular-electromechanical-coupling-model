cmake_minimum_required(VERSION 3.25)
project(MixVentricle CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(MixVentricle model.cuh model.cu)

set_target_properties(MixVentricle PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
