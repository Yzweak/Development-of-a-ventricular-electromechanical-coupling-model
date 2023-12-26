cmake_minimum_required(VERSION 3.25)
project(MixFiber CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(MixFiber model.cu model.cuh)

set_target_properties(MixFiber PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
