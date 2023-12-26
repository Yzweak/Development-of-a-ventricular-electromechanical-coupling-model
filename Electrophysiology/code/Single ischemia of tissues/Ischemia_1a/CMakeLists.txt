cmake_minimum_required(VERSION 3.25)
project(Ischemia_1a CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(Ischemia_1a model.cu model.cuh)

set_target_properties(Ischemia_1a PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
