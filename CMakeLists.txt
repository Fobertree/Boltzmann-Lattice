cmake_minimum_required(VERSION 3.27)
project(BoltzmannLattice)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_GENERATOR_PLATFORM "x64")

add_compile_options(-O3)

include_directories(include)
#include("C:/Users/Alexa/.vcpkg-clion/vcpkg/scripts/buildsystems/vcpkg.cmake")

find_package(OpenGL REQUIRED)
# below packages installed via vcpkg
find_package(glfw3 REQUIRED 3.4.0)
find_package(glad REQUIRED)
find_package(glm REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(BoltzmannLattice
        src/main.cpp
        include/math_util.h
        include/utils.h
        include/shader.h
        src/math_util.cpp
)


# link libraries after add_executable
target_link_libraries(BoltzmannLattice PRIVATE OpenGL::GL)
target_link_libraries(BoltzmannLattice PRIVATE glm::glm)
target_link_libraries(BoltzmannLattice PRIVATE glfw)
target_link_libraries(BoltzmannLattice PRIVATE glad::glad)
target_link_libraries(BoltzmannLattice PRIVATE Eigen3::Eigen)