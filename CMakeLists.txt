project(Heisenberg)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

set(SOURCES
    main.cpp)
set(HEADERS
    tdvp.h
    basisextension.h
    model.h
    interface.h
)


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC include)

target_link_libraries(${PROJECT_NAME} ${GTEST_LIBRARIES} -Wall -pedantic -m64 -std=c++17 -fconcepts -fopenmp -fPIC -I'/home/maciej/itensor' -O2 -DNDEBUG -Wall -Wno-unknown-pragmas -Wno-unused-variable -L'/home/maciej/itensor/lib' -litensor -larmadillo -lpthread -L/usr/lib -lblas -llapack)
