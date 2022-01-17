cmake_minimum_required(VERSION 3.2)
project(MetaMinimac2 VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions(-DVERSION="${PROJECT_VERSION}" -DUSER="${USER}" -DDATE="${DATE}")

find_package(savvy REQUIRED)
add_executable(MetaMinimac2
        src/Main.cpp
        src/MyVariables.h src/MarkovParameters.h src/simplex.h
        src/MetaMinimac.h src/MetaMinimac.cpp
        src/HaplotypeSet.h src/HaplotypeSet.cpp
        src/MarkovModel.h src/MarkovModel.cpp)
target_link_libraries(MetaMinimac2 savvy)

install(TARGETS MetaMinimac2 RUNTIME DESTINATION bin)
