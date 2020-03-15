cmake_minimum_required(VERSION 3.0.0)
project(lab-1 VERSION 0.1.0)

set(src_dir ${PROJECT_SOURCE_DIR}/src)
set(include_dir ${PROJECT_SOURCE_DIR}/include)
set(lib_dir ${PROJECT_SOURCE_DIR}/lib)
set(resource_dir ${PROJECT_SOURCE_DIR}/resource)
set(out_dir ${PROJECT_SOURCE_DIR}/out)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${out_dir})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${out_dir})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${out_dir})

get_filename_component(project_name ${PROJECT_SOURCE_DIR} NAME)

include_directories(${include_dir})
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

file(GLOB_RECURSE sources ${src_dir}/*.c)
add_executable(${project_name} ${sources})

set_property(TARGET ${project_name} PROPERTY C_STANDARD 11)
target_link_libraries(${project_name} gmp tcmalloc_minimal)

set (CMAKE_C_FLAGS_DEBUG      "${CMAKE_C_FLAGS_DEBUG}      -lto -O3")
set (CMAKE_CXX_FLAGS_DEBUG    "${CMAKE_CXX_FLAGS_DEBUG}    -lto -O3")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -lto -O3")

add_custom_command(
    TARGET ${project_name} PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${out_dir}
)

add_custom_command(
    TARGET ${project_name} PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E make_directory ${out_dir}
)

add_custom_command(
    TARGET ${project_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${resource_dir} ${out_dir}
)
