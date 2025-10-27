cmake_minimum_required(VERSION 4.0)
project(polynomial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

# macOS 特定配置
if(APPLE)
    execute_process(
        COMMAND xcrun --show-sdk-path
        OUTPUT_VARIABLE MACOS_SDK_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isysroot ${MACOS_SDK_PATH}")
endif()

# 设置输出目录
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# 收集所有的模块文件
file(GLOB_RECURSE MODULE_SOURCES "src/*.cppm")

# 创建动态链接库
add_library(polynomial SHARED)

target_sources(polynomial
    PUBLIC
        FILE_SET CXX_MODULES FILES
        ${MODULE_SOURCES}
)

# 设置库的输出目录
set_target_properties(polynomial PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
    CXX_SCAN_FOR_MODULES ON
)

# 添加examples子目录
add_subdirectory(examples)

add_custom_command(
    TARGET polynomial POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_BINARY_DIR}/compile_commands.json"
    "${PROJECT_SOURCE_DIR}/compile_commands.json"
    COMMENT "Copying compile_commands.json to project root"
)
