# 设置可执行文件输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/bin)

# 收集examples目录下所有的cpp文件
file(GLOB EXAMPLE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# 为每个example创建可执行文件
foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    # 获取文件名（不含扩展名）
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)

    # 创建可执行文件
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})

    # 链接动态库
    target_link_libraries(${EXAMPLE_NAME} PRIVATE polynomial)

    # 设置输出目录和模块扫描
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/bin
        CXX_SCAN_FOR_MODULES ON
    )
endforeach()
