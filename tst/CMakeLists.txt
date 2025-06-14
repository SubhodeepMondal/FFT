add_executable(${PROJECT_NAME}_test
    main.cpp
)
target_link_libraries(${PROJECT_NAME}_test
    ${CMAKE_BINARY_DIR}/lib ${PROJECT_NAME}
)
target_include_directories(${PROJECT_NAME}_test
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src
)