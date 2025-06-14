

add_library(${PROJECT_NAME} STATIC composite_radix_dft.cpp)


# set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/lib)


target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} m)



# set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/lib)
install(FILES composite_radix_dft.hpp DESTINATION include)