set(SRC_LIST
	source.cc
)

add_executable(simu ${SRC_LIST})
target_link_libraries(simu ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libplinkio.a ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libcsv.a ${Boost_LIBRARIES})
