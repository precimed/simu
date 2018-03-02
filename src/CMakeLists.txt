include_directories( ${PROJECT_SOURCE_DIR}/libplinkio/src )

set(SRC_LIST
	source.cc
)

add_executable(simu ${SRC_LIST})
set_target_properties(simu PROPERTIES LINK_FLAGS "-static" )
target_link_libraries(simu ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libplinkio.a ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libcsv.a ${Boost_LIBRARIES})

install( TARGETS simu DESTINATION bin )
