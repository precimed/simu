include_directories( ${PROJECT_SOURCE_DIR}/libplinkio/src )

set(SRC_LIST
	source.cc
)

add_executable(simu ${SRC_LIST})
target_link_libraries(simu libplinkio libcsv ${Boost_LIBRARIES})

install( TARGETS simu DESTINATION bin )
