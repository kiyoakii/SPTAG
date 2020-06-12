# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.

file(GLOB HDR_FILES ${PROJECT_SOURCE_DIR}/AnnService/inc/Core/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Core/Common/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Core/BKT/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Core/KDT/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Helper/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Helper/VectorSetReaders/*.h)
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/AnnService/src/Core/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Core/Common/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Core/BKT/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Core/KDT/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Helper/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Helper/VectorSetReaders/*.cpp)

include_directories(${PROJECT_SOURCE_DIR}/AnnService)

add_library (SPTAGLib SHARED ${SRC_FILES} ${HDR_FILES})
target_link_libraries (SPTAGLib)
add_library (SPTAGLibStatic STATIC ${SRC_FILES} ${HDR_FILES})
# set_target_properties(SPTAGLibStatic PROPERTIES OUTPUT_NAME SPTAGLib)

file(GLOB SERVER_HDR_FILES ${HDR_FILES} ${PROJECT_SOURCE_DIR}/AnnService/inc/Server/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Socket/*.h)
file(GLOB SERVER_FILES ${SRC_FILES} ${PROJECT_SOURCE_DIR}/AnnService/src/Server/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Socket/*.cpp)
add_executable (server ${SERVER_FILES} ${SERVER_HDR_FILES})
target_link_libraries(server ${Boost_LIBRARIES})

file(GLOB CLIENT_HDR_FILES ${HDR_FILES} ${PROJECT_SOURCE_DIR}/AnnService/inc/Client/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Socket/*.h)
file(GLOB CLIENT_FILES ${SRC_FILES} ${PROJECT_SOURCE_DIR}/AnnService/src/Client/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Socket/*.cpp)
add_executable (client ${CLIENT_FILES} ${CLIENT_HDR_FILES})
target_link_libraries(client ${Boost_LIBRARIES})

file(GLOB AGG_HDR_FILES ${HDR_FILES} ${PROJECT_SOURCE_DIR}/AnnService/inc/Aggregator/*.h ${PROJECT_SOURCE_DIR}/AnnService/inc/Socket/*.h)
file(GLOB AGG_FILES ${SRC_FILES} ${PROJECT_SOURCE_DIR}/AnnService/src/Aggregator/*.cpp ${PROJECT_SOURCE_DIR}/AnnService/src/Socket/*.cpp)
add_executable (aggregator ${AGG_FILES} ${AGG_HDR_FILES})
target_link_libraries(aggregator ${Boost_LIBRARIES})

file(GLOB BUILDER_HDR_FILES ${HDR_FILES} ${PROJECT_SOURCE_DIR}/AnnService/inc/IndexBuilder/*.h)
file(GLOB BUILDER_FILES ${SRC_FILES} ${PROJECT_SOURCE_DIR}/AnnService/src/IndexBuilder/*.cpp)
add_executable (indexbuilder ${BUILDER_FILES} ${BUILDER_HDR_FILES})
target_link_libraries(indexbuilder ${Boost_LIBRARIES})

file(GLOB SEARCHER_FILES ${SRC_FILES} ${PROJECT_SOURCE_DIR}/AnnService/src/IndexSearcher/*.cpp)
add_executable (indexsearcher ${SEARCHER_FILES} ${HDR_FILES})
target_link_libraries(indexsearcher ${Boost_LIBRARIES})

set(SSDSERVING_INC_DIR ${PROJECT_SOURCE_DIR}/AnnService/inc/SSDServing)
set(SSDSERVING_IMP_DIR ${PROJECT_SOURCE_DIR}/AnnService/src/SSDServing)
set(VECTORSEARCH_INC_DIR ${SSDSERVING_INC_DIR}/VectorSearch)
set(VECTORSEARCH_IMP_DIR ${SSDSERVING_IMP_DIR}/VectorSearch)
file(GLOB_RECURSE SSD_SERVING_HDR_FILES ${SSDSERVING_INC_DIR}/*.h)
file(GLOB_RECURSE SSD_SERVING_FILES ${SSDSERVING_IMP_DIR}/*.cpp)
if(NOT WIN32)
    list(REMOVE_ITEM SSD_SERVING_HDR_FILES 
            # prototype
            ${VECTORSEARCH_INC_DIR}/DiskAccessUtils.h
            ${VECTORSEARCH_INC_DIR}/DiskFileReader.h
            ${VECTORSEARCH_INC_DIR}/DiskListCommonUtils.h
            ${VECTORSEARCH_INC_DIR}/ExtraFullGraphSearcher.h
            ${VECTORSEARCH_INC_DIR}/IExtraSearcher.h
            ${VECTORSEARCH_INC_DIR}/PrioritizedDiskFileReader.h
            )
    list(REMOVE_ITEM SSD_SERVING_FILES
            # Implementation
            ${VECTORSEARCH_IMP_DIR}/DiskAccessUtils.cpp
            ${VECTORSEARCH_IMP_DIR}/DiskFileReader.cpp
            ${VECTORSEARCH_IMP_DIR}/PrioritizedDiskFileReader.cpp
            )
elseif(WIN32)
    list(REMOVE_ITEM SSD_SERVING_HDR_FILES 
            # prototype
            ${VECTORSEARCH_INC_DIR}/ExtraFullGraphSearcherLinux.h
            ${VECTORSEARCH_INC_DIR}/IExtraSearcherLinux.h
            )
endif()

add_executable(ssdserving ${SSD_SERVING_HDR_FILES} ${SSD_SERVING_FILES})
target_link_libraries(ssdserving SPTAGLibStatic ${Boost_LIBRARIES})
target_compile_definitions(ssdserving PRIVATE _exe)

# for Test
add_library(ssdservingLib ${SSD_SERVING_HDR_FILES} ${SSD_SERVING_FILES})
target_link_libraries(ssdservingLib SPTAGLibStatic ${Boost_LIBRARIES})

install(TARGETS SPTAGLib SPTAGLibStatic server client aggregator indexbuilder indexsearcher ssdserving
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
