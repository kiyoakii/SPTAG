# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.

set(AnnService ${PROJECT_SOURCE_DIR}/AnnService)

include_directories(${AnnService})

file(GLOB_RECURSE HDR_FILES ${AnnService}/inc/Core/*.h  ${AnnService}/inc/Helper/*.h)
file(GLOB_RECURSE SRC_FILES ${AnnService}/src/Core/*.cpp ${AnnService}/src/Helper/*.cpp)

list(REMOVE_ITEM HDR_FILES
    ${AnnService}/inc/Core/Common/DistanceUtils.h 
    ${AnnService}/inc/Core/Common/InstructionUtils.h
    ${AnnService}/inc/Core/Common/CommonUtils.h
    )

list(REMOVE_ITEM SRC_FILES
    ${AnnService}/src/Core/Common/DistanceUtils.cpp
    ${AnnService}/src/Core/Common/InstructionUtils.cpp
    )

add_library (DistanceUtils STATIC 
    inc/Core/Common/DistanceUtils.h 
    inc/Core/Common/InstructionUtils.h
    inc/Core/Common/CommonUtils.h
    inc/Core/Common/IQuantizer.h
    src/Core/Common/DistanceUtils.cpp
    src/Core/Common/InstructionUtils.cpp
    )

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(DistanceUtils PRIVATE -mavx2 -mavx -msse -msse2 -fPIC)
endif()

find_package(RocksDB CONFIG)
if((DEFINED RocksDB_DIR) AND RocksDB_DIR)
    list(APPEND RocksDB_LIBRARIES RocksDB::rocksdb)
    find_package(uring)
    message (STATUS "Found RocksDB ${RocksDB_VERSION}")
    message (STATUS "RocksDB: ${RocksDB_DIR}")
else()
    message (FATAL_ERROR "Could not find RocksDB!")
endif()

add_library (SPTAGLib SHARED ${SRC_FILES} ${HDR_FILES})
target_link_libraries (SPTAGLib DistanceUtils ${RocksDB_LIBRARIES} ${uring_LIBRARIES})
add_library (SPTAGLibStatic STATIC ${SRC_FILES} ${HDR_FILES})
target_link_libraries (SPTAGLibStatic DistanceUtils ${RocksDB_LIBRARIES} ${uring_LIBRARIES})
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(SPTAGLibStatic PRIVATE -fPIC)
endif()

install(TARGETS SPTAGLib SPTAGLibStatic
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

if (NOT LIBRARYONLY)
    file(GLOB SERVER_HDR_FILES ${AnnService}/inc/Server/*.h ${AnnService}/inc/Socket/*.h)
    file(GLOB SERVER_FILES ${AnnService}/src/Server/*.cpp ${AnnService}/src/Socket/*.cpp)
    add_executable (server ${SERVER_FILES} ${SERVER_HDR_FILES})
    target_link_libraries(server ${Boost_LIBRARIES} SPTAGLibStatic)

    file(GLOB CLIENT_HDR_FILES ${AnnService}/inc/Client/*.h ${AnnService}/inc/Socket/*.h)
    file(GLOB CLIENT_FILES ${AnnService}/src/Client/*.cpp ${AnnService}/src/Socket/*.cpp)
    add_executable (client ${CLIENT_FILES} ${CLIENT_HDR_FILES})
    target_link_libraries(client ${Boost_LIBRARIES} SPTAGLibStatic)

    file(GLOB AGG_HDR_FILES ${AnnService}/inc/Aggregator/*.h ${AnnService}/inc/Socket/*.h ${AnnService}/inc/Server/QueryParser.h)
    file(GLOB AGG_FILES ${AnnService}/src/Aggregator/*.cpp ${AnnService}/src/Socket/*.cpp ${AnnService}/src/Server/QueryParser.cpp)
    add_executable (aggregator ${AGG_FILES} ${AGG_HDR_FILES})
    target_link_libraries(aggregator ${Boost_LIBRARIES} SPTAGLibStatic)

    file(GLOB BUILDER_FILES ${AnnService}/src/IndexBuilder/*.cpp)
    add_executable (indexbuilder ${BUILDER_FILES})
    target_link_libraries(indexbuilder ${Boost_LIBRARIES} SPTAGLibStatic)

    file(GLOB SEARCHER_FILES ${AnnService}/src/IndexSearcher/*.cpp)
    add_executable (indexsearcher ${SEARCHER_FILES})
    target_link_libraries(indexsearcher ${Boost_LIBRARIES} SPTAGLibStatic)
    
    install(TARGETS server client aggregator indexbuilder indexsearcher
      RUNTIME DESTINATION bin
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib)
endif()

set(VECTORSEARCH_INC_DIR ${AnnService}/inc/SSDServing/VectorSearch)
set(VECTORSEARCH_IMP_DIR ${AnnService}/src/SSDServing/VectorSearch)
file(GLOB_RECURSE SSD_SERVING_HDR_FILES ${AnnService}/inc/SSDServing/*.h)
file(GLOB_RECURSE SSD_SERVING_FILES ${AnnService}/src/SSDServing/*.cpp)
if(NOT WIN32)
    list(REMOVE_ITEM SSD_SERVING_HDR_FILES 
            ${VECTORSEARCH_INC_DIR}/AsyncFileReader.h
            )
elseif(WIN32)
    list(REMOVE_ITEM SSD_SERVING_HDR_FILES 
            ${VECTORSEARCH_INC_DIR}/AsyncFileReaderLinux.h
            )
endif()

add_executable(ssdserving ${SSD_SERVING_HDR_FILES} ${SSD_SERVING_FILES})
target_link_libraries(ssdserving SPTAGLibStatic ${Boost_LIBRARIES} ${RocksDB_LIBRARIES})
target_compile_definitions(ssdserving PRIVATE _exe)

# for Test
add_library(ssdservingLib ${SSD_SERVING_HDR_FILES} ${SSD_SERVING_FILES})
target_link_libraries(ssdservingLib SPTAGLibStatic ${Boost_LIBRARIES})

find_package(MPI)
if (MPI_FOUND)
    message (STATUS "Found MPI.")
    message (STATUS "MPI Include Path: ${MPI_CXX_INCLUDE_PATH}")
    message (STATUS "MPI Libraries: ${MPI_CXX_LIBRARIES}")
    include_directories (${MPI_CXX_INCLUDE_PATH})
    file(GLOB PARTITION_FILES ${AnnService}/src/BalancedDataPartition/*.cpp)
    add_executable(balanceddatapartition ${PARTITION_FILES})
    target_link_libraries(balanceddatapartition SPTAGLibStatic ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()

install(TARGETS ssdservingLib ssdserving
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
