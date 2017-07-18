include_directories(${PROJECT_PATH}/extern/CompactNSearch/include)

add_executable(FDMSPH
    main.cpp
    json.hpp
    sphsolver.hpp
    sphsolver.cpp
    particleAttrib.hpp
    particleAttrib.cpp
    helpers.hpp
    filewriter.hpp
    poissonsample.hpp
)

target_link_libraries(FDMSPH CompactNSearch)

find_package(H5Part REQUIRED)
include_directories ( "${H5PART_INCLUDE_PATH}" )
target_link_libraries(FDMSPH "${H5PART_LIBRARIES}")

find_package(hdf5 REQUIRED)
include_directories ( "${HDF5_INCLUDE_PATH}" )
target_link_libraries(FDMSPH "${HDF5_LIBRARIES}")
