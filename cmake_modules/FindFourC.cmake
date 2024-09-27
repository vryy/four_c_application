IF(NOT FOUR_C_DIR)
    message(FATAL_ERROR "FOUR_C_DIR (containing 4CTargets.cmake) is not set")
ENDIF()

find_file(FOURC_TARGETS_PATH "4CTargets.cmake" PATHS ${FOUR_C_DIR})
if(FOURC_TARGETS_PATH)
    include(${FOURC_TARGETS_PATH})
    include_directories(${FOUR_C_INCLUDE_DIR})

    ## find the 4C library
    FIND_LIBRARY(FOUR_C_LIBRARY NAMES 4C PATHS ${FOUR_C_LIBRARY_DIR})

    ## look for dependencies
    set(FOUR_C_TPL_LIBRARIES "")
    # backtrace
    if(FOUR_C_WITH_BACKTRACE MATCHES ON)
        find_file(BACKTRACE_INFO_PATH "Backtrace.cmake" PATHS ${FOUR_C_DIR})
        include(${BACKTRACE_INFO_PATH})
        include_directories(${Backtrace_INCLUDE_DIRS})
        set(FOUR_C_TPL_LIBRARIES ${FOUR_C_TPL_LIBRARIES} ${Backtrace_LIBRARIES})
    endif()
    # HDF5
    if(FOUR_C_WITH_HDF5 MATCHES ON)
        find_file(HDF5_INFO_PATH "HDF5.cmake" PATHS ${FOUR_C_DIR})
        include(${HDF5_INFO_PATH})
        include_directories(${HDF5_INCLUDE_DIRS})
        set(FOUR_C_TPL_LIBRARIES ${FOUR_C_TPL_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    endif()
    # cln
    if(FOUR_C_WITH_CLN MATCHES ON)
        find_file(CLN_INFO_PATH "CLN.cmake" PATHS ${FOUR_C_DIR})
        include(${CLN_INFO_PATH})
        include_directories(${CLN_INCLUDE_DIRS})
        set(FOUR_C_TPL_LIBRARIES ${FOUR_C_TPL_LIBRARIES} ${CLN_LIBRARIES})
    endif()
    # fftw
    if(FOUR_C_WITH_FFTW MATCHES ON)
        find_file(FFTW_INFO_PATH "FFTW.cmake" PATHS ${FOUR_C_DIR})
        include(${FFTW_INFO_PATH})
        include_directories(${FFTW_INCLUDE_DIRS})
        set(FOUR_C_TPL_LIBRARIES ${FOUR_C_TPL_LIBRARIES} ${FFTW_LIBRARIES})
    endif()
    # qhull
    if(FOUR_C_WITH_QHULL MATCHES ON)
        find_file(QHULL_INFO_PATH "Qhull.cmake" PATHS ${FOUR_C_DIR})
        include(${QHULL_INFO_PATH})
        include_directories(${QHULL_INCLUDE_DIRS})
        set(FOUR_C_TPL_LIBRARIES ${FOUR_C_TPL_LIBRARIES} ${QHULL_LIBRARIES})
    endif()

    SET(FOUR_C_LIBRARIES ${FOUR_C_LIBRARY} ${FOUR_C_TPL_LIBRARIES})
    SET(FOUR_C_FOUND TRUE)

    message(STATUS "4C found at " ${FOUR_C_LIBRARIES})
else()
    message(FATAL_ERROR "${FOUR_C_DIR} does not contain 4CTargets.cmake")
endif()
