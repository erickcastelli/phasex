cmake_minimum_required (VERSION 2.6)

set (CMAKE_INSTALL_PREFIX "/usr/local/bin")
project (phasex)

find_package( Boost REQUIRED COMPONENTS program_options regex )
include_directories( ${Boost_INCLUDE_DIRS} )

set (CMAKE_CXX_STANDARD 11)

add_executable(phasex recode.cpp recode.hpp vcf.cpp vcf.hpp external.h functions.cpp functions.hpp main.cpp phase-ps.cpp phase-ps.hpp ThreadPool.h)
target_link_libraries(phasex pthread)
