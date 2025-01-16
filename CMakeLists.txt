cmake_minimum_required(VERSION 3.16...3.27)

project(Simulation)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

# Correctly identify source and header files
file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc) # Changed to lowercase .cc for sources
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh) # Should be .hh for headers, not .CC

add_executable(sim sim.cc ${sources} ${headers})
target_link_libraries(sim ${Geant4_LIBRARIES})

add_custom_target(Simulation DEPENDS sim)

