cmake_minimum_required(VERSION 3.16)
project(SEE_in_vacuum)

# Check for conda environment and add to CMAKE_PREFIX_PATH
if(DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{CONDA_PREFIX})
    message(STATUS "Found conda environment: $ENV{CONDA_PREFIX}")
endif()

# Find Geant4 package
find_package(Geant4 REQUIRED)

# Include Geant4 directories
include(${Geant4_USE_FILE})

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.cc
    src/DetectorConstruction.cc
    src/PrimaryGeneratorAction.cc
    src/PhysicsList.cc
    src/ActionInitialization.cc
    src/RunAction.cc
    src/SteppingAction.cc
)

# Header files
set(HEADERS
    include/DetectorConstruction.hh
    include/PrimaryGeneratorAction.hh
    include/PhysicsList.hh
    include/ActionInitialization.hh
    include/RunAction.hh
    include/SteppingAction.hh
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Link Geant4 libraries
target_link_libraries(${PROJECT_NAME} ${Geant4_LIBRARIES})
