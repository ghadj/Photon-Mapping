# Top-Level CMakeList.txt

cmake_minimum_required (VERSION 2.6)
project ( Photon_Mapping )

set(PROJECT_LIBRARIES_DIR libs)

set(PROJECT_SOURCE_DIR src)

set(PROJECT_DIR ./)

find_package (SDL)

if ( NOT SDL_FOUND )
   message ( FATAL_ERROR "SDL not found!" )
endif ( NOT SDL_FOUND )

add_executable( Photon_Mapping ${PROJECT_SOURCE_DIR}/photon_map.cpp )

if(SDL_FOUND)
	include_directories(
		${SDL_INCLUDE_DIR}
		${PROJECT_LIBRARIES_DIR}/glm
	    ${PROJECT_SOURCE_DIR}
        ${PROJECT_DIR}
    )
	#link_libraries(${SDL_LIBRARY})
    target_link_libraries(Photon_Mapping ${SDL_LIBRARY})
endif(SDL_FOUND)

