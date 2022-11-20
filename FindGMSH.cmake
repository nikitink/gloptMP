find_path(GMSH_INCLUDE NAMES gmsh.h)
find_library(GMSH_LIBRARY NAMES gmsh)
if(GMSH_LIBRARY AND NOT TARGET GMSH::GMSH)
    add_library(GMSH::GMSH INTERFACE IMPORTED)
    set_target_properties(GMSH::GMSH PROPERTIES
        INTERFACE_LINK_LIBRARIES "${GMSH_LIBRARY}"
    )
    if(GMSH_INCLUDE)
        set_target_properties(GMSH::GMSH PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE}"
        )
        file(STRINGS "${GMSH_INCLUDE}/gmsh.h" _gmsh_header_version REGEX "GMSH_API_VERSION.*\"[0-9.]\+\"" LIMIT_COUNT 1)
        string(REGEX MATCH "GMSH_API_VERSION \+\"([0-9.]\+)\"" _ "${_gmsh_header_version}")
        set(GMSH_VERSION ${CMAKE_MATCH_1})
    endif()
endif()

if(NOT GMSH_LIBRARY)
    message(WARNING "FindGMSH: Please set GMSH_ROOT pointing to your install prefix of GMSH.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMSH
    REQUIRED_VARS GMSH_LIBRARY
    VERSION_VAR GMSH_VERSION
)

mark_as_advanced(GMSH_INCLUDE GMSH_LIBRARY)
