# polyamri.c always gets touched so that its version information is up to date.
add_custom_target(update_polyamri_c ALL
                  cmake -E touch polyamri.c)

# Library.
add_library(polyamri polyamri.c)# tile.c tile_set.c grid_level.c 
#            grid_interpolator.c grid_hierarchy.c data_hierarchy.c 
#            linear_grid_interpolator.c write_silo_data_hierarchy.c
#            interpreter_register_amr_functions.c)
add_dependencies(polyamri update_polyamri_c)

set(POLYAMRI_LIBS polyamri;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYAMRI_LIBS polyamri;${POLYMEC_LIBS})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyamri DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polyamri)

