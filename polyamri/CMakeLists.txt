# polyamri.c always gets touched so that its version information is up to date.
add_custom_target(update_polyamri_c ALL
                  cmake -E touch polyamri.c)

# Library.
add_library(polyamri polyamri.c amr_patch.c amr_patch_set.c amr_grid.c 
            amr_grid_interpolator.c amr_grid_hierarchy.c amr_data_hierarchy.c 
            linear_amr_grid_interpolator.c
            silo_file_write_amr_data_hierarchy.c)
#            interpreter_register_amr_functions.c)
add_dependencies(polyamri update_polyamri_c)

set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyamri DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polyamri)

