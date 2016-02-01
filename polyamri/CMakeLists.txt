# Library.
add_library(polyamri polyamri.c 
            str_grid_patch.c str_grid_patch_filler.c str_grid.c
            str_grid_cell_data.c
            amr_patch.c amr_grid.c amr_grid_data.c
            amr_grid_interpolator.c amr_grid_hierarchy.c amr_data_hierarchy.c 
            amr_grid_assembly.c grid_to_bbox_coord_mapping.c
            partition_amr_grid.c
            linear_amr_grid_interpolator.c
            silo_file_amr_methods.c)
#            interpreter_register_amr_functions.c)

set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyamri DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/polyamri)

