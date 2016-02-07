# Library.
add_library(polyamri polyamri.c 
            str_grid_patch.c str_grid_patch_filler.c str_grid.c str_grid_assembly.c
            str_grid_cell_data.c str_grid_face_data.c str_grid_edge_data.c str_grid_node_data.c
            str_ode_integrator.c
            grid_to_bbox_coord_mapping.c
            silo_file_str_methods.c
            interpreter_register_polyamri_functions.c)

set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyamri DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/polyamri)

