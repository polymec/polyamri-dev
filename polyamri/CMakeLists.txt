include(add_polymec_library)

# Library.
add_polymec_library(polyamri polyamri.c 
                    str_grid_patch.c str_grid.c str_grid_factory.c 
                    str_grid_assembly.c str_grid_assembly_factory.c
                    str_grid_patch_filler.c str_grid_cell_filler.c str_grid_cell_filler_factory.c
                    str_grid_cell_data.c str_grid_face_data.c str_grid_edge_data.c str_grid_node_data.c
                    str_grid_cell_solver.c hypre_str_grid_cell_solvers.c
                    str_ode_integrator.c ark_str_ode_integrator.c
                    grid_to_bbox_coord_mapping.c
#                    cubed_sphere_coord_mappings.c
                    silo_file.c
                    interpreter_register_polyamri_functions.c)

set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYAMRI_LIBRARIES polyamri;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polyamri DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/polyamri)

