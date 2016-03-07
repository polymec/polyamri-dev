function(add_polyamri_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${POLYAMRI_LIBRARIES})
  endif()
endfunction(add_polyamri_library)

