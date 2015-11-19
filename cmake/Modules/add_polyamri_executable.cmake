function(add_polyamri_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${POLYAMRI_LIBRARIES})
endfunction(add_polyamri_executable)

