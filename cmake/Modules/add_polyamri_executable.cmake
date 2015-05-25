# This function adds a (serial) unit test executable to be built using cmockery.
function(add_polyamri_executable exe)
  add_executable(${exe} ${ARGN})
  set_target_properties(${exe} PROPERTIES OUTPUT_NAME ${exe} LINKER_LANGUAGE CXX)
  target_link_libraries(${exe} ${POLYAMRI_LIBS})
endfunction(add_polyamri_executable)

