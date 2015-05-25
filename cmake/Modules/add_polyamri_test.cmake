# This function adds a (serial) unit test executable to be built using cmockery.
function(add_polyamri_test exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} cmockery ${POLYAMRI_LIBS})
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
  add_test(${exe} ${exe})
endfunction()

# This function adds a parallel unit test executable to be built using gtest.
# The procs argument is a list of numbers of processes to be run.
# 1 test run will be generated for each processor number value.
function(add_mpi_polyamri_test exe procs)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${POLYAMRI_LIBS})
  foreach (proc ${procs})
    add_test(${exe}_${proc}_proc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${exe} ${MPIEXEC_POSTFLAGS})
  endforeach()
endfunction()

