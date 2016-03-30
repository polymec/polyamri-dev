// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <dlfcn.h>
#include "polyamri/hypre_str_grid_cell_solvers.h"

// HYPRE types and definitions.
#define HYPRE_PARCSR  5555
typedef double HYPRE_Real;
typedef double _Complex HYPRE_Complex;
typedef long long HYPRE_Int;
typedef void* HYPRE_StructGrid;
typedef void* HYPRE_StructSolver;
typedef void* HYPRE_StructMatrix;
typedef void* HYPRE_StructStencil;
typedef void* HYPRE_StructVector;
typedef void* HYPRE_ParCSRMatrix;
typedef void* HYPRE_ParVector;
//typedef HYPRE_Int (*HYPRE_PtrToSolverFcn)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector);
//typedef HYPRE_Int (*HYPRE_PtrToParSolverFcn)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
#if POLYMEC_HAVE_MPI
typedef MPI_Comm HYPRE_MPI_Comm;
#else
typedef HYPRE_Int HYPRE_MPI_Comm;
#endif

//------------------------------------------------------------------------
//                    Dynamically-loaded HYPRE library.
//------------------------------------------------------------------------

// Here's a table of function pointers for the HYPRE library.
typedef struct
{
  void* lib;

  HYPRE_Int (*HYPRE_GetError)();
  HYPRE_Int (*HYPRE_ClearAllErrors)();

  HYPRE_Int (*HYPRE_StructGridCreate)(HYPRE_MPI_Comm, int, HYPRE_StructGrid*);
  HYPRE_Int (*HYPRE_StructGridDestroy)(HYPRE_StructGrid);
  HYPRE_Int (*HYPRE_StructGridSetExtents)(HYPRE_StructGrid, int*, int*);
  HYPRE_Int (*HYPRE_StructGridAssemble)(HYPRE_StructGrid);
  HYPRE_Int (*HYPRE_StructGridSetPeriodic)(HYPRE_StructGrid, int*);
  HYPRE_Int (*HYPRE_StructGridSetNumGhost)(HYPRE_StructGrid, int*);

  HYPRE_Int (*HYPRE_StructMatrixCreate)(MPI_Comm, HYPRE_StructGrid, HYPRE_StructStencil, HYPRE_StructMatrix*);
  HYPRE_Int (*HYPRE_StructMatrixDestroy)(HYPRE_StructMatrix);
  HYPRE_Int (*HYPRE_StructMatrixInitialize)(HYPRE_StructMatrix);
  HYPRE_Int (*HYPRE_StructMatrixSetValues)(HYPRE_StructMatrix, int*, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixAddToValues)(HYPRE_StructMatrix, int*, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixSetConstantValues)(HYPRE_StructMatrix, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixAddToConstantValues)(HYPRE_StructMatrix, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixSetBoxValues)(HYPRE_StructMatrix, int*, int*, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixAddToBoxValues)(HYPRE_StructMatrix, int*, int*, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixAssemble)(HYPRE_StructMatrix);
  HYPRE_Int (*HYPRE_StructMatrixSetSymmetric)(HYPRE_StructMatrix, int);
  HYPRE_Int (*HYPRE_StructMatrixSetConstantEntries)(HYPRE_StructMatrix, int, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructMatrixSetNumGhost)(HYPRE_StructMatrix, int*);
  HYPRE_Int (*HYPRE_StructMatrixMatVec)(HYPRE_StructMatrix, HYPRE_StructVector, HYPRE_Real beta, HYPRE_StructVector);

  HYPRE_Int (*HYPRE_StructStencilCreate)(int, int, HYPRE_StructStencil*);
  HYPRE_Int (*HYPRE_StructStencilDestroy)(HYPRE_StructStencil);
  HYPRE_Int (*HYPRE_StructStencilSetElement)(HYPRE_StructStencil, int, int*);

  HYPRE_Int (*HYPRE_StructVectorCreate)(MPI_Comm, HYPRE_StructGrid, HYPRE_StructVector*);
  HYPRE_Int (*HYPRE_StructVectorDestroy)(HYPRE_StructVector);
  HYPRE_Int (*HYPRE_StructVectorInitialize)(HYPRE_StructVector);
  HYPRE_Int (*HYPRE_StructVectorSetValues)(HYPRE_StructVector, int*, HYPRE_Real);
  HYPRE_Int (*HYPRE_StructVectorAddToValues)(HYPRE_StructVector, int*, HYPRE_Real);
  HYPRE_Int (*HYPRE_StructVectorSetBoxValues)(HYPRE_StructVector, int*, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructVectorAddToBoxValues)(HYPRE_StructVector, int*, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructVectorAssemble)(HYPRE_StructVector);
  HYPRE_Int (*HYPRE_StructVectorGetValues)(HYPRE_StructVector, int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_StructVectorGetBoxValues)(HYPRE_StructVector, int*, int*, HYPRE_Real*);

  HYPRE_Int (*HYPRE_StructSMGCreate)(HYPRE_MPI_Comm,  HYPRE_StructSolver*);
  HYPRE_Int (*HYPRE_StructSMGDestroy)(HYPRE_StructSolver);
  HYPRE_Int (*HYPRE_StructSMGSetup)(HYPRE_StructSolver, HYPRE_StructMatrix, HYPRE_StructVector, HYPRE_StructVector);
  HYPRE_Int (*HYPRE_StructSMGSolve)(HYPRE_StructSolver, HYPRE_StructMatrix, HYPRE_StructVector, HYPRE_StructVector);
  HYPRE_Int (*HYPRE_StructSMGSetTol)(HYPRE_StructSolver, HYPRE_Real);
  HYPRE_Int (*HYPRE_StructSMGSetMaxIter)(HYPRE_StructSolver, int);
  HYPRE_Int (*HYPRE_StructSMGSetRelChange)(HYPRE_StructSolver, int);
  HYPRE_Int (*HYPRE_StructSMGSetZeroGuess)(HYPRE_StructSolver);
  HYPRE_Int (*HYPRE_StructSMGSetNonZeroGuess)(HYPRE_StructSolver);
  HYPRE_Int (*HYPRE_StructSMGSetNumPreRelax)(HYPRE_StructSolver, int);
  HYPRE_Int (*HYPRE_StructSMGSetNumPostRelax)(HYPRE_StructSolver, int);
  HYPRE_Int (*HYPRE_StructSMGGetNumIterations)(HYPRE_StructSolver, int*);
  HYPRE_Int (*HYPRE_StructSMGGetFinalRelativeResidualNorm)(HYPRE_StructSolver, HYPRE_Real*);
} hypre_lib_t;

// Use this to retrieve symbols from dynamically loaded libraries.
#define FETCH_SYMBOL(dylib, symbol_name, function_ptr, fail_label) \
  { \
    void* ptr = dlsym(dylib, symbol_name); \
    if (ptr == NULL) \
    { \
      log_urgent("%s: unable to find %s in dynamic library.", __func__, symbol_name); \
      goto fail_label; \
    } \
    *((void**)&(function_ptr)) = ptr; \
  } 

static hypre_lib_t* get_hypre(const char* hypre_dir)
{
  // Try to find HYPRE.
  char hypre_path[FILENAME_MAX+1];
  snprintf(hypre_path, FILENAME_MAX, "%s/libHYPRE%s", hypre_dir, SHARED_LIBRARY_SUFFIX);

  // Try to open libHYPRE and mine it for symbols.
  log_debug("get_hypre: Opening HYPRE library at %s.", hypre_path);
  void* hypre = dlopen(hypre_path, RTLD_NOW);
  if (hypre == NULL)
  {
    char* msg = dlerror();
    polymec_error("get_hypre: %s.", msg);
  }

#if APPLE
  // Mac-specific trick: 
  // If the HYPRE library is parallel, it will contain a reference to MPI symbols, 
  // even if those symbols are not defined within the library itself.
  bool hypre_is_parallel = (dlsym(hypre, "MPI_Abort") != NULL);
#if POLYMEC_HAVE_MPI
  if (!hypre_is_parallel)
  {
    log_urgent("get_hypre: Polymec is configured to use MPI, but HYPRE is not.\n"
               "  HYPRE must be built using -DHYPRE_SEQUENTIAL=OFF."); 
    goto failure;
  }
#else
  if (hypre_is_parallel)
  {
    log_urgent("get_hypre: Polymec is configured serially, but HYPRE is not.\n"
               "  HYPRE must be built using -DHYPRE_SEQUENTIAL=ON."); 
    goto failure;
  }
#endif
#endif
  log_debug("get_hypre: Succeeded.");

  // HYPRE defines a function called hypre_printf if it is configured for "big" (64-bit) integers.
  bool hypre_uses_64bit_indices = (dlsym(hypre, "hypre_printf") != NULL);
  if (sizeof(index_t) == sizeof(int64_t))
  {
    if (!hypre_uses_64bit_indices)
    {
      log_urgent("get_hypre: Since polymec is configured for 64-bit indices,\n"
                 "  HYPRE must be built using -DHYPRE_BIGINT=ON.");
      goto failure;
    }
  }
  else
  {
    if (hypre_uses_64bit_indices)
    {
      log_urgent("get_hypre: Since polymec is configured for 32-bit indices,\n"
                 "  HYPRE must be built using -DHYPRE_BIGINT=OFF.");
      goto failure;
    }
  }

  // Allocate our HYPRE library object.
  hypre_lib_t* hypre_lib = polymec_malloc(sizeof(hypre_lib_t));
  hypre_lib->lib = hypre;

  // Get the symbols.
#define FETCH_HYPRE_SYMBOL(symbol_name) \
  FETCH_SYMBOL(hypre, #symbol_name, hypre_lib->symbol_name, failure);

  FETCH_HYPRE_SYMBOL(HYPRE_GetError);
  FETCH_HYPRE_SYMBOL(HYPRE_ClearAllErrors);

  FETCH_HYPRE_SYMBOL(HYPRE_StructGridCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_StructGridDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_StructGridSetExtents);
  FETCH_HYPRE_SYMBOL(HYPRE_StructGridAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_StructGridSetPeriodic);
  FETCH_HYPRE_SYMBOL(HYPRE_StructGridSetNumGhost);

  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixSetConstantValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixAddToConstantValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixSetBoxValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixAddToBoxValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixSetSymmetric);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixSetConstantEntries);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixSetNumGhost);
  FETCH_HYPRE_SYMBOL(HYPRE_StructMatrixMatVec);

  FETCH_HYPRE_SYMBOL(HYPRE_StructStencilCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_StructStencilDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_StructStencilSetElement);

  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorSetBoxValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorAddToBoxValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorGetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_StructVectorGetBoxValues);

  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetZeroGuess);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetNonZeroGuess);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetNumPreRelax);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGSetNumPostRelax);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_StructSMGGetFinalRelativeResidualNorm);

#undef FETCH_HYPRE_SYMBOL

  log_debug("get_hypre: Got HYPRE symbols.");
  return hypre_lib;

failure:
  dlclose(hypre);
  polymec_free(hypre_lib);
  return NULL;
}

static void free_hypre(hypre_lib_t* hypre)
{
  dlclose(hypre->lib);
}

typedef struct
{
  hypre_lib_t* hypre;
  HYPRE_StructGrid grid;
  HYPRE_StructMatrix A;
  HYPRE_StructStencil stencil;
  HYPRE_StructVector X, B;
} smg_t;

static bool smg_solve(void* context, str_grid_cell_data_t* X)
{
  return false;
}

static void smg_free(void* context)
{
  smg_t* smg = context;
  smg->hypre->HYPRE_StructVectorDestroy(smg->B);
  smg->hypre->HYPRE_StructVectorDestroy(smg->X);
  smg->hypre->HYPRE_StructMatrixDestroy(smg->A);
  smg->hypre->HYPRE_StructStencilDestroy(smg->stencil);
  smg->hypre->HYPRE_StructGridDestroy(smg->grid);
  free_hypre(smg->hypre);
  polymec_free(smg);
}

str_grid_cell_solver_t* hypre_smg_helmholtz_str_grid_cell_solver_new(const char* hypre_dir,
                                                                     MPI_Comm comm,
                                                                     str_grid_t* grid,
                                                                     int num_comps)
{
  smg_t* smg = polymec_malloc(sizeof(smg_t));
  smg->hypre = get_hypre(hypre_dir);

  // Set up the grid.
  smg->hypre->HYPRE_StructGridCreate(comm, 3, &smg->grid);
  // FIXME
  smg->hypre->HYPRE_StructGridAssemble(smg->grid);

  // Set up the stencil.
  smg->hypre->HYPRE_StructStencilCreate(3, 7, &smg->stencil);

  // Now create the Laplacian matrix.
  smg->hypre->HYPRE_StructMatrixCreate(comm, smg->grid, smg->stencil, &smg->A);
  smg->hypre->HYPRE_StructMatrixInitialize(smg->A);
  // FIXME
  smg->hypre->HYPRE_StructMatrixAssemble(smg->A);

  // Finally, set up X and B vectors.
  smg->hypre->HYPRE_StructVectorCreate(comm, smg->grid, &smg->X);
  smg->hypre->HYPRE_StructVectorCreate(comm, smg->grid, &smg->B);
  smg->hypre->HYPRE_StructVectorInitialize(smg->X);
  smg->hypre->HYPRE_StructVectorInitialize(smg->B);
  // FIXME
  smg->hypre->HYPRE_StructVectorAssemble(smg->X);
  smg->hypre->HYPRE_StructVectorAssemble(smg->B);

  // Stash the library.
  str_grid_cell_solver_vtable vtable = {.solve = smg_solve, 
                                        .dtor = smg_free};
  return str_grid_cell_solver_new("HYPRE SMG Helmholtz", comm,
                                  smg, vtable, grid, num_comps);
}

void hypre_smg_helmholtz_str_grid_cell_solver_set_operator(str_grid_cell_solver_t* solver,
                                                           real_t time, 
                                                           real_t alpha, 
                                                           real_t beta, 
                                                           st_func_t* A)
{
}

void hypre_smg_helmholtz_str_grid_cell_solver_set_rhs(str_grid_cell_solver_t* solver, 
                                                      real_t time,
                                                      real_t gamma, 
                                                      real_t delta, 
                                                      st_func_t* k, 
                                                      str_grid_cell_data_t* B, 
                                                      real_t epsilon, 
                                                      st_func_t* C)
{
}

