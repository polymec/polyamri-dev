// Copyright (c) 2014-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file really only checks that the version of polyamri matches that 
// of polymec.

#include "polyamri/polyamri.h"
#include "polyamri/polyamri_version.h"

void polyamri_version_fprintf(const char* exe_name, FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "%s v%s\n", exe_name, POLYAMRI_VERSION);
  (void)POLYAMRI_GIT_DIFFS;
}

// We're not relying on MPI-3 just yet.
//// We use some fancy MPI-3 stuff in this here library.
//#if MPI_VERSION < 3
//#error "MPI-3 is required for polyamri."
//#endif

#if POLYAMRI_MAJOR_VERSION > 0
// The major version of polyamri must match that of polymec.
#include "core/polymec_version.h"
#if POLYMEC_MAJOR_VERSION != POLYAMRI_MAJOR_VERSION
#error "The installed major version of polymec differs from that of polyamri. Please make sure these versions match."
#endif
#endif

