// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYAMRI_H
#define POLYAMRI_H

#include "core/polymec.h"

// Returns true if the given number is a power of 2, false otherwise.
static inline bool is_power_of_2(int number)
{
  int x = 1;
  while (x < number)
  {
    x *= 2;
    if (x == number)
      return true;
  }
  return false;
}

#endif

