// Copyright (c) 2014-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "polyamri/amr_patch_set.h"

struct amr_patch_set_t
{
  ptr_array_t* patches;
  ptr_array_t* bboxes;
  ptr_array_t* data;
};

amr_patch_set_t* amr_patch_set_new()
{
  amr_patch_set_t* patches = polymec_malloc(sizeof(amr_patch_set_t));
  patches->patches = ptr_array_new();
  patches->bboxes = ptr_array_new();
  patches->data = ptr_array_new();
  return patches;
}

void amr_patch_set_add_with_data(amr_patch_set_t* patches, amr_patch_t* patch, bbox_t* bbox, void* data, void (*data_dtor)(void*))
{
  ptr_array_append_with_dtor(patches->patches, patch, DTOR(amr_patch_free));

  // Create our own copy of the bounding box.
  bbox_t* my_bbox = polymec_malloc(sizeof(bbox_t));
  my_bbox->x1 = bbox->x1;
  my_bbox->x2 = bbox->x2;
  my_bbox->y1 = bbox->y1;
  my_bbox->y2 = bbox->y2;
  my_bbox->z1 = bbox->z1;
  my_bbox->z2 = bbox->z2;
  ptr_array_append_with_dtor(patches->bboxes, my_bbox, DTOR(polymec_free));
  
  // Associate data.
  if (data != NULL)
    ptr_array_append_with_dtor(patches->data, data, data_dtor);
  else
    ptr_array_append_with_dtor(patches->data, NULL, NULL);
}

void amr_patch_set_add(amr_patch_set_t* patches, amr_patch_t* patch, bbox_t* bbox)
{
  amr_patch_set_add_with_data(patches, patch, bbox, NULL, NULL);
}

void amr_patch_set_free(amr_patch_set_t* patches)
{
  ptr_array_free(patches->patches);
  ptr_array_free(patches->bboxes);
  ptr_array_free(patches->data);
  polymec_free(patches);
}

int amr_patch_set_size(amr_patch_set_t* patches)
{
  return patches->patches->size;
}

bool amr_patch_set_next(amr_patch_set_t* patches, int* pos, amr_patch_t** patch, bbox_t** bbox)
{
  if (*pos >= patches->patches->size)
    return false;
  *patch = patches->patches->data[*pos];
  *bbox = patches->bboxes->data[*pos];
  ++(*pos);
  return true;
}

bool amr_patch_set_next_with_data(amr_patch_set_t* patches, int* pos, amr_patch_t** patch, bbox_t** bbox, void** data)
{
  if (*pos >= patches->patches->size)
    return false;
  *patch = patches->patches->data[*pos];
  *bbox = patches->bboxes->data[*pos];
  *data = patches->data->data[*pos];
  ++(*pos);
  return true;
}

