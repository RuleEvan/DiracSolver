#ifndef FILE_IO_H
#define FILE_IO_H
#include "av18.h"

typedef struct nuMatParams
{
  char *density_file, *out_file_base;
  int n_body;
  int spec_dep;
  int j_op, t_op;
  
} nuMatParams;

nuMatParams* read_parameter_file(char* parameter_file);
#endif
