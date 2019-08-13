#include "file_io.h"

nuMatParams* read_parameter_file(char *param_file) {
  FILE *in_file;
  nuMatParams *nmp = malloc(sizeof(*nmp));
  in_file = fopen(param_file, "r");
  nmp->density_file = malloc(sizeof(char)*100);
  fscanf(in_file, "%s\n", nmp->density_file);
  nmp->out_file_base = malloc(sizeof(char)*100);
  fscanf(in_file, "%s\n", nmp->out_file_base);
  fscanf(in_file, "%d\n", &nmp->n_body);
  if ((nmp->n_body != 1) && (nmp->n_body != 2)) {printf("Incorrect request for %d body operator, please type only 1 or 2 here.\n", nmp->n_body); exit(0);}
  fscanf(in_file, "%d,%d\n", &nmp->j_op, &nmp->t_op);
  fscanf(in_file, "%d\n", &nmp->spec_dep);
  if ((nmp->spec_dep != 0) && (nmp->spec_dep) != 1) {printf("Invalid flag for spectator dependence %d, please type only 0 or 1 here\n", nmp->spec_dep); exit(0);}
  return nmp;
}

