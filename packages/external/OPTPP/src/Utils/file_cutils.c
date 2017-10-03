
/*
 * C file utilities needed by PDS
 *
 * WARNING:  These are PDS-specific!!!
 *
 */

#include <stdio.h>
#include "pds.h"

int bin_open(char *filename, FILE **fp)
{
  /* open read/write, create for backward compatibility 
     (read likely not needed); don't specify permissions */
  *fp = fopen(filename, "w+");
  if (*fp == NULL)
    return -1;
  else
    return 0;
}

int bin_close(FILE *fp)
{
  int error;

  error = fclose(fp);
  return error;
}
