/* $Id: var_bnds_only.c 217 2013-11-25 21:59:49Z tplante $
/* $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/1-var-bnds-only/var_bnds_only.c $

/* @HEADER
 * ************************************************************************
 *
 *         HOPSPACK: Hybrid Optimization Parallel Search Package
 *                 Copyright 2009-2013 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This file is part of HOPSPACK.
 *
 * HOPSPACK is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.  If not, see http://www.gnu.org/licenses/.
 *
 * Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
 *                 or Todd Plantenga (tplante@sandia.gov)
 *
 * ************************************************************************
 * @HEADER
 */

/*!
  @file var_bnds_only.c
  @brief Evaluates function values for var_bnds_only example 1.

  This example provides a C language evaluation program for an optimization
  problem with one objective and simple bound constraints.  The problem is:

    minimize  f(x1,x2) = x1^2 + (2 * x2^2)

    subject to  -1 <= x1 <= 1
                -1 <= x2 <= 1

  The exact solution point is x1=0, x2=0, with objective f(0,0) = 0.

  The input file has format:
    F
    2
    1.000e-01
    2.000e-01
*/

#include <stdio.h>
#include <string.h>


/* --------------------------------------------------------------------
 *  Platform/build specific symbols
 * --------------------------------------------------------------------
 */
#if defined(_WIN32)
  #if (_MSC_VER >= 1400)
    //---- RECENT WINDOWS MSVC COMPILERS OFFER SECURE STRING FUNCTIONS.
    #define HAVE_MSVC_SECURE_STRING_FNS
  #endif
#endif


/* --------------------------------------------------------------------
 *  Internal Function evaluate_f
 * --------------------------------------------------------------------
 */
static double  evaluate_f (const double  x1,
                           const double  x2)
{

    double  f = (x1 * x1) + (2.0 * x2 * x2);
    return( f );
}


/* --------------------------------------------------------------------
 *  Internal Function read_input_file
 * --------------------------------------------------------------------
 */
static int  read_input_file (const char   * const  szExeName,
                             const char   * const  szFileName,
                                   double * const  pX1,
                                   double * const  pX2)
{
    #define REQBUFLEN  10
    char    reqBuf[REQBUFLEN];

    FILE *  fp;
    int     i;

    /* OPEN THE INPUT FILE. */
   #if defined(HAVE_MSVC_SECURE_STRING_FNS)
    i = fopen_s (&fp, szFileName, "r");
    if (i != 0)
   #else
    fp = fopen (szFileName, "r");
    if (fp == NULL)
   #endif
    {
        fprintf (stderr, "%s - ERROR opening input file '%s'.\n",
                 szExeName, szFileName);
        return( -1 );
    }

    /* READ AND VERIFY THE INPUT REQUEST TYPE. */
    fgets (reqBuf, REQBUFLEN - 1, fp);
    for (i = 0; i < REQBUFLEN; i++)
    {
        if ((reqBuf[i] == '\n') || (reqBuf[i] == '\r'))
        {
            reqBuf[i] = 0;
            break;
        }
    }
    if ((reqBuf[0] != 'F') || (reqBuf[1] != 0))
    {
        fprintf (stderr, "%s - ERROR reading request type from '%s'.\n",
                 szExeName, szFileName);
        fprintf (stderr, "  Read '%s'\n", reqBuf);
        fclose (fp);
        return( -2 );
    }

    /* READ THE LENGTH OF x. */
   #if defined(HAVE_MSVC_SECURE_STRING_FNS)
    if ((fscanf_s (fp, "%d", &i) != 1) || (i != 2))
   #else
    if ((fscanf (fp, "%d", &i) != 1) || (i != 2))
   #endif
    {
        fprintf (stderr, "%s - ERROR reading n from '%s'.\n",
                 szExeName, szFileName);
        fclose (fp);
        return( -3 );
    }

    /* READ x. */
   #if defined(HAVE_MSVC_SECURE_STRING_FNS)
    if ((fscanf_s (fp, "%le", pX1) != 1) || (fscanf_s (fp, "%le", pX2) != 1))
   #else
    if ((fscanf (fp, "%le", pX1) != 1) || (fscanf (fp, "%le", pX2) != 1))
   #endif
    {
        fprintf (stderr, "%s - ERROR reading x from '%s'.\n",
                 szExeName, szFileName);
        fclose (fp);
        return( -4 );
    }

    fclose (fp);

    /* RETURN SUCCESS. */
    return( 0 );
}


/* --------------------------------------------------------------------
 *  Internal Function write_output_file
 * --------------------------------------------------------------------
 */
static int  write_output_file (const char   * const  szExeName,
                               const char   * const  szFileName,
                               const double          f)
{
    FILE *  fp;

   #if defined(HAVE_MSVC_SECURE_STRING_FNS)
    int  i;
    i = fopen_s (&fp, szFileName, "w");
    if (i != 0)
   #else
    fp = fopen (szFileName, "w");
    if (fp == NULL)
   #endif
    {
        fprintf (stderr, "%s - ERROR opening output file '%s'.\n",
                 szExeName, szFileName);
        return( -10 );
    }

    /* WRITE THE NUMBER OF OBJECTIVES AND THEIR VALUES TO THE OUTPUT
     * WITH MAXIMUM PRECISION (16 SIGNIFICANT DIGITS).
     */
    fprintf (fp, "1\n");
    fprintf (fp, "%23.15e\n", f);

    fclose (fp);

    return( 0 );
}


/* --------------------------------------------------------------------
 *  Main routine for evaluation executable.
 *
 *  Each execution reads x1 and x2 from a file, evaluates the objective
 *  function, and writes the function value to an output file.  HOPSPACK
 *  generates the file names and decides when to invoke this executable.
 *  If there is an error, either create no output file, or create an
 *  output file with an error message.
 *
 *  @param argc  Number of command line arguments.
 *  @param argv  Command line arguments (input and output file names).
 *  @return      0 if successful.
 * --------------------------------------------------------------------
 */
int  main (const int           argc,
           const char * const  argv[])
{
    double  x1, x2;
    double  f;
    int     nRetStatus;


    /* CHECK THE COMMAND LINE ARGUMENTS.
     * IN THIS EXAMPLE ONLY THE INPUT AND OUTPUT FILE NAMES ARE USED.
     */
    if (argc != 5)
    {
        fprintf (stderr, "usage: %s <input file> <output file> <tag> <type>\n",
                 argv[0]);
        return( -100 );
    }

    nRetStatus = read_input_file (argv[0], argv[1], &x1, &x2);
    if (nRetStatus != 0)
        return( nRetStatus );

    f = evaluate_f (x1, x2);

    nRetStatus = write_output_file (argv[0], argv[2], f);
    return( nRetStatus );
}
