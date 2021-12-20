/*********************************************************************
Copyright 2008 Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
retains certain rights in this software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of Sandia Corporation nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***********************************************************************/

/* Dt.c */

/* Variant of ls providing YYYYMMDD for latest file, */
/* or complete date and time independent of age */
/* or hex or decimal Unix times. */
/* Author: David M. Gay */
/* Re-implementation of a command I wrote and used for years at Bell Labs. */
#if 0 /* usage of the original version: */

usage: Dt [-[dfmtux]{prefix} [-[dfmtux]{prefix} ...]] file [file...]
	-d = date first*
	-F = include file name first*
	-f = include file name last
	-s = include file size
	-m = month first*
	-t = include time*
	-u = show Unix time
	-# = show Unix time in hex
	-8 = yyyymmdd first*
	-x = postfix
	-y = year first*
 *Latest file only

#endif

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#if defined(_WIN32) || defined(_WIN64)
#include <io.h>
#else
#include <unistd.h>
#endif

 static char *progname;

 static int
usage(int rc)
{
	fprintf(rc ? stderr : stdout, "%s: [option] file [file...]\nOptions:\n\
	-a ==> show last-access time\n\
	-c ==> show last-change time\n\
	-M ==> show last-modification time (default)\n\
	-u ==> show Unix time in decimal\n\
	-# ==> show Unix time in hex\n\
	-8 ==> show YYYYMMDD for latest file\n\
	-f ==> show file name (with -8)\n\
	No option ==> show each file's complete date and time.\n",
		progname);
	return rc;
	}

 int
main(int argc, char **argv)
{
	char *fmt, *s, *txname;
	int rc, showname, unixtime, yyyymmdd;
	struct stat B;
	struct tm T, *tp;
	time_t *st_time, tx;

	progname = *argv++;
	unixtime = yyyymmdd = 0;
	fmt = "%d%02d%02d %02d:%02d:%02d  %s\n";
	st_time = &B.st_mtime;
	while(argc > 1 && *(s = *argv) == '-') {
		--argc; ++argv;
		while(*++s) {
		 switch(*s) {
		  case 'a':
			st_time = &B.st_atime;
			break;
		  case 'c':
			st_time = &B.st_ctime;
			break;
		  case 'f':
			unixtime = 0;
			yyyymmdd = 2;
			fmt = "%d%02d%02d  %s\n";
			break;
		  case 'M':
			st_time = &B.st_mtime;
			break;
		  case 'u':
			yyyymmdd = 0;
			unixtime = 1;
			fmt = "%ld  %s\n";
			break;
		  case '#':
			yyyymmdd = 0;
			unixtime = 1;
			fmt = "%lx  %s\n";
			break;
		  case '8':
			unixtime = 0;
			if (!yyyymmdd) {
				yyyymmdd = 1;
				fmt = "%d%02d%02d\n";
				}
			break;
		  case '?':
			return usage(s[2] != 0);
		  case '-':
			if (!s[1])
				goto optsdone;
			return usage(strcmp(s+1,"help"));
		  default:
			return usage(1);
		  }
		 }
		}
 optsdone:
	rc = 0;
	if (argc < 2)
		return usage(1);
	tx = 0;
	txname = "";
	while(s = *argv++) {
		if (stat(s,&B)) {
			fprintf(stderr, "Cannot stat \"%s\"\n", s);
			rc |= 1;
			continue;
			}
		if (yyyymmdd) {
			if (tx < B.st_mtime) {
				tx = B.st_mtime;
				txname = s;
				}
			}
		else if (unixtime)
			printf(fmt, *st_time, s);
		else if (tp = localtime_r(st_time, &T))
			printf(fmt, tp->tm_year + 1900, tp->tm_mon + 1,
				tp->tm_mday, tp->tm_hour, tp->tm_min, tp->tm_sec, s);
		else {
			fprintf(stderr, "localtime_r failure!\n");
			rc |= 2;
			}
		}
	if (yyyymmdd)
		if (tp = localtime_r(&tx, &T))
			printf(fmt, tp->tm_year + 1900, tp->tm_mon + 1,
				tp->tm_mday, txname);
		else {
			fprintf(stderr, "localtime_r failure!\n");
			rc |= 2;
			}
	return rc;
	}
