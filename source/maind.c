/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996                           *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.							              *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <stdlib.h>
#ifdef __TURBOC__
#include <dir.h>
#endif
#define min(x,y) ((x) < (y)? (x) : (y))

#include <piplib/piplib.h>

#ifdef UNIX
#include <sys/times.h>
struct tms chrono;
#endif

char version[]="Version E.2\n";

long int cross_product, limit;
int allocation, comptage;
int verbose = 0;
int profondeur = 0;
int compa_count;

FILE *dump = NULL;
char dump_name[] = "XXXXXX";

#define INLENGTH 1024

char inbuff[INLENGTH];
int inptr = 256;
int proviso = 0;


int dgetc(FILE *foo)
{
 char *p;
 if(inptr >= proviso)
   {p = fgets(inbuff, INLENGTH, foo);
    if(p == NULL) return EOF;
    proviso = min(INLENGTH, strlen(inbuff));
    inptr = 0;
    if(verbose > 0) fprintf(dump, "-- %s", inbuff);
  }
 return inbuff[inptr++];
}

int dscanf(FILE *foo, char *format, Entier *val)
{
 char * p;
 int c;
 for(;inptr < proviso; inptr++)
   if(inbuff[inptr] != ' ' && inbuff[inptr] != '\n' && inbuff[inptr] != '\t')
				break;
 while(inptr >= proviso)
   {p = fgets(inbuff, 256, foo);
    if(p == NULL) return EOF;
    proviso = strlen(inbuff);
    if(verbose > 0) {
      fprintf(dump, ".. %s", inbuff);
      fflush(dump);
    }
    for(inptr = 0; inptr < proviso; inptr++)
       if(inbuff[inptr] != ' '
       && inbuff[inptr] != '\n'
       && inbuff[inptr] != '\t') break;
  }
 if(sscanf(inbuff+inptr, FORMAT, val) != 1) return -1;
 
 for(; inptr < proviso; inptr++)
	if((c = inbuff[inptr]) != '-' && !isdigit(c)) break;
 return 0;
}

void balance(FILE *foo, FILE *bar)
{
 int level = 0;
 int c;
 while((c = dgetc(foo)) != EOF)
     {
      switch(c)
	  {case '(' : level++; break;
	   case ')' : if(--level == 0) return;
	  }
      putc(c, bar);
     }
}

void escape(FILE *foo, FILE *bar, int level)
{int c;
 while((c = dgetc(foo)) != EOF)
   switch(c)
     {case '(' : level ++; break;
     case ')' : if(--level == 0)
		     { fprintf(bar, "\nSyntax error\n)\n");
		       return;
		     }
     }
}

char * getenv();

int main(int argc, char *argv[])

{
 FILE *in, *out;
 Tableau *ineq, *context, *ctxt;
 int nvar, nparm, ni, nc, bigparm;
 int nq; char * g;
 int simple = 0;
 struct high_water_mark hq;
 int c, non_vide;
 int p, q, xq;
 long temps;
 char *date;
 Entier D, x;

 in = stdin; out = stdout;
 p = 1;
 if(argc > 1)
     if(strcmp(argv[1], "-s") == 0)
	 {verbose = -1;
	  p = 2;
	 }
     else if(strcmp(argv[1], "-v") == 0)
	 {verbose = 1;
	  p = 2;
	  g = getenv("DEBUG");
	  if(g && *g)
	    {dump = fopen(g, "w");
	     if(dump == NULL)
	       {fprintf(stderr, "%s unaccessible\n", g);
		verbose = 0;
	      }
	   }
	  else
	    {
	     mkstemp(dump_name);
	     dump = fopen(dump_name, "w");
	    }
	  }
 if(verbose >= 0) fprintf(stderr, version);
 if(argc > p)
     if(strcmp(argv[p], "-z") == 0) {
	  simple = 1;
	  p++;
	  }
 if(argc>p)
     {in = fopen(argv[p], "r");
      if(in == NULL)
	  {fprintf(stderr, "%s unaccessible\n", argv[p]);
	   exit(1);
	  }
     }
 p++;
 if(argc>p)
     {out = fopen(argv[p], "w");
      if(out == NULL)
	  {fprintf(stderr, "%s unaccessible\n", argv[p]);
	   exit(2);
	  }
     }
 limit = 0L;
 p++;
 if(argc > p) limit = atol(argv[p]);
 sol_init();
 tab_init();
 while((c = dgetc(in)) != EOF)
     {if(c != '(') continue;
      fprintf(out, "(");
      balance(in, out);
      if(dscanf(in, FORMAT, &x) < 0){escape(in, out, 1); continue;}
      else nvar = (int) x;
      if(dscanf(in, FORMAT, &x) < 0){escape(in, out, 1); continue; }
      else nparm = (int) x;
      if(dscanf(in, FORMAT, &x) < 0){escape(in, out, 1); continue; }
      else ni = (int) x;
      if(dscanf(in, FORMAT, &x) < 0){escape(in, out, 1); continue; }
      else nc = (int) x;
      if(dscanf(in, FORMAT, &x) < 0){escape(in, out, 1); continue; }
      else bigparm = (int) x;
      if(dscanf(in, FORMAT, &x) < 0){escape(in, out, 1); continue; }
      else nq = (int) x;
      if(verbose > 0) {fprintf(dump, "%d %d %d %d %d %d\n",nvar, nparm, ni, nc,
			      bigparm, nq);
		       fflush(dump);
		      }
      cross_product = 0;
      hq = tab_hwm();
      if(verbose > 0) {fprintf(dump, "hwm %x\n", g);
                       fflush(dump);
                      }
      ineq = tab_get(in, ni, nvar+nparm+1, nvar);
      if(ineq == NULL){escape(in, out, 2); continue;}
      context = tab_get(in, nc, nparm+1, 0);
      if(ineq == NULL){escape(in, out, 2); continue;}
      xq = p = sol_hwm();
/* verification de la non-vacuite' du contexte */
      if(nc)
	  {ctxt = expanser(context, nparm, nc, nparm+1, nparm, 0, 0);
	   traiter(ctxt, NULL, True, UN, nparm, 0, nc, 0, -1);
	   non_vide = is_not_Nil(p);
	   sol_reset(p);
	  }
      else non_vide = True;
      if(non_vide) {
	   compa_count = 0;
	   D = traiter(ineq, context, nq, UN, nvar, nparm, ni, nc, bigparm);
	   putc(' ',out);
	   fprintf(out, FORMAT, D);
	   fputs(" )",out);
	   if(simple) sol_simplify(xq);
	   q = sol_hwm();
	   while((xq = sol_edit(out, xq)) != q);
	   sol_reset(p);
	  }
      else fprintf(out, "void\n");
      tab_reset(hq);
      if(verbose > 0) fflush(dump);
      fprintf(out, ")\n");
      fflush(out);
      if(verbose >= 0) fprintf(stderr,"cross : %ld, alloc : %d, compa : %d\n\r",
				cross_product, allocation, compa_count);
      comptage++;
     }
#ifdef UNIX
 times(& chrono);
 fprintf(stderr, "n %d u %d''' s %d'''\r\n",
	comptage, chrono.tms_utime, chrono.tms_stime);
#endif
 exit(0);
}

