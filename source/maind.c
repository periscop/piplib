/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                  maind.h                                   *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
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

char version[]="Version E.2 $Revision: 1.3 $\n";


extern long int cross_product, limit ;
extern int allocation, comptage, verbose ;
extern FILE * dump ;
extern int compa_count ;
extern char dump_name[] ;


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
 #if defined(LINEAR_VALUE_IS_MP)
 int x ;
 #else
 Entier D, x;
 #endif
 
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_init_set_si(UN, 1);
 mpz_init_set_si(ZERO, 0);
 #else
 UN   = VAL_UN ;
 ZERO = VAL_ZERO ;
 #endif

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
      if(dscanf(in, &x) < 0){escape(in, out, 1); continue;}
      else 
        nvar = (int) x;
      if(dscanf(in, &x) < 0){escape(in, out, 1); continue; }
      else
        nparm = (int) x;
      if(dscanf(in, &x) < 0){escape(in, out, 1); continue; }
      else 
        ni = (int) x;
      if(dscanf(in, &x) < 0){escape(in, out, 1); continue; }
      else
        nc = (int) x;
      if(dscanf(in, &x) < 0){escape(in, out, 1); continue; }
      else 
        bigparm = (int) x;
      if(dscanf(in, &x) < 0){escape(in, out, 1); continue; }
      else 
        nq = (int) x;
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
           #if defined(LINEAR_VALUE_IS_MP)
	   traiter(ctxt, NULL, True, nparm, 0, nc, 0, -1);
	   #else
	   traiter(ctxt, NULL, True, UN, nparm, 0, nc, 0, -1);
           #endif
	   non_vide = is_not_Nil(p);
	   sol_reset(p);
	  }
      else non_vide = True;
      if(non_vide) {
	   compa_count = 0;
           #if defined(LINEAR_VALUE_IS_MP)
	   traiter(ineq, context, nq, nvar, nparm, ni, nc, bigparm);
	   putc(' ',out);
	   mpz_out_str(out, 10, ineq->determinant);
	   #else
	   D = traiter(ineq, context, nq, UN, nvar, nparm, ni, nc, bigparm);
	   putc(' ',out);
	   fprintf(out, FORMAT, D);
           #endif
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

