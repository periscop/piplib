/********************************************************/
/* Part of MultiPrecision PIP port by Zbigniew Chamski  */
/* and Paul Feautrier.                                  */
/* Based on straight PIP E.1 version by Paul Feautrier  */
/* <Paul.Feautrier@inria.fr>                            */
/*                                                      */
/* and a previous port (C) Zbigniew CHAMSKI, 1993.      */
/* <Zbigniew.Chamski@philips.com>                       */
/*                                                      */
/* Copying subject to the terms and conditions of the   */
/* GNU General Public License.                          */
/*                                                      */
/* Send questions, bug reports and fixes to:            */
/* <Paul.Feautrier@inria.fr>                            */
/********************************************************/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <stdlib.h>


#define min(x,y) ((x) < (y)? (x) : (y))

#include <piplib/piplib.h>

#include <sys/times.h>
struct tms chrono;

long int cross_product, limit;
int allocation, comptage;
int verbose = 0;
int profondeur = 0;
float the_clock = 0.0;
int compa_count;

FILE *dump = NULL;
char dump_name[] = "PipXXXXXX";

#define INLENGTH 1024

char inbuff[INLENGTH];
int inptr = 256;
int proviso = 0;

Entier UN;
Entier ZERO;

char version[] = "Version MP.1\n";

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

int dscanf(FILE *foo, char *format, int *val)
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
 int x;
 mpz_init_set_si(UN, 1);
 mpz_init_set_si(ZERO, 0);

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
	       {fprintf(stdout, "%s unaccessible\n", g);
		verbose = 0;
	      }
	   }
	  else
	    {
	     mkstemp(dump_name);
	     dump = fopen(dump_name, "w");
	    }
	  }
 if(verbose >= 0) fprintf(stdout, version);
 if(argc > p)
     if(strcmp(argv[p], "-z") == 0) {
	  simple = 1;
	  p++;
	  }
 if(argc>p)
     {in = fopen(argv[p], "r");
      if(in == NULL)
	  {fprintf(stdout, "%s unaccessible\n", argv[p]);
	   exit(1);
	  }
     }
 p++;
 if(argc>p)
     {out = fopen(argv[p], "w");
      if(out == NULL)
	  {fprintf(stdout, "%s unaccessible\n", argv[p]);
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
      if(dscanf(in, "d", &x) < 0){escape(in, out, 1); continue;}
      else nvar = x;
      if(dscanf(in, "d", &x) < 0){escape(in, out, 1); continue; }
      else nparm = x;
      if(dscanf(in, "d", &x) < 0){escape(in, out, 1); continue; }
      else ni = x;
      if(dscanf(in, "d", &x) < 0){escape(in, out, 1); continue; }
      else nc = x;
      if(dscanf(in, "d", &x) < 0){escape(in, out, 1); continue; }
      else bigparm = x;
      if(dscanf(in, "d", &x) < 0){escape(in, out, 1); continue; }
      else nq = x;
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
	   traiter(ctxt, NULL, True, nparm, 0, nc, 0, -1);
	   non_vide = is_not_Nil(p);
	   sol_reset(p);
	  }
      else non_vide = True;
      if(non_vide) {
	   compa_count = 0;
	   traiter(ineq, context, nq, nvar, nparm, ni, nc, bigparm);
	   putc(' ',out);
	   mpz_out_str(out, 10, ineq->determinant);
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
      if(verbose >= 0)fprintf(stdout,"cross : %ld, alloc : %d, compa : %d\n\r",
				cross_product, allocation, compa_count);
      comptage++;
     }
 if(verbose >= 0){
      times(& chrono);
      fprintf(stdout, "n %d u %d''' s %d'''\r\n",
	comptage, chrono.tms_utime, chrono.tms_stime);
      }
 exit(0);
}

