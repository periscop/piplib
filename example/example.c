/* This is a very simple example of how to use the PipLib inside your programs.
 * You should compile it by typing 'make' (after edition of the makefile), then
 * test it for instance by typing 'more FILE.pol | ./example'. Finally you can
 * compare results given by PIP by typing 'pip32 FILE.dat'
 */

# include <stdio.h>
#include <piplib/piplib.h>

int main()
{ int bignum ;
  PipMatrix  * domain, * context  ;
  PipQuast   * solution ;
  PipOptions * options ;
  
  printf("[PIP2-like future input] Please enter:\n- the context matrix,\n") ;
  context = pip_matrix_read(stdin) ;
  pip_matrix_print(stdout,context) ;

  printf("- the bignum column (start at 0, -1 if no bignum),\n") ;
  fscanf(stdin," %d",&bignum) ;
  printf("%d\n",bignum) ;

  printf("- the constraint matrix.\n") ;
  domain = pip_matrix_read(stdin) ;
  pip_matrix_print(stdout,domain) ;
  printf("\n") ;
  
  options = pip_options_init() ;
  /*options->Maximize = 1 ;*/

  /* The bignum in PIP1 is fixed on the constraint matrix, here is
   * the translation.
   */
  if (bignum > 0)
  bignum += domain->NbColumns - context->NbColumns ;
  
  solution = pip_solve(domain,context,bignum,options) ;

  pip_options_free(options) ;
  pip_matrix_free(domain) ;
  pip_matrix_free(context) ;

  pip_quast_print(stdout,solution,0) ;

  pip_quast_free(solution) ;
  return 0 ;
}
