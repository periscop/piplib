/* This is a very simple example of how to use the PipLib inside your programs.
 * You should compile it by typing 'make' (after edition of the makefile), then
 * test it for instance by typing 'more FILE.pol | ./example'. Finally you can
 * compare results given by PIP by typing 'pip32 FILE.dat'
 */

# include <stdio.h>
# include <piplib/piplib32.h>

int main()
{ PipMatrix  * domain, * context  ;
  PipQuast   * solution ;
  PipOptions * options ;
  
  printf("Enter the Matrices :\n") ;
  domain = pip_matrix_read(stdin) ;
  pip_matrix_print(stdout,domain) ;
  
  context = pip_matrix_read(stdin) ;
  pip_matrix_print(stdout,context) ;
  printf("\n") ;
  
  options = pip_options_init() ;

  solution = pip_solve(domain,context,-1,options) ;

  pip_options_free(options) ;
  pip_matrix_free(domain) ;
  pip_matrix_free(context) ;

  pip_quast_print(stdout,solution,0) ;

  pip_quast_free(solution) ;
  return 0 ;
}
