/********************************************************/
/* Part of MultiPrecision PIP port (C) Zbigniew Chamski */
/* and Paul Feautrier.                                  */
/* Based on straight PIP version E.1 by Paul Feautrier  */
/* (<Paul.Feautrier@inria.fr>)                          */
/*                                                      */
/* and a previous port (C) Zbigniew CHAMSKI, 1993.      */
/* (Zbigniew.Chamski@philips.com)                       */
/*                                                      */
/* Copying subject to the terms and conditions of the   */
/* GNU General Public License.                          */
/*                                                      */
/* Send questions, bug reports and fixes to:            */
/*   <Paul.Feautrier@inria.fr>                          */
/********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <piplib/piplib.h>

extern long int cross_product, limit;
extern int verbose;
extern FILE *dump;

struct S
    {int flags;
     Entier param1, param2;
    };

/* In some cases, param1 and param2 are used to store small integers.
   These integers are nevertheless converted to bignums in the
   interest of simplicity. Since the handling of the solution is 
   a negligible part of the algorithm, the corresponding overhead is
   deemed negligible.
*/


#define Free 0
#define Nil  1
#define If   2
#define List 3
#define Form 4
#define New  5
#define Div  6
#define Val  7
#define Error 8

struct S sol_space[SOL_SIZE];
static int sol_free;

/* pgcd has been removed and replaced by mpz_gcd
 */

void sol_init(void)
{
 sol_free = 0;
}

int sol_hwm()
{
 return(sol_free);
}

void sol_reset(p)
int p;
{int i;
 if(p<0 || p>=SOL_SIZE)
     {fprintf(stdout, "Syserr : sol_reset : Memory allocation error\n");
      exit(40);
     }
 for(i=p; i<sol_free; i++){
   mpz_clear(sol_space[i].param1);
   mpz_clear(sol_space[i].param2);
 }
 sol_free = p;
}

struct S *sol_alloc(void)
{struct S *r;
 r = sol_space + sol_free;
 r->flags = Free;
 mpz_init(r->param1);
 mpz_init(r->param2);
 sol_free++;
 if(sol_free >= SOL_SIZE)
     {fprintf(stdout, "The solution is too complex! : sol\n");
      exit(26);
     }
     return(r);
}

void sol_nil(void)
{
 struct S * r;
 r = sol_alloc();
 r -> flags = Nil;
 if(verbose > 0)
   {fprintf(dump, "\nNil");
    fflush(dump);
  }
}

void sol_error(int c)
{
 struct S *r;
 r = sol_alloc();
 r->flags = Nil;
 mpz_set_si(r->param1, c);
 if(verbose > 0) {
     fprintf(dump, "Erreur %d\n", c);
     fflush(dump);
     }
}

int is_not_Nil(p)
int p;
{
 return(sol_space[p].flags != Nil);
}

void sol_if(void)
{
 struct S *r;
 r = sol_alloc();
 r -> flags = If;
 if(verbose > 0) {
     fprintf(dump, "\nIf ");
     fflush(dump);
   }
}

void sol_list(n)
int n;
{struct S * r;
 r = sol_alloc();
 r->flags = List;
 mpz_set_si(r->param1, n);
 if(verbose > 0) {
     fprintf(dump, "\nList %d ", n);
     fflush(dump);
}
}

void sol_forme(l)
int l;
{
 struct S *r;
 r = sol_alloc();
 r -> flags = Form;
 mpz_set_ui(r -> param1, l);
 if(verbose > 0) {
     fprintf(dump, "\nForme %d ", l);
     fflush(dump);
   }
}

void sol_new(k)
int k;
{
 struct S *r;
 if(verbose > 0) {
     fprintf(dump, "New %d ", k);
     fflush(dump);
   }
 r = sol_alloc();
 r -> flags = New;
 mpz_set_ui(r -> param1, k);
}

void sol_div()
{
 struct S *r;
 r = sol_alloc();
 r -> flags = Div;
 if(verbose > 0) {
     fprintf(dump, "Div ");
     fflush(dump);
   }
}

void sol_val(n, d)
Entier n, d;
{
 struct S *r;

 r = sol_alloc();
 r -> flags = Val;
 mpz_set(r -> param1, n);
 mpz_set(r -> param2, d);
 if(verbose > 0) {
   fputs("val(", dump);
   mpz_out_str(dump, 10, n);
   putc('/', dump);
   mpz_out_str(dump, 10, d);
   fputs(") ", dump);
  }
}

int skip(int);

/* a` partir d'un point de la solution, sauter un objet
bien forme' ainsi qu'un e'ventuel New et pointer sur l'objet
suivant */

int skip_New (int i)
{
 if(sol_space[i].flags != New) return i;
 i = skip(i+1);      /* sauter le Div */
 return i;
}
/* au lancement, i indexe une cellule qui est la te^te d'un objet.
   la valeur retourne'e est la te^te de l'objet qui suit. Les
   objets de type New sont e'limine's                                */

int skip (int i)
{int n, f;
 while((f = sol_space[i].flags) == Free || f == Error) i++;
 switch (sol_space[i].flags) {
 case Nil : case Val : i++; break;
 case New : i = skip_New(i); break;
 case If : i = skip(i+1);        /* sauter le pre'dicat */
	   i = skip(i);          /* sauter le vrai */
	   i = skip(i); break;   /* sauter le faux */
 case List : case Form : n = mpz_get_si(sol_space[i].param1);
	   i++;
	   while(n--) i = skip(i);
	   break;
 case Div : i = skip(i+1);       /* sauter la forme */
	    i = skip(i);         /* sauter le diviseur */
	    break;
 default : fprintf(stdout,
	      "Syserr : skip : unknown %d\n", sol_space[i].flags);
 }
 return skip_New(i);
}
/* simplification de la solution : e'limination des constructions
   (if p () ()). N'est en service qu'en pre'sence de l'option -z */

void sol_simplify(int i)
{int j, k, l;
 if(sol_space[i].flags == If) {
     j = skip(i+1);        /* j : de'but de la partie vraie */
     k = skip(j);          /* k : de'but de la partie fausse */
     sol_simplify(k);
     sol_simplify(j);
     if(sol_space[j].flags == Nil && sol_space[k].flags == Nil) {
	 sol_space[i].flags = Nil;
	 if(k >= sol_free - 1) sol_free = i+1;
	 else for(l = i+1; l<=k; l++) sol_space[l].flags = Free;
       }
   }

}
/* e'dition de la solution */

int sol_edit(FILE *foo, int i)
{int j, n;
 struct S *p;
 Entier N, D, d;
 mpz_init(N);
 mpz_init(D);
 mpz_init(d);
 p = sol_space + i;
 for(;;) {
   if(p->flags == Free) {
     p++;
     i++;
     continue;
   }
   if(p->flags == New) {
     n = mpz_get_si(p->param1);
     fprintf(foo, "(newparm %d ", n);
     if(verbose>0)fprintf(dump, "(newparm %d ", n);
     i = sol_edit(foo, ++i);
     p = sol_space +i;
     fprintf(foo, ")\n");
     if(verbose>0)fprintf(dump, ")\n");
     continue;
   }
   break;
 }
 switch(p->flags)
   {
   case Nil : fprintf(foo, "()\n");
     if(verbose>0)fprintf(dump, "()\n");
     i++; break;
   case Error : n = mpz_get_si(p->param1);
     fprintf(foo, "Error %d\n", n);
     if(verbose>0)fprintf(dump, "Error %d\n", n);
     i++; break;
   case If  : fprintf(foo, "(if ");
     if(verbose>0)fprintf(dump, "(if ");
     i = sol_edit(foo, ++i);
     i = sol_edit(foo, i);
     i = sol_edit(foo, i);
     fprintf(foo, ")\n");
     if(verbose>0)fprintf(dump, ")\n");
     break;
   case List: fprintf(foo, "(list ");
     if(verbose>0)fprintf(dump, "(list ");
     n = mpz_get_si(p->param1);
     i++;
     while(n--) i = sol_edit(foo, i);
     fprintf(foo, ")\n");
     if(verbose>0)fprintf(dump, ")\n");
     break;
   case Form: fprintf(foo, "#[");
     if(verbose>0)fprintf(dump, "#[");
     n = mpz_get_si(p->param1);
     for(j = 0; j<n; j++){
       i++; p++;
       mpz_set(N, p->param1); mpz_set(D, p->param2);
       mpz_gcd(d, N, D);
       if(mpz_cmp(d, D) == 0){
	 putc(' ', foo);
	 mpz_divexact(N, N, d);
	 mpz_out_str(foo, 10, N);
         if(verbose>0){
	   putc(' ', dump);
	   mpz_out_str(dump, 10, N);
	 }
       }
       else{
	 mpz_divexact(N, N, d);
	 mpz_divexact(D, D, d);
	 putc(' ', foo);
	 mpz_out_str(foo, 10, N);
	 putc('/', foo);
	 mpz_out_str(foo, 10, D);
	 if(verbose>0){
	   putc(' ', dump);
	   mpz_out_str(dump, 10, N);
	   putc('/', dump);
	   mpz_out_str(dump, 10, D);
	 }
       }
     }
     fprintf(foo, "]\n");
     if(verbose>0)
       fputs("]\n", dump);
     i++;
     break;
   case Div : fprintf(foo, "(div ");
     if(verbose>0)fprintf(dump, "(div ");
     i = sol_edit(foo, ++i);
     i = sol_edit(foo, i);
     fprintf(foo, ")\n");
     if(verbose>0)fprintf(dump, ")\n");
     break;
   case Val : mpz_set(N, p->param1); mpz_set(D, p->param2);
     mpz_gcd(d, N, D);
     if(mpz_cmp(d, D) == 0){
       mpz_divexact(N, N, d);
       putc(' ', foo);
       mpz_out_str(foo, 10, N);
       if(verbose>0){
	 putc(' ', dump);
	 mpz_out_str(dump, 10, N);
       }
     }
     else{
       mpz_divexact(N, N, d);
       mpz_divexact(D, D, d);
       putc(' ', foo);
       mpz_out_str(foo, 10, N);
       fprintf(foo, "/");
       mpz_out_str(foo, 10, D);
       if(verbose>0){
	 putc(' ', dump);
	 mpz_out_str(dump, 10, N);
	 fprintf(dump, "/");
	 mpz_out_str(dump, 10, D);
       }
     }
     i++;
     break;
   default  : fprintf(foo, "Inconnu : sol\n");
     if(verbose>0)fprintf(dump, "Inconnu : sol\n");
   }
 mpz_clear(d);
 mpz_clear(D);
 mpz_clear(N);
 return(i);
}





/* Fonction sol_vector_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a un Vector dans la grammaire dans une structure de type PipVector. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipVector contenant les informations de ce Vector.
 * Premiere version : Ced. 20 juillet 2001. 
 * 24 octobre 2002 : Premiere version MP. 
 */
PipVector * sol_vector_edit(int * i)
{ int j, n ;
  struct S *p ;
  Entier N, D, d ;
  PipVector * vector ;
  
  mpz_init(N) ;
  mpz_init(D) ;
  mpz_init(d) ;

  vector = (PipVector *)malloc(sizeof(PipVector)) ;
  if (vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  p = sol_space + (*i) ;
  n = mpz_get_si(p->param1) ;
  vector->nb_elements = n ;
  vector->the_vector = (Entier *)malloc(sizeof(Entier)*n) ;
  if (vector->the_vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  vector->the_deno = (Entier *)malloc(sizeof(Entier)*n) ;
  if (vector->the_deno == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  
  for (j=0;j<n;j++)
  { (*i)++ ;
    p++ ;
    mpz_set(N,p->param1) ;
    mpz_set(D,p->param2) ;
    mpz_gcd(d, N, D);

    mpz_init(vector->the_vector[j]) ;
    mpz_divexact(vector->the_vector[j],N,d) ;

    mpz_init(vector->the_deno[j]) ;
    if (mpz_cmp(d, D) == 0)
    mpz_set(vector->the_deno[j],UN) ;
    else
    mpz_divexact(vector->the_deno[j],D,d) ;
  }
  (*i)++ ;

  mpz_clear(d);
  mpz_clear(D);
  mpz_clear(N);
  return(vector) ;
}


/* Fonction sol_newparm_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a un Newparm dans la grammaire dans une structure de type PipNewparm. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipNewparm contenant les informations de ce Newparm.
 * Premiere version : Ced. 18 octobre 2001. 
 * 24 octobre 2002 : Premiere version MP. 
 */
PipNewparm * sol_newparm_edit(int * i)
{ struct S * p ;
  PipNewparm * newparm, * newparm_new, * newparm_now ;

  /* On place p au lieu de lecture. */
  p = sol_space + (*i) ;
  /* On passe le New et le Div pour aller a Form et lire le VECTOR. */
  (*i) += 2 ;

  newparm = (PipNewparm *)malloc(sizeof(PipNewparm)) ;
  if (newparm == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  newparm->vector = sol_vector_edit(i) ;
  newparm->rank = mpz_get_si(p->param1) ;
  /* On met p a jour pour lire le denominateur (un Val de param2 UN). */
  p = sol_space + (*i) ;
  mpz_init_set(newparm->deno,p->param1) ;
  newparm->next = NULL ;

  newparm_now = newparm ;
  if (verbose)
  { fprintf(dump,"\n(newparm ") ;
    fprintf(dump,"%d",newparm->rank) ;
    fprintf(dump," (div ") ;
    pip_vector_print(dump,newparm->vector) ;
    fprintf(dump," ") ;
    mpz_out_str(dump,10,newparm->deno) ;
    fprintf(dump,"))") ;
  }
  
  /* On passe aux elements suivants. */
  (*i) ++ ;
  p = sol_space + (*i) ;
  while (p->flags == New)
  { (*i) += 2 ;
    newparm_new = (PipNewparm *)malloc(sizeof(PipNewparm)) ;
    if (newparm_new == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    newparm_new->vector = sol_vector_edit(i) ;
    newparm_new->rank = mpz_get_si(p->param1) ;
    p = sol_space + (*i) ;
    mpz_init_set(newparm_new->deno,p->param1) ;
    newparm_new->next = NULL ;
      
    newparm_now->next = newparm_new ;
    newparm_now = newparm_now->next ;
    if (verbose)
    { fprintf(dump,"\n(newparm ") ;
      fprintf(dump,"%d",newparm_new->rank) ;
      fprintf(dump," (div ") ;
      pip_vector_print(dump,newparm_new->vector) ;
      fprintf(dump," ") ;
      mpz_out_str(dump,10,newparm_new->deno) ;
      fprintf(dump,"))") ;
    }
    (*i) ++ ;
    p = sol_space + (*i) ;
  }
  return(newparm) ;
}


/* Fonction sol_list_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a une List dans la grammaire dans une structure de type PipList. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipList contenant les informations de cette List.
 * Premiere version : Ced. 18 octobre 2001. 
 */
PipList * sol_list_edit(int * i, int nb_elements)
{ PipList * list, * list_new, * list_now ;
  
  /* Pour le premier element. */
  list = (PipList *)malloc(sizeof(PipList)) ;
  if (list == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  list->vector = sol_vector_edit(i) ;
  list->next = NULL ;

  list_now = list ;
  if (verbose)
  { fprintf(dump,"\n(list ") ;
    pip_vector_print(dump,list->vector) ;
  }
  nb_elements-- ;

  /* Pour les elements suivants. */
  while (nb_elements--)
  { list_new = (PipList *)malloc(sizeof(PipList)) ;
    if (list_new == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    list_new->vector = sol_vector_edit(i) ;
    list_new->next = NULL ;
		    
    if (verbose)
    { fprintf(dump,"\n") ;
      pip_vector_print(dump,list_new->vector) ;
    }
    list_now->next = list_new ;
    list_now = list_now->next ;
  }
  if (verbose)
  fprintf(dump,"\n)") ;
  
  return(list) ;
}


/* Fonction sol_quast_edit :
 * Cette fonction a pour but de placer les informations de la solution
 * (qui sont contenues dans le tableau sol_space) dans une structure de
 * type PipQuast en vue d'une utilisation directe de la solution par une
 * application exterieure. Elle prend en parametre un pointeur vers une
 * case memoire contenant le numero de cellule du tableau sol_space
 * a partir de laquelle on doit commencer la lecture des informations. Elle
 * recoit aussi l'adresse du PipQuast qui l'a appelle (pour le champ father).
 * Elle retourne un pointeur vers une structure de type PipQuast qui
 * contient toutes les informations sur la solution (sous forme d'arbre).
 * Remarques : cette fonction lit les informations comme elles doivent
 * se presenter a la fin du traitement. Elle respecte scrupuleusement
 * la grammaire attendue et n'accepte de passer des cellules a Free
 * qu'entre une des trois grandes formes (if, list ou suite de newparm).
 * 20 juillet 2001 : Premiere version, Ced. 
 * 31 juillet 2001 : Ajout du traitement de l'option verbose = code*2 :0( 
 * 18 octobre 2001 : Grands changements dus a l'eclatement de la structure
 *                   PipVector en PipVector, PipNewparm et PipList, et
 *                   eclatement de la fonction avec sol_newparm_edit et
 *                   sol_list_edit.
 * 24 octobre 2002 : Premiere version MP. 
 */
PipQuast * sol_quast_edit(int * i, PipQuast * father)
{ int nb_elements ;
  struct S * p ;
  PipQuast * solution ;
  PipList * list_new, * list_now ;
  PipNewparm * newparm_new, * newparm_now ;
      
  /* On place p au lieu de lecture. */
  p = sol_space + (*i) ;
  /* En cas d'utilisation de l'option de simplification, une plage de
   * structures S peut avoir les flags a Free. On doit alors les passer.
   */
  while (p->flags == Free)
  { p ++ ;
    (*i) ++ ;
  }
  
  solution = (PipQuast *)malloc(sizeof(PipQuast)) ;
  if (solution == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  solution->newparm = NULL ;
  solution->list = NULL ;
  solution->condition = NULL ;
  solution->next_then = NULL ;
  solution->next_else = NULL ;
  solution->father = father ;
  
  /* On peut commencer par une chaine de nouveaux parametres... */
  if (p->flags == New)
  { solution->newparm = sol_newparm_edit(i) ;
    p = sol_space + (*i) ;
  }
  
  /* ...ensuite soit par une liste (vide ou non) soit par un if. */
  (*i)++ ; /* Factorise de List, Nil et If. */
  switch (p->flags)
  { case List : if ((nb_elements = mpz_get_si(p->param1)) != 0)
                solution->list = sol_list_edit(i,nb_elements) ;
		break ;
    case Nil  : if (verbose)
		fprintf(dump,"\n()") ;
                break ;
    case If   : solution->condition = sol_vector_edit(i) ;
                if (verbose)
		{ fprintf(dump,"\n(if ") ;
                  pip_vector_print(dump,solution->condition) ;
                }
		solution->next_then = sol_quast_edit(i,solution) ;
                solution->next_else = sol_quast_edit(i,solution) ;
                if (verbose)
		fprintf(dump,"\n)") ;
                break ;
    default   : fprintf(stderr,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                if (verbose)
		fprintf(dump,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                exit(1) ;
  }
  
  return(solution) ;
}

























