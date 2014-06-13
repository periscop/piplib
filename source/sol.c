/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                   sol.h                                    *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This library is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at    *
 * your option) any later version.                                            *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Paul Feautrier and Cedric Bastoul                               *
 *                                                                            *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pip.h"

extern int PIPLIB_NAME(verbose);
extern FILE *PIPLIB_NAME(dump);

struct PIPLIB_NAME(S)
    {int flags;
     PIPLIB_NAME(piplib_int_t) param1, param2;
    };

#define Free 0
#define Nil  1
#define If   2
#define List 3
#define Form 4
#define New  5
#define Div  6
#define Val  7
#define Error 8

struct PIPLIB_NAME(S) * PIPLIB_NAME(sol_space);
static int PIPLIB_NAME(sol_free);

void PIPLIB_NAME(sol_init)(void)
{
 PIPLIB_NAME(sol_free) = 0;
 PIPLIB_NAME(sol_space) = (struct PIPLIB_NAME(S) *)malloc(
                            SOL_SIZE*sizeof(struct PIPLIB_NAME(S))) ;
}

void PIPLIB_NAME(sol_close)(void)
{
 free(PIPLIB_NAME(sol_space)) ;
}

int PIPLIB_NAME(sol_hwm)()
{
 return(PIPLIB_NAME(sol_free));
}

void PIPLIB_NAME(sol_reset)(p)
int p;
{int i;
 if(p<0 || p>=SOL_SIZE)
     {fprintf(stderr,
              "Syserr : PIPLIB_NAME(sol_reset) : Memory allocation error\n");
      exit(40);
     }
 for(i=p; i<PIPLIB_NAME(sol_free); i++){
   piplib_int_clear(PIPLIB_NAME(sol_space)[i].param1);
   piplib_int_clear(PIPLIB_NAME(sol_space)[i].param2);
 }
 PIPLIB_NAME(sol_free) = p;
}

struct PIPLIB_NAME(S) *PIPLIB_NAME(sol_alloc)(void)
{struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_space) + PIPLIB_NAME(sol_free);
 r->flags = Free;
 piplib_int_init_set_si(r->param1,0);
 piplib_int_init_set_si(r->param2,0);
 PIPLIB_NAME(sol_free)++;
 if(PIPLIB_NAME(sol_free) >= SOL_SIZE)
     {fprintf(stderr, "The solution is too complex! : sol\n");
      exit(26);
     }
     return(r);
}

void PIPLIB_NAME(sol_nil)(void)
{
 struct PIPLIB_NAME(S) * r;
 r = PIPLIB_NAME(sol_alloc)();
 r -> flags = Nil;
 if(PIPLIB_NAME(verbose) > 0)
   {fprintf(PIPLIB_NAME(dump), "\nNil");
    fflush(PIPLIB_NAME(dump));
  }
}

void PIPLIB_NAME(sol_error)(int c)
{
 struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_alloc)();
 r->flags = Nil;
 piplib_int_set_si(r->param1, c);
 if(PIPLIB_NAME(verbose) > 0) {
     fprintf(PIPLIB_NAME(dump), "Erreur %d\n", c);
     fflush(PIPLIB_NAME(dump));
     }
}

int PIPLIB_NAME(is_not_Nil)(p)
int p;
{
 return(PIPLIB_NAME(sol_space)[p].flags != Nil);
}

void PIPLIB_NAME(sol_if)(void)
{
 struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_alloc)();
 r -> flags = If;
 if(PIPLIB_NAME(verbose) > 0) {
     fprintf(PIPLIB_NAME(dump), "\nIf ");
     fflush(PIPLIB_NAME(dump));
   }
}

void PIPLIB_NAME(sol_list)(n)
int n;
{struct PIPLIB_NAME(S) * r;
 r = PIPLIB_NAME(sol_alloc)();
 r->flags = List;
 piplib_int_set_si(r->param1, n);
 if(PIPLIB_NAME(verbose) > 0) {
     fprintf(PIPLIB_NAME(dump), "\nList %d ", n);
     fflush(PIPLIB_NAME(dump));
}
}

void PIPLIB_NAME(sol_forme)(l)
int l;
{
 struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_alloc)();
 r -> flags = Form;
 piplib_int_set_si(r -> param1, l);
 if(PIPLIB_NAME(verbose) > 0) {
     fprintf(PIPLIB_NAME(dump), "\nForme %d ", l);
     fflush(PIPLIB_NAME(dump));
   }
}

void PIPLIB_NAME(sol_new)(k)
int k;
{
 struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_alloc)();
 r -> flags = New;
 piplib_int_set_si(r -> param1, k);
 if(PIPLIB_NAME(verbose) > 0) {
     fprintf(PIPLIB_NAME(dump), "New %d ", k);
     fflush(PIPLIB_NAME(dump));
   }
}

void PIPLIB_NAME(sol_div)()
{
 struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_alloc)();
 r -> flags = Div;
 if(PIPLIB_NAME(verbose) > 0) {
     fprintf(PIPLIB_NAME(dump), "Div ");
     fflush(PIPLIB_NAME(dump));
   }
}

void PIPLIB_NAME(sol_val)(n, d)
PIPLIB_NAME(piplib_int_t) n, d;
{
 struct PIPLIB_NAME(S) *r;
 r = PIPLIB_NAME(sol_alloc)();
 r -> flags = Val;
 piplib_int_assign(r->param1, n);
 piplib_int_assign(r->param2, d);
 if(PIPLIB_NAME(verbose) > 0) {
   fprintf(PIPLIB_NAME(dump), "val(");
   piplib_int_print(PIPLIB_NAME(dump), n);
   fprintf(PIPLIB_NAME(dump), "/");
   piplib_int_print(PIPLIB_NAME(dump), d);
   fprintf(PIPLIB_NAME(dump), ") ");
   fflush(PIPLIB_NAME(dump));
  }
}

void PIPLIB_NAME(sol_val_one)(PIPLIB_NAME(piplib_int_t) n) {
  PIPLIB_NAME(piplib_int_t) one;
  piplib_int_init_set_si(one, 1);
  PIPLIB_NAME(sol_val)(n, one);
  piplib_int_clear(one);
}

void PIPLIB_NAME(sol_val_zero_one)() {
  PIPLIB_NAME(piplib_int_t) zero;
  PIPLIB_NAME(piplib_int_t) one;
  piplib_int_init_set_si(zero, 0);
  piplib_int_init_set_si(one, 1);
  PIPLIB_NAME(sol_val)(zero, one);
  piplib_int_clear(one);
  piplib_int_clear(zero);
}

int PIPLIB_NAME(skip)(int);

/* a` partir d'un point de la solution, sauter un objet
bien forme' ainsi qu'un e'ventuel New et pointer sur l'objet
suivant */

int PIPLIB_NAME(skip_New) (int i)
{
 if(PIPLIB_NAME(sol_space)[i].flags != New) return i;
 i = PIPLIB_NAME(skip)(i+1);      /* sauter le Div */
 return i;
}
/* au lancement, i indexe une cellule qui est la te^te d'un objet.
   la valeur retourne'e est la te^te de l'objet qui suit. Les
   objets de type New sont e'limine's                                */

int PIPLIB_NAME(skip) (int i)
{int n, f;
 while((f = PIPLIB_NAME(sol_space)[i].flags) == Free || f == Error) i++;
 switch (PIPLIB_NAME(sol_space)[i].flags) {
 case Nil : case Val : i++; break;
 case New : i = PIPLIB_NAME(skip_New)(i); break;
 case If : i = PIPLIB_NAME(skip)(i+1);        /* sauter le pre'dicat */
	   i = PIPLIB_NAME(skip)(i);          /* sauter le vrai */
	   i = PIPLIB_NAME(skip)(i); break;   /* sauter le faux */
 case List : case Form :
           n = piplib_int_get_si(PIPLIB_NAME(sol_space)[i].param1);
	   i++;
	   while(n--) i = PIPLIB_NAME(skip)(i);
	   break;
 case Div : i = PIPLIB_NAME(skip)(i+1);       /* sauter la forme */
	    i = PIPLIB_NAME(skip)(i);         /* sauter le diviseur */
	    break;
 default : fprintf(stderr,
                   "Syserr : PIPLIB_NAME(skip) : unknown %d\n",
                   PIPLIB_NAME(sol_space)[i].flags);
 }
 return PIPLIB_NAME(skip_New)(i);
}
/* simplification de la solution : e'limination des constructions
   (if p () ()). N'est en service qu'en pre'sence de l'option -z */

void PIPLIB_NAME(sol_simplify)(int i)
{int j, k, l;
 if(PIPLIB_NAME(sol_space)[i].flags == If) {
     j = PIPLIB_NAME(skip)(i+1);        /* j : debut de la partie vraie */
     k = PIPLIB_NAME(skip)(j);          /* k : debut de la partie fausse */
     PIPLIB_NAME(sol_simplify)(k);
     PIPLIB_NAME(sol_simplify)(j);
     if (PIPLIB_NAME(sol_space)[j].flags == Nil &&
         PIPLIB_NAME(sol_space)[k].flags == Nil) {
	 PIPLIB_NAME(sol_space)[i].flags = Nil;
	 if (k >= PIPLIB_NAME(sol_free) - 1) 
	    PIPLIB_NAME(sol_reset)(i+1);
	 else for(l = i+1; l<=k; l++) PIPLIB_NAME(sol_space)[l].flags = Free;
       }
   }

}
/* e'dition de la solution */

int PIPLIB_NAME(sol_edit)(FILE *foo, int i)
{int j, n;
 struct PIPLIB_NAME(S) *p;
 PIPLIB_NAME(piplib_int_t) N, D, d;
 piplib_int_init(N);
 piplib_int_init(D);
 piplib_int_init(d);
 
 p = PIPLIB_NAME(sol_space) + i;
 for(;;) {
   if(p->flags == Free) {
     p++;
     i++;
     continue;
   }
   if(p->flags == New) {
     n = piplib_int_get_si(p->param1);
     fprintf(foo, "(newparm %d ", n);
     if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "(newparm %d ", n);
     i = PIPLIB_NAME(sol_edit)(foo, ++i);
     p = PIPLIB_NAME(sol_space) +i;
     fprintf(foo, ")\n");
     if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), ")\n");
     continue;
   }
   break;
 }
 switch(p->flags){
 case Nil : fprintf(foo, "()\n");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "()\n");
   i++; break;
 case Error :
   fprintf(foo, "Error %d\n", piplib_int_get_si(p->param1));
   if(PIPLIB_NAME(verbose)>0)
   fprintf(PIPLIB_NAME(dump), "Error %d\n", piplib_int_get_si(p->param1));
   i++; break;
 case If  : fprintf(foo, "(if ");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "(if ");
   i = PIPLIB_NAME(sol_edit)(foo, ++i);
   i = PIPLIB_NAME(sol_edit)(foo, i);
   i = PIPLIB_NAME(sol_edit)(foo, i);
   fprintf(foo, ")\n");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), ")\n");
   break;
 case List: fprintf(foo, "(list ");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "(list ");
   n = piplib_int_get_si(p->param1);
   i++;
   while(n--) i = PIPLIB_NAME(sol_edit)(foo, i);
   fprintf(foo, ")\n");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), ")\n");
   break;
 case Form: fprintf(foo, "#[");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "#[");
   n = piplib_int_get_si(p->param1);
   for(j = 0; j<n; j++){
     i++; p++;
     piplib_int_assign(N, p->param1);
     piplib_int_assign(D, p->param2);
     piplib_int_gcd(d, N, D);
     if(piplib_int_eq(d, D)) {
       putc(' ', foo);
       piplib_int_div_exact(N, N, d);
       piplib_int_print(foo, N);
       if(PIPLIB_NAME(verbose)>0){
         putc(' ', PIPLIB_NAME(dump));
         piplib_int_print(PIPLIB_NAME(dump), N);
       }
     }
     else{
       piplib_int_div_exact(N, N, d);
       piplib_int_div_exact(D, D, d);
       putc(' ', foo);
       piplib_int_print(foo, N);
       putc('/', foo);
       piplib_int_print(foo, D);
       if(PIPLIB_NAME(verbose)>0){
         putc(' ', PIPLIB_NAME(dump));
         piplib_int_print(PIPLIB_NAME(dump), N);
         putc('/', PIPLIB_NAME(dump));
         piplib_int_print(PIPLIB_NAME(dump), D);
       }
     }
   }
   fprintf(foo, "]\n");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "]\n");
   i++;
   break;
 case Div : fprintf(foo, "(div ");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "(div ");
   i = PIPLIB_NAME(sol_edit)(foo, ++i);
   i = PIPLIB_NAME(sol_edit)(foo, i);
   fprintf(foo, ")\n");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), ")\n");
   break;
 case Val :
   piplib_int_assign(N, p->param1);
   piplib_int_assign(D, p->param2);
   piplib_int_gcd(d, N, D);
   if (piplib_int_eq(d, D)) {
     piplib_int_div_exact(N, N, d);
     putc(' ', foo);
     piplib_int_print(foo, N);
     if(PIPLIB_NAME(verbose)>0){
       putc(' ', PIPLIB_NAME(dump));
       piplib_int_print(PIPLIB_NAME(dump), N);
     }
   }
   else{
     piplib_int_div_exact(N, N, d);
     piplib_int_div_exact(D, D, d);
     putc(' ', foo);
     piplib_int_print(foo, N);
     fprintf(foo, "/");
     piplib_int_print(foo, D);
     if(PIPLIB_NAME(verbose)>0){
       putc(' ', PIPLIB_NAME(dump));
       piplib_int_print(PIPLIB_NAME(dump), N);
       fprintf(PIPLIB_NAME(dump), "/");
       piplib_int_print(PIPLIB_NAME(dump), D);
     }
   }
   i++;
   break;
 default  : fprintf(foo, "Inconnu : sol\n");
   if(PIPLIB_NAME(verbose)>0)fprintf(PIPLIB_NAME(dump), "Inconnu : sol\n");
 }
 piplib_int_clear(d);
 piplib_int_clear(D);
 piplib_int_clear(N);
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
 */
PIPLIB_NAME(PipVector) * PIPLIB_NAME(sol_vector_edit)(int *i, int Bg,
                                                      int Urs_p, int flags)
{ int j, k, n, unbounded  = 0, first_urs;
  struct PIPLIB_NAME(S) *p ;
  PIPLIB_NAME(piplib_int_t) N, D, d ;
  PIPLIB_NAME(PipVector) * vector ;

  piplib_int_init(N);
  piplib_int_init(D);
  piplib_int_init(d);
  
  vector = (PIPLIB_NAME(PipVector) *)malloc(sizeof(PIPLIB_NAME(PipVector))) ;
  if (vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  p = PIPLIB_NAME(sol_space) + (*i) ;
  n = piplib_int_get_si(p->param1);
  if (flags & SOL_REMOVE)
    --n;
  n -= Urs_p;
  first_urs = Urs_p + (Bg >= 0);
  vector->nb_elements = n ;
  vector->the_vector = (PIPLIB_NAME(piplib_int_t) *)malloc(
                         sizeof(PIPLIB_NAME(piplib_int_t))*n) ;
  if (vector->the_vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  vector->the_deno = (PIPLIB_NAME(piplib_int_t) *)malloc(
                       sizeof(PIPLIB_NAME(piplib_int_t))*n) ;
  if (vector->the_deno == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  
  for (j=0, k=0; k < n; j++) {
    (*i)++ ;
    p++ ;

    piplib_int_assign(N, p->param1);
    piplib_int_assign(D, p->param2);
    piplib_int_gcd(d, N, D);

    if ((flags & SOL_SHIFT) && j == Bg) {
      piplib_int_sub(N, N, D);   /* subtract 1 */
      if (piplib_int_zero(N) == 0)
	unbounded = 1;
    }

    if ((flags & SOL_REMOVE) && j == Bg)
      continue;

    if (first_urs <= j && j < first_urs+Urs_p)
      continue;

    piplib_int_init(vector->the_vector[k]);
    piplib_int_div_exact(vector->the_vector[k], N, d);
    if (flags & SOL_NEGATE)
      piplib_int_oppose(vector->the_vector[k], vector->the_vector[k]);
    piplib_int_init(vector->the_deno[k]);
    if (piplib_int_eq(d, D))
      piplib_int_set_si(vector->the_deno[k], 1);
    else
      piplib_int_div_exact(vector->the_deno[k], D, d);
    ++k;
  }
  if (unbounded)
    for (k=0; k < n; k++)
      piplib_int_set_si(vector->the_deno[k], 0);
  (*i)++ ;

  piplib_int_clear(d);
  piplib_int_clear(D);
  piplib_int_clear(N);

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
 */
PIPLIB_NAME(PipNewparm) * PIPLIB_NAME(sol_newparm_edit)(int *i, int Bg,
                                                        int Urs_p, int flags)
{ struct PIPLIB_NAME(S) * p ;
  PIPLIB_NAME(PipNewparm) * newparm,
                          * newparm_first = NULL,
                          * newparm_now = NULL;

  /* On place p au lieu de lecture. */
  p = PIPLIB_NAME(sol_space) + (*i) ;

  do {
    /* On passe le New et le Div pour aller a Form et lire le VECTOR. */
    (*i) += 2 ;

    newparm = (PIPLIB_NAME(PipNewparm) *)malloc(
                sizeof(PIPLIB_NAME(PipNewparm)));
    if (newparm == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    newparm->vector = PIPLIB_NAME(sol_vector_edit)(i, Bg, Urs_p, flags);
    newparm->rank = piplib_int_get_si(p->param1);
    /* On met p a jour pour lire le denominateur (un Val de param2 1). */
    p = PIPLIB_NAME(sol_space) + (*i) ;
    piplib_int_init(newparm->deno);
    piplib_int_assign(newparm->deno, p->param1);
    if (flags & SOL_REMOVE)
      newparm->rank--;
    newparm->rank -= Urs_p;
    newparm->next = NULL ;

    if (newparm_now)
      newparm_now->next = newparm;
    else
      newparm_first = newparm;
    newparm_now = newparm ;
    if (PIPLIB_NAME(verbose) > 0)
    { fprintf(PIPLIB_NAME(dump),"\n(newparm ") ;
      fprintf(PIPLIB_NAME(dump), "%i", newparm->rank) ;
      fprintf(PIPLIB_NAME(dump)," (div ") ;
      PIPLIB_NAME(pip_vector_print)(PIPLIB_NAME(dump),newparm->vector) ;
      fprintf(PIPLIB_NAME(dump)," ") ;
      piplib_int_print(PIPLIB_NAME(dump), newparm->deno);
      fprintf(PIPLIB_NAME(dump),"))") ;
    }
  
    /* On passe aux elements suivants. */
    (*i) ++ ;
    p = PIPLIB_NAME(sol_space) + (*i) ;
  } while (p->flags == New);

  return newparm_first;
}


/* Fonction sol_list_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a une List dans la grammaire dans une structure de type PipList. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau sol_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipList contenant les informations de cette List.
 * Premiere version : Ced. 18 octobre 2001. 
 * 16 novembre 2005 : Ced. Prise en compte du cas 0 éléments, avant impossible.
 */
PIPLIB_NAME(PipList) * PIPLIB_NAME(sol_list_edit)(int *i, int nb_elements,
                                                  int Bg, int Urs_p, int flags)
{ PIPLIB_NAME(PipList) * list, * list_new, * list_now ;
  
  /* Pour le premier element. */
  list = (PIPLIB_NAME(PipList) *)malloc(sizeof(PIPLIB_NAME(PipList))) ;
  if (list == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  list->next = NULL ;
  
  if (nb_elements == 0)
  { list->vector = NULL ;
    return(list) ;
  }
  
  list->vector = PIPLIB_NAME(sol_vector_edit)(i, Bg, Urs_p, flags);

  list_now = list ;
  if (PIPLIB_NAME(verbose) > 0)
  { fprintf(PIPLIB_NAME(dump),"\n(list ") ;
    PIPLIB_NAME(pip_vector_print)(PIPLIB_NAME(dump),list->vector) ;
  }
  nb_elements-- ;

  /* Pour les elements suivants. */
  while (nb_elements--)
  { list_new = (PIPLIB_NAME(PipList) *)malloc(sizeof(PIPLIB_NAME(PipList))) ;
    if (list_new == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    list_new->vector = PIPLIB_NAME(sol_vector_edit)(i, Bg, Urs_p, flags);
    list_new->next = NULL ;
		    
    if (PIPLIB_NAME(verbose) > 0)
    { fprintf(PIPLIB_NAME(dump),"\n") ;
      PIPLIB_NAME(pip_vector_print)(PIPLIB_NAME(dump),list_new->vector) ;
    }
    list_now->next = list_new ;
    list_now = list_now->next ;
  }
  if (PIPLIB_NAME(verbose) > 0)
  fprintf(PIPLIB_NAME(dump),"\n)") ;
  
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
 * 20  juillet 2001 : Premiere version, Ced. 
 * 31  juillet 2001 : Ajout du traitement de l'option verbose = code*2 :0( 
 * 18  octobre 2001 : Grands changements dus a l'eclatement de la structure
 *                    PipVector en PipVector, PipNewparm et PipList, et
 *                    eclatement de la fonction avec sol_newparm_edit et
 *                    sol_list_edit.
 * 16 novembre 2005 : (debug) Même si une liste est vide il faut la créer pour
 *                    afficher plus tard le (list), repéré par Sven Verdoolaege.
 */
PIPLIB_NAME(PipQuast) *PIPLIB_NAME(sol_quast_edit)(
  int *i, PIPLIB_NAME(PipQuast) *father, int Bg, int Urs_p, int flags)
{ int nb_elements ;
  struct PIPLIB_NAME(S) * p ;
  PIPLIB_NAME(PipQuast) * solution ;
    
  /* On place p au lieu de lecture. */
  p = PIPLIB_NAME(sol_space) + (*i) ;
  /* En cas d'utilisation de l'option de simplification, une plage de
   * structures S peut avoir les flags a Free. On doit alors les passer.
   */
  while (p->flags == Free)
  { p ++ ;
    (*i) ++ ;
  }
  
  solution = (PIPLIB_NAME(PipQuast) *)malloc(sizeof(PIPLIB_NAME(PipQuast))) ;
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
  { solution->newparm = PIPLIB_NAME(sol_newparm_edit)(i, Bg, Urs_p,
                                                      flags & SOL_REMOVE);
    p = PIPLIB_NAME(sol_space) + (*i) ;
  }
  
  /* ...ensuite soit par une liste (vide ou non) soit par un if. */
  (*i)++ ; /* Factorise de List, Nil et If. */
  switch (p->flags)
  { case List :
                nb_elements = piplib_int_get_si(p->param1) ;
                solution->list = PIPLIB_NAME(sol_list_edit)(i, nb_elements,
                                                            Bg, Urs_p, flags);
		if (flags & SOL_DUAL)
		    solution->next_then = PIPLIB_NAME(sol_quast_edit)(i, solution,
                                                              Bg, Urs_p, 0);
		break ;
    case Nil  : if (PIPLIB_NAME(verbose) > 0)
		fprintf(PIPLIB_NAME(dump),"\n()") ;
                break ;
    case If   : solution->condition = 
			    PIPLIB_NAME(sol_vector_edit)(i, Bg, Urs_p, flags & SOL_REMOVE);
                if (PIPLIB_NAME(verbose) > 0)
		{ fprintf(PIPLIB_NAME(dump),"\n(if ") ;
                  PIPLIB_NAME(pip_vector_print)(
                    PIPLIB_NAME(dump),solution->condition) ;
                }
		solution->next_then = PIPLIB_NAME(sol_quast_edit)(i, solution,
		                                                  Bg, Urs_p, flags);
		solution->next_else = PIPLIB_NAME(sol_quast_edit)(i, solution,
		                                                  Bg, Urs_p, flags);
                if (PIPLIB_NAME(verbose) > 0)
		fprintf(PIPLIB_NAME(dump),"\n)") ;
                break ;
    default   : fprintf(stderr,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                if (PIPLIB_NAME(verbose) > 0)
		fprintf(PIPLIB_NAME(dump),"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                exit(1) ;
  }
  
  return(solution) ;
}
