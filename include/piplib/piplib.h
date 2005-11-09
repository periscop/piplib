/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.h                                   *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
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
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

/* Premiere version du 18 septembre 2002. */

#if !defined(LINEAR_VALUE_IS_LONGLONG) && !defined(LINEAR_VALUE_IS_INT)
#if !defined(LINEAR_VALUE_IS_MP)
# error Please define LINEAR_VALUE_IS_* or #include polylib32.h or polylib64.h
#endif
#endif

#if defined(LINEAR_VALUE_IS_LONGLONG)
# define Entier   long long
# define FORMAT   "%lld"
# define VAL_UN   1LL
# define VAL_ZERO 0LL
#elif defined(LINEAR_VALUE_IS_INT) 
# define Entier   long int
# define FORMAT   "%ld"
# define VAL_UN   1L
# define VAL_ZERO 0L
#elif defined(LINEAR_VALUE_IS_MP) 
# include <gmp.h>
# define Entier   mpz_t
# define FORMAT   "%d"
# define GMP_INPUT_FORMAT   "%lZd"
#endif

#ifndef PIPLIB_H
#define PIPLIB_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

# include <piplib/type.h>
# include <piplib/sol.h>
# include <piplib/tab.h>
# include <piplib/funcall.h>


/* Structure PipMatrix :
 * Structure de matrice au format PolyLib. Le premier element d'une ligne
 * indique quand il vaut 1 que la ligne decrit une inequation de la forme
 * p(x)>=0 et quand il vaut 0, que la ligne decrit une egalite de la forme
 * p(x)=0. Le dernier element de chaque ligne correspond au coefficient
 * constant.
 */
struct pipmatrix
{ unsigned NbRows, NbColumns ;
  Entier **p ;
  Entier *p_Init ;
  int p_Init_size;	        /* Only for PolyLib compatibility under MP
                                 * version: PolyLib makes sometimes
				 * overestimates on the size of the matrices,
				 * in order to go faster. Thus
				 * NbRows*NbColumns is not the number of
				 * allocated elements. With MP version, we
				 * have to think to mpz_clear() all the
				 * initialized elements before freing, then
				 * we need to know the number of allocated
				 * elements: p_Init_size.
				 */
} ;
typedef struct pipmatrix PipMatrix ;


/* Structure PipVector :
 * Cette structure contient un Vector de 'nb_elements' la ieme composante de
 * ce vecteur vaut the_vector[i]/the_deno[i].
 */
struct pipvector
{ int nb_elements ;             /* Nombre d'elements du vecteur. */
  Entier * the_vector ;         /* Numerateurs du vecteur. */
  Entier * the_deno ;           /* Denominateurs du vecteur. */
} ;
typedef struct pipvector PipVector ;


/* Structure PipNewparm :
 * Liste chainee de Newparm, les informations d'un newparm etant son rang, un
 * vecteur de coefficients et un denominateur. Le newparm est egal a la division
 * du vecteur par le denominateur.
 */
struct pipnewparm
{ int rank ;                    /* Rang du 'newparm'. */
  PipVector * vector ;          /* Le vector decrivant le newparm. */
  Entier deno ;                 /* Denominateur du 'newparm'. */
  struct pipnewparm * next ;    /* Pointeur vers le newparm suivant. */
} ;
typedef struct pipnewparm PipNewparm ;


/* Structure PipList :
 * Liste chainee de Vector.
 */
struct piplist
{ PipVector * vector ;          /* Le vector contenant la partie de solution. */
  struct piplist * next ;       /* Pointeur vers l'element suivant. */
} ;
typedef struct piplist PipList ;


/* Structure pipquast :
 * Arbre binaire. Conformement a la grammaire de sortie (voir mode d'emploi), un
 * noeud de l'arbre des solutions debute par une liste de 'newparm'. Il continue
 * ensuite soit par une 'list' (alors condition vaut null), soit par un 'if'
 * (alors le champ condition contient la condition).
 */
struct pipquast
{ PipNewparm * newparm ;        /* Les 'newparm'. */
  PipList * list ;              /* La 'list' si pas de 'if'. */
  PipVector * condition ;       /* La condition si 'if'. */
  struct pipquast * next_then ; /* Noeud si condition et si verifiee. */
  struct pipquast * next_else ; /* Noeud si condition et si non verifiee. */
  struct pipquast * father ;    /* Pointeur vers le quast pere. */
} ;      
typedef struct pipquast PipQuast ;


/* Structure pipoptions:
 * This structure contains each option that can be set to change the PIP
 * behaviour.
 */
struct pipoptions
{ int Nq ;                      /* 1 if an integer solution is needed,
                                 * 0 otherwise.
				 */
  int Verbose ;                 /* -1 -> absolute silence,
                                 *  0 -> relative silence,
                                 *  1 -> information on cuts when an integer
				 *       solution is needed,
                                 *  2 -> information sur les pivots et les
				 *       déterminants,
                                 *  3 -> information on arrays,
                                 * Each option include the preceding.
				 */
  int Simplify ;                /* Set to 1 to eliminate some trivial
                                 * solutions, 0 otherwise.
				 */
  int Deepest_cut ;             /* Set to 1 to include deepest cut
                                 * algorithm.
				 */
  int Max ;                     /* Set to 0. */
} ;      
typedef struct pipoptions PipOptions ;


/* Prototypes des fonctions d'affichages des structures de la PipLib. */
void pip_matrix_print(FILE *, PipMatrix *) ;
void pip_vector_print(FILE *, PipVector *) ;
void pip_newparm_print(FILE * foo, PipNewparm *, int indent) ;
void pip_list_print(FILE * foo, PipList *, int indent) ;
void pip_quast_print(FILE *, PipQuast *, int) ;


/* Prototypes des fonctions de liberation memoire des structures de la PipLib.*/
void pip_matrix_free(PipMatrix *) ;
void pip_vector_free(PipVector *) ;
void pip_newparm_free(PipNewparm *) ;
void pip_list_free(PipList *) ;
void pip_quast_free(PipQuast *) ;
void pip_options_free(PipOptions *) ;


/* Prototypes des fonctions d'acquisition de matrices de contraintes et
 * options.
 */
PipMatrix * pip_matrix_alloc(unsigned, unsigned) ;
PipMatrix * pip_matrix_read(FILE *) ;
PipOptions * pip_options_init(void) ;
 

/* initialization of pip library */
void pip_init();
void pip_close();


/* Prototype de la fonction de resolution :
 * pip_solve resoud le probleme qu'on lui passe en parametre, suivant les
 * options elles aussi en parametre. Elle renvoie la solution sous forme
 * d'un arbre de PipQuast. Parametres :
 * - probleme :
 * 1 PipMatrix  : systeme des inequations definissant le domaine des inconnues,
 * 2 PipMatrix  : systeme des inequations satisfaites par les parametres,
 * 3 int        : column rank of the bignum, or negative value if there
 *                is no big parameter.
 * 4 PipOptions : options for PIP.
 */ 
PipQuast * pip_solve(PipMatrix *, PipMatrix *, int, PipOptions *) ;

/* Ced : ajouts specifiques a la PipLib pour funcall. */
Tableau * tab_Matrix2Tableau(PipMatrix *, int, int, int, int, int) ;
Tableau * tab_Matrix2TableauMax(PipMatrix *, int, int, int, int) ;
PipQuast * sol_quast_edit(int *, PipQuast *) ;

#if defined(__cplusplus)
  }
#endif 
#endif /* define PIPLIB_H */
