/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.h                                   *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
 *                                                                            *
 * This library is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at    *
 * your option) any later version.                                            *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

/* Premiere version du 18 septembre 2002. */


#if !defined(PIPLIB_INT_SP) && !defined(PIPLIB_INT_DP) && !defined(PIPLIB_INT_GMP) && !defined(PIPLIB_INT_OSL)
  #error "Define of PIPLIB_INT_SP or PIPLIB_INT_DP or PIPLIB_INT_GMP or \
          PIPLIB_INT_OSL not found"
#endif


#include <stdio.h>
#include <math.h>


#ifdef PIPLIB_INT_GMP

  #include <gmp.h>

  #define PIPLIB_NAME(name) name##_gmp

  #define PIPLIB_ONE_DETERMINANT

  typedef mpz_t PIPLIB_NAME(piplib_int_t);
  #define piplib_int_format "%d"

  #define piplib_int_init(i) (mpz_init(i))
  #define piplib_int_init_set(i, v) (mpz_init_set(i, v))
  #define piplib_int_init_set_si(i, v) (mpz_init_set_si(i, v))
  #define piplib_int_assign(r, i) (mpz_set(r, i))
  #define piplib_int_set_si(r, i) (mpz_set_si(r, i))
  #define piplib_int_clear(i) (mpz_clear(i))
  #define piplib_int_print(file, i) (mpz_out_str(file, 10, i))
  #define piplib_int_sscanf(string, i) (gmp_sscanf(string, "%lZd", i))

  #define piplib_int_get_si(i) ((int)mpz_get_si(i))
  #define piplib_int_get_d(i) (mpz_get_d(i))

  #define piplib_int_add(r, a, b) (mpz_add(r, a, b))
  #define piplib_int_sub(r, a, b) (mpz_sub(r, a, b))
  #define piplib_int_increment(r, i) (mpz_add_ui(r, i, 1))
  #define piplib_int_decrement(r, i) (mpz_sub_ui(r, i, 1))
  #define piplib_int_mul(r, a, b) (mpz_mul(r, a, b))
  #define piplib_int_div_exact(q, a, b) (mpz_divexact(q, a, b))
  #define piplib_int_floor_div_q(q, a, b) (mpz_fdiv_q(q, a, b))
  #define piplib_int_floor_div_r(r, a, b) (mpz_fdiv_r(r, a, b))
  #define piplib_int_floor_div_q_r(q, r, a, b) (mpz_fdiv_qr(q, r, a, b))
  #define piplib_int_mod(mod, a, b) (mpz_mod(mod, a, b))
  #define piplib_int_gcd(gcd, a, b) (mpz_gcd(gcd, a, b))
  #define piplib_int_oppose(r, i) (mpz_neg(r, i))
  #define piplib_int_size_in_base_2(i) (mpz_sizeinbase(i, 2))
  #define piplib_int_size_in_base_10(i) (mpz_sizeinbase(i, 10))

  #define piplib_int_eq(a, b) (mpz_cmp(a, b) == 0)
  #define piplib_int_ne(a, b) (mpz_cmp(a, b) != 0)
  #define piplib_int_zero(i) (mpz_sgn(i) == 0)
  #define piplib_int_one(i) (mpz_cmp_si(i, 1) == 0)
  #define piplib_int_pos(i) (mpz_sgn(i) > 0)
  #define piplib_int_neg(i) (mpz_sgn(i) < 0)

#endif

#ifdef PIPLIB_INT_OSL

  #include <osl/int.h>

  #define PIPLIB_NAME(name) name##_osl

  #define PIPLIB_ONE_DETERMINANT

  extern int PIPLIB_INT_OSL_PRECISION;

  typedef osl_int_t PIPLIB_NAME(piplib_int_t);
  #define piplib_int_format "%d"

  #define piplib_int_init(i) (osl_int_init(PIPLIB_INT_OSL_PRECISION, &(i)))
  #define piplib_int_init_set(i, v) (osl_int_init_set(PIPLIB_INT_OSL_PRECISION, &(i), v))
  #define piplib_int_init_set_si(i, v) (osl_int_init_set_si(PIPLIB_INT_OSL_PRECISION, &(i), v))
  #define piplib_int_assign(r, i) (osl_int_assign(PIPLIB_INT_OSL_PRECISION, &(r), i))
  #define piplib_int_set_si(r, i) (osl_int_set_si(PIPLIB_INT_OSL_PRECISION, &(r), i))
  #define piplib_int_clear(i) (osl_int_clear(PIPLIB_INT_OSL_PRECISION, &(i)))
  #define piplib_int_print(file, i) (osl_int_print(file, PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_sscanf(string, i) (osl_int_sscanf(string, PIPLIB_INT_OSL_PRECISION, &(i)))

  #define piplib_int_get_si(i) (osl_int_get_si(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_get_d(i) (osl_int_get_d(PIPLIB_INT_OSL_PRECISION, i))

  #define piplib_int_add(r, a, b) (osl_int_add(PIPLIB_INT_OSL_PRECISION, &(r), a, b))
  #define piplib_int_sub(r, a, b) (osl_int_sub(PIPLIB_INT_OSL_PRECISION, &(r), a, b))
  #define piplib_int_increment(r, i) (osl_int_increment(PIPLIB_INT_OSL_PRECISION, &(r), i))
  #define piplib_int_decrement(r, i) (osl_int_decrement(PIPLIB_INT_OSL_PRECISION, &(r), i))
  #define piplib_int_mul(r, a, b) (osl_int_mul(PIPLIB_INT_OSL_PRECISION, &(r), a, b))
  #define piplib_int_div_exact(q, a, b) (osl_int_div_exact(PIPLIB_INT_OSL_PRECISION, &(q), a, b))
  #define piplib_int_floor_div_q(q, a, b) (osl_int_floor_div_q(PIPLIB_INT_OSL_PRECISION, &(q), a, b))
  #define piplib_int_floor_div_r(r, a, b) (osl_int_floor_div_r(PIPLIB_INT_OSL_PRECISION, &(r), a, b))
  #define piplib_int_floor_div_q_r(q, r, a, b) (osl_int_floor_div_q_r(PIPLIB_INT_OSL_PRECISION, &(q), &(r), a, b))
  #define piplib_int_mod(mod, a, b) (osl_int_mod(PIPLIB_INT_OSL_PRECISION, &(mod), a, b))
  #define piplib_int_gcd(gcd, a, b) (osl_int_gcd(PIPLIB_INT_OSL_PRECISION, &(gcd), a, b))
  #define piplib_int_oppose(r, i) (osl_int_oppose(PIPLIB_INT_OSL_PRECISION, &(r), i))
  #define piplib_int_size_in_base_2(i) (osl_int_size_in_base_2(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_size_in_base_10(i) (osl_int_size_in_base_10(PIPLIB_INT_OSL_PRECISION, i))

  #define piplib_int_eq(a, b) (osl_int_eq(PIPLIB_INT_OSL_PRECISION, a, b))
  #define piplib_int_ne(a, b) (osl_int_ne(PIPLIB_INT_OSL_PRECISION, a, b))
  #define piplib_int_zero(i) (osl_int_zero(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_one(i) (osl_int_one(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_pos(i) (osl_int_pos(PIPLIB_INT_OSL_PRECISION, i))
  #define piplib_int_neg(i) (osl_int_neg(PIPLIB_INT_OSL_PRECISION, i))

#endif

#ifdef PIPLIB_INT_SP

  #define PIPLIB_NAME(name) name##_sp

  typedef long int PIPLIB_NAME(piplib_int_t);
  #define piplib_int_format "%ld"

#endif

#ifdef PIPLIB_INT_DP

  #define PIPLIB_NAME(name) name##_dp

  typedef long long int PIPLIB_NAME(piplib_int_t);
  #define piplib_int_format "%lld"

#endif


#if defined(PIPLIB_INT_SP) || defined(PIPLIB_INT_DP)

  // Copy from osl_int
  long long int PIPLIB_NAME(piplib_llgcd)(long long int const, long long int const);
  long long int PIPLIB_NAME(piplib_llgcd_llabs)(long long int const, long long int const);
  size_t PIPLIB_NAME(piplib_lllog2)(long long int);
  size_t PIPLIB_NAME(piplib_lllog10)(long long int);
  long long int PIPLIB_NAME(piplib_llmod)(long long int const, long long int const);

  #define piplib_int_init(i) (i = 0)
  #define piplib_int_init_set(i, v) (i = v)
  #define piplib_int_init_set_si(i, v) (i = v)
  #define piplib_int_assign(r, i) (r = i)
  #define piplib_int_set_si(r, i) (r = (PIPLIB_NAME(piplib_int_t))i)
  #define piplib_int_clear(i) do { } while (0)
  #define piplib_int_print(file, i) (fprintf(file, piplib_int_format, i))
  #define piplib_int_sscanf(string, i) (sscanf(string, piplib_int_format, &i))

  #define piplib_int_get_si(i) ((int)(i))
  #define piplib_int_get_d(i) ((double)(i))

  #define piplib_int_add(r, a, b) (r = a + b)
  #define piplib_int_sub(r, a, b) (r = a - b)
  #define piplib_int_increment(r, i) (r = i + 1)
  #define piplib_int_decrement(r, i) (r = i - 1)
  #define piplib_int_mul(r, a, b) (r = a * b)
  #define piplib_int_div_exact(q, a, b) (q = (a) / (b))
  //#define piplib_int_floor_div_q(q, a, b) (q = (PIPLIB_NAME(piplib_int_t))(piplib_ll_floor_div_q(a, b)))
  #define piplib_int_floor_div_q(q, a, b) (q = (PIPLIB_NAME(piplib_int_t))((a - PIPLIB_NAME(piplib_llmod)(a, b)) / (b)))
  #define piplib_int_floor_div_r(r, a, b) (r = (PIPLIB_NAME(piplib_int_t))PIPLIB_NAME(piplib_llmod)(a, b))
  //#define piplib_int_floor_div_r(r, a, b) (r = (PIPLIB_NAME(piplib_int_t))(piplib_ll_floor_div_r(a, b)))
  #define piplib_int_floor_div_q_r(q, r, a, b) do { piplib_int_floor_div_q(q, a, b); piplib_int_floor_div_r(r, a, b); } while (0)
  #define piplib_int_mod(mod, a, b) (mod = (PIPLIB_NAME(piplib_int_t))(PIPLIB_NAME(piplib_llmod)(a, b)))
  #define piplib_int_gcd(gcd, a, b) (gcd = (PIPLIB_NAME(piplib_int_t))(PIPLIB_NAME(piplib_llgcd_llabs)(a, b)))
  #define piplib_int_oppose(r, i) (r = - (i))
  #define piplib_int_size_in_base_2(i) (PIPLIB_NAME(piplib_lllog2)(i))
  #define piplib_int_size_in_base_10(i) (PIPLIB_NAME(piplib_lllog10)(i))

  #define piplib_int_eq(a, b) (a == b)
  #define piplib_int_ne(a, b) (a != b)
  #define piplib_int_zero(i) (i == 0)
  #define piplib_int_one(i) (i == 1)
  #define piplib_int_pos(i) (i > 0)
  #define piplib_int_neg(i) (i < 0)

#endif


#if defined(__cplusplus)
extern "C" 
  {
#endif 


/**
 * @brief Structure PipMatrix
 * 
 * Structure de matrice au format PolyLib. Le premier element d'une ligne
 * indique quand il vaut 1 que la ligne decrit une inequation de la forme
 * p(x)>=0 et quand il vaut 0, que la ligne decrit une egalite de la forme
 * p(x)=0. Le dernier element de chaque ligne correspond au coefficient
 * constant.
 */
struct PIPLIB_NAME(pipmatrix) {
  unsigned int NbRows;    /**< Number of rows */
  unsigned int NbColumns; /**< Number of columns */
  PIPLIB_NAME(piplib_int_t)** p;       /**< Data */
  PIPLIB_NAME(piplib_int_t)* p_Init;   /**< Init */
  int p_Init_size;        /**< Only for PolyLib compatibility under MP version:
                               PolyLib makes sometimes overestimates on the size
                               of the matrices, in order to go faster.
                               Thus NbRows*NbColumns is not the number of
                               allocated elements. With MP version, we have to
                               think to mpz_clear() all the initialized elements
                               before freing, then we need to know the number of
                               allocated elements: p_Init_size. */
};
typedef struct PIPLIB_NAME(pipmatrix) PIPLIB_NAME(PipMatrix);


/**
 * @brief Structure PipVector
 * 
 * Cette structure contient un Vector de 'nb_elements' la ieme composante de
 * ce vecteur vaut the_vector[i]/the_deno[i].
 */
struct PIPLIB_NAME(pipvector) {
  int nb_elements;          /**< Nombre d'elements du vecteur. */
  PIPLIB_NAME(piplib_int_t)* the_vector; /**< Numerateurs du vecteur. */
  PIPLIB_NAME(piplib_int_t)* the_deno;   /**< Denominateurs du vecteur. */
};
typedef struct PIPLIB_NAME(pipvector) PIPLIB_NAME(PipVector);


/**
 * @brief Structure PipNewparm
 * 
 * Liste chainee de Newparm, les informations d'un newparm etant son rang, un
 * vecteur de coefficients et un denominateur. Le newparm est egal a la division
 * du vecteur par le denominateur.
 */
struct PIPLIB_NAME(pipnewparm) {
  int rank;                /**< Rang du 'newparm'. */
  PIPLIB_NAME(PipVector)* vector;       /**< Le vector decrivant le newparm. */
  PIPLIB_NAME(piplib_int_t) deno;       /**< Denominateur du 'newparm'. */
  struct PIPLIB_NAME(pipnewparm)* next; /**< Pointeur vers le newparm suivant. */
};
typedef struct PIPLIB_NAME(pipnewparm) PIPLIB_NAME(PipNewparm);


/**
 * @brief Structure PipList
 * 
 * Liste chainee de Vector.
 */
struct PIPLIB_NAME(piplist) {
  PIPLIB_NAME(PipVector)* vector;    /**< Le vector contenant la partie de solution. */
  struct PIPLIB_NAME(piplist)* next; /**< Pointeur vers l'element suivant. */
};
typedef struct PIPLIB_NAME(piplist) PIPLIB_NAME(PipList);


/**
 * @brief Structure pipquast
 * 
 * Arbre binaire. Conformement a la grammaire de sortie (voir mode d'emploi), un
 * noeud de l'arbre des solutions debute par une liste de 'newparm'. Il continue
 * ensuite soit par une 'list' (alors condition vaut null), soit par un 'if'
 * (alors le champ condition contient la condition).
 */
struct PIPLIB_NAME(pipquast) {
  PIPLIB_NAME(PipNewparm)* newparm;        /**< Les 'newparm'. */
  PIPLIB_NAME(PipList)* list;              /**< La 'list' si pas de 'if'. */
  PIPLIB_NAME(PipVector)* condition;       /**< La condition si 'if'. */
  struct PIPLIB_NAME(pipquast)* next_then; /**< Noeud si condition et si verifiee. */
  struct PIPLIB_NAME(pipquast)* next_else; /**< Noeud si condition et si non verifiee. */
  struct PIPLIB_NAME(pipquast)* father;    /**< Pointeur vers le quast pere. */
};
typedef struct PIPLIB_NAME(pipquast) PIPLIB_NAME(PipQuast);


/**
 * @brief Structure pipoptions
 * 
 * This structure contains each option that can be set to change the PIP
 * behaviour.
 */
struct PIPLIB_NAME(pipoptions) {
  /** 1 if an integer solution is needed, 0 otherwise. */
  int Nq;
  
  /**
   * @brief Verbose
   * 
   * -1 -> absolute silence @n
   *  0 -> relative silence @n
   *  1 -> information on cuts when an integer solution is needed @n
   *  2 -> information sur les pivots et les déterminants @n
   *  3 -> information on arrays
   * 
   * Each option include the preceding
   */
  int Verbose;
  
  /** Set to 1 to eliminate some trivial solutions, 0 otherwise */
  int Simplify;
  
  /** Set to 1 to include deepest cut algorithm */
  int Deepest_cut;
  
  /** Set to 1 if maximum is needed */
  int Maximize;
  
  /**
   * @brief Signs of parameters
   * 
   * -1 -> all parameters may be negative @n
   *  0 -> all parameters are non-negative
   */
  int Urs_parms;
  
  /**
   * @brief Signs of unknowns
   * 
   * -1 -> all unknowns may be negative @n
   *  0 -> all unknowns are non-negative
   */
  int Urs_unknowns;
  
  /** To compute the dual */
  int Compute_dual;
};
typedef struct PIPLIB_NAME(pipoptions) PIPLIB_NAME(PipOptions);

void PIPLIB_NAME(pip_options_print)(FILE*, PIPLIB_NAME(PipOptions)*);


/* Fonctions d'affichages des structures de la PipLib. */
void PIPLIB_NAME(pip_matrix_print)(FILE*, PIPLIB_NAME(PipMatrix)*);
void PIPLIB_NAME(pip_vector_print)(FILE*, PIPLIB_NAME(PipVector)*);
void PIPLIB_NAME(pip_newparm_print)(FILE*, PIPLIB_NAME(PipNewparm)*, int);
void PIPLIB_NAME(pip_list_print)(FILE*, PIPLIB_NAME(PipList)*, int);
void PIPLIB_NAME(pip_quast_print)(FILE*, PIPLIB_NAME(PipQuast)*, int);


/* Fonctions de liberation memoire des structures de la PipLib.*/
void PIPLIB_NAME(pip_matrix_free)(PIPLIB_NAME(PipMatrix)*);
void PIPLIB_NAME(pip_vector_free)(PIPLIB_NAME(PipVector)*);
void PIPLIB_NAME(pip_newparm_free)(PIPLIB_NAME(PipNewparm)*);
void PIPLIB_NAME(pip_list_free)(PIPLIB_NAME(PipList)*);
void PIPLIB_NAME(pip_quast_free)(PIPLIB_NAME(PipQuast)*);
void PIPLIB_NAME(pip_options_free)(PIPLIB_NAME(PipOptions)*);


/* Fonctions d'acquisition de matrices de contraintes et options. */
PIPLIB_NAME(PipMatrix)* PIPLIB_NAME(pip_matrix_alloc)(unsigned int, unsigned int);
PIPLIB_NAME(PipMatrix)* PIPLIB_NAME(pip_matrix_read)(FILE*);
PIPLIB_NAME(PipOptions)* PIPLIB_NAME(pip_options_init)(void);


/* initialization of pip library */
void PIPLIB_NAME(pip_init)();
void PIPLIB_NAME(pip_close)();


/**
 * @brief Fonction de resolution
 * 
 * PIPLIB_NAME(pip_solve) resoud le probleme qu'on lui passe en parametre, suivant les
 * options elles aussi en parametre. Elle renvoie la solution sous forme
 * d'un arbre de PipQuast.
 * 
 * @param[in] domain     Inequations definissant le domaine des inconnues
 * @param[in] parameters Inequations satisfaites par les parametres
 * @param[in] bignum     Column rank of the bignum, or negative value if there is no big parameter
 * @param[in] options    PipLib options
 */ 
PIPLIB_NAME(PipQuast)* PIPLIB_NAME(pip_solve)(PIPLIB_NAME(PipMatrix)* domain, PIPLIB_NAME(PipMatrix)* parameters,
                    int bignum, PIPLIB_NAME(PipOptions)* options);


/** sol_quast_edit */
PIPLIB_NAME(PipQuast)* PIPLIB_NAME(sol_quast_edit)(int* i, PIPLIB_NAME(PipQuast)* father,
                         int Bg, int Urs_p, int flags);


#if defined(__cplusplus)
  }
#endif



/* Old names (undef) */

#undef piplib_int_t
#undef PipMatrix
#undef PipVector
#undef PipNewparm
#undef PipList
#undef PipQuast
#undef PipOptions
#undef PipMatrix

#undef pip_options_print

#undef pip_matrix_print
#undef pip_vector_print
#undef pip_newparm_print
#undef pip_list_print
#undef pip_quast_print

#undef pip_matrix_free
#undef pip_vector_free
#undef pip_newparm_free
#undef pip_list_free
#undef pip_quast_free
#undef pip_options_free

#undef pip_matrix_alloc
#undef pip_matrix_read
#undef pip_options_init

#undef pip_init
#undef pip_close

#undef pip_solve

#undef sol_quast_edit

#undef Entier

/* Old names (define) */

#define piplib_int_t piplib_int_t_dp
#define PipMatrix PipMatrix_dp
#define PipVector PipVector_dp
#define PipNewparm PipNewparm_dp
#define PipList PipList_dp
#define PipQuast PipQuast_dp
#define PipOptions PipOptions_dp
#define PipMatrix PipMatrix_dp

#define pip_options_print pip_options_print_dp

#define pip_matrix_print pip_matrix_print_dp
#define pip_vector_print pip_vector_print_dp
#define pip_newparm_print pip_newparm_print_dp
#define pip_list_print pip_list_print_dp
#define pip_quast_print pip_quast_print_dp

#define pip_matrix_free pip_matrix_free_dp
#define pip_vector_free pip_vector_free_dp
#define pip_newparm_free pip_newparm_free_dp
#define pip_list_free pip_list_free_dp
#define pip_quast_free pip_quast_free_dp
#define pip_options_free pip_options_free_dp

#define pip_matrix_alloc pip_matrix_alloc_dp
#define pip_matrix_read pip_matrix_read_dp
#define pip_options_init pip_options_init_dp

#define pip_init pip_init_dp
#define pip_close pip_close_dp

#define pip_solve pip_solve_dp

#define sol_quast_edit sol_quast_edit_dp

#define Entier piplib_int_t_dp
