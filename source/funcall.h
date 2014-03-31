/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 funcall.h                                  *
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
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 *****************************************************************************/

#ifndef FUNCALL_H
#define FUNCALL_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

#define TRAITER_INT		(1 << 0)	/* Compute integer optimum */
#define TRAITER_DUAL		(1 << 1)	/* Compute dual variables */
void PIPLIB_NAME(traiter)(PIPLIB_NAME(Tableau) *tp, PIPLIB_NAME(Tableau) *ctxt, int nvar, int nparm, int ni, int nc,
	     int bigparm, int flags);
int PIPLIB_NAME(integrer)(PIPLIB_NAME(Tableau) **, PIPLIB_NAME(Tableau) **, int *, int *, int *, int *, int);

int PIPLIB_NAME(dgetc)(FILE *foo);
FILE *PIPLIB_NAME(pip_create_dump_file)();
int PIPLIB_NAME(sol_hwm)(void);
void PIPLIB_NAME(sol_simplify)(int);
int PIPLIB_NAME(is_not_Nil)(int);
int PIPLIB_NAME(sol_edit)(FILE *, int);
void PIPLIB_NAME(tab_reset)(struct PIPLIB_NAME(high_water_mark));
void PIPLIB_NAME(sol_reset)(int);
struct PIPLIB_NAME(high_water_mark) PIPLIB_NAME(tab_hwm)(void);
PIPLIB_NAME(Tableau) *PIPLIB_NAME(tab_get)(FILE *, int,int,int);
int PIPLIB_NAME(tab_simplify)(PIPLIB_NAME(Tableau) *tp, int cst);
void PIPLIB_NAME(sol_init)(void);
void PIPLIB_NAME(sol_close)(void);
void PIPLIB_NAME(tab_init)(void);
void PIPLIB_NAME(tab_close)(void);
void PIPLIB_NAME(sol_if)(void);
void PIPLIB_NAME(sol_forme)(int);
void PIPLIB_NAME(sol_val)(PIPLIB_NAME(piplib_int_t), PIPLIB_NAME(piplib_int_t));
void PIPLIB_NAME(sol_val_one)(PIPLIB_NAME(piplib_int_t));
void PIPLIB_NAME(sol_val_zero_one)();
void PIPLIB_NAME(sol_nil)(void);
void PIPLIB_NAME(sol_error)(int);
PIPLIB_NAME(Tableau) * PIPLIB_NAME(tab_alloc)(int, int, int);
void PIPLIB_NAME(sol_list)(int);
void PIPLIB_NAME(tab_display)(PIPLIB_NAME(Tableau) *, FILE *);
PIPLIB_NAME(Tableau) * PIPLIB_NAME(expanser)(PIPLIB_NAME(Tableau) *, int, int, int, int, int, int);
void PIPLIB_NAME(sol_new)(int);
void PIPLIB_NAME(sol_div)(void);

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
