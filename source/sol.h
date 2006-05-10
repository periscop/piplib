/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                  sol.h                                     *
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

#ifndef SOL_H
#define SOL_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

void sol_init(void);
int sol_hwm(void);
void sol_reset(int);
void sol_nil(void);
void sol_if(void);
void sol_list(int);
void sol_form(int);
void sol_new(int);
void sol_div(void);
void sol_val(Entier, Entier);
int sol_edit(FILE *, int);
int is_not_Nil(int);
void sol_simplify(int);

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
