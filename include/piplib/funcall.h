/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 funcall.h                                  *
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

Entier pgcd(Entier, Entier);
int sol_hwm(void);
void sol_simplify(int);
int is_not_Nil(int);
int sol_edit(FILE *, int);
void tab_reset(struct high_water_mark);
void sol_reset(int);
struct high_water_mark tab_hwm(void);
Tableau *tab_get(FILE *, int,int,int);
void sol_init(void);
void tab_init(void);
void sol_if(void);
void sol_forme(int);
void sol_val(Entier, Entier);
void sol_nil(void);
void sol_error(int);
int integrer(Tableau **, Tableau **, Entier, int *, int *, int *, int *);
Tableau * tab_alloc(int, int, int);
void sol_list(int);
void tab_display(Tableau *, FILE *);
Entier traiter(Tableau *, Tableau *, int, Entier, int, int, int, int, int);
Tableau * expanser(Tableau *, int, int, int, int, int, int);
int llog(Entier);
Entier mod(Entier,Entier);
void sol_new(int);
void sol_div(void);
