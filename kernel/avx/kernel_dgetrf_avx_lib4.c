/**************************************************************************************************
*                                                                                                 *
* This file is part of HPMPC.                                                                     *
*                                                                                                 *
* HPMPC -- Library for High-Performance implementation of solvers for MPC.                        *
* Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.                *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, giaf (at) dtu.dk                                                       *
*                                                                                                 *
**************************************************************************************************/

#include <mmintrin.h>
#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <pmmintrin.h>  // SSE3
#include <smmintrin.h>  // SSE4
#include <immintrin.h>  // AVX



void kernel_dgetrf_r_nn_8x4_lib4(int kmax, double *A0, int sda, double *B, int sdb, int alg, double *C0, int sdc, double *D0, int sdd, double *inv_diag_D, double *E0, int sde)
	{

	double *A1 = A0 + 4*sda;
	double *C1 = C0 + 4*sdc;
	double *D1 = D0 + 4*sdd;
	double *E1 = E0 + 4*sde;

	__builtin_prefetch( B+0 );
	__builtin_prefetch( B+8 );

	const int bs = 4;

	const int B_next = bs*sdb;

	__builtin_prefetch( B+B_next+0 );
	__builtin_prefetch( B+B_next+8 );

	int k;

	__m128d
		ones,
		v_0;

	__m256d
		zeros,
		a_0, a_4, A_0, A_4,
		b_0,
		tmp,
		d_0, d_1, d_2, d_3,
		d_4, d_5, d_6, d_7,
		c_0, 
		e_0, t_0,
		u_0, u_1;
	

	// correction phase

	d_0 = _mm256_setzero_pd();
	d_1 = _mm256_setzero_pd();
	d_2 = _mm256_setzero_pd();
	d_3 = _mm256_setzero_pd();
	d_4 = _mm256_setzero_pd();
	d_5 = _mm256_setzero_pd();
	d_6 = _mm256_setzero_pd();
	d_7 = _mm256_setzero_pd();

	if(kmax<=0)
		goto add;

	// prefetch
	a_0 = _mm256_load_pd( &A0[0+bs*0] );
	a_4 = _mm256_load_pd( &A1[0+bs*0] );

	for(k=0; k<kmax-3; k+=4)
		{

		__builtin_prefetch( B+2*B_next+0 );

		b_0 = _mm256_broadcast_sd( &B[0+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		A_0 = _mm256_load_pd( &A0[0+bs*1] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		A_4 = _mm256_load_pd( &A1[0+bs*1] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		b_0 = _mm256_broadcast_sd( &B[1+bs*0] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		a_0 = _mm256_load_pd( &A0[0+bs*2] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*1] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		a_4 = _mm256_load_pd( &A1[0+bs*2] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*2] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*3] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );

		__builtin_prefetch( B+2*B_next+8 );

		b_0 = _mm256_broadcast_sd( &B[2+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		A_0 = _mm256_load_pd( &A0[0+bs*3] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		A_4 = _mm256_load_pd( &A1[0+bs*3] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		b_0 = _mm256_broadcast_sd( &B[3+bs*0] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		a_0 = _mm256_load_pd( &A0[0+bs*4] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*1] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		a_4 = _mm256_load_pd( &A1[0+bs*4] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*2] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*3] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		A0 += 4*bs;
		A1 += 4*bs;
		B  += B_next;

		}
	for(; k<kmax; k++)
		{

		a_0 = _mm256_load_pd( &A0[0+bs*0] );
		a_4 = _mm256_load_pd( &A1[0+bs*0] );
		b_0 = _mm256_broadcast_sd( &B[0+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		A0 += 1*bs;
		A1 += 1*bs;
		B  += 1;

		}
	
	add:

	if(alg!=0)
		{
		c_0 = _mm256_load_pd( &C0[0+bs*0] );
		d_0 = _mm256_add_pd( c_0, d_0 );
		c_0 = _mm256_load_pd( &C0[0+bs*1] );
		d_1 = _mm256_add_pd( c_0, d_1 );
		c_0 = _mm256_load_pd( &C0[0+bs*2] );
		d_2 = _mm256_add_pd( c_0, d_2 );
		c_0 = _mm256_load_pd( &C0[0+bs*3] );
		d_3 = _mm256_add_pd( c_0, d_3 );
		c_0 = _mm256_load_pd( &C1[0+bs*0] );
		d_4 = _mm256_add_pd( c_0, d_4 );
		c_0 = _mm256_load_pd( &C1[0+bs*1] );
		d_5 = _mm256_add_pd( c_0, d_5 );
		c_0 = _mm256_load_pd( &C1[0+bs*2] );
		d_6 = _mm256_add_pd( c_0, d_6 );
		c_0 = _mm256_load_pd( &C1[0+bs*3] );
		d_7 = _mm256_add_pd( c_0, d_7 );
		}
	


	// solution of top & correction of bottom

	zeros = _mm256_setzero_pd();

	e_0 = _mm256_load_pd( &E0[0+bs*0] );
	a_0 = _mm256_load_pd( &E1[0+bs*0] );
	e_0 = _mm256_blend_pd( e_0, zeros, 0x1 );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_0 ) );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_0 = _mm256_sub_pd( d_0, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_4 = _mm256_sub_pd( d_4, tmp );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_1 ) );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_1 = _mm256_sub_pd( d_1, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_5 = _mm256_sub_pd( d_5, tmp );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_2 ) );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_3 ) );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_7 = _mm256_sub_pd( d_7, tmp );

	e_0 = _mm256_load_pd( &E0[0+bs*1] );
	a_0 = _mm256_load_pd( &E1[0+bs*1] );
	e_0 = _mm256_blend_pd( e_0, zeros, 0x3 );
	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_0 ), 0x3 );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_0 = _mm256_sub_pd( d_0, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_4 = _mm256_sub_pd( d_4, tmp );
	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_1 ), 0x3 );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_1 = _mm256_sub_pd( d_1, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_5 = _mm256_sub_pd( d_5, tmp );
	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_2 ), 0x3 );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_3 ), 0x3 );
	t_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x00 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_7 = _mm256_sub_pd( d_7, tmp );

	e_0 = _mm256_load_pd( &E0[0+bs*2] );
	a_0 = _mm256_load_pd( &E1[0+bs*2] );
	e_0 = _mm256_blend_pd( e_0, zeros, 0x7 );
	t_0 = _mm256_permute2f128_pd( d_0, d_0, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0x0 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_0 = _mm256_sub_pd( d_0, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_4 = _mm256_sub_pd( d_4, tmp );
	t_0 = _mm256_permute2f128_pd( d_1, d_1, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0x0 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_1 = _mm256_sub_pd( d_1, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_5 = _mm256_sub_pd( d_5, tmp );
	t_0 = _mm256_permute2f128_pd( d_2, d_2, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0x0 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	t_0 = _mm256_permute2f128_pd( d_3, d_3, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0x0 );
	tmp = _mm256_mul_pd( e_0, t_0 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_7 = _mm256_sub_pd( d_7, tmp );

	a_0 = _mm256_load_pd( &E1[0+bs*3] );
	t_0 = _mm256_permute2f128_pd( d_0, d_0, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0xf );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_4 = _mm256_sub_pd( d_4, tmp );
	t_0 = _mm256_permute2f128_pd( d_1, d_1, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0xf );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_5 = _mm256_sub_pd( d_5, tmp );
	t_0 = _mm256_permute2f128_pd( d_2, d_2, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0xf );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	t_0 = _mm256_permute2f128_pd( d_3, d_3, 0x11 );
	t_0 = _mm256_permute_pd( t_0, 0xf );
	tmp = _mm256_mul_pd( a_0, t_0 );
	d_7 = _mm256_sub_pd( d_7, tmp );

	// store upper
	_mm256_store_pd( &D0[0+bs*0], d_0 );
	_mm256_store_pd( &D0[0+bs*1], d_1 );
	_mm256_store_pd( &D0[0+bs*2], d_2 );
	_mm256_store_pd( &D0[0+bs*3], d_3 );


	// factorization of bottom

	ones = _mm_set_pd( 1.0, 1.0 );


	// first column
	u_0 = _mm256_blend_pd( u_0, d_4, 0x1 );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_4 ) );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_4 = _mm256_mul_pd( d_4, u_1 );
	_mm_store_sd( &inv_diag_D[0], v_0 );
	d_4 = _mm256_blend_pd( d_4, u_0, 0x1 );
	_mm256_store_pd( &D1[0+bs*0], d_4 );


	// second column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_5 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_4, u_0 );
	d_5 = _mm256_sub_pd( d_5, tmp );
	u_0 = _mm256_blend_pd( u_0, d_5, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_5 ), 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_5 = _mm256_mul_pd( d_5, u_1 );
	_mm_store_sd( &inv_diag_D[1], v_0 );
	d_5 = _mm256_blend_pd( d_5, u_0, 0x3 );
	_mm256_store_pd( &D1[0+bs*1], d_5 );


	// third column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_6 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_4, u_0 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	u_0 = _mm256_blend_pd( u_0, d_6, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_6 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_5, u_1 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	u_0 = _mm256_blend_pd( u_0, d_6, 0x4 );

	v_0 = _mm256_extractf128_pd( d_6, 0x1 );
	v_0 = _mm_movedup_pd( v_0 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_6 = _mm256_mul_pd( d_6, u_1 );
	_mm_store_sd( &inv_diag_D[2], v_0 );
	d_6 = _mm256_blend_pd( d_6, u_0, 0x7 );
	_mm256_store_pd( &D1[0+bs*2], d_6 );


	// fourth column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_7 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_4, u_0 );
	d_7 = _mm256_sub_pd( d_7, tmp );
	u_0 = _mm256_blend_pd( u_0, d_7, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_7 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_5, u_1 );
	d_7 = _mm256_sub_pd( d_7, tmp );
	u_0 = _mm256_blend_pd( u_0, d_7, 0x4 );

	u_1 = _mm256_permute2f128_pd( d_7, d_7, 0x11 );
	u_1 = _mm256_permute_pd( u_1, 0x0 );
	tmp = _mm256_mul_pd( d_6, u_1 );
	d_7 = _mm256_sub_pd( d_7, tmp );
	u_0 = _mm256_blend_pd( u_0, d_7, 0x8 );
	_mm256_store_pd( &D1[0+bs*3], u_0 );

	v_0 = _mm256_extractf128_pd( d_7, 0x1 );
	v_0 = _mm_permute_pd( v_0, 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	_mm_store_sd( &inv_diag_D[3], v_0 );

	return;

	}



void kernel_dgetrf_l_nn_8x4_lib4(int kmax, double *A0, int sda, double *B, int sdb, int alg, double *C0, int sdc, double *D0, int sdd, double *inv_diag_D)
	{

	double *A1 = A0 + 4*sda;
	double *C1 = C0 + 4*sdc;
	double *D1 = D0 + 4*sdd;

	__builtin_prefetch( B+0 );
	__builtin_prefetch( B+8 );

	const int bs = 4;

	const int B_next = bs*sdb;

	__builtin_prefetch( B+B_next+0 );
	__builtin_prefetch( B+B_next+8 );

	int k;

	__m128d
		ones,
		v_0;

	__m256d
		a_0, a_4, A_0, A_4,
		b_0,
		tmp,
		d_0, d_1, d_2, d_3,
		d_4, d_5, d_6, d_7,
		c_0, 
		u_0, u_1;
	

	// correction phase

	d_0 = _mm256_setzero_pd();
	d_1 = _mm256_setzero_pd();
	d_2 = _mm256_setzero_pd();
	d_3 = _mm256_setzero_pd();
	d_4 = _mm256_setzero_pd();
	d_5 = _mm256_setzero_pd();
	d_6 = _mm256_setzero_pd();
	d_7 = _mm256_setzero_pd();

	if(kmax<=0)
		goto add;

	// prefetch
	a_0 = _mm256_load_pd( &A0[0+bs*0] );
	a_4 = _mm256_load_pd( &A1[0+bs*0] );

	for(k=0; k<kmax-3; k+=4)
		{

		__builtin_prefetch( B+2*B_next+0 );

		b_0 = _mm256_broadcast_sd( &B[0+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		A_0 = _mm256_load_pd( &A0[0+bs*1] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		A_4 = _mm256_load_pd( &A1[0+bs*1] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		b_0 = _mm256_broadcast_sd( &B[1+bs*0] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		a_0 = _mm256_load_pd( &A0[0+bs*2] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*1] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		a_4 = _mm256_load_pd( &A1[0+bs*2] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*2] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*3] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );

		__builtin_prefetch( B+2*B_next+8 );

		b_0 = _mm256_broadcast_sd( &B[2+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		A_0 = _mm256_load_pd( &A0[0+bs*3] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		A_4 = _mm256_load_pd( &A1[0+bs*3] );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		b_0 = _mm256_broadcast_sd( &B[3+bs*0] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		a_0 = _mm256_load_pd( &A0[0+bs*4] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*1] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		a_4 = _mm256_load_pd( &A1[0+bs*4] );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*2] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*3] );
		tmp = _mm256_mul_pd( A_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( A_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		A0 += 4*bs;
		A1 += 4*bs;
		B  += B_next;

		}
	for(; k<kmax; k++)
		{

		a_0 = _mm256_load_pd( &A0[0+bs*0] );
		a_4 = _mm256_load_pd( &A1[0+bs*0] );
		b_0 = _mm256_broadcast_sd( &B[0+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_4 = _mm256_sub_pd( d_4, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_5 = _mm256_sub_pd( d_5, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_6 = _mm256_sub_pd( d_6, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );
		tmp = _mm256_mul_pd( a_4, b_0 );
		d_7 = _mm256_sub_pd( d_7, tmp );


		A0 += 1*bs;
		A1 += 1*bs;
		B  += 1;

		}
	
	add:

	if(alg!=0)
		{
		c_0 = _mm256_load_pd( &C0[0+bs*0] );
		d_0 = _mm256_add_pd( c_0, d_0 );
		c_0 = _mm256_load_pd( &C0[0+bs*1] );
		d_1 = _mm256_add_pd( c_0, d_1 );
		c_0 = _mm256_load_pd( &C0[0+bs*2] );
		d_2 = _mm256_add_pd( c_0, d_2 );
		c_0 = _mm256_load_pd( &C0[0+bs*3] );
		d_3 = _mm256_add_pd( c_0, d_3 );
		c_0 = _mm256_load_pd( &C1[0+bs*0] );
		d_4 = _mm256_add_pd( c_0, d_4 );
		c_0 = _mm256_load_pd( &C1[0+bs*1] );
		d_5 = _mm256_add_pd( c_0, d_5 );
		c_0 = _mm256_load_pd( &C1[0+bs*2] );
		d_6 = _mm256_add_pd( c_0, d_6 );
		c_0 = _mm256_load_pd( &C1[0+bs*3] );
		d_7 = _mm256_add_pd( c_0, d_7 );
		}
	


	// factorization

	ones = _mm_set_pd( 1.0, 1.0 );


	// first column
	u_0 = _mm256_blend_pd( u_0, d_0, 0x1 );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_0 ) );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_0 = _mm256_mul_pd( d_0, u_1 );
	d_4 = _mm256_mul_pd( d_4, u_1 );
	_mm_store_sd( &inv_diag_D[0], v_0 );
	d_0 = _mm256_blend_pd( d_0, u_0, 0x1 );
	_mm256_store_pd( &D0[0+bs*0], d_0 );
	_mm256_store_pd( &D1[0+bs*0], d_4 );


	// second column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_1 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_0, u_0 );
	d_1 = _mm256_sub_pd( d_1, tmp );
	tmp = _mm256_mul_pd( d_4, u_0 );
	d_5 = _mm256_sub_pd( d_5, tmp );
	u_0 = _mm256_blend_pd( u_0, d_1, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_1 ), 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_1 = _mm256_mul_pd( d_1, u_1 );
	d_5 = _mm256_mul_pd( d_5, u_1 );
	_mm_store_sd( &inv_diag_D[1], v_0 );
	d_1 = _mm256_blend_pd( d_1, u_0, 0x3 );
	_mm256_store_pd( &D0[0+bs*1], d_1 );
	_mm256_store_pd( &D1[0+bs*1], d_5 );


	// third column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_2 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_0, u_0 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	tmp = _mm256_mul_pd( d_4, u_0 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	u_0 = _mm256_blend_pd( u_0, d_2, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_2 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_1, u_1 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	tmp = _mm256_mul_pd( d_5, u_1 );
	d_6 = _mm256_sub_pd( d_6, tmp );
	u_0 = _mm256_blend_pd( u_0, d_2, 0x4 );

	v_0 = _mm256_extractf128_pd( d_2, 0x1 );
	v_0 = _mm_movedup_pd( v_0 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_2 = _mm256_mul_pd( d_2, u_1 );
	d_6 = _mm256_mul_pd( d_6, u_1 );
	_mm_store_sd( &inv_diag_D[2], v_0 );
	d_2 = _mm256_blend_pd( d_2, u_0, 0x7 );
	_mm256_store_pd( &D0[0+bs*2], d_2 );
	_mm256_store_pd( &D1[0+bs*2], d_6 );


	// fourth column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_3 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_0, u_0 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	tmp = _mm256_mul_pd( d_4, u_0 );
	d_7 = _mm256_sub_pd( d_7, tmp );
	u_0 = _mm256_blend_pd( u_0, d_3, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_3 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( d_1, u_1 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	tmp = _mm256_mul_pd( d_5, u_1 );
	d_7 = _mm256_sub_pd( d_7, tmp );
	u_0 = _mm256_blend_pd( u_0, d_3, 0x4 );

	u_1 = _mm256_permute2f128_pd( d_3, d_3, 0x11 );
	u_1 = _mm256_permute_pd( u_1, 0x0 );
	tmp = _mm256_mul_pd( d_2, u_1 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	tmp = _mm256_mul_pd( d_6, u_1 );
	d_7 = _mm256_sub_pd( d_7, tmp );
	u_0 = _mm256_blend_pd( u_0, d_3, 0x8 );
	_mm256_store_pd( &D0[0+bs*3], u_0 );

	v_0 = _mm256_extractf128_pd( d_3, 0x1 );
	v_0 = _mm_permute_pd( v_0, 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	_mm_store_sd( &inv_diag_D[3], v_0 );
	d_7 = _mm256_mul_pd( d_7, u_1 );
	_mm256_store_pd( &D1[0+bs*3], d_7 );

	return;

	}



void kernel_dgetrf_nn_4x4_lib4(int kmax, double *A, double *B, int sdb, int alg, double *C, double *D, double *inv_diag_D)
	{

	const int bs = 4;

	int k;

	__m128d
		ones,
		v_0;

	__m256d
		a_0,
		b_0,
		tmp,
		d_0, d_1, d_2, d_3,
		c_0, 
		u_0, u_1;
	

	// correction phase

	d_0 = _mm256_setzero_pd();
	d_1 = _mm256_setzero_pd();
	d_2 = _mm256_setzero_pd();
	d_3 = _mm256_setzero_pd();

	for(k=0; k<kmax-3; k+=4)
		{

		a_0 = _mm256_load_pd( &A[0+bs*0] );
		b_0 = _mm256_broadcast_sd( &B[0+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );


		a_0 = _mm256_load_pd( &A[0+bs*1] );
		b_0 = _mm256_broadcast_sd( &B[1+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		b_0 = _mm256_broadcast_sd( &B[1+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );


		a_0 = _mm256_load_pd( &A[0+bs*2] );
		b_0 = _mm256_broadcast_sd( &B[2+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		b_0 = _mm256_broadcast_sd( &B[2+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );


		a_0 = _mm256_load_pd( &A[0+bs*3] );
		b_0 = _mm256_broadcast_sd( &B[3+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		b_0 = _mm256_broadcast_sd( &B[3+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );


		A += 4*bs;
		B += bs*sdb;

		}
	for(; k<kmax; k++)
		{

		a_0 = _mm256_load_pd( &A[0+bs*0] );
		b_0 = _mm256_broadcast_sd( &B[0+bs*0] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_0 = _mm256_sub_pd( d_0, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*1] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_1 = _mm256_sub_pd( d_1, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*2] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_2 = _mm256_sub_pd( d_2, tmp );
		b_0 = _mm256_broadcast_sd( &B[0+bs*3] );
		tmp = _mm256_mul_pd( a_0, b_0 );
		d_3 = _mm256_sub_pd( d_3, tmp );


		A += 1*bs;
		B += 1;

		}
	
	if(alg!=0)
		{
		c_0 = _mm256_load_pd( &C[0+bs*0] );
		d_0 = _mm256_add_pd( c_0, d_0 );
		c_0 = _mm256_load_pd( &C[0+bs*1] );
		d_1 = _mm256_add_pd( c_0, d_1 );
		c_0 = _mm256_load_pd( &C[0+bs*2] );
		d_2 = _mm256_add_pd( c_0, d_2 );
		c_0 = _mm256_load_pd( &C[0+bs*3] );
		d_3 = _mm256_add_pd( c_0, d_3 );
		}
	


	// factorization

	ones = _mm_set_pd( 1.0, 1.0 );


	// first column
	u_0 = _mm256_blend_pd( u_0, d_0, 0x1 );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_0 ) );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_0 = _mm256_mul_pd( d_0, u_1 );
	_mm_store_sd( &inv_diag_D[0], v_0 );
	d_0 = _mm256_blend_pd( d_0, u_0, 0x1 );
	_mm256_store_pd( &D[0+bs*0], d_0 );


	// second column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_1 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_0, d_0 );
	d_1 = _mm256_sub_pd( d_1, tmp );
	u_0 = _mm256_blend_pd( u_0, d_1, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_1 ), 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_1 = _mm256_mul_pd( d_1, u_1 );
	_mm_store_sd( &inv_diag_D[1], v_0 );
	d_1 = _mm256_blend_pd( d_1, u_0, 0x3 );
	_mm256_store_pd( &D[0+bs*1], d_1 );


	// third column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_2 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_0, d_0 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	u_0 = _mm256_blend_pd( u_0, d_2, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_2 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_1, d_1 );
	d_2 = _mm256_sub_pd( d_2, tmp );
	u_0 = _mm256_blend_pd( u_0, d_2, 0x4 );

	v_0 = _mm256_extractf128_pd( d_2, 0x1 );
	v_0 = _mm_movedup_pd( v_0 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	d_2 = _mm256_mul_pd( d_2, u_1 );
	_mm_store_sd( &inv_diag_D[2], v_0 );
	d_2 = _mm256_blend_pd( d_2, u_0, 0x7 );
	_mm256_store_pd( &D[0+bs*2], d_2 );


	// fourth column
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( d_3 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_0, d_0 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	u_0 = _mm256_blend_pd( u_0, d_3, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( d_3 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_1, d_1 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	u_0 = _mm256_blend_pd( u_0, d_3, 0x4 );

	u_1 = _mm256_permute2f128_pd( d_3, d_3, 0x11 );
	u_1 = _mm256_permute_pd( u_1, 0x0 );
	tmp = _mm256_mul_pd( u_1, d_2 );
	d_3 = _mm256_sub_pd( d_3, tmp );
	u_0 = _mm256_blend_pd( u_0, d_3, 0x8 );
	_mm256_store_pd( &D[0+bs*3], u_0 );

	v_0 = _mm256_extractf128_pd( d_3, 0x1 );
	v_0 = _mm_permute_pd( v_0, 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	_mm_store_sd( &inv_diag_D[3], v_0 );

	return;

	}



void corner_dgetrf_nn_4x4_lib4(double *C, double *LU, double *inv_diag_U)
	{

	const int bs = 4;

	__m128d
		ones,
		v_0;

	__m256d
		tmp,
		c_0, c_1, c_2, c_3,
		u_0, u_1;
	
	ones = _mm_set_pd( 1.0, 1.0 );



	// first column
	c_0 = _mm256_load_pd( &C[0+bs*0] );

	u_0 = _mm256_blend_pd( u_0, c_0, 0x1 );
	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( c_0 ) );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	c_0 = _mm256_mul_pd( c_0, u_1 );
	_mm_store_sd( &inv_diag_U[0], v_0 );
	c_0 = _mm256_blend_pd( c_0, u_0, 0x1 );
	_mm256_store_pd( &LU[0+bs*0], c_0 );



	// second column
	c_1 = _mm256_load_pd( &C[0+bs*1] );

	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( c_1 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_0, c_0 );
	c_1 = _mm256_sub_pd( c_1, tmp );
	u_0 = _mm256_blend_pd( u_0, c_1, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( c_1 ), 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	c_1 = _mm256_mul_pd( c_1, u_1 );
	_mm_store_sd( &inv_diag_U[1], v_0 );
	c_1 = _mm256_blend_pd( c_1, u_0, 0x3 );
	_mm256_store_pd( &LU[0+bs*1], c_1 );



	// third column
	c_2 = _mm256_load_pd( &C[0+bs*2] );

	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( c_2 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_0, c_0 );
	c_2 = _mm256_sub_pd( c_2, tmp );
	u_0 = _mm256_blend_pd( u_0, c_2, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( c_2 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_1, c_1 );
	c_2 = _mm256_sub_pd( c_2, tmp );
	u_0 = _mm256_blend_pd( u_0, c_2, 0x4 );

	v_0 = _mm256_extractf128_pd( c_2, 0x1 );
	v_0 = _mm_movedup_pd( v_0 );
	v_0 = _mm_div_pd( ones, v_0 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	c_2 = _mm256_mul_pd( c_2, u_1 );
	_mm_store_sd( &inv_diag_U[2], v_0 );
	c_2 = _mm256_blend_pd( c_2, u_0, 0x7 );
	_mm256_store_pd( &LU[0+bs*2], c_2 );



	// fourth column
	c_3 = _mm256_load_pd( &C[0+bs*3] );

	v_0 = _mm_movedup_pd( _mm256_castpd256_pd128( c_3 ) );
	u_0 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_0, c_0 );
	c_3 = _mm256_sub_pd( c_3, tmp );
	u_0 = _mm256_blend_pd( u_0, c_3, 0x2 );

	v_0 = _mm_permute_pd( _mm256_castpd256_pd128( c_3 ), 0x3 );
	u_1 = _mm256_permute2f128_pd( _mm256_castpd128_pd256( v_0 ), _mm256_castpd128_pd256( v_0 ), 0x0 );
	tmp = _mm256_mul_pd( u_1, c_1 );
	c_3 = _mm256_sub_pd( c_3, tmp );
	u_0 = _mm256_blend_pd( u_0, c_3, 0x4 );

	u_1 = _mm256_permute2f128_pd( c_3, c_3, 0x11 );
	u_1 = _mm256_permute_pd( u_1, 0x0 );
	tmp = _mm256_mul_pd( u_1, c_2 );
	c_3 = _mm256_sub_pd( c_3, tmp );
	u_0 = _mm256_blend_pd( u_0, c_3, 0x8 );
	_mm256_store_pd( &LU[0+bs*3], u_0 );

	v_0 = _mm256_extractf128_pd( c_3, 0x1 );
	v_0 = _mm_permute_pd( v_0, 0x3 );
	v_0 = _mm_div_pd( ones, v_0 );
	_mm_store_sd( &inv_diag_U[3], v_0 );


	}

