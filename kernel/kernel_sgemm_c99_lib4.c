/**************************************************************************************************
*                                                                                                 *
* This file is part of HPMPC.                                                                     *
*                                                                                                 *
* HPMPC -- Library for High-Performance implementation of solvers for MPC.                        *
* Copyright (C) 2014 by Technical University of Denmark. All rights reserved.                     *
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

// normal-transposed, 4x4 with data packed in 4
void kernel_sgemm_pp_nt_4x4_lib4(int kmax, float *A, float *B, float *C, int ldc, int alg)
	{
	
	if(kmax<=0)
		return;

/*	const int bs = 4;*/

	int k;
	
	float
		a_0, a_1, a_2, a_3,
		b_0, b_1, b_2, b_3,
		c_00=0, c_01=0, c_02=0, c_03=0,
		c_10=0, c_11=0, c_12=0, c_13=0,
		c_20=0, c_21=0, c_22=0, c_23=0,
		c_30=0, c_31=0, c_32=0, c_33=0;
	
	k = 0;
	for(; k<kmax-3; k+=4)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		b_1 = B[1];
		b_2 = B[2];
		b_3 = B[3];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;

		c_03 += a_0 * b_3;
		c_13 += a_1 * b_3;
		c_23 += a_2 * b_3;
		c_33 += a_3 * b_3;
		
		
		a_0 = A[4];
		a_1 = A[5];
		a_2 = A[6];
		a_3 = A[7];
		
		b_0 = B[4];
		b_1 = B[5];
		b_2 = B[6];
		b_3 = B[7];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;

		c_03 += a_0 * b_3;
		c_13 += a_1 * b_3;
		c_23 += a_2 * b_3;
		c_33 += a_3 * b_3;
		
		
		a_0 = A[8];
		a_1 = A[9];
		a_2 = A[10];
		a_3 = A[11];
		
		b_0 = B[8];
		b_1 = B[9];
		b_2 = B[10];
		b_3 = B[11];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;

		c_03 += a_0 * b_3;
		c_13 += a_1 * b_3;
		c_23 += a_2 * b_3;
		c_33 += a_3 * b_3;
		
		
		a_0 = A[12];
		a_1 = A[13];
		a_2 = A[14];
		a_3 = A[15];
		
		b_0 = B[12];
		b_1 = B[13];
		b_2 = B[14];
		b_3 = B[15];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;

		c_03 += a_0 * b_3;
		c_13 += a_1 * b_3;
		c_23 += a_2 * b_3;
		c_33 += a_3 * b_3;
		
		A += 16;
		B += 16;

		}
	
	for(; k<kmax; k++)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		b_1 = B[1];
		b_2 = B[2];
		b_3 = B[3];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;

		c_03 += a_0 * b_3;
		c_13 += a_1 * b_3;
		c_23 += a_2 * b_3;
		c_33 += a_3 * b_3;
		
		A += 4;
		B += 4;
		
		}

	if(alg==0)
		{
		C[0+ldc*0] = c_00;
		C[1+ldc*0] = c_10;
		C[2+ldc*0] = c_20;
		C[3+ldc*0] = c_30;

		C[0+ldc*1] = c_01;
		C[1+ldc*1] = c_11;
		C[2+ldc*1] = c_21;
		C[3+ldc*1] = c_31;

		C[0+ldc*2] = c_02;
		C[1+ldc*2] = c_12;
		C[2+ldc*2] = c_22;
		C[3+ldc*2] = c_32;

		C[0+ldc*3] = c_03;
		C[1+ldc*3] = c_13;
		C[2+ldc*3] = c_23;
		C[3+ldc*3] = c_33;
		}
	else if(alg==1)
		{
		C[0+ldc*0] += c_00;
		C[1+ldc*0] += c_10;
		C[2+ldc*0] += c_20;
		C[3+ldc*0] += c_30;

		C[0+ldc*1] += c_01;
		C[1+ldc*1] += c_11;
		C[2+ldc*1] += c_21;
		C[3+ldc*1] += c_31;

		C[0+ldc*2] += c_02;
		C[1+ldc*2] += c_12;
		C[2+ldc*2] += c_22;
		C[3+ldc*2] += c_32;

		C[0+ldc*3] += c_03;
		C[1+ldc*3] += c_13;
		C[2+ldc*3] += c_23;
		C[3+ldc*3] += c_33;
		}
	else
		{
		C[0+ldc*0] -= c_00;
		C[1+ldc*0] -= c_10;
		C[2+ldc*0] -= c_20;
		C[3+ldc*0] -= c_30;

		C[0+ldc*1] -= c_01;
		C[1+ldc*1] -= c_11;
		C[2+ldc*1] -= c_21;
		C[3+ldc*1] -= c_31;

		C[0+ldc*2] -= c_02;
		C[1+ldc*2] -= c_12;
		C[2+ldc*2] -= c_22;
		C[3+ldc*2] -= c_32;

		C[0+ldc*3] -= c_03;
		C[1+ldc*3] -= c_13;
		C[2+ldc*3] -= c_23;
		C[3+ldc*3] -= c_33;
		}

	}



// normal-transposed, 4x3 with data packed in 4
void kernel_sgemm_pp_nt_4x3_lib4(int kmax, float *A, float *B, float *C, int ldc, int alg)
	{
	
	if(kmax<=0)
		return;

/*	const int bs = 4;*/

	int k;
	
	float
		a_0, a_1, a_2, a_3,
		b_0, b_1, b_2,
		c_00=0, c_01=0, c_02=0,
		c_10=0, c_11=0, c_12=0,
		c_20=0, c_21=0, c_22=0,
		c_30=0, c_31=0, c_32=0;
	
	k = 0;
	for(; k<kmax-3; k+=4)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		b_1 = B[1];
		b_2 = B[2];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;
		
		
		a_0 = A[4];
		a_1 = A[5];
		a_2 = A[6];
		a_3 = A[7];
		
		b_0 = B[4];
		b_1 = B[5];
		b_2 = B[6];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;
		
		
		a_0 = A[8];
		a_1 = A[9];
		a_2 = A[10];
		a_3 = A[11];
		
		b_0 = B[8];
		b_1 = B[9];
		b_2 = B[10];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;
		
		
		a_0 = A[12];
		a_1 = A[13];
		a_2 = A[14];
		a_3 = A[15];
		
		b_0 = B[12];
		b_1 = B[13];
		b_2 = B[14];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;
		
		A += 16;
		B += 16;

		}
	
	for(; k<kmax; k++)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		b_1 = B[1];
		b_2 = B[2];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;

		c_02 += a_0 * b_2;
		c_12 += a_1 * b_2;
		c_22 += a_2 * b_2;
		c_32 += a_3 * b_2;
		
		A += 4;
		B += 4;
		
		}

	if(alg==0)
		{
		C[0+ldc*0] = c_00;
		C[1+ldc*0] = c_10;
		C[2+ldc*0] = c_20;
		C[3+ldc*0] = c_30;

		C[0+ldc*1] = c_01;
		C[1+ldc*1] = c_11;
		C[2+ldc*1] = c_21;
		C[3+ldc*1] = c_31;

		C[0+ldc*2] = c_02;
		C[1+ldc*2] = c_12;
		C[2+ldc*2] = c_22;
		C[3+ldc*2] = c_32;
		}
	else if(alg==1)
		{
		C[0+ldc*0] += c_00;
		C[1+ldc*0] += c_10;
		C[2+ldc*0] += c_20;
		C[3+ldc*0] += c_30;

		C[0+ldc*1] += c_01;
		C[1+ldc*1] += c_11;
		C[2+ldc*1] += c_21;
		C[3+ldc*1] += c_31;

		C[0+ldc*2] += c_02;
		C[1+ldc*2] += c_12;
		C[2+ldc*2] += c_22;
		C[3+ldc*2] += c_32;
		}
	else
		{
		C[0+ldc*0] -= c_00;
		C[1+ldc*0] -= c_10;
		C[2+ldc*0] -= c_20;
		C[3+ldc*0] -= c_30;

		C[0+ldc*1] -= c_01;
		C[1+ldc*1] -= c_11;
		C[2+ldc*1] -= c_21;
		C[3+ldc*1] -= c_31;

		C[0+ldc*2] -= c_02;
		C[1+ldc*2] -= c_12;
		C[2+ldc*2] -= c_22;
		C[3+ldc*2] -= c_32;
		}

	}



// normal-transposed, 4x2 with data packed in 4
void kernel_sgemm_pp_nt_4x2_lib4(int kmax, float *A, float *B, float *C, int ldc, int alg)
	{
	
	if(kmax<=0)
		return;

/*	const int ldc = 4;*/

	int k;
	
	float
		a_0, a_1, a_2, a_3,
		b_0, b_1,
		c_00=0, c_01=0,
		c_10=0, c_11=0,
		c_20=0, c_21=0,
		c_30=0, c_31=0;
	
	k = 0;
	for(; k<kmax-3; k+=4)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		b_1 = B[1];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;
		
		
		a_0 = A[4];
		a_1 = A[5];
		a_2 = A[6];
		a_3 = A[7];
		
		b_0 = B[4];
		b_1 = B[5];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;
		
		
		a_0 = A[8];
		a_1 = A[9];
		a_2 = A[10];
		a_3 = A[11];
		
		b_0 = B[8];
		b_1 = B[9];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;
		
		
		a_0 = A[12];
		a_1 = A[13];
		a_2 = A[14];
		a_3 = A[15];
		
		b_0 = B[12];
		b_1 = B[13];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;
		
		A += 16;
		B += 16;

		}
	
	for(; k<kmax; k++)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		b_1 = B[1];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		c_01 += a_0 * b_1;
		c_11 += a_1 * b_1;
		c_21 += a_2 * b_1;
		c_31 += a_3 * b_1;
		
		A += 4;
		B += 4;
		
		}

	if(alg==0)
		{
		C[0+ldc*0] = c_00;
		C[1+ldc*0] = c_10;
		C[2+ldc*0] = c_20;
		C[3+ldc*0] = c_30;

		C[0+ldc*1] = c_01;
		C[1+ldc*1] = c_11;
		C[2+ldc*1] = c_21;
		C[3+ldc*1] = c_31;
		}
	else if(alg==1)
		{
		C[0+ldc*0] += c_00;
		C[1+ldc*0] += c_10;
		C[2+ldc*0] += c_20;
		C[3+ldc*0] += c_30;

		C[0+ldc*1] += c_01;
		C[1+ldc*1] += c_11;
		C[2+ldc*1] += c_21;
		C[3+ldc*1] += c_31;
		}
	else
		{
		C[0+ldc*0] -= c_00;
		C[1+ldc*0] -= c_10;
		C[2+ldc*0] -= c_20;
		C[3+ldc*0] -= c_30;

		C[0+ldc*1] -= c_01;
		C[1+ldc*1] -= c_11;
		C[2+ldc*1] -= c_21;
		C[3+ldc*1] -= c_31;
		}

	}



// normal-transposed, 4x1 with data packed in 4
void kernel_sgemm_pp_nt_4x1_lib4(int kmax, float *A, float *B, float *C, int ldc, int alg)
	{
	
	if(kmax<=0)
		return;

/*	const int ldc = 4;*/

	int k;
	
	float
		a_0, a_1, a_2, a_3,
		b_0,
		c_00=0,
		c_10=0,
		c_20=0,
		c_30=0;
	
	k = 0;
	for(; k<kmax-3; k+=4)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		
		a_0 = A[4];
		a_1 = A[5];
		a_2 = A[6];
		a_3 = A[7];
		
		b_0 = B[4];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		
		a_0 = A[8];
		a_1 = A[9];
		a_2 = A[10];
		a_3 = A[11];
		
		b_0 = B[8];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		
		a_0 = A[12];
		a_1 = A[13];
		a_2 = A[14];
		a_3 = A[15];
		
		b_0 = B[12];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		A += 16;
		B += 16;

		}
	
	for(; k<kmax; k++)
		{

		a_0 = A[0];
		a_1 = A[1];
		a_2 = A[2];
		a_3 = A[3];
		
		b_0 = B[0];
		
		c_00 += a_0 * b_0;
		c_10 += a_1 * b_0;
		c_20 += a_2 * b_0;
		c_30 += a_3 * b_0;
		
		A += 4;
		B += 4;
		
		}

	if(alg==0)
		{
		C[0+ldc*0] = c_00;
		C[1+ldc*0] = c_10;
		C[2+ldc*0] = c_20;
		C[3+ldc*0] = c_30;
		}
	else if(alg==1)
		{
		C[0+ldc*0] += c_00;
		C[1+ldc*0] += c_10;
		C[2+ldc*0] += c_20;
		C[3+ldc*0] += c_30;
		}
	else
		{
		C[0+ldc*0] -= c_00;
		C[1+ldc*0] -= c_10;
		C[2+ldc*0] -= c_20;
		C[3+ldc*0] -= c_30;
		}

	}

