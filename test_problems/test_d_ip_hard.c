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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) || defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || defined(TARGET_AMD_SSE3)
#include <xmmintrin.h> // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

#include "../include/aux_d.h"
#include "../include/aux_s.h"
#include "../include/lqcp_solvers.h"
#include "../include/mpc_solvers.h"
#include "../problem_size.h"
#include "../include/block_size.h"
#include "tools.h"
#include "test_param.h"



/************************************************ 
Mass-spring system: nx/2 masses connected each other with springs (in a row), and the first and the last one to walls. nu (<=nx) controls act on the first nu masses. The system is sampled with sampling time Ts. 
************************************************/
void mass_spring_system(double Ts, int nx, int nu, int N, double *A, double *B, double *b, double *x0)
	{

	int nx2 = nx*nx;

	int info = 0;

	int pp = nx/2; // number of masses
	
/************************************************
* build the continuous time system 
************************************************/
	
	double *T; d_zeros(&T, pp, pp);
	int ii;
	for(ii=0; ii<pp; ii++) T[ii*(pp+1)] = -2;
	for(ii=0; ii<pp-1; ii++) T[ii*(pp+1)+1] = 1;
	for(ii=1; ii<pp; ii++) T[ii*(pp+1)-1] = 1;

	double *Z; d_zeros(&Z, pp, pp);
	double *I; d_zeros(&I, pp, pp); for(ii=0; ii<pp; ii++) I[ii*(pp+1)]=1.0; // = eye(pp);
	double *Ac; d_zeros(&Ac, nx, nx);
	dmcopy(pp, pp, Z, pp, Ac, nx);
	dmcopy(pp, pp, T, pp, Ac+pp, nx);
	dmcopy(pp, pp, I, pp, Ac+pp*nx, nx);
	dmcopy(pp, pp, Z, pp, Ac+pp*(nx+1), nx); 
	free(T);
	free(Z);
	free(I);
	
	d_zeros(&I, nu, nu); for(ii=0; ii<nu; ii++) I[ii*(nu+1)]=1.0; //I = eye(nu);
	double *Bc; d_zeros(&Bc, nx, nu);
	dmcopy(nu, nu, I, nu, Bc+pp, nx);
	free(I);
	
/************************************************
* compute the discrete time system 
************************************************/

	double *bb; d_zeros(&bb, nx, 1);
	dmcopy(nx, 1, bb, nx, b, nx);
		
	dmcopy(nx, nx, Ac, nx, A, nx);
	dscal_3l(nx2, Ts, A);
	expm(nx, A);
	
	d_zeros(&T, nx, nx);
	d_zeros(&I, nx, nx); for(ii=0; ii<nx; ii++) I[ii*(nx+1)]=1.0; //I = eye(nx);
	dmcopy(nx, nx, A, nx, T, nx);
	daxpy_3l(nx2, -1.0, I, T);
	dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);
	
	int *ipiv = (int *) malloc(nx*sizeof(int));
	dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
	free(ipiv);

	free(Ac);
	free(Bc);
	free(bb);
	
			
/************************************************
* initial state 
************************************************/
	
	if(nx==4)
		{
		x0[0] = 5;
		x0[1] = 10;
		x0[2] = 15;
		x0[3] = 20;
		}
	else
		{
		int jj;
		for(jj=0; jj<nx; jj++)
			x0[jj] = 1;
		}

	}



int main()
	{
	
	printf("\n");
	printf("\n");
	printf("\n");
	printf(" HPMPC -- Library for High-Performance implementation of solvers for MPC.\n");
	printf(" Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.\n");
	printf("\n");
	printf(" HPMPC is distributed in the hope that it will be useful,\n");
	printf(" but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf(" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
	printf(" See the GNU Lesser General Public License for more details.\n");
	printf("\n");
	printf("\n");
	printf("\n");

#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) || defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || defined(TARGET_AMD_SSE3)
/*	printf("\nflush subnormals to zero\n");*/
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON); // flush to zero subnormals !!! works only with one thread !!!
#endif

	int ii, jj, idx;
	
	int rep, nrep=NREP;

	int nx = NX; // number of states (it has to be even for the mass-spring system test problem)
	int nu = NU; // number of inputs (controllers) (it has to be at least 1 and at most nx/2 for the mass-spring system test problem)
	int N  = NN; // horizon lenght
	int nb = NB; // number of box constrained inputs and states
	int ng = 0; //4;  // number of general constraints

	int nbu = nu<nb ? nu : nb ;
	int ngu = nu<ng ? nu : ng ; // TODO remove when not needed any longer in tests

	printf(" Test problem: mass-spring system with %d masses and %d controls.\n", nx/2, nu);
	printf("\n");
	printf(" MPC problem size: %d states, %d inputs, %d horizon length, %d two-sided box constraints, %d two-sided general constraints.\n", nx, nu, N, nb, ng);
	printf("\n");
#if IP == 1
	printf(" IP method parameters: primal-dual IP, double precision, %d maximum iterations, %5.1e exit tolerance in duality measure (edit file test_param.c to change them).\n", K_MAX, MU_TOL);
#elif IP == 2
	printf(" IP method parameters: predictor-corrector IP, double precision, %d maximum iterations, %5.1e exit tolerance in duality measure (edit file test_param.c to change them).\n", K_MAX, MU_TOL);
#else
	printf(" Wrong value for IP solver choice: %d\n", IP);
#endif

	int info = 0;
		
	const int bs = D_MR; //d_get_mr();
	const int ncl = D_NCL;
	const int nal = bs*ncl; // number of doubles per cache line
	
	const int nz = nx+nu+1;
	const int pnz = bs*((nz+bs-1)/bs);
	const int pnx = bs*((nx+bs-1)/bs);
	const int cnz = ncl*((nx+nu+1+ncl-1)/ncl);
	const int cnx = ncl*((nx+ncl-1)/ncl);
	const int cng = ncl*((ng+ncl-1)/ncl);
	//const int pnb = bs*((2*nb+bs-1)/bs); // packed number of box constraints
	const int pnb = bs*((nb+bs-1)/bs); // SIMD aligned number of one-sided box constraints !!!!!!!!!!!!
	const int anz = nal*((nz+nal-1)/nal);
	const int anx = nal*((nx+nal-1)/nal);
//	const int anb = nal*((2*nb+nal-1)/nal); // cache aligned number of box constraints
	const int anb = nal*((nb+nal-1)/nal); // cache aligned number of one-sided box constraints !!!!!!!!!!!! // TODO ALIGN TO SIMD-ALIGNED BOUNDARIES
	const int ang = nal*((ng+nal-1)/nal); // cache aligned number of one-sided box constraints !!!!!!!!!!!! // TODO ALIGN TO SIMD-ALIGNED BOUNDARIES

//	const int pad = (ncl-nx%ncl)%ncl; // packing between BAbtL & P
//	const int cnl = cnz<cnx+ncl ? nx+pad+cnx+ncl : nx+pad+cnz;
	const int cnl = cnz<cnx+ncl ? cnx+ncl : cnz;
	
/************************************************
* dynamical system
************************************************/	

	double *A; d_zeros(&A, nx, nx); // states update matrix

	double *B; d_zeros(&B, nx, nu); // inputs matrix

	double *b; d_zeros(&b, nx, 1); // states offset
	double *x0; d_zeros(&x0, nx, 1); // initial state

	double Ts = 0.5; // sampling time
	mass_spring_system(Ts, nx, nu, N, A, B, b, x0);
	
	for(jj=0; jj<nx; jj++)
		b[jj] = 0.1;
	
	for(jj=0; jj<nx; jj++)
		x0[jj] = 0;
	x0[0] = 3.5;
	x0[1] = 3.5;
	
//	d_print_mat(nx, nx, A, nx);
//	d_print_mat(nx, nu, B, nx);
//	d_print_mat(nx, 1, b, nx);
//	d_print_mat(nx, 1, x0, nx);
	
	/* packed */
/*	double *BAb; d_zeros(&BAb, nx, nz);*/

/*	dmcopy(nx, nu, B, nx, BAb, nx);*/
/*	dmcopy(nx, nx, A, nx, BAb+nu*nx, nx);*/
/*	dmcopy(nx, 1 , b, nx, BAb+(nu+nx)*nx, nx);*/
	
	/* transposed */
/*	double *BAbt; d_zeros_align(&BAbt, pnz, pnz);*/
/*	for(ii=0; ii<nx; ii++)*/
/*		for(jj=0; jj<nz; jj++)*/
/*			{*/
/*			BAbt[jj+pnz*ii] = BAb[ii+nx*jj];*/
/*			}*/

	/* packed into contiguous memory */
	double *pBAbt; d_zeros_align(&pBAbt, pnz, cnx);
/*	d_cvt_mat2pmat(nz, nx, 0, bs, BAbt, pnz, pBAbt, cnx);*/
/*	d_cvt_tran_mat2pmat(nx, nz, 0, bs, BAb, nx, pBAbt, cnx);*/

	d_cvt_tran_mat2pmat(nx, nu, 0, bs, B, nx, pBAbt, cnx);
	d_cvt_tran_mat2pmat(nx, nx, nu, bs, A, nx, pBAbt+nu/bs*cnx*bs+nu%bs, cnx);
	for (jj = 0; jj<nx; jj++)
		pBAbt[(nx+nu)/bs*cnx*bs+(nx+nu)%bs+jj*bs] = b[jj];

/*	d_print_pmat (nz, nx, bs, pBAbt, cnx);*/
/*	exit(1);*/

/************************************************
* box constraints
************************************************/	

	double *db; d_zeros_align(&db, 2*pnb, 1);
	for(jj=0; jj<nbu; jj++)
		{
		db[jj]     = - 0.5;   //   umin
		db[pnb+jj] = - 0.5;   // - umax
		}
	for(; jj<nb; jj++)
		{
		db[jj]     = - 4.0;   //   xmin
		db[pnb+jj] = - 4.0;   // - xmax
		}

/************************************************
* general constraints
************************************************/	

	double *dg; d_zeros_align(&dg, 2*ng, 1);
	for(jj=0; jj<2*ngu; jj++)
		dg[jj] = - 0.5;   // umin
	for(; jj<2*ng; jj++)
		dg[jj] = - 4.0;   // xmin
	
	double *D; d_zeros(&D, ng, nu+nx);
	for(jj=0; jj<ng; jj++)
		D[jj*(ng+1)] = 1.0;
	//d_print_mat(ng, nu+nx, D, ng);

	double *pD; d_zeros_align(&pD, pnz, cng);
	d_cvt_tran_mat2pmat(ng, nu+nx, 0, bs, D, ng, pD, cng);
	//d_print_pmat(nu+nx, ng, bs, pD, cng);
	//exit(1);
	// TODO arrived here working on general constraints

/************************************************
* cost function
************************************************/	

	double *Q; d_zeros_align(&Q, pnz, pnz);
	for(ii=0; ii<nu; ii++) Q[ii*(pnz+1)] = 2.0;
	for(; ii<pnz; ii++) Q[ii*(pnz+1)] = 1.0;
	for(ii=0; ii<nu+nx; ii++) Q[nx+nu+ii*pnz] = 0.1;
/*	Q[(nx+nu)*(pnz+1)] = 1e35; // large enough (not needed any longer) */
	
	/* packed into contiguous memory */
	double *pQ; d_zeros_align(&pQ, pnz, cnz);
	d_cvt_mat2pmat(nz, nz, 0, bs, Q, pnz, pQ, cnz);

	// linear part copied on another vector
	double *q; d_zeros_align(&q, anz, 1);
	for(ii=0; ii<nu+nx; ii++) q[ii] = Q[nx+nu+ii*pnz];

	// maximum element in cost functions
	double mu0 = 1.0;
	for(ii=0; ii<nu+nx; ii++)
		for(jj=0; jj<nu+nx; jj++)
			mu0 = fmax(mu0, Q[jj+nz*ii]);
	//printf("\n mu0 = %f\n", mu0);

/************************************************
* matrices series
************************************************/	

	double *(hpQ[N+1]);
	double *(hq[N+1]);
	double *(hux[N+1]);
	double *(hpi[N+1]);
	double *(hlam[N+1]);
	double *(ht[N+1]);
	double *(hpBAbt[N]);
	double *(hdb[N+1]);
	double *(hrb[N]);
	double *(hrq[N+1]);
	double *(hrd[N+1]);

	for(jj=0; jj<N; jj++)
		{
		//d_zeros_align(&hpQ[jj], pnz, cnz);
		hpQ[jj] = pQ;
		//d_zeros_align(&hq[jj], anz, 1);
		hq[jj] = q;
		d_zeros_align(&hux[jj], anz, 1);
		d_zeros_align(&hpi[jj], anx, 1);
		d_zeros_align(&hlam[jj],2*pnb, 1); // TODO pnb
		d_zeros_align(&ht[jj], 2*pnb, 1); // TODO pnb
		hpBAbt[jj] = pBAbt;
		hdb[jj] = db;
		d_zeros_align(&hrb[jj], anx, 1);
		d_zeros_align(&hrq[jj], anz, 1);
		d_zeros_align(&hrd[jj], 2*pnb, 1); // TODO pnb
		}
	//d_zeros_align(&hpQ[N], pnz, cnz);
	hpQ[N] = pQ;
	//d_zeros_align(&hq[N], anz, 1);
	hq[N] = q;
	d_zeros_align(&hux[N], anz, 1);
	d_zeros_align(&hpi[N], anx, 1);
	d_zeros_align(&hlam[N], 2*pnb, 1); // TODO pnb
	d_zeros_align(&ht[N], 2*pnb, 1); // TODO pnb
	hdb[N] = db;
	d_zeros_align(&hrq[N], anz, 1);
	d_zeros_align(&hrd[N], 2*pnb, 1); // TODO pnb
	
	// starting guess
	for(jj=0; jj<nx; jj++) hux[0][nu+jj]=x0[jj];

/************************************************
* riccati-like iteration
************************************************/

	//double *work; d_zeros_align(&work, (N+1)*(pnz*cnl + 4*anz + 4*anb + 2*anx) + 3*anz, 1); // work space
	//double *work; d_zeros_align(&work, (N+1)*(pnz*cnl + 5*anz + 4*anb + 2*anx) + 3*anz, 1); // work space TODO change work space on other files !!!!!!!!!!!!!
	//double *work; d_zeros_align(&work, (N+1)*(pnz*cnl + 5*anz + 4*anb + 2*anx) + anz + pnz*cnx, 1); // work space TODO change work space on other files !!!!!!!!!!!!!
	double *work; d_zeros_align(&work, (N+1)*(pnz*cnl + 5*anz + 8*pnb + 2*anx) + anz + pnz*cnx, 1); // work space TODO change work space on other files !!!!!!!!!!!!!
/*	for(jj=0; jj<( (N+1)*(pnz*cnl + 4*anz + 4*anb + 2*anx) + 3*anz ); jj++) work[jj] = -1.0;*/
	int kk = 0; // acutal number of iterations
/*	char prec = PREC; // double/single precision*/
/*	double sp_thr = SP_THR; // threshold to switch between double and single precision*/
	int k_max = K_MAX; // maximum number of iterations in the IP method
	double mu_tol = MU_TOL; // tolerance in the duality measure
	double alpha_min = ALPHA_MIN; // minimum accepted step length
	double sigma[] = {0.4, 0.3, 0.01}; // control primal-dual IP behaviour
	double *stat; d_zeros(&stat, 5, k_max); // stats from the IP routine
	int compute_mult = COMPUTE_MULT;
	int warm_start = WARM_START;
	double mu = -1.0;
	int hpmpc_status;
	


	/* initizile the cost function */
//	for(ii=0; ii<N; ii++)
//		{
//		for(jj=0; jj<pnz*cnz; jj++) hpQ[ii][jj]=pQ[jj];
//		}
//	for(jj=0; jj<pnz*cnz; jj++) hpQ[N][jj]=pQ[jj];



	// solution of unconstrained problem
//	double *(hpL[N+1]);
//	double *(hPb[N]);
//	for(jj=0; jj<=N; jj++)
//		{
//		d_zeros_align(&hpL[jj], pnz, cnl);
//		d_zeros_align(&hPb[jj], anx, 1);
//		}
//	double *diag; d_zeros_align(&diag, pnz, 1);
//	int update_hessian = 0;
//	double **Qd;
//	double **Ql;
//	d_ric_sv_mpc(nx, nu, N, hpBAbt, hpQ, update_hessian, Qd, Ql, hux, hpL, work, diag, COMPUTE_MULT, hpi);



	// initial states
	double xx0[] = {3.5, 3.5, 3.66465, 2.15833, 1.81327, -0.94207, 1.86531, -2.35760, 2.91534, 1.79890, -1.49600, -0.76600, -2.60268, 1.92456, 1.66630, -2.28522, 3.12038, 1.83830, 1.93519, -1.87113};



	/* warm up */

	// initialize states and inputs
	for(ii=0; ii<=N; ii++)
		for(jj=0; jj<nx+nu; jj++)
			hux[ii][jj] = 0;

	hux[0][nu+0] = xx0[0];
	hux[0][nu+1] = xx0[1];

//	// solution of unconstrained problem as warm start
//	warm_start = 1;
//	d_ric_trs_mpc(nx, nu, N, hpBAbt, hpL, hq, hux, work, 1, hPb, COMPUTE_MULT, hpi);
//	//d_ric_sv_mpc(nx, nu, N, hpBAbt, hpQ, update_hessian, Qd, Ql, hux, hpL, work, diag, COMPUTE_MULT, hpi);

	// call the IP solver
//	if(FREE_X0==0)
//		{
		if(IP==1)
			hpmpc_status = d_ip_hard_mpc(&kk, k_max, mu0, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, ng, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
		else
			hpmpc_status = d_ip2_hard_mpc(&kk, k_max, mu0, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, ng, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
//		}
//	else
//		{
//		if(IP==1)
//			hpmpc_status = d_ip_box_mhe_old(&kk, k_max, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
//		else
//			hpmpc_status = d_ip2_box_mhe_old(&kk, k_max, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
//		}



	int kk_avg = 0;

	/* timing */
	struct timeval tv0, tv1;
	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{

		idx = rep%10;
		x0[0] = xx0[2*idx];
		x0[1] = xx0[2*idx+1];

		// initialize states and inputs
		for(ii=0; ii<=N; ii++)
			for(jj=0; jj<nx+nu; jj++)
				hux[ii][jj] = 0;

		hux[0][nu+0] = xx0[2*idx];
		hux[0][nu+1] = xx0[2*idx+1];

//		// solution of unconstrained problem as warm start
//		warm_start = 1;
//		d_ric_trs_mpc(nx, nu, N, hpBAbt, hpL, hq, hux, work, 1, hPb, COMPUTE_MULT, hpi);
//		//d_ric_sv_mpc(nx, nu, N, hpBAbt, hpQ, update_hessian, Qd, Ql, hux, hpL, work, diag, COMPUTE_MULT, hpi);

		// call the IP solver
//		if(FREE_X0==0)
//			{
			if(IP==1)
				hpmpc_status = d_ip_hard_mpc(&kk, k_max, mu0, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, ng, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
			else
				hpmpc_status = d_ip2_hard_mpc(&kk, k_max, mu0, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, ng, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
//			}
//		else
//			{
//			if(IP==1)
//				hpmpc_status = d_ip_box_mhe_old(&kk, k_max, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
//			else
//				hpmpc_status = d_ip2_box_mhe_old(&kk, k_max, mu_tol, alpha_min, warm_start, sigma, stat, nx, nu, N, nb, hpBAbt, hpQ, hdb, hux, compute_mult, hpi, hlam, ht, work);
//			}

		kk_avg += kk;

		}
	
	gettimeofday(&tv1, NULL); // stop
	


	double time = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);
	
/*	printf("\nnx\tnu\tN\tkernel\n\n");*/
/*	printf("\n%d\t%d\t%d\t%e\n\n", nx, nu, N, time);*/
	
	printf("\n");
	printf(" Average number of iterations over %d runs: %5.1f\n", nrep, kk_avg / (double) nrep);
	printf("\n");
	printf(" Average solution time over %d runs: %5.2e seconds\n", nrep, time);
	printf("\n");



	// restore linear part of cost function 
	for(ii=0; ii<N; ii++)
		{
		for(jj=0; jj<nx+nu; jj++) hq[ii][jj] = Q[nx+nu+pnz*jj];
		}
	for(jj=0; jj<nx+nu; jj++) hq[N][jj] = Q[nx+nu+pnz*jj];

	// residuals computation
//	if(FREE_X0==0)
		d_res_ip_hard_mpc(nx, nu, N, nb, ng, hpBAbt, hpQ, hq, hux, hdb, hpi, hlam, ht, hrq, hrb, hrd, &mu);
//	else
//		d_res_ip_box_mhe_old(nx, nu, N, nb, hpBAbt, hpQ, hq, hux, hdb, hpi, hlam, ht, hrq, hrb, hrd, &mu);


	if(PRINTSTAT==1)
		{

		printf("\n");
		printf("\n");
		printf(" Print IP statistics of the last run\n");
		printf("\n");

		for(jj=0; jj<kk; jj++)
			printf("k = %d\tsigma = %f\talpha = %f\tmu = %f\t\tmu = %e\talpha = %f\tmu = %f\tmu = %e\n", jj, stat[5*jj], stat[5*jj+1], stat[5*jj+2], stat[5*jj+2], stat[5*jj+3], stat[5*jj+4], stat[5*jj+4]);
		printf("\n");
		
		}

	if(PRINTRES==1)
		{

		printf("\n");
		printf("\n");
		printf(" Print solution\n");
		printf("\n");

		printf("\nu = \n\n");
		for(ii=0; ii<N; ii++)
			d_print_mat(1, nu, hux[ii], 1);
		
		printf("\nx = \n\n");
		for(ii=0; ii<=N; ii++)
			d_print_mat(1, nx, hux[ii]+nu, 1);
		
		printf("\nlam = \n\n");
		for(ii=0; ii<=N; ii++)
			d_print_mat(1, 2*pnb, hlam[ii], 1);
		
		}

	if(PRINTRES==1 && COMPUTE_MULT==1)
		{
		// print result 
		// print result 
		printf("\n");
		printf("\n");
		printf(" Print residuals\n\n");
		printf("\n");
		printf("\n");
		printf("rq = \n\n");
//		if(FREE_X0==0)
//			{
			d_print_mat(1, nu, hrq[0], 1);
			for(ii=1; ii<=N; ii++)
/*				d_print_mat_e(1, nx+nu, hrq[ii], 1);*/
				d_print_mat(1, nx+nu, hrq[ii], 1);
//			}
//		else
//			{
//			for(ii=0; ii<=N; ii++)
///*				d_print_mat_e(1, nx+nu, hrq[ii], 1);*/
//				d_print_mat(1, nx+nu, hrq[ii], 1);
//			}
		printf("\n");
		printf("\n");
		printf("rb = \n\n");
		for(ii=0; ii<N; ii++)
/*			d_print_mat_e(1, nx, hrb[ii], 1);*/
			d_print_mat(1, nx, hrb[ii], 1);
		printf("\n");
		printf("\n");
		printf("rd = \n\n");
		for(ii=0; ii<=N; ii++)
/*			d_print_mat_e(1, 2*nb, hrd[ii], 1);*/
			d_print_mat(1, 2*nb, hrd[ii], 1);
		printf("\n");
		printf("\n");
		printf("mu = %e\n\n", mu);
		
		}

/*	printf("\nnx\tnu\tN\tkernel\n\n");*/
/*	printf("\n%d\t%d\t%d\t%e\n\n", nx, nu, N, time);*/
	
/************************************************
* free memory and return
************************************************/

	free(A);
	free(B);
	free(b);
	free(x0);
/*	free(BAb);*/
/*	free(BAbt);*/
	free(pBAbt);
	free(db);
	free(Q);
	free(pQ);
	free(q);
	free(work);
	free(stat);
	for(jj=0; jj<N; jj++)
		{
		//free(hpQ[jj]);
		//free(hq[jj]);
		free(hux[jj]);
		free(hpi[jj]);
		free(hlam[jj]);
		free(ht[jj]);
		free(hrb[jj]);
		free(hrq[jj]);
		free(hrd[jj]);
		}
	//free(hpQ[N]);
	//free(hq[N]);
	free(hux[N]);
	free(hpi[N]);
	free(hlam[N]);
	free(ht[N]);
	free(hrq[N]);
	free(hrd[N]);



	return 0;

	}



