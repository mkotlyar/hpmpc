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

#include <math.h>

#include "../include/aux_d.h"
#include "../include/aux_s.h"
#include "../include/lqcp_solvers.h"
#include "../include/block_size.h"
#include "../include/mpc_aux.h"






/* primal-dual interior-point method, box constraints, time invariant matrices (mpc version) */
int d_ip_soft_mpc(int *kk, int k_max, double mu_tol, double alpha_min, int warm_start, double *sigma_par, double *stat, int nx, int nu, int N, int nb, double **pBAbt, double **pQ, double **Z, double **z, double **db, double **ux, int compute_mult, double **pi, double **lam, double **t, double *work_memory)
	{
	
	// TODO Z and z are supposed to be aligned with the corresponding components of t and lam (i.e. have 2*nu dummy elements at the beginning)

/*printf("\ncazzo\n");*/

/*	int nbx = nb - nu;*/
/*	if(nbx<0)*/
/*		nbx = 0;*/
	int nbu = nu<nb ? nu : nb ;
	int nbx = nb-nu>0 ? nb-nu : 0 ;

	// indeces
	int jj, ll, ii, bs0;


	// constants
	const int bs = D_MR; //d_get_mr();
	const int ncl = D_NCL;
	const int nal = bs*ncl; // number of doubles per cache line

	const int nz = nx+nu+1;
	const int nxu = nx+nu;
	const int pnz = bs*((nz+bs-1)/bs);
	const int pnx = bs*((nx+bs-1)/bs);
	const int cnz = ncl*((nz+ncl-1)/ncl);
	const int cnx = ncl*((nx+ncl-1)/ncl);
	const int anz = nal*((nz+nal-1)/nal);
	const int anx = nal*((nx+nal-1)/nal);
	const int anb = nal*((2*nb+nal-1)/nal); // cache aligned number of box and soft constraints

	const int pad = (ncl-nx%ncl)%ncl; // packing between BAbtL & P
	const int cnl = cnz<cnx+ncl ? nx+pad+cnx+ncl : nx+pad+cnz;
	
	// initialize work space
	double *ptr;
	ptr = work_memory;

	double *(dux[N+1]);
	double *(dpi[N+1]);
	double *(pL[N+1]);
	double *(pd[N+1]); // pointer to diagonal of Hessian
	double *(pl[N+1]); // pointer to linear part of Hessian
	double *(pl2[N+1]); // pointer to linear part of Hessian (backup)
	double *(bd[N+1]); // backup diagonal of Hessian
	double *(bl[N+1]); // backup linear part of Hessian
	double *work;
	double *diag;
	double *(dlam[N+1]);
	double *(dt[N+1]);
	double *(lamt[N+1]);
	double *(t_inv[N+1]);
	double *(Zl[N+1]); // inverse of the diagonal of the matrix of the cost funciton of the soft constraint slack variables as updated in the IP
	double *(zl[N+1]); // linear part of the cost funciton of the soft constraint slack variables as updated in the IP

//	ptr += (N+1)*(pnz + pnx + pnz*cnl + 8*pnz) + 3*pnz;

	// inputs and states
	for(jj=0; jj<=N; jj++)
		{
		dux[jj] = ptr;
		ptr += anz;
		}

	// equality constr multipliers
	for(jj=0; jj<=N; jj++)
		{
		dpi[jj] = ptr;
		ptr += anx;
		}
	
	// Hessian
	for(jj=0; jj<=N; jj++)
		{
		pd[jj] = pQ[jj];
		pl[jj] = pQ[jj] + ((nu+nx)/bs)*bs*cnz + (nu+nx)%bs;
		bd[jj] = ptr;
		bl[jj] = ptr + anz;
		ptr += 2*anz;
		// backup
		for(ll=0; ll<nx+nu; ll++)
			{
			bd[jj][ll] = pQ[jj][(ll/bs)*bs*cnz+ll%bs+ll*bs];
			bl[jj][ll] = pQ[jj][((nx+nu)/bs)*bs*cnz+(nx+nu)%bs+ll*bs];
			}
		}

	// work space
	for(jj=0; jj<=N; jj++)
		{
		pL[jj] = ptr;
		ptr += pnz*cnl;
		}
	
	for(jj=0; jj<=N; jj++)
		{
		pl2[jj] = ptr;
		ptr += anz;
		}

	work = ptr;
	ptr += 2*anz;

	diag = ptr;
	ptr += anz;

	// slack variables, Lagrangian multipliers for inequality constraints and work space (assume # box constraints <= 2*(nx+nu) < 2*pnz)
	for(jj=0; jj<=N; jj++)
		{
		dlam[jj] = ptr;
		dt[jj]   = ptr + 2*anb;;
		ptr += 4*anb;
		}
	for(jj=0; jj<=N; jj++)
		{
		lamt[jj] = ptr;
		ptr += 2*anb;
		}
	for(jj=0; jj<=N; jj++)
		{
		t_inv[jj] = ptr;
		ptr += 2*anb;
		}

	// updated cost function of soft constraint slack variables
	for(jj=0; jj<=N; jj++)
		{
		Zl[jj] = ptr;
		zl[jj] = ptr + anb;;
		ptr += 2*anb;
		}


	double temp0, temp1;
	double alpha, mu;
	double mu_scal = 1.0/(N*2*(nb+nbx));
	double sigma, sigma_decay, sigma_min;

	sigma = sigma_par[0]; //0.4;
	sigma_decay = sigma_par[1]; //0.3;
	sigma_min = sigma_par[2]; //0.01;
	


	// initialize ux & t>0 (slack variable)
	d_init_var_soft_mpc(N, nx, nu, nb, ux, pi, db, t, lam, warm_start);



	// initialize pi
	for(jj=0; jj<=N; jj++)
		for(ll=0; ll<nx; ll++)
			dpi[jj][ll] = 0.0;



	// initialize dux
	for(ll=0; ll<nx; ll++)
		dux[0][nu+ll] = ux[0][nu+ll];



	// compute the duality gap
	alpha = 0.0; // needed to compute mu !!!!!
	d_compute_mu_soft_mpc(N, nx, nu, nb, &mu, mu_scal, alpha, lam, dlam, t, dt);

	// set to zero iteration count
	*kk = 0;	

	// larger than minimum accepted step size
	alpha = 1.0;
	


	// IP loop		
	while( *kk<k_max && mu>mu_tol && alpha>=alpha_min )
		{
						
/*printf("\niteration %d, mu %f\n", *kk, mu);*/
		//update cost function matrices and vectors (box constraints)

		// box constraints

		// update hessian
		d_update_hessian_soft_mpc(N, nx, nu, nb, cnz, sigma*mu, t, t_inv, lam, lamt, dlam, bd, bl, pd, pl, pl2, db, Z, z, Zl, zl);

/*return;*/


		// compute the search direction: factorize and solve the KKT system
		d_ric_sv_mpc(nx, nu, N, pBAbt, pQ, dux, pL, work, diag, compute_mult, dpi);



		// compute t_aff & dlam_aff & dt_aff & alpha
		alpha = 1.0;
		d_compute_alpha_soft_mpc(N, nx, nu, nb, &alpha, t, dt, lam, dlam, lamt, dux, db, Zl, zl);

		stat[5*(*kk)] = sigma;
		stat[5*(*kk)+1] = alpha;
			
		alpha *= 0.995;




		// update x, u, lam, t & compute the duality gap mu
		d_update_var_soft_mpc(N, nx, nu, nb, &mu, mu_scal, alpha, ux, dux, t, dt, lam, dlam, pi, dpi);
		
		stat[5*(*kk)+2] = mu;
		


		// update sigma
		sigma *= sigma_decay;
		if(sigma<sigma_min)
			sigma = sigma_min;
		
		if(alpha<0.3)
			sigma = sigma_par[0];



		// increment loop index
		(*kk)++;



		} // end of IP loop
	

	// restore Hessian
	for(jj=0; jj<=N; jj++)
		{
		for(ll=0; ll<nx+nu; ll++)
			{
			pQ[jj][(ll/bs)*bs*cnz+ll%bs+ll*bs] = bd[jj][ll];
			pQ[jj][((nx+nu)/bs)*bs*cnz+(nx+nu)%bs+ll*bs] = bl[jj][ll];
			}
		}

/*printf("\nfinal iteration %d, mu %f\n", *kk, mu);*/



	// successful exit
	if(mu<=mu_tol)
		return 0;
	
	// max number of iterations reached
	if(*kk>=k_max)
		return 1;
	
	// no improvement
	if(alpha<alpha_min)
		return 2;
	
	// impossible
	return -1;

	} // end of ipsolver

