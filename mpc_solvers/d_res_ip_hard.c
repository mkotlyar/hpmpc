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


#include "../include/blas_d.h"
#include "../include/block_size.h"



/* supports the problem size to change stage-wise */
void d_res_mpc_hard_tv(int N, int const * nx, int const * nu, int const * nb, int **idxb, int const * ng, double **hpBAbt, double **hb, double **hpQ, double **hq, double **hux, double **hpDCt, double **hd, double **hpi, double **hlam, double **ht, double **hrq, double **hrb, double **hrd, double *mu)
	{

	const int bs = D_MR;
	const int ncl = D_NCL;

	static double temp[D_MR] = {};

	int ii, jj;
	
	int nu0, nu1, cnux0, nx0, nx1, nxm, cnx0, cnx1, nb0, pnb, ng0, png, cng, nb_tot;


	// initialize mu
	nb_tot = 0;
	mu[0] = 0;



	// first block
	ii = 0;
	nu0 = nu[ii];
	nu1 = nu[ii+1];
	nx0 = nx[ii]; // nx1;
	nx1 = nx[ii+1];
	cnx1  = (nx1+ncl-1)/ncl*ncl;
	cnux0 = (nu0+nx0+ncl-1)/ncl*ncl;
	nb0 = nb[ii];
	pnb = (nb0+bs-1)/bs*bs;
	ng0 = ng[ii];
	png = (ng0+bs-1)/bs*bs;
	cng = (ng0+ncl-1)/ncl*ncl;
	nb_tot += nb0 + ng0;

	for(jj=0; jj<nb0; jj++)
		mu[0] += hlam[ii][jj] * ht[ii][jj] + hlam[ii][pnb+jj] * ht[ii][pnb+jj];
	for(jj=0; jj<ng0; jj++) 
		mu[0] += hlam[ii][2*pnb+jj] * ht[ii][2*pnb+jj] + hlam[ii][2*pnb+png+jj] * ht[ii][2*pnb+png+jj];

	for(jj=0; jj<nb0; jj++)
		{
		hrd[ii][jj]     =   hux[ii][idxb[ii][jj]] - hd[ii][jj]     - ht[ii][jj];
		hrd[ii][pnb+jj] = - hux[ii][idxb[ii][jj]] + hd[ii][pnb+jj] - ht[ii][pnb+jj];
		}
	if(ng0>0)
		{
		dgemv_t_lib(nu0+nx0, ng0, hpDCt[ii], cng, hux[ii], 0, hrd[ii]+2*pnb, hrd[ii]+2*pnb);
		for(jj=0; jj<ng0; jj++)
			{
			hrd[ii][2*pnb+png+jj] = - hrd[ii][2*pnb+jj];
			hrd[ii][2*pnb+jj] += - hd[ii][2*pnb+jj] - ht[ii][2*pnb+jj];
			hrd[ii][2*pnb+png+jj] += hd[ii][2*pnb+png+jj] - ht[ii][2*pnb+png+jj];
			}
		}

//	for(jj=0; jj<nu0; jj++) 
//		hrq[ii][jj] = - hq[ii][jj];
//	for(jj=0; jj<nx0; jj++) 
//		hrq[ii][nu0+jj] = - hq[ii][nu0+jj]; // + hpi[ii-1][jj];
	for(jj=0; jj<nu0+nx0; jj++) 
		hrq[ii][jj] = - hq[ii][jj];
	for(jj=0; jj<nb0; jj++) 
		hrq[ii][idxb[ii][jj]] += hlam[ii][jj] - hlam[ii][pnb+jj];
	dsymv_lib(nu0+nx0, nu0+nx0, hpQ[ii], cnux0, hux[ii], -1, hrq[ii], hrq[ii]);
	if(ng0>0)
		{
		// TODO work space + one dgemv call
		dgemv_n_lib(nu0+nx0, ng0, hpDCt[ii], cng, hlam[ii]+2*pnb, 1, hrq[ii], hrq[ii]);
		dgemv_n_lib(nu0+nx0, ng0, hpDCt[ii], cng, hlam[ii]+2*pnb+png, -1, hrq[ii], hrq[ii]);
		}

	for(jj=0; jj<nx1; jj++) 
		hrb[ii][jj] = hux[ii+1][nu1+jj] - hb[ii][jj];
	dgemv_nt_lib(nu0+nx0, nx1, hpBAbt[ii], cnx1, hpi[ii], hux[ii], -1, -1, hrq[ii], hrb[ii], hrq[ii], hrb[ii]);



	// middle blocks
	for(ii=1; ii<N; ii++)
		{
		nu0 = nu1;
		nu1 = nu[ii+1];
		nx0 = nx1;
		nx1 = nx[ii+1];
		cnx0 = cnx1;
		cnx1  = (nx1+ncl-1)/ncl*ncl;
		cnux0  = (nu0+nx0+ncl-1)/ncl*ncl;
		nb0 = nb[ii];
		pnb = (nb0+bs-1)/bs*bs;
		ng0 = ng[ii];
		png = (ng0+bs-1)/bs*bs;
		cng = (ng0+ncl-1)/ncl*ncl;
		nb_tot += nb0 + ng0;

		for(jj=0; jj<nb0; jj++)
			mu[0] += hlam[ii][jj] * ht[ii][jj] + hlam[ii][pnb+jj] * ht[ii][pnb+jj];
		for(jj=0; jj<ng0; jj++) 
			mu[0] += hlam[ii][2*pnb+jj] * ht[ii][2*pnb+jj] + hlam[ii][2*pnb+png+jj] * ht[ii][2*pnb+png+jj];

		for(jj=0; jj<nb0; jj++)
			{
			hrd[ii][jj]     =   hux[ii][idxb[ii][jj]] - hd[ii][jj]     - ht[ii][jj];
			hrd[ii][pnb+jj] = - hux[ii][idxb[ii][jj]] + hd[ii][pnb+jj] - ht[ii][pnb+jj];
			}
		if(ng0>0)
			{
			dgemv_t_lib(nu0+nx0, ng0, hpDCt[ii], cng, hux[ii], 0, hrd[ii]+2*pnb, hrd[ii]+2*pnb);
			for(jj=0; jj<ng0; jj++)
				{
				hrd[ii][2*pnb+png+jj] = - hrd[ii][2*pnb+jj];
				hrd[ii][2*pnb+jj] += - hd[ii][2*pnb+jj] - ht[ii][2*pnb+jj];
				hrd[ii][2*pnb+png+jj] += hd[ii][2*pnb+png+jj] - ht[ii][2*pnb+png+jj];
				}
			}

		for(jj=0; jj<nu0; jj++) 
			hrq[ii][jj] = - hq[ii][jj];
		for(jj=0; jj<nx0; jj++) 
			hrq[ii][nu0+jj] = - hq[ii][nu0+jj] + hpi[ii-1][jj];
		for(jj=0; jj<nb0; jj++) 
			hrq[ii][idxb[ii][jj]] += hlam[ii][jj] - hlam[ii][pnb+jj];
		dsymv_lib(nu0+nx0, nu0+nx0, hpQ[ii], cnux0, hux[ii], -1, hrq[ii], hrq[ii]);
		if(ng0>0)
			{
			// TODO work space + one dgemv call
			dgemv_n_lib(nu0+nx0, ng0, hpDCt[ii], cng, hlam[ii]+2*pnb, 1, hrq[ii], hrq[ii]);
			dgemv_n_lib(nu0+nx0, ng0, hpDCt[ii], cng, hlam[ii]+2*pnb+png, -1, hrq[ii], hrq[ii]);
			}

		for(jj=0; jj<nx1; jj++) 
			hrb[ii][jj] = hux[ii+1][nu1+jj] - hb[ii][jj];
		dgemv_nt_lib(nu0+nx0, nx1, hpBAbt[ii], cnx1, hpi[ii], hux[ii], -1, -1, hrq[ii], hrb[ii], hrq[ii], hrb[ii]);

		}
	


	// last block
	ii = N;
	nu0 = nu1;
	nx0 = nx1;
	cnux0  = (nu0+nx0+ncl-1)/ncl*ncl;
	nb0 = nb[ii];
	pnb = (nb0+bs-1)/bs*bs;
	ng0 = ng[ii];
	png = (ng0+bs-1)/bs*bs;
	cng = (ng0+ncl-1)/ncl*ncl;
	nb_tot += nb0 + ng0;

	for(jj=0; jj<nb0; jj++)
		mu[0] += hlam[ii][jj] * ht[ii][jj] + hlam[ii][pnb+jj] * ht[ii][pnb+jj];
	for(jj=0; jj<ng0; jj++) 
		mu[0] += hlam[ii][2*pnb+jj] * ht[ii][2*pnb+jj] + hlam[ii][2*pnb+png+jj] * ht[ii][2*pnb+png+jj];

	for(jj=0; jj<nb0; jj++)
		{
		hrd[ii][jj]     =   hux[ii][idxb[ii][jj]] - hd[ii][jj]     - ht[ii][jj];
		hrd[ii][pnb+jj] = - hux[ii][idxb[ii][jj]] + hd[ii][pnb+jj] - ht[ii][pnb+jj];
		}
	if(ng0>0)
		{
		dgemv_t_lib(nu0+nx0, ng0, hpDCt[ii], cng, hux[ii], 0, hrd[ii]+2*pnb, hrd[ii]+2*pnb);
		for(jj=0; jj<ng0; jj++)
			{
			hrd[ii][2*pnb+png+jj] = - hrd[ii][2*pnb+jj];
			hrd[ii][2*pnb+jj] += - hd[ii][2*pnb+jj] - ht[ii][2*pnb+jj];
			hrd[ii][2*pnb+png+jj] += hd[ii][2*pnb+png+jj] - ht[ii][2*pnb+png+jj];
			}
		}

	for(jj=0; jj<nx0; jj++) 
		hrq[ii][nu0+jj] = hpi[ii-1][jj] - hq[ii][nu0+jj];
	for(jj=0; jj<nb0; jj++) 
		hrq[ii][idxb[ii][jj]] += hlam[ii][jj] - hlam[ii][pnb+jj];
	dsymv_lib(nx0+nu0%bs, nx0+nu0%bs, hpQ[ii]+nu0/bs*bs*cnux0+nu0/bs*bs*bs, cnux0, hux[ii]+nu0/bs*bs, -1, hrq[ii]+nu0/bs*bs, hrq[ii]+nu0/bs*bs);
	if(ng0>0)
		{
		// TODO work space + one dgemv call
		dgemv_n_lib(nu0+nx0, ng0, hpDCt[ii], cng, hlam[ii]+2*pnb, 1, hrq[ii], hrq[ii]);
		dgemv_n_lib(nu0+nx0, ng0, hpDCt[ii], cng, hlam[ii]+2*pnb+png, -1, hrq[ii], hrq[ii]);
		}
	


	// normalize mu
	if(nb_tot!=0)
		mu[0] /= 2.0*nb_tot;



#if 1
	// change sign of residuals
	for(ii=0; ii<=N; ii++)
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			hrq[ii][jj] = - hrq[ii][jj];
	for(ii=0; ii<N; ii++)
		for(jj=0; jj<nx[ii+1]; jj++)
			hrb[ii][jj] = - hrb[ii][jj];
	for(ii=0; ii<=N; ii++)
		{
		pnb = (nb[ii]+bs-1)/bs*bs;
		png = (ng[ii]+bs-1)/bs*bs;
		for(jj=0; jj<nb[ii]; jj++)
			{
			hrd[ii][jj]     = - hrd[ii][jj];
			hrd[ii][pnb+jj] = - hrd[ii][pnb+jj];
			}
		for(jj=0; jj<ng[ii]; jj++)
			{
			hrd[ii][2*pnb+jj]     = - hrd[ii][2*pnb+jj];
			hrd[ii][2*pnb+png+jj] = - hrd[ii][2*pnb+png+jj];
			}
		}
#endif



	return;

	}


