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
* Author: Mikhail Katliar, mikhai.katliar (at) tuebingen.mpg.de                                   *
*                                                                                                 *
**************************************************************************************************/

#include "mex_helper/MexHelper.hpp"

#include "/opt/hpmpc/include/c_interface.h"

#include <stdexcept>
#include <array>
#include <cmath>
#include <algorithm>
#include <vector>

namespace hpmpc_wrapper
{
    static int fortran_order_d_ip_ocp_hard_tv(
		int *kk, int k_max, double mu0, double mu_tol,
		int N, int const *nx, int const *nu, int const *nb, int const * const *hidxb, int const *ng, int N2,
		int warm_start,
		double const * const *A, double const * const *B, double const * const *b,
		double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
		double const * const *lb, double const * const *ub,
		double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
		double * const *x, double * const *u, double * const *pi, double * const *lam,
		double *inf_norm_res,
		void *work0,
		double *stat)
	{
		return ::fortran_order_d_ip_ocp_hard_tv(
			kk, k_max, mu0, mu_tol,
			N, const_cast<int*>(nx), const_cast<int*>(nu), const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2,
			warm_start,
			const_cast<double**>(A), const_cast<double**>(B), const_cast<double**>(b),
			const_cast<double**>(Q), const_cast<double**>(S), const_cast<double**>(R), const_cast<double**>(q), const_cast<double**>(r),
			const_cast<double**>(lb), const_cast<double**>(ub),
			const_cast<double**>(C), const_cast<double**>(D), const_cast<double**>(lg), const_cast<double**>(ug),
			const_cast<double**>(x), const_cast<double**>(u), const_cast<double**>(pi), const_cast<double**>(lam), 
			inf_norm_res,
			work0,
			stat);
	}

	int hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2)
	{
		return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
			N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2);
	}
}

std::invalid_argument errorInconsistentMatrixSize(std::string const& name, int stage)
{
    return std::invalid_argument("Inconsistent size of " + name + " at stage " + std::to_string(stage));
}

/*
class Bounds
{
public:
    Bounds(mxArray const * lbx, mxArray const * ubx, mxArray const * lbu, mxArray const * ubu)
    {
        // Check lbx
        mxArray const * const lbx = mx::GetField(problem, k, "lbx");
        auto const size_lbx(lbx);
        auto const nx = size_lbx[0];

        if (size_lbx != mx::ArraySize(nx, 1))
            throw errorInconsistentMatrixSize("lbx");

        // Check ubx
        mxArray const * const ubx = mx::GetField(problem, k, "ubx");
        if (mx::ArraySize(ubx) != mx::ArraySize(nx, 1))
            throw errorInconsistentMatrixSize("ubx");
        
        // Check lbu
        mxArray const * const lbu = mx::GetField(problem, k, "lbu");
        auto const size_lbu(lbu);
        auto const nu = size_lbu[0];

        if (size_lbu != mx::ArraySize(nu, 1))
            throw errorInconsistentMatrixSize("lbu");

        // Check ubu
        mxArray const * const ubu = mx::GetField(problem, k, "ubu");
        if (mx::ArraySize(ubu) != mx::ArraySize(nu, 1))
            throw errorInconsistentMatrixSize("ubu");
        
        int idx_bound = 0;
        double const * lbx_data = mx::GetPr(lbx);
        double const * ubx_data = mx::GetPr(ubx);
        double const * lbu_data = mx::GetPr(lbu);
        double const * ubu_data = mx::GetPr(ubu);
        
        for (int i = 0; i < nu; ++i)
        {
            if (std::isfinite(lbu_data[i]) && std::isfinite(ubu_data[i]))
                
                
            if (std::isinf(lbu_data[i]) && lbu_data[i] < 0)
        }
    }
    
private:
    std::vector<double> lb_;
    std::vector<double> ub_;
    std::vector<int> idxb_;
};
*/

class ProblemStage
{
public:
    ProblemStage(mxArray const * problem, mwIndex k)
    {
        // Check Q
        mxArray const * const Q = mx::GetField(problem, k, "Q");
        nx_ = mx::GetM(Q);

        if (!mx::GetSize(Q).isEqual({nx_, nx_}))
            throw errorInconsistentMatrixSize("Q", k);

        // Check R
        mxArray const * const R = mx::GetField(problem, k, "R");
        nu_ = mx::GetM(R);

        if (!mx::GetSize(R).isEqual({nu_, nu_}))
            throw errorInconsistentMatrixSize("R", k);

        // Check S
        mxArray const * const S = mx::GetField(problem, k, "S");
        if (!mx::GetSize(S).isEqual({nu_, nx_}))
            throw errorInconsistentMatrixSize("S", k);

        // Check q
        mxArray const * const q = mx::GetField(problem, k, "q");
        if (!mx::GetSize(q).isEqual({nx_, 1}))
            throw errorInconsistentMatrixSize("q", k);

        // Check r
        mxArray const * const r = mx::GetField(problem, k, "r");
        if (!mx::GetSize(r).isEqual({nu_, 1}))
            throw errorInconsistentMatrixSize("r", k);

        // Check A
        mxArray const * const A = mx::GetField(problem, k, "A");
        nxNext_ = mx::GetM(A);

        if (!mx::GetSize(A).isEqual({nxNext_, nx_}))
            throw errorInconsistentMatrixSize("A", k);

        // Check B
        mxArray const * const B = mx::GetField(problem, k, "B");
        if (!mx::GetSize(B).isEqual({nxNext_, nu_}))
            throw errorInconsistentMatrixSize("B", k);

        // Check b
        mxArray const * const b = mx::GetField(problem, k, "b");
        if (!mx::GetSize(b).isEqual({nxNext_, 1}))
            throw errorInconsistentMatrixSize("b", k);

        // Check C
        mxArray const * const C = mx::GetField(problem, k, "C");
        ng_ = mx::GetM(C);

        if (!mx::GetSize(C).isEqual({ng_, nx_}))
            throw errorInconsistentMatrixSize("C", k);

        // Check D
        mxArray const * const D = mx::GetField(problem, k, "D");
        if (!mx::GetSize(D).isEqual({ng_, nu_}))
            throw errorInconsistentMatrixSize("D", k);

        // Check lg
        mxArray const * const lg = mx::GetField(problem, k, "lg");
        if (!mx::GetSize(lg).isEqual({ng_, 1}))
            throw errorInconsistentMatrixSize("lg", k);

        // Check ug
        mxArray const * const ug = mx::GetField(problem, k, "ug");
        if (!mx::GetSize(ug).isEqual({ng_, 1}))
            throw errorInconsistentMatrixSize("ug", k);

        // Check lb
        mxArray const * const lb = mx::GetField(problem, k, "lb");
        nb_ = mx::GetM(lb);

        if (!mx::GetSize(lb).isEqual({nb_, 1}))
            throw errorInconsistentMatrixSize("lb", k);

        // Check ub
        mxArray const * const ub = mx::GetField(problem, k, "ub");
        if (!mx::GetSize(ub).isEqual({nb_, 1}))
            throw errorInconsistentMatrixSize("ub", k);

        // Check idxb
        mxArray const * const idxb = mx::GetField(problem, k, "idxb");
        if (!mx::GetSize(idxb).isEqual({1, nb_}))
            throw errorInconsistentMatrixSize("idxb", k);

        if (!mxIsInt32(idxb))
            throw std::invalid_argument("idxb must be an array of int32");

        A_ = mx::GetPr(A);
        B_ = mx::GetPr(B);
        b_ = mx::GetPr(b);
        Q_ = mx::GetPr(Q);
        R_ = mx::GetPr(R);
        S_ = mx::GetPr(S);
        q_ = mx::GetPr(q);
        r_ = mx::GetPr(r);
        C_ = mx::GetPr(C);
        D_ = mx::GetPr(D);
        lg_ = mx::GetPr(lg);
        ug_ = mx::GetPr(ug);

        lb_ = mx::GetPr(lb);
        ub_ = mx::GetPr(ub);
        idxb_ = mx::GetData<std::int32_t>(idxb);
    }

    int nx() const { return nx_; }
    int nu() const { return nu_; }
    int nb() const { return nb_; }
    int ng() const { return ng_; }
    int const * idxb() const { return idxb_; }
    int nxNext() const { return nxNext_; }
    double const * Q() const { return Q_; }
    double const * R() const { return R_; }
    double const * S() const { return S_; }
    double const * q() const { return q_; }
    double const * r() const { return r_; }
    double const * A() const { return A_; }
    double const * B() const { return B_; }
    double const * b() const { return b_; }
    double const * lb() const { return lb_; }
    double const * ub() const { return ub_; }
    double const * C() const { return C_; }
    double const * D() const { return D_; }
    double const * lg() const { return lg_; }
    double const * ug() const { return ug_; }

private:
    int nx_ = 0;
    int nu_ = 0;
    int nb_ = 0;
    int ng_ = 0;
    int const * idxb_ = 0;
    int nxNext_ = 0;

    double const * Q_ = nullptr;
    double const * R_ = nullptr;
    double const * S_ = nullptr;
    double const * q_ = nullptr;
    double const * r_ = nullptr;
    double const * A_ = nullptr;
    double const * B_ = nullptr;
    double const * b_ = nullptr;
    double const * lb_ = nullptr;
    double const * ub_ = nullptr;
    double const * C_ = nullptr;
    double const * D_ = nullptr;
    double const * lg_ = nullptr;
    double const * ug_ = nullptr;
};


class Problem
:   public std::vector<ProblemStage>
{
public:
    Problem(mxArray const * problem)
    {
        if (!mx::IsStruct(problem))
            throw std::invalid_argument("The 'problem' argument must be a structure");

        if (!mx::IsVector(problem))
            throw std::invalid_argument("The 'problem' argument must be a vector");

        auto const N = mx::GetNumberOfElements(problem);

        if (N < 2)
            throw std::invalid_argument("The 'problem' argument must be at least of size 2");
        
        reserve(N);

        for (int k = 0; k < N; ++k)
            emplace_back(problem, k);

        for (int k = 0; k + 1 < N; ++k)
        {
            if (at(k).nxNext() != at(k + 1).nx())
                throw std::invalid_argument("Inconsistent A, B matrix sizes between stages " 
                    + std::to_string(k) + " and " + std::to_string(k + 1));
        }
    }
};


class Options
{
public:
    Options(mxArray const * options)
    {
        if (!mxIsStruct(options))
        {
            throw std::invalid_argument("The 'options' argument must be a structure");
        }

        if (!mxIsScalar(options))
        {
            throw std::invalid_argument("The 'options' argument must be a scalar");
        }
        
        kMax_ = static_cast<int>(mx::GetScalar<double>(mx::GetField(options, 0, "k_max")));
        mu0_ = mx::GetScalar<double>(mx::GetField(options, 0, "mu0"));
        tol_ = mx::GetScalar<double>(mx::GetField(options, 0, "tol"));
        
        auto const warm_start = mx::GetField(options, 0, "warm_start");
        warmStart_ = warm_start && mxIsLogicalScalar(warm_start) && *mxGetLogicals(warm_start);
    }

    int k_max() const
    {
        return kMax_;
    }

    double mu0() const
    {
        return mu0_;
    }

    double tol() const
    {
        return tol_;
    }

    bool warmStart() const
    {
        return warmStart_;
    }
    
private:
    int kMax_;
    double mu0_;
    double tol_;
    bool warmStart_;
};


class SolutionStage
{
public:
    SolutionStage(int nx, int nu, int nb, int ng, int nx_next)
    :   x_(mx::CreateDoubleMatrix(nx, 1))
    ,   u_(mx::CreateDoubleMatrix(nu, 1))
    ,   pi_(mx::CreateDoubleMatrix(nx_next, 1))
    ,   lam_(mx::CreateDoubleMatrix(2 * nb + 2 * ng, 1))
    {        
    }
    
    mxArray * x()
    {
        return x_;
    }
    
    mxArray * u()
    {
        return u_;
    }
    
    mxArray * pi()
    {
        return pi_;
    }
    
    mxArray * lam()
    {
        return lam_;
    }
    
    double * x_data()
    {
        return mx::GetPr(x_);
    }
    
    double * u_data()
    {
        return mx::GetPr(u_);
    }
    
    double * pi_data()
    {
        mx::GetPr(pi_);
    }
    
    double * lam_data()
    {
        return mx::GetPr(lam_);
    }
    
private:
    mxArray * x_;
    mxArray * u_;
    mxArray * pi_;
    mxArray * lam_;
};


class Result
{
public:
    Result(std::vector<ProblemStage> const& problem, Options const& options)
    {
        // Number of QP stages, including the terminal one.
        auto const N = problem.size();
        
        // HPMPC input variables
        std::vector<double const *> A;
        std::vector<double const *> B;
        std::vector<double const *> b;
        std::vector<double const *> Q;
        std::vector<double const *> S;
        std::vector<double const *> R;
        std::vector<double const *> q;
        std::vector<double const *> r;
        std::vector<double const *> lb;
        std::vector<double const *> ub;
        std::vector<double const *> C;
        std::vector<double const *> D;
        std::vector<double const *> lg;
        std::vector<double const *> ug;
        std::vector<int> nx;
        std::vector<int> nu;
        std::vector<int> nb;
        std::vector<int> ng;
        std::vector<int const *> idxb;
        
        nx.reserve(N);
        nu.reserve(N);
        nb.reserve(N);
        ng.reserve(N);
        idxb.reserve(N);

        A .reserve(N);
        B .reserve(N);
        b .reserve(N);
        Q .reserve(N);
        S .reserve(N);
        R .reserve(N);
        q .reserve(N);
        r .reserve(N);
        lb.reserve(N);
        ub.reserve(N);
        C .reserve(N);
        D .reserve(N);
        lg.reserve(N);
        ug.reserve(N);

        // Fill HPMPC input arrays
        for (auto const& stage : problem)
        {
            nx.push_back(stage.nx());
            nu.push_back(stage.nu());
            ng.push_back(stage.ng());
            nb.push_back(stage.nb());

            A.push_back(stage.A());
            B.push_back(stage.B());
            b.push_back(stage.b());
            Q.push_back(stage.Q());
            R.push_back(stage.R());
            S.push_back(stage.S());
            q.push_back(stage.q());
            r.push_back(stage.r());
            C.push_back(stage.C());
            D.push_back(stage.D());
            lg.push_back(stage.lg());
            ug.push_back(stage.ug());

            lb.push_back(stage.lb());
            ub.push_back(stage.ub());
            idxb.push_back(stage.idxb());
        }
        
        // Prepare output variables
        std::vector<SolutionStage> solution;
        solution.reserve(N);
        
        for (auto const& stage : problem)
            solution.emplace_back(stage.nx(), stage.nu(), stage.nb(), stage.ng(), stage.nxNext());
        
        // HPMPC output variables
        std::vector<double *> x;
        std::vector<double *> u;
        std::vector<double *> pi;
        std::vector<double *> lam;
        
        x.reserve(N);
        u.reserve(N);
        pi.reserve(N);
        lam.reserve(N);
        
        for (auto& stage : solution)
        {
            x.push_back(stage.x_data());
            u.push_back(stage.u_data());
            pi.push_back(stage.pi_data());
            lam.push_back(stage.lam_data());
        }
        
        infNormRes_ = mx::CreateDoubleMatrix(1, 4);
        stat_ = mx::CreateDoubleMatrix(5, options.k_max());
        
        // Allocate solver workspace
        int const N_hpmpc = N - 1;

        std::vector<char> solver_workspace(hpmpc_wrapper::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
            N_hpmpc, nx.data(), nu.data(), nb.data(), idxb.data(), ng.data(), N_hpmpc));
        
        // Call the solver
        int iter_count = 0;
        auto const ret = hpmpc_wrapper::fortran_order_d_ip_ocp_hard_tv(&iter_count, options.k_max(), options.mu0(), options.tol(), N_hpmpc,
            nx.data(), nu.data(), nb.data(), idxb.data(), ng.data(), N_hpmpc, options.warmStart() ? 1 : 0, A.data(), B.data(), b.data(),
            Q.data(), S.data(), R.data(), q.data(), r.data(), lb.data(), ub.data(), C.data(), D.data(),
            lg.data(), ug.data(), x.data(), u.data(), pi.data(), lam.data(), mx::GetPr(infNormRes_),
            solver_workspace.data(), mx::GetPr(stat_));

        // Set the number of cols in the stat matrix equal to number of actual iterations
        mx::SetN(stat_, iter_count);

        // Create the solution structure
        solution_ = mx::CreateStructMatrix(1, N, {"x", "u", "pi", "lam"});

        // Copy solution mx arrays to the output structure.
        for (mwIndex k = 0; k < N; ++k)
        {
            mx::SetField(solution_, k, "x", solution[k].x());
            mx::SetField(solution_, k, "u", solution[k].u());
            mx::SetField(solution_, k, "pi", solution[k].pi());
            mx::SetField(solution_, k, "lam", solution[k].lam());
        }

        // Store the return code
        returnCode_ = mx::CreateNumericScalar(ret);
    }

    mxArray * solution()
    {
        return solution_;
    }

    mxArray * returnCode()
    {
        return returnCode_;
    }

    mxArray * stat()
    {
        return stat_;
    }

    mxArray * infNormRes()
    {
        return infNormRes_;
    }

private:
    mxArray * solution_;
    mxArray * returnCode_;
    mxArray * stat_;
    mxArray * infNormRes_;
};


// The entry function
extern "C" void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{		
    // Check arguments
    if (nrhs != 2)
        throw std::invalid_argument("HPMPC_ip_mpc_hard_new() expects 2 input arguments");

    Result result { Problem { prhs[0] }, Options { prhs[1] } };
    
    if (nlhs > 0)
        plhs[0] = result.solution();

    if (nlhs > 1)
        plhs[1] = result.returnCode();

    if (nlhs > 2)
        plhs[2] = result.stat();

    if (nlhs > 3)
        plhs[3] = result.infNormRes();
}

