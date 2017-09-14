% Test problem

%% Set dimensions
clear;

nx = 3;			% number of states
nu = 2;				% number of inputs (controls)
N = 1;				% horizon length

nb  = nx + nu;		% number of two-sided box constraints
ng  = 0;            % number of two-sided general constraints
ngN = 0;           % number of two-sided general constraints on last stage

%% Stage 1

qp(1).A = zeros(0, nx);
qp(1).B = zeros(0, nu);
qp(1).b = zeros(0, 1);

qp(1).Q = [ ...
    2.73449304081811,	1.88589647487061,	2.07850896302612;
	1.88589647487061,	2.23398400982993,	2.04607087035366;
	2.07850896302612,	2.04607087035366,	2.75909914758841;
];
qp(1).R = [ ...
    2.58480403401811,	2.27683785836085;
	2.27683785836085,	2.48531656389865;
];
qp(1).S = [ ...
    1.94423942432824,	1.95671439063179;
	2.31644542710892,	2.08748866537729;
	2.46061457306194,	1.94728938270016;
].';
qp(1).q = [ ...
    0.276025076998578;
    0.679702676853675;
    0.655098003973841;
];
qp(1).r = [ ...
    0.162611735194631;
	0.118997681558377;
];

qp(1).C = zeros(ng, nx);
qp(1).D  = zeros(ng, nu);
qp(1).lg = zeros(ng, 1);
qp(1).ug = zeros(ng, 1);

qp(1).lb = repmat(-10000, nb, 1);
qp(1).ub = repmat(10000, nb, 1);
qp(1).idxb = int32(0 : nb - 1);

%% Stage 2

qp(2).A = [];
qp(2).B = [];
qp(2).b = zeros(0, 1);

qp(2).Q = [];
qp(2).R = [];
qp(2).S = [];
qp(2).q = zeros(0, 1);
qp(2).r = zeros(0, 1);

qp(2).C = zeros(ngN, 0);
qp(2).D = zeros(ngN, 0);
qp(2).lg = zeros(ngN, 1);
qp(2).ug = zeros(ngN, 1);

qp(2).lb = zeros(0, 1);
qp(2).ub = zeros(0, 1);
qp(2).idxb = int32(0 : -1);

%% Set options and run solver
options.warm_start = false; % read initial guess from x and u
options.mu0 = 2;        % max element in cost function as estimate of max multiplier
options.k_max = 20;		% maximim number of iterations
options.tol = 1e-8;		% tolerance in the duality measure

[sol, ret, stat, inf_norm_res] = HPMPC_ip_mpc_hard_new(qp, options);
