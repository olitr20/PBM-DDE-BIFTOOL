%% Continuation and stability of steady states (equilibria)
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(47), 19/04/2014
% </html>
%
% Once the user-defined functions are prepared, DDE-Biftool can compute and
% continue equilibria of the DDE, and compute their linearized stability,
% thus detecting local bifurcations. This demo requires <demo1_funcs.html> to
% have run beforehand.
%%
%#ok<*ASGLU,*NOPTS,*NASGU>
%
%% Initial guess for steady state
% It is clear that the neuron DDE has a steady state solution
% $(x_1^*,x_2^*)=(0,0)$ for all values of the parameters. We define a first
% steady state solution using the parameter values $\kappa=0.5$,
% $\beta=-1$, $a_{12}=1$, $a_{21}=2.34$, $\tau_1=\tau_2=0.2$ and
% $\tau_s=1.5$. Remember that we chose |par| $=[\kappa, \beta, a_{12},
% a_{21},\tau_1,\tau_2, \tau_s]$.
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4; p.c = -1; p.d = 0;
p.theta_u = 0.7; p.theta_v = 0.5;
p.tau_1 = 0.2; p.tau_2 = p.tau_1;

stst.kind='stst';
% stst.parameter=[p.alpha, p.beta, p.a, p.b, p.c, p.d, p.theta_u, p.theta_v, p.tau_1, p.tau_2];
stst.parameter=[p.alpha, p.beta, p.a, p.b, p.c, p.d, p.theta_u, p.theta_v, p.tau_1];
% stst.x=[0.103596478880629;1] % theta_u = 0.5
% stst.x=[0.599328884047436;0] % theta_u = 0.6
stst.x=[0.698598941633828;0] % theta_u = 0.7
% stst.x=[0.797713228634354;0] % theta_u = 0.8
% stst.x=[0.896403521119372;0] % theta_u = 0.9

%% Linear stability of initial equilibrium
% We get default point method parameters and correct the point, which,
% being already a correct solution, remains unchanged. Computing and
% plotting stability of the corrected point reveals it has one unstable
% real mode, see figure.

flag_newhheur=1; % flag_newhheur=1 is the default choice if this argument is omitted
method=df_mthod(funcs,'stst',flag_newhheur);
method.stability.minimal_real_part=-1
[stst,success]=p_correc(funcs,stst,[],[],method.point)
% compute its stability:
stst.stability=p_stabil(funcs,stst,method.stability)
figure(1); clf;
p_splot(stst); % plot its stability:
%% Ask for roots with more negative real part
% In both figures, approximations $(\times)$ and corrections $(*)$ are
% nearly indistinguishable.
method.stability.minimal_real_part=-2; 
stst.stability=p_stabil(funcs,stst,method.stability); % recompute stability:
figure(2); clf;
p_splot(stst); % replot stability
%% Figures: Spectrum of equilibrium
% Approximated $(\times)$ and corrected $(*)$ roots of the characteristic
% equation of neuron system at its steady state solution
% $(x_1^*,x_2^*)=(0,0)$. Real parts computed up to $\Re(\lambda)\geq
% -\frac{1}{\tau}$ (top), $\Re(\lambda)\geq -2$ (bottom).

%% Initialize branch of trivial equilibria
% We will use this point as a first point to compute a branch of steady
% state solutions. First, we obtain an empty branch with free parameter
% $a_{21}$, limited by $a_{21}\in[0,5]$ and $\Delta a_{21}\leq 0.2$ between
% points.

% get an empty branch with ind_taus as a free parameter:
branch1=df_brnch(funcs,ind_theta_u,'stst');
% set bounds for continuation parameter
branch1.parameter.min_bound(1,:)=[ind_theta_u 0.4];
branch1.parameter.max_bound(1,:)=[ind_theta_u 1];
branch1.parameter.max_step(1,:)=[ind_theta_u 0.005];
% use stst as a first branch point:
branch1.point=stst;

branch1.parameter
branch1.parameter.min_bound

%%  Extend and continue branch of trivial equilibria
% To obtain a second starting point we change  parameter value $a_{21}$
% slightly and correct again.

figure(3); clf;

stst.parameter(ind_theta_u)=stst.parameter(ind_theta_u)+0.005;
[stst,success]=p_correc(funcs,stst,[],[],method.point)
% use as a second branch point:
branch1.point(2)=stst;
branch1.method.continuation.plot=1;

% continue in one direction:
[branch1,s,f,r]=br_contn(funcs,branch1,200)
% turn the branch around:
branch1=br_rvers(branch1);
% continue in the other direction:
[branch1,s,f,r]=br_contn(funcs,branch1,200)
xlabel('\theta_{u}');ylabel('u');

%% Stability of branch of equilibria
% During continuation, sixteen points were successfully computed ($s=16$)
% before the right boundary $a_{21}=5$ was hit (signalled by a warning). No
% corrections failed ($f=0$) and no computed points were later rejected
% ($r=0$). Reversing the order of the branch points allows to continue to
% the left.
%
% After obtaining suitable measure structures we plot the real part of the
% approximated and corrected roots of the characteristic equation along the
% branch, (see figure). Notice the strange behaviour (coinciding of several
% complex pairs of roots) at $a_{21}=0$. At this parameter value one of the
% couplings between the neurons is broken. In fact, for $a_{21}=0$, the
% evolution of the second component is independent of the evolution of the
% first.
branch1.method.stability.minimal_real_part=-2;
branch1=br_stabl(funcs,branch1,0,0);

% obtain suitable scalar measures to plot stability along branch:
[xm,ym]=df_measr(1,branch1)

figure(4); clf;
ym.subfield='l0';
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
plot([0 250],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');
%% Figures: Stability of equilibria
%
% <html><a name=stststability></a>
% </html>
%
% Real parts of the approximated (top) and corrected (top,bottom) roots of
% the characteristic equation versus $a_{21}$ (top) respectively the point
% number along the branch (bottom).
%% Save and continue with Hopf bifurcations: <demo1_hopf.html>
save('demo1_stst_results.mat');
