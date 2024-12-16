%% Hopf bifurcation
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(73), 31/12/2014
% </html>
%
% The eigenvalues of the linearized system along branches of equilibria
% indicate potential bifurcations. In this demo complex conjugate pairs of
% eigenvalues cross the imaginary axis, corresponding to Hopf bifurcations.
% The demo will proceed to continue two of these Hopf bifurcations in two
% system parameters $a_{21}$ and $\tau_s$. This part requires to run
% <demo1_stst.html> first
%%
%#ok<*ASGLU,*NOPTS,*NASGU>
%
%% Locating the first Hopf point
% Where eigenvalue curves in the stability plot (see <demo1_stst.html#stststability>)
% cross the zero line, bifurcations occur. If we want to compute the Hopf
% bifurcation near $a_{21}\approx0.8$ we need its point number. This is
% most easily obtained by plotting the stability versus the point numbers
% along the branch. We select the last point with positive eigenvalues and
% turn it into an (approximate) Hopf bifurcation point. We correct the Hopf
% point using appropriate method parameters and one free parameter
% ($a_{21}$). We then copy the corrected point to keep it for later use.
% Computing and plotting stability of the Hopf point clearly reveals the
% pair of purely imaginary eigenvalues.
ind_hopf=find(arrayfun(@(x)real(x.stability.l0(1))>0,branch1.point),1,'last');
hopf=p_tohopf(funcs,branch1.point(ind_hopf));
method=df_mthod(funcs,'hopf',flag_newhheur); % get hopf calculation method parameters:
method.stability.minimal_real_part=-1;
[hopf,success]=p_correc(funcs,hopf,ind_theta_u,[],method.point) % correct hopf
first_hopf=hopf;                    % store hopf point in other variable for later use
hopf.stability=p_stabil(funcs,hopf,method.stability); % compute stability of hopf point
figure(5); clf;
p_splot(hopf);                     % plot stability of hopf point
%% Figure: eigenvalues at Hopf point
% Characteristic roots at Hopf point: a pair of pure imaginary eigenvalues
% is clearly visible.

%% Initialize and continue first Hopf bifurcation
% In order to follow a branch of Hopf bifurcations in the two parameter
% space $(a_{21},\tau_s)$ we again need two starting points. Hence we use
% the Hopf point already found and one perturbed in $\tau_s$ and corrected
% in $a_{21}$, to start on a branch of Hopf bifurcations. For the free
% parameters, $a_{21}$ and $\tau_s$, we provide suitable intervals,
% $a_{21}\in[0,4]$ and $\tau_s\in[0,10]$, and maximal stepsizes, $0.2$ for
% $a_{21}$ and $0.5$ for $\tau_s$. We continue the branch on both sides by
% an intermediate order reversal and a second call to |br_contn|.
branch2=df_brnch(funcs,[ind_theta_u,ind_taus],'hopf'); % use hopf point as first point of hopf branch:
branch2.parameter.min_bound(1,:)=[ind_theta_u 0.4];
branch2.parameter.max_bound(1:2,:)=[[ind_theta_u 1]' [ind_taus 0.5]']';
branch2.parameter.max_step(1:2,:)=[[ind_theta_u 0.005]' [ind_taus 0.005]']';
branch2.point=hopf;

hopf.parameter(ind_theta_u)=hopf.parameter(ind_theta_u)-0.0001; % perturb hopf point
[hopf,success]=p_correc(funcs,hopf,ind_taus,[],method.point); % correct hopf point, recompute stability
branch2.point(2)=hopf;                                 % use as second point of hopf branch:
figure(6); clf;
[branch2,s,f,r]=br_contn(funcs,branch2,1000);            % continue with plotting hopf branch:
branch2=br_rvers(branch2);                             % reverse Hopf branch
[branch2,s,f,r]=br_contn(funcs,branch2,1000);            % continue in other direction
xlabel('\theta_{u}');ylabel('\tau');

%% Save and continue 
% Continue with with periodic orbits <demo1_psol.html> or normal forms
% <demo1_normalforms.html>.
save('demo1_hopf_results.mat');
