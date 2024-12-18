%% Continuation of folds of periodic orbits
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
% </html>
%
% The extension ddebiftool_extra_psol is able to continue local
% bifurcations of periodic orbits in two parameters. This demo shows how
% one can continue folds (saddle-nodes) of periodic orbits i nthe neuron
% problem. (requires running of <demo1_psol.html> first).

%% Add extension folder to path
% The extension is installed in a separate folder. In addition, a folder
% with utilities for user convenience is loaded.
addpath('../ddebiftool_extra_psol/');
addpath('../ddebiftool_utilities/');
%#ok<*ASGLU,*NOPTS,*NASGU>
%% Speed up computations by vectorization
% The functions |neuron_sys_rhs| and |neuron_sys_deri| are
% not vectorized. In order to speed up computations we re-define
% |neuron_sys_rhs|, and replace |neuron_sys_seri| with the default
% finite-difference approximation
% neuron_sys_rhs=@(xx,par)[...
%     -par(1)*xx(1,1,:)+par(2)*tanh(xx(1,4,:))+par(3)*tanh(xx(2,3,:));....
%     -par(1)*xx(2,1,:)+par(2)*tanh(xx(2,4,:))+par(4)*tanh(xx(1,2,:))];
neuron_sys_rhs=@(xx,par)[...
    -xx(1,1,:)+1/(1+exp(-par(2)*(par(7)+par(3)*xx(1,2,:)+par(4)*xx(2,2,:))));....
    par(1)*(-xx(2,1,:)+1/(1+exp(-par(2)*(par(8)+par(5)*xx(1,2,:)+par(6)*xx(2,2,:)))))];
vfuncs=set_funcs(...
    'sys_rhs',neuron_sys_rhs,...
    'sys_tau',@()9,...
    'x_vectorized',true);
%% Find initial guess
% For the branch computed in <demo1_psol.html> one fold occured at the
% maximal parameter value. So we extract its index |indmax|. 
[~,indmax1]=max(arrayfun(@(x)x.parameter(ind_theta_u),branch5.point));
[~,indmax2]=min(arrayfun(@(x)x.parameter(ind_theta_u),branch5.point));

%% Initialize branch and set up extended system
% Then we call |SetupPOfold| to initialize the branch and to create the
% functions for the extended system. For the core DDE-Biftool routines the
% fold of periodic orbits is of the same type as a standard periodic orbit.
% However, the user-provided right-hand side has been extended (eg,
% foldfuncs.sys_rhs is different from funcs.sys_rhs). SetupPOfold has three
% mandatory parameters:
%
% * |funcs| containing the user-problem functions (in the demo: |vfuncs|),
% * |branch|, the branch on which the fold occurred (in the demo:
% |branch5|), and
% * |ind|, the index in |branch.point| near which the fold occurs (in the
% demo |indmax|).
%
% Output (to be fed into |br_contn|):
% 
% * |foldfuncs|: functions for the extended system (derived from |funcs|)
% * |FPObranch|: branch of folds of periodic orbits with two points
%
% The other parameters are optional name-value pairs. Important are
%
% * |'contpar'|: typically two integers ,the two continuation parameters,
% * |'dir'|: in which direction the first step is taken), and 
% * |'step'|: the length of the first step. 
% Note that name-value pairs are also passed on to fields of the |branch|
% structure: the output |FPObranch| inherits all fields from the input
% |branch|. Additional name-value argument pairs canbe used to change
% selected fields.
branch5.method.point.newton_max_iterations=16;
[foldfuncs,branch7]=SetupPOfold(vfuncs,branch5,indmax1,'contpar',[ind_theta_u,ind_taus],...
    'dir',ind_taus,'print_residual_info',1,'step',0.01,'plot_measure',[],...
    'min_bound',[ind_theta_u,0.572535; ind_taus,0.102212], ...
    'max_bound',[ind_theta_u,0.827663; ind_taus,0.5],...
    'max_step',[ind_theta_u,0.01; ind_taus,0.01]);
%% Branch continuation
% The output of |SetupPOfold| can be fed directly into |br_contn| to
% perform a branch continuation. Note that the computation is typically
% slower than a standard continuation of periodic orbits because the system
% dimension is twice as large and additional delays have been introduced.
figure(13); clf;
branch7.method.point.print_residual_info=0;
branch7=br_contn(foldfuncs,branch7,100);
branch7=br_rvers(branch7);
branch7=br_contn(foldfuncs,branch7,100);
xlabel('\theta_{u}');ylabel('\tau');
title('Continuation of fold of periodic orbits');
%% Initilaise and continue second periodic orbit branch
[foldfuncs,branch8]=SetupPOfold(vfuncs,branch5,indmax2,'contpar',[ind_theta_u,ind_taus],...
    'dir',ind_taus,'print_residual_info',1,'step',0.01,'plot_measure',[],...
    'min_bound',[ind_theta_u,0.4; ind_taus,0.0], ...
    'max_bound',[ind_theta_u,1.0; ind_taus,0.5],...
    'max_step',[ind_theta_u,0.01; ind_taus,0.01]);
branch8.method.point.print_residual_info=0;
branch8=br_contn(foldfuncs,branch8,100);
branch8=br_rvers(branch8);
branch8=br_contn(foldfuncs,branch8,100);
%% Extracting solution components
% Since the fold continuation solves an extended system (with additional
% components) it does not make sense to compute the linear stability of the
% resulting orbits directly. However, |foldfuncs| has an additional field
% |'get_comp'|. This is a function that extacts the components of the
% extended system that correspond to the solution of the original system.
% The output of |get_comp| is an array of periodic orbits that are located
% precisely on the fold. 
pf_orbits1=foldfuncs.get_comp(branch7.point,'solution');
pf_orbits2=foldfuncs.get_comp(branch8.point,'solution');

%% Stability of fold orbits
% We compute Floquet multipliers for these orbits using the utility
% function |GetStability(psol_array,...)|. Its optional inputs instruct it
% to ignore the two eigenvalues closest to 1 when computing stability. It
% also returns |triv_defect| which measures the distance of the supposedly
% trivial eigenvalues (as computed) to their correct values. |GetStability|
% needs the optional input |'funcs'| if its first argument does not yet
% have a stability structure.
[~,~,~,pf_orbits1]=GetStability(pf_orbits1,'funcs',vfuncs);
[~,~,~,pf_orbits2]=GetStability(pf_orbits2,'funcs',vfuncs);

theta_u_pfold1=arrayfun(@(x)x.parameter(7),branch7.point);
taus_pfold1=arrayfun(@(x)x.parameter(9),branch7.point);
stability_pfold1=arrayfun(@(x)real(x.stability.mu(1)),pf_orbits1);

theta_u_pfold2=arrayfun(@(x)x.parameter(7),branch8.point);
taus_pfold2=arrayfun(@(x)x.parameter(9),branch8.point);
stability_pfold2=arrayfun(@(x)real(x.stability.mu(1)),pf_orbits2);

figure(14);
subplot(2,1,1)
plot(theta_u_pfold1,stability_pfold1)
xlabel('\theta_{u}');
ylabel('Stability');
subplot(2,1,2)
plot(theta_u_pfold2,stability_pfold2)
xlabel('\theta_{u}');
ylabel('Stability');
%% Save results (end of tutorial demo, but try also <demo1_hcli.html>)
save('demo1_POfold_results.mat')
