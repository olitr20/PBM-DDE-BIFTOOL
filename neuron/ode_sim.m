%% Simulation
% Select Parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4;
p.c = -1; p.d = 0;

% Find u* and v* such that θu = 0.4 & θv = 0.5
p.theta_u = 0.7;
p.theta_v = 0.5;

biases = @(vars) [
        f_inv(vars(1),p) - (p.a .* vars(1)) - (p.b .* vars(2)) - p.theta_u;
        f_inv(vars(2),p) - (p.c .* vars(1)) - (p.d .* vars(2)) - p.theta_v;
    ];

sol = fsolve(biases, [0.5, 0.5], optimoptions('fsolve', 'Display', 'off'));

p.u = sol(1);
p.v = sol(2);

clear sol

% Define simulation parameters
xmax = 30; xstep = 0.1; % x range
tmax = 1500; % max simulation time
smpl = 5; % initial value sample rate
[X,Y] = meshgrid(0:smpl:xmax);
% x0 = [X(:), Y(:)]; clear X Y
x0 = [1; 0];
tspan = [0 tmax];

[t,xode] = ode23s(@(t, x) odefun(t, x, p), tspan, x0);

fprintf('u = %.15f\n', xode(end,1));
fprintf('v = %.15f\n', xode(end,2));
%% --------------------------------------------------------------------- %%
% ------------------------------- f(x,p) -------------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function f = f(z)
        % f = heaviside(z);
        f = 1 ./ (1 + exp(-600 * z));
    end

% ----------------------------- f_inv(x,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ---------------------------- odefun(t,y,Z) ---------------------------- %
    function d = odefun(~,x,p)
        u = x(1);
        v = x(2);

        dudt = -u + f(p.theta_u + p.a .* u + p.b .* v);
        dvdt = p.alpha .* (-v + f(p.theta_v + p.c .* u + p.d .* v));
    
        d = [dudt; dvdt];
    end
