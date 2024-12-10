function J = neuron_sys_deri(xx,par,nx,np,v)
%% Partial derivatives of r.h.s wrt xx and par
%
% (c) DDE-BIFTOOL v. 3.1.1(69), 24/12/2014
%
%%

% zu = par(7)+par(3)*xx(1,2)+par(4)*xx(2,3);
% zv = par(8)+par(5)*xx(1,3)+par(6)*xx(2,2);
zu = par(7)+par(3)*xx(1,2)+par(4)*xx(2,2);
zv = par(8)+par(5)*xx(1,2)+par(6)*xx(2,2);
sigmau = 1/(1+exp(-par(2)*zu));
sigmav = 1/(1+exp(-par(2)*zv));

J = [];

if isscalar(nx) && isempty(np) && isempty(v)
    switch nx
    case 0 % wrt u(t) and v(t)
        J = [-1, 0;
             0, -par(1)];
    case 1 % wrt u(t-\tau_1) and v(t-\tau_1)
        J = [par(2)*par(3)*sigmau*(1-sigmau), 0;
             0, par(1)*par(2)*par(6)*sigmav*(1-sigmav)];
    case 2 % wrt u(t-\tau_2) and v(t-\tau_2)
        J = [0, par(2)*par(4)*sigmau*(1-sigmau);
             par(1)*par(2)*par(5)*sigmav*(1-sigmav), 0];
    end
elseif isempty(nx) && isscalar(np) && isempty(v)
    switch np
    case 1 % wrt \alpha
        J = [0;
             -xx(2,1)+sigmav];
    case 2 % wrt \beta
        J = [zu*(sigmau-sigmau^2);
             par(1)*zv*(sigmav-sigmav^2)];
    case 3 % wrt a
        J = [par(2)*xx(1,2)*(sigmau-sigmau^2);
             0];
    case 4 % wrt b
        % J = [par(2)*xx(2,3)*(sigmau-sigmau^2);
        %      0];
        J = [par(2)*xx(2,2)*(sigmau-sigmau^2);
             0];
    case 5 % wrt c
        % J = [0;
        %      par(1)*par(2)*xx(1,3)*(sigmav-sigmav^2)];
        J = [0;
             par(1)*par(2)*xx(1,2)*(sigmav-sigmav^2)];
    case 6 % wrt d
        J = [0;
             par(1)*par(2)*xx(2,2)*(sigmav-sigmav^2)];
    case 7 % wrt \theta_u
        J = [par(2)*(sigmau-sigmau^2);
             0];
    case 8 % wrt \theta_v
        J = [0;
             par(1)*par(2)*(sigmav-sigmav^2)];
    case 9 % wrt \tau_1
        J = [0;
             0];
    case 10 % wrt \tau_2
        J = [0;
             0];
    end
elseif isscalar(nx) && isscalar(np) && isempty(v)
    switch nx
    case 0 % wrt u(t) and v(t)
        switch np
        case 1 % wrt \alpha
            J = [0, 0;
                 0, -1];
        case 2 % wrt \beta
            J = [0, 0;
                 0, 0];
        case 3 % wrt a
            J = [0, 0;
                 0, 0];
        case 4 % wrt b
            J = [0, 0;
                 0, 0];
        case 5 % wrt c
            J = [0, 0;
                 0, 0];
        case 6 % wrt d
            J = [0, 0;
                 0, 0];
        case 7 % wrt \theta_u
            J = [0, 0;
                 0, 0];
        case 8 % wrt \theta_v
            J = [0, 0;
                 0, 0];
        case 9 % wrt \tau_1
            J = [0, 0;
                 0, 0];
        case 10 % wrt \tau_2
            J = [0, 0;
                 0, 0];
        end
    case 1 % wrt u(t-\tau_1) and v(t-\tau1)
        switch np
        case 1 % wrt \alpha
            J = [0, 0;
                 0, par(2)*par(6)*sigmav*(1-sigmav)];
            case 2 % wrt \beta
            J = [par(3)*sigmau*(1-sigmau)*(par(2)*zu*(1-2*sigmau)+1), 0;
                 0, par(1)*par(6)*sigmav*(1-sigmav)*(par(2)*zv*(1-2*sigmav)+1)];
        case 3 % wrt a
            J = [par(2)*sigmau*(1-sigmau)*(par(2)*par(3)*xx(1,2)*(1-2*sigmau)+1), 0;
                 0, 0];
        case 4 % wrt b
            % J = [par(2)^2*par(3)*xx(2,3)*sigmau*(1-sigmau)*(1-2*sigmau), 0;
            %      0, 0];
            J = [par(2)^2*par(3)*xx(2,2)*sigmau*(1-sigmau)*(1-2*sigmau), 0;
                 0, 0];
        case 5 % wrt c
            % J = [0, 0;
            %      0, par(1)*par(2)^2*par(6)*xx(1,3)*sigmav*(1-sigmav)*(1-2*sigmav)];
            J = [0, 0;
                 0, par(1)*par(2)^2*par(6)*xx(1,2)*sigmav*(1-sigmav)*(1-2*sigmav)];
        case 6 % wrt d
            J = [0, 0;
                 0, par(1)*par(2)^2*par(6)*xx(2,2)*sigmav*(1-sigmav)*(1-2*sigmav)+par(1)*par(2)*sigmav*(1-sigmav)];
        case 7 % wrt \theta_u
            J = [par(2)^2*par(3)*sigmau*(1-sigmau)*(1-2*sigmau), 0;
                 0, 0];
        case 8 % wrt \theta_v
            J = [0, 0;
                 0, par(1)*par(2)^2*par(6)*sigmav*(1-sigmav)*(1-2*sigmav)];
        case 9 % wrt \tau_1
            J = [0, 0;
                 0, 0];
        case 10 % wrt \tau_2
            J = [0, 0;
                 0, 0];
        end
	case 2 % wrt u(t-\tau2) and v(t-\tau2)
        switch np
        case 1 % wrt \alpha
            J = [0, 0;
                 par(2)*par(5)*sigmav*(1-sigmav), 0];
        case 2 % wrt \beta
            J = [0, par(4)*sigmau*(1-sigmau)*(par(2)*zu*(1-2*sigmau)+1);
                 par(1)*par(5)*sigmav*(1-sigmav)*(par(2)*zv*(1-2*sigmav)+1), 0];
        case 3 % wrt a
            J = [0, par(2)^2*par(4)*xx(1,2)*sigmau*(1-sigmau)*(1-2*sigmau);
                 0, 0];
        case 4 % wrt b
            % J = [0, par(2)^2*par(4)*xx(2,3)*sigmau*(1-sigmau)*(1-2*sigmau)+par(2)*sigmau*(1-sigmau);
            %      0, 0];
            J = [0, par(2)^2*par(4)*xx(2,2)*sigmau*(1-sigmau)*(1-2*sigmau)+par(2)*sigmau*(1-sigmau);
                 0, 0];
        case 5 % wrt c
            % J = [0, 0;
            %      par(1)*par(2)^2*par(5)*xx(1,3)*sigmav*(1-sigmav)*(1-2*sigmav)+par(1)*par(2)*sigmav*(1-sigmav), 0];
            J = [0, 0;
                 par(1)*par(2)^2*par(5)*xx(1,2)*sigmav*(1-sigmav)*(1-2*sigmav)+par(1)*par(2)*sigmav*(1-sigmav), 0];

        case 6 % wrt d
            J = [0, 0;
                 par(1)*par(2)^2*par(5)*xx(2,2)*sigmav*(1-sigmav)*(1-2*sigmav), 0];
        case 7 % wrt \theta_u
            J = [0, par(2)^2*par(4)*sigmau*(1-sigmau)*(1-2*sigmau);
                 0, 0];
        case 8 % wrt \theta_v
            J = [0, 0;
                 par(1)*par(2)^2*par(5)*sigmav*(1-sigmav)*(1-2*sigmav), 0];
        case 9 % wrt \tau_1
            J = [0, 0;
                 0, 0];
        case 10 % wrt \tau2
            J = [0, 0;
                 0, 0];
        end
    end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
    nx1 = nx(1); 
    nx2 = nx(2);
    switch nx1
    case 0 % wrt u(t) and v(t)
        switch nx2
        case 0 % wrt u(t) and v(t)
            J = [0, 0;
                 0, 0];
        case 1 % wrt u(t-\tau_1) and v(t-\tau_1)
            J = [0, 0;
                 0, 0];
        case 2 % wrt u(t-\tau2) and v(t-\tau_2)
            J = [0, 0;
                 0, 0];
        end
    case 1 % wrt u(t-\tau_1) and v(t-\tau_1)
        switch nx2
        case 0 % wrt u(t) and v(t)
            J = [0, 0;
                 0, 0];
        case 1 % wrt u(t-\tau_1) and v(t-\tau_1)
            J = [par(2)^2*par(3)^2*sigmau*(1-sigmau)*(1-2*sigmau), 0;
                 0, par(1)*par(2)^2*par(6)^2*sigmav*(1-sigmav)*(1-2*sigmav)];
        case 2 % wrt u(t-\tau2) and v(t-\tau_2)
            J = [0, 0;
                 0, 0];
        end
    case 2 % wrt u(t-\tau2) and v(t-\tau_2)
        switch nx2
        case 0 % wrt u(t) and v(t)
            J = [0, 0;
                 0, 0];
        case 1 % wrt u(t-\tau_1) and v(t-\tau_1)
            J = [0, 0; 
                 0, 0];
        case 2 % wrt u(t-\tau2) and v(t-\tau_2)
            J = [0, par(2)^2*par(4)^2*sigmau*(1-sigmau)*(1-2*sigmau);
                 par(1)*par(2)^2*par(5)^2*sigmav*(1-sigmav)*(1-2*sigmav),0];
        end
    end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end
end
