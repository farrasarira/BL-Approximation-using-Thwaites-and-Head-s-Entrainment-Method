function [theta_turb,H_turb,Cf_turb,disp_turb,vtran,sep_u,sep_l] = TURBULENT(vtan,X,i_s,nu,trans_u,trans_l,theta_thwaites,H_thwaites,Cf_thwaites,disp_thwaites)
%TURBULENT Summary of this function goes here
%   Detailed explanation goes here
    numPan = size(X,1);
    vtan = abs(vtan);
    theta_turb = theta_thwaites;
    H_turb = H_thwaites;
    H_1 = H_turb;
    Cf_turb = Cf_thwaites;
    disp_turb = disp_thwaites;
    vtran = size(X,1);
    stag_u = i_s;
    stag_l = i_s-1;

    for i=trans_u:numPan
        if i == numPan
            Re_theta = vtan(i)*theta_turb(i)/nu;
            Cf_turb(i) = 0.246/(10^(0.678*H_turb(i)) * Re_theta^0.268);
            disp_turb(i) = theta_turb(i) * H_turb(i);
            break;
        end
        % Calculate Re_theta
        Re_theta = vtan(i)*theta_turb(i)/nu;
        Cf_turb(i) = 0.246/(10^(0.678*H_turb(i)) * Re_theta^0.268);
        
        % Calculate theta
        dVin_dx = (vtan(i+1) - vtan(i))/(X(i+1) - X(i));
        dtheta_dx = Cf_turb(i)/2 - (2+H_turb(i))*theta_turb(i)/vtan(i) * dVin_dx;
        theta_turb(i+1) = theta_turb(i) + dtheta_dx*(X(i+1) - X(i));
          
        % Calculate H and H_1
%         if H_turb(i) <= 1.6
%             H_1(i) = 0.8234*(H_turb(i)-1.1)^-1.287 + 3.3;
%         else
%             H_1(i) = 1.5501*(H_turb(i)-0.6778)^-3.064 + 3.3;
%         end
        H_1(i) = 3.0445 + 0.8702/(H_turb(i)-1.1)^1.2721;
        
        H_1(i+1) = (vtan(i+1)*theta_turb(i+1))^-1  *  (vtan(i)*theta_turb(i)*H_1(i) + (0.0306*vtan(i))/(H_1(i)-3.0)^-0.6169 * (X(i+1)-X(i)) );
        H_turb(i+1) = 1.1 + ( 0.8702 / (H_1(i+1) - 3.0445) )^(1/1.2721);

        % Calculate disp_t and 
        disp_turb(i) = theta_turb(i) * H_turb(i);
    end

    for i=trans_l:-1:1
        if i == 1
            Re_theta = vtan(i)*theta_turb(i)/nu;
            Cf_turb(i) = 0.246/(10^(0.678*H_turb(i)) * Re_theta^0.268);
            disp_turb(i) = theta_turb(i) * H_turb(i);
            break;
        end
        % Calculate Re_theta
        Re_theta = vtan(i)*theta_turb(i)/nu;
        Cf_turb(i) = 0.246/(10^(0.678*H_turb(i)) * Re_theta^0.268);
        
        % Calculate theta
        dVin_dx = (vtan(i-1) - vtan(i))/(X(i-1) - X(i));
        dtheta_dx = Cf_turb(i)/2 - (2+H_turb(i))*theta_turb(i)/vtan(i) * dVin_dx;
        theta_turb(i-1) = theta_turb(i) + dtheta_dx*(X(i-1) - X(i));
          
        % Calculate H and H_1
%         if H_turb(i) <= 1.6
%             H_1(i) = 0.8234*(H_turb(i)-1.1)^-1.287 + 3.3;
%         else
%             H_1(i) = 1.5501*(H_turb(i)-0.6778)^-3.064 + 3.3;
%         end
        H_1(i) = 3.0445 + 0.8702/(H_turb(i)-1.1)^1.2721;
        
        H_1(i-1) = (vtan(i-1)*theta_turb(i-1))^-1  *  (vtan(i)*theta_turb(i)*H_1(i) + (0.0306*vtan(i))/(H_1(i)-3.0)^-0.6169 * (X(i-1)-X(i)) );
        H_turb(i-1) = 1.1 + ( 0.8702 / (H_1(i-1) - 3.0445) )^(1/1.2721);

        % Calculate disp_t and 
        disp_turb(i) = theta_turb(i) * H_turb(i);
    end

    %% Calculate Transpirational Velocity
    vtran(stag_u) = (vtan(stag_u+1)*disp_turb(stag_u+1)-vtan(stag_u)*disp_turb(stag_u))/(X(stag_u+1)-X(stag_u));
    vtran(stag_l) = (vtan(stag_l-1)*disp_turb(stag_l-1)-vtan(stag_l)*disp_turb(stag_l))/(X(stag_l-1)-X(stag_l));

    for i=stag_u+1:numPan
        vtran(i) = (vtan(i)*disp_turb(i)-vtan(i-1)*disp_turb(i-1))/(X(i)-X(i-1));
    end
    for i=stag_l-1:-1:1
        vtran(i) = (vtan(i)*disp_turb(i)-vtan(i+1)*disp_turb(i+1))/(X(i)-X(i+1));
    end

    %% Determine Separation point
    %Seperation occurs when H = 0.35
    sep_u = numPan;
    sep_l = 1;
    for i=stag_u:numPan
        if(H_turb(i)) >= 3.5
            sep_u = i;
            break;
        end
    end

    for i=stag_l:-1:1
        if(H_turb(i)) >= 3.5
            sep_l = i;
            break;
        end
    end
    
end

