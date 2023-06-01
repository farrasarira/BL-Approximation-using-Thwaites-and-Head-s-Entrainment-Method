function [sep_u,sep_l,theta,disp_t,H,Cf,trans_u,trans_l,vtran] = THWAITES(vtan,X,i_s,nu)
    numPan = size(X,1);
    vtan = abs(vtan);
    theta = zeros(numPan,1);
    stag_u = i_s;
    stag_l = i_s-1;
    m = zeros(numPan,1);
    l = zeros(numPan,1);
    H = zeros(numPan,1);
    disp_t = zeros(numPan,1);
    Cf = zeros(numPan,1);
    vtran = zeros(numPan,1);

    

    %% Calculate Theta
    theta(stag_u) = sqrt(0.0735*nu/(vtan(stag_u+1)-vtan(stag_u))*abs(X(stag_u+1)-X(stag_u)));
    for i=stag_u+1:numPan
        int1 = 0;
        for j=stag_u+1:i
            int1 = int1 + (vtan(j)^5+vtan(j-1)^5)*abs(X(j)-X(j-1))/2;
        end
        theta(i) = sqrt(0.441*nu*int1/vtan(i)^6);
    end

    theta(stag_l) = sqrt(0.0735*nu/(vtan(stag_l-1)-vtan(stag_l))*(X(stag_l-1)-X(stag_l)));
    for i=stag_l-1:-1:1
        int1 = 0;
        for j=stag_l-1:-1:i
            int1 = int1 + (vtan(j)^5+vtan(j+1)^5) * abs(X(j)-X(j+1))/2;
        end
        theta(i) = sqrt(0.441*nu*int1/vtan(i)^6);
    end

    %% Calculate paramter m
    for i=stag_u+1:numPan
        m(i) = 1 * theta(i)^2/nu*(vtan(i)-vtan(i-1))/(X(i)-X(i-1));
    end
    for i=stag_l-1:-1:1
        m(i) = 1 * theta(i)^2/nu*(vtan(i)-vtan(i+1))/(X(i)-X(i+1));
    end
   
    m(stag_u) = 0.075;
    m(stag_l) = 0.075;
    % figure(15);
    % plot(X,m,'*')
    %% Calculate l(m), H(m), displacement thickeness, and separation point
    
    m_data = [-0.250, -0.200, -0.140, -0.120, -0.100, -0.080, -0.064, -0.048, -0.032, -0.016, 0.000, 0.016, 0.032, 0.040, 0.048, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080, 0.084, 0.086, 0.088, 0.090];
    l_data = [ 0.500,  0.463,  0.404,  0.382,  0.359,  0.333,  0.313,  0.291,  0.268,  0.244, 0.220, 0.195, 0.168, 0.153, 0.138, 0.122, 0.113, 0.104, 0.095, 0.085, 0.072, 0.056, 0.038, 0.027, 0.015, 0.000];
    H_data = [ 2.000,  2.070,  2.180,  2.230,  2.280,  2.340,  2.390,  2.440,  2.490,  2.550, 2.610, 2.670, 2.750, 2.810, 2.870, 2.940, 2.990, 3.040, 3.090, 3.150, 3.220, 3.300, 3.390, 3.440, 3.490, 3.550];

    
    
    for i=stag_u:numPan
        if (m(i) >= -0.1) && (m(i) < 0)
            l(i) = 0.22 + 1.420*m(i)+ 0.018*m(i)/(0.107+m(i));
            H(i) = 0.0731/(0.14+m(i))+2.088;
        elseif (m(i) >= 0) && (m(i) <= 0.25)
            l(i) = 0.22 + 1.57*m(i) - 1.8*m(i)^2;
            H(i) = 2.61 - 3.75*m(i) + 5.24*m(i)^2;
        else
            l(i) = l(i-1);
            H(i) = H(i-1);
        end

%         l(i) = (m(i) + 0.09)^0.62;
%         H(i) = 2 +4.14*(0.25-m(i)) -83.5*(0.25-m(i))^2 +854*(0.25-m(i))^3 -3337*(0.25-m(i))^4 +4576*(0.25-m(i))^5;      

%         l(i) = interp1(m_data, l_data, m(i), 'spline');
%         H(i) = interp1(m_data, H_data, m(i), 'spline');
        
        disp_t(i) = H(i)*theta(i);
        if theta(i) == 0
            Cf(i) = 0;
        else
            Cf(i) = 2*nu/(vtan(i)*theta(i))*l(i);
        end
    end
    
    for i=stag_l:-1:1
        if (m(i) >= -0.1) && (m(i) < 0)
            l(i) = 0.22 + 1.420*m(i)+ 0.018*m(i)/(0.107+m(i));
            H(i) = 0.0731/(0.14+m(i))+2.088;
        elseif (m(i) >= 0) && (m(i) <= 0.25)
            l(i) = 0.22 + 1.57*m(i) - 1.8*m(i)^2;
            H(i) = 2.61 - 3.75*m(i) + 5.24*m(i)^2;
        else
            l(i) = l(i+1);
            H(i) = H(i+1);
        end
        
%         l(i) = (m(i) + 0.09)^0.62;
%         H(i) = 2 +4.14*(0.25-m(i)) -83.5*(0.25-m(i))^2 +854*(0.25-m(i))^3 -3337*(0.25-m(i))^4 +4576*(0.25-m(i))^5;
        
        disp_t(i) = H(i)*theta(i);
        if theta(i) == 0
            Cf(i) = 0;
        else
            Cf(i) = 2*nu/(vtan(i)*theta(i))*l(i);
        end
    end

    %% Determine Transition Point using Michel Method
    %Transition occurs when Re_theta >= R_transX
    dtrans_u=0;
    for i=stag_u:numPan
        Re_theta = vtan(i)*theta(i)/nu;
        Re_X = vtan(i)*X(i)/nu;
        Re_transX = 2.9*Re_X^0.4;

        if Re_theta > Re_transX
            trans_u = i+dtrans_u;
            break;
        end
    end
    
    dtrans_l=0;
    for i=stag_l:-1:1
        Re_theta = vtan(i)*theta(i)/nu;
        Re_X = vtan(i)*X(i)/nu;
        Re_transX = 2.9*Re_X^0.4;

        if Re_theta > Re_transX
            trans_l = i-dtrans_l;
            break;
        end
    end

    %% Determine Separation point
    %Seperation occurs when H = 0.35
    sep_u = numPan;
    sep_l = 1;
    for i=stag_u:numPan
        if(H(i)) >= 3.5
            sep_u = i;
            break;
        end
    end

    for i=stag_l:-1:1
        if(H(i)) >= 3.5
            sep_l = i;
            break;
        end
    end
    
    %% Calculate Transpirational Velocity
    vtran(stag_u) = (vtan(stag_u+1)*disp_t(stag_u+1)-vtan(stag_u)*disp_t(stag_u))/(X(stag_u+1)-X(stag_u));
    vtran(stag_l) = (vtan(stag_l-1)*disp_t(stag_l-1)-vtan(stag_l)*disp_t(stag_l))/(X(stag_l-1)-X(stag_l));

    for i=stag_u+1:numPan
        vtran(i) = (vtan(i)*disp_t(i)-vtan(i-1)*disp_t(i-1))/(X(i)-X(i-1));
    end
    for i=stag_l-1:-1:1
        vtran(i) = (vtan(i)*disp_t(i)-vtan(i+1)*disp_t(i+1))/(X(i)-X(i+1));
    end
end

