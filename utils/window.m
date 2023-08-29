function w = window(r, omega, j, type, shape)

    if strcmp(type, 'rs')
        w=ws_meyer(r*2^(-j));
    elseif strcmp(type, 'radial')
        w=wr_meyer(r*2^(-j));
    elseif strcmp(type, 'angular')
        if strcmp(shape, 'meyer_adapted')
            w=wa_meyer_adapted(omega);
        elseif strcmp(shape, 'meyer_inv')
            w=wa_meyer_adapted_inv(omega);
        elseif strcmp(shape, 'meyer')
            w=wa_meyer(omega);
        elseif strcmp(shape, 'shannon')
            w=wa_shannon(omega);
        else
            error('not implemented')
        end
    else
        disp('This type of window is unknown')
        w=zeros(size(r));
    end      
        
end

% Myer low-pass
function u = ws_meyer(R)
    u=zeros(size(R));
    u(R<4/3) = 1;
    case4=(R>=4/3) & (R<5/3); u(case4) = cos(pi/2*ny(3*R(case4)-4));
end

% Meyer W(r)
function u = wr_meyer(R)
    u=zeros(size(R));
    case2=(R>=2/3) & (R<=5/6); u(case2) = cos(pi/2*ny(5-6*R(case2)));
    case3=(R>=5/6) & (R<=4/3); u(case3) = 1;
    case4=(R>=4/3) & (R<=5/3); u(case4) = cos(pi/2*ny(3*R(case4)-4));
end

% Meyer V(omega)
function v = wa_meyer(OM)
    v =zeros(size(OM));
    case2=(abs(OM)>=1/3) & (abs(OM)<=2/3); v(case2) = cos(pi/2*ny(3*abs(OM(case2))-1));
    case3=(abs(OM)<=1/3); v(case3) = 1;
end

function v = wa_meyer_adapted(OM)
    v =zeros(size(OM));
    case2=(abs(OM)>=1/3) & (abs(OM)<=2/3); v(case2) = cos(pi/2*ny(3*abs(OM(case2))-1));
    case3=(abs(OM)<=1/3); v(case3) = 1;
end

function v = wa_meyer_adapted_inv(OM)
    v =zeros(size(OM));
    case3=(abs(OM*2)<=1/3); v(case3) = 1;
    case2=(abs(OM-1/3)>=1/3) & (abs(OM-1/3)<=2/3); v(case2) = cos(pi/2*ny(3*abs(OM(case2))-1));
end

% Shannon low-pass
function u = ws_shannon(R)
    u=zeros(size(R));
    u(R<=1.5) = 1;
end

% Shannon W(r)
function u = wr_shannon(R)
    u=zeros(size(R));
    case2=(R>0.75) & (R<=1.5); u(case2) = 1;
end

% Shannon V(omega)
function v = wa_shannon(OM)
    v =zeros(size(OM));
    case2=(OM<1/2 & OM>=-1/2); v(case2) = 1;
end
 
function y = ny(x)
    y=zeros(size(x));
    case2=(x>0)&(x<1); y(case2)= s(x(case2)-1)./(s(x(case2)-1)+s(x(case2)));
    case3=(x>=1); y(case3)=1;
    
    y = 3*x.^2 - 2*x.^3;
    %     y = x;
end

function y = s(x)
    y=exp(-((1+x).^(-2)+(1-x).^(-2)));
end



