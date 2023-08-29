function v = window_adapted(OM, denom, shift, b, Phi)
    v = left((OM+shift)/denom*b) + middle(OM, Phi) + right((OM-shift)/denom*b);
end

function r = right(OM)
    r = zeros(size(OM));
    case1=(OM>=1/3) & (OM<=2/3); r(case1) = cos(pi/2*ny(3*abs(OM(case1))-1));
end					    

function l = left(OM)
    l = zeros(size(OM));
    case1=(OM<=-1/3) & (OM>=-2/3); l(case1) = cos(pi/2*ny(3*abs(OM(case1))-1));
end

function m = middle(OM, Phi)
    m = zeros(size(OM));
    if isempty(Phi)
        case3=(abs(OM)<=1/3);
    else
        case3=(abs(OM)<=Phi);
    end
    m(case3) = 1;
end

function y = ny(x)
    y=zeros(size(x));
    case2=(x>0)&(x<1); y(case2)= s(x(case2)-1)./(s(x(case2)-1)+s(x(case2)));
    case3=(x>=1); y(case3)=1;
end

function y = s(x)
    y=exp(-((1+x).^(-2)+(1-x).^(-2)));
end
