function [r,omega] = rpc(k1,k2)
% Transformation from cartesian to recto-polar coordinates

    % radius
    r = max(abs(k1),abs(k2));

    % polar angle
    omega1=atan2(k1,k2);

    omega=zeros(size(omega1));
    case1=abs(omega1)<=pi/4; omega(case1)=tan(omega1(case1));
    case2=(omega1>=pi/4)&(omega1<3*pi/4); omega(case2)=2-cot(omega1(case2));
    case3=abs(omega1)>=3*pi/4; omega(case3) = 4+tan(omega1(case3));
    case4=(omega1>=-3*pi/4)&(omega1<-pi/4); omega(case4)=-2-cot(omega1(case4));

end
