function [W, COORD] = create_tiling(N, nscales, Phi, shape)

% This function create the Curvelet tiling in frequency domain.
%
% INPUT:
%   N       - Image size.
%   nscales - Number of scales in the Curvelet decomposition.
%   Phi     - Limited angular range parameter. If Phi is non empty, the tiling
%             will be adapted to the to the visibile and invisible cone.
%   shape   - To be removed.
%
% OUTPUT:
%   W       - Cell array containing the Curvelet tiles at given scales and
%             rotations.
%   COORD   - Cell array containing the support of the Curvelets in frequency
%             domain.


J = log2(N);
jmax = J-1; jmin = J-nscales;
k = linspace(-N/2, N/2, N);
[k1, k2] = ndgrid(k, k);

[omega, r] = cart2pol(k1, k2);
omega = omega';

% limited view adapted case
if ~isempty(Phi)
    
    Phi_vis = Phi;              % angular range for the visible cone
    Phi_inv = pi/2-Phi_vis;     % angular range for the invisible cone
    
    omega_vis_right = omega;
    omega_vis_left = fliplr(omega);

    omega_inv_top = omega - pi/2;
    omega_inv_bot = omega + pi/2;

    W_vis = {};
    W_inv = {};
    COORD_vis = {};
    COORD_inv = {};
    
    wr = window(r, [], jmin-1, 'rs', shape);
    w  = wr.*window([], omega_vis_right/(2*Phi_vis), jmin-1, 'angular', shape);
    w  = w + wr.*window([], omega_vis_left/(2*Phi_vis), jmin-1, 'angular', shape);

    W_vis{1,1}{1,1} = w;
    COORD_vis{1,1}{1,1} = find_coord(w, jmin);
    
    w  = wr.*window_adapted(omega_inv_top, 2*Phi_inv, (Phi_inv-Phi_vis), Phi_inv/Phi_vis, Phi_inv-1/3*Phi_vis);
    w  = w + wr.*window_adapted(omega_inv_bot, 2*Phi_inv, (Phi_inv-Phi_vis), Phi_inv/Phi_vis, Phi_inv-1/3*Phi_vis);
   
    W_inv{1,1}{1,1} = w;  
    COORD_inv{1,1}{1,1} = find_coord(w, jmin);
    
    scale = 2;
    for j = jmin:jmax
        wr = window(r, [], j, 'radial', shape);

        Ntheta = 2^(ceil((j-jmin)/2)+1);
        
        % Right and left visible part
        L = -Ntheta/2:Ntheta/2-1;
        rot = 1;
        for ell = L
            w = wr.*window([], omega_vis_right/(2*Phi_vis)*Ntheta-ell-0.5, [], 'angular', shape);
            W_vis{1,scale}{1,rot} = w;
            COORD_vis{1,scale}{1,rot} = find_coord(w, j+1);

            rot = rot + 1;
            
            w = wr.*window([], omega_vis_left/(2*Phi_vis)*Ntheta-ell-0.5, [], 'angular', shape);
            W_vis{1,scale}{1,rot} = w;
            COORD_vis{1,scale}{1,rot} = find_coord(w, j+1);
            
            rot = rot + 1;
        end
    
        % Top and bottom invisible part
        Ntheta2 = Ntheta;
        L = -Ntheta2/2:Ntheta2/2-1;
        rot = 1;
        for ell = L
            w = wr.*window_adapted(omega_inv_top*Ntheta2-(0.5+ell)*2*Phi_inv, 2*Phi_inv, (Phi_inv-Phi_vis), Phi_inv/Phi_vis, Phi_inv-1/3*Phi_vis);
          

            W_inv{1,scale}{1,rot} = w;
            COORD_inv{1,scale}{1,rot} = find_coord(w, j+1);

            rot = rot + 1;

            w = wr.*window_adapted(omega_inv_bot*Ntheta2-(0.5+ell)*2*Phi_inv, 2*Phi_inv, (Phi_inv-Phi_vis), Phi_inv/Phi_vis, Phi_inv-1/3*Phi_vis);
          
            W_inv{1,scale}{1,rot} = w;
            COORD_inv{1,scale}{1,rot} = find_coord(w, j+1);
    
            rot = rot + 1;
        end
    
        scale = scale + 1;
    end
    W = {W_vis, W_inv};
    COORD = {COORD_vis, COORD_inv};

% Regular full tiling case
else
    W  = {};
    w = window(r, [], jmin-1, 'rs', shape);
    COORD{1,1}{1,1} = find_coord(w, jmin);
    W{1,1}{1,1} = w;
    scale = 2;
    for j = jmin:jmax
        wr = window(r, [], j, 'radial', shape);
        Ntheta = 2^(ceil((j-jmin)/2));
    
        % Right and top visible/invisible part
        rot = 1;
        for ell = -4*Ntheta:4*Ntheta-1
            OM = mod(omega/(2*pi)*Ntheta*8-ell+1.5, 8*Ntheta)-1;
            w = wr.*window([], OM, [], 'angular', shape);
            W{1,scale}{1,rot} = w;
            COORD{1,scale}{1,rot} = find_coord(w, j+1);
            rot = rot + 1;
        end    
        scale = scale + 1;
    end

end

end
