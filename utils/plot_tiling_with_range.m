function plot_tiling_with_range(W,Phi,fig_title)

figure,
if isa(W, 'cell')
    if length(W) ~= 2
        H = zeros(size(W{1,1}{1,1}));
        V = zeros(size(W{1,1}{1,1}));
        for scale=1:length(W)
            for ell=1:length(W{1,scale})
                w = W{1,scale}{1,ell};
                H = H + w;
                V = V + w.^2;

%                 subplot(121), imagesc(V)
%                 subplot(122), imagesc(H)
%                 drawnow
%                 pause(0.4)
            end
        end
    else
        H = zeros(size(W{1,1}{1,1}{1,1}));
        V = zeros(size(W{1,1}{1,1}{1,1}));
        for vi=1:length(W)
            for scale=1:length(W{1,vi})
                for ell=1:length(W{1,vi}{1,scale})
                    w = W{1,vi}{1,scale}{1,ell};
                    H = H + w;
                    V = V + w.^2;
% 
%                     subplot(121), imagesc(V)
%                     subplot(122), imagesc(H)
%                     drawnow
%                     pause(0.4)
                end
            end
        end
    end
axis off
axis equal
imagesc(H)

NN = size(H,1);
N = NN/2;
[x1, y1] = pol2cart(Phi, N);
[x2, y2] = pol2cart(-Phi, N);
line([NN/2+1, x1*NN], [NN/2+1, y1*NN], 'Color', 'red')
line([NN/2+1, -x1*NN], [NN/2+1, -y1*NN], 'Color', 'red')
line([NN/2+1, x2*NN], [NN/2+1, y2*NN], 'Color', 'red')
line([NN/2+1, -x2*NN], [NN/2+1, -y2*NN], 'Color', 'red')
else
    NN = size(W,1);
    N = NN/2;
    
    [x1, y1] = pol2cart(Phi, N);
    [x2, y2] = pol2cart(-Phi, N);
    
    axis off
    axis equal
    imagesc(W)
    line([NN/2+1, x1*NN], [NN/2+1, y1*NN], 'Color', 'red')
    line([NN/2+1, -x1*NN], [NN/2+1, -y1*NN], 'Color', 'red')
    line([NN/2+1, x2*NN], [NN/2+1, y2*NN], 'Color', 'red')
    line([NN/2+1, -x2*NN], [NN/2+1, -y2*NN], 'Color', 'red')
    % % imagesc(abs(fftshift(fft2(frec))))
    % L = sqrt(N^2+N^2);
    % x=NN/2+1;
    % y=NN/2+1;
    % x2=x+(L*cos(Phi));
    % y2=y+(L*sin(Phi));
    % plot([x x2],[y y2], 'r')
    % x2=x-(L*cos(Phi));
    % y2=y-(L*sin(Phi));
    % plot([x x2],[y y2], 'r')
    % x2=x+(L*cos(-Phi));
    % y2=y+(L*sin(-Phi));
    % plot([x x2],[y y2], 'r')
    % x2=x-(L*cos(-Phi));
    % y2=y-(L*sin(-Phi));
    % plot([x x2],[y y2], 'r')
    if ~isempty(fig_title)
        title(fig_title)
    end
end
end

