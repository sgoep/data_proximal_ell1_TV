function R = transform(x, W, COORD, N, subsamp, direction)

if strcmp(direction, 'analysis_loc')
    x = x(1:end-1, 1:end-1);
end
func = strcat(direction, '(x, W, COORD, N, subsamp)');
R = eval(func);

end

function C = analysis_loc(x, W, COORD, N, subsamp)
    H = zeros(size(x));
%     figure
    if subsamp
        fhat = fftshift(fft2(ifftshift(x))); 
    else
        fhat = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)));
    end
%     fhat = x/sqrt(prod(size(x)));
    if length(W) ~= 2
        for scale=1:length(W)
            for ell=1:length(W{1,scale})
                w = W{1,scale}{1,ell};
             
                c = w.*fhat;
                if subsamp
                    coord = COORD{1,scale}{1,ell};
                    ca = ((coord(3) - coord(1))*(coord(4) - coord(2)))/N^2;
                    cnew = c(coord(2):coord(4),coord(1):coord(3))*ca;
                    C{1,scale}{1,ell} = fftshift(ifft2(ifftshift(cnew)));
                else
                    cnew = c;
                    C{1,scale}{1,ell} = (fftshift(ifft2(ifftshift(cnew)))*sqrt(prod(size(cnew))));
                end

                
%                 H = H + cnew;
%                 imagesc(abs(H)), drawnow, pause(0.5)
            end
        end
    else
        for vi=1:length(W)
            for scale=1:length(W{1,vi})
                for ell=1:length(W{1,vi}{1,scale})
                    w = W{1,vi}{1,scale}{1,ell};

                    c = w.*fhat;
                    if subsamp
                        coord = COORD{1,vi}{1,scale}{1,ell};
                        ca = ((coord(3) - coord(1))*(coord(4) - coord(2)))/N^2;
                        cnew = c(coord(2):coord(4),coord(1):coord(3))*ca;
                        C{1,vi}{1,scale}{1,ell} = fftshift(ifft2(ifftshift(cnew)));
                    else
                        cnew = c;
                        C{1,vi}{1,scale}{1,ell} = (fftshift(ifft2(ifftshift(cnew)))*sqrt(prod(size(cnew))));
                    end

                    
%                     H = H + w;
%                     imagesc(H), drawnow, pause(0.2)
                end
            end
        end
    end
end

function f = synthesis_loc(x, W, COORD, N, subsamp)
%     fhat = zeros(size(W{1,1}{1,1}));
%     N = size(fhat, 1);
%     figure
    fhat = zeros(N);
    if length(W) ~= 2
        for scale=1:length(W)
            for ell=1:length(W{1,scale})
                if subsamp
                    c = zeros(N-1);
                    coord = COORD{1,scale}{1,ell};
                    calt = fftshift(fft2(x{1,scale}{1,ell}));
                    c(coord(2):coord(4),coord(1):coord(3)) = calt;
                    ca = N^2/prod(size(calt));
                    c = c*ca;
                else
                    c = fftshift(fft2(ifftshift((x{1,scale}{1,ell}))))/sqrt(prod(size(x{1,scale}{1,ell})));
                end
                w = W{1,scale}{1,ell};
%                 figure, imagesc(abs(c))
                fhat = fhat + w.*c;
            end
        end
    else
        for vi=1:length(W)
            for scale=1:length(W{1,vi})
                for ell=1:length(W{1,vi}{1,scale})
                    if subsamp
                        c = zeros(N);
                        coord = COORD{1,vi}{1,scale}{1,ell};
                        calt = fftshift(fft2(x{1,vi}{1,scale}{1,ell}));
                        c(coord(2):coord(4),coord(1):coord(3)) = calt;
                        ca = N^2/prod(size(calt));
                        c = c*ca;
                    else
                        c =  fftshift(fft2(ifftshift((x{1,vi}{1,scale}{1,ell}))))/sqrt(prod(size(x{1,vi}{1,scale}{1,ell})));
                    end
                    w = W{1,vi}{1,scale}{1,ell};
                    fhat = fhat + w.*c;
%                     imagesc(abs(fhat)), drawnow
                end
            end
        end
    end   
    if subsamp
        f = real((ifft2(ifftshift(fhat))));%*sqrt(prod(size(fhat)));
    else
        f = real(ifftshift(ifft2(ifftshift(fhat))))*sqrt(prod(size(fhat)));
    end
    f = padarray(f, [1,1], 0, 'post');
end