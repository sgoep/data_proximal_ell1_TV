function plot_tiling(W)

if length(W) ~= 2
    img = zeros(size(W{1,1}{1,1}));
    for scale=1:length(W)
        for ell=1:length(W{1,scale})
            img = img + W{1,scale}{1,ell};
        end
    end
else
    img = zeros(size(W{1,1}{1,1}{1,1}));
    for vi=1:length(W)
        for scale=1:length(W{1,vi})
            for ell=1:length(W{1,vi}{1,scale})
                img = img + W{1,vi}{1,scale}{1,ell};
            end
        end
    end
end

figure, imagesc(img)

end

