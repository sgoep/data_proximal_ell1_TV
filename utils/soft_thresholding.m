function y = SoftThresh(x,alpha)
    if isa(x, 'cell')
        if length(x) ~= 2
            for scale=1:length(x)
                if length(alpha) ~= 1
                    if scale > 1
                        lambda = alpha(scale-1);
                    else
                        lambda = 0;
                    end
                else
                    if scale > 1
                        lambda = alpha;
                    else
                        lambda = 0;
                    end
                end

                for ell=1:length(x{1,scale})
                    tmp = x{1,scale}{1,ell};
                    y{1,scale}{1,ell} = sign(tmp).*(max(0, abs(tmp) - lambda));
                end
            end
        else
            for vi=1:length(x)
                for scale=1:length(x{1,vi})
                    if length(alpha) ~= 1 
                        if scale > 1
                            lambda = alpha(scale-1);
                        else
                            lambda = 0;
                        end
                    else
                        if scale > 1
                            lambda = alpha;
                        else
                            lambda = 0;
                        end
                    end

                    for ell=1:length(x{1,vi}{1,scale})
                        tmp = x{1,vi}{1,scale}{1,ell};
                        y{1,vi}{1,scale}{1,ell} = sign(tmp).*(max(0, abs(tmp) - lambda));
                    end
                end
            end
        end
    else
        y = sign(x).*(max(0,abs(x)-alpha));
    end

end