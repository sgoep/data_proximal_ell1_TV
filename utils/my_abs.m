function y = my_abs(x)
    if isa(x, 'cell')
        if length(x) ~= 2
            for scale=1:length(x)
                for ell=1:length(x{1,scale})
                    tmp = x{1,scale}{1,ell};
                    y{1,scale}{1,ell} = abs(tmp);
                end
            end
        else
            for vi=1:length(x)
                for scale=1:length(x{1,vi})
                    for ell=1:length(x{1,vi}{1,scale})
                        tmp = x{1,vi}{1,scale}{1,ell};
                        y{1,vi}{1,scale}{1,ell} = abs(tmp);
                    end
                end
            end
        end
    else
        y = abs(x);
    end

end
