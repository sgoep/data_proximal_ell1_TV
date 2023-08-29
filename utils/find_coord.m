function coord = find_coord(w, j)

N = size(w, 1);

Ka = min(N, 2^(j+1));

% Ka = min(N/2, 2^(j+1));
ymin = N/2-Ka/2+1;
ymax = N/2+Ka/2;
xmin = N/2-Ka/2+1;
xmax = N/2+Ka/2;


% [yId, xId] = find(w);
% xmin = min(xId);
% ymin = min(yId);
% xmax = max(xId);
% ymax = max(yId);
% % coord = w(ymin:ymax,xmin:xmax);
coord = [xmin, ymin, xmax, ymax];

end

