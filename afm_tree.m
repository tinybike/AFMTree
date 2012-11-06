% afm_tree.m
% "Anti-ferromagnetic" model for tree growth

% Define an NxN grid (T), with each grid point big enough to hold one tree.
% Each element T(i,j) = 1 if there a tree has randomly started to grow
% there, and T(i,j) = 0 otherwise.
N = 100;
T = round(rand(N));

% Z is an NxN grid that has the number of neighboring trees of each site
% in the T grid
% F is a 3x3 "filter" matrix that contains the neighbors to be counted
% F = toeplitz([0 1 0]);     % no diagonal neighbors
F = [1 1 1;1 0 1;1 1 1];    % diagonals included!
Z = conv2(T,F,'same');

% Trees don't like to be next to one another.  Suppose the default growth
% rate (for a tree with no neighbors) is set at 1.  Further suppose the
% growth rate is decreased as exp(-z), where z is the number of neighbors
% that particular tree has.  K is a "growth rate matrix":
K = exp(-Z);
% K = 1./(1 + exp(Z - 2));

% To model the actual tree sizes, suppose they obey a logistic equation:
% S(t) = 2/(1 + exp(-k*t)) - 1
% where S(t) is the size of any given tree at time t, as a function of its
% growth rate k, given by the appropriate element of the matrix K
% Note: time (t) is stored as a matrix, since new trees will need to have 
% their time set to zero
tMax = 1000;
S = cell(tMax);
t = ones(N);
for i = 1:tMax
    S{i} = 2*(1./(1 + exp(-K.*t)) - 1/2).*T;
    %S{i} = 400./(1 + exp(-K.*(t/6 + 7))).*T;
    % Every time step, there's a 1% chance for a seed to land on each site,
    % so that a new tree grows there (starting at size 0)
    newTrees = rand(N) < 0.005;
    % Reset time for new trees
    t = t.*((~T + newTrees) < 2);
    % Update the tree matrix
    T = T + newTrees;
    T(T > 1) = 1;
    % Get new growth rate matrix
    Z = conv2(T,F,'same');
    K = exp(-Z);
    %K = 1./(1 + exp(Z - 2));
    % Every 10 time steps, each tree has a 10% chance to be killed
    if mod(i,10) == 0
        T = (rand(N) > 0.1).*T;
    end
    t = t + 1;
end

% Let's see what these look like...
% (for really big N, use shading interp to smooth out plots!)
nRows = 2;
nCols = 2;
figure(1)
clf

subplot(nRows,nCols,1)
surf(S{1})
zlabel('size')
title('t = 1')
view(2)
colorbar

subplot(nRows,nCols,2)
surf(S{50})
zlabel('size')
title('t = 50')
view(2)
colorbar

subplot(nRows,nCols,3)
surf(S{100})
zlabel('size')
title('t = 100')
view(2)
colorbar

subplot(nRows,nCols,4)
surf(S{1000})
zlabel('size')
title('t = 1000')
view(2)
colorbar

% Histograms
bins = 50;
[y50,x50] = hist(S{50}(:),bins);
[y100,x100] = hist(S{100}(:),bins);
[y1000,x1000] = hist(S{1000}(:),bins);

figure(2)
clf
plot(x50,y50,'r-')
hold on
plot(x100,y100,'b-')
plot(x1000,y1000,'k-')
plot(x50,y50,'r.')
plot(x100,y100,'b.')
plot(x1000,y1000,'k.')
xlabel('size')
ylabel('count')
legend('t=50','t=100','t=1000')
