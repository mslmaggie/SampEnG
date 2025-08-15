clear all

%% Values for logistic map

% The logistic map becomes chaotic for 3.57 < r < 4.
% This code sweeps r from 3.55 to 4.00 and computes complexity measures.
% r is the growth rate parameter in the logistic map.

a1 = 3.55;       % lower bound of r sweep
a2 = 4.00;       % upper bound of r sweep
def = 10e-5;     % step size for r (here: 1e-4)

auxi = length(a1:def:a2);  % number of r values to test
res = zeros(auxi, 3);      % matrix to store results [SampEn, SEG_undir, SEG_dir]

in = 2;        % embedding dimension for analysis (m in Sample Entropy)

siz = 2^12;    % total length of generated time series (4096 points)

j = 1;         % index for storing results

k = 2^10;      % number of transient iterations to discard (1024 points)

N = siz - k + 1;  % length of time series after discarding transient

% Construct adjacency-like matrices for graph-based complexity analysis
DPath = sparse(diag(ones(1, N-1), 1));          % path connection
UPath_undir = sparse(transpose(DPath) + DPath); % undirected path adjacency
UPath_dir = sparse([zeros(N-1,1), diag(ones(N-1,1)); zeros(1,N)]); % directed path adjacency

%% Loop over r values
for r = a1:def:a2
    % Generate logistic map series
    xe = zeros(siz, 1);        % preallocate series
    xe(1, 1) = 0.65;           % initial condition
    for i = 1:siz-1
        xe(i+1, 1) = r * xe(i, 1) * (1 - xe(i, 1)); % logistic map equation
    end

    % Compute complexity metrics on the time series after transient
    res(j, 1) = sampen(xe(k:end, 1)', in, 0.2);                  % Sample Entropy
    res(j, 2) = SEG(xe(k:end, 1)', UPath_undir, in, 0.2);     % SEG with undirected path
    res(j, 3) = SEG(xe(k:end, 1)', UPath_dir, in, 0.2);       % SEG with directed path

    j = j + 1;     % increment result index

    disp(r)        % display current r in console
end

% Save results to file for plotting
save('Log2.mat', 'res')
