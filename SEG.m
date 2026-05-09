function Out_SampEn = SEG(x, Adj, m, r, dist_type)

% SEG: Sample Entropy on graphs
% - x: a vector (length N) representing a signal on graph nodes
% - Adj: N x N adjacency (or weighted adjacency) matrix
% - m: embedding dimension (length of patterns)
% - r: tolerance for matching (scaled by std(x))
% - dist_type: distance metric for pdist (e.g., 'euclidean', 'chebychev', etc.)

% Error detection and defaults
if nargin < 4, error('Not enough parameters.'); end
if nargin < 5
    dist_type = 'chebychev';
    %fprintf('[WARNING] Using default distance method: chebychev.\n');
end
if ~isvector(x)
    error('The signal parameter must be a vector.');
end
if ~ischar(dist_type)
    error('Distance must be a string.');
end
if m > length(x)
    error('Embedding dimension must be smaller than the signal length (m<N).');
end

N = length(x); % number of vertices

if or(size(Adj, 1) ~= N, size(Adj, 1) ~= N)
    error('not square matrix');
else
    % ----------------------------------------------
    % Generate m+1 pattern matrix 
    % ----------------------------------------------
    L = 1;
    m1 = m+1;
   
    % Step 1: Matrix Lap where first column is the signal, and the columns the j-average
    Lap = zeros(N, m1);
    Lap(:, 1) = x';
    Aux = eye(N);
    t = 2;
    for j = L:L:(m1-1)*L
        Aux = Aux * (Adj^L);
        Lap(:, t) = sparse(diag(1./sum(Aux, 2))) * (Aux * x');
        t = t + 1;
    end

    % Tolerance
    sigma = std(x);
    tol = r * sigma;

    % ---------------------------------------------------
    % Count for
    % m = m and m = m+1 pattern matrices
    %---------------------------------------------------
    % m embedding 
    Xm_tmp = Lap(:, 1:m1);
    
    % Rows valid according to Xm
    valid_rows_m = ~any(isnan(Xm_tmp), 2);

    X_m  = Lap(valid_rows_m, 1:m);
    X_m1 = Lap(valid_rows_m, 1:m1);

    d_m = pdist(X_m(:,1:m), dist_type);             % Compute pairwise of m dimension vectors in rows

    Nm = size(X_m, 1);
        
 
    if isempty(d_m)
        % If B = 0, SampEn is not defined: no regularity detected
        %   Note: Upper bound is returned
        Out_SampEn = Inf;
    else      

        % Compute pairwise of m+1 dimension vectors in rows
        d_m1 = pdist(X_m1, dist_type);

        % Count for number of distances <= tolerance
        count_m = sum(d_m  <= tol);
        count_m1 = sum(d_m1 <= tol);

        % Output 0 if SampEn is undefined
        if count_m1 == 0
            Out_SampEn = Inf;
        else
            % Compute sample entropy
            Out_SampEn = -log(count_m1 / count_m);
        end
    end

    if isinf(Out_SampEn)
        % Note: SampEn has the following limits:
        %       - Upper bound: log(N-m)+log(N-m-1)-log(2)
        % Exactly one match in m+1
        % But dependent on graph size
        
        Out_SampEn = log(Nm)+log(Nm-1)-log(2);

    end
    
end
