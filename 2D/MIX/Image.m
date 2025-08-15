clear all

%% Define parameters
count = 1;
rep = 20;                  % Number of repetitions for averaging

% List of image sizes to test (square images: n x n)
vector = 10:10:150;         % Main case

const = size(vector, 2);    % Number of different sizes tested
RESU = {};                  % Cell array to store results for each size

% Probabilities for the p
p_list = 0.1:0.1:0.9;       % p values from 0.1 to 0.9
p_length = size(p_list, 2);

% Parameters for complexity calculation
r = 0.2;  % tolerance
m = 2;    % embedding dimension

% Preallocate results for repetitions × p values
res = zeros(rep, p_length);

%% Loop over different image sizes
for n = vector
    RES{count} = zeros(rep, 3);  % Placeholder
    X = zeros(n, n);             % Image matrix
    
    % Get adjacency matrix for an n x n grid with 8-neighbour connectivity
    A = Grid2D8nei(n, n);  
    B = triu(A);                 % Upper triangular part (avoids duplication)
    
    % Repeat experiment `rep` times for statistical reliability
    for kk = 1:rep
        t = 1;  % Counter
        
        % Loop over different p values
        for i = 1:p_length
            p = p_list(i);  % Current mixing probability
            
            % Binary mask: Z(i,j) = 1 with prob p, else 0
            Z = binornd(1, p * ones(n));
            
            % Noise matrix Y: uniform random values in [-√3, √3]
            Y = 2 * rand(n, n) * sqrt(3) - sqrt(3);
            
            % Construct MIX process image:
            % If Z(i,j)=0 → deterministic sinusoidal pattern
            % If Z(i,j)=1 → noise value from Y
            for ii = 1:n
                for jj = 1:n
                    X(ii, jj) = (1 - Z(ii, jj)) * (sin(2 * pi * ii / 12) + sin(2 * pi * jj / 12)) ...
                                + Z(ii, jj) * Y(ii, jj);
                end
            end
             
            % Convert image to a 1D vector
            matcol = reshape(X, 1, []);
            
            % Compute SEG average using sparse graph structure B
            res1(kk, t) = SEG(matcol, sparse(B), m, r);
                    
            t = t + 1;
        end
        
        % Display progress
        {kk, t}
    end
       
    % Store all results for this image size
    RESU{count} = res1; 
    count = count + 1;
    disp(n)  % Show current image size in console
end    

%% Resize results so RES1{p} is size × repetition
for kk = 1:p_length
    for mm = 1:const
        RES1{kk}(:, mm) = RESU{mm}(:, kk);
    end    
end 

%% Save data for plot
save('m_2_2D_MIX_averaged.mat', 'RES1');
