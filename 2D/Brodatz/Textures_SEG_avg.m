clear all
%% ------------------------------------------------------------------------
% 1. Read and scale textures
d = 'textures';
files = {'a_05.gif','b_15.gif','c_30.gif','d_36.gif','e_45.gif', ...
         'f_75.gif','g_93.gif','h_95.gif','i_102.gif'};
Nimg  = numel(files);

I = cell(1,Nimg);
for k = 1:Nimg
    I{k} = mat2gray(imread(fullfile(d,files{k})));   % 0‑1 double
end

%% ------------------------------------------------------------------------
% 2. Parameters for the SEG calculation
patchSz = 128;              % size of each square patch
m       = 2;                % embedding dimension
r_list  = 0.05:0.05:0.40;   % tolerance values

% adjacency for one 128×128 grid is reused for every patch
A   = Grid2D8nei(patchSz,patchSz);
B   = sparse(triu(A));

% result array: rows = image, cols = r values
res = zeros(Nimg,numel(r_list));

%% ------------------------------------------------------------------------
% 3. Main loop
for k = 1:Nimg
    img = I{k};
    [h,w] = size(img);

    % ensure the image is at least 640×640 so we can fit 5×5 patches
    assert(h >= patchSz*5 && w >= patchSz*5, ...
        'Image %s is too small for 25 non‑overlapping %dx%d patches.', ...
        files{k}, patchSz, patchSz);

    % coordinates of the 5×5 grid
    rowIdx = 1:patchSz:patchSz*5;
    colIdx = 1:patchSz:patchSz*5;

    for rIdx = 1:numel(r_list)
        r   = r_list(rIdx);
        val = zeros(25,1);      % holds SampEn for the 25 patches
        p   = 1;                % patch counter

        for rr = rowIdx
            for cc = colIdx
                patch      = img(rr:rr+patchSz-1, cc:cc+patchSz-1);
                x          = reshape(patch,1,[]);
                val(p)     = SEG(x,B,m,r);   % Sample Entropy on patch
                p = p + 1;
            end
        end

        res(k,rIdx) = mean(val);   % ***average over 25 patches***

        disp({k,rIdx})
    end
end


%% ------------------------------------------------------------------------
save('SEG_Brodatz_texture_m2_25patch.mat','res')
