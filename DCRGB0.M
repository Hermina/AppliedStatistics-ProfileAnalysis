function res = dcrgb0(n)

% - utility program for selecting plot colors

while (n>10)
   n = n - 10;
end

% build cell-array of colors to choose from:
colors{1}  = [0 0 1];          % blue
colors{2}  = [1 0 0];          % red
colors{3}  = [0 1 0];          % green
colors{4}  = [0.5 0 0.5];      % purple
colors{5}  = [0 0.5 0];        % dark green
colors{6}  = [0.28 0.73 0.94]; % sky
colors{7}  = [0 0 0];          % black
colors{8}  = [0 0 0.5];        % navy
colors{9}  = [1 0.5 0];        % orange
colors{10} = [1 0 1];          % magenta

res = colors{n};

