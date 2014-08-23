function res = dcsymb0(n)
% - utility program for selecting plot symbols

while (n>10)
   n = n - 10;
end

% symbols to choose from:
symbols = {'+','*','d','.','x','o','p','s','^','v','<','h','>'};;

res = symbols{n};