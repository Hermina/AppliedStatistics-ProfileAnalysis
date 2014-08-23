function PROFANA(X,alpha)
%PROFANA Profile Analysis of Multivariate Repeated Measures.
% Profile analysis is a special application of multivariate analysis of
% variance (MANOVA) in which several dependent variables are measured and
% they are all measured on the same scale. Is also used in research where
% subjects are measured repeatedly on the same dependent variable. In this
% case, profile analysis is an alternative to univariate repeated measures
% analysis of variance. The major question to be answered by profile
% analysis is whether profiles of groups differ on a set of measures. To
% apply profile analysis, all measures must have the same range of possible
% scores, with the same score value having the same meaning on all the 
% measures. The restriction on scaling of the measures is present because
% in two of the major tests of profile analysis (parallelism and flatness)
% the numbers that are actually tested are difference scores between
% dependent variables measured on adjacent occasions. Difference scores are
% called segments in profile analysis.
% There are three main questions to be answered on profile analysis:
% --1. Do different groups have parallel profiles?: test of parallelism.
%    When profile analysis is used as a substitute for univariate repeated
%    measures ANOVA, the parallelism test is the test of interaction.
% --2. Whether or not groups produce parallel profiles, does one group, on
%    average, score higher on the collected set of measures than another?:
%    test of difference in levels. In regular ANOVA, this question is
%    answered by test of the groups hypothesis. In repeated-measures ANOVA,
%    it address as the between-subjects main effects quastion.
% --3. Do the dependent variables all elicit the same average response?:
%    test of flatness. This is relevant only if the profiles are parallel.
%    If not, the at least one of them is necessarily not flat. This test 
%    evaluates the within-subjects main effect hypothesis in repeated-
%    measures ANOVA.
%
% Syntax: function PROFANA(X,alpha) 
%      
% Inputs:
%      X - data matrix (Size of matrix must be n-by-(1+p); sample=column 1,
%          variables=column 2:p). 
%  alpha - significance level (default = 0.05).
%
% Outputs:
%        - Complete analysis of variance for difference levels test.
%        - Complete analysis of variance for parallelism test.
%        - Complete Hotelling's T-squared analysis for flatness test.
%        - Figure of profiles of dependent variables for the studied groups.
%
% Example: From the data Table 10.3 of Tabachnick and Fidell (1996, p. 448).
%        Three groups (independent variables) whose profiles are compared are
%        belly dancers (1), politicians (2), and administrators (3). The
%        five respondents in each of these occupational groups participate in
%        four leisure activities (dependent variables) and, during each, are
%        asked to rate their satisfaction on a 10-points scale. Dependent 
%        variables are read(1), dance(2), TV(3), and ski(4). Profile analysis
%        of difference in levels, parallelism, and flatness are tested with a
%        significance of 0.10.
%
%                                      Activity
%                 ------------------------------------------
%                  Group       1       2       3       4
%                 ------------------------------------------
%                    1         7      10       6       5
%                              8       9       5       7  
%                              5      10       5       8
%                              6      10       6       8
%                              7       8       7       9
%                    2         4       4       4       4
%                              6       4       5       3
%                              5       5       5       6 
%                              6       6       6       7
%                              4       5       6       5
%                    3         3       1       1       2
%                              5       3       1       5
%                              4       2       2       5
%                              7       1       2       4
%                              6       3       3       3
%                 ------------------------------------------
%                                       
% Data matrix must be:
%   X=[1 7 10 6 5;1 8 9 5 7;1 5 10 5 8;1 6 10 6 8;1 7 8 7 9;2 4 4 4 4;2 6 4 5 3;
%   2 5 5 5 6;2 6 6 6 7;2 4 5 6 5;3 3 1 1 2;3 5 3 1 5;3 4 2 2 5;3 7 1 2 4;3 6 3 3 3];
%
% Calling on Matlab the function: 
%   PROFANA(X,0.1)
%
% Answer is:
%
% The number of groups are: 3
% The number of dependent variables are: 4
%
% Difference in Levels Test.
% 
% Analysis of Variance Summary Table for Test of Levels Effect.
% ----------------------------------------------------------------------------
% SOV                        SS          df         MS          F        P
% ----------------------------------------------------------------------------
% Between Groups          172.900         2       86.450      44.145   0.0000
% Within Groups            23.500        12        1.958
% Total                   196.400        14
% ----------------------------------------------------------------------------
% The associated probability for the F test is smaller than 0.10
% So, the assumption that differences between means of the groups over the dependent
% variables are equal was met.
% The strength of association between groups and averaged dependent variables are 0.88
% 
% Parallelism Test.
% It is considering as a small sample problem (n < 25).
% 
% Multivariate Analysis of Variance Table.
% --------------------------------------------
% No. data    Samples     Variables       L
% --------------------------------------------
%    15          3            4        0.0763
% --------------------------------------------
% 
% ------------------------------------------------------------------------------
% Test                 Statistic     df1     df2         F       P    Conclusion
% ------------------------------------------------------------------------------
% Rao                    0.076         6      20       8.74   0.0001       S
% Pillai                 1.433         6      22       9.28   0.0000       S
% Lawley-Hotelling       5.428         6      18       8.14   0.0002       S
% Roy                    3.541       3.0    11.0      12.98   0.0006       S
% ------------------------------------------------------------------------------
% With a given significance of: 0.10
% According to the P-value, the sample mean vectors could be significant (S) or
% not significant (NS).
% The strength of association between groups and averaged dependent variables are 0.72
% Here, if significant, this leading to rejection of the hypothesis of parallelism.
% 
% Flatness Test.
% As n is less than 50: 15
% We use the F approximation.
% 
% -------------------------------------------------------------------------------------
% Sample-size    Variables       T2         F         df1           df2           P
% -------------------------------------------------------------------------------------
%      15            3         2.5825     8.6082        3            10        0.0040
% -------------------------------------------------------------------------------------
% The associated probability for the F test is smaller than 0.10
% There is significant deviation from flatness.
% The strength of association between groups and averaged dependent variables are 0.72
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, A. Castro-Perez
%             and K. Barba-Rojo
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
% Copyright (C) December 24, 2006.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez and K. Barba-Rojo. (2006).
%   PROFANA:Profile Analysis of Multivariate Repeated Measures. A MATLAB file.
%   [WWW document]. URL http://www.mathworks.com/matlabcentral/fileexchange/
%   loadFile.do?objectId=13483
%
%  Reference:
%  Tabachnick, B. G. and Fidell, L. S. (1996), Using Multivariate Statistics.
%           3rd. ed. NY: HarperCollins College Pub, p. 441-505.
%

if nargin < 2,
   alpha = 0.05;
end 

if (alpha <= 0 | alpha >= 1)
   fprintf('Warning: significance level must be between 0 and 1\n');
   return;
end;

if nargin < 1, 
   error('Requires at least one input argument.');
   return;
end;

g = max(X(:,1));
p = size(X(:,2:end),2);
fprintf('The number of groups are:%2i\n', g);
fprintf('The number of dependent variables are:%2i\n\n', p);

%%%Difference in Levels.
disp('Difference in Levels Test.')
MGL = [];M = [];n = [];
indice=X(:,1);
for i = 1:g
   Xe=find(indice==i);
   eval(['XG' num2str(i) '= X(Xe,2:end);']);
   eval(['mG' num2str(i) '= mean(XG' num2str(i) ',2);'])
   eval(['mGl' num2str(i) '= mean(XG' num2str(i) ',1);'])
   eval(['mmG' num2str(i) '= mean(mG' num2str(i) ');'])
   eval(['nG' num2str(i) '= length(XG' num2str(i) ') ;'])
   eval(['xm = mmG' num2str(i) ';'])
   eval(['xmgl = mGl' num2str(i) ';'])
   eval(['xn = nG' num2str(i) ';'])
   MGL = [MGL;xmgl];M = [M;xm];n = [n;xn];
end

S1 = [];
for i = 1:g
    eval(['S1' num2str(i) '= nG' num2str(i) '*(mmG' num2str(i) ' - mean(M)).^2;'])
    eval(['xS1 = S1' num2str(i) ';'])
    S1 = [S1;xS1];
end
S1
SSBG = p*sum(S1);  %sum of squares between groups
v1 = g-1;
MSBG = SSBG/v1;

S2 = [];
for i = 1:g
    eval(['S2' num2str(i) '= (mG' num2str(i) ' - mmG' num2str(i) ').^2;'])
    eval(['xS2 = S2' num2str(i) ';'])
    S2 = [S2;xS2];
end
S2
SSWG = p*sum(S2);  %sum of squares within groups
v2 = sum(n)-g;
MSWG = SSWG/v2;

SST = SSBG+SSWG;
v3 = v1+v2;
FL = MSBG/MSWG; 
P = 1 - fcdf(FL,v1,v2);  %probability associated to the F-statistic.  

disp(' ')
disp('Analysis of Variance Summary Table for Test of Levels Effect.')
fprintf('----------------------------------------------------------------------------\n');
disp('SOV                        SS          df         MS          F        P')
fprintf('----------------------------------------------------------------------------\n');
fprintf('Between Groups  %15.3f%10i%13.3f%12.3f%9.4f\n',SSBG,v1,MSBG,FL,P);
fprintf('Within Groups%18.3f%10i%13.3f\n',SSWG,v2,MSWG);
fprintf('Total%26.3f%10i\n',SST,v3);
fprintf('----------------------------------------------------------------------------\n');
if P >= alpha;
    fprintf('The associated probability for the F test is equal or larger than% 3.2f\n', alpha);
    disp('So, the assumption that differences between means of the groups over the dependent')
    disp('variables are equal was met.');
else
    fprintf('The associated probability for the F test is smaller than% 3.2f\n', alpha);
    disp('So, the assumption that differences between means of the groups over the dependent')
    disp('variables are equal was met.');
end
n2 = SSBG/SST;
fprintf('The strength of association between groups and averaged dependent variables are% 3.2f\n', n2);

%%%Parallelism.
disp(' ')
disp('Parallelism Test.')

Xi = X(:,2:p);
Xf = X(:,3:p+1);

X = [X(:,1) [Xi-Xf]];

MAOV1(X,alpha);
disp('Here, if significant, this leading to rejection of the hypothesis of parallelism.')

%%%Flatness.
disp(' ')
disp('Flatness Test.')

g = max(X(:,1));

N = [];
indice = X(:,1);
for i = 1:g
   Xe = find(indice==i);
   eval(['X' num2str(i) '= X(Xe,2);']);
   eval(['n' num2str(i) '= length(X' num2str(i) ') ;'])
   eval(['xn= n' num2str(i) ';'])
   N = [N,xn];
end

[f,c] = size(X);
X = X(:,2:c);

[n,p] = size(X);
r = 1;
r2 = N(1);

for k = 1:g
   eval(['M' num2str(k) '= mean(X(r:r2,:));']);
   if k < g
      r = r+N(k);
      r2 = r2+N(k+1);
   end
end
M = mean(X);

dT = [];
for k = 1:p
   eval(['dT  = [dT,(X(:,k) - mean(X(:,k)))];']);
end
T = dT'*dT;  %total sum of squares

r = 1;
r2 = N(1);
g = length(N);
for k = 1:g
   Md(k,:) = (mean(X(r:r2,:)) - mean(X));
   if k < g
      r = r+N(k);
      r2 = r2+N(k+1);
   end
end

H = [];  %between samples sum of squares
for k = 1:g
   h = N(k)*(Md(k,:)'*Md(k,:));
   if k == 1
      H = h;
   else
      H = H+h;
   end
end

E = T-H;  %within samples sum of squares

T2 = sum(N)*M*inv(E)*M';  %observed Hotelling's T-squared

LW = 1/(1+T2);

if n >= 50  %Chi-square approximation.   
    X2 = T2;
    v = p+1; %Degrees of freedom.
    P = 1-chi2cdf(X2,v); %Probability that null Ho: is true.
    fprintf('As n is equal or larger than 50:% i\n', sum(n));
    disp('We use the Chi-square approximation.')
    disp(' ')
    fprintf('----------------------------------------------------------------------------\n');
    disp(' Sample-size    Variables      T2          Chi-sqr.         df          P')
    fprintf('----------------------------------------------------------------------------\n');
    fprintf('%8.i%13.i%15.4f%14.4f%11.i%14.4f\n\n',n,p,T2,X2,v,P);
    fprintf('----------------------------------------------------------------------------\n');
    if P >= alpha;
        fprintf('The associated probability for the Chi-squared test is equal or larger than% 3.2f\n', alpha);
        disp('There is not significant deviation from flatness.');
    else
        fprintf('The associated probability for the Chi-squared test is smaller than% 3.2f\n', alpha);
        disp('There is significant deviation from flatness.');
    end 
    n2 = 1-LW;
    fprintf('The strength of association between groups and averaged dependent variables are% 3.2f\n', n2); 
else  %F approximation.
    F = T2*(sum(N)-g-p+1)/p;  %F approximation
    v1 = p;  %Numerator degrees of freedom.
    v2 = sum(N)-g-p+1;  %Denominator degrees of freedom.
    P = 1-fcdf(F,v1,v2);  %Probability that null Ho: is true.
    fprintf('As n is less than 50:% i\n', sum(n));
    disp('We use the F approximation.')
    disp(' ')
    fprintf('-------------------------------------------------------------------------------------\n');
    disp(' Sample-size    Variables       T2         F         df1           df2           P')
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('%8.i%13.i%15.4f%11.4f%9.i%14.i%14.4f\n',n,p,T2,F,v1,v2,P);
    fprintf('-------------------------------------------------------------------------------------\n');
    if P >= alpha;
        fprintf('The associated probability for the F test is equal or larger than% 3.2f\n', alpha);
        disp('There is not significant deviation from flatness.');
    else
        fprintf('The associated probability for the F test is smaller than% 3.2f\n', alpha);
        disp('There is significant deviation from flatness.');
    end
    n2 = 1-LW;
    fprintf('The strength of association between groups and averaged dependent variables are% 3.2f\n', n2); 
end

G = [];
for k = 1:g
    G = [G;k];
end

sc = [G MGL];

figure;
hold on;

lg = [];
for k = 1:g
    plot(1:p+1,sc(k,2:end),DCSYMB0(k),'LineStyle','--','Color',DCRGB0(k),'MarkerFaceColor',DCRGB0(k),'MarkerEdgeColor',DCRGB0(k));
    lg = [lg,['''Group ' num2str(k) ''',']];
    axis([0 p+2 min(min((sc(:,2:end))))-1 max(max((sc(:,2:end))))+1])
    %line('LineStyle','--','Color',DCRGB0(k));
end
lg(end)=' ';
eval(['legend(' lg ')']);

title('Profiles of Dependent Variables for the Studied Groups.');
xlabel('Dependent  Variables');
ylabel('Mean  of  Dependent  Variable');
set(gca,'XTick',1:1:p+1)

hold off;

return,