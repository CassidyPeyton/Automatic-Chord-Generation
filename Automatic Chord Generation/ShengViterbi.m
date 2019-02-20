function [T1,T2,path] = ShengViterbi(transMat, loglikeMat, initProb)
% Implementation of the Viterbi algorithm to find the path of states that has the
% highest posterior probability in a finite-state hidden Markov model
% (HMM).
%
% Input
%   - transMat      : the transition probability matrix, P(S_n | S_n-1).
%                       Rows correspond to the starting states; columns
%                       corresond to the ending states. Each row sums to 1.
%                       (size: nState * nState)
%   - loglikeMat    : the log-likelihood matrix of each state to explain
%                       each observation, log( P(O_n | S_n) ) Rows
%                       correspond to states; columns correspond to
%                       observations. (size: nState * nObserv)
%   - initProb      : the initial probability of all states, P(S_0). A
%                       column vector with nState elements.
%
% Output
%   - path          : the best path of states (S_1, ..., S_N) that gives
%                       the maximal posterior probability, P(S_1, ..., S_N
%                       | O_1, ... O_N). A column vector with nObserv elements.
%
% Author: XXX
% Created: XXX
% Last modified: XXX

if nargin<3     % use a uniform distribution if the initial probability is not given
    initProb = ones(size(transMat,1),1)/size(transMat,1); 
end
if size(transMat,1) ~= size(loglikeMat, 1) || size(transMat,1) ~= length(initProb)
    error('The number of states is not consistent in the transition matrix, the likelihood matrix, and the initial probability!');
end

% Your implementation starts here...
% the Viterbi algorithm idea from wikipedia
% https://en.wikipedia.org/wiki/Viterbi_algorithm
[K,T] = size(loglikeMat);
T1 = zeros(K,T);
T2 = zeros(K,T);
S = 1:7;
initProb = log(initProb);
T1(:,1) = initProb;
transMat = log(transMat);
for i = 2:T
    for j = 1:K
         [M,I] = max(T1(:,i-1) + transMat(:,j));
         T1(j,i) = M + loglikeMat(j,i);    
         T2(j,i) = I;     
    end
end

z = zeros(1,T);
x = zeros(1,T);
for i = 1:T
      [~,I] = max(T1(:,i));
      z(i) = I;
      x(i) = S(z(i));
end

for i = T:-1:2
    z(i-1) = T2(z(i),i);
    x(i-1) = S(z(i-1));
end
path = x;

end