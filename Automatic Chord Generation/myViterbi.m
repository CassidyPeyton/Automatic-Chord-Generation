function path = myViterbi(transMat, loglikeMat, initProb)
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
initProb=initProb;
transMat = transMat;
[K , N] = size(loglikeMat);
T1=zeros(K,N);
T2=zeros(K,N);
for i=1:K
    T1(i,1)=initProb(i)+loglikeMat(i,1);
    T2(i,1)=0;
end
for i=2:N
    for j=1:K
        [M,I]=max(T1(:,i-1)+transMat(:,j)+loglikeMat(j,i));
        T1(j,i)=M;
        T2(j,i)=I;
    end
end
z=zeros(1,N);
X=zeros(1,N);
[M,I]=max(T1(:,N));
z(N)=I;

X(N)=(I-2);
for i=1:N-1
    z(N-i)=T2(z(N-i+1),N-i+1);
    X(N-i)=(z(N-i)-2);
end
path=X;
end