function [opinion, consensusTime, nPositiveAgent, deltaOpinion] = ...
 vectorized_simulate(meetingA, meetingB, opinion, mu, u, selfOpinionFactor, ...
 otherOpinionFactor, opinionBias, normalAgentFlag, nPositiveAgent, ...
 consensusFraction)
% "Parallel" simulation of our "employed" agent-model for a fixed size of agents
% and precalculated meeting times.
% Performing several simulations in parallel has the great advantage that one
% can simulate a single time step of the simulation using vectorization
% technique.
% m: [1]: number of simulations done in parallel
% t: [1]: number of time steps
% n: [1]: number of agents
%
% INPUT
% meetingA, meetingB: [m t]: agent meetingA(i, j) and meetingB(i, j) meet in
% simulation i at time j
% opinion: [m n]: opinion(i, j) is the initial opinion of agent j in 
%  simulation i
% mu: [m 1]: mu(i) is the convergence factor in simulation i
% u: [m 1]: u(i) is the uncertainty level in simulation i
% selfOpinionFactor, otherOpinionFactor, opinionBias: [m n]: opinion of agent j
%  in simulation i when meeting other agent is selfOpinionFactor(i, j) *
%  selfOpinion + otherOpinionFactor(i, j) * otherOpinion + opinionBias(i, j)
% normalAgentFlag: [m n]: normalAgentFlag(i, j) is 1 iff agent j in simulation i
%  is a "normal" agent i.e. not an emplyed agent
% nPositiveAgent: [m 1] number of agents with positive opinion (at initial time)
% consensusFraction: [1] fraction of (normal) agents with positive opinion s.t.
%  consensus is reached 
%
% OUTPUT
% opinion: [m n]: opinion(i,j) is the opinion of agent j in simulation i after
%  performing the simulation
% consensusTime: [m 1]: consensusTime(i) gives the number of time steps T taken
%  until in simulation i the fraction of agents with positive opinion is greater
%  or equal consensusFraction, i.e., 1 <= T <= t, if no consensus is reached
%  then T = t+1
% nPositiveAgent: [m 1]: number of agents with positive opinion after performing
%  the simulation
% deltaOpinion: [m 1]: deltaOpinion(i) gives the sum of the absolute values of
%  changes of agents' opinion in simulation i

m = size(meetingA, 1);
t = size(meetingA, 2);
n = size(opinion, 2);

% minOpinion, maxOpinion: [m 1]: minOpinion(i) resp. maxOpinion(i) define
%  the maximal and minimal opinion an agent can hold in simulation i
minOpinion = -ones(m, 1);
maxOpinion = ones(m, 1);
deltaOpinion = zeros(m, 1);
consensusTime = repmat(t+1, m, 1);

% calculate one-dimensional indices into [m n]-matrix,
% i.e. opinion(i, meetingA(i, j)) == opinion(m*(meetingA(i, j)-1)+1)
meetingA = m*(meetingA-1)+repmat((1:m)', 1, t);
meetingB = m*(meetingB-1)+repmat((1:m)', 1, t);

for k = 1:t
    oldOpinionA = opinion(meetingA(:, k));
    oldOpinionB = opinion(meetingB(:, k));
    
    % calculate opinion held of the interacting agents at time step k
    opinionA = selfOpinionFactor(meetingA(:, k)) .* oldOpinionA ...
             + otherOpinionFactor(meetingA(:, k)) .* oldOpinionB ...
             + opinionBias(meetingA(:,k));
    opinionB = selfOpinionFactor(meetingB(:, k)) .* oldOpinionB ...
             + otherOpinionFactor(meetingB(:, k)) .* oldOpinionA ...
             + opinionBias(meetingB(:,k));
    opinionA = min(opinionA, maxOpinion);
    opinionA = max(opinionA, minOpinion);
    opinionB = min(opinionB, maxOpinion);
    opinionB = max(opinionB, minOpinion);
    
    % normalAgentA, normalAgentB: [m 1]: 0-1-vector. An entry is 1 iff the
    % corresponding agent is a normal agent.
    normalAgentA = normalAgentFlag(meetingA(:, k));
    normalAgentB = normalAgentFlag(meetingB(:, k));

    diffOpinion = opinionA-opinionB;
    
    % changeOpinion is primarily (mu .* diffOpinion) but has entries equal to 0
    % where the absolute value of diffOpinion is smaller than the uncertainty
    % level (plus a small epsilon for comparing floating point numbers)
    changeOpinion = mu .* (abs(diffOpinion)<u+1e-9) .* diffOpinion;

    % changeOpinionA, changeOpinionB describe the effective change of agents'
    % opinion at time step k, i.e. employed agents don't change their opinion
    changeOpinionA = normalAgentA .* changeOpinion;
    changeOpinionB = normalAgentB .* changeOpinion;

    newOpinionA = oldOpinionA - changeOpinionA;
    newOpinionB = oldOpinionB + changeOpinionB;

    % accumulate opinion change
    deltaOpinion = deltaOpinion + abs(changeOpinionA) + abs(changeOpinionB);
    
    % update number of agents with positive opinion
    nPositiveAgent = nPositiveAgent-(oldOpinionA>0)+(newOpinionA>0) ...
                                   -(oldOpinionB>0)+(newOpinionB>0);
    % update consensusTime, i.e. set it if consensus is reached for the first
    % time
    consensusTime = min(consensusTime, ...
                        (nPositiveAgent<consensusFraction*n)*(t+1)+k);
    
    opinion(meetingA(:, k)) = newOpinionA;
    opinion(meetingB(:, k)) = newOpinionB;
end

end % vectorized_simulate(...)
