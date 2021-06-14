function selected_cells = selectCandidates(assemblie,min_prob)

if min_prob==1
    selected_cells = assemblie;
    return
end

numCandidates = length(assemblie);
prob_vector = 1:-(1-min_prob)/numCandidates:min_prob;
selected_cells = [];
permutation = randperm(numCandidates);

for i = 1:numCandidates
    coinflip = rand(1);
    if coinflip<prob_vector(i)
        selected_cells = [selected_cells assemblie(permutation(i))];
    end
end