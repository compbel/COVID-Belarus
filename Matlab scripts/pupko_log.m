function [L,C,labels] = pupko_log(AM,WM,traits,traitLeafs,Q,traitRoot,traitFreq)

nNodes = size(AM,1);
nTraits = length(traits);
L = zeros(nNodes,nTraits);
C = zeros(nNodes,nTraits);
root = find(sum(AM,1) == 0);
G = digraph(AM);
dfsorder = flip(dfsearch(G,root))';
outdeg = sum(AM,2);
indeg = sum(AM,1);
leafs = (find(outdeg == 0))';

% L(x,j) is the likelihood of the best reconstruction of this subtree on the condition that the father of 
%node x is assigned character state j. C(x,j) is the character state assigned to node x in this optimal 
%conditional reconstruction.


for i = 1:length(leafs)
    l = leafs(i);
    par = find(AM(:,l));
    P = expm(Q*WM(par,l));
    trait_leaf = find(strcmp(traits,traitLeafs{i}));
    C(l,:) = trait_leaf;
    L(l,:) = log(P(:,trait_leaf))';
end

for i = dfsorder
    if (outdeg(i) == 0)
        continue;
    end
    child = find(AM(i,:));
    if i ~= root
        par = find(AM(:,i));
        P = expm(Q*WM(par,i));
        for j = 1:nTraits
            [L(i,j),C(i,j)] = max(log(P(j,:)) + L(child(1),:) + L(child(2),:));
        end
    else
        L(i,:) = log(traitFreq) + L(child(1),:) + L(child(2),:);
    end
end

R = zeros(1,nNodes);
if ~isempty(traitRoot)
    R(root) = traitRoot;
else
    R(root) = find(L(root,:) == max(L(root,:)));
end
dfsorder = flip(dfsorder);
for i = dfsorder
    if i == root
        continue;
    end
    par = find(AM(:,i));
    R(i) = C(i,R(par));
end

labels = cell(1,nNodes);
for i = 1:length(labels)
    labels{i} = traits{R(i)};
end

