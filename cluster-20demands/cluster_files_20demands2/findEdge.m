function [edgeList, edgeIdxList, networkIncidenceMatrix, linkLength] = findEdge(adj, adjw)
% generate edgeList and edgeIdxList

% create edgeList
N = size(adj, 1);
adj = spdiags(sparse(N, 1), 0, adj);
[rowidx, colidx, ~] = find(adj);
edgeList = [rowidx, colidx];
edgeList = sortrows(edgeList);

% create edgeIdxList
edgeAll = ones(N)-diag(ones(N, 1));
[rowidx, colidx] = find(edgeAll);
edgeAll = [rowidx, colidx];
edgeAll = sortrows(edgeAll);
[~, edgeIdxList] = ismember(edgeList, edgeAll, 'rows');

rowidx = [(1:length(edgeIdxList))', (1:length(edgeIdxList))']; 
rowidx = rowidx.'; 
rowidx = rowidx(:);
colidx = edgeList.';
colidx = colidx(:);
nodeNum = size(adj, 1);
networkIncidenceMatrix = full(sparse(rowidx, colidx, 1, ...
    length(edgeIdxList), nodeNum));

for i=1:size(edgeList, 1)
    networkIncidenceMatrix(i, edgeList(i, 2)) = -1;
end

linkLength = adjw(sub2ind(size(adjw), edgeList(:, 1), edgeList(:, 2)));
end