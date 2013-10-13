d = 0;
% create new list of lists with same dimensions as domain
% used to keep track of which edges have been computed
done = ;

for i in range(len(domain)):
	for j in range(len(domain[i)):
			if !done[i][j]:
				k = find_edge(domain, i, j);
				d += tensor(edge_length(domain,i,j), angle(domain, i, k));

return d;

