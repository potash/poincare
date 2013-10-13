# attempts to generate dirichlet domain
# init is the initial word lenghth to use and max is the maximum number of subsequent length two runs
# the defaults work for all of the exceptional manifolds
def get_domain(gens, rels=[], init=5, max=4):

    print 'Determining automatic structure...',; sys.stdout.flush();
    gap.load_package('kbmag');
    G = gap_group(len(gens), rels);
    rws = G.KBMAGRewritingSystem();
    
    o = rws.OptionsRecordOfKBMAGRewritingSystem();r
    #o.maxstoredlen = gap([10,10]);
    
    auto = rws.AutomaticStructure(gap(true)); # large = true
    print auto; 
    if not auto: return None;

    print 'Estimating domain...',; sys.stdout.flush();
    domain, words = gen_domain(gens, init);
    print `len(domain)` + ' sides';
    
    ngens = gens; n = 0;
    
    while not check_poincare(domain, gap_word(words, G), rws) and n < max:
        print 'Estimating domain...',; sys.stdout.flush();
    
        words = clean_words(words, ngens);     
        ngens = to_matrices(words, ngens);

        domain, nwords = gen_domain(ngens, 2);
        
        words = expand_word(nwords, words);
        
        print `len(domain)` + ' sides';
        n = n + 1;

    return domain;
    
def check_poincare(domain, words, rws):
    print 'Checking hypotheses...',; sys.stdout.flush();
    
    # (estimate) domain is a polyhedron 
    if not check_polyhedron(domain): print false; return false;
    
    # three edges at every vertex
    if not check_vertices(domain): print false; return false;
    

    """    
    # (estimate) paired sides match
    for pair in pairs:
        side0 = to_upper(domain[pair[0]]);
        side1 = to_upper(domain[pair[1]]);
        side1 = transform( side1, to_matrix(words[pair[0]], gens));
        if not eq_side(side0, side1): print false; return false;
    """
    
    indices = get_indices(domain);
    
    # side pairing
    pairs = get_pairs(words, rws);
    if pairs == false: print false; return false;
    
    # cycle length 3
    for pair in pairs:
        g = words[pair[0]];
        for e in indices[pair[0]]:
            h = words[e];
            found = false;
            for e2 in indices[pair[1]]:
                h2 = words[e2];
                if eq_word(h, g*h2, rws): found = true; break;
            if not found: print false; return false;
    
    print true; return true;
    

def check_cycle(pairs, indices, words, rws):
    for pair in pairs:
        g = words[pair[0]];
        for e in indices[pair[0]]:
            h = words[e];
            found = false;
            for e2 in indices[pair[1]]:
                h2 = words[e2];
                if eq_word(h, g*h2, rws): found = true; break;
            if not found: print false; return false;
    return true;
    
# check if each vertex is common to exactly three sides
def check_vertices(domain):
    vert = [];
    for side in domain:
        for edge in side:
            vert.append(edge[0]);

    while len(vert) > 0:
        j = 1;
        count = 1;
        while j < len(vert):
            if eq_vec(vert[0], vert[j]): 
                vert.pop(j); count = count + 1; j = j - 1;
            j = j + 1;

        if count != 3: return false;

        vert.pop(0);

    return true;

# check if side is properly ordered, etc
def check_side(side):
    for i in range(len(side)) :
        if not eq_vec(side[i][1], side[ (i+1) % len(side)][0] ) : 
            return false;
    return true;

# check if each edge appears exactly twice
def check_edges(domain):
    edges = [];
    for side in domain:
        edges.extend(side);

    i = 0;
    while i < len(edges):
        for j in range(i + 1, len(edges)):
            if eq_edge(edges[i], edges[j]):
                edges.pop(j); edges.pop(i);
                i = i - 1; break;
        i = i + 1;

    return len(edges) == 0;

def check_polyhedron(domain):
    for i in range(len(domain)) :
        if check_side (domain[i]) == false :
            print 'check_domain: side ' + `i` + ' broken';
            return false;

    return check_edges(domain);
    
def get_pairs(words, rws):
    # find side pairing
    pairs = [];
    sides = range(0, len(words));   # indices to sides

    i = 0;
    while i < len(sides):
        j = i + 1;
        while j < len(sides):                   # search for its inverse
            if eq_word(words[sides[i]], words[sides[j]].Inverse(), rws):
                break;
            j = j + 1;
        
        if j == len(sides): # sides[i] is unpaired
            return false;
        else:
            pairs.append((sides[i], sides[j]));  # add pair to the list
            sides.pop(j); sides.pop(i);          # remove pair from the search
    return pairs;


# returns index of the side which shares the edge domain[k][l], or -1 if the edge is not shared
def find_edge(domain, k, l):
    edge = domain[k][l];
    
    for i in range(len(domain)):
        if i != k:
            for j in range(len(domain[i])):
                if eq_edge(edge, domain[i][j]): return i;
    return -1;
                  

# returns a list with the same structure as the domain, but each edge is a pair of indices of sides that create that edge    
def get_indices(domain):
    indices = [ [] for side in domain ];

    for i in range(len(domain)):
        for j in range(len(domain[i])):
            indices[i].append( find_edge(domain, i, j) );
                                 
    return indices;

def gap_group(n, rels):
    F = gap('FreeGroup(' + str(n) + ')');
    rels = gap(gap_word(rels, F));
    
    return F/rels;
    
def gap_word(word, group):
        
    if len(word) != 0 and isinstance(word[0], list):
        return [gap_word(w, group) for w in word];
    
    gens = group.GeneratorsOfGroup();
    w = group.Identity();
    
    for i in range(0, len(word)):
        if word[i] > 0: w = w*gens[word[i]];
        elif word[i] < 0: w = w*gens[-word[i]].Inverse();
        
    return w;

# removes inverses from a list of words (should use rewriting system? check speed)
def clean_words(words, gens=false):
    nwords = list(words);
    
    for i in range(len(nwords) - 1):
        j = i + 1;
        while j < len(nwords):
            if are_inverses(nwords[i], nwords[j], gens): 
                nwords.pop(j); j = j - 1;
            j = j + 1;
            
    return nwords;

# inverts the list representation of a word
def invert_word(word):
    w = []; n = len(word);
    for i in range(n):
        w.append(-word[n-i-1]);
    return w;

# expand a word whose letters are indices of old words
def expand_word(new, old):
    if len(new) != 0 and isinstance(new[0], list):
        return [expand_word(w, old) for w in new];
    w = [];
    
    for i in range(0, len(new)):
        if new[i] > 0: w.extend(old[new[i]-1]);
        elif new[i] < 0: w.extend(invert_word(old[-new[i]-1]));
        
    return w;
    
def eq_word(a, b, rws):
    return rws.ReducedForm((a*b.Inverse()).UnderlyingElement()).IsOne();


def read_rws(filename):
    print filename;
    return gap("ReadRWS(\"" + filename + "\")");

"""
kbmag_root
def automatic_structure
  kbprog, etc.
def read_rws(filename)
  rws = gap("ReadRWS(\" + filename + \")");
def write_rws(rws, filename)
  rws.WriteRWS(filename);
"""

"""
G = gap_group(n, rels);
rws = G.KBMAGRewritingSystem();

write_rws(rws, 'kbmag/x0');
kbprog('kbmag/x0');
gpmakefsa('kbmag/x0');
gpaxioms('kbmag/x0');

rws = read_rws(
"""

