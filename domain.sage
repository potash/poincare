from sage.modules.free_module_element import FreeModuleElement;
from sage.rings.real_mpfi import RealIntervalFieldElement;
from sage.rings.complex_interval import ComplexIntervalFieldElement
from sage.plot.plot3d.index_face_set import IndexFaceSet;
from sage.matrix.matrix import Matrix;
#from sage.algebras.quaternion_algebra_element import QuaternionAlgebraElement_abstract;

# initialize everything to the specified precision;
def init(n):
    global R, C, U3, i, j, k, R3, M, I, pi, R_0, R_1, R_2, R_4, R3_0;
    R = RealIntervalField(n);
    C = R.complex_field();
    U3.<i,j,k> = QuaternionAlgebra(R, -1,-1);
    R3 = VectorSpace(R,3);
    M = MatrixSpace(C,2,2);

    # constants
    I = C.gen();
    pi = R.pi();
    R_0 = R(0);
    R_2 = R(2);
    R_4 = R(4);
    R_1 = R(1);
    R3_0 = R3([0,0,0]);

init(256);

def eq_rif(a, b):
    return not (a != b);

def eq_cif(a, b):
    return not (a.real() != b.real()) and not (a.imag() != b.imag());

def eq_vec(a, b):
    return eq_rif(a[0], b[0]) and eq_rif(a[1], b[1]) and eq_rif(a[2], b[2]);

def eq_mat(a, b):
    for i in range(a.nrows()):
        for j in range(a.ncols()):
            if not eq_cif(a[i,j], b[i,j]):
                return false;
    return true;

# check if two edges (i.e. unordered pair of elements of R3) are the same:
def eq_edge(a, b):
    return ( eq_vec(a[0], b[0]) and eq_vec(a[1], b[1]) ) or ( eq_vec(a[0], b[1]) and eq_vec(a[1], b[0]) );

# check if two sides are the same
def eq_side(a, b):
    if len(a) != len(b): return false;

    for e1 in a:
        found = false;
        for e2 in b:
            if eq_edge(e1, e2): found = true; break;
        if not found: return false;

    return true;

# return the center of a vector or list with interval components
def center(v):
    if isinstance(v, list) or isinstance(v, Matrix):
        return [center(x) for x in v];
        
    elif isinstance(v, FreeModuleElement):
        if len(v) == 2: return (v[0].center(), v[1].center());
        if len(v) == 3: return (RR(v[0].center()), RR(v[1].center()), RR(v[2].center()) );

    elif isinstance(v, RealIntervalFieldElement):
        return v.center();

    else: return v;

def invert_matrix(m):
    return R_1/m.det() * M ( [ [m[1,1], -m[0,1]], [-m[1,0], m[0,0]] ] );
    
def invert_quat(q):
    x = q;
    a = x[0]^2 + x[1]^2 + x[2]^2 + x[3]^2;
    b = x[0] - x[1]*i -x[2]*j - x[3]*k;
    return 1/a * b;
        
def complex_to_quat(z):
    return z.real() + z.imag()*i;
    
def vector_to_quat(z):
    return z[0] + z[1]*i + z[2]*j;
    
def quat_to_R3(q):
    return R3([q[0], q[1], q[2]]);

# transform a point of U^3 by an element of SL(2, C)
def transform(v, g):
    if isinstance(v, list):
        return [transform(x, g) for x in v];
    else:
        w = vector_to_quat(v);
        a = complex_to_quat(g[0,0])*w + complex_to_quat(g[0,1]);
        b = complex_to_quat(g[1,0])*w + complex_to_quat(g[1,1]);

        return quat_to_R3( a * invert_quat(b) );
      
# go from U^3 to D^3
def to_klein(v):
    if isinstance(v, list):
        return [to_klein(x) for x in v];

    return R3((2*v[0], 2*v[1], v[0]^2 + v[1]^2 + v[2]^2 - 1))/(1+ v[0]^2 + v[1]^2 + v[2]^2);

# go from D^3 to U^3
def to_upper(v):
    if isinstance(v, list):
        return [to_upper(x) for x in v];

    return R3([v[0], v[1], sqrt(1 - v[0]^2 - v[1]^2 - v[2]^2)])/(1 - v[2]);
    
# generate the list of words on num generators up to the specified length
# each word is a list of integers where n is the nth generator and -n is its inverse
def gen_words(num, length):
    
    words = [[i] for i in range(1, num + 1)];
    words.extend([[-i] for i in range(1, num + 1)]);

    last = words;    

    for i in range(length - 1):
        new = [];
        for word in last:
            last = word[len(word) - 1];
            
            for i in range(1, num + 1) :
                if last != -i : 
                    w = list(word); w.append(i);
                    new.append(w);
                if last != i : 
                    w = list(word); w.append(-i);
                    new.append(w);

        last = new;
        words.extend(new);

    return words;

# compute the matrix of a word encoded as above
def to_matrix(word, gens, inv=[]):
    m = M(1);
    
    if len(inv) != len(gens):
        inv = [];
        for gen in gens:
            inv.append(invert_matrix(gen));

    for i in range(0, len(word)):
        if abs(word[i]) > len(gens): print word; print len(gens); print len(inv);
        if word[i] > 0 : m = m*gens[word[i] - 1];
        if word[i] < 0 : m = m*inv[- word[i] - 1];
    
    return m;
    
# compute the matrix of a word encoded as above
def to_matrices(word, gens):
    inv = [];
    for gen in gens:
        inv.append(invert_matrix(gen));
    
    return [to_matrix(w, gens, inv) for w in word];
    
def conjugate(m, g):
    if isinstance(m, list):
        return [conjugate(n, g) for n in m];
    
    return g*m*invert_matrix(g);

# compute generators specified by l,d,r as in the paper
def x_gens(n):
    a = sqrt(C(X_GENS[n][0])); b = sqrt(C(X_GENS[n][1])); c = sqrt(C(X_GENS[n][2]));

    return [M( [ [a, 0], [0, 1/a] ] ),
            M( [ [b*(c+1/c)/2, b*(c-1/c)/2], [(c-1/c)/(2*b), (c+1/c)/(2*b)] ] ) ];

def x_rels(n):
    return X_RELS[n];

# checks if two words are inverses of each other
# if gens is specified then will check numerically
def are_inverses(word1, word2, gens=false):

    exact = true;
    if not len(word1) == len(word2): exact = false;
    else:
        for i in range(len(word1)):
            if word1[i] != -word2[len(word2) - 1 - i]: exact = false; break;
    if exact or gens == false: return exact;

    return eq_mat(to_matrix(word1, gens)*to_matrix(word2, gens), M(1));
            
# add extra edges to keep side connected
def fix_side(side):
    i = 0;
    while i < len(side) :
        if not eq_vec(side[i][1], side[(i + 1) % len(side)][0]) : 
            side.insert(i + 1, [side[i][1], side[(i + 1) % len(side)][0]]); i = i + 1;
        i = i + 1;
    return side;
            
# order the edges of a side so they are connected (assumes each vertex appears exactly twice)
def order_side(side):
    if side == [] : return side;

    for i in range(len(side) - 1) :
        j = i + 1;
        while j < len(side) and not eq_vec(side[i][1], side[i+1][0]):
            if eq_vec(side[j][0], side[i][1]):
                side.insert(i+1, side.pop(j));
                break;
            if eq_vec(side[j][1], side[i][1]):
                edge = side.pop(j); edge.reverse();
                side.insert(i+1, edge);
                break;
            j = j + 1;

    if eq_vec(side[len(side) - 1][0], side[0][0]):
        edge = side.pop(); edge.reverse();
        side.append(edge);

    return side;

# removes duplicates and single vertex of elements in the list
def clean_side(side):

    for edge in side:
        if eq_vec(edge[0], edge[1]): side.remove(edge);

    i = 0;
    while i < len(side) - 1:
        j = i + 1;
        while j < len(side):
            if eq_edge(side[i], side[j]): 
                side.pop(j); j = j - 1;
            j = j + 1;
        i = i + 1;
    return side;

# get the distinct vertices in a list
def clean_edge(edge):
    for i in range(len(edge) - 1):
        j = i + 1;
        while j < len(edge):
            if eq_vec(edge[i], edge[j]): 
                edge.pop(j); j = j - 1;
            j = j + 1;
    return edge;

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

    if len(edges) != 0: print 'check_edges: ' + `len(edges)` + ' unpaired edges';
    return len(edges) == 0;

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

# check each side of the domain
def check_domain(domain):
    for i in range(len(domain)) :
        if check_side (domain[i]) == false :
            print 'check_domain: side ' + `i` + ' broken';
            return false;

    return check_edges(domain) and check_vertices(domain);
 

    
# dirichlet domain (in the disc model) is represented by a list of sides; each side is a list of vertices.

# initialize the domain to the unit box
def init_domain():
    v = [ 
        R3([1,1,-1]), R3([1,-1,-1]), R3([-1,-1,-1]), R3([-1,1,-1]),
        R3([1,1,1]), R3([1,-1,1]), R3([-1,-1,1]), R3([-1,1,1]) 
    ];

    domain = [ 
        [[v[0], v[1]], [v[1], v[2]], [v[2], v[3]], [v[3], v[0]]],    # back
        [[v[4], v[5]], [v[5], v[6]], [v[6], v[7]], [v[7], v[4]]],    # front
        [[v[1], v[5]], [v[5], v[6]], [v[6], v[2]], [v[2], v[1]]],    # bottom
        [[v[0], v[4]], [v[4], v[7]], [v[7], v[3]], [v[3], v[0]]],    # top
        [[v[4], v[5]], [v[5], v[1]], [v[1], v[0]], [v[0], v[4]]],    # right
        [[v[7], v[6]], [v[6], v[2]], [v[2], v[3]], [v[3], v[7]]]     # left
    ];

    return domain;

# returns distance from vertex v to the plane P = {x | n.x = -p} (Hessian normal form)
def distance(v, n, p):
    return n*v + p;
    


# intersect the edge by the plane defined by the orbit; returns [new edge, intersection data]
def intersect_edge(edge, orbit):
    dist = [distance(edge[0], orbit, sqrt(1-orbit*orbit) - 1), distance(edge[1], orbit, sqrt(1-orbit*orbit) - 1)]
    
    data = [false, eq_rif(dist[0], R_0), eq_rif(dist[1], R_0)];    # [new vertex?, e0 on plane?, e1 on plane?]
    
    # if both are above, discard the edge
    if dist[0] > 0 and dist[1] > 0:
        return [];

    # if the first is above and second is below, replace the first with the intersection point
    elif dist[1] < 0 and dist[0] > 0:
        data[0] = true; data[1] = true;
        t = dist[1] / (orbit*(edge[1]-edge[0]));
        return [ [t*edge[0] + (1-t)*edge[1], edge[1]], data];

    # if the second is above and first is below, replace the second with the intersection point
    elif dist[0] < 0 and dist[1] > 0:
        data[0] = true; data[2] = true;
        t = dist[1] / (orbit*(edge[1]-edge[0]));
        return [ [edge[0], t*edge[0] + (1-t)*edge[1]], data];

    elif dist[1] > 0 and eq_rif(dist[0], 0):
        data[1] = true;
        return [[edge[0], edge[0]], data];

    elif dist[0] > 0 and eq_rif(dist[1], 0):
        data[2] = true;
        return [[edge[1], edge[1]], data];

    # otherwise, keep the edge
    return [edge, data];

# return the result of intersecting the domain with the plane defined by the orbit
# words is a list of words corresponding to each side, word is the word corresponding to the orbit
# these will be updated to account for any sides or edges that are culled or created
def intersect_domain(domain, g, words=false, word=false):
    orbit = g;
    
    bad = false;    # is one of the sides contained in the plane?
    f1 = [];        # potential side created by plane

    d = [];         # the new domain

    culled = [];    # original index of sides which have been culled

    for side in domain:
        f = [];     # the new side
        e1 = [];    # potential edge created by the plane

        for edge in side:
            e = intersect_edge(edge, orbit);
            if e != []: 
                # make sure the edge is not a single vertex before adding it
                if not eq_vec(e[0][0], e[0][1]):
                    f.append(e[0]);
               
                if e[1][1]: e1.append(e[0][0]);
                if e[1][2]: e1.append(e[0][1]);

        # if the intersection left more than one edge, add it
        if len(f) > 1 : d.append(fix_side(f));
        elif words != false: culled.append(domain.index(side)); 

        if e1 != []:
            e1 = clean_edge(e1);
            if len(e1) > 2: bad = true;
            elif len(e1) == 2: f1.append(e1);

    for i in range(len(culled)):
        words.pop(culled[i] - i);
    
    # TODO: verify the assumption that a new plane is created iff. 
    # no side is contained in the intersection plane and len(f1) >= 3
    if not bad:
       f1 = order_side(clean_side(f1));
       if len(f1) >= 3: # TODO: can len(f1) == 2?
            d.append(f1);
            if words != false: words.append(word);

    return d;

# generate the domain defined by words in a, b up to the specified length
def gen_domain(gens, length):
    all_words = gen_words(len(gens), length);

    domain = init_domain();
    words = [ [] for side in domain ];    # the words associated to the sides of the domain

    inv = [];
    for gen in gens:
        inv.append(invert_matrix(gen));

    x = R3([R_0, R_0, R_1]);

    for i in range(len(all_words)):
        word = all_words[i];
        g = to_matrix(word, gens, inv)
        orbit = to_klein(transform(x, g));
        if not eq_vec(orbit, R3([R_0, R_0, R_0])):    # ignore identity
            domain = intersect_domain(domain, orbit, words, word);

    return domain, words;
    
# return a graphics plot of the domain
def get_plot(domain, solid=true, color='blue', highlight=false):
    if solid:
        vertices = [ [ center(edge[0]) for edge in side] for side in domain ];
        plot = IndexFaceSet(vertices, color=color);

        if highlight != false:
            vlast = [ [ center(edge[0]) for edge in domain[len(domain) - 1] ]];
            plot += IndexFaceSet(vlast, color=highlight);

        return plot;
    
    else:
        plot = line3d([(0,0,0), (0,0,0)]);
        for side in domain:
            for edge in side:
                plot += line3d( [center(edge[0]), center(edge[1])], color=color);

        if highlight != false:
            for edge in domain[len(domain) - 1]:
                plot += line3d( [center(edge[0]), center(edge[1])], color=highlight);

    return plot;
