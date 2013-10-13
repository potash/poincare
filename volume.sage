# precomputes the specified number of coefficients of the power series for the lobachevsky function
def init_lob(num):
    global lob_c;
    
    try: lob_c
    except NameError:
        lob_c = [0]*(num+1);
        for n in range(1, len(lob_c)):
            lob_c[n] = abs(R(bernoulli(2*n)))*R_2^(R_2*n+R_1)/((R_4*n)*R(factorial(2*n+1)));

init_lob(500);

# find the center of a side
def side_center(side):
    sum = R3([R_0,R_0,R_0]);
    for edge in side:
        sum = sum + edge[0];
    return R_1/R(len(side))*sum;

# iterate a sequence in reverse 
def reverse(sequence):
    i = len(sequence)
    while i > 0:
        i = i - 1
        yield sequence[i]

# compute lob using power series of pre-computed coefficients, Horner's method, and dynamic table lookup
def lob(t):
    global lob_c, lob_data;

    # power series is only defined on 0 < t < pi
    if eq_rif(t, pi) or eq_rif(t, R_0): return R_0;
    if not (t - pi < R_0): t = t - pi;
    
    # try to look up value
    for x in lob_data:
        if eq_rif(x[0], t):
            return x[1];
    
    sum = lob_c[len(lob_c) - 1];
    t2 = t*t;
    
    for n in range(1, len(lob_c)):
        sum = lob_c[len(lob_c)-1-n] + sum*t2;
        
    sum = (sum + R_1 - log(R_2*t))*t;
    lob_data.append((t, sum));
    
    return sum;

def cross_product(v1, v2):
    return R3([ v1[1]*v2[2]-v1[2]*v2[1] , v1[2]*v2[0]-v1[0]*v2[2] , v1[0]*v2[1]-v1[1]*v2[0] ]);


# inner product in klein model
def klein_product(a, b, x):

    k = x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
    l = (R_1 - k)*(R_1 - k);
    e11 = (R_1 - k + x[0]*x[0]) / l;
    e22 = (R_1 - k + x[1]*x[1]) / l;
    e33 = (R_1 - k + x[2]*x[2]) / l;
    e12 = x[0]*x[1]/l;
    e23 = x[1]*x[2]/l;
    e31 = x[2]*x[0]/l;

    return a[0]*b[0]*e11 + a[1]*b[1]*e22 + a[2]*b[2]*e33 +(a[0]*b[1]+a[1]*b[0])*e12 + (a[0]*b[2]+a[2]*b[0])*e31 + (a[1]*b[2]+a[2]*b[1])*e23;

# normalize a matrix
def normalize(m):
    return m / sqrt(m.det());

# returns the angle between tangent vectors a,b at x
def klein_angle(a,b,x):
    return acos(klein_product(a,b,x)/(sqrt(klein_product(a,a,x))*sqrt(klein_product(b,b,x))));

# find the hyperbolic volume of a tetrahedron with one ideal point
def volIT(a):
    return (lob((a[0]-a[1]-a[2]+pi)/R_2)+lob((-a[0]+a[1]-a[2]+pi)/R_2)+lob((-a[0]-a[1]+a[2]+pi)/R_2)-lob((a[0]+a[1]+a[2]+pi)/R_2)+lob((a[0]-a[4]-a[5]+pi)/R_2)+lob((-a[0]+a[4]-a[5]+pi)/R_2) +lob((-a[0]-a[4]+a[5]+pi)/R_2)-lob((a[0]+a[4]+a[5]+pi)/R_2)+lob((a[3]-a[1]-a[5]+pi)/R_2)+lob((-a[3]+a[1]-a[5]+pi)/R_2)+lob((-a[3]-a[1]+a[5]+pi)/R_2)-lob((a[3]+a[1]+a[5]+pi)/R_2) +lob((a[3]-a[4]-a[2]+pi)/R_2)+lob((-a[3]+a[4]-a[2]+pi)/R_2)+lob((-a[3]-a[4]+a[2]+pi)/R_2)-lob((a[3]+a[4]+a[2]+pi)/R_2))/R_2+lob((a[0]+a[3]+a[1]+a[4])/R_2)+lob((a[0]+a[3]+a[2]+a[5])/R_2) +lob((a[1]+a[4]+a[2]+a[5])/R_2);

# rotate v around the z-axis by theta
def rotate(v,theta):
    a=cos(theta)
    b=sin(theta)
    return R3([ a*v[0]+b*v[1], -b*v[0]+a*v[1], v[2]]);

# returns the angle on x-y plane
def get_angle(v):
    t = atan(v[1]/v[0]);
    if v[0] < R_0:
        t= t+pi;
    return t;

# rotate v[2], v[3] in x-y plane so they are in first and second quadrants, respectively
def reposition2(v):

    avg = (v[2]+v[3])/R_2;
    t = get_angle(avg);

    v[2] = rotate (v[2],t-pi/R_2);
    v[3] = rotate (v[3],t-pi/R_2);

    if get_angle(v[2]) > get_angle(v[3]):
        temp=v[2];
        v[2]=v[3];
        v[3]=temp;

    return v;


# transforms a list of 4 vertices
def transform_tetra(v, g):
    return [ transform(v[0], g), transform(v[1], g), transform(v[2], g), transform(v[3], g) ];

# reposition v so that v[0] and v[1] are on the z-axis with v[0] below v[1]
def reposition(v):
    # put v[0] at (0, 0, 1)
    if not (eq_rif(v[0][0], R_0) and eq_rif(v[0][1], R_0)):
        A = M( [[R_1, -v[0][0] - v[0][1]*I ], [R_0, v[0][2]] ] );
        v = transform_tetra(v, A);
    
    # if v[1] is already on the z-axis
    if eq_rif(v[1][0], R_0) and eq_rif(v[1][1], R_0):
        # and v[1] is above v[0]
        if v[1][2] > v[0][2]:
            return v;
        # if v[1] is below, invert
        else:
            A = M([ [R_0, R_1], [-R_1, R_0] ]);
            return transform_tetra(v, A);

    x = sqrt(v[1][0]*v[1][0] + v[1][1]*v[1][1]);
    n = v[1]*v[1];

    A = M([ [v[1][0]/x - v[1][1]/x * I, (R_1-n+sqrt(R_4*x*x+(R_1-n)*(R_1-n)))/(R_2*x)],
          [-(v[1][0]/x-I*v[1][1]/x)*(R_1-n+sqrt(R_4*x*x+(R_1-n)*(R_1-n)))/(R_2*x), R_1] ]);

    A = A / sqrt(A.det());  # normalize
    return transform_tetra(v, A);

# get normal by getting center of the plane not including the origin
def last_normal(v,u,w):
    t = ((w*w)*(u*v)-(u*w)*(v*w))/((u*w)*(u*w)-(u*u)*(w*w));
    s = ((u*u)*(w*v)-(u*w)*(v*u))/((u*w)*(u*w)-(u*u)*(w*w));
    return v + t*u + s*w;

def volT(v):
    v = reposition(v);
    v = reposition2(v);
    
    v = to_klein(v);
    z = R3([R_0,R_0,R_1]);
    sum = R3([R_0,R_0,R_0]);

    n1 = cross_product(v[3]-v[0],v[1]-v[0]);
    n2 = cross_product(v[1]-v[0],v[2]-v[0]);
    n3 = cross_product(v[2]-v[0],v[3]-v[0]);
    n4 = last_normal(v[3],v[2]-v[3],v[1]-v[3]);
    n5 = last_normal(z,v[3]-z,v[2]-z);

    a = n4[0]*n4[0] + n4[1]*n4[1] + n4[2]*n4[2];
    p4 = R_1/a * n4;

    a = n5[0]*n5[0] + n5[1]*n5[1] + n5[2]*n5[2];
    p5 = R_1/a * n5;

    vol1 = volIT([klein_angle(-n1,n2,v[0]),
              klein_angle(-n1,n3,v[0]),
              klein_angle(-n2,n3,v[0]),
              klein_angle(-n3,-p5+v[2],v[2]),
              klein_angle(-n2,-p5+v[2],v[2]),
              klein_angle(-n1,-p5+v[3],v[3])]);


    vol2 = volIT([klein_angle(-n1,n2,v[1]),
              klein_angle(-n1,p4-v[1],v[1]),
              klein_angle(-n2,p4-v[1],v[1]),
              klein_angle(-p4+v[2],-p5+v[2],v[2]),
              klein_angle(-n2,-p5+v[2],v[2]),
              klein_angle(-n1,-p5+v[3],v[3])]);
    
    return vol1 - vol2;

# takes domain vertices in klein model (in order to find centers of sides)!
def volume(domain):
    global lob_data; lob_data = [];
    total = R_0;
    origin = R3([R_0,R_0,R_1]);
  
    for side in domain:
        center = to_upper(side_center(side));
        side = to_upper(side);
        for edge in side:
            vol = volT([origin, center, edge[0], edge[1]]);
            total = total + vol;

    return total;
