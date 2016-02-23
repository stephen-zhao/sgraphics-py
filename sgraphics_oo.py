import math


################################################################################
################################################################################
################################################################################
# CLASS Vector: Represents vectors. Provides corresponding vector operations
#               and related methods.

class Vector:
  
    #debug
    def get_v(self):
        return self.__v

    
    ############################################################################
    # CONSTRUCTORS:
    
    # __init(*v) initializes self with a tuple of its component
    #       coordinates.
    def __init__(self, *v):
        self.__v = list(v)
        self.dim = len(v)


    # build_0(n) builds the R^n zero-vector
    @classmethod
    def build_0(cls, n):
        return cls(*(0,)*n)


    # build_rep(n, x) builds the R^n vector with all components as x
    @classmethod
    def build_rep(cls, n, x):
        return cls(*(x,)*n)


    # build_seq(a, b, by=1) builds the vector with components starting at a,
    #       ending at one before b, stepping up 'by' per component
    @classmethod
    def build_seq(cls, a, b, by=1):
        assert(by!=0)
        assert(a!=b)
        assert((a-b)*by < 0)
        return cls(*tuple(a+by*i for i in range(int((b-a)/by)+1)))


    # build_count(n, start=0, by=1) builds the vector in R^n with components
    #       starting at 'start', counting up 'by' per component
    @classmethod
    def build_count(cls, n, start=0, by=1):
        return cls(*tuple(start+by*i for i in range(n)))
    


    # build_ei(n, i) builds the e_i basis vector in R^n (zero-indexing)
    @classmethod
    def build_ei(cls, n, i):
        return cls(*(1 if j == i else 0 for j in range(n)))



    ############################################################################
    # GETTERS: Functions to get a coordinate (projection onto basis), the vector
    #          itself in list form, the norm of the vector, the distance to
    #          another vector, and the iterator function.
    

    # get(i) gets self's ith coordinate.
    def get(self, i):
        return self.__v[i]


    # get_list() gets self's vector as a list
    def get_list(self):
        return self.__v
    

    # get_norm() gets the norm of self.
    def get_norm(self):
        return math.sqrt(Vector.dot(self, self))


    # get_dist(u) gets the distance from self to u.
    def get_dist(self, u):
        v = Vector.sub(self, u)
        return math.sqrt(Vector.dot(v, v))


    # __iter__() returns the list holding the vector coordinates.
    def __iter__(self):
        return iter(self.__v)


    
    ############################################################################
    # SETTER OPERATIONS: Methods to operate and mutate self. Includes addition,
    #                    subtraction, scalar multiplcation & division,
    #                    projection, and perpendicular. These methods
    #                    operate in the same memory location and mutate the 
    #                    operand. They return self to allow chaining of
    #                    operations.

    # set_add(u) sets self to self + u, returns self.
    def set_add(self, u):
        assert(u.dim == self.dim)
        for i in range(self.dim):
            self.__v[i] += u.get(i)
        return self
    

    # set_sub(u) sets self to self - u, returns self.
    def set_sub(self, u):
        assert(u.dim == self.dim)
        for i in range(self.dim):
            self.__v[i] -= u.get(i)
        return self
    

    # set_smult(s) sets self to (scalar s) * self, returns self.
    def set_smult(self, s):
        for i in range(self.dim):
            self.__v[i] *= s
        return self


    # set_sdiv(s) sets self to self / (scalar s), returns self.
    def set_sdiv(self, s):
        assert(s != 0)
        for i in range(self.dim):
            self.__v[i] /= s
        return self


    # set_proj(a) sets self to the projection of self onto a, returns self.
    def set_proj(self, a):
        assert(a.dim == self.dim)
        assert(any(a))
        sm = Vector.dot(self, a)/Vector.dot(a, a)
        for i in range(self.dim):
            self.__v[i] = sm * a.get(i)
        return self


    # set_perp(a) sets self to the perpendicular of self onto a, returns self.
    def set_perp(self, a):
        assert(a.dim == self.dim)
        assert(any(a))
        sm = Vector.dot(self, a)/Vector.dot(a, a)
        self.set_sub(Vector.smult(sm, a))
        return self


    
    ############################################################################
    # OPERATIONS: Class functions that operate on real numbers and Vector 
    #             objects and return a reference to a new Vector / a real
    #             number. These include vector functions for addition, 
    #             subtraction, scalar multiplcation & division, projection, 
    #             and perpendicular, As well as scalar functions for dot
    #              product and norm/distance.

    # add(*vs) adds one or more vectors in vs together, returns the new sum
    #       vector.
    @classmethod
    def add(cls, *vs):
        assert(v.dim == vs[0].dim for v in vs)
        u = [0 for i in range(vs[0].dim)]
        for i in range(vs[0].dim):
            u[i] = sum(map(lambda v: v.get(i), vs))
        return Vector(*tuple(u))


    # sub(*vs) subtracts one or more vectors in vs[1:] from vs[0], returns the
    #       new difference vector.
    @classmethod
    def sub(cls, *vs):
        assert(v.dim == vs[0].dim for v in vs)
        u = [vs[0].get(i) for i in range(vs[0].dim)]
        for i in range(vs[0].dim):
            u[i] -= sum(map(lambda v: v.get(i), vs[1:]))
        return Vector(*tuple(u))


    # smult(s, v) multiplies v by the scalar s, returns the new scalar multiple
    #       vector.
    @classmethod
    def smult(cls, s, v):
        return Vector(*tuple(map(lambda x: x*s, v)))


    # sdiv(s, b) divides v by the scalar s, returns the new scalar scalar
    #       multiple vector.
    @classmethod
    def sdiv(cls, s, v):
        assert(s != 0)
        return Vector(*tuple(map(lambda x: x/s, v)))


    # dot(v, u) takes the dot product of v and u, returns a real number.
    @classmethod
    def dot(cls, v, u):
        assert(v.dim == u.dim)
        return sum([v.get(i)*u.get(i) for i in range(v.dim)])


    # proj(v, a) projects v onto a, returns the new projected vector.
    @classmethod
    def proj(cls, v, a):
        u = Vector(*v)
        u.set_proj(a)
        return u


    # perp(v, a) takes the perpendicular of v onto a, returns the new
    #       perpendicular vector.
    @classmethod
    def perp(cls, v, a):
        u = Vector(*v)
        u.set_perp(a)
        return u


    # norm(v) takes the norm of v, returns a real number.
    @classmethod
    def norm(cls, v):
        return math.sqrt(Vector.dot(v, v))


    # dist(v, u) finds the distance from v to u, returns a real number.
    @classmethod
    def dist(cls, v, u):
        d = Vector.sub(v, u)
        return math.sqrt(Vector.dot(d, d))



    ###########################################################################
    # DISPLAY / FORMATTING: Methods/functions that format Vectors for output
    def toString(self, col=False):
        if not col:
            s = "[  "+"{:>7.6G}  "*self.dim+"]"
        else:
            s = "/ {:>7.6G} \\\n"+"| {:>7.6G} |\n"*(self.dim-2)+"\\ {:>7.6G} /"
        return s.format(*tuple(self.__v))

    def print(self, col=False):
        print(self.toString(col))






        


################################################################################
################################################################################
################################################################################
# CLASS Vector3D: Represents vectors in R^3. Provides corresponding vector
#                 operations and related methods.

# TODO: - create superclass "Vector" for general vectors in R^n
#       - refactor Vector3D to be subclass of Vector; many methods can be
#           migrated and generalized

class Vector3D(Vector):
    
    ############################################################################
    # CLASS VARIABLES:
    
    dim = 3
    

    ############################################################################
    # CONSTRUCTORS:
    
    # __init(*v) initializes self with 3 dimensions and a tuple of its
    #       coordinates.
    def __init__(self, *v):
        self.__v = list(v)
        self.dim = 3


    # build_0() builds the R^3 zero-vector
    @classmethod
    def build_0(cls):
        return cls(0,0,0)

    # build_1() builds a vector with all components at 1
    @classmethod
    def build_1(cls):
        return cls(1,1,1)


    # build_i() builds the i unit vector
    @classmethod
    def build_i(cls):
        return cls(1,0,0)


    # build_j() builds the j unit vector
    @classmethod
    def build_j(cls):
        return cls(0,1,0)


    # build_k() builds the k unit vector
    @classmethod
    def build_k(cls):
        return cls(0,0,1)


    # build_ei(i) builds the e_i basis vector in R^3
    @classmethod
    def build_ei(cls, i):
        return cls(*(1 if j == i else 0 for j in range(cls.dim)))


    # build_fromSpherical(r, theta, phi) builds a vector of r length with theta
    #       azimuth and phi inclination.
    @classmethod
    def build_fromSpherical(cls, r, theta, phi):
        return cls(r*math.sin(phi)*math.cos(theta),
                   r*math.sin(phi)*math.sin(theta),
                   r*math.cos(phi))


    # build_fromEulerAngles(r, alpha, beta, gamma) builds a vector of r length
    #       with Euler angles alpha, beta, gamma

    # TODO more
    


    ############################################################################
    # GETTERS: Functions to get x,y,z coordinates, the norm of the vector, the
    #          distance to another vector and the iterator function.

    # get_x() gets self's x-coordinate.
    def get_x(self):
        return self.__v[0]


    # get_y() gets self's y-coordinate.
    def get_y(self):
        return self.__v[1]


    # get_z() gets self's z-coordinate.
    def get_z(self):
        return self.__v[2]
    

    # get(i) gets self's ith coordinate.
    def get(self, i):
        return self.__v[i]


    # get_norm() gets the norm of self.
    def get_norm(self):
        return math.sqrt(Vector3D.dot(self, self))


    # get_dist(u) gets the distance from self to u.
    def get_dist(self, u):
        return math.sqrt(Vector3D.dot(self, u))


    # __iter__() returns the list holding the vector coordinates.
    def __iter__(self):
        return iter(self.__v)


    ############################################################################
    # SETTER OPERATIONS: Methods to operate and mutate self. Includes addition,
    #                    subtraction, scalar multiplcation & division, cross
    #                    product, projection, and perpendicular. These methods
    #                    operate in the same memory location and mutate the 
    #                    operand. They return self to allow chaining of
    #                    operations.

    # set_add(u) sets self to self + u, returns self.
    def set_add(self, u):
        assert(u.dim == self.dim)
        for i in range(self.dim):
            self.__v[i] += u.get(i)
        return self
    

    # set_sub(u) sets self to self - u, returns self.
    def set_sub(self, u):
        assert(u.dim == self.dim)
        for i in range(self.dim):
            self.__v[i] -= u.get(i)
        return self
    

    # set_smult(s) sets self to (scalar s) * self, returns self.
    def set_smult(self, s):
        for i in range(self.dim):
            self.__v[i] *= s
        return self


    # set_sdiv(s) sets self to self / (scalar s), returns self.
    def set_sdiv(self, s):
        assert(s != 0)
        for i in range(self.dim):
            self.__v[i] /= s
        return self


    # set_cross(u) sets self to self cross multiplied with u, returns self.
    def set_cross(self, u):
        assert(u.dim == self.dim == cls.dim)
        temp = self.__v[:]
        self.__v[0] = temp[1]*u.get_z() - temp[2]*u.get_y()
        self.__v[1] = temp[2]*u.get_x() - temp[0]*u.get_z()
        self.__v[2] = temp[0]*u.get_y() - temp[1]*u.get_x()
        return self


    # set_proj(a) sets self to the projection of self onto a, returns self.
    def set_proj(self, a):
        assert(a.dim == self.dim)
        assert(any(a))
        sm = Vector3D.dot(self, a)/Vector3D.dot(a, a)
        for i in range(self.dim):
            self.__v[i] = sm * a.get(i)
        return self


    # set_perp(a) sets self to the perpendicular of self onto a, returns self.
    def set_perp(self, a):
        assert(a.dim == self.dim)
        assert(any(a))
        sm = Vector3D.dot(self, a)/Vector3D.dot(a, a)
        self.set_sub(Vector3D.smult(sm, a))
        return self
        
    
    ############################################################################
    # OPERATIONS: Class functions that operate on real numbers and Vector3D 
    #             objects and return a reference to a new Vector3D / a real
    #             number. These include vector functions for addition, 
    #             subtraction, scalar multiplcation & division, cross product,
    #             projection, and perpendicular, As well as scalar functions
    #             for dot product and norm/distance.

    # add(*vs) adds one or more vectors in vs together, returns the new sum
    #       vector.
    @classmethod
    def add(cls, *vs):
        assert(v.dim == cls.dim for v in vs)
        u = [0, 0, 0]
        for i in range(cls.dim):
            u[i] = sum(map(lambda v: v.get(i), vs))
        return Vector3D(*tuple(u))


    # sub(*vs) subtracts one or more vectors in vs[1:] from vs[0], returns the
    #       new difference vector.
    @classmethod
    def sub(cls, *vs):
        assert(v.dim == cls.dim for v in vs)
        u = [vs[0].get(i) for i in range(cls.dim)]
        for i in range(cls.dim):
            u[i] -= sum(map(lambda v: v.get(i), vs[1:]))
        return Vector3D(*tuple(u))


    # smult(s, v) multiplies v by the scalar s, returns the new scalar multiple
    #       vector.
    @classmethod
    def smult(cls, s, v):
        return Vector3D(*tuple(map(lambda x: x*s, v)))


    # sdiv(s, b) divides v by the scalar s, returns the new scalar scalar
    #       multiple vector.
    @classmethod
    def sdiv(cls, s, v):
        assert(s != 0)
        return Vector3D(*tuple(map(lambda x: x/s, v)))


    # dot(v, u) takes the dot product of v and u, returns a real number.
    @classmethod
    def dot(cls, v, u):
        assert(v.dim == u.dim == cls.dim)
        return sum([v.get(i)*u.get(i) for i in range(cls.dim)])


    # cross(*vs) takes the sequential cross product of all vectors in vs,
    #       returns the new resulting product vector.
    @classmethod
    def cross(cls, *vs):
        u = Vector3D(vs[0])
        for v in vs[1:]:
            u.set_cross(v)
        return u


    # proj(v, a) projects v onto a, returns the new projected vector.
    @classmethod
    def proj(cls, v, a):
        u = Vector3D(*v)
        u.set_proj(a)
        return u


    # perp(v, a) takes the perpendicular of v onto a, returns the new
    #       perpendicular vector.
    @classmethod
    def perp(cls, v, a):
        u = Vector3D(*v)
        u.set_perp(a)
        return u


    # norm(v) takes the norm of v, returns a real number.
    @classmethod
    def norm(cls, v):
        return v.get_norm()


    # dist(v, u) finds the distance from v to u, returns a real number.
    @classmethod
    def dist(cls, v, u):
        return math.sqrt(Vector3D.dot(v, u))

        
    ###########################################################################
    # DISPLAY / FORMATTING: Methods/functions that format Vector3Ds for output
    def toString(self, col=False):
        if not col:
            return "[  {:>9.8G}  {:>9.8G}  {:>9.8G}  ]".format(*tuple(self.__v))
        else:
            return "/ {:>9.8G} \\\n| {:>9.8G} |\n\\ {:>9.8G} /".format(*tuple(self.__v))

    def print(self, col=False):
        print(self.toString(col))





# CLASS Matrix: Represents matrices of arbitrary size.
# TODO: - entire thing
class Matrix:
     
    ############################################################################
    # CLASS VARIABLES:

    # TODO: add class vars

    ############################################################################
    # CONSTRUCTORS:
    def __init__(self, nRows, nCols, *vs):
        self.rows = rows
        self.cols = cols
        self.M = [v for v in vs]
        
        
        return something
        # TODO


    def __iter__(self):
        return iter(self.M)
    
    def set_add(self, M):
        return self


    ###########################################################################
    # DISPLAY / FORMATTING: Methods/functions that format Matrices for output
    def toString(self):
        "/ {:>9.8G} \\\n| {:>9.8G} |\n\\ {:>9.8G} /".format(*tuple(self.__v))

    def print(self):
        print(self.toString())



#################### TESTS #####################################################

a = Vector3D(0,1,2)
b = Vector3D(3,4,5)
c = Vector3D.add(a,b)
n = Vector3D(0, 0, 1)

c.print()
c.print(True) #-> [3, 5, 7]
c.set_add(a)
c.print(True) #-> [3, 6, 9]
c.set_add(a).set_add(a)
c.print(True) #-> [3, 8, 13]
c.set_sub(Vector3D(3,7,11))
c.print(True) #-> [0, 1, 2]
a.set_add(c).set_add(c).set_add(c).set_sdiv(4)
a.print(True) #-> [0,1,2]
d = Vector3D.sub(b,a)
d.print(True) #-> [3, 3, 3]
