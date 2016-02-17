import math

################################################################################
################################################################################
################################################################################
# CLASS Vector3D: Represents vectors in R^3. Provides corresponding vector
#                 operations and related methods.

# TODO: - create superclass "Vector" for general vectors in R^n
#       - refactor Vector3D to be subclass of Vector; many methods can be
#           migrated and generalized

class Vector3D:
    
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


    # build_fromSpherical(r, theta, phi) builds a vector of r length with theta
    #       azimuth and phi inclination.
    @classmethod
    def build_fromSpherical(cls, r, theta, phi):
        return cls(r*math.sin(phi)*math.cos(theta),
                   r*math.sin(phi)*math.sin(theta),
                   r*math.cos(phi))


    # build_fromEulerAngles(r, alpha, beta, gamma) builds a vector of r length
    #       with Euler angles alpha, beta, gamma

    # TODO
    


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
        return sqrt(Vector3D.dot(self, self))


    # get_dist(u) gets the distance from self to u.
    def get_dist(self, u):
        return sqrt(Vector3D.dot(self, u))


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


    # set_proj(n) sets self to the projection of self onto n, returns self.
    def set_proj(self, n):
        assert(n.dim == self.dim)
        assert(any(n))
        sm = Vector3D.dot(self, n)/Vector3D.dot(n, n)
        for i in range(self.dim):
            self.__v[i] = sm * n.get(i)
        return self


    # set_perp(n) sets self to the perpendicular of self onto n, returns self.
    def set_perp(self, n):
        assert(n.dim == self.dim)
        assert(any(n))
        sm = Vector3D.dot(self, n)/Vector3D.dot(n, n)
        self.set_sub(Vector3D.smult(sm, n))
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
        u = [vs[0].get(i) for i in range(3)]
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


    # proj(v, n) projects v onto n, returns the new projected vector.
    @classmethod
    def proj(cls, v, n):
        u = Vector3D(*v)
        u.set_proj(n)
        return u


    # perp(v, n) takes the perpendicular of v onto n, returns the new
    #       perpendicular vector.
    @classmethod
    def perp(cls, v, n):
        u = Vector3D(*v)
        u.set_perp(n)
        return u


    # norm(v) takes the norm of v, returns a real number.
    @classmethod
    def norm(cls, v):
        return v.get_norm()


    # dist(v, u) finds the distance from v to u, returns a real number.
    @classmethod
    def dist(cls, v, u):
        return sqrt(Vector3D.dot(v, u))

        
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
