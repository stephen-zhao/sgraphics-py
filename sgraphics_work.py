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


    # build_1(n) builds the R^n vector with all components at 1
    @classmethod
    def build_1(cls, n):
        return cls(*(1,)*n)


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
        return sqrt(Vector.dot(self, self))


    # get_dist(u) gets the distance from self to u.
    def get_dist(self, u):
        return sqrt(Vector.dot(self, u))


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
        return v.get_norm()


    # dist(v, u) finds the distance from v to u, returns a real number.
    @classmethod
    def dist(cls, v, u):
        return sqrt(Vector.dot(v, u))



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
