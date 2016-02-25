import math


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# CLASS Vector: Represents vectors. Provides corresponding vector operations
#               and related methods.
################################################################################
################################################################################
################################################################################

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
        self.__dim = len(v)


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
    #          another vector, the dimension of vector, and the iterator
    #          function.
    

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


    # get_dim() gets self's dimension (number of components).
    def get_dim(self):
        return self.__dim


    # __iter__() returns the list holding the vector coordinates.
    def __iter__(self):
        return iter(self.__v)


    
    ############################################################################
    # SETTER OPERATIONS: Methods to operate and mutate self. Includes addition,
    #                    subtraction, scalar multiplcation & division,
    #                    projection, and perpendicular, as well as adding and
    #                    removing dimensions. These methods operate in the same
    #                    memory location and mutate the operand. They return
    #                    self to allow chaining of operations.

    # set_add(u) sets self to self + u, returns self.
    def set_add(self, u):
        assert(u.get_dim() == self.__dim)
        for i in range(self.__dim):
            self.__v[i] += u.get(i)
        return self
    

    # set_sub(u) sets self to self - u, returns self.
    def set_sub(self, u):
        assert(u.get_dim() == self.__dim)
        for i in range(self.__dim):
            self.__v[i] -= u.get(i)
        return self
    

    # set_smult(s) sets self to (scalar s) * self, returns self.
    def set_smult(self, s):
        for i in range(self.__dim):
            self.__v[i] *= s
        return self


    # set_sdiv(s) sets self to self / (scalar s), returns self.
    def set_sdiv(self, s):
        assert(s != 0)
        for i in range(self.__dim):
            self.__v[i] /= s
        return self


    # set_proj(a) sets self to the projection of self onto a, returns self.
    def set_proj(self, a):
        assert(a.get_dim() == self.__dim)
        assert(any(a))
        sm = Vector.dot(self, a)/Vector.dot(a, a)
        for i in range(self.__dim):
            self.__v[i] = sm * a.get(i)
        return self


    # set_perp(a) sets self to the perpendicular of self onto a, returns self.
    def set_perp(self, a):
        assert(a.get_dim() == self.__dim)
        assert(any(a))
        sm = Vector.dot(self, a)/Vector.dot(a, a)
        self.set_sub(Vector.smult(sm, a))
        return self


    # set_addDim(init=0) sets self to self with an additional appended dimension
    #       initialized to init.
    def set_addDim(self, init=0):
        self.__v += [init]
        self.__dim += 1
        return self


    # set_remDim(n=-1) sets self to self with dimension n removed; default
    #       dimension is the last dimension.
    def set_remDim(self, n=-1):
        self.__v.pop(n)
        self.__dim -= 1
        return self


        
    
    ############################################################################
    # OPERATIONS: Class functions that operate on real numbers and Vector 
    #             objects and return references to new Vectors / real
    #             numbers. These include functions for addition, subtraction, 
    #             scalar multiplcation & division, projection, perpendicular, 
    #             dot product and norm/distance.

    # add(*vs) adds one or more vectors in vs together, returns the new sum
    #       vector.
    @classmethod
    def add(cls, *vs):
        v_dim = vs[0].get_dim()
        assert(v.get_dim() == v_dim for v in vs)
        u = [0 for i in range(v_dim)]
        for i in range(v_dim):
            u[i] = sum(map(lambda v: v.get(i), vs))
        return Vector(*tuple(u))


    # sub(*vs) subtracts one or more vectors in vs[1:] from vs[0], returns the
    #       new difference vector.
    @classmethod
    def sub(cls, *vs):
        v_dim = vs[0].get_dim()
        assert(v.get_dim() == v_dim for v in vs)
        u = [vs[0].get(i) for i in range(v_dim)]
        for i in range(v_dim):
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
        v_dim = v.get_dim()
        assert(v_dim == u.get_dim())
        return sum([v.get(i)*u.get(i) for i in range(v_dim)])


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
            s = "[  "+"{:>7.6G}  "*self.__dim+"]"
        else:
            s = "/ {:>7.6G} \\\n" + \
                "| {:>7.6G} |\n"*(self.__dim-2) + \
                "\\ {:>7.6G} /"
        return s.format(*tuple(self.__v))

    def print(self, col=False):
        print(self.toString(col))






        


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# CLASS Vector3D: Represents vectors in R^3. Provides corresponding vector
#                 operations and related methods.
################################################################################
################################################################################
################################################################################

# TODO: - refactor Vector3D to be subclass of Vector; many methods can be
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

        
    ############################################################################
    # DISPLAY / FORMATTING: Methods/functions that format Vector3Ds for output
    def toString(self, col=False):
        if not col:
            return "[  {:>9.8G}  {:>9.8G}  {:>9.8G}  ]"\
                   .format(*tuple(self.__v))
        else:
            return "/ {:>9.8G} \\\n| {:>9.8G} |\n\\ {:>9.8G} /"\
                   .format(*tuple(self.__v))

    def print(self, col=False):
        print(self.toString(col))





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# CLASS Matrix: Represents matrices of arbitrary size.
# TODO: - entire thing
################################################################################
################################################################################
################################################################################
class Matrix:

    ############################################################################
    # CLASS VARIABLES:
    
    __RESTRICTION = 0x1F8EAA2

    # TODO: add class vars


    ############################################################################
    # CONSTRUCTORS:
    def __init__(self, token, nRows, nCols, *rs):
        assert(token == self.__RESTRICTION)
        self.__nRows = nRows
        self.__nCols = nCols
        self.Ma = [r for r in rs]
        
                
    # build_fromRows(*rs) builds the matrix with rows from the lists in rs
    @classmethod
    def build_fromRows(cls, *rs):
        r0_len = len(rs[0])
        assert(len(r) == r0_len for r in rs)
        return cls(cls.__RESTRICTION, len(rs), r0_len, *rs)
        

    # build_fromCols(*cs) builds the matrix with columns from the lists in cs
    @classmethod
    def build_fromCols(cls, *cs):
        c0_len = len(cs[0])
        assert(len(c) == c0_len for c in cs)
        return cls(cls.__RESTRICTION, c0_len, len(cs), \
                   *(cls.array_transpose(cs)))


    # build_fromVRows(*rs) builds the matrix with rows from row Vectors in rs
    @classmethod
    def build_fromVRows(cls, *rs):
        r0_dim = rs[0].get_dim()
        assert(r.get_dim() == r0_dim for r in rs)
        return cls(cls.__RESTRICTION, len(rs), r0_dim, \
                   *(r.get_list() for r in rs))


    # build_fromVCols(*cs) builds the matrix with cols from col Vectors in cs
    @classmethod
    def build_fromVCols(cls, *cs):
        c0_dim = cs[0].get_dim()
        assert(c.get_dim() == c0_dim for c in cs)
        rs = [[] for i in range(c0_dim)]
        for i in range(c0_dim):
            for j in range(len(cs)):
                rs[i].append(cs[j].get(i))
        return cls(cls.__RESTRICTION, c0_dim, len(cs), *rs)

    
    ############################################################################
    # GETTERS: Functions to get the matrix as an array (list of row lists or
    #          list of column lists), a list of column vectors, a list of row
    #          vectors, a flattened list (row or column prioritized), and the
    #          dimensions of the matrix.


    # get_array(prioritizeRow=True) gets the matrix as a list of lists,
    #       prioritizing either row or col, row by default
    def get_array(self, prioritizeRow=True):
        if prioritizeRow:
            return self.Ma
        else:
            return self.array_transpose(self.Ma)


    # get_lVCols() gets the matrix as a list of column vectors
    def get_lVCols(self):
        return [Vector(*c) for c in self.array_transpose(self.Ma)]



    # get_lVRows() gets the matrix as a list of row vectors
    def get_lVRows(self):
        return [Vector(*r) for r in self.Ma]



    # get_flattened(prioritizeRow=True) gets the matrix as a flattened list,
    #       prioritizing either row or col, row by default
    def get_flattened(self, prioritizeRow=True):
        if prioritizeRow:
            return [a for r in self.Ma for a in r]
        else:
            return [a for r in self.array_transpose(self.Ma) for a in r]


    # get_nRows() gets the number of rows of self
    def get_nRows(self):
        return self.__nRows
    

    # get_nCols() gets the number of columns of self
    def get_nCols(self):
        return self.__nCols


    # get_dim() gets a tuple containing (nRows, nCols)
    def get_dim(self):
        return (self.__nRows, self.__nCols)



    ############################################################################
    # OPERATIONS: Instance operations at apply operations on self and return
    #             references to new matrices.

    # add(A) adds matrix A with same dimensions as self to self, returns a new
    #       matrix.
    def add(self, A):
        assert(self.get_dim() == A.get_dim())
        A_arr = A.get_array()
        rs = [[self.Ma[i][j] + A_arr[i][j] for j in range(self.__nCols)] \
              for i in range(self.__nRows)]
        return Matrix.build_fromRows(*rs)

    
    # sub(A) subtracts matrix A with same dimensions as self from self, returns
    #       a new matrix.
    def sub(self, A):
        assert(self.get_dim() == A.get_dim())
        A_arr = A.get_array()
        rs = [[self.Ma[i][j] - A_arr[i][j] for j in range(self.__nCols)] \
              for i in range(self.__nRows)]
        return Matrix.build_fromRows(*rs)


    # smult(s) multiplies self by a scalar s, returns a new matrix.
    def smult(self, s):
        rs = [[s*self.Ma[i][j] for j in range(self.__nCols)] \
              for i in range(self.__nRows)]
        return Matrix.build_fromRows(*rs)

    
    # sdiv(s) divides self by a scalar s, returns a new matrix.
    def sdiv(self, s):
        assert(s != 0)
        rs = [[self.Ma[i][j]/s for j in range(self.__nCols)] \
              for i in range(self.__nRows)]
        return Matrix.build_fromRows(*rs)


    # transpose() returns self's transpose matrix.
    def transpose(self):
        return Matrix.build_fromCols(*self.Ma)


    # postmult(A) multiplies self by A (with appropriate dimensions), and
    #       returns the new matrix
    def postmult(self, A):
        assert(self.__nCols == A.get_nRows())
        selfRows = self.get_lVRows()
        ACols = A.get_lVCols()
        rs = [[Vector.dot(row, col) for col in ACols] for row in selfRows]
        return Matrix.build_fromRows(*rs)

    
    # premult(A) multiplies A (with appropriate dimensions) by self, and
    #       returns the new matrix
    def premult(self, A):
        assert(A.get_nCols() == self.__nRows)
        ARows = A.get_lVRows()
        selfCols = self.get_lVCols()
        rs = [[Vector.dot(row, col) for col in selfCols] for row in ARows]
        return Matrix.build_fromRows(*rs)


    # rowReduce() determines the row-reduced echelon form of self, and returns
    #       the new matrix
    def rowReduce(self):
        lRows = self.get_lVRows()
        # i, j is the coordinate of a candidate pivot
        j = 0
        for i in range(self.__nRows):
            while ((j < self.__nCols) and (lRows[i].get(j) == 0)):
                is0 = self.__findAndSwapRowWithNo0(lRows, i, j)
                if is0:
                    j+=1

            if j == self.__nCols:
                break

            lRows[i].set_sdiv(lRows[i].get(j))

            for i1 in range(self.__nRows):
                if i1 != i:
                    lRows[i1].set_sub(Vector.smult(lRows[i1].get(j), lRows[i]))

            j+=1

        return Matrix.build_fromVRows(*lRows)
            
        

    ###########################################################################
    # DISPLAY / FORMATTING: Methods/functions that format Matrices for output
    def toString(self, padding=8, precision=1):
        if self.__nRows == 1:
            s = "[ "+"{:>7.6G} "*self.__nCols+"]"
        else:
            s = "/ "+"{:>7.6G} "*self.__nCols+"\\\n" + \
                ("| "+"{:>7.6G} "*self.__nCols+"|\n")*(self.__nRows - 2) + \
                "\\ "+"{:>7.6G} "*self.__nCols+"/"
        return s.format(*tuple(self.get_flattened()))


    def print(self):
        print(self.toString())


    ############################################################################
    # UTILITY: miscellaneous functions

    
    # array_transpose(arr) tranposes an 'array' in the form of a list of lists,
    #       swapping rows with columns
    @classmethod
    def array_transpose(cls, arr):
        arr_cols = len(arr[0])
        arr_rows = len(arr)
        arrT = [[] for i in range(arr_cols)]
        for i in range(arr_cols):
            for j in range(arr_rows):
                arrT[i].append(arr[j][i])
        return arrT


    # __findAndSwapRowWithNo0(rs, i, j) finds the first row after the ith row 
    #       that has a non-zero jth component and swaps that row with the ith
    #       row; returns 0 on success or 1 on failure
    @classmethod
    def __findAndSwapRowWithNo0(cls, rs, i, j):
        if i != (len(rs)-1):
            for k in range(i+1, len(rs)):
                if rs[k].get(j) != 0:
                    temp = rs[k]
                    rs[k] = rs[i]
                    rs[i] = temp
                    return 0
        return 1
        



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


A = Matrix.build_fromRows([2,4,6,8],[1,2,1,0])
B = Matrix.build_fromRows([1,2,3],[1,4,5],[1,2,0],[0,0,1])
C = Matrix.build_fromRows([0,1,2,-2],[0,3,5,1],[1,2,5,0])

