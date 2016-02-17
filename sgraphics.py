from PIL import Image, ImageFilter
import copy
import math

################################################################################
#####                                          #################################
#####   sGraphics , for 3-dimensional spaces   #################################
#####                                          #################################
################################################################################


####### Vector Functions #######################################################

def v_add(v, u):
    assert(len(v) == len(u))
    for i in range(len(v)):
        v[i] = v[i] + u[i]
    return v


def v_sub(v, u):
    assert(len(v) == len(u))
    for i in range(len(v)):
        v[i] = v[i] - u[i]
    return v


def v_dot(v1, v2):
    assert(len(v1) == len(v2))
    sum = 0
    for i in range(len(v1)):
        sum += v1[i] * v2[i]
    return sum


def v_smult(s, v):
    v[:] = map(lambda x: x*s, v)
    return v


def v_sdiv(s, v):
    assert(s)
    v[:] = list(map(lambda x: x/s, v))
    return v


def v_norm(v):
    return sqrt(v_dot(v, v))


def v_proj(n, v):
    assert(any(n))
    c = v_dot(n, v)/v_dot(n, n)
    v[:] = list(map(lambda x: x*c, n))
    return v


def v_perp(n, v):
    assert(any(n))
    c = v_dot(n, v)/v_dot(n, n)
    u = list(map(lambda x: x*c, n))
    v_sub(v, u)
    return v



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


####### Matrix Functions #######################################################

def m_add(M1, M2):
    M = [0]*len(M1)
    for i in range(len(M1)):
        M[i] = v_add(M1[i], M2[i])
    return M


def m_sub(M1, M2):
    M = [0]*len(M1)
    for i in range(len(M1)):
        M[i] = v_sub(M1[i], M2[i])
    return M


def m_smult(s, M):
    return list(map(lambda r: v_smult(s, r), M))


def m_sdiv(s, M):
    return list(map(lambda r: v_sdiv(s, r), M))


def m_transpose(A):
    B = [[] for i in range(len(A[0]))]
    for i in range(len(B)):
        for j in range(len(A)):
            B[i].append(A[j][i])
    return B


def m_mult(A, B1):
    assert(len(A[0]) == len(B1))
    B = m_transpose(B1)
    C = []
    for a in A:
        C.append([])
        for b in B:
            C[-1].append(v_dot(a, b))
    return C


def m_multAll(li_M):
    M = li_M.pop(0)
    while (li_M):
        M = m_mult(M, li_M.pop(0))
    return M


def m_getColFactor(r1, r2, j):
    return r1[j]/r2[j]


def m_findAndSwapRowNo0(A, row, col):
    if row != (len(A)-1):
        for i in range(row+1, len(A)):
            if A[i][col] != 0:
                temp = A[i]
                A[i] = A[row]
                A[row] = temp
                return 0
    return 1


def m_rowReduce(A):
    M = copy.deepcopy(A)
    # i, j is the coordinate of a pivot
    j = 0
    for i in range(len(M)):
        while ((j < len(M[0])) and (M[i][j] == 0)):
            is0 = m_findAndSwapRowNo0(M, i, j)
            if is0:
                j+=1

        if j == len(M[0]):
            break
        
        M[i] = v_sdiv(M[i][j], M[i])

        for i1 in range(0, len(M)):
            if i1 != i:
                M[i1] = v_sub(M[i1], v_smult(m_getColFactor(M[i1], M[i], j) , M[i]))
        j+=1
        
    return M


def m_display(M, padding=   8):
    entry = "%"+str(padding)+".1f  "
    if len(M) == 1:
        print("[  ", end='')
        for a in M[0]:
            print(entry % a, end='')
        print("]")
    else:
        print("/  ", end='')
        for a in M[0]:
            print(entry % a, end='')
        print("\\")
        for i in range(1, len(M)-1):
            print("|  ", end='')
            for a in M[i]:
                print(entry % a, end='')
            print("|")
        print("\\  ", end='')
        for a in M[-1]:
            print(entry % a, end='')
        print("/")

        

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


####### Linear Transformation Standard Matrices ################################

DIREC_X = 0
DIREC_Y = 1
DIREC_Z = 2


def vtoRow(v):
    return [v]


def vtoCol(v):
    return [[x] for x in v]


def v_makeTransAble(v):
    return v + [1]


def m_makeTransAble(M):
    M += [[0 for x in M[0]]]
    for v in M:
        v += [0]
    M[-1][-1] = 1
    return M


def sm_identity(transAble=True, dim=3):
    return [[1 if i == j else 0 for j in range(dim+1)] for i in range(dim+1)]


def sm_proj(n, transAble=True, dim=3):
    M = m_transpose([v_smult(ni/v_dot(n,n),n) for ni in n])
    if transAble:
        return m_makeTransAble(M)
    else:
        return M
    

def sm_scale(v, transAble=True, dim=3):
    if transAble:
        v = v_makeTransAble(v)
        M = [[v[i] if j == i else 0 for j in range(dim+1)] for i in range(dim+1)]
        M[-1][-1] = 1
        return M
    else:
        return [[v[i] if j == i else 0 for j in range(dim)] for i in range(dim)]


def sm_translation(v, dim=3):
    M = [[1 if j == i else 0 for j in range(dim+1)] for i in range(dim+1)]
    for i in range(len(v)):
        M[i][-1] = v[i]
    return M


def sm_rotationd(alpha, beta, gamma, transAble=True, dim=3):
    assert(dim == 3)
    return sm_rotationr(math.radians(alpha), math.radians(beta), math.radians(gamma), transAble, dim)

        
def sm_rotationr(alpha, beta, gamma, transAble=True, dim=3):
    assert(dim == 3)
    Rz = [[math.cos(alpha), -math.sin(alpha), 0],
          [math.sin(alpha),  math.cos(alpha), 0],
          [              0,                0, 1]]
    Ry = [[ math.cos(beta), 0,   math.sin(beta)],
          [              0, 1,                0],
          [-math.sin(beta), 0,   math.cos(beta)]]
    Rx = [[1,               0,                0],
          [0, math.cos(gamma), -math.sin(gamma)],
          [0, math.sin(gamma),  math.cos(gamma)]]
    M = m_multAll([Rz, Ry, Rx])
    if transAble:
        return m_makeTransAble(M)
    else:
        return M


def sm_applyTransform(M, v, transAble=True):
    assert(len(M[0]) == len(v))
    u = copy.deepcopy(v)
    for i in range(len(v)):
        v[i] = v_dot(M[i], u)
    return v







class Mesh:
    def __init__(self, triangles, offset=0):
        self.triangles = triangles
        self.offset = offset



class Cube(Mesh):
    def __init__(self, edgeLength=1, offset=0):
        self.edgeLength = edgeLength
        self.vertices = cls.getVertices(edgeLength)
        self.edges = cls.getEdges(self.vertices)

    def getVertices(cls, s):
        return [[s/2,s/2,s/2],[-s/2,s/2,s/2],[s/2,-s/2,s/2],[s/2,s/2,-s/2],
                [s/2,-s/2,-s/2],[-s/2,s/2,-s/2],[-s/2,-s/2,s/2],[-s/2,-s/2,-s/2]]

    def getEdges(cls, vs):
        return 0
    
        
        
    

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

B = [[2, 3, 5],[6, 7, 8], [1, 1, 2]]
C = [[2, 4, -2, 2],[3, 5, -5, 1],[0, 2, 1, -2]]
D = [[1, 1, 2, 1, 3],[1, 2, 4, 1, 7], [1, 0, 0, 1, -21]]
E = [[1, 4, 2, 0],[2, 5, 1, 0], [3, 6, 0, 0]]
F = [[14, 10000, 1],[9, 1, -9], [10, 10, 10]]
G = [[1, 2, 3], [4, 4, 2]]
H = [[1, 0],[1, 3], [1, 2]]



transformations = [[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1]]



img = Image.new('RGB', (255, 255), "black")

point = [0, 0, 0]
camera_normal = [1, 0, 0]
camera_position = [-8, 0, 0]

a = [1, 2, 12]
n = [2, 2, 2]

print(v_proj(n, a))




##pixels = img.load()
##
##for i in range(img.size[0]):
##    for j in range(img.size[1]):
##        pixels[i, j] = (i, j, j)
##
##
##img.show()
