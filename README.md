# sgraphics-py



## Brief Summary
sgraphics-py is my own attempt at creating a set of classes and methods to perform matrix operations in 3D vectorspaces on meshes. Python, being a versatile and adaptive language, was the perfect choice to do the prototyping. Once complete, there will be provided native functionality for vectors and their operations, matrices and their operations, the creation of mesh objects, and an API for rendering them to the screen (using OpenGL). I may switch to C++ to complete the project once prototyping is complete. Below is a list of the provided functionality so far.




## Classes
The documentation below indicates all implemented and some planned functionality. TODO's are marked italicized and marked with @TODO

_@TODO: Vector and Vector3D might be renamed to V and V3d, respectively, for shorter names. This should actually benefit readability and allow faster code development._


### class Vector:
Vector represents vectors in R^n, providing related vector operations and methods.

#### CONSTRUCTORS: (class methods)
| Function |  Description |
|---------:|--------------|
| `Vector(*v)` | takes a variable number of positional arguments to instantiate the vector with components `*v` |
| `build_0(n)` | builds the R^`n` zero-vector |
| `build_rep(n, x)` | builds the R^`n` vector with all components filled with `x` |
| `build_seq(a, b, by=1)` | builds a vector containing the sequence from `a` up to, but not including, `b`, incrementing by `by` each subsequent component |
| `build_count(n, start=0, by=1)` | builds the R^`n` vector containing the sequence starting with `start`, incrementing by `by` each subsequent component |
| `build_ei(n, i)` | builds the R^`n` basis vector (e_i) in direction `i` |

#### GETTERS: (instance methods)
| Function |  Description |
|---------:|--------------|
| `get(i)` | gets `self`'s `i`th coordinate |
| `get_list()` | gets `self` as a list |
| `get_norm()` | gets `self`'s norm (square root of dot product with itself) (distance from `self` to origin) |
| `get_dist(u)` | gets the distance from `self` to `u` |

#### SETTER OPERATORS: (instance methods)
These setter operators work inplace (in the instance's memory location) and return self, to allow chaining of operations. These will be useful when performing transformations on vertices of meshes while ensuring all references to the vertices are changed as well.

For example, `Vector(1,1,1).set_add(Vector(1,3,5)).set_sub(Vector(4,2,9)).set_smult(2) -> [-4  4 -6]`

| Function |  Description |
|---------:|--------------|
| `set_add(u)` | sets `self` to `self`+`u`, returns `self` |
| `set_sub(u)` | sets `self` to `self`-`u`, returns `self` |
| `set_smult(s)` | sets `self` to `self`*(scalar `s`), returns `self` |
| `set_sdiv(s)` | sets `self` to `self`/(scalar `s`), returns `self` |
| `set_proj(a)` | sets `self` to `self` projected onto `a`, returns `self` |
| `set_perp(a)` | sets `self` to `self` taken perpendicular to `a`, returns `self` |

#### OPERATIONS: (class methods)
| Function |  Description |
|---------:|--------------|
| `add(*vs)` | adds one or more vectors in `vs` together, produces the new sum vector |
| `sub(*vs)` | subtracts one or more vectors in `vs[1:]` from `v[0]`, produces the new difference vector |
| `smult(s, v)` | multiplies `v` by scalar `s`, returns the new multiple vector |
| `sdiv(s, v)` | divides `v` by scalar `s`, returns the new quotient vector |
| `dot(v, u)` | takes the dot product of `v` and `u`, returns a float |
| `proj(v, a)` | projects `v` onto `a`, returns the new projected vector |
| `perp(v, a)` | takes the perpendicular of `v` to `a`, returns the new perpendicular vector |
| `norm(v)` | takes the norm of `v`, returns a float |
| `dist(v, u)` | takes the distance from `v` to `u`, returns a float |

#### DISPLAY / FORMATTING: (instance methods)
| Function |  Description |
|---------:|--------------|
| `toString(col=False)` | returns a string representation of `self`, either as a column vector (`col=True`) or as a row vector (`col=False`) |
| `print(col=False)` | prints `self` to the screen, either as a column vector (`col=True`) or as a row vector (`col=False`) |


### class Vector3D:
Vector3D represents vectors in R^3, providing related vector operations and methods. _@TODO: Refactor Vector3D to be a subclass of Vector, or remove it entirely and move functionality over to class Vector_

#### CONSTRUCTORS: (class methods)
| Function |  Description |
|---------:|--------------|
| `Vector3D(*v)` | takes 3 positional arguments to instantiate the vector with *v x,y,z-coordinates |
| `build_0()` | builds the R^3 zero-vector [0 0 0] |
| `build_1()` | builds the R^3 vector of all 1's [1 1 1] |
| `build_i()` | builds the i unit vector [1 0 0] |
| `build_j()` | builds the j unit vector [0 1 0] |
| `build_k()` | builds the k unit vector [0 0 1] |
| `build_fromSpherical(r, theta, phi)` | builds the vector with norm r with theta azimuth and phi inclination. |
| _`build_fromEulerAngles(r, alpha, beta, gamma)` | @TODO: builds the vector with norm r with alpha, beta, gamma euler angles_ |

#### GETTERS: (instance methods)
| Function |  Description |
|---------:|--------------|
| `get_x()` | gets self's x-coordinate |
| `get_y()` | gets self's y-coordinate |
| `get_z()` | gets self's z-coordinate |
| `get(i)` | gets self's ith coordinate |
| `get_norm()` | gets self's norm (square root of dot product with itself) (distance from self to origin) |
| `get_dist(u)` | gets the distance from self to u |

#### SETTER OPERATORS: (instance methods)
These setter operators work inplace (in the instance's memory location) and return self, to allow chaining of operations. These will be useful when performing transformations on vertices of meshes while ensuring all references to the vertices are changed as well.

For example, `Vector3D.build_1().set_add(Vector3D(1,3,5)).set_sub(Vector3D(4,2,9)).set_smult(2) -> [-4  4 -6]`

| Function |  Description |
|---------:|--------------|
| `set_add(u)` | sets self to self+u, returns self |
| `set_sub(u)` | sets self to self-u, returns self |
| `set_smult(s)` | sets self to self*(scalar s), returns self |
| `set_sdiv(s)` | sets self to self/(scalar s), returns self |
| `set_cross(u)` | sets self to self cross multiplied with u, returns self |
| `set_proj(n)` | sets self to self projected onto n, returns self |
| `set_perp(n)` | sets self to self taken perpendicular to n, returns self |

#### OPERATIONS: (class methods)
| Function |  Description |
|---------:|--------------|
| `add(*vs)` | adds one or more vectors in vs together, produces the new sum vector |
| `sub(*vs)` | subtracts one or more vectors in vs[1:] from v[0], produces the new difference vector |
| `smult(s, v)` | multiplies v by scalar s, returns the new multiple vector |
| `sdiv(s, v)` | divdes v by scalar s, returns the new quotient vector |
| `dot(v, u)` | takes the dot product of v and u, returns a float |
| `cross(*vs)` | takes the sequential cross products of all vectors in vs, returns the new resulting product vector |
| `proj(v, n)` | projects v onto n, returns the new projected vector |
| `perp(v, n)` | takes the perpendicular of v to n, returns the new perpendicular vector |
| `norm(v)` | takes the norm of v, returns a float |
| `dist(v, u)` | takes the distance from v to u, returns a float |

#### DISPLAY / FORMATTING: (instance methods)
| Function |  Description |
|---------:|--------------|
| `toString(col=False)` | returns a string representation of self, either as a column vector (`col=True`) or as a row vector (`col=False`) |
| `print(col=False)` | prints self to the screen, either as a column vector (`col=True`) or as a row vector (`col=False`) |


### class Matrix:
_@TODO: functions for determining consistency, rank, nullity, solving and returning solutions_ (Many have already been implemented procedurally that require rewriting to object-oriented)

#### CONSTRUCTORS: (class methods)
The default constructor is private. Use the following functions to make a matrix.

| Function |  Description |
|---------:|--------------|
| `build_fromRows(*rs)` | builds a matrix from lists representing the rows of the matrix. |
| `build_fromCols(*cs)` | builds a matrix from lists representing the columns of the matrix. |
| `build_fromVRows(*rs)` | builds a matrix from vectors representing the rows of the matrix. |
| `build_fromVCols(*rs)` | builds a matrix from vectors representing the columns of the matrix. |

#### GETTERS: (instance methods)
| Function |  Description |
|---------:|--------------|
| `get(i, j)` | gets `self`'s (`i`, `j`)th entry |
| `get_array(prioritizeRow=True)` | gets `self` as a list of lists, either prioritizing rows or columns; row by default |
| `get_lVCols()` | gets `self` as a list of column vectors |
| `get_lVRows()` | gets `self` as a list of row vectors |
| `get_flattened(prioritizeRow=True)` | gets `self` as a flattened list, prioritizing rows or columns; row by default |
| `get_nRows()` | gets the number of rows in `self` |
| `get_nCols()` | gets the number of columns in `self` |
| `get_dim()` | gets the number of rows and number of columns as a tuple |

#### OPERATIONS: (instance methods)
| Function |  Description |
|---------:|--------------|
| `add(A)` | adds matrix `A` to `self`; returns the new sum matrix |
| `sub(A)` | subtracts matrix `A` from `self`; returns the new difference matrix |
| `smult(s)` | multiplies `self` by scalar `s`; returns the new scaled matrix |
| `sdiv(s)` | divides `self` by scalar `s`; returns the new quotient vector |
| `transpose()` | returns `self`'s transpose matrix |
| `postmult(A)` | multiplies `self` by `A`; returns the new product matrix |
| `premult(A)` | multiplies `A` by `self`; returns the new product matrix |
| `rowReduce()` | returns the Row Reduced Echelon Form of self |

#### DISPLAY / FORMATTING: (instance methods)
| Function |  Description |
|---------:|--------------|
| `toString()` | returns a string representation of `self` |
| `print(col=False)` | prints `self` to screen |



### class Transformation:
I'm unsure whether this will stay a class or be merged into the general matrix class
_@TODO: functions to generate standard matrices for scaling, rotation, translation for R^3 (4x4 matrices), projection, perpendicular, identity, etc._ (Many have already been implemented procedurally that require rewriting to object-oriented)



### class Mesh:
A class to superclass all specific mesh classes such as Cube, Cuboid, Parallelepiped, Sphere, Cylinder, etc. Methods to generate lists of edges, triangles, vertices, from given parameters in construction. No plans for textures yet. I'm still reading into that.

