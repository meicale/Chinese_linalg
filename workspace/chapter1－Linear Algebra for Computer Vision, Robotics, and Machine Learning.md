Chapter 1 

Introduction 

As we explaffined in the preface, this first volume covers "classical" linear algebra, up to and including the primary decomposition and the Jordan form. Besides covering the standard topics, we discuss a few topics that are important for applications. These include: 
1. 
Haar bases and the corresponding Haar wavelets, a fundamental tool in signal process-ing and computer graphics. 

2. 
Hadamard matrices which have applications in error correcting codes, signal processing, and low rank approximation. 

3. 
Affine maps (see Section 4.4). 	These are usually ignored or treated in a somewhat obscure fashion. Yet they play an important role in computer vision and robotics. There is a clean and elegant way to define affine maps. One simply has to define affine combinations. Linear maps preserve linear combinations, and similarly affine maps preserve affine combinations. 

4. 
Norms and matrix norms (Chapter 8). These are used extensively in optimization theory. 

5. 
Convergence of sequences and series in a normed vector space. Banach spaces (see Section 8.7). The matrix exponential $e^A$ and its basic properties (see Section 8.8). In particular, we prove the Rodrigues formula for rotations in SO(3) and discuss the surjectivity of the exponential map exp: $so(3) \to SO(3)$, where so(3) is the real vector space of 3x3 skew symmetric matrices (see Section 11.7). We also show that $det(e^A) = e^{tr(A)}$ (see Section 14.5). 


6. 
The group of unit quaternions, 	SU(2), and the representation of rotations in SO(3) by unit quaternions (Chapter 15). We define a homomorphism r : $SU(2) \to SO(3)$ and prove that it is surjective and that its kernel is {-I,I}. We compute the rota-tion matrix $R_q$ associated with a unit quaternion q, and give an algorithm to construct a quaternion from a rotation matrix. We also show that the exponential map exp: $su(2) \to SU(2)$ is surjective, where su(2) is the real vector space of skew-Hermitian 2x2 matrices with zero trace. We discuss quaternion interpolation and prove the famous slerp interpolation formula due to Ken Shoemake. 

7. 
An introduction to algebraic and spectral graph theory. We define the graph Laplacian and prove some of its basic properties (see Chapter 18). In Chapter 19, we explain how the eigenvectors of the graph Laplacian can be used for graph drawing. 

8. 
Applications of SVD and pseudo-inverses, in particular, principal component analysis, for short PCA (Chapter 21). 

9. 
Methods for computing eigenvalues and eigenvectors are discussed in Chapter 17. We first focus on the QR algorithm due to Rutishauser, Francis, and Kublanovskaya. See Sections 17.1 and 17.3. We then discuss how to use an Arnoldi iteration, in combination with the QR algorithm, to approximate eigenvalues for a matrix A of large dimension. See Section 17.4. The special case where A is a symmetric (or Hermitian) tridiagonal matrix, involves a Lanczos iteration, and is discussed in Section 17.6. In Section 17.7, we present power iterations and inverse (power) iterations. 

Five topics are covered in more detail than usual. These are 
1. 
Matrix factorizations such as LU, PA = LU, Cholesky, and reduced row echelon form (rref). Deciding the solvablity of a linear system Ax = b, and describing the space of solutions when a solution exists. See Chapter 7. 

2. 
Duality (Chapter 10). 

3. 
Dual norms (Section 13.7). 

4. 
The geometry of the orthogonal groups O(n) and SO(n), and of the unitary groups U(n) and SU(n). 

5. 
The spectral theorems (Chapter 16). 

Most texts omit the proof that the PA = LU factorization can be obtaffined by a simple modi.cation of Gaussian elimination. We give a complete proof of Theorem 7.5 in Section 7.6. We also prove the uniqueness of the rref of a matrix; see Proposition 7.19. 

At the most basic level, duality corresponds to transposition. But duality is really the bijection between subspaces of a vector space E (say .nite-dimensional) and subspaces of linear forms (subspaces of the dual space E.) established by two maps: the first map assigns to a subspace V of E the subspace V 0 of linear forms that vanish on V ; the second map assigns to a subspace U of linear forms the subspace U0 consisting of the vectors in E on which all linear forms in U vanish. The above maps define a bijection such that dim(V ) + dim(V 0)= dim(E), dim(U) + dim(U0) = dim(E), V 00 = V , and U00 = U. 

Another important fact is that if E is a .nite-dimensional space with an inner product $u, v \longmapsto\langle u, v\rangle$ (or a Hermitian inner product if E is a complex vector space), then there is a canonical isomorphism between E and its dual $E^*$. . This means that every linear form $f\in E^*$ is uniquely represented by some vector $u\in E$, in the sense that $f(v) = <v,u>$ for all $v\in E$. As a consequence, every linear map f has an adjoint $f^*$ such that <f(u),v)> = <u, f^*(v)> for all $u,v \in E$. 

Dual norms show up in convex optimization; see Boyd and Vandenberghe [11]. 

Because of their importance in robotics and computer vision, we discuss in some detail the groups of isometries O(E) and SO(E) of a vector space with an inner product. The isometries in O(E) are the linear maps such that $f \circ f^{*}=f^{*} \circ f=\mathrm{id}$, and the direct isometries in SO(E), also called rotations, are the isometries in O(E) whose determinant is equal to +1. We also discuss the hermitian counterparts U(E) and SU(E). 

We prove the spectral theorems not only for real symmetric matrices, but also for real and complex normal matrices. 

We stress the importance of linear maps. Matrices are of course invaluable for computing and one needs to develop skills for manipulating them. But matrices are used to represent a linear map over a basis (or two bases), and the same linear map has di.erent matrix representations. In fact, we can view the various normal forms of a matrix (Schur, SVD, Jordan) as a suitably convenient choice of bases. 

We have listed most of the Matlab functions relevant to numerical linear algebra and have included Matlab programs implementing most of the algorithms discussed in this book. 