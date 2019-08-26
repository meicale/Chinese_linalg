Chapter 7 
Gaussian Elimination, LU-Factorization, Cholesky Factorization, Reduced Row Echelon Form 
In this chapter we assume that all vector spaces are over the .eld R. All results that do not rely on the ordering on R or on taking square roots hold for arbitrary .elds. 
7.1 Motivating Example: Curve Interpolation 
Curve interpolation is a problem that arises frequently in computer graphics and in robotics (path planning). There are many ways of tackling this problem and in this section we will describe a solution using cubic splines. Such splines consist of cubic B¡äezier curves. They are often used because they are cheap to implement and give more .exibility than quadratic B¡äezier curves. 
A cubic B¡äezier curve C(t) (in R2 or R3) is speci.ed by a list of four control points (b0,b2,b2,b3) and is given parametrically by the equation 
C(t) = (1 . t)3 b0 + 3(1 . t)2tb1 + 3(1 . t)t2 b2 + t3 b3. 
Clearly, C(0) = b0, C(1) = b3, and for t ¡Ê [0, 1], the point C(t) belongs to the convex hull of the control points b0,b1,b2,b3. The polynomials 
(1 . t)3 , 3(1 . t)2t, 3(1 . t)t2 ,t3 
are the Bernstein polynomials of degree 3. 
Typically, we are only interested in the curve segment corresponding to the values of t in the interval [0, 1]. Still, the placement of the control points drastically a.ects the shape of the curve segment, which can even have a self-intersection; See Figures 7.1, 7.2, 7.3 illustrating various con.gurations. 
191 


Figure 7.1: A ¡°standard¡± B¡äezier curve. 


Figure 7.2: A B¡äezier curve with an in.ection point. 


Figure 7.3: A self-intersecting B¡äezier curve. 
Interpolation problems require .nding curves passing through some given data points and possibly satisfying some extra constraints. 
A B¡äezier spline curve F is a curve which is made up of curve segments which are B¡äezier curves, say C1,...,Cm (m ¡Ý 2). We will assume that F de.ned on [0,m], so that for i =1,...,m, 
F (t)= Ci(t . i + 1),i . 1 ¡Ü t ¡Ü i. 
7.1. MOTIVATING EXAMPLE: CURVE INTERPOLATION 
Typically, some smoothness is required between any two junction points, that is, between any two points Ci(1) and Ci+1(0), for i =1,...,m . 1. We require that Ci(1) = Ci+1(0) (C0-continuity), and typically that the derivatives of Ci at1 andof Ci+1 at 0 agree up to second order derivatives. This is called C2-continuity, and it ensures that the tangents agree as well as the curvatures. 
There are a number of interpolation problems, and we consider one of the most common problems which can be stated as follows: 
Problem: Given N + 1 data points x0,...,xN , .nd a C2 cubic spline curve F such that F (i)= xi for all i,0 ¡Ü i ¡Ü N (N ¡Ý 2). 
A way to solve this problem is to .nd N + 3 auxiliary points d.1,...,dN+1, called de Boor control points, from which N B¡äezier curves can be found. Actually, 
d.1 = x0 and dN+1 = xN 
so we only need to .nd N + 1 points d0,...,dN . 
It turns out that the C2-continuity constraints on the N B¡äezier curves yield only N . 1 equations, so d0 and dN can be chosen arbitrarily. In practice, d0 and dN are chosen according to various end conditions, such as prescribed velocities at x0 and xN . For the time being, we will assume that d0 and dN are given. 
Figure 7.4 illustrates an interpolation problem involving N +1=7+1=8 data points. The control points d0 and d7 were chosen arbitrarily. 

Figure 7.4: A C2 cubic interpolation spline curve passing through the points x0,x1, x2,x3, x4,x5, x6,x7. 
It can be shown that d1,...,dN.1 are given by the linear system 
.
..
..
. 
7 16x1 . 3 
2 d12 d0 
141 0 

d2 
6x2 
...... 

...
...
...... 
...... 

. 

. 

. 

dN.2 
...... 

= 

...... 

. 

. 

. 

6xN.2 
...... 

.

. 

. 

. 

0 141 

1 7 dN.1 6xN.1 . 3 dN
22 
We will show later that the above matrix is invertible because it is strictly diagonally dominant. 
Once the above system is solved, the B¡äezier cubics C1,..., CN are determined as follows (we assume N ¡Ý 2): For 2 ¡Ü i ¡Ü N . 1, the control points (bi 0,b1i ,bi 2,b3i ) of Ci are given by 
bi 
= xi.1
0 
21 
bi 
1 = di.1 + di
33 
12 
bi 
2 = di.1 + di
33 
bi = xi.
3 
The control points (b01,b11,b21,b31) of C1 are given by b1 
= x0
0 b1 1 = d0 
11 
b1 
2 = d0 + d1
22 
b1 
= x1,
3 
and the control points (bN ,bN ,bN ,bN ) of CN are given by 
0123 
bN 
= xN.1
0 
11 
bN 
1 = dN.1 + dN
22 
bN 
2 = dN bN 
= xN .
3 
Figure 7.5 illustrates this process spline interpolation for N = 7. 
We will now describe various methods for solving linear systems. Since the matrix of the above system is tridiagonal, there are specialized methods which are more e.cient than the general methods. We will discuss a few of these methods. 


7.2 Gaussian Elimination 
Let A be an n ¡Á n matrix, let b ¡Ê Rn be an n-dimensional vector and assume that A is invertible. Our goal is to solve the system Ax = b. Since A is assumed to be invertible, 
7.2. GAUSSIAN ELIMINATION 

Figure 7.5: A C2 cubic interpolation of x0,x1, x2,x3, x4,x5, x6,x7 with associated color coded B¡äezier cubics. 
we know that this system has a unique solution x = A.1b. Experience shows that two counter-intuitive facts are revealed: 
(1) 
One should avoid computing the inverse A.1 of A explicitly. This is ine.cient since it would amount to solving the n linear systems Au(j) = ej for j =1,...,n, where ej = (0,..., 1,..., 0) is the jth canonical basis vector of Rn (with a 1 is the jth slot). By doing so, we would replace the resolution of a single system by the resolution of n systems, and we would still have to multiply A.1 by b. 

(2) 
One does not solve (large) linear systems by computing determinants (using Cramer¡¯s formulae) since this method requires a number of additions (resp. multiplications) proportional to (n + 1)! (resp. (n + 2)!). 


The key idea on which most direct methods (as opposed to iterative methods, that look for an approximation of the solution) are based is that if A is an upper-triangular matrix, which means that aij = 0for1 ¡Ü j<i ¡Ü n (resp. lower-triangular, which means that aij =0 for 1 ¡Ü i<j ¡Ü n), then computing the solution x is trivial. Indeed, say A is an upper-triangular matrix 
.
. 
A = 

........ 

a11 a12 ¡¤¡¤¡¤ a1 n.2 a1 n.1 a1 n 
0 ¡¤¡¤¡¤
a22 a2 n.2 a2 n.1 a2 n 
...
.
... .
00... . 
.. .  . . .  . . .  
0  0  ¡¤ ¡¤ ¡¤  0  an.1 n.1  an.1 n  
0  0  ¡¤ ¡¤ ¡¤  0  0  an n  

........ 

. 

Then det(A)= a11a22 ¡¤¡¤¡¤ ann = 0, which implies that aii 0 for i =1,...,n, and we can 
= solve the system Ax = b from bottom-up by back-substitution. That is, .rst we compute xn from the last equation, next plug this value of xn into the next to the last equation and compute xn.1 from it, etc. This yields 
.1 
= 	a
xn nnbn xn.1 = a .1 (bn.1 . an.1 nxn)
n.1 n.1
. 
. 
. 
.1 
x1 = a11 (b1 . a12x2 .¡¤ ¡¤¡¤. a1 nxn). 
Note that the use of determinants can be avoided to prove that if A is invertible then aii 0 for i =1,...,n.
= 	Indeed, it can be shown directly (by induction) that an upper (or 
lower) triangular matrix is invertible i. all its diagonal entries are nonzero. If A is lower-triangular, we solve the system from top-down by forward-substitution. Thus, what we need is a method for transforming a matrix to an equivalent one in upper-
triangular form. This can be done by elimination. Let us illustrate this method on the following example: 
2x + y + z =5 4x . 6y = .2 
.2x +7y +2z =9. We can eliminate the variable x from the second and the third equation as follows: Subtract twice the .rst equation from the second and add the .rst equation to the third. We get the new system 
2x 	+ y + z =5 . 8y . 2z = .12 
8y +3z = 14. This time we can eliminate the variable y from the third equation by adding the second equation to the third: 
2x 	+ y + z =5 . 8y . 2z = .12 
z =2. This last system is upper-triangular. Using back-substitution, we .nd the solution: z = 2, y = 1, x = 1. 
Observe that we have performed only row operations. The general method is to iteratively eliminate variables using simple row operations (namely, adding or subtracting a multiple of a row to another row of the matrix) while simultaneously applying these operations to the vector b, to obtain a system, MAx = Mb, where MA is upper-triangular. Such a method is called Gaussian elimination. However, one extra twist is needed for the method to work in all cases: It may be necessary to permute rows, as illustrated by the following example: 
x + y + z =1 x + y +3z =1 2x +5y +8z =1. 
7.2. GAUSSIAN ELIMINATION 
In order to eliminate x from the second and third row, we subtract the .rst row from the second and we subtract twice the .rst row from the third: 
x + y + z =1 2z =0 3y +6z = .1. 
Now the trouble is that y does not occur in the second row; so, we can¡¯t eliminate y from the third row by adding or subtracting a multiple of the second row to it. The remedy is simple: Permute the second and the third row! We get the system: 
x + y + z =1 3y +6z = .1 2z =0, 
which is already in triangular form. Another example where some permutations are needed is: 
z =1 .2x +7y +2z =1 4x . 6y = .1. 
First we permute the .rst and the second row, obtaining 
.2x +7y +2z =1 z =1 4x . 6y = .1, 
and then we add twice the .rst row to the third, obtaining: 
.2x +7y +2z =1 z =1 8y +4z =1. 
Again we permute the second and the third row, getting 
.2x +7y +2z =1 8y +4z =1 z =1, 
an upper-triangular system. Of course, in this example, z is already solved and we could have eliminated it .rst, but for the general method, we need to proceed in a systematic fashion. 
We now describe the method of Gaussian elimination applied to a linear system Ax = b, where A is assumed to be invertible. We use the variable k to keep track of the stages of elimination. Initially, k = 1. 
(1) 
The .rst step is to pick some nonzero entry 	ai 1 in the .rst column of A. Such an entry must exist, since A is invertible (otherwise, the .rst column of A would be the zero vector, and the columns of A would not be linearly independent. Equivalently, we would have det(A) = 0). The actual choice of such an element has some impact on the numerical stability of the method, but this will be examined later. For the time being, we assume that some arbitrary choice is made. This chosen element is called the pivot of the elimination step and is denoted ¦Ð1 (so, in this .rst step, ¦Ð1 = ai 1). 

(2) 
Next we permute the row (i) corresponding to the pivot with the .rst row. Such a step is called pivoting. So after this permutation, the .rst element of the .rst row is nonzero. 

(3) 
We now eliminate the variable 	x1 from all rows except the .rst by adding suitable multiples of the .rst row to these rows. More precisely we add .ai 1/¦Ð1 times the .rst row to the ith row for i =2,...,n. At the end of this step, all entries in the .rst column are zero except the .rst. 

(4) 
Increment k by1. If k = n, stop. Otherwise, k<n, and then iteratively repeat Steps (1), (2), (3) on the (n . k + 1) ¡Á (n . k + 1) subsystem obtained by deleting the .rst k . 1 rows and k . 1 columns from the current system. 


If we let A1 = A and Ak =(aij (k)) be the matrix obtained after k . 1 elimination steps (2 ¡Ü k ¡Ü n), then the kth elimination step is applied to the matrix Ak of the form 
. 

a(k) 1 1  a(k) 1 2  ¡¤ ¡¤ ¡¤  ¡¤ ¡¤ ¡¤  ¡¤ ¡¤ ¡¤  a(k) 1 n  
0  a(k) 2 2  ¡¤ ¡¤ ¡¤  ¡¤ ¡¤ ¡¤  ¡¤ ¡¤ ¡¤  a(k) 2 n  
. . .  .. .  ...  . . .  . . .  
0  0  0  a(k) k k  ¡¤ ¡¤ ¡¤  a(k) k n  
.  .  .  .  .  
.  .  .  .  .  
.  .  .  .  .  
0  0  0  a(k) n k  ¡¤ ¡¤ ¡¤  a(k) n n  

. 

......... 

......... 

Ak = 
. 

Actually, note that 
(k)(i)
a= a
ij ij 
for all i, j with 1 ¡Ü i ¡Ü k . 2 and i ¡Ü j ¡Ü n, since the .rst k . 1 rows remain unchanged after the (k . 1)th step. 
We will prove later that det(Ak)= ¡À det(A). Consequently, Ak is invertible. The fact that Ak is invertible i. A is invertible can also be shown without determinants from the fact that there is some invertible matrix Mk such that Ak = MkA, as we will see shortly. 
Since Ak is invertible, some entry aik (k) with k ¡Ü i ¡Ü n is nonzero. Otherwise, the last n . k + 1 entries in the .rst k columns of Ak would be zero, and the .rst k columns of Ak would yield k vectors in Rk.1 . But then the .rst k columns of Ak would be linearly 

7.3. ELEMENTARY MATRICES AND ROW OPERATIONS 
dependent and Ak would not be invertible, a contradiction. This situation is illustrated by the following matrix for n = 5 and k = 3: 
.

.

(3) (3) (3) (3) (3)
aaaaa
11 12 13 13 15 
(3) (3) (3) (3)
0 aaaa
22 23 24 25 
...... 

...... 

(3) (3)
000 aa
34 35 
. 

(3) (3)
000 aa
44 4 n 
(3) (3)
000 aa
54 55 
The .rst three columns of the above matrix are linearly dependent. So one of the entries aik (k) with k ¡Ü i ¡Ü n can be chosen as pivot, and we permute the kth 
(k)(k)
row with the ith row, obtaining the matrix ¦Á(k) =(¦Á). The new pivot is ¦Ðk = ¦Á, and 
jl kk we zero the entries i = k +1,...,n in column k by adding .¦Á(i k k)/¦Ðk times row k to row i. At the end of this step, we have Ak+1. Observe that the .rst k . 1 rows of Ak are identical to the .rst k . 1 rows of Ak+1. 
The process of Gaussian elimination is illustrated in schematic form below: 
.
..
..
..
. 
¡Á¡Á¡Á¡Á ¡Á¡Á¡Á¡Á ¡Á¡Á¡Á¡Á ¡Á¡Á¡Á¡Á 

... 

¡Á¡Á¡Á¡Á 

¡Á¡Á¡Á¡Á 

... 

=. 

... 

0 

¡Á¡Á¡Á

... 

=. 

... 

0 

¡Á¡Á¡Á 

0 0 

¡Á¡Á 

... 

=. 

... 

0 

¡Á¡Á¡Á 

00 

¡Á¡Á 

...

. 

0 

¡Á¡Á¡Á 

¡Á¡Á¡Á¡Á 0 ¡Á¡Á¡Á 0 0 ¡Á¡Á 00 0 ¡Á 


7.3 Elementary Matrices and Row Operations 
It is easy to .gure out what kind of matrices perform the elementary row operations used during Gaussian elimination. The key point is that if A = PB, where A, B are m ¡Á n matrices and P is a square matrix of dimension m, if (as usual) we denote the rows of A and B by A1,...,Am and B1,...,Bm, then the formula 
 
m
aij = pikbkj k=1 
giving the (i, j)th entry in A shows that the ith row of A is a linear combination of the rows of B: 
Ai = pi1B1 + ¡¤¡¤¡¤ + pimBm. 
Therefore, multiplication of a matrix on the left by a square matrix performs row opera-tions. Similarly, multiplication of a matrix on the right by a square matrix performs column operations 
The permutation of the kth row with the ith row is achieved by multiplying A on the left by the transposition matrix P (i, k), which is the matrix obtained from the identity matrix 
.
. 
1  
1  
0  1  
1  
...  
1  

.............. 

.............. 

.

P (i, k)= 

1  0  
1  
1  
. . 

For example, if m = 3, 

. 

001 
010

.

P (1, 3) = 

, 

100 
then 
.
..
..
. 
001 b11 b12 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b1n b31 b32 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b3n 
P (1, 3)B =

.

010

.
.

b21 b22 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b2n 
.

=

.

b21 b22 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b2n 
.

. 

100 b31 b32 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b3n b11 b12 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b1n 
Observe that det(P (i, k)) = .1. Furthermore, P (i, k) is symmetric (P (i, k)T = P (i, k)), and 
P (i, k).1 = P (i, k). 
During the permutation Step (2), if row k and row i need to be permuted, the matrix A is multiplied on the left by the matrix Pk such that Pk = P (i, k), else we set Pk = I. Adding ¦Â times row j to row i (with i = j) is achieved by multiplying A on the left by the elementary matrix, Ei,j;¦Â = I + ¦Âeij, 
where

 

1 if k = i and l = j
(eij)kl =
0 if k = i or l = j, 
i.e., 
.
.
.
. 
11 

Ei,j;¦Â 
= 

.............. 

1 

1 

.............. 

or 

Ei,j;¦Â 
= 

.............. 

1 

1 

.............. 

, 

1 

1 ¦Â 

1

1 

.

.

.

. 

.

. 

1 

¦Â 1 

1 

1 

11 

7.3. ELEMENTARY MATRICES AND ROW OPERATIONS 
on the left, i>j, and on the right, i<j. The index i is the index of the row that is changed by the multiplication. For example, if m = 3 and we want to add twice row 1 to row 3, since ¦Â = 2, j = 1 and i = 3, we form 
. .. .. . 
100 000 100 E3,1;2 = I +2e31 = .010. + .000. = .010., 001 200 201 
and calculate 
... . 
100 b11 b12 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b1n 
E3,1;2B = .010..b21 b22 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b2n . 
201 b31 b32 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b3n
.. 
b11 b12 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b1n 
..
= b21 b22 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ b2n . 2b11 + b31 2b12 + b32 ¡¤¡¤¡¤ ¡¤¡¤¡¤ ¡¤¡¤¡¤ 2b1n + b3n 
Observe that the inverse of Ei,j;¦Â = I + ¦Âeij is Ei,j;.¦Â = I . ¦Âeij and that det(Ei,j;¦Â) = 1. Therefore, during Step 3 (the elimination step), the matrix A is multiplied on the left by a product Ek of matrices of the form Ei,k;¦Âi,k , with i>k. 
Consequently, we see that 
Ak+1 = EkPkAk, 
and then 
Ak = Ek.1Pk.1 ¡¤¡¤¡¤ E1P1A. 
This justi.es the claim made earlier that Ak = MkA for some invertible matrix Mk; we can pick 
Mk = Ek.1Pk.1 ¡¤¡¤¡¤ E1P1, 
a product of invertible matrices. 
The fact that det(P (i, k)) = .1 and that det(Ei,j;¦Â) = 1 implies immediately the fact claimed above: We always have 
det(Ak)= ¡À det(A). 
Furthermore, since 
Ak = Ek.1Pk.1 ¡¤¡¤¡¤ E1P1A 
and since Gaussian elimination stops for k = n, the matrix 
An = En.1Pn.1 ¡¤¡¤¡¤ E2P2E1P1A 
is upper-triangular. Also note that if we let M = En.1Pn.1 ¡¤¡¤¡¤ E2P2E1P1, then det(M)= ¡À1, and 
det(A)= ¡À det(An). 
The matrices P (i, k) and Ei,j;¦Â are called elementary matrices. We can summarize the above discussion in the following theorem: 
Theorem 7.1. (Gaussian elimination) Let A be an n ¡Á n matrix (invertible or not). Then there is some invertible matrix M so that U = MA is upper-triangular. The pivots are all nonzero i. A is invertible. 
Proof. We already proved the theorem when A is invertible, as well as the last assertion. Now A is singular i. some pivot is zero, say at Stage k of the elimination. If so, we must have a(i k k) = 0 for i = k, . . . , n; but in this case, Ak+1 = Ak and we may pick Pk = Ek = I. 
Remark: Obviously, the matrix M can be computed as 
M = En.1Pn.1 ¡¤¡¤¡¤ E2P2E1P1, 
but this expression is of no use. Indeed, what we need is M.1; when no permutations are needed, it turns out that M.1 can be obtained immediately from the matrices Ek¡¯s, in fact, from their inverses, and no multiplications are necessary. 
Remark: Instead of looking for an invertible matrix M so that MA is upper-triangular, we can look for an invertible matrix M so that MA is a diagonal matrix. Only a simple change to Gaussian elimination is needed. At every Stage k, after the pivot has been found and pivoting been performed, if necessary, in addition to adding suitable multiples of the kth row to the rows below row k in order to zero the entries in column k for i = k +1,...,n, also add suitable multiples of the kth row to the rows above row k in order to zero the entries in column k for i =1,...,k . 1. Such steps are also achieved by multiplying on the left by elementary matrices Ei,k;¦Âi,k , except that i<k, so that these matrices are not lower-triangular matrices. Nevertheless, at the end of the process, we .nd that An = MA, is a diagonal matrix. 
This method is called the Gauss-Jordan factorization. Because it is more expensive than Gaussian elimination, this method is not used much in practice. However, Gauss-Jordan factorization can be used to compute the inverse of a matrix A. Indeed, we .nd the jth column of A.1 by solving the system Ax(j) = ej (where ej is the jth canonical basis vector of Rn). By applying Gauss-Jordan, we are led to a system of the form Djx(j) = Mjej, where Dj is a diagonal matrix, and we can immediately compute x(j). 
It remains to discuss the choice of the pivot, and also conditions that guarantee that no permutations are needed during the Gaussian elimination process. We begin by stating a necessary and su.cient condition for an invertible matrix to have an LU-factorization (i.e., Gaussian elimination does not require pivoting). 


7.4 LU-Factorization 
De.nition 7.1. We say that an invertible matrix A has an LU-factorization if it can be written as A = LU, where U is upper-triangular invertible and L is lower-triangular, with Lii = 1 for i =1,...,n. 
7.4. LU-FACTORIZATION 
A lower-triangular matrix with diagonal entries equal to 1 is called a unit lower-triangular matrix. Given an n ¡Á n matrix A =(aij), for any k with 1 ¡Ü k ¡Ü n, let A(1 : k, 1: k) denote the submatrix of A whose entries are aij, where 1 ¡Ü i, j ¡Ü k. 1 For example, if A is the 5 ¡Á 5 
matrix

.
. 
a11 a12 a13 a14 a15 a21 a22 a23 a24 a25 
a31 a32 a33 a34 a35 
a41 a42 a43 a44 a45 a51 a52 a53 a54 a55 
..... 

..... 

A = 

, 

then

.
. 
a11 a12 a13 
.

.

A(1 : 3, 1 : 3) = 

a21 a22 a23 
. 

a31 a32 a33 Proposition 7.2. Let A be an invertible n ¡Á n-matrix. Then A has an LU-factorization A = LU i. every matrix A(1 : k, 1: k) is invertible for k =1,...,n. Furthermore, when A has an LU-factorization, we have det(A(1 : k, 1: k)) = ¦Ð1¦Ð2 ¡¤¡¤¡¤ ¦Ðk,k =1, . . . , n, where ¦Ðk is the pivot obtained after k . 1 elimination steps. Therefore, the kth pivot is given 
by 

¦Ðk 
= 

.. . 

a11 = det(A(1 : 1, 1 : 1)) if k =1 
det(A(1 : k, 1: k)) 

if k =2, . . . , n. 

det(A(1 : k . 1, 1: k . 1)) 

Proof. First assume that A = LU is an LU-factorization of A. We can write 

  
 
  

 
  
A(1 : k, 1: k) A2 L1 0 U1 U2 L1U1 L1U2
A ===,
A3 A4L3 L40 U4L3U1 L3U2 + L4U4
where L1,L4 are unit lower-triangular and U1,U4 are upper-triangular. (Note, A(1 : k, 1: k), L1, and U1 are k ¡Á k matrices; A2 and U2 are k ¡Á (n . k) matrices; A3 and L3 are (n . k) ¡Á k matrices; A4, L4, and U4 are (n . k) ¡Á (n . k) matrices.) Thus, 
A(1 : k, 1: k)= L1U1, 
and since U is invertible, U1 is also invertible (the determinant of U is the product of the diagonal entries in U, which is the product of the diagonal entries in U1 and U4). As L1 is invertible (since its diagonal entries are equal to 1), we see that A(1 : k, 1: k) is invertible for k =1,...,n. 
Conversely, assume that A(1 : k, 1: k) is invertible for k =1,...,n. We just need to show that Gaussian elimination does not need pivoting. We prove by induction on k that the kth step does not need pivoting. 
1We are using Matlab¡¯s notation. 
This holds for k = 1, since A(1 :1, 1 :1) = (a11), so a11 = 0. Assume that no pivoting was necessary for the .rst k . 1 steps (2 ¡Ü k ¡Ü n . 1). In this case, we have 
Ek.1 ¡¤¡¤¡¤ E2E1A = Ak, 
where L = Ek.1 ¡¤¡¤¡¤ E2E1 is a unit lower-triangular matrix and Ak(1 : k, 1: k) is upper-triangular, so that LA = Ak can be written as 
L1 0 A(1 : k, 1: k) A2 U1 B2 
= ,
L3 L4 A3 A4 0 B4 
where L1 is unit lower-triangular and U1 is upper-triangular. (Once again A(1 : k, 1: k), L1, and U1 are k ¡Á k matrices; A2 and B2 are k ¡Á (n . k) matrices; A3 and L3 are (n . k) ¡Á k matrices; A4, L4, and B4 are (n . k) ¡Á (n . k) matrices.) But then, 
L1A(1 : k, 1: k)) = U1, 
where L1 is invertible (in fact, det(L1) = 1), and since by hypothesis A(1 : k, 1: k) is invertible, U1 is also invertible, which implies that (U1)kk = 0, since U1 is upper-triangular. Therefore, no pivoting is needed in Step k, establishing the induction step. Since det(L1) = 1, we also have 
det(U1) = det(L1A(1 : k, 1: k)) = det(L1) det(A(1 : k, 1: k)) = det(A(1 : k, 1: k)), 
and since U1 is upper-triangular and has the pivots ¦Ð1,...,¦Ðk on its diagonal, we get 
det(A(1 : k, 1: k)) = ¦Ð1¦Ð2 ¡¤¡¤¡¤ ¦Ðk,k =1, . . . , n, 
as claimed. 
Remark: The use of determinants in the .rst part of the proof of Proposition 7.2 can be avoided if we use the fact that a triangular matrix is invertible i. all its diagonal entries are nonzero. 
Corollary 7.3. (LU-Factorization) Let A be an invertible n ¡Á n-matrix. If every matrix A(1 : k, 1: k) is invertible for k =1,...,n, then Gaussian elimination requires no pivoting and yields an LU-factorization A = LU. 
Proof. We proved in Proposition 7.2 that in this case Gaussian elimination requires no pivoting. Then since every elementary matrix Ei,k;¦Â is lower-triangular (since we always arrange that the pivot ¦Ðk occurs above the rows that it operates on), since E.1 = 
i,k;¦Â Ei,k;.¦Â and the Eks are products of Ei,k;¦Âi,k s, from 
En.1 ¡¤¡¤¡¤ E2E1A = U, 

7.4. LU-FACTORIZATION 
where U is an upper-triangular matrix, we get A = LU, 
E.1E.1 E.1
where L = ¡¤¡¤¡¤ is a lower-triangular matrix. Furthermore, as the diagonal 
12 n.1 
entries of each Ei,k;¦Â are 1, the diagonal entries of each Ek are also 1. Example 7.1. The reader should verify that 
.
..
..
. 
2110 1000 2110 0
... 
... 

... 

... 

... 

... 

4331 

210 

0111 

= 

0 6798 3411 0002 
8795 

431 

0022 

is an LU-factorization. 

One of the main reasons why the existence of an LU-factorization for a matrix A is interesting is that if we need to solve several linear systems Ax = b corresponding to the same matrix A, we can do this cheaply by solving the two triangular systems 
Lw = b, and Ux = w. 
There is a certain asymmetry in the LU-decomposition A = LU of an invertible matrix A. Indeed, the diagonal entries of L are all 1, but this is generally false for U. This asymmetry can be eliminated as follows: if 
D = diag(u11,u22,...,unn) 
is the diagonal matrix consisting of the diagonal entries in U (the pivots), then we if let UÊ¬ 
= D.1U, we can write A = LDUÊ¬, 
where L is lower-triangular, UÊ¬ is upper-triangular, all diagonal entries of both L and UÊ¬ 
are 1, and D is a diagonal matrix of pivots. Such a decomposition leads to the following de.nition. 
De.nition 7.2. We say that an invertible n¡Án matrix A has an LDU-factorization if it can be written as A = LDUÊ¬, where L is lower-triangular, UÊ¬ is upper-triangular, all diagonal 
entries of both L and UÊ¬ are 1, and D is a diagonal matrix. 
We will see shortly than if A is real symmetric, then UÊ¬ = LT . 
As we will see a bit later, real symmetric positive de.nite matrices satisfy the condition of Proposition 7.2. Therefore, linear systems involving real symmetric positive de.nite matrices can be solved by Gaussian elimination without pivoting. Actually, it is possible to do better: this is the Cholesky factorization. 
If a square invertible matrix A has an LU-factorization, then it is possible to .nd L and U while performing Gaussian elimination. Recall that at Step k, we pick a pivot ¦Ðk = aik (k) =0 in the portion consisting of the entries of index j ¡Ý k of the k-th column of the matrix Ak obtained so far, we swap rows i and k if necessary (the pivoting step), and then we zero the entries of index j = k +1,...,n in column k. Schematically, we have the following steps: 
.
..
..
. 
¡Á ¡Á ¡Á¡Á¡Á ¡Á ¡Á ¡Á¡Á¡Á ¡Á¡Á¡Á¡Á¡Á 
¡Á 
¡Á 
..... 
..... 

..... 

..... 

..... 

..... 

(k)
0 ¡Á ¡Á¡Á 

0

0 

¡Á¡Á¡Á¡Á

¡Á¡Á¡Á

a

pivot 
. 
ik 
elim
.
0 

¡Á ¡Á¡Á 

0 0 ¡Á¡Á¡Á

0 ¡Á ¡Á¡Á¡Á

=

=

. 

(k) 
0 0 ¡Á¡Á¡Á

0 

0

¡Á 
0 ¡Á ¡Á¡Á¡Á 0 ¡Á ¡Á¡Á¡Á 0 0 ¡Á¡Á¡Á 
More precisely, after permuting row k and row i (the pivoting step), if the entries in column k below row k are ¦Ák+1k,...,¦Ánk, then we add .¦Ájk/¦Ðk times row k to row j; this process is illustrated below: 
¡Á¡Á 

¡Á ¡Á¡Á¡Á

a

ik 
.

.
.

..
.
..
(k) (k) 
¦Ðk ¦Ðk
a

a

kk ik 
......... 

......... 

......... 

......... 

= 

........ 

(k) (k) 
¦Ák+1k 
. 

. 

. 

¦Áik 
. 

. 

........ 

........ 

0 

. 

. 

. 

0 

. 

. 

........ 

a

a

k+1k k+1k 
. 

. 

. 

(k) 
. 

. 

. 

(k) 
pivot 
. 
elim
.
=

=

. 

a

a

ik kk 
. 

. 

. 

. 

..

. 

. 

(k)(k)
aa¦Ánk 0 
nk nk 
Then if we write ¦Éjk = ¦Ájk/¦Ðk for j = k +1,...,n, the kth column of L is 
.
. 

.......... 

0 . 
. 
. 0 1 
¦Ék+1k 
. 
. 
. 
¦Énk 
.......... 

. 

Observe that the signs of the multipliers .¦Ájk/¦Ðk have been .ipped. Thus, we obtain the unit lower triangular matrix 
.
. 
L = 

...... 

100 ¡¤¡¤¡¤ 0 ¦É21 10 ¡¤¡¤¡¤ 0 ¦É31 ¦É32 1 ¡¤¡¤¡¤ 0 
... 
.
....
... .0 
¦Én1 ¦Én2 ¦Én3 ¡¤¡¤¡¤ 1 
...... 

. 

It is easy to see (and this is proven in Theorem 7.5) that the inverse of L is obtained from 

7.4. LU-FACTORIZATION L by .ipping the signs of the ¦Éij: 
.
. 
L.1 
= 

...... 

100 ¡¤¡¤¡¤ 0 .¦É21 10 ¡¤¡¤¡¤ 0 
.¦É31 .¦É32 1 ¡¤¡¤¡¤ 0 
... 
.
. ...
. .. .0 
.¦Én1 .¦Én2 .¦Én3 ¡¤¡¤¡¤ 1 
...... 

. 

Furthermore, if the result of Gaussian elimination (without pivoting) is U = En.1 ¡¤¡¤¡¤ E1A, then 
.
.
.
. 
1 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 01 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 

Ek 
= 

........ 

........ 

and 

E.1 
k 
= 

........

........ 

, 

. . ...
.
. . ...
.
. 

.. ... 

. 

.. ...

. 

.

. 

. ... 

. 

. ... 

0 ¡¤¡¤¡¤ 10 ¡¤¡¤¡¤ 0 

0 ¡¤¡¤¡¤ .¦Ék+1k 1 ¡¤¡¤¡¤ 0 
0 ¡¤¡¤¡¤ 10 ¡¤¡¤¡¤ 0 

0 ¡¤¡¤¡¤ ¦Ék+1k 1 ¡¤¡¤¡¤ 0 
.... .
.
.... .
.
.. .. 

.. 

.. .. 

..

. 

.

.. .. 

. .. .. 

. 

0 ¡¤¡¤¡¤ .¦Énk 0 ¡¤¡¤¡¤ 10 ¡¤¡¤¡¤ ¦Énk 0 ¡¤¡¤¡¤ 1 
so the kth column of Ek is the kth column of L.1 . Here is an example illustrating the method. 
Example 7.2. Given 

.
. 
... 

11 1 0 
1 .10 1 

11 .10 
1 .10 .1 

...

A = A1 = 
, 

we have the following sequence of steps: The .rst pivot is ¦Ð1 = 1 in row 1, and we substract row 1 from rows 2, 3, and 4. We get 
.
.
.
. 
1110 1000 

A2 
= 

... 

0 .2 .11 

00 .20 

...

L1 
= 

... 

1100 

1010 

...

. 

0 .2 .1 .1 1001 
The next pivot is ¦Ð2 = .2 in row 2, and we subtract row 2 from row 4 (and add 0 times row 2 to row 3). We get 
.
.
.
. 
1110 1000 

A3 
= 

... 

0 .2 .11 

00 .20 

...

L2 
= 

... 

1100 

1010 

...

. 

00 0 .2 1101 

The next pivot is ¦Ð3 = .2 in row 3, and since the fourth entry in column 3 is already a zero, we add 0 times row 3 to row 4. We get 
.
.
.
. 
1110 1000 

A4 
= 

... 

0 .2 .11 

00 .20 

...

L3 
= 

... 

1100 

1010 

...

. 

00 0 .2 1101 The procedure is .nished, and we have 
.
.
.
. 
1000 1110 0
...
... 

... 

...

.2 .1

110 

0 

1 

L = L3 U = A4 
= 

= 

.

.2

0 1101 00 0 .2 
It is easy to check that indeed 
101 

00 

0 

.
..
..
. 
1000 1110 1110 

... 

1100 

1010 

... 
... 

... 

= 

... 

... 

.2 .1 

1 

.1

0 

1 

01 

LU 

= 

= A.

.2 

.1 

0

00 

0 

11 

1101 00 0 .21 .10 .1 
We now show how to extend the above method to deal with pivoting e.ciently. This is the PA = LU factorization. 
7.5 PA = LU Factorization 
The following easy proposition shows that, in principle, A can be premultiplied by some permutation matrix P , so that PA can be converted to upper-triangular form without using any pivoting. Permutations are discussed in some detail in Section 6.1, but for now we just need this de.nition. For the precise connection between the notion of permutation (as discussed in Section 6.1) and permutation matrices, see Problem 7.16. 
De.nition 7.3. A permutation matrix is a square matrix that has a single 1 in every row and every column and zeros everywhere else. 
It is shown in Section 6.1 that every permutation matrix is a product of transposition matrices (the P (i, k)s), and that P is invertible with inverse P T . 
Proposition 7.4. Let A be an invertible n ¡Á n-matrix. There is some permutation matrix P so that (PA)(1 : k, 1: k) is invertible for k =1,...,n. 
Proof. The case n = 1 is trivial, and so is the case n = 2 (we swap the rows if necessary). If n ¡Ý 3, we proceed by induction. Since A is invertible, its columns are linearly independent; in particular, its .rst n . 1 columns are also linearly independent. Delete the last column of 
7.5. PA = LU FACTORIZATION 
A. Since the remaining n . 1 columns are linearly independent, there are also n . 1 linearly independent rows in the corresponding n ¡Á (n . 1) matrix. Thus, there is a permutation of these n rows so that the (n . 1) ¡Á (n . 1) matrix consisting of the .rst n . 1 rows is invertible. But then there is a corresponding permutation matrix P1, so that the .rst n . 1 
rows and columns of P1A form an invertible matrix AÊ¬ . Applying the induction hypothesis to the (n . 1) ¡Á (n . 1) matrix AÊ¬, we see that there some permutation matrix P2 (leaving 
the nth row .xed), so that (P2P1A)(1 : k, 1: k) is invertible, for k =1,...,n . 1. Since A is invertible in the .rst place and P1 and P2 are invertible, P1P2A is also invertible, and we are done. 
Remark: One can also prove Proposition 7.4 using a clever reordering of the Gaussian elimination steps suggested by Trefethen and Bau [68] (Lecture 21). Indeed, we know that if A is invertible, then there are permutation matrices Pi and products of elementary matrices Ei, so that 
An = En.1Pn.1 ¡¤¡¤¡¤ E2P2E1P1A, 
where U = An is upper-triangular. For example, when n = 4, we have E3P3E2P2E1P1A = U. We can de.ne new matrices EÊ¬ ,EÊ¬ ,EÊ¬ which are still products of elementary matrices so 
123 
that we have EÊ¬ EÊ¬ EÊ¬ P3P2P1A = U. 
321Indeed, if we let EÊ¬ = E3, EÊ¬ = P3E2P .1, and E1 Ê¬ = P3P2E1P2 .1P3 .1, we easily verify that each EÊ¬ isaproduc3tofelemen2tarymatr3ices and that 
k 
E3Ê¬ E2Ê¬ E1Ê¬ P3P2P1 = E3(P3E2P3 .1)(P3P2E1P .1P3 .1)P3P2P1
2 
= E3P3E2P2E1P1. 
It can also be proven that E1Ê¬ ,E2Ê¬ ,E3 Ê¬ are lower triangular (see Theorem 7.5). 
In general, we let 
Ek Ê¬ = Pn.1 ¡¤¡¤¡¤ Pk+1EkPk.+11 ¡¤¡¤¡¤ Pn..11, 
and we have EnÊ¬.1 ¡¤¡¤¡¤ E1Ê¬ Pn.1 ¡¤¡¤¡¤ P1A = U, 
where each Ej Ê¬ is a lower triangular matrix (see Theorem 7.5). 
It is remarkable that if pivoting steps are necessary during Gaussian elimination, a very simple modi.cation of the algorithm for .nding an LU-factorization yields the matrices L, U, and P , such that PA = LU. To describe this new method, since the diagonal entries of L are 1s, it is convenient to write 
L = I +¦«. 
Then in assembling the matrix ¦« while performing Gaussian elimination with pivoting, we make the same transposition on the rows of ¦« (really ¦«k.1) that we make on the rows of A (really Ak) during a pivoting step involving row k and row i. We also assemble P by starting with the identity matrix and applying to P the same row transpositions that we apply to A and ¦«. Here is an example illustrating this method. 
Example 7.3. Given 

.
. 
11 1 0 
11 .10 

A = A1 = 
... 

...

,

1 .10 1 
1 .10 .1 

we have the following sequence of steps: We initialize ¦«0 = 0 and P0 = I4. The .rst pivot is ¦Ð1 = 1 in row 1, and we subtract row 1 from rows 2, 3, and 4. We get 
.
.
.
.
.
. 
1 1 1 0 0000 1000 

A2 
= 

... 

00 .20 

0 .2 .11 

...

¦«1 
= 

... 

1000 

1000 

...

P1 
= 

... 

0100 

0010 

...

. 

0 .2 .1 .1 1000 0001 
The next pivot is ¦Ð2 = .2 in row 3, so we permute row 2 and 3; we also apply this permutation to ¦« and P : 
.
.
.
.
.
. 
1 1 1 0 0000 1000 

A

Ê¬ 
3 
= 

... 

0 .2 .11 

00 .20 

...

¦«

Ê¬ 
2 
= 

... 

1000 

1000 

...

P2 
= 

... 

0010 

0100 

...

. 

0 .2 .1 .1 1000 0001 Next we subtract row 2 from row 4 (and add 0 times row 2 to row 3). We get 
.
.
.
.
.
. 
1 1 1 0 0000 1000 

A3 
= 

... 

0 .2 .11 

00 .20 

...

¦«2 
= 

... 

1000 

1000 

...

P2 
= 

... 

0010 

0100 

...

. 

00 0 .2 1100 0001 
The next pivot is ¦Ð3 = .2 in row 3, and since the fourth entry in column 3 is already a zero, we add 0 times row 3 to row 4. We get 
.
.
.
.
.
. 
1 1 1 0 0000 1000 

A4 
= 

... 

0 .2 .11 

00 .20 

...

¦«3 
= 

... 

1000 

1000 

...

P3 
= 

... 

0010 

0100 

...

. 

00 0 .2 1100 0001 The procedure is .nished, and we have 
.
.
.
. 
1000 1110 

L =¦«3 + I 
= 

... 

1100 

1010 

...

U = A4 
= 

... 

... 

0 

.2 .1 

1 

00 

.2 

0 

1101 00 0 .2 

.
. 
... 

1000 
0010 

0100 
0001 

...

P = P3 = 
. 

7.5. PA = LU FACTORIZATION It is easy to check that indeed 
.
..
..
. 
1000 1110 1110 

... 

1100 

1010 

... 
... 

... 

= 

... 

... 

0 

.2 .1 

1 

1 

.1 

01 

LU = 

00 

.2 

0 

11 

.1 

0 

1101 00 0 .21 .10 .1 

and

.
..
..
. 
1000 1110 1110 

... 

0010 

0100 

... 
... 

... 

= 

... 

...

11 

.1 

0 

1 

.1 

01 

PA = 

. 

1 

.1 

01 

11 

.1 

0 

0001 1 .10 .11 .10 .1 
Using the idea in the remark before the above example, we can prove the theorem below which shows the correctness of the algorithm for computing P, L and U using a simple adaptation of Gaussian elimination. 
We are not aware of a detailed proof of Theorem 7.5 in the standard texts. Although Golub and Van Loan [30] state a version of this theorem as their Theorem 3.1.4, they say that ¡°The proof is a messy subscripting argument.¡± Meyer [48] also provides a sketch of proof (see the end of Section 3.10). In view of this situation, we o.er a complete proof. It does involve a lot of subscripts and superscripts, but in our opinion, it contains some techniques that go far beyond symbol manipulation. 
Theorem 7.5. For every invertible n ¡Á n-matrix A, the following hold: 
(1) 
There is some permutation matrix P , some upper-triangular matrix U, and some unit lower-triangular matrix L, so that PA = LU (recall, Lii =1 for i =1,...,n). Fur-thermore, if P = I, then L and U are unique and they are produced as a result of Gaussian elimination without pivoting. 

(2) 
If En.1 ...E1A = U is the result of Gaussian elimination without pivoting, write as 


(k)(k)(k)
usual Ak = Ek.1 ...E1A (with Ak =(a)), and let ¦Éik = a/a, with 1 ¡Ü k ¡Ü n . 1
ij ikkk 
and k +1 ¡Ü i ¡Ü n. Then 
.
. 
L = 

...... 

100 ¡¤¡¤¡¤ 0 ¦É21 10 ¡¤¡¤¡¤ 0 ¦É31 ¦É32 1 ¡¤¡¤¡¤ 0 
... 
.
....
... . 0 
¦Én1 ¦Én2 ¦Én3 ¡¤¡¤¡¤ 1 
...... 

, 

where the kth column of L is the kth column of Ek .1, for k =1,...,n . 1. 
(3) If En.1Pn.1 ¡¤¡¤¡¤ E1P1A = U is the result of Gaussian elimination with some pivoting, write Ak = Ek.1Pk.1 ¡¤¡¤¡¤ E1P1A, and de.ne Ejk, with 1 ¡Ü j ¡Ü n . 1 and j ¡Ü k ¡Ü n . 1, such that, for j =1,...,n . 2, 
Ej 
j = Ej Ejk = PkEjk.1Pk, for k = j +1,...,n . 1, 
and 
En.1 
n.1 = En.1. 
Then, 
Ek = PkPk.1 ¡¤¡¤¡¤ Pj+1EjPj+1 ¡¤¡¤¡¤ Pk.1Pk
j = En.1 En.1
U ¡¤¡¤¡¤ Pn.1 ¡¤¡¤¡¤ P1A,
n.11 
and if we set 
P = Pn.1 ¡¤¡¤¡¤ P1 ).1 (En.1
L =(E1 n.1 ¡¤¡¤¡¤ n.1 ).1 , 
then 
PA = LU. (.1) 
Furthermore, 
(Ejk).1 = I + Ejk , 1 ¡Ü j ¡Ü n . 1,j ¡Ü k ¡Ü n . 1, where Ejk is a lower triangular matrix of the form 
.
. 
Ek 
= 
j 
......... 

0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 
. . ...
.
......
.
. . ... 
0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 
0 ¡¤¡¤¡¤ ¦É

(k) 
j+1j 
0 ¡¤¡¤¡¤ 0 

.... .
.
.. ....
.
.... . 
(k)
0 ¡¤¡¤¡¤ ¦É0 ¡¤¡¤¡¤ 0
nj 
......... 

, 

we have  
Ek j = I . Ek j ,  
and  
Ek j = PkEk.1 j ,  1 ¡Ü j ¡Ü n . 2, j + 1 ¡Ü k ¡Ü n . 1,  

where Pk = I or else Pk = P (k, i) for some i such that k +1 ¡Ü i ¡Ü n; if Pk = I, this means that (Ejk).1 is obtained from (Ejk.1).1 by permuting the entries on rows i and k in column j. Because the matrices (Ejk).1 are all lower triangular, the matrix L is also lower triangular. 
7.5. PA = LU FACTORIZATION In order to .nd L, de.ne lower triangular n ¡Á n matrices ¦«k of the form 
. 

0 0000 ¡¤¡¤¡¤ ¡¤¡¤¡¤ 	0 0 
. ............... 
(k) 	.. .. 
¦Ë0 0 00 .
21 
............... 
. 

(k)(k) ..	....
¦Ë¦Ë. 00 . 
0 
. . ..	.
.
. ....	. 
... 00 . .. 0 
.

31 32 
¦«k = 
(k)(k)(k)
¦Ë¦Ë¡¤¡¤¡¤ ¦Ë0 
¡¤¡¤¡¤ ¡¤¡¤¡¤ 

k+11 k+12 
k+1k 
(k)(k)(k) ..
¦Ë¦Ë¡¤¡¤¡¤ ¦Ë0 
0 
.. ....
.	.
. . ......
.	.
.. .... 
. ¡¤¡¤¡¤ 

k+21 k+22 
k+2k 
(k)(k)(k)
¦Ë¦Ë¡¤¡¤¡¤ ¦Ë0 ¡¤¡¤¡¤ ¡¤¡¤¡¤ 0
n1 n2 nk 
to assemble the columns of L iteratively as follows: let (.¦Ék+1k,..., .¦Énk )
(k)(k) 
be the last n . k elements of the kth column of Ek, and de.ne ¦«k inductively by setting 
.
. .
¦«1 = .... 0 ¦É(1) 21 . . . ¦É(1) n1  0 0 . . . 0  ¡¤ ¡¤ ¡¤ ¡¤ ¡¤ ¡¤ ... ¡¤ ¡¤ ¡¤  0 0 . . . 0 ....,  
then for k = 2, . . . , n . 1, de.ne  
¦«Ê¬ k = Pk¦«k.1,  (.2)  
and ¦«k = (I + ¦«Ê¬ k)E.1 k . I, with  

. ............... 
............... 
0 000 0 ¡¤¡¤¡¤ ¡¤¡¤¡¤ 0 
£§(k.1) 	.. ..
¦Ë000 0 .. 0
21 £§(k.1) £§(k.1) ......¦Ë¦Ë. 00 .. 0
31 32 
. . 	...
.
. ..... ¦«k = £§(k. .1) £§(k. .1) . £§(k0 .1) 0 . .. ,
¦Ë¦Ë¡¤¡¤¡¤ ¦Ë0 ¡¤¡¤¡¤ ¡¤¡¤¡¤ 0
k1 k2 k (k.1) ¦Ë£§(k.1) ¦Ë£§(k.1) ¡¤¡¤¡¤ ¦Ë£§(k.1) ¦É(k) ... ¡¤¡¤¡¤ 0
k+11 k+12 k+1 (k.1) k+1k 
.. ..	..
.	.
. ... .... £§ . £§ .. £§ . .... 
(k.1) (k.1) (k.1) (k)
¦Ë¦Ë¡¤¡¤¡¤ ¦Ë¦É¡¤¡¤¡¤ ¡¤¡¤¡¤ 0
n1 n2 nk.1 nk 
with Pk = I or Pk = P (k, i) for some i>k. This means that in assembling L, row k and row i of ¦«k.1 need to be permuted when a pivoting step permuting row k and row i of Ak is required. Then 
I +¦«k =(E1 k).1 ¡¤¡¤¡¤ (Ekk).1 
¦«k = E1 k + ¡¤¡¤¡¤ + Ekk , 
for k =1,...,n . 1, and therefore 
L = I +¦«n.1. 
The proof of Theorem 7.5, which is very technical, is given in Section 7.6. 
We emphasize again that Part (3) of Theorem 7.5 shows the remarkable fact that in assembling the matrix L while performing Gaussian elimination with pivoting, the only change to the algorithm is to make the same transposition on the rows of ¦«k.1 that we make on the rows of A (really Ak) during a pivoting step involving row k and row i. We can also assemble P by starting with the identity matrix and applying to P the same row transpositions that we apply to A and ¦«. Here is an example illustrating this method. 
Example 7.4. Consider the matrix 
.
. 
1  2  .3  4  
A = ... 4 2  8 3  12 2  .8 1 ....  
.3  .1  1  .4  

We set P0 = I4, and we can also set ¦«0 = 0. The .rst step is to permute row 1 and row 2, using the pivot 4. We also apply this permutation to P0: 
.
.
.
. 
4  8  12  .8  0  1  0  0  
AÊ¬ 1 = ... 1 2  2 3  .3 2  4 1 ... P1 = ... 1 0  0 0  0 1  0 0 ....  
.3  .1  1  .4  0  0  0  1  

Next we subtract 1/4 times row 1 from row 2, 1/2 times row 1 from row 3, and add 3/4 times row 1 to row 4, and start assembling ¦«: 
.
.
.
.
.
. 
...
... 
...
... 
...
48 12 .8 0 000 0100 
..00661/40001000..
A2 =¦«1 = P1 = . 
0 .1 .45 1/2 000 0010 05 10 .10 .3/4000 0001 
Next we permute row 2 and row 4, using the pivot 5. We also apply this permutation to ¦« and P : 
. . . . . .  
4  8  12  .8  0  0  0  0  0  1  0  0  
AÊ¬ 3 = ... 0 0  5 .1  10 .4  .10 5 ... ¦«Ê¬ 2 = ... .3/4 1/2  0 0  0 0  0 0 ... P2 = ... 0 0  0 0  0 1  1 0 ....  
0  0  .6  6  1/4  0  0  0  1  0  0  0  
Next we add 1/5 times row 2 to row 3, and update ¦«Ê¬ 2:  

.
.
..
.
. 
...
48 12 .8 0 000 0100 0
... 
...
... 
05 10 .10 
¦«2 = 
00 .23 
...
... 
0001 0010 
.3/4 

00 

A3 
= 

P2 
= 

. 

1/2 

.1/5 

0 

0 00 .66 1/4 0 00 1000 
7.5. PA = LU FACTORIZATION 
Next we permute row 3 and row 4, using the pivot .6. We also apply this permutation to ¦« and P : 
.
.
.
.
.
. 
48 12 .8 0 000 0100 0
...
A

Ê¬ 
4 
= 

... 

...

... 

... 

...

. 

.10 

.3/4

05 10 

00 

0001

Ê¬
¦«

P3 
= 

= 

3
00 .66 

0 00 .23 1/2 .1/500 0010 
1/4 00 

1000 

Finally we subtract 1/3 times row 3 from row 4, and update ¦«

.
. 
3:
Ê¬ 
.
.
.
. ...
48 12 .8 0 000 0100 1
... 

...

... 

...

... 

.10 

.3/4

05 10 

0 00 

000 

A4 
¦«3 
P3 
= 

= 

= 

. 

00 .66 

0 0001 1/2 .1/51/30 0010 
Consequently, adding the identity to ¦«3, we obtain 
1/40 00 

100 

.
.
.
.
.
. 
1 0 00 4812 .8 0100 

... 

.3/4 

1 00 

1/40 10 

...

, 

U 

= 

... 

...

, 

P 

= 

... 

0001 

1000 

...

.10

05 10 

L = 

. 

00 

.6 

6 

1/2  .1/5  1/3  1  0  0  0  1  0  0  1  0  
We check that  
. .. . . .  

0100 12 .34 4812 .8 

PA = 

... 

0001 

1000 

... 
... 

4 812 

.8 

2321 

... 

= 

... 

.3 .11 .4 

12 

.3 

4 

...

, 

0010 .3 .11 .4 2321 
and that 
.
..
..
. 
1 0 00 4812 .8 4 812 .8 

LU 

= 

... 

.3/4 

1 00 

1/40 10 

... 
... 

... 

= 

... 

.3 .11 .4 

12 

.3 

4 

... 

= P A. 

.10

05 10 

00 

.6 

6 

1/2 .1/51/31 000 1 2321 
Note that if one willing to overwrite the lower triangular part of the evolving matrix A, one can store the evolving ¦« there, since these entries will eventually be zero anyway! There is also no need to save explicitly the permutation matrix P . One could instead record the permutation steps in an extra column (record the vector (¦Ð(1),...,¦Ð(n)) corresponding to the permutation ¦Ð applied to the rows). We let the reader write such a bold and space-e.cient version of LU-decomposition! 
Remark: In Matlab the function lu returns the matrices P, L, U involved in the PA = LU factorization using the call [L, U, P ] = lu(A). 
As a corollary of Theorem 7.5(1), we can show the following result. 
Proposition 7.6. If an invertible real symmetric matrix A has an LU-decomposition, then A has a factorization of the form 
A = LDLT , 
where L is a lower-triangular matrix whose diagonal entries are equal to 1, and where D consists of the pivots. Furthermore, such a decomposition is unique. 
Proof. If A has an LU-factorization, then it has an LDU factorization 
A = LDU, 
where L is lower-triangular, U is upper-triangular, and the diagonal entries of both L and U are equal to 1. Since A is symmetric, we have 
LDU = A = AT = UTDLT , 
with UT lower-triangular and DLT upper-triangular. By the uniqueness of LU-factorization (Part (1) of Theorem 7.5), we must have L = UT (and DU = DLT), thus U = LT , as claimed. 
Remark: It can be shown that Gaussian elimination plus back-substitution requires n3/3+ O(n2) additions, n3/3+ O(n2) multiplications and n2/2+ O(n) divisions. 
7.6 Proof of Theorem 7.5 . 
Proof. (1) The only part that has not been proven is the uniqueness part (when P = I). Assume that A is invertible and that A = L1U1 = L2U2, with L1,L2 unit lower-triangular and U1,U2 upper-triangular. Then we have 
L.1 = U2U.1 
.
2 L11 
However, it is obvious that L.21 is lower-triangular and that U1 .1 is upper-triangular, and so L.1 
2 L1 is lower-triangular and U2U1 .1 is upper-triangular. Since the diagonal entries of L1 and L2 are 1, the above equality is only possible if U2U1 .1 = I, that is, U1 = U2, and so L1 = L2. 
E.1E.1 E.1
(2) When P = I, we have L = 2 ¡¤¡¤¡¤ n.1, where Ek is the product of n . k elementary matrices of the form Ei,k;.¿Ûi , where Ei,k;.¿Ûi subtracts ¦Éi times row k from row i, 
(k)(k)
with ¦Éik = aik /akk ,1 ¡Ü k ¡Ü n . 1, and k +1 ¡Ü i ¡Ü n. Then it is immediately veri.ed that 
.
. 
Ek = 
........ 

1 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 . . ...
.
......
. 0 ¡¤¡¤¡¤ 10 ¡¤¡¤¡¤ 0 
. . ... 
0 ¡¤¡¤¡¤ .¦Ék+1k 1 ¡¤¡¤¡¤ 0 
.... .
.
.. ....
.
.... . 
........ 

, 

0 ¡¤¡¤¡¤ .¦Énk 0 ¡¤¡¤¡¤ 1 

7.6. PROOF OF THEOREM 7.5 ¢à 
and that 

.
. 
1 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 

E.1 
= 
k 
........ 

. . ...
.
......
.
. . ... 
0  ¡¤ ¡¤ ¡¤  1  0  ¡¤ ¡¤ ¡¤  0  
0  ¡¤ ¡¤ ¡¤  ¦Ék+1k  1  ¡¤ ¡¤ ¡¤  0  
. . .  . . .  . . .  . . .  ...  . . .  
0  ¡¤ ¡¤ ¡¤  ¦Énk  0  ¡¤ ¡¤ ¡¤  1  

........ 

. 

.
. 
If we de.ne Lk by 
. 
.
1 0 0 00.0 . 
.
¦É21 1 0 00.0 .
. 
..
¦É31 ¦É32 . 0 0.0 
............ 

............ 

Lk = 
.. .
.
.... 
. . . 10.0 
¦Ék+11 ¦Ék+12 ¡¤¡¤¡¤ ¦Ék+1k 1 ¡¤¡¤¡¤ 0 
.. ..
.
. ... . 
. . . .0.0 
¦Én1 ¦Én2 ¡¤¡¤¡¤ ¦Énk 0 ¡¤¡¤¡¤ 1 
for k =1,...,n . 1, we easily check that L1 = E1 .1, and that 
Lk = Lk.1Ek .1 , 2 ¡Ü k ¡Ü n . 1, 
because multiplication on the right by Ek .1 adds ¦Éi times column i to column k (of the matrix Lk.1) with i>k, and column i of Lk.1 has only the nonzero entry 1 as its ith element. Since 
= E.1 E.1
Lk 1 ¡¤¡¤¡¤ k , 1 ¡Ü k ¡Ü n . 1, 
we conclude that L = Ln.1, proving our claim about the shape of L. 
(3) 
Step 1. Prove (.1). 
First we prove by induction on k that 

Ak+1 = Ekk ¡¤¡¤¡¤ E1 kPk ¡¤¡¤¡¤ P1A, k =1,...,n . 2. 
For k = 1, we have A2 = E1P1A = E11P1A, since E11 = E1, so our assertion holds trivially. Now if k ¡Ý 2, Ak+1 = EkPkAk, 
and by the induction hypothesis, 
= Ek.1 Ek.1Ek.1
Ak ¡¤¡¤¡¤ Pk.1 ¡¤¡¤¡¤ P1A.
k.1 21 
Because Pk is either the identity or a transposition, Pk 2 = I, so by inserting occurrences of PkPk as indicated below we can write 
Ak+1 = EkPkAk Ek.1Ek.1 
= EkPkEk.1 ¡¤¡¤¡¤ Pk.1 ¡¤¡¤¡¤ P1A
k.1 21 = EkPkEk.1(PkPk) ¡¤¡¤¡¤ (PkPk)Ek.1(PkPk)Ek.1(PkPk)Pk.1 ¡¤¡¤¡¤ P1A
k.1 21 = Ek(PkEk.1Pk) ¡¤¡¤¡¤ (PkEk.1Pk)(PkEk.1Pk)PkPk.1 ¡¤¡¤¡¤ P1A.
k.1 21 
Observe that Pk  has been ¡°moved¡± to the right of the elimination steps.  However, by  
de.nition,  
Ek j = PkEk.1 j Pk,  j = 1, . . . , k . 1  
Ek k = Ek,  

so we get 
Ak+1 = EkkEkk .1 ¡¤¡¤¡¤ E2 kE1 kPk ¡¤¡¤¡¤ P1A, establishing the induction hypothesis. For k = n . 2, we get = En.1 En.1
U = An.1 ¡¤¡¤¡¤ Pn.1 ¡¤¡¤¡¤ P1A,
n.11 
as claimed, and the factorization PA = LU with 
P = Pn.1 ¡¤¡¤¡¤ P1 ).1 (En.1).1
L =(En.1 ¡¤¡¤¡¤ 
1 n.1 
is clear. Step 2. Prove that the matrices (Ejk).1 are lower-triangular. To achieve this, we prove that the matrices Ejk are strictly lower triangular matrices of a very special form. Since for j =1,...,n . 2, we have Ejj = Ej, 
Ejk = PkEjk.1Pk,k = j +1,...,n . 1, 
since En.1 = En.1 and P .1 = Pk, we get (Ej).1 = E.1 for j =1,..., n . 1, and for 
n.1 k jj 
j =1,...,n . 2, we have (Ejk).1 = Pk(Ejk.1).1Pk,k = j +1,...,n . 1. Since (Ejk.1).1 = I + Ejk.1 and Pk = P (k, i) is a transposition or Pk = I, so Pk 2 = I, and we get (Ek).1 = Pk(Ek.1).1Pk = Pk(I + Ek.1 = P 2 + Pk Ek.1 
jj j )Pkk j Pk = I + Pk Ejk.1 Pk. 

7.6. PROOF OF THEOREM 7.5 ¢à 
Therefore, we have (Ejk).1 = I + Pk Ejk.1 Pk, 1 ¡Ü j ¡Ü n . 2,j +1 ¡Ü k ¡Ü n . 1. We prove for j =1,...,n . 1, that for k = j, . . . , n . 1, each Ejk is a lower triangular matrix 
of the form 

.
. 
Ek 
= 
j 
......... 

0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 . . ...
.
......
. 0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 
. . ... 
0 ¡¤¡¤¡¤ ¦É

(k) 
j+1j 
0 ¡¤¡¤¡¤ 0 

.... .
.
.. ....
.
.... . 
(k)
0 ¡¤¡¤¡¤ ¦É0 ¡¤¡¤¡¤ 0
nj 
......... 

, 

and that Ejk = Pk Ejk.1 , 1 ¡Ü j ¡Ü n . 2,j +1 ¡Ü k ¡Ü n . 1, 
with Pk = I or Pk = P (k, i) for some i such that k +1 ¡Ü i ¡Ü n. 
For each j (1 ¡Ü j ¡Ü n . 1) we proceed by induction on k = j, . . . , n . 1. Since (Ejj).1 = Ej .1 and since Ej .1 is of the above form, the base case holds. 
For the induction step, we only need to consider the case where Pk = P (k, i) is a trans-position, since the case where Pk = I is trivial. We have to .gure out what Pk Ejk.1 Pk = P (k, i) Ejk.1 P (k, i) is. However, since 
.
. 
Ek.1 
= 
j 
......... 

0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 . . ...
.
......
. 0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 
. . ... 
0 ¡¤¡¤¡¤ ¦É

(k.1) 
j+1j 
0 ¡¤¡¤¡¤ 0 

.... .
.
.. ....
.
.... . 
(k.1)
0 ¡¤¡¤¡¤ ¦É0 ¡¤¡¤¡¤ 0
nj 
......... 

, 

and because k +1 ¡Ü i ¡Ü n and j ¡Ü k . 1, multiplying Ejk.1 on the right by P (k, i) will permute columns i and k, which are columns of zeros, so 
P (k, i) Ejk.1 P (k, i)= P (k, i) Ejk.1 , 
and thus, (Ejk).1 = I + P (k, i) Ejk.1 . 
But since (Ek).1 = I + Ek 
jj , 
we deduce that Ejk = P (k, i) Ejk.1 . 
We also know that multiplying Ejk.1 on the left by P (k, i) will permute rows i and k, which shows that Ejk has the desired form, as claimed. Since all Ejk are strictly lower triangular, all (Ejk).1 = I + Ejk are lower triangular, so the product 
).1 (En.1
L =(E1 n.1 ¡¤¡¤¡¤ n.1 ).1 
is also lower triangular. Step 3. Express L as L = I +¦«n.1, with ¦«n.1 = E11 + ¡¤¡¤¡¤ + Enn..11 . From Step 1 of Part (3), we know that 
).1 (En.1
L =(E1 n.1 ¡¤¡¤¡¤ n.1 ).1 . 
We prove by induction on k that ).1 ).1
I +¦«k =(E1 k ¡¤¡¤¡¤ (Ekk ¦«k = E1 k + ¡¤¡¤¡¤ + Ekk , 
for k =1,...,n . 1. If k = 1, we have E11 = E1 and 
.
. 
E1 = 
....

1  0  ¡¤ ¡¤ ¡¤  0  
.¦É(1) 21  1  ¡¤ ¡¤ ¡¤  0  
. . .  . . .  ...  . . .  
.¦É(1) n1  0  ¡¤ ¡¤ ¡¤  1  

.... 

. 

We also get 

.
. 
(E.1).1 
1 
= 

.... 

1  0  ¡¤ ¡¤ ¡¤  0  
¦É(1) 21  1  ¡¤ ¡¤ ¡¤  0  
. . .  . . .  ...  . . .  
¦É(1) n1  0  ¡¤ ¡¤ ¡¤  1  

.... 

= I +¦«1. 
Since (E.1).1 = I + E1, we .nd that we get ¦«1 = E1, and the base step holds. 
11 1 
Since (Ejk).1 = I + Ejk with 
.
. 
Ek 
= 
j 
......... 

0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 . . ...
.
......
. 0 ¡¤¡¤¡¤ 00 ¡¤¡¤¡¤ 0 
. . ... 
0 ¡¤¡¤¡¤ ¦É

(k) 
j+1j 
0 ¡¤¡¤¡¤ 0 

.... .
.
.. ....
.
.... . 
(k)
0 ¡¤¡¤¡¤ ¦É0 ¡¤¡¤¡¤ 0
nj 
......... 

and EikEjk =0 if i<j, as in part (2) for the computation involving the products of Lk¡¯s, we get 
(Ek.1).1 ¡¤¡¤¡¤ (Ek.1).1 = I + Ek.1 + ¡¤¡¤¡¤ + Ek.1 2 ¡Ü k ¡Ü n. (.)
1 k.11 k.1 , 

7.7. DEALING WITH ROUNDOFF ERRORS; PIVOTING STRATEGIES 221 
Similarly, from the fact that Ejk.1 P (k, i)= Ejk.1 if i ¡Ý k + 1 and j ¡Ü k . 1 and since 
(Ejk).1 = I + PkEjk.1 , 1 ¡Ü j ¡Ü n . 2,j +1 ¡Ü k ¡Ü n . 1, 
we get 
).1 ¡¤ + Ek.1
(Ek).1 ¡¤¡¤¡¤ (Ek = I + Pk(Ek.1 + ¡¤¡¤ ), 2 ¡Ü k ¡Ü n . 1. (..)
1 k.11 k.1 
By the induction hypothesis, ).1
I +¦«k.1 =(Ek.1).1 ¡¤¡¤¡¤ (Ek.1 ,
1 k.1 
and from (.), we get 
= Ek.1 ¡¤ + Ek.1
¦«k.1 + ¡¤¡¤ 
1 k.1 . 
Using (..), we deduce that 
(E1 k).1 ¡¤¡¤¡¤ (Ekk .1).1 = I + Pk¦«k.1. 
Since Ekk = Ek, we obtain ).1
(Ek ¡¤¡¤¡¤ (Ek ).1(Ek).1 =(I + Pk¦«k.1)E.1 .
1 k.1kk 
However, by de.nition I +¦«k =(I + Pk¦«k.1)Ek .1 , 
which proves that ).1 (Ek ).1(Ek).1
I +¦«k =(E1 k ¡¤¡¤¡¤ k.1k , (.) 
and .nishes the induction step for the proof of this formula. If we apply Equation (.) again with k + 1 in place of k, we have 
(Ek).1 ¡¤¡¤¡¤ (Ek).1 = I + Ek + ¡¤¡¤¡¤ + Ekk ,
1 k 1 
and together with (.), we obtain, 
¦«k = E1 k + ¡¤¡¤¡¤ + Ekk , 
also .nishing the induction step for the proof of this formula. For k = n .1 in (.), we obtain the desired equation: L = I +¦«n.1. 
7.7 Dealing with Roundo. Errors; Pivoting Strategies 
Let us now brie.y comment on the choice of a pivot. Although theoretically, any pivot can be chosen, the possibility of roundo. errors implies that it is not a good idea to pick very small pivots. The following example illustrates this point. Consider the linear system 
10.4x + y =1 x + y =2. 
Since 10.4 is nonzero, it can be taken as pivot, and we get 
10.4x + y =1 (1 . 104)y =2 . 104 . 
Thus, the exact solution is 
104 104 . 2 
x = ,y = . 
104 . 1104 . 1 However, if roundo. takes place on the fourth digit, then 104 . 1 = 9999 and 104 . 2 = 9998 will be rounded o. both to 9990, and then the solution is x = 0 and y = 1, very far from the exact solution where x ¡Ö 1 and y ¡Ö 1. The problem is that we picked a very small pivot. If instead we permute the equations, the pivot is 1, and after elimination we get the system 
x + y =2 (1 . 10.4)y =1 . 2 ¡Á 10.4 . 
This time, 1 . 10.4 =0.9999 and 1 . 2 ¡Á 10.4 =0.9998 are rounded o. to 0.999 and the solution is x =1,y = 1, much closer to the exact solution. To remedy this problem, one may use the strategy of partial pivoting. This consists of choosing during Step k (1 ¡Ü k ¡Ü n . 1) one of the entries aik (k) such that 
(k)(k)
|a| = max |a|.
ik pk 
k¡Üp¡Ün 
By maximizing the value of the pivot, we avoid dividing by undesirably small pivots. 
Remark: A matrix, A, is called strictly column diagonally dominant i. 
n |ajj| > |aij|, for j =1,...,n i=1,iJ
=j 
(resp. strictly row diagonally dominant i. 
n |aii| > |aij|, for i =1, . . . , n.) j=1,j=i
J
For example, the matrix

.
. 
...... 

7 
1
2 
141 0 ...
...
... 0 141 
7
1 
2 
...... 

of the curve interpolation problem discussed in Section 7.1 is strictly column (and row) diagonally dominant. 

7.8. GAUSSIAN ELIMINATION OF TRIDIAGONAL MATRICES 
It has been known for a long time (before 1900, say by Hadamard) that if a matrix A is strictly column diagonally dominant (resp. strictly row diagonally dominant), then it is invertible. It can also be shown that if A is strictly column diagonally dominant, then Gaussian elimination with partial pivoting does not actually require pivoting (see Problem 7.12). 
(k)
Another strategy, called complete pivoting, consists in choosing some entry a, where 
ij 
k ¡Ü i, j ¡Ü n, such that 
(k)(k)
|a| = max |a|.
ij pq 
k¡Üp,q¡Ün 
However, in this method, if the chosen pivot is not in column k, it is also necessary to permute columns. This is achieved by multiplying on the right by a permutation matrix. However, complete pivoting tends to be too expensive in practice, and partial pivoting is the method of choice. 
A special case where the LU-factorization is particularly e.cient is the case of tridiagonal matrices, which we now consider. 


7.8 Gaussian Elimination of Tridiagonal Matrices 
Consider the tridiagonal matrix 
.
. 
A = 

.......... 

b1 c1 a2 b2 c2 a3 b3 c3 
.. .
.

.......... 

.

. 

.

. 

.

. 

an.2  bn.2  cn.2  
an.1  bn.1  cn.1  
an  bn  

De.ne the sequence ¦Ä0 =1,¦Ä1 = b1,¦Äk = bk¦Äk.1 . akck.1¦Äk.2, 2 ¡Ü k ¡Ü n. 
Proposition 7.7. If A is the tridiagonal matrix above, then ¦Äk = det(A(1 : k, 1: k)) for k =1,...,n. 
Proof. By expanding det(A(1 : k, 1: k)) with respect to its last row, the proposition follows by induction on k. 
Theorem 7.8. If A is the tridiagonal matrix above and ¦Äk =0 for k =1,...,n, then A has the following LU-factorization: 
.

.

¦Ä1
.
. 
1 

c1
................ 

................ 

¦Ä0¦Ä0
.............. 

a2 
.............. 

¦Ä2
1 

¦Ä1 c2¦Ä1¦Ä1 ¦Ä3
1

a3 ¦Ä2 c3¦Ä2
A = 

.

..
..
.. 
..
..
.. 
¦Än.3 ¦Än.1
1

an.1 ¦Än.2 cn.1¦Än.2¦Än.2 
¦Än 
¦Än.1 
¦Än.1 
1

an 
Proof. Since ¦Äk = det(A(1 : k, 1: k)) = 0 for k =1,...,n, by Theorem 7.5 (and Proposition 7.2), we know that A has a unique LU-factorization. Therefore, it su.ces to check that the proposed factorization works. We easily check that 
(LU)kk+1 = ck, 1 ¡Ü k ¡Ü n . 1 
(LU)kk.1 = ak, 2 ¡Ü k ¡Ü n (LU)kl =0, |k . l|¡Ý 2 ¦Ä1
(LU)11 == b1
¦Ä0 akck.1¦Äk.2 + ¦Äk
(LU)kk == bk, 2 ¡Ü k ¡Ü n,¦Äk.1 
since ¦Äk = bk¦Äk.1 . akck.1¦Äk.2. 
It follows that there is a simple method to solve a linear system Ax = d where A is tridiagonal (and ¦Äk = 0 for k =1,...,n). For this, it is convenient to ¡°squeeze¡± the diagonal matrix ¦¤ de.ned such that ¦¤kk = ¦Äk/¦Äk.1 into the factorization so that A =(L¦¤)(¦¤.1U), and if we let 
c1 ¦Äk.1 ¦Än 
z1 = ,zk = ck , 2 ¡Ü k ¡Ü n . 1,zn == bn . anzn.1,b1 ¦Äk ¦Än.1 
7.9. SPD MATRICES AND THE CHOLESKY DECOMPOSITION 
A =(L¦¤)(¦¤.1U) is written as 
.
. 
...................... 

1 z1 1 
...................... 

. 

.
. 
c1 z2 
1 z3 
............. 

............. 

z1  
c2  
a2  
z2  c3  
a3  

z3 
..

A = 

.

.

. 

. 

...  ...  
cn.1  
an.1  
zn.1  
an  zn  

1 zn.2 
1 zn.1 1 
As a consequence, the system Ax = d can be solved by constructing three sequences: First, the sequence 
c1 ck 
z1 = ,zk = ,k =2,...,n . 1,zn = bn . anzn.1,
b1 bk . akzk.1 
corresponding to the recurrence ¦Äk = bk¦Äk.1 . akck.1¦Äk.2 and obtained by dividing both sides of this equation by ¦Äk.1, next 
d1 dk . akwk.1 
w1 = ,wk = ,k =2, . . . , n, 
b1 bk . akzk.1 
corresponding to solving the system L¦¤w = d, and .nally 
xn = wn,xk = wk . zkxk+1,k = n . 1,n . 2,..., 1, 
corresponding to solving the system ¦¤.1Ux = w. 
Remark: It can be veri.ed that this requires 3(n . 1) additions, 3(n . 1) multiplications, and 2n divisions, a total of 8n . 6 operations, which is much less that the O(2n3/3) required by Gaussian elimination in general. 
We now consider the special case of symmetric positive de.nite matrices (SPD matrices). 


7.9 SPD Matrices and the Cholesky Decomposition 
Recall that an n ¡Á n real symmetric matrix A is positive de.nite i. 
x TAx > 0 for all x ¡Ê Rn with x =0. 
Equivalently, A is symmetric positive de.nite i. all its eigenvalues are strictly positive. The following facts about a symmetric positive de.nite matrix A are easily established (some left as an exercise): 
(1) 
The matrix A is invertible. (Indeed, if Ax = 0, then xTAx = 0, which implies x = 0.) 

(2) 
We have aii > 0 for i =1,...,n. (Just observe that for x = ei, the ith canonical basis 


T
vector of Rn, we have ei Aei = aii > 0.) 
(3) 
For every n ¡Á n real invertible matrix Z, the matrix ZTAZ is real symmetric positive de.nite i. A is real symmetric positive de.nite. 

(4) 
The set of n ¡Á n real symmetric positive de.nite matrices is convex. This means that if A and B are two n ¡Á n symmetric positive de.nite matrices, then for any ¦Ë ¡Ê R such that 0 ¡Ü ¦Ë ¡Ü 1, the matrix (1 . ¦Ë)A + ¦ËB is also symmetric positive de.nite. Clearly since A and B are symmetric, (1 . ¦Ë)A + ¦ËB is also symmetric. For any nonzero x ¡Ê Rn, we have xTAx > 0 and xTBx > 0, so 


x T((1 . ¦Ë)A + ¦ËB)x = (1 . ¦Ë)x TAx + ¦ËxTBx > 0, 
because 0 ¡Ü ¦Ë ¡Ü 1, so1.¦Ë ¡Ý 0 and ¦Ë ¡Ý 0, and1.¦Ë and ¦Ë can¡¯t be zero simultaneously. 
(5) The set of n ¡Á n real symmetric positive de.nite matrices is a cone. This means that if A is symmetric positive de.nite and if ¦Ë> 0 is any real, then ¦ËA is symmetric positive de.nite. Clearly ¦ËA is symmetric, and for nonzero x ¡Ê Rn, we have xTAx > 0, and since ¦Ë> 0, we have xT¦ËAx = ¦ËxTAx > 0. 
Remark: Given a complex m ¡Á n matrix A, we de.ne the matrix A as the m ¡Á n matrix A =(aij). Then we de.ne A. as the n ¡Á m matrix A. =(A)T =(AT). The n ¡Á n complex matrix A is Hermitian if A. = A. This is the complex analog of the notion of a real symmetric matrix. A Hermitian matrix A is positive de.nite if 
z . Az > 0 for all z ¡Ê Cn with z =0. 
It is easily veri.ed that Properties (1)-(5) hold for Hermitian positive de.nite matrices; replace T by .. 
It is instructive to characterize when a 2 ¡Á 2 real symmetric matrix A is positive de.nite. 
Write  
A =  a c  c b  .  
Then we have 
 x  y  a c  c b  x y  = ax 2 + 2cxy + by2 .     

If the above expression is strictly positive for all nonzero vectorsxy, then for x =1,y =0 we get a> 0 and for x =0,y = 1 we get b> 0. Then we can write 
2
¡Ì c 2 c2 
ax 2 +2cxy + by2 = ax + ¡Ì y + by2 . y 
aa 
2   
¡Ì c 1 22 
= ax + ¡Ì y + ab . c y. (.) 
aa
7.9. SPD MATRICES AND THE CHOLESKY DECOMPOSITION 
Since a> 0, if ab . c2 ¡Ü 0, then we can choose y> 0 so that the second term is negative or zero, and we can set x = .(c/a)y to make the .rst term zero, in which case ax2 +2cxy+by2 ¡Ü 0, so we must have ab . c2 > 0. 
Conversely, if a> 0,b > 0 and ab > c2, then for any (x, y) = (0, 0), if y = 0, then x =0 and the .rst term of (.) is positive, and if y = 0, then the second term of (.) is positive. Therefore, the symmetric matrix A is positive de.nite i. 
a> 0, b> 0, ab>c2 . (.) 
Note that ab . c2 = det(A), so the third condition says that det(A) > 0. 
Observe that the condition b> 0 is redundant, since if a> 0 and ab > c2, then we must have b> 0 (and similarly b> 0 and ab > c2 implies that a> 0). 
We can try to visualize the space of 2 ¡Á 2 real symmetric positive de.nite matrices in R3, by viewing (a, b, c) as the coordinates along the x, y, z axes. Then the locus determined by the strict inequalities in (.) corresponds to the region on the side of the cone of equation xy = z2 that does not contain the origin and for which x> 0 and y> 0. For z = ¦Ä .xed, the equation xy = ¦Ä2 de.ne a hyperbola in the plane z = ¦Ä. The cone of equation xy = z2 consists of the lines through the origin that touch the hyperbola xy = 1 in the plane z = 1. We only consider the branch of this hyperbola for which x> 0 and y> 0. See Figure 7.6. 
It is not hard to show that the inverse of a real symmetric positive de.nite matrix is also real symmetric positive de.nite, but the product of two real symmetric positive de.nite matrices may not be symmetric positive de.nite, as the following example shows: 
¡Ì¡Ì ¡Ì 
11 1/ 2 .12 02/ 2
¡Ì¡Ì = ¡Ì¡Ì . 
12 .1/ 23/ 2 .1/ 25/ 2 
According to the above criterion, the two matrices on the left-hand side are real symmetric positive de.nite, but the matrix on the right-hand side is not even symmetric, and 
¡Ì¡Ì ¡Ì
02/ 2 .62/ 2 
.61 ¡Ì¡Ì = .61 ¡Ì = .1/ 5,
.1/ 25/ 21 11/ 2 
even though its eigenvalues are both real and positive. 
Next we show that a real symmetric positive de.nite matrix has a special LU-factorization of the form A = BBT, where B is a lower-triangular matrix whose diagonal elements are strictly positive. This is the Cholesky factorization. 
First we note that a symmetric positive de.nite matrix satis.es the condition of Propo-sition 7.2. 
Proposition 7.9. If A is a real symmetric positive de.nite matrix, then A(1 : k, 1: k) is symmetric positive de.nite and thus invertible for k =1,...,n. 
Proof. Since A is symmetric, each A(1 : k, 1: k) is also symmetric. If w ¡Ê Rk , with 1 ¡Ü k ¡Ü n, we let x ¡Ê Rn be the vector with xi = wi for i =1,...,k and xi = 0 for Figure 7.6: Two views of the surface xy = z2 in R3 . The intersection of the surface with a constant z plane results in a hyperbola. The region associated with the 2 ¡Á 2 symmetric positive de.nite matrices lies in ¡±front¡± of the green side. 

i = k +1,...,n. Now since A is symmetric positive de.nite, we have xTAx > 0 for all x ¡Ê Rn with x = 0. This holds in particular for all vectors x obtained from nonzero vectors w ¡Ê Rk as de.ned earlier, and clearly 
x TAx = w TA(1 : k, 1: k) w, 
which implies that A(1 : k, 1: k) is positive de.nite. Thus, by Fact 1 above, A(1 : k, 1: k) is also invertible. 
Proposition 7.9 also holds for a complex Hermitian positive de.nite matrix. Proposition 
7.9 can be strengthened as follows: A real symmetric (or complex Hermitian) matrix A is positive de.nite i. det(A(1 : k, 1: k)) > 0 for k =1,...,n. The above fact is known as Sylvester¡¯s criterion. We will prove it after establishing the 
Cholesky factorization. Let A be an n ¡Á n real symmetric positive de.nite matrix and write 
W T
a11 
A = ,
WC 

7.9. SPD MATRICES AND THE CHOLESKY DECOMPOSITION 
where C is an(n . 1) ¡Á (n . 1) symmetric matrix and W is an (n . 1) ¡Á 1 matrix. Since A
¡Ì 
is symmetric positive de.nite, a11 > 0, and we can compute ¦Á = a11. The trick is that we can factor A uniquely as 
a11 W T ¦Á 01 0 ¦ÁW T/¦Á
A == ,
W C W/¦Á I 0 C . WW T/a11 0 I 
i.e., as A = B1A1B1 T, where B1 is lower-triangular with positive diagonal entries. Thus, B1 is invertible, and by Fact (3) above, A1 is also symmetric positive de.nite. 
Remark: The matrix C .WW T/a11 is known as the Schur complement of the matrix (a11). 
Theorem 7.10. (Cholesky factorization) Let A be a real symmetric positive de.nite matrix. Then there is some real lower-triangular matrix B so that A = BBT . Furthermore, B can be chosen so that its diagonal elements are strictly positive, in which case B is unique. 
Proof. We proceed by induction on the dimension n of A. For n = 1, we must have a11 > 0,
¡Ì 
and if we let ¦Á = a11 and B =(¦Á), the theorem holds trivially. If n ¡Ý 2, as we explained above, again we must have a11 > 0, and we can write 
a11 W T ¦Á 01 0 ¦ÁW T/¦Á
A == = B1A1BT 
W C W/¦Á I 0 C . WW T/a11 0 I 1 , 
¡Ì 
where ¦Á = a11, the matrix B1 is invertible and 
10 
A1 = 0 C . WW T/a11 
is symmetric positive de.nite. However, this implies that C . WW T/a11 is also symmetric positive de.nite (consider xTA1x for every x ¡Ê Rn with x = 0 and x1 = 0). Thus, we can apply the induction hypothesis to C . WW T/a11 (which is an (n . 1) ¡Á (n . 1) matrix), and we .nd a unique lower-triangular matrix L with positive diagonal entries so that 
C . WW T/a11 = LLT . 
But then we get 
¦Á 01 0 ¦ÁW T/¦ÁA = W/¦Á I 0 C . WW T/a11 0 I ¦Á 0 10 ¦ÁW T/¦Á
= 
W/¦Á I 0 LLT 0 I ¦Á 0 10 10 ¦ÁW T/¦Á
= 
W/¦Á I 0 L 0 LT 0 I ¦Á 0 ¦ÁW T/¦Á
= . 
W/¦Á L 0 LT 
Therefore, if we let  
B =  ¦Á W/¦Á  0 L  ,  

we have a unique lower-triangular matrix with positive diagonal entries and A = BBT . 
Remark: The uniqueness of the Cholesky decomposition can also be established using the uniqueness of an LU-decomposition. Indeed, if A = B1B1 T = B2B2 T where B1 and B2 are lower triangular with positive diagonal entries, if we let ¦¤1 (resp. ¦¤2) be the diagonal matrix consisting of the diagonal entries of B1 (resp. B2) so that (¦¤k)ii =(Bk)ii for k =1, 2, then we have two LU-decompositions 
A =(B1¦¤.11)(¦¤1B1 T)=(B2¦¤.21)(¦¤2B2 T) 
with B1¦¤.112 1 , ¦¤2B2 T upper triangular. By uniquenes 
,B2¦¤.1 unit lower triangular, and ¦¤1BT of LU-factorization (Theorem 7.5(1)), we have 
B1¦¤.1 = B2¦¤.1 , ¦¤1BT =¦¤2B2 T ,
12 1 
and the second equation yields B1¦¤1 = B2¦¤2. (.) 
The diagonal entries of B1¦¤1 are (B1)2 and similarly the diagonal entries of B2¦¤2 are (B2)2 
ii ii, 
so the above equation implies that 
(B1)2 =(B2)2 i =1, . . . , n. 
ii ii, 
Since the diagonal entries of both B1 and B2 are assumed to be positive, we must have 
(B1)ii =(B2)ii,i =1,...,n; 
that is, ¦¤1 =¦¤2, and since both are invertible, we conclude from (.) that B1 = B2. 
Theorem 7.10 also holds for complex Hermitian positive de.nite matrices. In this case, we have A = BB. for some unique lower triangular matrix B with positive diagonal entries. 
The proof of Theorem 7.10 immediately yields an algorithm to compute B from A by solving for a lower triangular matrix B such that A = BBT (where both A and B are real matrices). For j =1,...,n, 
 j.1  1/2 b2
bjj =ajj . jk, k=1 
and for i = j +1,...,n (and j =1,...,n . 1) 
  
j.1 
bij =aij . bikbjk/bjj. k=1 

7.9. SPD MATRICES AND THE CHOLESKY DECOMPOSITION 
The above formulae are used to compute the jth column of B from top-down, using the .rst j . 1 columns of B previously computed, and the matrix A. In the case of n = 3, A = BBT yields 
.
..
..
. 
a11 a12 a31 b11 00 b11 b21 b31 
.

a21 a22 a32
.

=

.

b21 b22 0
.
.

0 b22 b32
. 

a31 a32 a33 b31 b32 b33 00 b33 
.
. 
b2 
11 b11b21 b11b31 
=

.

b2
b11b21 21 + b2 b21b31 + b22b32
22 
.

. 

b2
b11b31 b21b31 + b22b32 31 + b2 33
32 + b2 
We work down the .rst column of A, compare entries, and discover that 
¡Ì 
= b2 
a11 11 b11 = a11 
a21 a21 = b11b21 b21 = b11 a31 a31 = b11b31 b31 = . b11 
Next we work down the second column of A using previously calculated expressions for b21 and b31 to .nd that 
= b2 
a22 22 b22 = a22 . b2 
21 + b2 21 
1 
2 
a32 . b21b31 a32 = b21b31 + b22b32 b32 = . b22 
Finally, we use the third column of A and the previously calculated expressions for b31 and b32 to determine b33 as 
1 
2
a33 = b312 + b2 33 b33 = a33 . b2 . b322 
32 + b2 31 
. 

For another example, if 

.
. 
....... 

1  1  1  1  1  1  
1  2  2  2  2  2  
1  2  3  3  3  3  

1  2  3  4  4  4  
1  2  3  4  5  5  
1  2  3  4  5  6  

....... 

A = 

, 

we .nd that 

.
. 
....... 

1  0  0  0  0  0  
1  1  0  0  0  0  
1  1  1  0  0  0  

1  1  1  1  0  0  
1  1  1  1  1  0  
1  1  1  1  1  1  

....... 

B = 

. 

We leave it as an exercise to .nd similar formulae (involving conjugation) to factor a complex Hermitian positive de.nite matrix A as A = BB. . The following Matlab program implements the Cholesky factorization. 
function B = Cholesky(A) n = size(A,1); B = zeros(n,n); for j = 1:n-1; 
if 	j == 1 
B(1,1) = sqrt(A(1,1)); 
fori =2:n 

B(i,1) = A(i,1)/B(1,1); 
end 

else 
B(j,j) = sqrt(A(j,j) -B(j,1:j-1)*B(j,1:j-1)¡¯); 
for i = j+1:n 

B(i,j) = (A(i,j) -B(i,1:j-1)*B(j,1:j-1)¡¯)/B(j,j); 
end 

end end B(n,n) = sqrt(A(n,n) -B(n,1:n-1)*B(n,1:n-1)¡¯); end 
If we run the above algorithm on the following matrix 
.
. 
A = 

..... 

41000 
14100 

01410 

00141 
00014 

..... 

, 

we obtain 

.
. 
B = 

..... 

2.00000 0 0 0 0.5000 1.93650 0 0 00.5164 1.93220 0 
0 00.5175 1.9319 0 0 0 00.5176 1.9319 
..... 

. 

The Cholesky factorization can be used to solve linear systems Ax = b where A is symmetric positive de.nite: Solve the two systems Bw = b and BTx = w. 
Remark: It can be shown that this methods requires n3/6+ O(n2) additions, n3/6+ O(n2) multiplications, n2/2+O(n) divisions, and O(n) square root extractions. Thus, the Cholesky method requires half of the number of operations required by Gaussian elimination (since Gaussian elimination requires n3/3+ O(n2) additions, n3/3+ O(n2) multiplications, and 

7.9. SPD MATRICES AND THE CHOLESKY DECOMPOSITION 
n2/2+ O(n) divisions). It also requires half of the space (only B is needed, as opposed to both L and U). Furthermore, it can be shown that Cholesky¡¯s method is numerically stable (see Trefethen and Bau [68], Lecture 23). In Matlab the function chol returns the lower-triangular matrix B such that A = BBT using the call B = chol(A, ¡®lower¡¯). 
Remark: If A = BBT, where B is any invertible matrix, then A is symmetric positive de.nite. 
Proof. Obviously, BBT is symmetric, and since B is invertible, BT is invertible, and from 
x TAx = x TBBT x =(BT x)TBT x, 
it is clear that xTAx > 0 if x = 0. 
We now give three more criteria for a symmetric matrix to be positive de.nite. 
Proposition 7.11. Let A be any n ¡Á n real symmetric matrix. The following conditions are equivalent: 
(a) 
A is positive de.nite. 

(b) 
All principal minors of A are positive; that is: det(A(1 : k, 1: k)) > 0 for k =1,...,n (Sylvester¡¯s criterion). 

(c) 
A has an LU-factorization and all pivots are positive. 

(d) 
A has an LDLT-factorization and all pivots in D are positive. 


Proof. By Proposition 7.9, if A is symmetric positive de.nite, then each matrix A(1 : k, 1: k) is symmetric positive de.nite for k =1,...,n. By the Cholsesky decomposition, A(1 : k, 1: k)= QTQ for some invertible matrix Q, so det(A(1 : k, 1: k)) = det(Q)2 > 0. This shows that (a) implies (b). 
If det(A(1 : k, 1: k)) > 0 for k =1,...,n, then each A(1 : k, 1: k) is invertible. By Proposition 7.2, the matrix A has an LU-factorization, and since the pivots ¦Ðk are given by 
. 
.a11 = det(A(1 : 1, 1 : 1)) if k =1 ¦Ðk = det(A(1 : k, 1: k))
. if k =2, . . . , n, 
det(A(1 : k . 1, 1: k . 1)) 
we see that ¦Ðk > 0 for k =1,...,n. Thus (b) implies (c). 
Assume A has an LU-factorization and that the pivots are all positive. Since A is symmetric, this implies that A has a factorization of the form 
A = LDLT , 
with L lower-triangular with 1s on its diagonal, and where D is a diagonal matrix with positive entries on the diagonal (the pivots). This shows that (c) implies (d). 
Given a factorization A = LDLT with all pivots in D positive, if we form the diagonal matrix 
¡Ì ¡Ì¡Ì 
D = diag( ¦Ð1,..., ¦Ðn) 
¡Ì 
and if we let B = LD, then we have 
A = BBT , 
with B lower-triangular and invertible. By the remark before Proposition 7.11, A is positive de.nite. Hence, (d) implies (a). 
Criterion (c) yields a simple computational test to check whether a symmetric matrix is positive de.nite. There is one more criterion for a symmetric matrix to be positive de.nite: its eigenvalues must be positive. We will have to learn about the spectral theorem for symmetric matrices to establish this criterion. 
Proposition 7.11 also holds for complex Hermitian positive de.nite matrices, where in (d), the factorization LDLT is replaced by LDL. . 
For more on the stability analysis and e.cient implementation methods of Gaussian elimination, LU-factoring and Cholesky factoring, see Demmel [16], Trefethen and Bau [68], Ciarlet [14], Golub and Van Loan [30], Meyer [48], Strang [63, 64], and Kincaid and Cheney [39]. 


7.10 Reduced Row Echelon Form (RREF) 
Gaussian elimination described in Section 7.2 can also be applied to rectangular matrices. This yields a method for determining whether a system Ax = b is solvable and a description of all the solutions when the system is solvable, for any rectangular m ¡Á n matrix A. 
It turns out that the discussion is simpler if we rescale all pivots to be 1, and for this we need a third kind of elementary matrix. For any ¦Ë = 0, let Ei,¦Ë be the n ¡Á n diagonal matrix 
.
. 
Ei,¦Ë 
= 

.......... 

1  
...  
1  

¦Ë 

1  
...  
1  

.......... 

, 

with (Ei,¦Ë)ii = ¦Ë (1 ¡Ü i ¡Ü n). Note that Ei,¦Ë is also given by Ei,¦Ë = I +(¦Ë . 1)eii, 
and that Ei,¦Ë is invertible with 
E.1 
i,¦Ë = Ei,¦Ë.1 . 
7.10. REDUCED ROW ECHELON FORM 
Now after k . 1 elimination steps, if the bottom portion 
(k)(k)(k)
(akk ,ak+1k,...,amk) 
of the kth column of the current matrix Ak is nonzero so that a pivot ¦Ðk can be chosen, after a permutation of rows if necessary, we also divide row k by ¦Ðk to obtain the pivot 1, and not only do we zero all the entries i = k +1,...,m in column k, but also all the entries i =1,...,k . 1, so that the only nonzero entry in column k is a 1 in row k. These row operations are achieved by multiplication on the left by elementary matrices. 
(k)(k)(k)
If a= a= ¡¤¡¤¡¤ = a= 0, we move on to column k + 1. 
kk k+1k mk 
When the kth column contains a pivot, the kth stage of the procedure for converting a matrix to rref consists of the following three steps illustrated below: 
.
..
. 
1 ¡Á 0 ¡Á ¡Á¡Á¡Á 1 ¡Á 0 ¡Á ¡Á¡Á¡Á 

....... 

pivot 
=. 
....... 
000 ¡Á ¡Á¡Á¡Á 000 ¡Á ¡Á¡Á¡Á 
....... 

....... 

001 

001

¡Á ¡Á¡Á¡Á 

¡Á ¡Á¡Á¡Á 
(k)
000 

¡Á ¡Á¡Á¡Á 

rescale
.
000 

¡Á¡Á¡Á

a

ik 
=

000 ¡Á ¡Á¡Á¡Á 

000 ¡Á ¡Á¡Á¡Á 

(k)
000 a
000

¡Á¡Á¡Á 

¡Á ¡Á¡Á¡Á

ik 
.
..
. 
1 ¡Á 0 ¡Á¡Á¡Á¡Á 1 ¡Á 0 0 ¡Á¡Á¡Á 

....... 
elim 
=. 
000 ¡Á¡Á¡Á¡Á 000 0 ¡Á¡Á¡Á 
If the kth column does not contain a pivot, we simply move on to the next column. 
The result is that after performing such elimination steps, we obtain a matrix that has a special shape known as a reduced row echelon matrix, for short rref. 
Here is an example illustrating this process: Starting from the matrix 
....... 

....... 

....... 

001 

001 0

¡Á¡Á¡Á¡Á 

¡Á¡Á¡Á 

0001 

000 1

¡Á¡Á¡Á 

¡Á¡Á¡Á 

. 

000 

000 0 ¡Á¡Á¡Á

¡Á¡Á¡Á¡Á 

000 

000 0 ¡Á¡Á¡Á

¡Á¡Á¡Á¡Á 

.
. 
. 

1021 5 
1152 7

.

A1 = 
, 

1 2 8 4 12 
we perform the following steps 
.
. 
10215 

A1 .¡ú A2 
=

.

01312

.

, 

02637 
by subtracting row 1 from row 2 and row 3; 
.
.
.
..
. 
10215 1021 5 102 1 
5 
A2 .¡ú 
.

02637

. .¡ú 

.

0133/27/2

. .¡ú A3 
=

.

013 3/27/2

.

, 

01312 013 1 2 000 .1/2 .3/2 

after choosing the pivot 2 and permuting row 2 and row 3, dividing row 2 by 2, and sub-
tracting row 2 from row 3;  
.  .  .  .  
1  0  2  1  5  1  0  2  0  2  
A3 .¡ú  .0  1  3  3/2  7/2. .¡ú A4 = .0  1  3  0  .1.,  
0  0  0  1  3  0  0  0  1  3  

after dividing row 3 by .1/2, subtracting row 3 from row 1, and subtracting (3/2) ¡Á row 3 from row 2. 
It is clear that columns 1, 2 and 4 are linearly independent, that column 3 is a linear combination of columns 1 and 2, and that column 5 is a linear combination of columns 1, 2, 4. 
In general, the sequence of steps leading to a reduced echelon matrix is not unique. For example, we could have chosen 1 instead of 2 as the second pivot in matrix A2. Nevertheless, the reduced row echelon matrix obtained from any given matrix is unique; that is, it does not depend on the the sequence of steps that are followed during the reduction process. This fact is not so easy to prove rigorously, but we will do it later. 
If we want to solve a linear system of equations of the form Ax = b, we apply elementary row operations to both the matrix A and the right-hand side b. To do this conveniently, we form the augmented matrix (A, b), which is the m ¡Á (n + 1) matrix obtained by adding b as an extra column to the matrix A. For example if 
.. .. 
1021 	5 
.. ..
A = 	1152 and b =7 , 1284 12 
then the augmented matrix is 
.	. 
10	21 5 
.	.
(A, b)= 	1152 7 . 1 2 8 4 12 
Now for any matrix M, since M(A, b)=(MA, Mb), 
performing elementary row operations on (A, b) is equivalent to simultaneously performing operations on both A and b. For example, consider the system 
x1 +2x3 + x4 =5 
x1 + x2 +5x3 +2x4 =7 
x1 +2x2 +8x3 +4x4 = 12. 
Its augmented matrix is the matrix 
.	. 
1021 5 (A, b)= .1152 7 . 
1 2 8 4 12 


7.10. REDUCED ROW ECHELON FORM 
considered above, so the reduction steps applied to this matrix yield the system 
x1 +2x3 =2 x2 +3x3 = .1 x4 =3. 
This reduced system has the same set of solutions as the original, and obviously x3 can be chosen arbitrarily. Therefore, our system has in.nitely many solutions given by 
x1 =2 . 2x3,x2 = .1 . 3x3,x4 =3, 
where x3 is arbitrary. 
The following proposition shows that the set of solutions of a system Ax = b is preserved by any sequence of row operations. 
Proposition 7.12. Given any m ¡Á n matrix A and any vector b ¡Ê Rm, for any sequence of elementary row operations E1,...,Ek, if P = Ek ¡¤¡¤¡¤ E1 and (AÊ¬,bÊ¬)= P (A, b), then the 
solutions of Ax = b are the same as the solutions of AÊ¬x = bÊ¬ . 
Proof. Since each elementary row operation Ei is invertible, so is P , and since (AÊ¬,bÊ¬)= P (A, b), then AÊ¬ = PA and bÊ¬ = Pb. If x is a solution of the original system Ax = b, then multiplying both sides by P we get P Ax = Pb; that is, AÊ¬x = bÊ¬, so x is a solution of the new system. Conversely, assume that x is a solution of the new system, that is AÊ¬x = bÊ¬ . Then because AÊ¬ = PA, bÊ¬ = Pb, and P is invertible, we get 
Ax = P .1AÊ¬ x = P .1bÊ¬ = b, 
so x is a solution of the original system Ax = b. 
Another important fact is this: 
Proposition 7.13. Given an m¡Án matrix A, for any sequence of row operations E1,...,Ek, if P = Ek ¡¤¡¤¡¤ E1 and B = PA, then the subspaces spanned by the rows of A and the rows of B are identical. Therefore, A and B have the same row rank. Furthermore, the matrices A and B also have the same (column) rank. 
Proof. Since B = PA, from a previous observation, the rows of B are linear combinations of the rows of A, so the span of the rows of B is a subspace of the span of the rows of A. Since P is invertible, A = P .1B, so by the same reasoning the span of the rows of A is a subspace of the span of the rows of B. Therefore, the subspaces spanned by the rows of A and the rows of B are identical, which implies that A and B have the same row rank. 
Proposition 7.12 implies that the systems Ax = 0 and Bx = 0 have the same solutions. Since Ax is a linear combinations of the columns of A and Bx is a linear combinations of the columns of B, the maximum number of linearly independent columns in A is equal to the maximum number of linearly independent columns in B; that is, A and B have the same rank. 
Remark: The subspaces spanned by the columns of A and B can be di.erent! However, their dimension must be the same. 
We will show in Section 7.14 that the row rank is equal to the column rank. This will also be proven in Proposition 10.13 Let us now de.ne precisely what is a reduced row echelon matrix. 
De.nition 7.4. An m ¡Á n matrix A is a reduced row echelon matrix i. the following con-ditions hold: 
(a) 
The .rst nonzero entry in every row is 1. This entry is called a pivot. 

(b) 
The .rst nonzero entry of row i + 1 is to the right of the .rst nonzero entry of row i. 


(c) The entries above a pivot are zero. If a matrix satis.es the above conditions, we also say that it is in reduced row echelon form, 
for short rref . Note that Condition (b) implies that the entries below a pivot are also zero. For example, 
the matrix 

.
. 
. 

1601 
0012

.

A = 

0000 
is a reduced row echelon matrix. In general, a matrix in rref has the following shape: 

.
. 
......... 

100 ¡Á¡Á 00 ¡Á 010 ¡Á¡Á 00 ¡Á 001 ¡Á¡Á 00 ¡Á 0000 010 ¡Á 0000 001 ¡Á 0000 0000 0000 0000 
......... 

if the last row consists of zeros, or 

.
. 
....... 

100 ¡Á¡Á 00 ¡Á 0 ¡Á 010 ¡Á¡Á 00 ¡Á 0 ¡Á 001 ¡Á¡Á 00 ¡Á 0 ¡Á 0000 010 ¡Á 0 ¡Á 0000 001 ¡Á¡Á 0 0000 0000 1 ¡Á 
....... 

if the last row contains a pivot. 
The following proposition shows that every matrix can be converted to a reduced row echelon form using row operations. 

7.10. REDUCED ROW ECHELON FORM 
Proposition 7.14. Given any m ¡Á n matrix A, there is a sequence of row operations E1,...,Ek such that if P = Ek ¡¤¡¤¡¤ E1, then U = PA is a reduced row echelon matrix. 
Proof. We proceed by induction on m. If m = 1, then either all entries on this row are zero, so A = 0, or if aj is the .rst nonzero entry in A, let P =(a .1) (a 1 ¡Á 1 matrix); clearly, PA 
j 
is a reduced row echelon matrix. 
Let us now assume that m ¡Ý 2. If A = 0, we are done, so let us assume that A = 0. Since A = 0, there is a leftmost column j which is nonzero, so pick any pivot ¦Ð = aij in the jth column, permute row i and row 1 if necessary, multiply the new .rst row by ¦Ð.1, and clear out the other entries in column j by subtracting suitable multiples of row 1. At the end of this process, we have a matrix A1 that has the following shape: 
....
... 
0 ¡¤¡¤¡¤ 01 . ¡¤¡¤¡¤ . 0 ¡¤¡¤¡¤ 00 . ¡¤¡¤¡¤ . 
... 
A1 = ,
... 
0 ¡¤¡¤¡¤ 00 . ¡¤¡¤¡¤ . 
... 
where . stands for an arbitrary scalar, or more concisely 
... 
01 B 
.... 
A1 = ,
00 D 
where D is a (m . 1) ¡Á (n . j) matrix (and B is a 1 ¡Á n . j matrix). If j = n, we are done. 
Otherwise, by the induction hypothesis applied to D, there is a sequence of row operations that converts D to a reduced row echelon matrix RÊ¬, and these row operations do not a.ect 
the .rst row of A1, which means that A1 is reduced to a matrix of the form 
01 B R = 00 RÊ¬ . 
Because RÊ¬ is a reduced row echelon matrix, the matrix R satis.es Conditions (a) and (b) of the reduced row echelon form. Finally, the entries above all pivots in RÊ¬ can be cleared out by subtracting suitable multiples of the rows of RÊ¬ containing a pivot. The resulting matrix 
also satis.es Condition (c), and the induction step is complete. 
Remark: There is a Matlab function named rref that converts any matrix to its reduced row echelon form. 
If A is any matrix and if R is a reduced row echelon form of A, the second part of Proposition 7.13 can be sharpened a little, since the structure of a reduced row echelon matrix makes it clear that its rank is equal to the number of pivots. 
Proposition 7.15. The rank of a matrix A is equal to the number of pivots in its rref R. 
.
. 


7.11 	RREF, Free Variables, and Homogenous Linear Systems 
Given a system of the form Ax = b, we can apply the reduction procedure to the augmented matrix (A, b) to obtain a reduced row echelon matrix (AÊ¬,bÊ¬) such that the system AÊ¬x = bÊ¬ 
has the same solutions as the original system Ax = b. The advantage of the reduced system AÊ¬x = bÊ¬ is that there is a simple test to check whether this system is solvable, and to .nd 
its solutions if it is solvable. Indeed, if any row of the matrix AÊ¬ is zero and if the corresponding entry in bÊ¬ is nonzero, 
then it is a pivot and we have the ¡°equation¡± 
0=1, 
which means that the system AÊ¬x = bÊ¬ has no solution. On the other hand, if there is no pivot in bÊ¬, then for every row i in which bÊ¬ = 0, there is some column j in AÊ¬ where the entry on row i is 1 (a pivot). Consequently,iwe can assign arbitrary values to the variable xk if column k does not contain a pivot, and then solve for the pivot variables. 
For example, if we consider the reduced row echelon matrix 
.. 
16010 
(AÊ¬,bÊ¬)= .00120. , 00001 there is no solution to AÊ¬x = bÊ¬ because the third equation is 0 = 1. On the other hand, the 
reduced system  .  .  
1  6  0  1  1  
(AÊ¬, bÊ¬) = .0  0  1  2  3.  
0  0  0  0  0  

has solutions. We can pick the variables x2,x4 corresponding to nonpivot columns arbitrarily, and then solve for x3 (using the second equation) and x1 (using the .rst equation). 
The above reasoning proves the following theorem: 
Theorem 7.16. Given any system Ax = b where A is a m ¡Á n matrix, if the augmented matrix (A, b) is a reduced row echelon matrix, then the system Ax = b has a solution i. there is no pivot in b. In that case, an arbitrary value can be assigned to the variable xj if column j does not contain a pivot. 
De.nition 7.5. Nonpivot variables are often called free variables. 
Putting Proposition 7.14 and Theorem 7.16 together we obtain a criterion to decide 
whether a system Ax = b has a solution: Convert the augmented system (A, b) to a row reduced echelon matrix (AÊ¬,bÊ¬) and check whether bÊ¬ has no pivot. 
Remark: When writing a program implementing row reduction, we may stop when the last column of the matrix A is reached. In this case, the test whether the system Ax = b is 
7.11. RREF, FREE VARIABLES, HOMOGENEOUS SYSTEMS 
solvable is that the row-reduced matrix AÊ¬ has no zero row of index i>r such that bÊ¬ =0
i 
(where r is the number of pivots, and bÊ¬ is the row-reduced right-hand side). 
If wehavea homogeneous system Ax = 0, which means that b = 0, of course x = 0 is always a solution, but Theorem 7.16 implies that if the system Ax = 0 has more variables than equations, then it has some nonzero solution (we call it a nontrivial solution). 
Proposition 7.17. Given any homogeneous system Ax =0 of m equations in n variables, 
if m<n, then there is a nonzero vector x ¡Ê Rn such that Ax =0. Proof. Convert the matrix A to a reduced row echelon matrix AÊ¬ . We know that Ax = 0 i. AÊ¬x = 0. If r is the number of pivots of AÊ¬, we must have r ¡Ü m, so by Theorem 7.16 we may 
assign arbitrary values to n . r> 0 nonpivot variables and we get nontrivial solutions. 
Theorem 7.16 can also be used to characterize when a square matrix is invertible. First, note the following simple but important fact: 
If a square n ¡Á n matrix A is a row reduced echelon matrix, then either A is the identity or the bottom row of A is zero. 
Proposition 7.18. Let A be a square matrix of dimension n. The following conditions are equivalent: 
(a) The matrix A can be reduced to the identity by a sequence of elementary row operations. 

(b) The matrix A is a product of elementary matrices. 

(c) The matrix A is invertible. 

(d) The system of homogeneous equations Ax =0 has only the trivial solution x =0. 


Proof. First we prove that (a) implies (b). If (a) can be reduced to the identity by a sequence of row operations E1,...,Ep, this means that Ep ¡¤¡¤¡¤ E1A = I. Since each Ei is invertible, we get 
A = E1 .1 ¡¤¡¤¡¤ Ep .1 , 
where each Ei .1 is also an elementary row operation, so (b) holds. Now if (b) holds, since elementary row operations are invertible, A is invertible and (c) holds. If A is invertible, we already observed that the homogeneous system Ax = 0 has only the trivial solution x = 0, because from Ax = 0, we get A.1Ax = A.10; that is, x = 0. It remains to prove that (d) implies (a) and for this we prove the contrapositive: if (a) does not hold, then (d) does not hold. 
Using our basic observation about reducing square matrices, if A does not reduce to the identity, then A reduces to a row echelon matrix AÊ¬ whose bottom row is zero. Say AÊ¬ = PA, where P is a product of elementary row operations. Because the bottom row of AÊ¬ is zero, the system AÊ¬x = 0 has at most n . 1 nontrivial equations, and by Proposition 7.17, this system has a nontrivial solution x. But then, Ax = P .1AÊ¬x = 0 with x = 0, contradicting 
the fact that the system Ax = 0 is assumed to have only the trivial solution. Therefore, (d) implies (a) and the proof is complete. 
Proposition 7.18 yields a method for computing the inverse of an invertible matrix A: reduce A to the identity using elementary row operations, obtaining Ep ¡¤¡¤¡¤ E1A = I. Multiplying both sides by A.1 we get A.1 
= Ep ¡¤¡¤¡¤ E1. 
From a practical point of view, we can build up the product Ep ¡¤¡¤¡¤ E1 by reducing to row echelon form the augmented n ¡Á 2n matrix (A, In) obtained by adding the n columns of the identity matrix to A. This is just another way of performing the Gauss¨CJordan procedure. 
Here is an example: let us .nd the inverse of the matrix 
54 
A = . 
65 We form the 2 ¡Á 4 block matrix 
5410 
(A, I)= 
6501 and apply elementary row operations to reduce A to the identity. For example: 5410 5410 
(A, I)= .¡ú 
6501 11 .11 by subtracting row 1 from row 2, 5410 105 .4 
.¡ú 
11 .11 11 .11 by subtracting 4 ¡Á row 2 from row 1, 10 5 .4 105 .4 
.¡ú =(I,A.1),
11 .11 01 .65 by subtracting row 1 from row 2. Thus 5 .4 
A.1 
= .
.65 
Proposition 7.18 can also be used to give an elementary proof of the fact that if a square matrix A has a left inverse B (resp. a right inverse B), so that BA = I (resp. AB = I), then A is invertible and A.1 = B. This is an interesting exercise, try it! 

7.12. UNIQUENESS OF RREF 


7.12 Uniqueness of RREF Form 
For the sake of completeness, we prove that the reduced row echelon form of a matrix is unique. The neat proof given below is borrowed and adapted from W. Kahan. 
Proposition 7.19. Let A be any m ¡Á n matrix. If U and V are two reduced row echelon matrices obtained from A by applying two sequences of elementary row operations E1,...,Ep and F1,...,Fq, so that 
U = Ep ¡¤¡¤¡¤ E1A and V = Fq ¡¤¡¤¡¤ F1A, 
then U = V and Ep ¡¤¡¤¡¤ E1 = Fq ¡¤¡¤¡¤ F1. In other words, the reduced row echelon form of any matrix is unique. 
Proof. Let 
F .1
C = Ep ¡¤¡¤¡¤ E1F .1 ¡¤¡¤¡¤ 
1  q  
so that  
U = CV  and  V  = C.1U.  

We prove by induction on n that U = V (and C = I). 
Let ¦Éj denote the jth column of the identity matrix In, and let uj = U¦Éj, vj = V¦Éj, cj = C¦Éj, and aj = A¦Éj, be the jth column of U, V , C, and A respectively. 
First I claim that uj = 0 i. vj = 0 i. aj = 0. 
Indeed, if vj = 0, then (because U = CV ) uj = Cvj = 0, and if uj = 0, then vj = 
C.1
uj = 0. Since U = Ep ¡¤¡¤¡¤ E1A, we also get aj = 0 i. uj = 0. 
Therefore, we may simplify our task by striking out columns of zeros from U, V , and A, since they will have corresponding indices. We still use n to denote the number of columns of 
A. Observe that because U and V are reduced row echelon matrices with no zero columns, we must have u1 = v1 = ¦É1. 
Claim. If U and V are reduced row echelon matrices without zero columns such that U = CV , for all k ¡Ý 1, if k ¡Ü n, then ¦Ék occurs in U i. ¦Ék occurs in V , and if ¦Ék does occur in U, then 
1. 
¦Ék occurs for the same column index jk in both U and V ; 

2. 
the .rst jk columns of U and V match; 

3. 
the subsequent columns in U and V (of column index >jk) whose coordinates of index k + 1 through m are all equal to 0 also match. Let nk be the rightmost index of such a column, with nk = jk if there is none. 

4. 
the .rst nk columns of C match the .rst nk columns of In. 


We prove this claim by induction on k. 
For the base case k = 1, we already know that u1 = v1 = ¦É1. We also have 
c1 = C¦É1 = Cv1 = u1 = ¦É1. 
If vj = ¦Ë¦É1 for some ¦Ë ¡Ê R, then 
uj = U¦Éj = CV ¦Éj = Cvj = ¦ËC¦É1 = ¦Ëc1 = ¦Ë¦É1 = vj. 
A similar argument using C.1 shows that if uj = ¦Ë¦É1, then vj = uj. Therefore, all the columns of U and V proportional to ¦É1 match, which establishes the base case. Observe that if ¦É2 appears in U, then it must appear in both U and V for the same index, and if not then n1 = n and U = V . 
Next us now prove the induction step. If nk = n, then U = V and we are done. Otherwise, ¦Ék+1 appears in both U and V , in which case, by (2) and (3) of the induction hypothesis, it appears in both U and V for the same index, say jk+1. Thus, ujk+1 It follows 
= vjk+1 = ¦Ék+1. 
that 
ck+1 = C¦Ék+1 = Cvjk+1 = ujk+1 = ¦Ék+1, so the .rst jk+1 columns of C match the .rst jk+1 columns of In. Consider any subsequent column vj (with j>jk+1) whose elements beyond the (k + 1)th all vanish. Then vj is a linear combination of columns of V to the left of vj, so 
uj = Cvj = vj. 
because the .rst k + 1 columns of C match the .rst column of In. Similarly, any subsequent column uj (with j>jk+1) whose elements beyond the (k + 1)th all vanish is equal to vj. Therefore, all the subsequent columns in U and V (of index >jk+1) whose elements beyond the (k + 1)th all vanish also match, so the .rst nk+1 columns of C match the .rst nk+1 columns of C, which completes the induction hypothesis. 
We can now prove that U = V (recall that we may assume that U and V have no zero columns). We noted earlier that u1 = v1 = ¦É1, so there is a largest k ¡Ü n such that ¦Ék occurs in U. Then the previous claim implies that all the columns of U and V match, which means that U = V . 
The reduction to row echelon form also provides a method to describe the set of solutions of a linear system of the form Ax = b. 

7.13 Solving Linear Systems Using RREF 
First we have the following simple result. 
Proposition 7.20. Let A be any m ¡Á n matrix and let b ¡Ê Rm be any vector. If the system Ax = b has a solution, then the set Z of all solutions of this system is the set 
Z = x0 + Ker (A)= {x0 + x | Ax =0}, 
7.13. SOLVING LINEAR SYSTEMS USING RREF 
where x0 ¡Ê Rn is any solution of the system Ax = b, which means that Ax0 = b (x0 is called a special solution), and where Ker (A)= {x ¡Ê Rn | Ax =0}, the set of solutions of the homogeneous system associated with Ax = b. 
Proof. Assume that the system Ax = b is solvable and let x0 and x1 be any two solutions so that Ax0 = b and Ax1 = b. Subtracting the .rst equation from the second, we get 
A(x1 . x0)=0, 
which means that x1 . x0 ¡Ê Ker (A). Therefore, Z . x0 + Ker (A), where x0 is a special solution of Ax = b. Conversely, if Ax0 = b, then for any z ¡Ê Ker (A), we have Az = 0, and so 
A(x0 + z)= Ax0 + Az = b +0= b, 
which shows that x0 + Ker (A) . Z. Therefore, Z = x0 + Ker (A). 
Given a linear system Ax = b, reduce the augmented matrix (A, b) to its row echelon form (AÊ¬,bÊ¬). As we showed before, the system Ax = b has a solution i. bÊ¬ contains no pivot. Assume that this is the case. Then, if (AÊ¬,bÊ¬) has r pivots, which means that AÊ¬ has r pivots since bÊ¬ has no pivot, we know that the .rst r columns of Im appear in AÊ¬ . 
We can permute the columns of AÊ¬ and renumber the variables in x correspondingly so that the .rst r columns of Im match the .rst r columns of AÊ¬, and then our reduced echelon matrix is of the form (R, bÊ¬) with 
Ir F 
R = 0m.r,r 0m.r,n.r 
and bÊ¬ = d,0m.r 
where F is a r ¡Á (n . r) matrix and d ¡Ê Rr . Note that R has m . r zero rows. 
Then because  
Ir 0m.r,r  F 0m.r,n.r  d 0n.r  =  d 0m.r  = bÊ¬ ,  
we see that  
d  
x0 =  

0n.r 
is a special solution of Rx = bÊ¬, and thus to Ax = b. In other words, we get a special solution by assigning the .rst r components of bÊ¬ to the pivot variables and setting the nonpivot variables (the free variables) to zero. 
Here is an example of the preceding construction taken from Kumpel and Thorpe [40]. The linear system 
x1 . x2 + x3 + x4 . 2x5 = .1 .2x1 +2x2 . x3 + x5 =2 x1 . x2 +2x3 +3x4 . 5x5 = .1, 
is represented by the augmented matrix 
.. 
1 .111 .2 .1 (A, b)= ..22 .101 2 ., 1 .123 .5 .1 
where A is a 3 ¡Á 5 matrix. The reader should .nd that the row echelon form of this system is 
.. 
(AÊ¬,bÊ¬)= .01 .01 10 .21 .13 .01 . . 
0000 0 0 The 3 ¡Á 5 matrix AÊ¬ has rank 2. We permute the second and third columns (which is 
equivalent to interchanging variables x2 and x3) to form 
I2 F .1 .11 
R = ,F = . 
01,2 01,3 02 .3 
Then a special solution to this linear system is given by 
.. 
.1 
d 
..
x0 = =0 . 
03 03 
We can also .nd a basis of the kernel (nullspace) of A using F . If x =(u, v) is in the kernel of A, with u ¡Ê Rr and v ¡Ê Rn.r, then x is also in the kernel of R, which means that Rx = 0; that is, 
Ir F uu + Fv 0r 
== . 0m.r,r 0m.r,n.r v 0m.r 0m.r 
Therefore, u = .Fv, and Ker (A) consists of all vectors of the form 
.Fv .F 
= v, 
vIn.r 
for any arbitrary v ¡Ê Rn.r . It follows that the n . r columns of the matrix 
.F 
N = 
In.r 
form a basis of the kernel of A. This is because N contains the identity matrix In.r as a submatrix, so the columns of N are linearly independent. In summary, if N1,...,Nn.r are the columns of N, then the general solution of the equation Ax = b is given by 
Nn.r 
x = d + xr+1N1 + ¡¤¡¤¡¤ + xn,0n.r 

7.13. SOLVING LINEAR SYSTEMS USING RREF 
where xr+1,...,xn are the free variables; that is, the nonpivot variables. Going back to our example from Kumpel and Thorpe [40], we see that 
.
. 
11 .1 0 .2 .3 
.F 
N = =10 0 ,
I3 
01 0 00 1 
..... 
and that the general solution is given by 
..... 
..
..
..
.. 
..... 
..... 
..... 
..... 
.11 1 .1 
..... 
00 .2 .3
..... 
x =0+ x3 1+ x4 0+ x5 0 . 
001 0 
000 1 
..... 
..... 
In the general case where the columns corresponding to pivots are mixed with the columns 
corresponding to free variables, we .nd the special solution as follows. Let i1 < ¡¤¡¤¡¤ <ir be the indices of the columns corresponding to pivots. Assign bÊ¬ to the pivot variable 
xik for k =1,...,r, and set all other variables to 0. To .nd a bkasis of the kernel, we form the n . r vectors Nk obtained as follows. Let j1 < ¡¤¡¤¡¤ <jn.r be the indices of the columns corresponding to free variables. For every column jk corresponding to a free variable (1 ¡Ü k ¡Ü n . r), form the vector Nk de.ned so that the entries Nik 1 ,...,Nk are equal to the negatives of the .rst r entries in column jk (.ip the sign of these entries);irlet Njk k = 1, and set all other entries to zero. Schematically, if the column of index jk (corresponding to the 
free variable xjk ) is 
.
. 

........ 
¦Á1 
. 
. 
. 
¦Ár 
,
0 
. 
. 
. 
........ 

0 

is given by 

.
. 

10 

. 

. 

. 

. 

i1 . 1 
i1 
i1 +1 
. 

. 

. 

ir . 1 
ir 
ir +1 
. 

. 

. 

jk . 1 
jk 
jk +1 
. 

. 

............................ 

. 

. 

0 

.¦Á1 
0 

. 

. 

. 

0 

.¦Ár 
0 

. 

. 

. 

0 

1 

0 

. 

. 

............................ 

. 

. 

. 

n 0 
The presence of the 1 in position jk guarantees that N1,...,Nn.r are linearly indepen-dent. As an illustration of the above method, consider the problem of .nding a basis of the subspace V of n ¡Á n matrices A ¡Ê Mn(R) satisfying the following properties: 
1. 
The sum of the entries in every row has the same value (say c1); 

2. 
The sum of the entries in every column has the same value (say c2). 


It turns out that c1 = c2 and that the 2n.2 equations corresponding to the above conditions are linearly independent. We leave the proof of these facts as an interesting exercise. It can be shown using the duality theorem (Theorem 10.4) that the dimension of the space V of matrices satisying the above equations is n2 . (2n . 2). Let us consider the case n = 4. There are 6 equations, and the space V has dimension 10. The equations are 
a11 + a12 + a13 + a14 . a21 . a22 . a23 . a24 =0 a21 + a22 + a23 + a24 . a31 . a32 . a33 . a34 =0 a31 + a32 + a33 + a34 . a41 . a42 . a43 . a44 =0 a11 + a21 + a31 + a41 . a12 . a22 . a32 . a42 =0 a12 + a22 + a32 + a42 . a13 . a23 . a33 . a43 =0 a13 + a23 + a33 + a43 . a14 . a24 . a34 . a44 =0, 

7.13. SOLVING LINEAR SYSTEMS USING RREF 
and the corresponding matrix is 
.
. 
A = 

....... 

11 1 1 .1 .1 .1 .10 0 0 0 0 0 0 0 
00 0 0 1 1 1 1 .1 .1 .1 .10 0 0 0 
00 0 0 0 0 0 0 1 1 1 1 .1 .1 .1 .1 

1  .1  0  0  1  .1  0  0  1  .1  0  0  1  .1  0  0  
0  1  .1  0  0  1  .1  0  0  1  .1  0  0  1  .1  0  
0  0  1  .1  0  0  1  .1  0  0  1  .1  0  0  1  .1  

....... 

. 

The result of performing the reduction to row echelon form yields the following matrix in rref: 
.
. 
U = 

....... 

1  0  0  0  0  .1  .1  .1  0  .1  .1  .1  2  1  1  1  
0  1  0  0  0  1  0  0  0  1  0  0  .1  0  .1  .1  
0  0  1  0  0  0  1  0  0  0  1  0  .1  .1  0  .1  

0  0  0  1  0  0  0  1  0  0  0  1  .1  .1  .1  0  
0  0  0  0  1  1  1  1  0  0  0  0  .1  .1  .1  .1  
0  0  0  0  0  0  0  0  1  1  1  1  .1  .1  .1  .1  

....... 

The list pivlist of indices of the pivot variables and the list freelist of indices of the free variables is given by 
pivlist = (1, 2, 3, 4, 5, 9), freelist =(6, 7, 8, 10, 11, 12, 13, 14, 15, 16). 
After applying the algorithm to .nd a basis of the kernel of U, we .nd the following 16 ¡Á 10 
matrix

.
. 
BK = 

........................... 

1  1  1  1  1  1  .2  .1  .1  .1  
.1  0  0  .1  0  0  1  0  1  1  
0  .1  0  0  .1  0  1  1  0  1  
0  0  .1  0  0  .1  1  1  1  0  
.1  .1  .1  0  0  0  1  1  1  1  
1  0  0  0  0  0  0  0  0  0  
0  1  0  0  0  0  0  0  0  0  
0  0  1  0  0  0  0  0  0  0  

0  0  0  .1  .1  .1  1  1  1  1  
0  0  0  1  0  0  0  0  0  0  
0  0  0  0  1  0  0  0  0  0  
0  0  0  0  0  1  0  0  0  0  
0  0  0  0  0  0  1  0  0  0  
0  0  0  0  0  0  0  1  0  0  
0  0  0  0  0  0  0  0  1  0  
0  0  0  0  0  0  0  0  0  1  

........................... 

. 

The reader should check that that in each column j of BK, the lowest bold 1 belongs to the row whose index is the jth element in freelist, and that in each column j of BK, the signs of the entries whose indices belong to pivlist are the .ipped signs of the 6 entries in the column U corresponding to the jth index in freelist. We can now read o. from BK the 4 ¡Á 4 matrices that form a basis of V : every column of BK corresponds to a matrix whose rows have been concatenated. We get the following 10 matrices: 
.
.
.
.
.
. 
1 .100 10 .10 100 .1 

... 

.1 

1 00 

0 000 

...

, 

M2 
= 

... 

.1 

010 

0000 

...

, 

M3 
= 

... 

.1 

00 1 

0 00 0 

...

M1 
= 

, 

0 000 0000 0000 

.
.
.
.
.
. 
1 .100 10 .10 100 .1 

... 

0 000 

.1 

1 00 

...

, 

M5 
= 

... 

0000 

.1 

010 

...

, 

M6 
= 

... 

0 00 0 

.1 

00 1 

...

M4 
= 

, 

0 000 0000 0000 

.
.
.
.
.
. 
...
.2111 .1011 .1101 0
... 

...

... 

...

... 

1 000 

1 000 

1 00 

M7 
M8 M9 
= 

= 

=

, 

, 

,

1 000 

0 1 000 0 100 0 010 
1 000 

1 00 

.
. 
... 

.1110 
1 000 

1 000 
0 001 

...

M10 = 
. 

Recall that a magic square is a square matrix that satis.es the two conditions about the sum of the entries in each row and in each column to be the same number, and also the additional two constraints that the main descending and the main ascending diagonals add up to this common number. Furthermore, the entries are also required to be positive integers. For n = 4, the additional two equations are 
a22 + a33 + a44 . a12 . a13 . a14 =0 a41 + a32 + a23 . a11 . a12 . a13 =0, 
and the 8 equations stating that a matrix is a magic square are linearly independent. Again, by running row elimination, we get a basis of the ¡°generalized magic squares¡± whose entries are not restricted to be positive integers. We .nd a basis of 8 matrices. For n = 3, we .nd a basis of 3 matrices. 
A magic square is said to be normal if its entries are precisely the integers 1, 2 ...,n2 . Then since the sum of these entries is 
n2(n2 + 1) 
1+2+3+ ¡¤¡¤¡¤ + n 2 = ,
2 

7.14. ELEMENTARY MATRICES AND COLUMNS OPERATIONS 
and since each row (and column) sums to the same number, this common value (the magic sum) is 
n(n2 + 1) 
. 
2 It is easy to see that there are no normal magic squares for n = 2. For n = 3, the magic sum is 15, for n = 4, it is 34, and for n = 5, it is 65. In the case n = 3, we have the additional condition that the rows and columns add up to 15, so we end up with a solution parametrized by two numbers x1,x2; namely, 
.  .  
x1 + x2 . 5  10 . x2  10 . x1  
.20 . 2x1 . x2  5  2x1 + x2 . 10. .  
x1  x2  15 . x1 . x2  

Thus, in order to .nd a normal magic square, we have the additional inequality constraints 
x1 + x2 > 5 
x1 < 10 
x2 < 10 
2x1 + x2 < 20 
2x1 + x2 > 10 
x1 > 0 
x2 > 0 
x1 + x2 < 15, 
and all 9 entries in the matrix must be distinct. After a tedious case analysis, we discover the remarkable fact that there is a unique normal magic square (up to rotations and re.ections): 
.. 
276 
..
951 . 438 
It turns out that there are 880 di.erent normal magic squares for n = 4, and 275, 305, 224 normal magic squares for n = 5 (up to rotations and re.ections). Even for n = 4, it takes a fair amount of work to enumerate them all! Finding the number of magic squares for n> 5 is an open problem! 


7.14 Elementary Matrices and Columns Operations 
Instead of performing elementary row operations on a matrix A, we can perform elementary columns operations, which means that we multiply A by elementary matrices on the right. As elementary row and column operations, P (i, k), Ei,j;¦Â, Ei,¦Ë perform the following actions: 
1. 
As a row operation, P (i, k) permutes row i and row k. 

2. 
As a column operation, P (i, k) permutes column i and column k. 

3. 
The inverse of P (i, k) is P (i, k) itself. 

4. 
As a row operation, Ei,j;¦Â adds ¦Â times row j to row i. 

5. 
As a column operation, Ei,j;¦Â adds ¦Â times column i to column j (note the switch in the indices). 

6. 
The inverse of Ei,j;¦Â is Ei,j;.¦Â. 

7. 
As a row operation, Ei,¦Ë multiplies row i by ¦Ë. 

8. 
As a column operation, Ei,¦Ë multiplies column i by ¦Ë. 

9. 
The inverse of Ei,¦Ë is Ei,¦Ë.1 . 


We can de.ne the notion of a reduced column echelon matrix and show that every matrix can be reduced to a unique reduced column echelon form. Now given any m ¡Á n matrix A, if we .rst convert A to its reduced row echelon form R, it is easy to see that we can apply elementary column operations that will reduce R to a matrix of the form 
Ir 0r,n.r 
,
0m.r,r 0m.r,n.r 
where r is the number of pivots (obtained during the row reduction). Therefore, for every m ¡Á n matrix A, there exist two sequences of elementary matrices E1,...,Ep and F1,...,Fq, such that 
Ir 0r,n.r
Ep ¡¤¡¤¡¤ E1AF1 ¡¤¡¤¡¤ Fq = . 
0m.r,r 0m.r,n.r 
The matrix on the right-hand side is called the rank normal form of A. Clearly, r is the rank of A. As a corollary we obtain the following important result whose proof is immediate. 
Proposition 7.21. A matrix A and its transpose AT have the same rank. 

7.15 Transvections and Dilatations ¢à 
In this section we characterize the linear isomorphisms of a vector space E that leave every vector in some hyperplane .xed. These maps turn out to be the linear maps that are represented in some suitable basis by elementary matrices of the form Ei,j;¦Â (transvections) or Ei,¦Ë (dilatations). Furthermore, the transvections generate the group SL(E), and the dilatations generate the group GL(E). 
7.15. TRANSVECTIONS AND DILATATIONS ¢à 
Let H be any hyperplane in E, and pick some (nonzero) vector v ¡Ê E such that v/¡Ê H, so that E = H ¨’ Kv. 
Assume that f : E ¡ú E is a linear isomorphism such that f(u)= u for all u ¡Ê H, and that f is not the identity. We have 
f(v)= h + ¦Áv, for some h ¡Ê H and some ¦Á ¡Ê K, 
with ¦Á = 0, because otherwise we would have f(v)= h = f(h) since h ¡Ê H, contradicting the injectivity of f (v = h since v/¡Ê H). For any x ¡Ê E, if we write 
x = y + tv, for some y ¡Ê H and some t ¡Ê K, 
then f(x)= f(y)+ f(tv)= y + tf(v)= y + th + t¦Áv, 
and since ¦Áx = ¦Áy + t¦Áv, we get 
f(x) . ¦Áx = (1 . ¦Á)y + th f(x) . x = t(h +(¦Á . 1)v). 
Observe that if E is .nite-dimensional, by picking a basis of E consisting of v and basis vectors of H, then the matrix of f is a lower triangular matrix whose diagonal entries are all 1 except the .rst entry which is equal to ¦Á. Therefore, det(f)= ¦Á. 
Case 1 . ¦Á = 1. 
We have f(x)= ¦Áx i. (1 . ¦Á)y + th = 0 i. 
t 
y = h. 
¦Á . 1 
Then if we let w = h +(¦Á . 1)v, for y =(t/(¦Á . 1))h, we have 
tt t 
x = y + tv = h + tv =(h +(¦Á . 1)v)= w, 
¦Á . 1 ¦Á . 1¦Á . 1
which shows that f(x)= ¦Áx i. x ¡Ê Kw. Note that w ¡Ê/H, since ¦Á = 1 and v ¡Ê/H. Therefore, 
E = H ¨’ Kw, 
and f is the identity on H and a magni.cation by ¦Á on the line D = Kw. 
De.nition 7.6. Given a vector space E, for any hyperplane H in E, any nonzero vector u ¡Ê E such that u ¡Ê H, and any scalar ¦Á =0, 1, a linear map f such that f(x)= x for all x ¡Ê H and f(x)= ¦Áx for every x ¡Ê D = Ku is called a dilatation of hyperplane H, direction D, and scale factor ¦Á. 
If ¦ÐH and ¦ÐD are the projections of E onto H and D, then we have 
f(x)= ¦ÐH (x)+ ¦Á¦ÐD(x). 
The inverse of f is given by 
f.1(x)= ¦ÐH (x)+ ¦Á.1¦ÐD(x). 
When ¦Á = .1, we have f2 = id, and f is a symmetry about the hyperplane H in the 
direction D. This situation includes orthogonal re.ections about H. 
Case 2 . ¦Á = 1. 
In this case, 
f(x) . x = th, 
that is, f(x) . x ¡Ê Kh for all x ¡Ê E. Assume that the hyperplane H is given as the kernel of some linear form ., and let a = .(v). We have a = 0, since v ¡Ê/H. For any x ¡Ê E, we have 
.(x . a .1.(x)v)= .(x) . a .1.(x).(v)= .(x) . .(x)=0, 
which shows that x . a.1.(x)v ¡Ê H for all x ¡Ê E. Since every vector in H is .xed by f, we get 
x . a .1.(x)v = f(x . a .1.(x)v) = f(x) . a .1.(x)f(v), 
so 
.1 .1
f(x)= x + .(x)(f(av) . av). Since f(z) . z ¡Ê Kh for all z ¡Ê E, we conclude that u = f(a.1v) . a.1v = ¦Âh for some ¦Â ¡Ê K, so .(u) = 0, and we have 
f(x)= x + .(x)u, .(u)=0. (.) 
A linear map de.ned as above is denoted by ¦Ó.,u. 
Conversely for any linear map f = ¦Ó.,u given by Equation (.), where . is a nonzero linear form and u is some vector u ¡Ê E such that .(u)=0, if u = 0 , then f is the identity, so assume that u = 0. If so, we have f(x)= x i. .(x) = 0, that is, i. x ¡Ê H. We also claim that the inverse of f is obtained by changing u to .u. Actually, we check the slightly more general fact that 
¦Ó.,u . ¦Ó.,w = ¦Ó.,u+w. Indeed, using the fact that .(w) = 0, we have 
¦Ó.,u(¦Ó.,w(x)) = ¦Ó.,w(x)+ .(¦Ó.,w(x))u 
= ¦Ó.,w(x)+(.(x)+ .(x).(w))u 
= ¦Ó.,w(x)+ .(x)u 
= x + .(x)w + .(x)u 
= x + .(x)(u + w). 

7.15. TRANSVECTIONS AND DILATATIONS ¢à 
For v = .u, we have ¦Ó.,u+v = ..,0 = id, so ¦Ó.1 = ¦Ó.,.u, as claimed. 
.,u 
Therefore, we proved that every linear isomorphism of E that leaves every vector in some hyperplane H .xed and has the property that f(x) . x ¡Ê H for all x ¡Ê E is given by a map ¦Ó.,u as de.ned by Equation (.), where . is some nonzero linear form de.ning H and u is some vector in H. We have ¦Ó.,u = id i. u = 0. 
De.nition 7.7. Given any hyperplane H in E, for any nonzero nonlinear form . ¡Ê E. de.ning H (which means that H = Ker (.)) and any nonzero vector u ¡Ê H, the linear map f = ¦Ó.,u given by 
¦Ó.,u(x)= x + .(x)u, .(u)=0, 
for all x ¡Ê E is called a transvection of hyperplane H and direction u. The map f = ¦Ó.,u leaves every vector in H .xed, and f(x) . x ¡Ê Ku for all x ¡Ê E. 
The above arguments show the following result. 
Proposition 7.22. Let f : E ¡ú E be a bijective linear map and assume that f = id and that f(x)= x for all x ¡Ê H, where H is some hyperplane in E. If there is some nonzero vector u ¡Ê E such that u/¡Ê H and f(u) . u ¡Ê H, then f is a transvection of hyperplane H; otherwise, f is a dilatation of hyperplane H. 
Proof. Using the notation as above, for some v ¡Ê/H, we have f(v)= h + ¦Áv with ¦Á = 0, and write u = y + tv with y ¡Ê H and t = 0 since u/¡Ê H. If f(u) . u ¡Ê H, from 
f(u) . u = t(h +(¦Á . 1)v), 
we get (¦Á . 1)v ¡Ê H, and since v ¡Ê/H, we must have ¦Á = 1, and we proved that f is a transvection. Otherwise, ¦Á =0, 1, and we proved that f is a dilatation. 
If E is .nite-dimensional, then ¦Á = det(f), so we also have the following result. 
Proposition 7.23. Let f : E ¡ú E be a bijective linear map of a .nite-dimensional vector space E and assume that f = id and that f(x)= x for all x ¡Ê H, where H is some hyperplane in E. If det(f)=1, then f is a transvection of hyperplane H; otherwise, f is a dilatation of hyperplane H. 
Suppose that f is a dilatation of hyperplane H and direction u, and say det(f)= ¦Á =0, 1. Pick a basis (u, e2,...,en) of E where (e2,...,en) is a basis of H. Then the matrix of f is 
of the form

.
. 
.... 

¦Á 0 ¡¤¡¤¡¤ 0 
01 0 

..
.
. ..
. 00 ¡¤¡¤¡¤ 1 .. 
....

, 

which is an elementary matrix of the form E1,¦Á. Conversely, it is clear that every elementary matrix of the form Ei,¦Á with ¦Á =0, 1 is a dilatation. 
Now, assume that f is a transvection of hyperplane H and direction u ¡Ê H. Pick some v/¡Ê H, and pick some basis (u, e3,...,en) of H, so that (v, u, e3,...,en) is a basis of E. Since f(v) . v ¡Ê Ku, the matrix of f is of the form 
.
. 
.... 

10 ¡¤¡¤¡¤ 0 
¦Á 10 

..
.
. ..
. 00 ¡¤¡¤¡¤ 1 .. 
....

, 

which is an elementary matrix of the form E2,1;¦Á. Conversely, it is clear that every elementary matrix of the form Ei,j;¦Á (¦Á = 0) is a transvection. 
The following proposition is an interesting exercise that requires good mastery of the elementary row operations Ei,j;¦Â; see Problems 7.10 and 7.11. 
Proposition 7.24. Given any invertible n ¡Á n matrix A, there is a matrix S such that 
In.1 0 
SA == En,¦Á,
0 ¦Á 
with ¦Á = det(A), and where S is a product of elementary matrices of the form Ei,j;¦Â; that is, S is a composition of transvections. 
Surprisingly, every transvection is the composition of two dilatations! 
Proposition 7.25. If the .eld K is not of characteristic 2, then every transvection f of hyperplane H can be written as f = d2 . d1, where d1,d2 are dilatations of hyperplane H, where the direction of d1 can be chosen arbitrarily. 
Proof. Pick some dilatation d1 of hyperplane H and scale factor ¦Á =0, 1. Then, d2 = f .d.11 leaves every vector in H .xed, and det(d2)= ¦Á.1 = 1. By Proposition 7.23, the linear map d2 is a dilatation of hyperplane H, and we have f = d2 . d1, as claimed. 
Observe that in Proposition 7.25, we can pick ¦Á = .1; that is, every transvection of hyperplane H is the compositions of two symmetries about the hyperplane H, one of which can be picked arbitrarily. 
Remark: Proposition 7.25 holds as long as K = {0, 1}. The following important result is now obtained. 
Theorem 7.26. Let E be any .nite-dimensional vector space over a .eld K of characteristic not equal to 2. Then the group SL(E) is generated by the transvections, and the group GL(E) is generated by the dilatations. 

7.15. TRANSVECTIONS AND DILATATIONS ¢à 
Proof. Consider any f ¡Ê SL(E), and let A be its matrix in any basis. By Proposition 7.24, there is a matrix S such that 
In.1 0 
SA == En,¦Á,
0 ¦Á 
with ¦Á = det(A), and where S is a product of elementary matrices of the form Ei,j;¦Â. Since det(A) = 1, we have ¦Á = 1, and the result is proven. Otherwise, if f is invertible but f/¡Ê SL(E), the above equation shows En,¦Á is a dilatation, S is a product of transvections, and by Proposition 7.25, every transvection is the composition of two dilatations. Thus, the second result is also proven. 
We conclude this section by proving that any two transvections are conjugate in GL(E). Let ¦Ó.,u (u = 0) be a transvection and let g ¡Ê GL(E) be any invertible linear map. We have 
(g . ¦Ó.,u . g .1)(x)= g(g .1(x)+ .(g .1(x))u) 
= x + .(g .1(x))g(u). 
Let us .nd the hyperplane determined by the linear form x ¡ú .(g.1(x)). This is the set of 
vectors x ¡Ê E such that .(g.1(x)) = 0, which holds i. g.1(x) ¡Ê H i. x ¡Ê g(H). Therefore, Ker (..g.1)= g(H)= HÊ¬, and we have g(u) ¡Ê g(H)= HÊ¬, so g.¦Ó.,u .g.1 is the transvection of hyperplane HÊ¬ = g(H) and direction uÊ¬ = g(u) (with uÊ¬ ¡Ê HÊ¬). 
Conversely, let ¦Ó¦×,u£§ be some transvection (uÊ¬ = 0). Pick some vectors v, vÊ¬ such that .(v)= ¦×(vÊ¬) = 1, so that E = H ¨’ Kv = HÊ¬ ¨’ KvÊ¬ . 
There is a linear map g ¡Ê GL(E) such that g(u)= uÊ¬ , g(v)= vÊ¬, and g(H)= HÊ¬ . To 
de.ne g, pick a basis (v, u, e2,...,en.1) where (u, e2,...,en.1) is a basis of H and pick a basis (vÊ¬,uÊ¬,e2Ê¬ ,...,eÊ¬ n.1) where (uÊ¬,e2Ê¬ ,...,eÊ¬ n.1) is a basis of HÊ¬; then g is de.ned so that 
g(v)= vÊ¬ , g(u)= uÊ¬, and g(ei)= g(eiÊ¬ ), for i =2,...,n . 1. If n = 2, then ei and ei Ê¬ are 
missing. Then, we have 
(g . ¦Ó.,u . g .1)(x)= x + .(g .1(x))uÊ¬ . 
Now . . g.1 also determines the hyperplane HÊ¬ = g(H), so we have . . g.1 = ¦Ë¦× for some nonzero ¦Ë in K. Since vÊ¬ = g(v), we get 
.(v)= . . g .1(vÊ¬)= ¦Ë¦×(vÊ¬), 
and since .(v)= ¦×(vÊ¬) = 1, we must have ¦Ë = 1. It follows that 
(g . ¦Ó.,u . g .1)(x)= x + ¦×(x)uÊ¬ = ¦Ó¦×,u£§ (x). 
In summary, we proved almost all parts the following result. 
Proposition 7.27. Let E be any .nite-dimensional vector space. For every transvection ¦Ó.,u (u =0) and every linear map g ¡Ê GL(E), the map g . ¦Ó.,u . g.1 is the transvection of hyperplane g(H) and direction g(u) (that is, g . ¦Ó.,u . g.1 = ¦Ó..g.1,g(u)). For every other 
transvection ¦Ó¦×,u£§ (uÊ¬ =0), there is some g ¡Ê GL(E) such ¦Ó¦×,u£§ = g . ¦Ó.,u . g.1; in other words any two transvections (= id) are conjugate in GL(E). Moreover, if n ¡Ý 3, then the linear isomorphism g as above can be chosen so that g ¡Ê SL(E). 
Proof. We just need to prove that if n ¡Ý 3, then for any two transvections ¦Ó.,u and ¦Ó¦×,u£§ (u, uÊ¬ = 0), there is some g ¡Ê SL(E) such that ¦Ó¦×,u£§ = g .¦Ó.,u .g.1 . As before, we pick a basis (v, u, e2,...,en.1) where (u, e2,...,en.1) is a basis of H, we pick a basis (vÊ¬,uÊ¬,eÊ¬ 2,...,eÊ¬ )
n.1
where (uÊ¬,e2Ê¬ ,...,enÊ¬.1) is a basis of HÊ¬, and we de.ne g as the unique linear map such that g(v)= vÊ¬ , g(u)= uÊ¬, and g(ei)= eÊ¬ , for i =1,...,n . 1. But in this case, both H and HÊ¬ = g(H) have dimension at least 2,iso in any basis of HÊ¬ including uÊ¬, there is some basis vector eÊ¬ independent of uÊ¬, and we can rescale eÊ¬ in such a way that the matrix of g over 
2	2 
the two bases has determinant +1. 


7.16 Summary 
The main concepts and results of this chapter are listed below: 
. 	One does not solve (large) linear systems by computing determinants. 

. 	Upper-triangular (lower-triangular) matrices. 

. 	Solving by back-substitution (forward-substitution). 

. 	Gaussian elimination. 

. 	Permuting rows. 

. 	The pivot of an elimination step; pivoting. 

. 	Transposition matrix; elementary matrix. 

. 	The Gaussian elimination theorem (Theorem 7.1). 

. 	Gauss-Jordan factorization. 

. 	LU-factorization; Necessary and su.cient condition for the existence of an 
LU-factorization (Proposition 7.2). 


. 	LDU-factorization. 

. 	¡°PA = LU theorem¡± (Theorem 7.5). 

. 	LDLT-factorization of a symmetric matrix. 

7.16. SUMMARY 

. 	
Avoiding small pivots: partial pivoting; complete pivoting. 

. 	
Gaussian elimination of tridiagonal matrices. 

. 	
LU-factorization of tridiagonal matrices. 

. 	
Symmetric positive de.nite matrices (SPD matrices). 

. 	
Cholesky factorization (Theorem 7.10). 

. 	
Criteria for a symmetric matrix to be positive de.nite; Sylvester¡¯s criterion. 

. 	
Reduced row echelon form. 

. 	
Reduction of a rectangular matrix to its row echelon form. 

. 	
Using the reduction to row echelon form to decide whether a system Ax = b is solvable, and to .nd its solutions, using a special solution and a basis of the homogeneous system Ax = 0. 

. 	
Magic squares. 

. 	
Transvections and dilatations. 


7.17 Problems 
Problem 7.1. Solve the following linear systems by Gaussian elimination: 
..
...
......
. 
231 x 6 111 x 6 

.

12 .1

.
.

y

.

=

.

2

.

,

.

112

.
.

y

.

=

.

9

.. 

.3 .51 z .7 123 z 14 Problem 7.2. Solve the following linear system by Gaussian elimination: 
...
..
. 
1 211 x1 7 
... 

2 323 

.101 .1 

... 
... 

x2 
x3 
... 

= 

... 

14 

.1 

...

. 

.2 .14 0 x4 2 Problem 7.3. Consider the matrix 
.
. 
A =

. 

1 c 0 
241

.

. 

351 
When applying Gaussian elimination, which value of c yields zero in the second pivot posi-tion? Which value of c yields zero in the third pivot position? In this case, what can your say about the matrix A? 
...
.. 
Problem7.4. Solvethesystem . 
2110 x1 1 
... 

4331 

8795 

... 
... 

x2 
x3 
... 

= 

... 

.1 

.1 

... 

6798 x4 1 using the LU-factorization of Example 7.1. Problem 7.5. Apply rref to the matrix 
.
. 
A2 = 
... 

1 211 
2 323 

.101 .1 
.2 .13 0 

...

. 

Problem 7.6. Apply rref to the matrix 

.
. 
... 

1 4 916 
4 9 1625 

9 16 25 36 
16 25 36 49 

...

. 

7.17. PROBLEMS 
Problem 7.7. (1) Prove that the dimension of the subspace of 2 ¡Á 2 matrices A, such that the sum of the entries of every row is the same (say c1) and the sum of entries of every column is the same (say c2) is 2. 
(2) Prove that the dimension of the subspace of 2 ¡Á 2 matrices A, such that the sum of the entries of every row is the same (say c1), the sum of entries of every column is the same (say c2), and c1 = c2 is also 2. Prove that every such matrix is of the form 
ab 
,
ba 
and give a basis for this subspace. 
(3) Prove that the dimension of the subspace of 3 ¡Á 3 matrices A, such that the sum of the entries of every row is the same (say c1), the sum of entries of every column is the same (say c2), and c1 = c2 is 5. Begin by showing that the above constraints are given by the set of equations 
.
. 

............. 

a11 a12 a13 
a21 
a22 
a23 
a31 a32 a33 
............. 

= 

.
. 
11 1 .1 .1 .10 0 0 

.
. 

0 

..... 

00 0 1 1 1 .1 .1 .1 

1 .10 1 .10 1 .10 

01 .10 1 .10 1 .1 

..... 

..... 

0 

0 

0 

..... 

. 

01 1 .10 0 .10 0 

0 

Prove that every matrix satisfying the above constraints is of the form 

.
. 
. 

a + b . c .a + c + e .b + c + d .a . b + c + d + ea b
.

, 

c de 
with a, b, c, d, e ¡Ê R. Find a basis for this subspace. (Use the method to .nd a basis for the 
kernel of a matrix). 
Problem 7.8. If A is an n ¡Á n symmetric matrix and B is any n ¡Á n invertible matrix, 

prove that A is positive de.nite i. BTAB is positive de.nite. Problem 7.9. (1) Consider the matrix 
.
. 
A4 = 
... 

2 .10 0 
.12 .10 

0 .12 .1 
00 .12 

...

. 

Find three matrices of the form E2,1;¦Â1 ,E3,2;¦Â2 ,E4,3;¦Â3 , such that E4,3;¦Â3 E3,2;¦Â2 E2,1;¦Â1 A4 = U4 where U4 is an upper triangular matrix. Compute M = E4,3;¦Â3 E3,2;¦Â2 E2,1;¦Â1 
and check that 

.
. 
... 

2 .10 0 
03/2 .10 

004/3 .1 
00 05/4 

...

MA4 = U4 = 
. 

(2) Now consider the matrix 

.
. 
..... 

2  .1  0  0  0  
.1  2  .1  0  0  
0  .1  2  .1  0  

00 .12 .1 
000 .12 

..... 

A5 = 
. 

Find four matrices of the form E2,1;¦Â1 ,E3,2;¦Â2 ,E4,3;¦Â3 ,E5,4;¦Â4 , such that E5,4;¦Â4 E4,3;¦Â3 E3,2;¦Â2 E2,1;¦Â1 A5 = U5 where U5 is an upper triangular matrix. Compute M = E5,4;¦Â4 E4,3;¦Â3 E3,2;¦Â2 E2,1;¦Â1 
and check that 

.
. 
..... 

2 .10 0 0 
03/2 .10 0 

0  0  4/3  .1  0  
0  0  0  5/4  .1  
0  0  0  0  6/5  

..... 

MA5 = U5 = 
. 

(3) Write a Matlab program de.ning the function Ematrix(n, i, j, b) which is the n ¡Á n matrix that adds b times row j to row i. Also write some Matlab code that produces an n ¡Á n matrix An generalizing the matrices A4 and A5. 
Use your program to .gure out which .ve matrices Ei,j;¦Â reduce A6 to the upper triangular 
matrix

.
. 
....... 

2 .10 0 0 0 
03/2 .10 0 0 
004/3 .10 0 

0  0  0  5/4  .1  0  
0  0  0  0  6/5  .1  
0  0  0  0  0  7/6  

....... 

U6 = 
. 

7.17. PROBLEMS 
Also use your program to .gure out which six matrices Ei,j;¦Â reduce A7 to the upper trian-
gular matrix 

.
. 
U7 = 
......... 

2 .10 0 0 0 0 
03/2 .10 0 0 0 
004/3 .10 0 0 

00 05/4 .10 0 

0  0  0  0  6/5  .1  0  
0  0  0  0  0  7/6  .1  
0  0  0  0  0  0  8/7  

......... 

. 

(4) Find the lower triangular matrices L6 and L7 such that 
L6U6 = A6 

and L7U7 = A7. 
(5) It is natural to conjecture that there are n . 1 matrices of the form Ei,j;¦Â that reduce An to the upper triangular matrix 
.
. 
Un = 
.......... 

2 .10000 0 03/2 .1000 0 004/3 .100 0 00 05/4 .10 0 
0  0  0  0  6/5  ...  . . .  
. . .  . . .  . . .  . . .  ...  ...  .1  
0  0  0  0  ¡¤ ¡¤ ¡¤  0  (n + 1)/n  

.......... 

, 

namely, 
E2,1;1/2,E3,2;2/3,E4,3;3/4, ¡¤¡¤¡¤ ,En,n.1;(n.1)/n. It is also natural to conjecture that the lower triangular matrix Ln such that LnUn = An is given by 
Ln = E2,1;.1/2E3,2;.2/3E4,3;.3/4 ¡¤¡¤¡¤ En,n.1;.(n.1)/n, 
that is, 

.
. 
10000 0 	0 0 0 
.......... 

.1/21000 0 

0 .2/3100 0 

0 .
. 
..
000 .4/51 .. .... 
..
......
0 
00 .3/410 0

Ln = 
. 

. 

.

.... 

.......... 
0000 ¡¤¡¤¡¤ .(n . 1)/n 1 
264 CHAPTER 7. GAUSSIAN ELIMINATION, LU, CHOLESKY, ECHELON FORM Prove the above conjectures. 
(6) Prove that the last column of A.n 1 is 
.
. 
.... 

1/(n + 1) 
2/(n + 1) 

. 
. 
. 
n/(n + 1) 

.... 

. 

Problem 7.10. (1) Let A be any invertible 2 ¡Á 2 matrix 
ab 
A = . 
cd Prove that there is an invertible matrix S such that 10 
SA = ,
0 ad . bc 
where S is the product of at most four elementary matrices of the form Ei,j;¦Â. Conclude that every matrix A in SL(2) (the group of invertible 2 ¡Á 2 matrices A with det(A) = +1) is the product of at most four elementary matrices of the form Ei,j;¦Â. For any a =0, 1, give an explicit factorization as above for 
a 0 
A = .
.1
0 aWhat is this decomposition for a = .1? 
(2) Recall that a rotation matrix R (a member of the group SO(2)) is a matrix of the form 
cos ¦È . sin ¦È 
R = . 
sin ¦È cos ¦È Prove that if ¦È = k¦Ð (with k ¡Ê Z), any rotation matrix can be written as a product R = ULU, where U is upper triangular and L is lower triangular of the form 1 u 10 
U = ,L = . 
01 v 1 
Therefore, every plane rotation (except a .ip about the origin when ¦È = ¦Ð) can be written as the composition of three shear transformations! 
7.17. PROBLEMS 
Problem 7.11. (1) Recall that Ei,d is the diagonal matrix 
Ei,d = diag(1,..., 1, d, 1,..., 1), 
whose diagonal entries are all +1, except the (i, i)th entry which is equal to d. 
Given any n ¡Á n matrix A, for any pair (i, j) of distinct row indices (1 ¡Ü i, j ¡Ü n), prove that there exist two elementary matrices E1(i, j) and E2(i, j) of the form Ek,¿Û ;¦Â, such that 
Ej,.1E1(i, j)E2(i, j)E1(i, j)A = P (i, j)A, 
the matrix obtained from the matrix A by permuting row i and row j. Equivalently, we have 
E1(i, j)E2(i, j)E1(i, j)A = Ej,.1P (i, j)A, 
the matrix obtained from A by permuting row i and row j and multiplying row j by .1. 
Prove that for every i =2,...,n, there exist four elementary matrices E3(i, d),E4(i, d), E5(i, d),E6(i, d) of the form Ek,¿Û ;¦Â, such that 
E6(i, d)E5(i, d)E4(i, d)E3(i, d)En,d = Ei,d. 
What happens when d = .1, that is, what kind of simpli.cations occur? 
Prove that all permutation matrices can be written as products of elementary operations of the form Ek,¿Û ;¦Â and the operation En,.1. 
(2) Prove that for every invertible n ¡Á n matrix A, there is a matrix S such that 
In.1 0 
SA == En,d,
0 d 
with d = det(A), and where S is a product of elementary matrices of the form Ek,¿Û ;¦Â. 
In particular, every matrix in SL(n) (the group of invertible n ¡Á n matrices A with det(A) = +1) can be written as a product of elementary matrices of the form Ek,¿Û ;¦Â. Prove that at most n(n + 1) . 2 such transformations are needed. 
(3) Prove that every matrix in SL(n) can be written as a product of at most (n . 
1)(max{n, 3} + 1) elementary matrices of the form Ek,¿Û ;¦Â. 
Problem 7.12. A matrix A is called strictly column diagonally dominant i. 
n 
|ajj| > |aij|, for j =1,...,n i=1,iJ
=j 
Prove that if A is strictly column diagonally dominant, then Gaussian elimination with partial pivoting does not require pivoting, and A is invertible. 
.
..
. 
1000 1000 0 
0 1331 0121 
(2) What is the e.ect of the product (on the left) with 
E4,3;.1E3,2;.1E4,3;.1E2,1;.1E3,2;.1E4,3;.1 
...
... 

... 

... 

1100 

010 

E 

= 

. 

1210 

011 

on the matrix 

.
. 
1000 0
...
0 1331 
(3) Find the inverse of the matrix Pa3. 
(4) Consider the (n + 1) ¡Á (n + 1) Pascal matrix Pan whose ith row is given by the binomial coe.cients 
i . 1 ,
j . 1 with 1 ¡Ü i ¡Ü n +1, 1 ¡Ü j ¡Ü n + 1, and with the usual convention that 
0 i 
=1, =0 if j >i. 
0 j 
The matrix Pa3 is shown in Question (c) and Pa4 is shown below: 
... 

110 

Pa3 = 
. 

121 

.
. 
..... 

1  0  0  0  0  
1  1  0  0  0  
1  2  1  0  0  

13310 
14641 

..... 

Pa4 = 
. 

Find n elementary matrices Eik,jk;¦Âk such that 
10 

Ein,jn;¦Ân ¡¤¡¤¡¤ Ei1,j1;¦Â1 Pan = . 
0 Pan.1 Use the above to prove that the inverse of Pan is the lower triangular matrix whose ith row is given by the signed binomial coe.cients i . 1 
(.1)i+j.2 ,
j . 1 
7.17. PROBLEMS 
with 1 ¡Ü i ¡Ü n +1, 1 ¡Ü j ¡Ü n + 1. For example, 
.
. 
Pa.41 
= 

..... 

1 0 0 00 
.11 0 00 

1 .21 00 

.13 .310 
1 .46 .41 

..... 

. 

Hint. Given any n ¡Á n matrix A, multiplying A by the elementary matrix Ei,j;¦Â on the right yields the matrix AEi,j;¦Â in which ¦Â times the ith column is added to the jth column. 
Problem 7.14. (1) Implement the method for converting a rectangular matrix to reduced row echelon form in Matlab. 
(2) 
Use the above method to .nd the inverse of an invertible n ¡Á n matrix A by applying it to the the n ¡Á 2n matrix [AI] obtained by adding the n columns of the identity matrix to 

A. 


(3) Consider the matrix 
.
. 
A = 

...... 

12 3 4 ¡¤¡¤¡¤ n 23 4 5 ¡¤¡¤¡¤ n +1 34 5 6 ¡¤¡¤¡¤ n +2 
... .
.
.. . ..
. nn +1 n +2 n +3 ¡¤¡¤¡¤ 2n . 1 ... . 
...... 

. 

Using your program, .nd the row reduced echelon form of A for n =4,..., 20. 
Also run the Matlab rref function and compare results. 
Your program probably disagrees with rref even for small values of n. The problem is that some pivots are very small and the normalization step (to make the pivot 1) causes roundo. errors. Use a tolerance parameter to .x this problem. 
What can you conjecture about the rank of A? 
(4) Prove that the matrix A has the following row reduced form: 
.
. 
R = 

...... 

10 .1 .2 ¡¤¡¤¡¤ .(n . 2) 012 3 ¡¤¡¤¡¤ n . 1 000 0 ¡¤¡¤¡¤ 0 
... .
.
... ..
. 000 0 ¡¤¡¤¡¤ 0 ... . 
...... 

. 

Deduce from the above that A has rank 2. Hint. Some well chosen sequence of row operations. 
(5) Use your program to show that if you add any number greater than or equal to (2/25)n2 to every diagonal entry of A you get an invertible matrix! In fact, running the Matlab fuction chol should tell you that these matrices are SPD (symmetric, positive de.-nite). 
Problem 7.15. Let A be an n ¡Á n complex Hermitian positive de.nite matrix. Prove that the lower-triangular matrix B with positive diagonal entries such that A = BB. is given by the following formulae: For j =1,...,n, 
1/2
j.1 
bjj = ajj .|bjk|2 , k=1 
and for i = j +1,...,n (and j =1,...,n . 1) 
j.1 
bij = aij . bikbjk /bjj. k=1 
Problem 7.16. (Permutations and permutation matrices) A permutation can be viewed as an operation permuting the rows of a matrix. For example, the permutation 
1234 3421 
corresponds to the matrix 

.
. 
... 

0001 
0010 

1000 
0100 

...

P¦Ð = 
. 

Observe that the matrix P¦Ð has a single 1 on every row and every column, all other entries being zero, and that if we multiply any 4 ¡Á 4 matrix A by P¦Ð on the left, then the rows of A are permuted according to the permutation ¦Ð; that is, the ¦Ð(i)th row of P¦ÐA is the ith row of A. For example, 
.
..
..
. 
... 
0001 
a11 a12 a13 a14 a41 a42 a43 a44 
0
... 

... 

... 

... 

...

001 

a21 a22 a23 a24 a31 a32 a33 a34
P¦ÐA = 
= 

. 

0 
0100 

a41 a42 a43 a44 a21 a22 a23 a24 
Equivalently, the ith row of P¦ÐA is the ¦Ð.1(i)th row of A. In order for the matrix P¦Ð to move the ith row of A to the ¦Ð(i)th row, the ¦Ð(i)th row of P¦Ð must have a 1 in column i and zeros everywhere else; this means that the ith column of P¦Ð contains the basis vector e¦Ð(i), the vector that has a 1 in position ¦Ð(i) and zeros everywhere else. 
This is the general situation and it leads to the following de.nition. 
De.nition 7.8. Given any permutation ¦Ð :[n] ¡ú [n], the permutation matrix P¦Ð =(pij) representing ¦Ð is the matrix given by 
100 

a31 a32 a33 a34 a11 a12 a13 a14 
 

1 if i = ¦Ð(j) 
pij =
0 if i = ¦Ð(j); 
7.17. PROBLEMS 
equivalently, the jth column of P¦Ð is the basis vector e¦Ð(j).A permutation matrix P is any matrix of the form P¦Ð (where P is an n ¡Á n matrix, and ¦Ð :[n] ¡ú [n] is a permutation, for some n ¡Ý 1). 
Remark: There is a confusing point about the notation for permutation matrices. A per-mutation matrix P acts on a matrix A by multiplication on the left by permuting the rows of A. As we said before, this means that the ¦Ð(i)th row of P¦ÐA is the ith row of A, or equivalently that the ith row of P¦ÐA is the ¦Ð.1(i)th row of A. But then observe that the row index of the entries of the ith row of PA is ¦Ð.1(i), and not ¦Ð(i)! See the following example: 
.
..
..
. 
0001 
a11 a12 a13 a14 a41 a42 a43 a44 
... 

0010 

1000 

... 
... 

a21 a22 a23 a24 
a31 a32 a33 a34 
... 

= 

... 

a31 a32 a33 a34 
a11 a12 a13 a14 
...

, 

0100 
a41 a42 a43 a44 a21 a22 a23 a24 
where 

¦Ð.1(1) = 4 ¦Ð.1(2) = 3 ¦Ð.1(3) = 1 ¦Ð.1(4) = 2. 
Prove the following results 
(1) 
Given any two permutations ¦Ð1,¦Ð2 :[n] ¡ú [n], the permutation matrix P¦Ð2.¦Ð1 repre-senting the composition of ¦Ð1 and ¦Ð2 is equal to the product P¦Ð2 P¦Ð1 of the permutation matrices P¦Ð1 and P¦Ð2 representing ¦Ð1 and ¦Ð2; that is, 

P¦Ð2.¦Ð1 = P¦Ð2 P¦Ð1 . 

(2) 
The matrix P¦Ð.1 representing the inverse of the permutation ¦Ð1 is the inverse P¦Ð.1 1 of 


1 
the matrix P¦Ð1 representing the permutation ¦Ð1; that is, 
P¦Ð.1 = P¦Ð.1 1 . 

1 
Furthermore, 
P .1 )T 
¦Ð1 =(P¦Ð1 . 
(3) 
Prove that if P is the matrix associated with a transposition, then det(P )= .1. 

(4) 
Prove that if P is a permutation matrix, then det(P )= ¡À1. 

(5) 
Use permutation matrices to give 	another proof of the fact that the parity of the number of transpositions used to express a permutation ¦Ð depends only on ¦Ð. 




