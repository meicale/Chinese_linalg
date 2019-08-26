Chapter 2 
Vector Spaces, Bases, Linear Maps 
2.1 	Motivations: Linear Combinations, Linear Inde-pendence and Rank 
In linear optimization problems, we often encounter systems of linear equations. For example, consider the problem of solving the following system of three linear equations in the three variables x1,x2,x3 ∈ R: 
x1 +2x2 . x3 =1 
2x1 + x2 + x3 =2 
x1 . 2x2 . 2x3 =3. 
One way to approach this problem is introduce the “vectors” u, v, w, and b, given by 
.... 	.... 
12 .11 u = .... ...2.
2 v =1 w =1 b = 1 .2 .23 
and to write our linear system as 
x1u + x2v + x3w = b. 
In the above equation, we used implicitly the fact that a vector z can be multiplied by a scalar λ ∈ R, where 
... . 
z1 λz1 
... .
λz = λz2 = λz2 , z3 λz3 
and two vectors y and and z can be added, where 
..... . 
y1 z1 y1 + z1 y + z = .y2. + .z2. = .y2 + z2. . y3 z3 y3 + z3 
15 

Also, given a vector 
.. 
x1 
..
x = x2 , 
x3 
we de.ne the additive inverse .x of x (pronounced minus x) as 
.. 
.x1 
..
.x = 	.x2 . .x3 
Observe that .x =(.1)x, the scalar multiplication of x by .1. 
The set of all vectors with three components is denoted by R3×1 . The reason for using the notation R3×1 rather than the more conventional notation R3 is that the elements of R3×1 
are column vectors; they consist of three rows and a single column, which explains the superscript 3 × 1. On the other hand, R3 = R × R × R consists of all triples of the form (x1,x2,x3), with x1,x2,x3 ∈ R, and these are row vectors. However, there is an obvious bijection between R3×1 and R3 and they are usually identi.ed. For the sake of clarity, in this introduction, we will denote the set of column vectors with n components by Rn×1 . 
An expression such as 
x1u + x2v + x3w 
where u, v, w are vectors and the xis are scalars (in R) is called a linear combination. Using this notion, the problem of solving our linear system 
x1u + x2v + x3w = b. 
is equivalent to determining whether b can be expressed as a linear combination of u, v, w. 
Now if the vectors u, v, w are linearly independent, which means that there is no triple (x1,x2,x3) 
= (0, 0, 0) such that 
x1u + x2v + x3w =03, 
it can be shown that every vector in R3×1 can be written as a linear combination of u, v, w. Here, 03 is the zero vector 
.. 
0 
..
03 =0 . 0 
It is customary to abuse notation and to write 0 instead of 03. This rarely causes a problem because in most cases, whether 0 denotes the scalar zero or the zero vector can be inferred from the context. 
In fact, every vector z ∈ R3×1 can be written in a unique way as a linear combination 
z = x1u + x2v + x3w. 
2.1. MOTIVATIONS: LINEAR COMBINATIONS, LINEAR INDEPENDENCE, RANK17 

This is because if z = x1u + x2v + x3w = y1u + y2v + y3w, 
then by using our (linear!) operations on vectors, we get 
(y1 . x1)u +(y2 . x2)v +(y3 . x3)w =0, 
which implies that y1 . x1 = y2 . x2 = y3 . x3 =0, 
by linear independence. Thus, 
y1 = x1,y2 = x2,y3 = x3, 
which shows that z has a unique expression as a linear combination, as claimed. Then our equation x1u + x2v + x3w = b 
has a unique solution, and indeed, we can check that 
x1 =1.4 x2 = .0.4 x3 = .0.4 
is the solution. But then, how do we determine that some vectors are linearly independent? One answer is to compute a numerical quantity det(u, v, w), called the determinant of 
(u, v, w), and to check that it is nonzero. In our case, it turns out that 
   

12 .1 
21 1 

1 .2 .2

      

= 15,

   

det(u, v, w)=

which con.rms that u, v, w are linearly independent. 
Other methods, which are much better for systems with a large number of variables, consist of computing an LU-decomposition or a QR-decomposition, or an SVD of the matrix consisting of the three columns u, v, w, 
.
. 
12 .1 
21 1 

、
／ 
A =uvw=

.

.

. 

1 .2 .2 

If we form the vector of unknowns 

.
. 

x1 
.

.

x = 

x2 
, 

x3 
then our linear combination x1u + x2v + x3w can be written in matrix form as 
. ... 
12 .1 x1 
. ...
x1u + x2v + x3w =21 1 x2 , 1 .2 .2 x3 
so our linear system is expressed by 
. ... .. 
12 .1 x1 1 
. ... ..
21 1 x2 =2 , 1 .2 .2 x3 3 
or more concisely as Ax = b. 
Now what if the vectors u, v, w are linearly dependent? For example, if we consider the 
vectors  
.  .  .  .  .  .  
1  2  .1  
u = .2.  v = . 1 .  w = . 1 . ,  
1  .1  2  

we see that u . v = w, 
a nontrivial linear dependence. It can be veri.ed that u and v are still linearly independent. Now for our problem x1u + x2v + x3w = b 
it must be the case that b can be expressed as linear combination of u and v. However, it turns out that u, v, b are linearly independent (one way to see this is to compute the determinant det(u, v, b)= .6), so b cannot be expressed as a linear combination of u and v and thus, our system has no solution. 
If we change the vector b to 
.. 
3 
..
b =3 , 0 
then b = u + v, 
and so the system x1u + x2v + x3w = b 
has the solution x1 =1,x2 =1,x3 =0. 
Actually, since w = u . v, the above system is equivalent to 
(x1 + x3)u +(x2 . x3)v = b, 

2.1. MOTIVATIONS: LINEAR COMBINATIONS, LINEAR INDEPENDENCE, RANK19 

and because u and v are linearly independent, the unique solution in x1 + x3 and x2 . x3 is 
x1 + x3 =1 
x2 . x3 =1, 
which yields an in.nite number of solutions parameterized by x3, namely 
x1 =1 . x3 
x2 =1+ x3. 
In summary, a 3 × 3 linear system may have a unique solution, no solution, or an in.nite number of solutions, depending on the linear independence (and dependence) or the vectors u, v, w, b. This situation can be generalized to any n × n system, and even to any n × m system (n equations in m variables), as we will see later. 
The point of view where our linear system is expressed in matrix form as Ax = b stresses the fact that the map x → Ax is a linear transformation. This means that 
A(λx)= λ(Ax) 
for all x ∈ R3×1 and all λ ∈ R and that 
A(u + v)= Au + Av, 
for all u, v ∈ R3×1 . We can view the matrix A as a way of expressing a linear map from R3×1 to R3×1 and solving the system Ax = b amounts to determining whether b belongs to the 
image of this linear map. Given a 3 × 3 matrix  A = . . a11 a21  a12 a22  a13 a23 . . ,  
a31  a32  a33  

whose columns are three vectors denoted A1,A2,A3, and given any vector x =(x1,x2,x3), we de.ned the product Ax as the linear combination 
.. 
a11x1 + a12x2 + a13x3 Ax = x1A1 + x2A2 + x3A3 = .a21x1 + a22x2 + a23x3. . 
a31x1 + a32x2 + a33x3 
The common pattern is that the ith coordinate of Ax is given by a certain kind of product called an inner product, of a row vector, the ith row of A, times the column vector x:
.. 
／ 、. x1 .
ai1 ai2 ai3· x2 = ai1x1 + ai2x2 + ai3x3. 
x3 
More generally, given any two vectors x =(x1,...,xn) and y =(y1,...,yn) ∈ Rn, their inner product denoted x · y, or（x, y）, is the number 
.
. 

x · y =

／
x1 x2 
··· 

xn
、 
· 

y1 y2 n
= 
n....
... 
.... 
xiyi. 
i=1 
yn 
Inner products play a very important role. First, we quantity 
√ 
22 )1/2
lxl2 = x · x =(x1 + ··· + xn
is a generalization of the length of a vector, called the Euclidean norm, or �2-norm. Second, it can be shown that we have the inequality 
|x · y|≤lxllyl , 
so if x, y = 0, the ratio (x · y)/(lxllyl) can be viewed as the cosine of an angle, the angle between x and y. In particular, if x · y = 0 then the vectors x and y make the angle π/2, that is, they are orthogonal. The (square) matrices Q that preserve the inner product, in the sense that（Qx, Qy） =（x, y） for all x, y ∈ Rn, also play a very important role. They can be thought of as generalized rotations. 
Returning to matrices, if A is an m × n matrix consisting of n columns A1,...,An (in Rm), and B is a n × p matrix consisting of p columns B1,...,Bp (in Rn) we can form the p vectors (in Rm) 
AB1, . . . , ABp. 
These p vectors constitute the m × p matrix denoted AB, whose jth column is ABj. But we know that the ith coordinate of ABj is the inner product of the ith row of A by the jth 
column of B,

.
. 

n....
... 
.... 
b1j b2jn
=
／
、 
ai1 ai2 ··· ain· aikbkj. 
k=1 
bnj 
Thus we have de.ned a multiplication operation on matrices, namely if A =(aik) is a m × n matrix and if B =(bjk) if n × p matrix, then their product AB is the m × n matrix whose entry on the ith row and the jth column is given by the inner product of the ith row of A by the jth column of B, 
n
n(AB)ij = aikbkj. k=1 
Beware that unlike the multiplication of real (or complex) numbers, if A and B are two n×n matrices, in general, AB = BA. 
2.1. 	MOTIVATIONS: LINEAR COMBINATIONS, LINEAR INDEPENDENCE, RANK21 Suppose that A is an n × n matrix and that we are trying to solve the linear system Ax = b, with b ∈ Rn . Suppose we can .nd an n × n matrix B such that 
BAi = ei,i =1, . . . , n, with ei = (0,..., 0, 1, 0 ..., 0), where the only nonzero entry is 1 in the ith slot. If we form 
the n × n matrix

.
. 
In = 
........ 

100 ··· 00 
010 ··· 00 
001 ··· 00 

. . .  . . .  . . .  ...  . . .  . . .  
0  0  0  · · ·  1  0  
0  0  0  · · ·  0  1  

........ 

, 

called the identity matrix, whose ith column is ei, then the above is equivalent to 
BA = In. 
If Ax = b, then multiplying both sides on the left by B, we get 
B(Ax)= Bb. 
But is is easy to see that B(Ax)=(BA)x = Inx = x, so we must have 
x = Bb. 
We can verify that x = Bb is indeed a solution, because it can be shown that 
A(Bb)=(AB)b = Inb = b. 
What is not obvious is that BA = In implies AB = In, but this is indeed provable. The matrix B is usually denoted A.1 and called the inverse of A. It can be shown that it is the unique matrix such that 
AA.1 = A.1A = In. 
If a square matrix A has an inverse, then we say that it is invertible or nonsingular, otherwise we say that it is singular. We will show later that a square matrix is invertible i. its columns are linearly independent i. its determinant is nonzero. 
In summary, if A is a square invertible matrix, then the linear system Ax = b has the unique solution x = A.1b. In practice, this is not a good way to solve a linear system because computing A.1 is too expensive. A practical method for solving a linear system is Gaussian elimination, discussed in Chapter 7. Other practical methods for solving a linear system Ax = b make use of a factorization of A (QR decomposition, SVD decomposition), using orthogonal matrices de.ned next. 
T
Given an m × n matrix A =(akl), the n × m matrix AT =(aij) whose ith row is the ith column of A, which means that aT = aji for i =1,...,n and j =1,...,m, is called the 
ij 
transpose of A. An n × n matrix Q such that 
QQT = QTQ = In 
is called an orthogonal matrix. Equivalently, the inverse Q.1 of an orthogonal matrix Q is equal to its transpose QT. Orthogonal matrices play an important role. Geometrically, they correspond to linear transformation that preserve length. A major result of linear algebra states that every m × n matrix A can be written as 
A = V ΣUT, 
where V is an m × m orthogonal matrix, U is an n × n orthogonal matrix, and Σ is an m × n matrix whose only nonzero entries are nonnegative diagonal entries σ1 ≥ σ2 ≥ ··· ≥ σp, where p = min(m, n), called the singular values of A. The factorization A = V ΣUT is called a singular decomposition of A, or SVD. 
The SVD can be used to “solve” a linear system Ax = b where A is an m × n matrix, even when this system has no solution. This may happen when there are more equations that variables (m>n) , in which case the system is overdetermined. 
Of course, there is no miracle, an unsolvable system has no solution. But we can look for a good approximate solution, namely a vector x that minimizes some measure of the error Ax . b. Legendre and Gauss used lAx . bl22, which is the squared Euclidean norm of the error. This quantity is di.erentiable, and it turns out that there is a unique vector x+ of minimum Euclidean norm that minimizes lAx . bl22. Furthermore, x+ is given by the expression x+ = A+b, where A+ is the pseudo-inverse of A, and A+ can be computed from an SVD A = V ΣUT of A. Indeed, A+ = UΣ+V T, where Σ+ is the matrix obtained from Σ by replacing every positive singular value σi by its inverse σi .1, leaving all zero entries intact, and transposing. 
Instead of searching for the vector of least Euclidean norm minimizing lAx . bl22, we can add the penalty term K lxl22 (for some positive K> 0) to lAx . bl22 and minimize the quantity lAx . bl22 + K lxl22 . This approach is called ridge regression. It turns out that 
+
there is a unique minimizer xgiven by x+ =(ATA + KIn).1ATb, as shown in the second volume. 
Another approach is to replace the penalty term K lxl22 by K lxl1, where lxl1 = |x1| + ··· + |xn| (the �1-norm of x). The remarkable fact is that the minimizers x of lAx . bl2 +
2 
K lxl1 tend to be sparse, which means that many components of x are equal to zero. This approach known as lasso is popular in machine learning and will be discussed in the second volume. 
Another important application of the SVD is principal component analysis (or PCA), an important tool in data analysis. 

2.1. MOTIVATIONS: LINEAR COMBINATIONS, LINEAR INDEPENDENCE, RANK23 

Yet another fruitful way of interpreting the resolution of the system Ax = b is to view this problem as an intersection problem. Indeed, each of the equations x1 +2x2 . x3 =1 2x1 + x2 + x3 =2 x1 . 2x2 . 2x3 =3 de.nes a subset of R3 which is actually a plane. The .rst equation 
x1 +2x2 . x3 =1 de.nes the plane H1 passing through the three points (1, 0, 0), (0, 1/2, 0), (0, 0, .1), on the coordinate axes, the second equation 
2x1 + x2 + x3 =2 de.nes the plane H2 passing through the three points (1, 0, 0), (0, 2, 0), (0, 0, 2), on the coordinate axes, and the third equation x1 . 2x2 . 2x3 =3 de.nes the plane H3 passing through the three points (3, 0, 0), (0, .3/2, 0), (0, 0, .3/2), on the coordinate axes. See Figure 2.1. 

The intersection Hi∩Hj of any two distinct planes Hi and Hj is a line, and the intersection H1 ∩ H2 ∩ H3 of the three planes consists of the single point (1.4, .0.4, .0.4), as illustrated in Figure 2.2. 

The planes corresponding to the system 
x1 +2x2 . x3 =1 2x1 + x2 + x3 =2 x1 . x2 +2x3 =3, 
are illustrated in Figure 2.3. This system has no solution since there is no point simultane-

ously contained in all three planes; see Figure 2.4. 

2.1. MOTIVATIONS: LINEAR COMBINATIONS, LINEAR INDEPENDENCE, RANK25 


Finally, the planes corresponding to the system 
x1 +2x2 . x3 =3 2x1 + x2 + x3 =3 x1 . x2 +2x3 =0, 
are illustrated in Figure 2.5. 

This system has in.nitely many solutions, given parametrically by (1 . x3, 1+ x3,x3). Geometrically, this is a line common to all three planes; see Figure 2.6. 
Under the above interpretation, observe that we are focusing on the rows of the matrix A, rather than on its columns, as in the previous interpretations. 

theredlinecommontoallthreeplanes. Typicallythedatasetisrepresentedasan matrix A whereeachrowcorresponds ×mn ≥-dimensionaldatapointandtypically, Inmostapplications,thedataarenot toan nmn. {}independentsotherankof A isalotsmallerthanmin,andthethegoalof low-rank m,n《{}《matrixand C isa k matrix,with kmin(here,means“muchsmallerthan”): × nm,n. 
Another great example of a real-world problem where linear algebra proves to be very e.ective is the problem of data compression, that is, of representing a very large data set 
. 
using a much smaller amount of storage. 
. 
decomposition is to factor A as the product of two matrices B and C, where B is a m × k 
. 
......... 

A 

m × n 

......... 

= 

......... 

B 

m × k 

. ......... 
. 

.

C

. 

k × n 

Now it is generally too costly to .nd an exact factorization as above, so we look for a low-rank matrix A' which is a “good” approximation of A. In order to make this statement precise, we need to de.ne a mechanism to determine how close two matrices are. This can be done using matrix norms, a notion discussed in Chapter 8. The norm of a matrix A is a nonnegative real number lAl which behaves a lot like the absolute value |x| of a real number 
x. Then our goal is to .nd some low-rank matrix A' that minimizes the norm 
lA . A'l2 , 
over all matrices A' of rank at most k, for some given k《 min{m, n}. Some advantages of a low-rank approximation are: 

2.2. VECTOR SPACES 
1. 
Fewer elements are required to represent 	A; namely, k(m + n) instead of mn. Thus less storage and fewer operations are needed to reconstruct A. 

2. 
Often, the process for obtaining the decomposition exposes the underlying structure of the data. Thus, it may turn out that “most” of the signi.cant data are concentrated along some directions called principal directions. 


Low-rank decompositions of a set of data have a multitude of applications in engineering, including computer science (especially computer vision), statistics, and machine learning. As we will see later in Chapter 21, the singular value decomposition (SVD) provides a very satisfactory solution to the low-rank approximation problem. Still, in many cases, the data sets are so large that another ingredient is needed: randomization. However, as a .rst step, linear algebra often yields a good initial solution. 
We will now be more precise as to what kinds of operations are allowed on vectors. In the early 1900, the notion of a vector space emerged as a convenient and unifying framework for working with “linear” objects and we will discuss this notion in the next few sections. 


2.2 Vector Spaces 
A (real) vector space is a set E together with two operations, +: E ×E → E and ·: R×E → E, called addition and scalar multiplication, that satisfy some simple properties. First of all, E under addition has to be a commutative (or abelian) group, a notion that we review next. 
However, keep in mind that vector spaces are not just algebraic objects; they are also geometric objects. 
De.nition 2.1. A group is a set G equipped with a binary operation ·: G × G → G that associates an element a · b ∈ G to every pair of elements a, b ∈ G, and having the following properties: · is associative, has an identity element e ∈ G, and every element in G is invertible 
(w.r.t. ·). More explicitly, this means that the following equations hold for all a, b, c ∈ G: 
(G1) a · (b · c)=(a · b) · c. 	(associativity); 
(G2) a · e = e · a = a. 	(identity); 
(G3) For every a ∈ G, there is some a.1 ∈ G such that 
.1 .1
a · a= a· a = e. 	(inverse). 
A group G is abelian (or commutative) if 
a · b = b · a for all a, b ∈ G. 

A set M together with an operation ·: M × M → M and an element e satisfying only Conditions (G1) and (G2) is called a monoid. For example, the set N = {0, 1, . . . , n, . . .}of natural numbers is a (commutative) monoid under addition with identity element 0. However, it is not a group. 
Some examples of groups are given below. 
Example 2.1. 
1. 
The set 	Z = {..., .n, . . . , .1, 0, 1, . . . , n, . . .} of integers is an abelian group under addition, with identity element 0. However, Z. = Z .{0} is not a group under multiplication; it is a commutative monoid with identity element 1. 

2. 
The set Q of rational numbers (fractions p/q with p, q ∈ Z and q = 0) is an abelian group under addition, with identity element 0. The set Q. = Q.{0} is also an abelian group under multiplication, with identity element 1. 

3. 
Similarly, the sets R of real numbers and C of complex numbers are abelian groups under addition (with identity element 0), and R. = R .{0} and C. = C .{0} are abelian groups under multiplication (with identity element 1). 

4. 
The sets Rn and Cn of n-tuples of real or complex numbers are abelian groups under componentwise addition: 


(x1,...,xn)+(y1,...,yn)=(x1 + y1,...,xn + yn), 
with identity element (0,..., 0). 
5. 
Given any nonempty set S, the set of bijections f : S → S, also called permutations of S, is a group under function composition (i.e., the multiplication of f and g is the composition g . f), with identity element the identity function idS. This group is not abelian as soon as S has more than two elements. 

6. 
The set of n × n matrices with real (or complex) coe.cients is an abelian group under addition of matrices, with identity element the null matrix. It is denoted by Mn(R) (or Mn(C)). 

7. 
The set R[X] of all polynomials in one variable X with real coe.cients, 


P (X)= anXn + an.1Xn.1 + ··· + a1X + a0, 
(with ai ∈ R), is an abelian group under addition of polynomials. The identity element is the zero polynomial. 
8. 
The set of n ×n invertible matrices with real (or complex) coe.cients is a group under matrix multiplication, with identity element the identity matrix In. This group is called the general linear group and is usually denoted by GL(n, R) (or GL(n, C)). 

2.2. VECTOR SPACES 

9. 
The set of n×n invertible matrices with real (or complex) coe.cients and determinant +1 is a group under matrix multiplication, with identity element the identity matrix In. This group is called the special linear group and is usually denoted by SL(n, R) (or SL(n, C)). 

10. 
The set of n × n invertible matrices with real coe.cients such that RRT = RTR = In and of determinant +1 is a group (under matrix multiplication) called the special orthogonal group and is usually denoted by SO(n) (where RT is the transpose of the matrix R, i.e., the rows of RT are the columns of R). It corresponds to the rotations in Rn . 

11. 
Given an open interval (a, b), the set C(a, b) of continuous functions f :(a, b) → R is an abelian group under the operation f + g de.ned such that 


(f + g)(x)= f(x)+ g(x) 
for all x ∈ (a, b). 
It is customary to denote the operation of an abelian group G by +, in which case the inverse a.1 of an element a ∈ G is denoted by .a. 
The identity element of a group is unique. In fact, we can prove a more general fact: 
Proposition 2.1. If a binary operation ·: M × M → M is associative and if e ' ∈ M is a left identity and e '' ∈ M is a right identity, which means that 
e ' · a = a for all a ∈ M (G2l) 
and 
a · e '' = a for all a ∈ M, (G2r) then e ' = e '' . 
Proof. If we let a = e '' in equation (G2l), we get 
' '' '' 
e · e = e, 
and if we let a = e ' in equation (G2r), we get 
''' ' 
e · e = e, 
and thus 
' ''' '' 
e = e · e = e, as claimed. 
Proposition 2.1 implies that the identity element of a monoid is unique, and since every group is a monoid, the identity element of a group is unique. Furthermore, every element in a group has a unique inverse. This is a consequence of a slightly more general fact: 
Proposition 2.2. In a monoid M with identity element e, if some element a ∈ M has some left inverse a ' ∈ M and some right inverse a '' ∈ M, which means that 
a ' · a = e (G3l) 
and 
a · a '' = e, (G3r) 
then a ' = a '' . 
Proof. Using (G3l) and the fact that e is an identity element, we have 
' '' '' '' 
(a · a) · a = e · a = a. 
Similarly, Using (G3r) and the fact that e is an identity element, we have 
' '' 
a · (a · a '' )= a · e = a. 
However, since M is monoid, the operation · is associative, so 
'' ' '''' 
a = a · (a · a '' )=(a · a) · a = a, 
as claimed. 
Remark: Axioms (G2) and (G3) can be weakened a bit by requiring only (G2r) (the exis-tence of a right identity) and (G3r) (the existence of a right inverse for every element) (or (G2l) and (G3l)). It is a good exercise to prove that the group axioms (G2) and (G3) follow from (G2r) and (G3r). 
A vector space is an abelian group E with an additional operation ·: K × E → E called scalar multiplication that allows rescaling a vector in E by an element in K. The set K itself is an algebraic structure called a .eld. A .eld is a special kind of stucture called a ring. These notions are de.ned below. We begin with rings. 
De.nition 2.2. A ring is a set A equipped with two operations +: A × A → A (called addition) and .: A × A → A (called multiplication) having the following properties: 
(R1) A is an abelian group w.r.t. +; 
(R2) . is associative and has an identity element 1 ∈ A; 
(R3) . is distributive w.r.t. +. 
2.2. VECTOR SPACES 
The identity element for addition is denoted 0, and the additive inverse of a ∈ A is denoted by .a. More explicitly, the axioms of a ring are the following equations which hold for all a, b, c ∈ A: 
a +(b + c)=(a + b)+ c (associativity of +) (2.1) a + b = b + a (commutativity of +) (2.2) a +0=0+ a = a (zero) (2.3) a +(.a)=(.a)+ a = 0 (additive inverse) (2.4) a . (b . c)=(a . b) . c (associativity of .) (2.5) a . 1=1 . a = a (identity for .) (2.6) (a + b) . c =(a . c)+(b . c) (distributivity) (2.7) a . (b + c)=(a . b)+(a . c) (distributivity) (2.8) 
The ring A is commutative if 
a . b = b . a for all a, b ∈ A. 
From (2.7) and (2.8), we easily obtain 
a . 0=0 . a = 0 (2.9) a . (.b)=(.a) . b = .(a . b). (2.10) 
Note that (2.9) implies that if 1 = 0, then a = 0 for all a ∈ A, and thus, A = {0}. The ring A = {0} is called the trivial ring. A ring for which 1 = 0 is called nontrivial. The multiplication a . b of two elements a, b ∈ A is often denoted by ab. 
The abelian group Z is a commutative ring (with unit 1), and for any .eld K, the abelian group K[X] of polynomials is also a commutative ring (also with unit 1). The set Z/mZ of residues modulo m where m is a positive integer is a commutative ring. 
A .eld is a commutative ring K for which K .{0} is a group under multiplication. 
De.nition 2.3. A set K is a .eld if it is a ring and the following properties hold: 
(F1) 0=1; 
(F2) K. = K .{0} is a group w.r.t. . (i.e., every a = 0 has an inverse w.r.t. .); 
(F3) . is commutative. 
If . is not commutative but (F1) and (F2) hold, we say that we have a skew .eld (or noncommutative .eld). 
Note that we are assuming that the operation . of a .eld is commutative. This convention is not universally adopted, but since . will be commutative for most .elds we will encounter, we may as well include this condition in the de.nition. 
Example 2.2. 
1. 
The rings Q, R, and C are .elds. 

2. 
The set Z/pZ of residues modulo p where p is a prime number is .eld. 

3. 
The set of (formal) fractions f(X)/g(X) of polynomials f(X),g(X) ∈ R[X], where g(X) is not the zero polynomial, is a .eld. 


Vector spaces are de.ned as follows. 
De.nition 2.4. A real vector space isa set E (of vectors) together with two operations +: E × E → E (called vector addition)1 and ·: R × E → E (called scalar multiplication) satisfying the following conditions for all α, β ∈ R and all u, v ∈ E; 
(V0) E is an abelian group w.r.t. +, with identity element 0;2 
(V1) α · (u + v)=(α · u)+(α · v); 
(V2) (α + β) · u =(α · u)+(β · u); 
(V3) (α . β) · u = α · (β · u); 
(V4) 1 · u = u. 
In (V3), . denotes multiplication in R. 
Given α ∈ R and v ∈ E, the element α · v is also denoted by αv. The .eld R is often called the .eld of scalars. 
In De.nition 2.4, the .eld R may be replaced by the .eld of complex numbers C, in which case we have a complex vector space. It is even possible to replace R by the .eld of rational numbers Q or by any arbitrary .eld K (for example Z/pZ, where p is a prime number), in which case we have a K-vector space (in (V3), . denotes multiplication in the .eld K). In most cases, the .eld K will be the .eld R of reals, but all results in this chapter hold for vector spaces over an arbitrary .eld. 
From (V0), a vector space always contains the null vector 0, and thus is nonempty. From (V1), we get α · 0 = 0, and α · (.v)= .(α · v). From (V2), we get 0 · v = 0, and (.α) · v = .(α · v). 
Another important consequence of the axioms is the following fact: 
Proposition 2.3. For any u ∈ E and any λ ∈ R, if λ =0 and λ · u =0, then u =0. 
1The symbol + is overloaded, since it denotes both addition in the .eld R and addition of vectors in E. It is usually clear from the context which + is intended. 
2The symbol 0 is also overloaded, since it represents both the zero in R (a scalar) and the identity element of E (the zero vector). Confusion rarely arises, but one may prefer using 0 for the zero vector. 
2.2. VECTOR SPACES 
Proof. Indeed, since λ = 0, it has a multiplicative inverse λ.1, so from λ · u = 0, we get 
λ.1 
· (λ · u)= λ.1 · 0. 
However, we just observed that λ.1 · 0 = 0, and from (V3) and (V4), we have 
λ.1 
· (λ · u)=(λ.1λ) · u =1 · u = u, 
and we deduce that u = 0. 
Remark: One may wonder whether axiom (V4) is really needed. Could it be derived from the other axioms? The answer is no. For example, one can take E = Rn and de.ne ·: R× Rn → Rn by 
λ · (x1,...,xn) = (0,..., 0) 
for all (x1,...,xn) ∈ Rn and all λ ∈ R. Axioms (V0)–(V3) are all satis.ed, but (V4) fails. Less trivial examples can be given using the notion of a basis, which has not been de.ned yet. 
The .eld R itself can be viewed as a vector space over itself, addition of vectors being addition in the .eld, and multiplication by a scalar being multiplication in the .eld. 
Example 2.3. 
1. 
The .elds R and C are vector spaces over R. 

2. 
The groups Rn and Cn are vector spaces over R, with scalar multiplication given by 


λ(x1,...,xn)=(λx1, . . . , λxn), 
for any λ ∈ R and with (x1,...,xn) ∈ Rn or (x1,...,xn) ∈ Cn, and Cn is a vector space over C with scalar multiplication as above, but with λ ∈ C. 
3. The ring R[X]n of polynomials of degree at most n with real coe.cients is a vector space over R, and the ring C[X]n of polynomials of degree at most n with complex coe.cients is a vector space over C, with scalar multiplication λ·P (X) of a polynomial 
P (X)= amXm + am.1Xm.1 + ··· + a1X + a0 
(with ai ∈ R or ai ∈ C) by the scalar λ (in R or C), with m ≤ n, given by 
λ · P (X)= λamXm + λam.1Xm.1 + ··· + λa1X + λa0. 
4. 
The ring R[X] of all polynomials with real coe.cients is a vector space over R, and the ring C[X] of all polynomials with complex coe.cients is a vector space over C, with the same scalar multiplication as above. 

5. 
The ring of n × n matrices Mn(R) is a vector space over R. 

6. 
The ring of m × n matrices Mm,n(R) is a vector space over R. 

7. 
The ring C(a, b) of continuous functions f :(a, b) → R is a vector space over R, with the scalar multiplication λf of a function f :(a, b) → R by a scalar λ ∈ R given by 

(λf)(x)= λf(x), for all x ∈ (a, b). 

8. 
A very important example of vector space is the set of linear maps between two vector spaces to be de.ned in Section 2.7. Here is an example that will prepare us for the vector space of linear maps. Let X be any nonempty set and let E be a vector space. The set of all functions f : X → E can be made into a vector space as follows: Given any two functions f : X → E and g : X → E, let (f + g): X → E be de.ned such that 


(f + g)(x)= f(x)+ g(x) 
for all x ∈ X, and for every λ ∈ R, let λf : X → E be de.ned such that 
(λf)(x)= λf(x) 
for all x ∈ X. The axioms of a vector space are easily veri.ed. 
Let E be a vector space. We would like to de.ne the important notions of linear combi-nation and linear independence. 
Before de.ning these notions, we need to discuss a strategic choice which, depending how it is settled, may reduce or increase headaches in dealing with notions such as linear combinations and linear dependence (or independence). The issue has to do with using sets of vectors versus sequences of vectors. 
 
2.3 Indexed Families; the Sum Notation
i∈I ai 
Our experience tells us that it is preferable to use sequences of vectors; even better, indexed families of vectors. (We are not alone in having opted for sequences over sets, and we are in good company; for example, Artin [3], Axler [4], and Lang [41] use sequences. Nevertheless, some prominent authors such as Lax [44] use sets. We leave it to the reader to conduct a survey on this issue.) 
Given a set A, recall that a sequence is an ordered n-tuple (a1,...,an) ∈ An of elements from A, for some natural number n. The elements of a sequence need not be distinct and the order is important. For example, (a1,a2,a1) and (a2,a1,a1) are two distinct sequences in A3 . Their underlying set is {a1,a2}. 
What we just de.ned are .nite sequences, which can also be viewed as functions from {1, 2,...,n} to the set A; the ith element of the sequence (a1,...,an) is the image of i under the function. This viewpoint is fruitful, because it allows us to de.ne (countably) in.nite 
2.3. INDEXED FAMILIES; THE SUM NOTATION I∈I AI 35 
sequences as functions s: N → A. But then, why limit ourselves to ordered sets such as {1,...,n} or N as index sets? 
The main role of the index set is to tag each element uniquely, and the order of the tags is not crucial, although convenient. Thus, it is natural to de.ne the notion of indexed family. 
De.nition 2.5. Given a set A, an I-indexed family of elements of A, for short a family, is a function a: I → A where I is any set viewed as an index set. Since the function a is determined by its graph 
{(i, a(i)) | i ∈ I}, 
the family a can be viewed as the set of pairs a = {(i, a(i)) | i ∈ I}. For notational simplicity, we write ai instead of a(i), and denote the family a = {(i, a(i)) | i ∈ I} by (ai)i∈I . 
For example, if I = {r, g, b, y} and A = N, the set of pairs 
a = {(r, 2), (g, 3), (b, 2), (y, 11)} 
is an indexed family. The element 2 appears twice in the family with the two distinct tags r and b. 
When the indexed set I is totally ordered, a family (ai)i∈I is often called an I-sequence. Interestingly, sets can be viewed as special cases of families. Indeed, a set A can be viewed as the A-indexed family {(a, a) | a ∈ I} corresponding to the identity function. 
Remark: An indexed family should not be confused with a multiset. Given any set A,a multiset is a similar to a set, except that elements of A may occur more than once. For example, if A = {a, b, c, d}, then {a, a, a, b, c, c, d, d} is a multiset. Each element appears with a certain multiplicity, but the order of the elements does not matter. For example, a has multiplicity 3. Formally, a multiset is a function s : A → N, or equivalently a set of pairs {(a, i) | a ∈ A}. Thus, a multiset is an A-indexed family of elements from N, but not a N-indexed family, since distinct elements may have the same multiplicity (such as c an d in the example above). An indexed family is a generalization of a sequence, but a multiset is a generalization of a set. 
We also need to take care of an annoying technicality, which is to de.ne sums of the form i∈I ai, where I is any .nite index set and (ai)i∈I is a family of elements in some set A equiped with a binary operation +: A × A → A which is associative (Axiom (G1)) and commutative. This will come up when we de.ne linear combinations. 
The issue is that the binary operation + only tells us how to compute a1 + a2 for two elements of A, but it does not tell us what is the sum of three of more elements. For example, how should a1 + a2 + a3 be de.ned? 
What we have to do is to de.ne a1 +a2 +a3 by using a sequence of steps each involving two elements, and there are two possible ways to do this: a1 +(a2 + a3) and (a1 + a2)+ a3. If our operation + is not associative, these are di.erent values. If it associative, then a1+(a2+a3)= (a1 + a2)+ a3, but then there are still six possible permutations of the indices 1, 2, 3, and if 
+ is not commutative, these values are generally di.erent. If our operation is commutative, 
then all six permutations have the same value. Thus, if + is associative and commutative, it seems intuitively clear that a sum of the form i∈I ai does not depend on the order of the operations used to compute it. 
This is indeed the case, but a rigorous proof requires induction, and such a proof is surprisingly involved. Readers may accept without proof the fact that sums of the form i∈I ai are indeed well de.ned, and jump directly to De.nition 2.6. For those who want to 
see the gory details, here we go. 
First, we de.ne sums i∈I ai, where I is a .nite sequence of distinct natural numbers, say I =(i1,...,im). If I =(i1,...,im) with m ≥ 2, we denote the sequence (i2,...,im) by I .{i1}. We proceed by induction on the size m of I. Let 
n 
ai = ai1 , if m =1, i∈I .. nn 
ai = ai1 + ai , if m> 1. 
i∈Ii∈I.{i1} 
For example, if I = (1, 2, 3, 4), we have 
n 
ai = a1 +(a2 +(a3 + a4)). 
i∈I 
If the operation + is not associative, the grouping of the terms matters. For instance, in general 
a1 +(a2 +(a3 + a4)) = (a1 + a2)+(a3 + a4). 
However, if the operation + is associative, the sum i∈I ai should not depend on the grouping of the elements in I, as long as their order is preserved. For example, if I = (1, 2, 3, 4, 5), J1 = (1, 2), and J2 = (3, 4, 5), we expect that 
n  n  n  
ai =  aj  +  aj  .  
i∈I  j∈J1  j∈J2  

This indeed the case, as we have the following proposition. 
Proposition 2.4. Given any nonempty set A equipped with an associative binary operation +: A × A → A, for any nonempty .nite sequence I of distinct natural numbers and for any partition of I into p nonempty sequences Ik1 ,...,Ikp , for some nonempty sequence K = (k1,...,kp) of distinct natural numbers such that ki <kj implies that α<β for all α ∈ Iki and all β ∈ Ikj , for every sequence (ai)i∈I of elements in A, we have 
n nn 
aα = aα . 
α∈Ik∈Kα∈Ik 
2.3. INDEXED FAMILIES; THE SUM NOTATION I∈I AI 37 
Proof. We proceed by induction on the size n of I. 
If n = 1, then we must have p = 1 and Ik1 = I, so the proposition holds trivially. 
Next, assume n> 1. If p = 1, then Ik1 = I and the formula is trivial, so assume that p ≥ 2 and write J =(k2,...,kp). There are two cases. 
Case 1. The sequence Ik1 has a single element, say β, which is the .rst element of I. In this case, write C for the sequence obtained from I by deleting its .rst element β. By de.nition, 
nn 
aα = aβ + aα , 
α∈Iα∈C 
and 
nn nn 
aα = aβ + aα . 
k∈Kα∈Ik j∈Jα∈Ij 
Since |C| = n . 1, by the induction hypothesis, we have 
n nn 
aα = aα , 
α∈Cj∈Jα∈Ij 
which yields our identity. 
Case 2. The sequence Ik1 has at least two elements. In this case, let β be the .rst element of I (and thus of Ik1 ), let I ' be the sequence obtained from I by deleting its .rst element β, let I ' be the sequence obtained from Ik1 by deleting its .rst element β, and let I ' = Iki for
k1 ki 
i =2,...,p. Recall that J =(k2,...,kp) and K =(k1,...,kp). The sequence I ' has n . 1 elements, so by the induction hypothesis applied to I ' and the Ik' i , we get 
nnn nnn 
= aα = aα + aα . α∈Ik∈K 三三 α∈Ij
三 aαα∈Ik α∈Ik1 j∈J 
If we add the lefthand side to aβ, by de.nition we get 
n 
aα. α∈I 
If we add the righthand side to aβ, using associativity and the de.nition of an indexed sum, we get 
nn 
aβ + aα + aα
三
α∈Ij∈Jα∈Ij
k1 
nn 
= aβ + aα + aα 
α三∈Ij∈Jα∈Ijk1 
n nn nn 
= aα + aα = aα , 
α∈Ik1 j∈Jα∈Ij k∈Kα∈Ik 
as claimed. 

If I = (1,...,n), we also write n ai instead of ai. Since + is associative, Propo-
i=1 i∈I sition 2.4 shows that the sum in =1 ai is independent of the grouping of its elements, which justi.es the use the notation a1 + ··· + an (without any parentheses). If we also assume that our associative binary operation on A is commutative, then we can show that the sum i∈I ai does not depend on the ordering of the index set I. 
Proposition 2.5. Given any nonempty set A equipped with an associative and commutative binary operation +: A × A → A, for any two nonempty .nite sequences I and J of distinct natural numbers such that J is a permutation of I (in other words, the underlying sets of I and J are identical), for every sequence (ai)i∈I of elements in A, we have 
nn 
aα = aα. 
α∈Iα∈J 
Proof. We proceed by induction on the number p of elements in I. If p = 1, we have I = J and the proposition holds trivially. 
If p> 1, to simplify notation, assume that I = (1,...,p) and that J is a permutation (i1,...,ip) of I. First, assume that 2 ≤ i1 ≤ p . 1, let J ' be the sequence obtained from J by deleting i1, I ' be the sequence obtained from I by deleting i1, and let P = (1, 2,...,i1.1) and Q =(i1 +1,...,p . 1,p). Observe that the sequence I ' is the concatenation of the sequences P and Q. By the induction hypothesis applied to J ' and I ' , and then by Proposition 2.4 applied to I ' and its partition (P, Q), we have 
i1.1p
nnn n 
aα = aα = ai + ai . α∈J三 α∈I三 i=1 i=i1+1 
If we add the lefthand side to ai1 , by de.nition we get 
n 
aα. α∈J 
If we add the righthand side to ai1 , we get 
i1.1p
nn 
ai1 + ai + ai . 
i=1 i=i1+1 
Using associativity, we get 
i1.1pi1.1p
nn nn 
ai1 + ai + ai = ai1 + ai + ai , 
i=1 i=i1+1 i=1 i=i1+1 
then using associativity and commutativity several times (more rigorously, using induction 

2.3. INDEXED FAMILIES; THE SUM NOTATION I∈I AI 39 
on i1 . 1), we get 
i1.1pi1.1p
nnn n 
ai1 + ai + ai = ai + ai1 + ai 
i=1 i=i1+1 i=1 i=i1+1 p
n 
= ai, i=1 
as claimed. The cases where i1 = 1 or i1 = p are treated similarly, but in a simpler manner since either P = () or Q = () (where () denotes the empty sequence). 
Having done all this, we can now make sense of sums of the form i∈I ai, for any .nite indexed set I and any family a =(ai)i∈I of elements in A, where A is a set equipped with a binary operation + which is associative and commutative. 
Indeed, since I is .nite, it is in bijection with the set {1,...,n} for some n ∈ N, and any total ordering主 on I corresponds to a permutation I of {1,...,n} (where we identify a permutation with its image). For any total ordering主 on I, we de.ne i∈I, ai as 
nn 
ai = aj. 
i∈I, j∈I. 
Then for any other total ordering主' on I, we have 
nn 
ai = aj, i∈I, 三 j∈I.三 
and since I and I 三 are di.erent permutations of {1,...,n}, by Proposition 2.5, we have 
nn 
aj = aj. 
j∈I. j∈I.三 
Therefore, the sum i∈I, ai does not depend on the total ordering on I. We de.ne the sum i∈I ai as the common value i∈I, ai for all total orderings主 of I. Here are some examples with A = R: 
√ √√ 
1. If I = {1, 2, 3}, a = {(1, 2), (2, .3), (3, 2)}, then ai =2 . 3+ 2= .1+ 2. 
i∈I 
√ √√ 
2. If I = {2, 5, 7}, a = {(2, 2), (5, .3), (7, 2)}, then ai =2 . 3+ 2= .1+ 2. 

3. If I = {r, g, b}, a = {(r, 2), (g, .3), (b, 1)}, then ai =2 . 3+1 = 0. 


i∈I 
i∈I 
2.4 Linear Independence, Subspaces 
One of the most useful properties of vector spaces is that they possess bases. What this means is that in every vector space E, there is some set of vectors, {e1,...,en}, such that every vector v ∈ E can be written as a linear combination, 
v = λ1e1 + ··· + λnen, 
of the ei, for some scalars, λ1,...,λn ∈ R. Furthermore, the n-tuple, (λ1,...,λn), as above is unique. 
This description is .ne when E has a .nite basis, {e1,...,en}, but this is not always the case! For example, the vector space of real polynomials, R[X], does not have a .nite basis but instead it has an in.nite basis, namely 
1, X,X2, ...,Xn , ... 
Given a set A, recall that an I-indexed family (ai)i∈I of elements of A (for short, a family) is a function a: I → A, or equivalently a set of pairs {(i, ai) | i ∈ I}. We agree that when I = .,(ai)i∈I = .. A family (ai)i∈I is .nite if I is .nite. 
Remark: When considering a family (ai)i∈I , there is no reason to assume that I is ordered. The crucial point is that every element of the family is uniquely indexed by an element of 
I. Thus, unless speci.ed otherwise, we do not assume that the elements of an index set are ordered. 
Given two disjoint sets I and J, the union of two families (ui)i∈I and (vj)j∈J , denoted as (ui)i∈I ∪ (vj)j∈J , is the family (wk)k∈(I∪J) de.ned such that wk = uk if k ∈ I, and wk = vk if k ∈ J. Given a family (ui)i∈I and any element v, we denote by (ui)i∈I ∪k (v) the family (wi)i∈I∪{k} de.ned such that, wi = ui if i ∈ I, and wk = v, where k is any index such that k/∈ I. Given a family (ui)i∈I ,a subfamily of (ui)i∈I is a family (uj)j∈J where J is any subset of I. 
In this chapter, unless speci.ed otherwise, it is assumed that all families of scalars are .nite (i.e., their index set is .nite). 
De.nition 2.6. Let E be a vector space. A vector v ∈ E is a linear combination of a family (ui)i∈I of elements of E i. there is a family (λi)i∈I of scalars in R such that 
n 
v = λiui. i∈I 
When I = ., we stipulate that v = 0. (By Proposition 2.5, sums of the form i∈I λiui are well de.ned.) We say that a family (ui)i∈I is linearly independent i. for every family (λi)i∈I of scalars in R, 
n 
λiui = 0 implies that λi = 0 for all i ∈ I. i∈I 
2.4. LINEAR INDEPENDENCE, SUBSPACES 
Equivalently, a family (ui)i∈I is linearly dependent i. there is some family (λi)i∈I of scalars in R such that n λiui = 0 and λj = 0 for some j ∈ I. i∈I 
We agree that when I = ., the family . is linearly independent. 
Observe that de.ning linear combinations for families of vectors rather than for sets of vectors has the advantage that the vectors being combined need not be distinct. For example, for I = {1, 2, 3} and the families (u, v, u) and (λ1,λ2,λ1), the linear combination 
n 
λiui = λ1u + λ2v + λ1u 
i∈I 
makes sense. Using sets of vectors in the de.nition of a linear combination does not allow such linear combinations; this is too restrictive. 
Unravelling De.nition 2.6, a family (ui)i∈I is linearly dependent i. either I consists of a single element, say i, and ui = 0, or |I|≥ 2 and some uj in the family can be expressed as a linear combination of the other vectors in the family. Indeed, in the second case, there is some family (λi)i∈I of scalars in R such that 
n 
λiui = 0 and λj = 0 for some j ∈ I, i∈I 
and since |I|≥ 2, the set I .{j} is nonempty and we get 
n 
uj = .λ.1λiui.
j i∈(I.{j}) 
Observe that one of the reasons for de.ning linear dependence for families of vectors rather than for sets of vectors is that our de.nition allows multiple occurrences of a vector. This is important because a matrix may contain identical columns, and we would like to say that these columns are linearly dependent. The de.nition of linear dependence for sets does not allow us to do that. 
The above also shows that a family (ui)i∈I is linearly independent i. either I = ., or I consists of a single element i and ui = 0, or |I|≥ 2 and no vector uj in the family can be expressed as a linear combination of the other vectors in the family. 
When I is nonempty, if the family (ui)i∈I is linearly independent, note that ui = 0 for all i ∈ I. Otherwise, if ui = 0 for some i ∈ I, then we get a nontrivial linear dependence 
i∈I λiui = 0 by picking any nonzero λi and letting λk = 0 for all k ∈ I with k = i, since λi0=0. If |I|≥ 2, we must also have ui = uj for all i, j ∈ I with i = j, since otherwise we get a nontrivial linear dependence by picking λi = λ and λj = .λ for any nonzero λ, and letting λk = 0 for all k ∈ I with k = i, j. 
Thus, the de.nition of linear independence implies that a nontrivial linearly independent family is actually a set. This explains why certain authors choose to de.ne linear indepen-dence for sets of vectors. The problem with this approach is that linear dependence, which is the logical negation of linear independence, is then only de.ned for sets of vectors. How-ever, as we pointed out earlier, it is really desirable to de.ne linear dependence for families allowing multiple occurrences of the same vector. 
Example 2.4. 
1. 
Any two distinct scalars λ, μ =0 in R are linearly dependent. 

2. 
In R3, the vectors (1, 0, 0), (0, 1, 0), and (0, 0, 1) are linearly independent. See Figure 

3. 
In R4, the vectors (1, 1, 1, 1), (0, 1, 1, 1), (0, 0, 1, 1), and (0, 0, 0, 1) are linearly indepen-dent. 

4. 
In R2, the vectors u = (1, 1), v = (0, 1) and w = (2, 3) are linearly dependent, since 


2.7. 

w =2u + v. 
See Figure 2.8. 
When I is .nite, we often assume that it is the set I = {1, 2,...,n}. In this case, we denote the family (ui)i∈I as (u1,...,un). 
The notion of a subspace of a vector space is de.ned as follows. 
De.nition 2.7. Given a vector space E, a subset F of E is a linear subspace (or subspace) of E i. F is nonempty and λu + μv ∈ F for all u, v ∈ F , and all λ, μ ∈ R. 
It is easy to see that a subspace F of E is indeed a vector space, since the restriction of +: E × E → E to F × F is indeed a function +: F × F → F , and the restriction of ·: R× E → E to R× F is indeed a function ·: R× F → F . 
Since a subspace F is nonempty, if we pick any vector u ∈ F and if we let λ = μ = 0, then λu + μu =0u +0u = 0, so every subspace contains the vector 0. 
The following facts also hold. The proof is left as an exercise. 
2.4. LINEAR INDEPENDENCE, SUBSPACES 

Proposition 2.6. 
(1) 
The intersection of any family (even in.nite) of subspaces of a vector space 	E is a subspace. 

(2) 
Let F be any subspace of a vector space E. For any nonempty .nite index set I, 


if 	(ui)i∈I is any family of vectors ui ∈ F and (λi)i∈I is any family of scalars, then i∈I λiui ∈ F . 
The subspace {0} will be denoted by (0), or even 0 (with a mild abuse of notation). 
Example 2.5. 
1. In R2, the set of vectors u =(x, y) such that 
x + y =0 
is the subspace illustrated by Figure 2.9. 
2. In R3, the set of vectors u =(x, y, z) such that 
x + y + z =0 
is the subspace illustrated by Figure 2.10. 
3. For any n ≥ 0, the set of polynomials f(X) ∈ R[X] of degree at most n is a subspace of R[X]. 

Figure 2.9: The subspace x + y = 0 is the line through the origin with slope .1. It consists of all vectors of the form λ(.1, 1). 

Figure 2.10: The subspace x + y + z = 0 is the plane through the origin with normal (1, 1, 1). 
4. The set of upper triangular n × n matrices is a subspace of the space of n × n matrices. 
Proposition 2.7. Given any vector space E, if S is any nonempty subset of E, then the smallest subspace（S） (or Span(S)) of E containing S is the set of all (.nite) linear combi-nations of elements from S. 
Proof. We prove that the set Span(S) of all linear combinations of elements of S is a subspace of E, leaving as an exercise the veri.cation that every subspace containing S also contains Span(S). 
First, Span(S) is nonempty since it contains S (which is nonempty). If u = i∈I λiui 
2.4. LINEAR INDEPENDENCE, SUBSPACES 
and v = j∈J μjvj are any two linear combinations in Span(S), for any two scalars λ, μ ∈ R, 
nn 
λu + μv = λλiui + μμjvj 
i∈Ij∈J 

nn 
= λλiui + μμjvj 
i∈Ij∈J 
nn n 
= λλiui +(λλi + μμi)ui + μμjvj, 
i∈I.Ji∈I∩Jj∈J.I 
which is a linear combination with index set I ∪ J, and thus λu + μv ∈ Span(S), which proves that Span(S) is a subspace. 
One might wonder what happens if we add extra conditions to the coe.cients involved 
in forming linear combinations. Here are three natural restrictions which turn out to be 
important (as usual, we assume that our index sets are .nite): 
(1) Consider combinations i∈I λiui for which 
n 
λi =1. i∈I 
These are called a.ne combinations. One should realize that every linear combination i∈I λiui can be viewed as an a.ne combination. For example, if k is an index not in I, if we let J = I ∪{k}, uk = 0, and λk =1 . λi, then λjuj is an a.ne 
i∈Ij∈J 
combination and 
nn 
λiui = λjuj. 
i∈Ij∈J 
However, we get new spaces. For example, in R3, the set of all a.ne combinations of the three vectors e1 = (1, 0, 0),e2 = (0, 1, 0), and e3 = (0, 0, 1), is the plane passing through these three points. Since it does not contain 0 = (0, 0, 0), it is not a linear subspace. 
(2) Consider combinations i∈I λiui for which 
λi ≥ 0, for all i ∈ I. 
These are called positive (or conic) combinations. It turns out that positive combina-tions of families of vectors are cones. They show up naturally in convex optimization. 
(3) Consider combinations i∈I λiui for which we require (1) and (2), that is 
n 
λi =1, and λi ≥ 0 for all i ∈ I. 
i∈I 
These are called convex combinations. Given any .nite family of vectors, the set of all convex combinations of these vectors is a convex polyhedron. Convex polyhedra play a very important role in convex optimization. 
Remark: The notion of linear combination can also be de.ned for in.nite index sets I. To ensure that a sum i∈I λiui makes sense, we restrict our attention to families of .nite support. 
De.nition 2.8. Given any .eld K, a family of scalars (λi)i∈I has .nite support if λi =0 for all i ∈ I . J, for some .nite subset J of I. 
If (λi)i∈I is a family of scalars of .nite support, for any vector space E over K, for any (possibly in.nite) family (ui)i∈I of vectors ui ∈ E, we de.ne the linear combination i∈I λiui as the .nite linear combination λjuj, where J is any .nite subset of I such that λi =0 
j∈J 
for all i ∈ I . J. In general, results stated for .nite families also hold for families of .nite support. 
2.5 Bases of a Vector Space 
Given a vector space E, given a family (vi)i∈I , the subset V of E consisting of the null vector 0 and of all linear combinations of (vi)i∈I is easily seen to be a subspace of E. The family (vi)i∈I is an economical way of representing the entire subspace V , but such a family would be even nicer if it was not redundant. Subspaces having such an “e.cient” generating family (called a basis) play an important role and motivate the following de.nition. 
De.nition 2.9. Given a vector space E and a subspace V of E, a family (vi)i∈I of vectors vi ∈ V spans V or generates V i. for every v ∈ V , there is some family (λi)i∈I of scalars in R such that n 
v = λivi. i∈I 
We also say that the elements of (vi)i∈I are generators of V and that V is spanned by (vi)i∈I , or generated by (vi)i∈I . If a subspace V of E is generated by a .nite family (vi)i∈I , we say that V is .nitely generated. A family (ui)i∈I that spans V and is linearly independent is called a basis of V . 
Example 2.6. 
1. 
In R3, the vectors (1, 0, 0), (0, 1, 0), and (0, 0, 1), illustrated in Figure 2.9, form a basis. 

2. 
The vectors (1, 1, 1, 1), (1, 1, .1, .1), (1, .1, 0, 0), (0, 0, 1, .1) form a basis of R4 known as the Haar basis. This basis and its generalization to dimension 2n are crucial in wavelet theory. 

3. 
In the subspace of polynomials in R[X] of degree at most n, the polynomials 1, X, X2 , ...,Xn form a basis. 

4. 
The Bernstein polynomials (1 . X)n.kXk for k =0,...,n, also form a basis of 


n 
k 
that space. These polynomials play a major role in the theory of spline curves. 
2.5. BASES OF A VECTOR SPACE 
The .rst key result of linear algebra is that every vector space E has a basis. We begin with a crucial lemma which formalizes the mechanism for building a basis incrementally. 
Lemma 2.8. Given a linearly independent family (ui)i∈I of elements of a vector space E, if v ∈ E is not a linear combination of (ui)i∈I , then the family (ui)i∈I ∪k (v) obtained by adding v to the family (ui)i∈I is linearly independent (where k/∈ I). 
Proof. Assume that μv + i∈I λiui = 0, for any family (λi)i∈I of scalars in R. If μ = 0, then μ has an inverse (because R is a .eld), and thus we have v = . i∈I (μ.1λi)ui, showing that v is a linear combination of (ui)i∈I and contradicting the hypothesis. Thus, μ = 0. But then, we have i∈I λiui = 0, and since the family (ui)i∈I is linearly independent, we have λi =0 for all i ∈ I. 
The next theorem holds in general, but the proof is more sophisticated for vector spaces that do not have a .nite set of generators. Thus, in this chapter, we only prove the theorem for .nitely generated vector spaces. 
Theorem 2.9. Given any .nite family S =(ui)i∈I generating a vector space E and any linearly independent subfamily L =(uj)j∈J of S (where J . I), there is a basis B of E such that L . B . S. 
Proof. Consider the set of linearly independent families B such that L . B . S. Since this set is nonempty and .nite, it has some maximal element (that is, a subfamily B =(uh)h∈H of S with H . I of maximum cardinality), say B =(uh)h∈H . We claim that B generates 
E. Indeed, if B does not generate E, then there is some up ∈ S that is not a linear combination of vectors in B (since S generates E), with p/∈ H. Then by Lemma 2.8, the family B ' =(uh)h∈H∪{p} is linearly independent, and since L . B . B ' . S, this contradicts the maximality of B. Thus, B is a basis of E such that L . B . S. 
Remark: Theorem 2.9 also holds for vector spaces that are not .nitely generated. In this case, the problem is to guarantee the existence of a maximal linearly independent family B such that L . B . S. The existence of such a maximal family can be shown using Zorn’s lemma; see Lang [41] (Theorem 5.1). 
A situation where the full generality of Theorem 2.9 is needed is the case of the vector 
√ 
space R over the .eld of coe.cients Q. The numbers 1 and 2 are linearly independent 
√ 
over Q, so according to Theorem 2.9, the linearly independent family L = (1, 2) can be extended to a basis B of R. Since R is uncountable and Q is countable, such a basis must be uncountable! 
The notion of a basis can also be de.ned in terms of the notion of maximal linearly independent family and minimal generating family. 
De.nition 2.10. Let (vi)i∈I be a family of vectors in a vector space E. We say that (vi)i∈I a maximal linearly independent family of E if it is linearly independent, and if for any vector w ∈ E, the family (vi)i∈I ∪k {w} obtained by adding w to the family (vi)i∈I is linearly dependent. We say that (vi)i∈I a minimal generating family of E if it spans E, and if for any index p ∈ I, the family (vi)i∈I.{p} obtained by removing vp from the family (vi)i∈I does not span E. 
The following proposition giving useful properties characterizing a basis is an immediate consequence of Lemma 2.8. 
Proposition 2.10. Given a vector space E, for any family B =(vi)i∈I of vectors of E, the following properties are equivalent: 
(1) 
B is a basis of E. 

(2) 
B is a maximal linearly independent family of E. 

(3) 
B is a minimal generating family of E. 


Proof. We will .rst prove the equivalence of (1) and (2). Assume (1). Since B is a basis, it is a linearly independent family. We claim that B is a maximal linearly independent family. If B is not a maximal linearly independent family, then there is some vector w ∈ E such that the family B ' obtained by adding w to B is linearly independent. However, since B is a basis of E, the vector w can be expressed as a linear combination of vectors in B, contradicting the fact that B ' is linearly independent. 
Conversely, assume (2). We claim that B spans E. If B does not span E, then there is some vector w ∈ E which is not a linear combination of vectors in B. By Lemma 2.8, the family B ' obtained by adding w to B is linearly independent. Since B is a proper subfamily of B ' , this contradicts the assumption that B is a maximal linearly independent family. Therefore, B must span E, and since B is also linearly independent, it is a basis of E. 
Now we will prove the equivalence of (1) and (3). Again, assume (1). Since B is a basis, it is a generating family of E. We claim that B is a minimal generating family. If B is not a minimal generating family, then there is a proper subfamily B ' of B that spans E. Then, every w ∈ B . B ' can be expressed as a linear combination of vectors from B ' , contradicting the fact that B is linearly independent. 
Conversely, assume (3). We claim that B is linearly independent. If B is not linearly independent, then some vector w ∈ B can be expressed as a linear combination of vectors in B ' = B .{w}. Since B generates E, the family B ' also generates E, but B ' is a proper subfamily of B, contradicting the minimality of B. Since B spans E and is linearly independent, it is a basis of E. 
The second key result of linear algebra is that for any two bases (ui)i∈I and (vj)j∈J of a vector space E, the index sets I and J have the same cardinality. In particular, if E has a .nite basis of n elements, every basis of E has n elements, and the integer n is called the dimension of the vector space E. 
To prove the second key result, we can use the following replacement lemma due to Steinitz. This result shows the relationship between .nite linearly independent families and .nite families of generators of a vector space. We begin with a version of the lemma which is 
2.5. BASES OF A VECTOR SPACE 
a bit informal, but easier to understand than the precise and more formal formulation given in Proposition 2.12. The technical di.culty has to do with the fact that some of the indices need to be renamed. 
Proposition 2.11. (Replacement lemma, version 1) Given a vector space E, let (u1,...,um) be any .nite linearly independent family in E, and let (v1,...,vn) be any .nite family such that every ui is a linear combination of (v1,...,vn). Then we must have m ≤ n, and there is a replacement of m of the vectors vj by (u1,...,um), such that after renaming some of the indices of the vjs, the families (u1,...,um,vm+1,...,vn) and (v1,...,vn) generate the same subspace of E. 
Proof. We proceed by induction on m. When m = 0, the family (u1,...,um) is empty, and the proposition holds trivially. For the induction step, we have a linearly independent family (u1,...,um,um+1). Consider the linearly independent family (u1,...,um). By the induction hypothesis, m ≤ n, and there is a replacement of m of the vectors vj by (u1,...,um), such that after renaming some of the indices of the vs, the families (u1,...,um,vm+1,...,vn) and (v1,...,vn) generate the same subspace of E. The vector um+1 can also be expressed as a lin-ear combination of (v1,...,vn), and since (u1,...,um,vm+1,...,vn) and (v1,...,vn) generate the same subspace, um+1 can be expressed as a linear combination of (u1,...,um,vm+1,..., vn), say 
n
n
mnum+1 = λiui + λjvj. i=1 j=m+1 
We claim that λj = 0 for some j with m +1 ≤ j ≤ n, which implies that m +1 ≤ n. Otherwise, we would have 
n
mum+1 = λiui, 
i=1 a nontrivial linear dependence of the ui, which is impossible since (u1,...,um+1) are linearly independent. 
Therefore, m +1 ≤ n, and after renaming indices if necessary, we may assume that λm+1 = 0, so we get 
n
n
mn
= . (λ.1 m+1um+1 . (λ.1 i=1 j=m+2 
vm+1 m+1λi)ui . λ.1 m+1λj)vj. 
Observe that the families (u1,...,um,vm+1,...,vn) and (u1,...,um+1,vm+2,...,vn) generate the same subspace, since um+1 is a linear combination of (u1,...,um,vm+1,...,vn) and vm+1 is a linear combination of (u1,...,um+1,vm+2,...,vn). Since (u1,...,um,vm+1,...,vn) and (v1,...,vn) generate the same subspace, we conclude that (u1,...,um+1,vm+2,...,vn) and and (v1,...,vn) generate the same subspace, which concludes the induction hypothesis. 
Here is an example illustrating the replacement lemma. Consider sequences (u1,u2,u3) and (v1,v2,v3,v4,v5), where (u1,u2,u3) is a linearly independent family and with the uis expressed in terms of the vjs as follows: 
u1 = v4 + v5 
u2 = v3 + v4 . v5 
u3 = v1 + v2 + v3. 
From the .rst equation we get 
v4 = u1 . v5, 
and by substituting in the second equation we have 
u2 = v3 + v4 . v5 = v3 + u1 . v5 . v5 = u1 + v3 . 2v5. 
From the above equation we get 
v3 = .u1 + u2 +2v5, 
and so u3 = v1 + v2 + v3 = v1 + v2 . u1 + u2 +2v5. 
Finally, we get 
v1 = u1 . u2 + u3 . v2 . 2v5 Therefore we have 
v1 = u1 . u2 + u3 . v2 . 2v5 v3 = .u1 + u2 +2v5 v4 = u1 . v5, 
which shows that (u1,u2,u3,v2,v5) spans the same subspace as (v1,v2,v3,v4,v5). The vectors (v1,v3,v4) have been replaced by (u1,u2,u3), and the vectors left over are (v2,v5). We can rename them (v4,v5). 
For the sake of completeness, here is a more formal statement of the replacement lemma (and its proof). 
Proposition 2.12. (Replacement lemma, version 2) Given a vector space E, let (ui)i∈I be any .nite linearly independent family in E, where |I| = m, and let (vj)j∈J be any .nite family such that every ui is a linear combination of (vj)j∈J , where |J| = n. Then there exists a set L and an injection ρ: L → J (a relabeling function) such that L ∩ I = ., |L| = n . m, and the families (ui)i∈I ∪ (vρ(l))l∈L and (vj)j∈J generate the same subspace of E. In particular, m ≤ n. 
Proof. We proceed by induction on |I| = m. When m = 0, the family (ui)i∈I is empty, and the proposition holds trivially with L = J (ρ is the identity). Assume |I| = m + 1. Consider the linearly independent family (ui)i∈(I.{p}), where p is any member of I. By the induction hypothesis, there exists a set L and an injection ρ: L → J such that L ∩ (I .{p})= ., 
2.5. BASES OF A VECTOR SPACE 
|L| = n . m, and the families (ui)i∈(I.{p}) ∪ (vρ(l))l∈L and (vj)j∈J generate the same subspace of E. If p ∈ L, we can replace L by (L .{p}) ∪{p ' } where p ' does not belong to I ∪ L, and replace ρ by the injection ρ ' which agrees with ρ on L .{p} and such that ρ ' (p ' )= ρ(p). Thus, we can always assume that L ∩ I = .. Since up is a linear combination of (vj)j∈J and the families (ui)i∈(I.{p}) ∪ (vρ(l))l∈L and (vj)j∈J generate the same subspace of E, up is a linear combination of (ui)i∈(I.{p}) ∪ (vρ(l))l∈L. Let 
nn 
up = λiui + λlvρ(l). (1) 
i∈(I.{p}) l∈L 
If λl = 0 for all l ∈ L, we have 
n 
λiui . up =0, 
i∈(I.{p}) 
contradicting the fact that (ui)i∈I is linearly independent. Thus, λl = 0 for some l ∈ L, say l = q. Since λq = 0, we have 
nn 
(.λ.1λi)ui + λ.1 (.λ.1 
vρ(q) = qq up + λl)vρ(l). (2)
q i∈(I.{p}) l∈(L.{q}) 
We claim that the families (ui)i∈(I.{p}) ∪ (vρ(l))l∈L and (ui)i∈I ∪ (vρ(l))l∈(L.{q}) generate the same subset of E. Indeed, the second family is obtained from the .rst by replacing vρ(q) by up, and vice-versa, and up is a linear combination of (ui)i∈(I.{p}) ∪ (vρ(l))l∈L, by (1), and vρ(q) is a linear combination of (ui)i∈I ∪(vρ(l))l∈(L.{q}), by (2). Thus, the families (ui)i∈I ∪(vρ(l))l∈(L.{q}) and (vj)j∈J generate the same subspace of E, and the proposition holds for L .{q} and the restriction of the injection ρ: L → J to L .{q}, since L ∩ I = . and |L| = n . m imply that (L .{q}) ∩ I = . and |L .{q}| = n . (m + 1). 
The idea is that m of the vectors vj can be replaced by the linearly independent uis in such a way that the same subspace is still generated. The purpose of the function ρ: L → J is to pick n . m elements j1,...,jn.m of J and to relabel them l1,...,ln.m in such a way that these new indices do not clash with the indices in I; this way, the vectors vj1 ,...,vjn.m who “survive” (i.e. are not replaced) are relabeled vl1 ,...,vln.m , and the other m vectors vj with j ∈ J .{j1,...,jn.m} are replaced by the ui. The index set of this new family is I ∪ L. 
Actually, one can prove that Proposition 2.12 implies Theorem 2.9 when the vector space is .nitely generated. Putting Theorem 2.9 and Proposition 2.12 together, we obtain the following fundamental theorem. 
Theorem 2.13. Let E be a .nitely generated vector space. Any family (ui)i∈I generating E contains a subfamily (uj)j∈J which is a basis of E. Any linearly independent family (ui)i∈I can be extended to a family (uj)j∈J which is a basis of E (with I . J). Furthermore, for every two bases (ui)i∈I and (vj)j∈J of E, we have |I| = |J| = n for some .xed integer n ≥ 0. 
Proof. The .rst part follows immediately by applying Theorem 2.9 with L = . and S = (ui)i∈I . For the second part, consider the family S ' =(ui)i∈I ∪ (vh)h∈H , where (vh)h∈H is any .nitely generated family generating E, and with I ∩ H = .. Then apply Theorem 2.9 to L =(ui)i∈I and to S ' . For the last statement, assume that (ui)i∈I and (vj)j∈J are bases of 
E. Since (ui)i∈I is linearly independent and (vj)j∈J spans E, Proposition 2.12 implies that |I|≤|J|. A symmetric argument yields |J|≤|I|. 
Remark: Theorem 2.13 also holds for vector spaces that are not .nitely generated. 
De.nition 2.11. When a vector space E is not .nitely generated, we say that E is of in.nite dimension. The dimension of a .nitely generated vector space E is the common dimension n of all of its bases and is denoted by dim(E). 
Clearly, if the .eld R itself is viewed as a vector space, then every family (a) where a ∈ R and a = 0 is a basis. Thus dim(R) = 1. Note that dim({0}) = 0. 
De.nition 2.12. If E is a vector space of dimension n ≥ 1, for any subspace U of E, if dim(U) = 1, then U is called a line; if dim(U) = 2, then U is called a plane; if dim(U)= n.1, then U is called a hyperplane. If dim(U)= k, then U is sometimes called a k-plane. 
Let (ui)i∈I be a basis of a vector space E. For any vector v ∈ E, since the family (ui)i∈I generates E, there is a family (λi)i∈I of scalars in R, such that 
n 
v = λiui. i∈I 
A very important fact is that the family (λi)i∈I is unique. 
Proposition 2.14. Given a vector space E, let (ui)i∈I be a family of vectors in E. Let v ∈ E, and assume that v = i∈I λiui. Then the family (λi)i∈I of scalars such that v = i∈I λiui is unique i. (ui)i∈I is linearly independent. 
Proof. First, assume that (ui)i∈I is linearly independent. If (μi)i∈I is another family of scalars in R such that v = i∈I μiui, then we have 
n 
(λi . μi)ui =0, i∈I 
and since (ui)i∈I is linearly independent, we must have λi.μi = 0 for all i ∈ I, that is, λi = μi for all i ∈ I. The converse is shown by contradiction. If (ui)i∈I was linearly dependent, there would be a family (μi)i∈I of scalars not all null such that 
n 
μiui =0 
i∈I 
2.6. MATRICES 
and μj = 0 for some j ∈ I. But then, 
n
n
n
n 
v = λiui +0= λiui + μiui =(λi + μi)ui, 
i∈Ii∈Ii∈Ii∈I 
with λj = λj +μj since μj = 0, contradicting the assumption that (λi)i∈I is the unique family such that v = i∈I λiui. 
De.nition 2.13. If (ui)i∈I is a basis of a vector space E, for any vector v ∈ E, if (xi)i∈I is the unique family of scalars in R such that 
n 

v = xiui, i∈I 
each xi is called the component (or coordinate) of index i of v with respect to the basis (ui)i∈I . 
2.6 Matrices 
In Section 2.1 we introduced informally the notion of a matrix. In this section we de.ne matrices precisely, and also introduce some operations on matrices. It turns out that matri-ces form a vector space equipped with a multiplication operation which is associative, but noncommutative. We will explain in Section 3.1 how matrices can be used to represent linear maps, de.ned in the next section. 
De.nition 2.14. If K = R or K = C, an m × n-matrix over K is a family (aij)1≤i≤m, 1≤j≤n of scalars in K, represented by an array 
.
. 
.... 

a11 a12 ... a1 n a21 a22 ... a2 n 
.. .
.
. ...
.
.. . 
am 1 am 2 ... amn 
.... 

. 

In the special case where m = 1, we have a row vector, represented by (a11 ··· a1 n) and in the special case where n = 1, we have a column vector, represented by 
.
. 

.. 

a11 
. 
. 
. 
am 1 
..

. 

In these last two cases, we usually omit the constant index 1 (.rst index in case of a row, second index in case of a column). The set of all m × n-matrices is denoted by Mm,n(K) or Mm,n. An n × n-matrix is called a square matrix of dimension n. The set of all square matrices of dimension n is denoted by Mn(K), or Mn. 
Remark: As de.ned, a matrix A =(aij)1≤i≤m, 1≤j≤n is a family, that is, a function from {1, 2,...,m}×{1, 2,...,n} to K. As such, there is no reason to assume an ordering on the indices. Thus, the matrix A can be represented in many di.erent ways as an array, by adopting di.erent orders for the rows or the columns. However, it is customary (and usually convenient) to assume the natural ordering on the sets {1, 2,...,m} and {1, 2,...,n}, and to represent A as an array according to this ordering of the rows and columns. 
We de.ne some operations on matrices as follows. 
De.nition 2.15. Given two m × n matrices A =(aij) and B =(bij), we de.ne their sum A + B as the matrix C =(cij) such that cij = aij + bij; that is, 
.
..
. 
a11 a12 ... a1 n b11 b12 ... b1 n 
.... 

.... 

+ 

.... 

.... 

b21 b22 ... b2 na21 a22 ... a2 n 
.. .
.
.. .
.
.. 

.. 

.. 

..

. 

.

.. 

. 

.. 

. 

am 1 am 2 ... amn bm 1 bm 2 ... bmn 
.
. 
.... 

a11 + b11 a12 + b12 ... a1 n + b1 n a21 + b21 a22 + b22 ... a2 n + b2 n 
.. .
.
. ...
.
.. . 
am 1 + bm 1 am 2 + bm 2 ... amn + bmn 
.... 

= 

. 

For any matrix A =(aij), we let .A be the matrix (.aij). Given a scalar λ ∈ K, we de.ne the matrix λA as the matrix C =(cij) such that cij = λaij; that is 
.
..
. 
a11 a12 ... a1 n λa11 λa12 ... λa1 n 
λ 

.... 

.... 

= 

.... 

.... 

. 

λa21 λa22 ... λa2 na21 a22 ... a2 n 
.. .
.
.. .
.
.. 

.. 

.. 

..

. 

.

.. 

. 

.. 

. 

am 1 am 2 ... amn λam 1 λam 2 ... λamn 
Given an m × n matrices A =(aik) and an n × p matrices B =(bkj), we de.ne their product AB as the m × p matrix C =(cij) such that 
n
n

cij = aikbkj, k=1 
for 1 ≤ i ≤ m, and 1 ≤ j ≤ p. In the product AB = C shown below 
.
..
..
. 
b11 b12 ... b1 pa11 a12 ... a1 n c11 c12 ... c1 p
.... 

.... 
.... 

.... 

= 

.... 

....

, 

b21 b22 ... b2 pa21 a22 ... a2 n c21 c22 ... c2 p 
.. .
.
.. .
.
.. .
.
.. 

.. 

.. 

.. 

.. 

..

. 

. 

.

.. 

. 

.. 

. 

.. 

. 

am 1 am 2 ... amn bn 1 bn 2 ... bnp cm 1 cm 2 ... cmp 
note that the entry of index i and j of the matrix AB obtained by multiplying the matrices A and B can be identi.ed with the product of the row matrix corresponding to the i-th row 
2.6. MATRICES 
of A with the column matrix corresponding to the j-column of B: 
.
. 

b1 j 
nk=1
bnj 
De.nition 2.16. The square matrix In of dimension n containing 1 on the diagonal and 0 everywhere else is called the identity matrix. It is denoted by 
n
.. 

..

.

(ai 1 ··· ain) aikbkj.
. 

= 

. 

.
. 
.... 

10 ... 0 

01 ... 0 


.. .
.
.. ..
.
.. . 
00 ... 1 
.... 

In = 
De.nition 2.17. Given an m × n matrix A =(aij), its transpose AT =(aT ), is the 
ji
n × m-matrix such that aT = aij, for all i,1 ≤ i ≤ m, and all j,1 ≤ j ≤ n.
ji 
The transpose of a matrix A is sometimes denoted by At, or even by tA. Note that the transpose AT of a matrix A has the property that the j-th row of AT is the j-th column of 
A. In other words, transposition exchanges the rows and the columns of a matrix. Here is an example. If A is the 5 × 6 matrix 
.
. 
1 23456 
7 12345 

8 71234 

9 87123 
10 9 8 7 1 2 

..... 

A = 

..... 

, 

then AT is the 6 × 5 matrix 
.
. 
....... 

1  7  8  9  10  
2  1  7  8  9  
3  2  1  7  8  

4  3  2  1  7  
5  4  3  2  1  
6  5  4  3  2  

....... 

AT 
= 
. 

The following observation will be useful later on when we discuss the SVD. Given any m × n matrix A and any n × p matrix B, if we denote the columns of A by A1,...,An and the rows of B by B1,...,Bn, then we have 
AB = A1B1 + ··· + AnBn. 
For every square matrix A of dimension n, it is immediately veri.ed that AIn = InA = A. 
De.nition 2.18. For any square matrix A of dimension n, if a matrix B such that AB = BA = In exists, then it is unique, and it is called the inverse of A. The matrix B is also denoted by A.1 . An invertible matrix is also called a nonsingular matrix, and a matrix that is not invertible is called a singular matrix. 
Using Proposition 2.19 and the fact that matrices represent linear maps, it can be shown that if a square matrix A has a left inverse, that is a matrix B such that BA = I, or a right inverse, that is a matrix C such that AC = I, then A is actually invertible; so B = A.1 and C = A.1 . These facts also follow from Proposition 5.11. 
It is immediately veri.ed that the set Mm,n(K) of m × n matrices is a vector space under addition of matrices and multiplication of a matrix by a scalar. 
De.nition 2.19. The m × n-matrices Eij =(ehk), are de.ned such that eij = 1, and ehk = 0, if h = i or k = j; in other words, the (i, j)-entry is equal to 1 and all other entries are 0. 
Here are the Eij matrices for m = 2 and n = 3: 
100 010 001 
E11 = ,E12 = ,E13 = 
000 000 000 
000 000 000 
E21 = ,E22 = ,E23 = . 
100 010 001 
It is clear that every matrix A =(aij) ∈ Mm,n(K) can be written in a unique way as 
nmn
n 
A = aijEij. 
i=1 j=1 
Thus, the family (Eij)1≤i≤m,1≤j≤n is a basis of the vector space Mm,n(K), which has dimension mn. 
Remark: De.nition 2.14 and De.nition 2.15 also make perfect sense when K is a (com-mutative) ring rather than a .eld. In this more general setting, the framework of vector spaces is too narrow, but we can consider structures over a commutative ring A satisfying all the axioms of De.nition 2.4. Such structures are called modules. The theory of modules is (much) more complicated than that of vector spaces. For example, modules do not always have a basis, and other properties holding for vector spaces usually fail for modules. When a module has a basis, it is called a free module. For example, when A is a commutative ring, the structure An is a module such that the vectors ei, with (ei)i =1 and (ei)j = 0 for j = i, form a basis of An . Many properties of vector spaces still hold for An . Thus, An is a free module. As another example, when A is a commutative ring, Mm,n(A) is a free module with basis (Ei,j)1≤i≤m,1≤j≤n. Polynomials over a commutative ring also form a free module of in.nite dimension. 
The properties listed in Proposition 2.15 are easily veri.ed, although some of the com-putations are a bit tedious. A more conceptual proof is given in Proposition 3.1. 
2.7. LINEAR MAPS 
Proposition 2.15. (1) Given any matrices A ∈ Mm,n(K), B ∈ Mn,p(K), and C ∈ Mp,q(K), we have 
(AB)C = A(BC); 
that is, matrix multiplication is associative. 
(2) Given any matrices A, B ∈ Mm,n(K), and C, D ∈ Mn,p(K), for all λ ∈ K, we have 
(A + B)C = AC + BC A(C + D)= AC + AD (λA)C = λ(AC) A(λC)= λ(AC), 
so that matrix multiplication ·:Mm,n(K) × Mn,p(K) → Mm,p(K) is bilinear. 
The properties of Proposition 2.15 together with the fact that AIn = InA = A for all square n×n matrices show that Mn(K) is a ring with unit In (in fact, an associative algebra). This is a noncommutative ring with zero divisors, as shown by the following example. 
Example 2.7. For example, letting A, B be the 2 × 2-matrices 
10 00 
A = ,B = ,
00 10 
then 1000 00 
AB == ,
0010 00 
and 0010 00 
BA == . 
1000 10 
Thus AB = BA, and AB = 0, even though both A, B = 0. 
2.7 Linear Maps 
Now that we understand vector spaces and how to generate them, we would like to be able to transform one vector space E into another vector space F . A function between two vector spaces that preserves the vector space structure is called a homomorphism of vector spaces, or linear map. Linear maps formalize the concept of linearity of a function. 
Keep in mind that linear maps, which are transformations of space, are usually far more important than the spaces themselves. 
In the rest of this section, we assume that all vector spaces are real vector spaces, but all results hold for vector spaces over an arbitrary .eld. 
De.nition 2.20. Given two vector spaces E and F ,a linear map between E and F is a function f : E → F satisfying the following two conditions: 
f(x + y)= f(x)+ f(y) for all x, y ∈ E; f(λx)= λf(x) for all λ ∈ R,x ∈ E. 
Setting x = y = 0 in the .rst identity, we get f(0) = 0. The basic property of linear maps is that they transform linear combinations into linear combinations. Given any .nite family (ui)i∈I of vectors in E, given any family (λi)i∈I of scalars in R, we have 
nn 
f( λiui)= λif(ui). i∈Ii∈I 
The above identity is shown by induction on |I| using the properties of De.nition 2.20. 
Example 2.8. 
1. The map f : R2 → R2 de.ned such that 
x ' = x . y 

y ' = x + y is a linear map. The reader should check that it is the composition of a rotation by 
√ 
π/4 with a magni.cation of ratio 2. 
2. 
For any vector space E, the identity map id: E → E given by 

id(u)= u for all u ∈ E is a linear map. When we want to be more precise, we write idE instead of id. 

3. 
The map D : R[X] → R[X] de.ned such that 

D(f(X)) = f ' (X), 
where f ' (X) is the derivative of the polynomial f(X), is a linear map. 


4. 
The map Φ: C([a, b]) → R given by 


 b 
Φ(f)=f(t)dt, 
a 
where C([a, b]) is the set of continuous functions de.ned on the interval [a, b], is a linear map. 
2.7. LINEAR MAPS 
5. The function（.,.）: C([a, b]) ×C([a, b]) → R given by
b 
（f, g） = f(t)g(t)dt, 
a 
is linear in each of the variable f, g. It also satis.es the properties（f, g） =（g, f） and（f, f） = 0 i. f = 0. It is an example of an inner product. 
De.nition 2.21. Given a linear map f : E → F , we de.ne its image (or range) Im f = f(E), as the set 
Im f = {y ∈ F | (.x ∈ E)(y = f(x))}, and its Kernel (or nullspace) Ker f = f.1(0), as the set 
Ker f = {x ∈ E | f(x)=0}. 
The derivative map D : R[X] → R[X] from Example 2.8(3) has kernel the constant polynomials, so Ker D = R. If we consider the second derivative D . D : R[X] → R[X], then the kernel of D .D consists of all polynomials of degree ≤ 1. The image of D : R[X] → R[X] is actually R[X] itself, because every polynomial P (X)= a0Xn + ··· + an.1X + an of degree n is the derivative of the polynomial Q(X) of degree n + 1 given by 
Xn+1 X2 Q(X)= a0 + ··· + an.1 + anX. 
n +1 2 On the other hand, if we consider the restriction of D to the vector space R[X]n of polyno-mials of degree ≤ n, then the kernel of D is still R, but the image of D is the R[X]n.1, the vector space of polynomials of degree ≤ n . 1. 
Proposition 2.16. Given a linear map f : E → F , the set Im f is a subspace of F and the set Ker f is a subspace of E. The linear map f : E → F is injective i. Ker f = (0) (where 
(0) is the trivial subspace {0}). Proof. Given any x, y ∈ Im f, there are some u, v ∈ E such that x = f(u) and y = f(v), and for all λ, μ ∈ R, we have f(λu + μv)= λf(u)+ μf(v)= λx + μy, 
and thus, λx + μy ∈ Im f, showing that Im f is a subspace of F . Given any x, y ∈ Ker f, we have f(x)=0 and f(y) = 0, and thus, 
f(λx + μy)= λf(x)+ μf(y)=0, 
that is, λx + μy ∈ Ker f, showing that Ker f is a subspace of E. 
First, assume that Ker f = (0). We need to prove that f(x)= f(y) implies that x = y. However, if f(x)= f(y), then f(x) . f(y) = 0, and by linearity of f we get f(x . y) = 0. Because Ker f = (0), we must have x . y = 0, that is x = y, so f is injective. Conversely, assume that f is injective. If x ∈ Ker f, that is f(x) = 0, since f(0) = 0 wehave f(x)= f(0), and by injectivity, x = 0, which proves that Ker f = (0). Therefore, f is injective i. Ker f = (0). 
Since by Proposition 2.16, the image Im f of a linear map f is a subspace of F , we can de.ne the rank rk(f) of f as the dimension of Im f. 
De.nition 2.22. Given a linear map f : E → F , the rank rk(f) of f is the dimension of the image Im f of f. 
A fundamental property of bases in a vector space is that they allow the de.nition of linear maps as unique homomorphic extensions, as shown in the following proposition. 
Proposition 2.17. Given any two vector spaces E and F , given any basis (ui)i∈I of E, given any other family of vectors (vi)i∈I in F , there is a unique linear map f : E → F such that f(ui)= vi for all i ∈ I. Furthermore, f is injective i. (vi)i∈I is linearly independent, and f is surjective i. (vi)i∈I generates F . 
Proof. If such a linear map f : E → F exists, since (ui)i∈I is a basis of E, every vector x ∈ E can written uniquely as a linear combination 
n 
x = xiui, i∈I 
and by linearity, we must have 
nn 
f(x)= xif(ui)= xivi. i∈Ii∈I 
De.ne the function f : E → F , by letting 
n 
f(x)= xivi i∈I 
for every x = i∈I xiui. It is easy to verify that f is indeed linear, it is unique by the previous reasoning, and obviously, f(ui)= vi. Now assume that f is injective. Let (λi)i∈I be any family of scalars, and assume that 
n 
λivi =0. i∈I 
Since vi = f(ui) for every i ∈ I, we have 
nn n 
f( λiui)= λif(ui)= λivi =0. i∈Ii∈Ii∈I 
Since f is injective i. Ker f = (0), we have 
n 
λiui =0, i∈I 
2.7. LINEAR MAPS 
and since (ui)i∈I is a basis, we have λi = 0 for all i ∈ I, which shows that (vi)i∈I is linearly independent. Conversely, assume that (vi)i∈I is linearly independent. Since (ui)i∈I is a basis of E, every vector x ∈ E is a linear combination x = i∈I λiui of (ui)i∈I . If 
n 
f(x)= f( λiui)=0, i∈I 
then 
nn n 
λivi = λif(ui)= f( λiui)=0, i∈Ii∈Ii∈I 
and λi = 0 for all i ∈ I because (vi)i∈I is linearly independent, which means that x = 0. Therefore, Ker f = (0), which implies that f is injective. The part where f is surjective is left as a simple exercise. 
Figure 2.11 provides an illustration of Proposition 2.17 when E = R3 and V = R2 

By the second part of Proposition 2.17, an injective linear map f : E → F sends a basis (ui)i∈I to a linearly independent family (f(ui))i∈I of F , which is also a basis when f is bijective. Also, when E and F have the same .nite dimension n,(ui)i∈I is a basis of E, and 
f : E → F is injective, then (f(ui))i∈I is a basis of F (by Proposition 2.10). The following simple proposition is also useful. 
Proposition 2.18. Given any two vector spaces E and F , with F nontrivial, given any family (ui)i∈I of vectors in E, the following properties hold: 
(1) 
The family (ui)i∈I generates E i. for every family of vectors (vi)i∈I in F , there is at most one linear map f : E → F such that f(ui)= vi for all i ∈ I. 

(2) 
The family (ui)i∈I is linearly independent i. for every family of vectors (vi)i∈I in F , there is some linear map f : E → F such that f(ui)= vi for all i ∈ I. 


Proof. (1) If there is any linear map f : E → F such that f(ui)= vi for all i ∈ I, since (ui)i∈I generates E, every vector x ∈ E can be written as some linear combination 
n 
x = xiui, i∈I 
and by linearity, we must have 
nn 
f(x)= xif(ui)= xivi. i∈Ii∈I 
This shows that f is unique if it exists. Conversely, assume that (ui)i∈I does not generate E. Since F is nontrivial, there is some some vector y ∈ F such that y = 0. Since (ui)i∈I does not generate E, there is some vector w ∈ E that is not in the subspace generated by (ui)i∈I . By Theorem 2.13, there is a linearly independent subfamily (ui)i∈I0 of (ui)i∈I generating the same subspace. Since by hypothesis, w ∈ E is not in the subspace generated by (ui)i∈I0 , by Lemma 2.8 and by Theorem 2.13 again, there is a basis (ej)j∈I0∪J of E, such that ei = ui for all i ∈ I0, and w = ej0 for some j0 ∈ J. Letting (vi)i∈I be the family in F such that vi = 0 for all i ∈ I, de.ning f : E → F to be the constant linear map with value 0, we have a linear map such that f(ui)=0 for all i ∈ I. By Proposition 2.17, there is a unique linear map g : E → F such that g(w)= y, and g(ej)=0 for all j ∈ (I0 ∪ J) .{j0}. By de.nition of the basis (ej)j∈I0∪J of E, we have g(ui) = 0 for all i ∈ I, and since f = g, this contradicts the fact that there is at most one such map. See Figure 2.12. 
(2) If the family (ui)i∈I is linearly independent, then by Theorem 2.13, (ui)i∈I can be extended to a basis of E, and the conclusion follows by Proposition 2.17. Conversely, assume that (ui)i∈I is linearly dependent. Then there is some family (λi)i∈I of scalars (not all zero) such that 
n 
λiui =0. i∈I 
By the assumption, for any nonzero vector y ∈ F , for every i ∈ I, there is some linear map fi : E → F , such that fi(ui)= y, and fi(uj) = 0, for j ∈ I .{i}. Then we would get 
nn 
0= fi( λiui)= λifi(ui)= λiy, i∈Ii∈I 
and since y = 0, this implies λi = 0 for every i ∈ I. Thus, (ui)i∈I is linearly independent. 
Given vector spaces E, F , and G, and linear maps f : E → F and g : F → G, it is easily veri.ed that the composition g . f : E → G of f and g is a linear map. 
2.7. LINEAR MAPS 

De.nition 2.23. A linear map f : E → F is an isomorphism i. there is a linear map 
g : F → E, such that g . f = idE and f . g = idF . (.) 
The map g in De.nition 2.23 is unique. This is because if g and h both satisfy g.f = idE, f . g = idF , h . f = idE, and f . h = idF , then 
g = g . idF = g . (f . h)=(g . f) . h = idE . h = h. 
The map g satisfying (.) above is called the inverse of f and it is also denoted by f.1 . 
Observe that Proposition 2.17 shows that if F = Rn, then we get an isomorphism between any vector space E of dimension |J| = n and Rn . Proposition 2.17 also implies that if E and F are two vector spaces, (ui)i∈I is a basis of E, and f : E → F is a linear map which is an isomorphism, then the family (f(ui))i∈I is a basis of F . 
One can verify that if f : E → F is a bijective linear map, then its inverse f.1 : F → E, as a function, is also a linear map, and thus f is an isomorphism. 
Another useful corollary of Proposition 2.17 is this: 
Proposition 2.19. Let E be a vector space of .nite dimension n ≥ 1 and let f : E → E be any linear map. The following properties hold: 
(1) 
If f has a left inverse g, that is, if g is a linear map such that g . f = id, then f is an isomorphism and f.1 = g. 

(2) 
If f has a right inverse h, that is, if h is a linear map such that f . h = id, then f is an isomorphism and f.1 = h. 


Proof. (1) The equation g . f = id implies that f is injective; this is a standard result about functions (if f(x)= f(y), then g(f(x)) = g(f(y)), which implies that x = y since g . f = id). Let (u1,...,un) be any basis of E. By Proposition 2.17, since f is injective, (f(u1),...,f(un)) is linearly independent, and since E has dimension n, it is a basis of E (if (f(u1),...,f(un)) doesn’t span E, then it can be extended to a basis of dimension strictly greater than n, contradicting Theorem 2.13). Then f is bijective, and by a previous observation its inverse is a linear map. We also have 
g = g . id = g . (f . f.1)=(g . f) . f.1 = id . f.1 = f.1 . 
(2) The equation f . h = id implies that f is surjective; this is a standard result about functions (for any y ∈ E, we have f(h(y)) = y). Let(u1,...,un) be any basis of E. By Proposition 2.17, since f is surjective, (f(u1),...,f(un)) spans E, and since E has dimension n, it is a basis of E (if (f(u1),...,f(un)) is not linearly independent, then because it spans E, it contains a basis of dimension strictly smaller than n, contradicting Theorem 2.13). Then f is bijective, and by a previous observation its inverse is a linear map. We also have 
h = id . h =(f.1 . f) . h = f.1 . (f . h)= f.1 . id = f.1 . 
This completes the proof. 
De.nition 2.24. The set of all linear maps between two vector spaces E and F is denoted by Hom(E, F ) or by L(E; F ) (the notation L(E; F ) is usually reserved to the set of continuous linear maps, where E and F are normed vector spaces). When we wish to be more precise and specify the .eld K over which the vector spaces E and F are de.ned we write HomK (E, F ). 
The set Hom(E, F ) is a vector space under the operations de.ned in Example 2.3, namely 
(f + g)(x)= f(x)+ g(x) 
for all x ∈ E, and (λf)(x)= λf(x) 
for all x ∈ E. The point worth checking carefully is that λf is indeed a linear map, which uses the commutativity of . in the .eld K (typically, K = R or K = C). Indeed, we have 
(λf)(μx)= λf(μx)= λμf(x)= μλf(x)= μ(λf)(x). 
When E and F have .nite dimensions, the vector space Hom(E, F ) also has .nite di-mension, as we shall see shortly. 
2.8. LINEAR FORMS AND THE DUAL SPACE 
De.nition 2.25. When E = F , a linear map f : E → E is also called an endomorphism. The space Hom(E, E) is also denoted by End(E). 
It is also important to note that composition confers to Hom(E, E) a ring structure. Indeed, composition is an operation .: Hom(E, E) × Hom(E, E) → Hom(E, E), which is associative and has an identity idE, and the distributivity properties hold: 
(g1 + g2) . f = g1 . f + g2 . f; 
g . (f1 + f2)= g . f1 + g . f2. 
The ring Hom(E, E) is an example of a noncommutative ring. 
It is easily seen that the set of bijective linear maps f : E → E is a group under compo-sition. 
De.nition 2.26. Bijective linear maps f : E → E are also called automorphisms. The group of automorphisms of E is called the general linear group (of E), and it is denoted by GL(E), or by Aut(E), or when E = Rn, by GL(n, R), or even by GL(n). 
2.8 Linear Forms and the Dual Space 
We already observed that the .eld K itself (K = R or K = C) is a vector space (over itself). The vector space Hom(E, K) of linear maps from E to the .eld K, the linear forms, plays a particular role. In this section, we only de.ne linear forms and show that every .nite-dimensional vector space has a dual basis. A more advanced presentation of dual spaces and duality is given in Chapter 10. 
De.nition 2.27. Given a vector space E, the vector space Hom(E, K) of linear maps from E to the .eld K is called the dual space (or dual) of E. The space Hom(E, K) is also denoted by E., and the linear maps in E. are called the linear forms, or covectors. The dual space 
E.. 
of the space E. is called the bidual of E. 
As a matter of notation, linear forms f : E → K will also be denoted by starred symbol, such as u . , x . , etc. 
If E is a vector space of .nite dimension n and (u1,...,un) is a basis of E, for any linear form f. ∈ E., for every x = x1u1 + ··· + xnun ∈ E, by linearity we have 
f . (x)= f . (u1)x1 + ··· + f . (un)xn 
= λ1x1 + ··· + λnxn, 
with λi = f.(ui) ∈ K for every i,1 ≤ i ≤ n. Thus, with respect to the basis (u1,...,un), the linear form f. is represented by the row vector 
(λ1 ··· λn), 
we have 

.
. 

x1
、.. 
..

. 

.

λ1 ··· λn
f . (x)=
／
. 
xn 
, 

a linear combination of the coordinates of x, and we can view the linear form f. as a linear equation. If we decide to use a column vector of coe.cients 
. 

c1 
. 
c = .. 
cn instead of a row vector, then the linear form f. 
.. 
. .. 

is de.ned by 

f . (x)= c T x. 
The above notation is often used in machine learning. 
Example 2.9. Given any di.erentiable function f : Rn → R, by de.nition, for any x ∈ Rn , the total derivative dfx of f at x is the linear form dfx : Rn → R de.ned so that for all u =(u1,...,un) ∈ Rn , 
.
. 

n
.... 
u1 n
.f .f .f 
.
dfx(u)= (x) ··· (x).= (x) ui.
.
.x1 .xn .xi
i=1
un 
Example 2.10. Let C([0, 1]) be the vector space of continuous functions f : [0, 1] → R. The map I : C([0, 1]) → R given by 
1 
I(f)= f(x)dx for any f ∈C([0, 1]) 
0 
is a linear form (integration). 
Example 2.11. Consider the vector space Mn(R) of real n×n matrices. Let tr: Mn(R) → R be the function given by 
tr(A)= a11 + a22 + ··· + ann, 
called the trace of A. It is a linear form. Let s:Mn(R) → R be the function given by 
n
n
s(A)= aij, 
i,j=1 
where A =(aij). It is immediately veri.ed that s is a linear form. 
2.8. LINEAR FORMS AND THE DUAL SPACE 
Given a vector space E and any basis (ui)i∈I for E, we can associate to each ui a linear form ui . ∈ E., and the ui . have some remarkable properties. 
De.nition 2.28. Given a vector space E and any basis (ui)i∈I for E, by Proposition 2.17, for every i ∈ I, there is a unique linear form ui . such that 
 
1 if i = jui . (uj)=0 if i = j, 
for every j ∈ I. The linear form ui . is called the coordinate form of index i w.r.t. the basis (ui)i∈I . 
Remark: Given an index set I, authors often de.ne the so called “Kronecker symbol” δij such that  1 if i = jδij =
0 if i = j, 
for all i, j ∈ I. Then, ui .(uj)= δij. 
The reason for the terminology coordinate form is as follows: If E has .nite dimension and if (u1,...,un) is a basis of E, for any vector 
v = λ1u1 + ··· + λnun, 
we have 
u . i (v)= u . i (λ1u1 + ··· + λnun) 
.. . 
(u1)+ ··· + λiu (ui)+ ··· + λnu )
= λ1uiii (un
= λi, 
since ui .(uj)= δij. Therefore, ui . is the linear function that returns the ith coordinate of a vector expressed over the basis (u1,...,un). 
The following theorem shows that in .nite-dimension, every basis (u1,...,un) of a vector space E yields a basis (u1. ,...,u . n) of the dual space E., called a dual basis. 
Theorem 2.20. (Existence of dual bases) Let E be a vector space of dimension n. The following properties hold: For every basis (u1,...,un) of E, the family of coordinate forms (u1. ,...,u n. ) is a basis of E. (called the dual basis of (u1,...,un)). 
. ∈ E.
Proof. (a) If v is any linear form, consider the linear form 
f .. (u1)u .. (un. 
= v 1 + ··· + v )un. 
Observe that because ui .(uj)= δij, 
. (u1)u .. (un. 
f . (ui)=(v 1 + ··· + v )un)(ui) 
.. . 
= v . (u1)u1(ui)+ ··· + v . (ui)ui (ui)+ ··· + v . (un)un(ui) 
= v . (ui), 
and so f. and v . agree on the basis (u1,...,un), which implies that 
..	. 
v = f . = v . (u1)u1 + ··· + v . (un)un. 
.. 	..
Therefore, (u1,...,u n) spans E. . We claim that the covectors u1,...,u n are linearly inde-pendent. If not, we have a nontrivial linear dependence λ1u1 . + ··· + λnun . =0, and if we apply the above linear form to each ui, using a familar computation, we get 0= λiui . (ui)= λi, 
.. 	..
proving that u1,...,u n are indeed linearly independent. Therefore, (u1,...,u n) is a basis of E. . 
In particular, Theorem 2.20 shows a .nite-dimensional vector space and its dual E. have the same dimension. 
2.9 Summary 
The main concepts and results of this chapter are listed below: 
. 	
The notion of a vector space. 

. 	
Families of vectors. 

. 	
Linear combinations of vectors; linear dependence and linear independence of a family of vectors. 

. 	
Linear subspaces. 

. 	
Spanning (or generating) family; generators, .nitely generated subspace; basis of a subspace. 

. 	
Every linearly independent family can be extended to a basis (Theorem 2.9). 

. 	
A family B of vectors is a basis i. it is a maximal linearly independent family i. it is a minimal generating family (Proposition 2.10). 

. 	
The replacement lemma (Proposition 2.12). 

. 	
Any two bases in a .nitely generated vector space E have the same number of elements; this is the dimension of E (Theorem 2.13). 

. 	
Hyperplanes. 

2.10. PROBLEMS 

. 	
Every vector has a unique representation over a basis (in terms of its coordinates). 

. 	
Matrices 

. 	
Column vectors, row vectors. 

. 	
Matrix operations: addition, scalar multiplication, multiplication. 

. 	
The vector space Mm,n(K) of m × n matrices over the .eld K; The ring Mn(K) of n × n matrices over the .eld K. 

. 	
The notion of a linear map. 

. 	
The image Im f (or range) of a linear map f. 

. 	
The kernel Ker f (or nullspace) of a linear map f. 

. 	
The rank rk(f) of a linear map f. 

. 	
The image and the kernel of a linear map are subspaces. A linear map is injective i. its kernel is the trivial space (0) (Proposition 2.16). 

. 	
The unique homomorphic extension property of linear maps with respect to bases (Proposition 2.17 ). 

. 	
The vector space of linear maps HomK (E, F ). 

. 	
Linear forms (covectors) and the dual space E. . 

. 	
Coordinate forms. 

. 	
The existence of dual bases (in .nite dimension). 


2.10 Problems 
Problem 2.1. Let H be the set of 3 × 3 upper triangular matrices given by 
..  .  .  
.  1  a  b  .  
H = . .0 0  1 0  c 1 . | a, b, c ∈ R . .  

(1) Prove that H with the binary operation of matrix multiplication is a group; .nd explicitly the inverse of every matrix in H. Is H abelian (commutative)? 
(2) Given two groups G1 and G2, recall that a homomorphism if a function .: G1 → G2 
such that .(ab)= .(a).(b), a,b ∈ G1. 
Prove that .(e1)= e2 (where ei is the identity element of Gi) and that .(a .1)=(.(a)).1 ,a ∈ G1. 
(3) Let S1 be the unit circle, that is 
S1 iθ 
= {e = cos θ + i sin θ | 0 ≤ θ< 2π}, 
and let . be the function given by  
.  .  
1  a  b  
.  .0  1  c. = (a, c, e ib).  
0  0  1  

Prove that . is a surjective function onto G = R × R × S1 , and that if we de.ne multiplication on this set by 
ix1y2
(x1,y1,u1) · (x2,y2,u2)=(x1 + x2,y1 + y2,e u1u2), 
then G is a group and . is a group homomorphism from H onto G. 
(4) The kernel of a homomorphism .: G1 → G2 is de.ned as 
Ker (.)= {a ∈ G1 | .(a)= e2}. 
Find explicitly the kernel of . and show that it is a subgroup of H. 
Problem 2.2. For any m ∈ Z with m> 0, the subset mZ = {mk | k ∈ Z} is an abelian subgroup of Z. Check this. 
(1) Give a group isomorphism (an invertible homomorphism) from mZ to Z. 
(2) Check that the inclusion map i: mZ → Z given by i(mk)= mk is a group homomor-phism. Prove that if m ≥ 2 then there is no group homomorphism p: Z → mZ such that p . i = id. 
Remark: The above shows that abelian groups fail to have some of the properties of vector spaces. We will show later that a linear map satisfying the condition p . i = id always exists. 
Problem 2.3. Let E = R× R, and de.ne the addition operation 
(x1,y1)+(x2,y2)=(x1 + x2,y1, +y2),x1,x2,y1,y2 ∈ R, 
and the multiplication operation ·: R× E → E by 
λ · (x, y)=(λx, y), λ,x,y ∈ R. 
Show that E with the above operations + and · is not a vector space. Which of the axioms is violated? 
2.10. PROBLEMS 
Problem 2.4. (1) Prove that the axioms of vector spaces imply that 
α · 0=0 
0 · v =0 α · (.v)= .(α · v) (.α) · v = .(α · v), 
for all v ∈ E and all α ∈ K, where E is a vector space over K. 
(2) For every λ ∈ R and every x =(x1,...,xn) ∈ Rn, de.ne λx by 
λx = λ(x1,...,xn)=(λx1, . . . , λxn). 
Recall that every vector x =(x1,...,xn) ∈ Rn can be written uniquely as 
x = x1e1 + ··· + xnen, 
where ei = (0,..., 0, 1, 0,..., 0), with a single 1 in position i. For any operation ·: R× Rn → Rn, if · satis.es the Axiom (V1) of a vector space, then prove that for any α ∈ R, we have 
α · x = α · (x1e1 + ··· + xnen)= α · (x1e1)+ ··· + α · (xnen). 
Conclude that · is completely determined by its action on the one-dimensional subspaces of Rn spanned by e1,...,en. 
(3) 
Use (2) to de.ne operations ·: R× Rn → Rn that satisfy the Axioms (V1–V3), but for which Axiom V4 fails. 

(4) 
For any operation ·: R× Rn → Rn, prove that if · satis.es the Axioms (V2–V3), then for every rational number r ∈ Q and every vector x ∈ Rn, we have 


r · x = r(1 · x). 
In the above equation, 1 · x is some vector (y1,...,yn) ∈ Rn not necessarily equal to x = (x1,...,xn), and 
r(1 · x)=(ry1, . . . , ryn), 
as in Part (2). 
Use (4) to conclude that any operation ·: Q×Rn → Rn that satis.es the Axioms (V1–V3) is completely determined by the action of 1 on the one-dimensional subspaces of Rn spanned by e1,...,en. 
Problem 2.5. Let A1 be the following matrix: 
.  .  
2  3  1  
A1 = . 1  2  .1..  
.3  .5  1  

Prove that the columns of A1 are linearly independent. Find the coordinates of the vector x = (6, 2, .7) over the basis consisting of the column vectors of A1. 
Problem 2.6. Let A2 be the following matrix: 
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

Express the fourth column of A2 as a linear combination of the .rst three columns of A2. Is the vector x = (7, 14, .1, 2) a linear combination of the columns of A2? 
Problem 2.7. Let A3 be the following matrix: 
.
. 
A3 =
. 

111 
112

.

. 

123 
Prove that the columns of A1 are linearly independent. Find the coordinates of the vector 
x = (6, 9, 14) over the basis consisting of the column vectors of A3. 
Problem 2.8. Let A4 be the following matrix: 

.
. 
A4 = 
... 

1 211 
2 323 

.101 .1 
.2 .14 0 

...

. 

Prove that the columns of A4 are linearly independent. Find the coordinates of the vector x = (7, 14, .1, 2) over the basis consisting of the column vectors of A4. 
Problem 2.9. Consider the following Haar matrix 
.
. 
H = 

... 

11 1 0 
11 .10 

1 .10 1 
1 .10 .1 

...

. 

Prove that the columns of H are linearly independent. Hint. Compute the product HTH. 
Problem 2.10. Consider the following Hadamard matrix 
.
. 
H4 = 
... 

11 1 1 
1 .11 .1 

11 .1 .1 
1 .1 .11 

...

. 

Prove that the columns of H4 are linearly independent. Hint. Compute the product H4 TH4. 
2.10. PROBLEMS 
Problem 2.11. In solving this problem, do not use determinants. 
(1) Let (u1,...,um) and (v1,...,vm) be two families of vectors in some vector space E. Assume that each vi is a linear combination of the ujs, so that 
vi = ai 1u1 + ··· + aimum, 1 ≤ i ≤ m, 
and that the matrix A =(aij) is an upper-triangular matrix, which means that if 1 ≤ j< i ≤ m, then aij = 0. Prove that if (u1,...,um) are linearly independent and if all the diagonal entries of A are nonzero, then (v1,...,vm) are also linearly independent. Hint. Use induction on m. 
(2) 
Let A =(aij) be an upper-triangular matrix. Prove that if all the diagonal entries of A are nonzero, then A is invertible and the inverse A.1 of A is also upper-triangular. Hint. Use induction on m. 

Prove that if A is invertible, then all the diagonal entries of A are nonzero. 

(3) 
Prove that if the families (u1,...,um) and (v1,...,vm) are related as in (1), then (u1,...,um) are linearly independent i. (v1,...,vm) are linearly independent. 


Problem 2.12. In solving this problem, do not use determinants. Consider the n × n 
matrix

.
. 
120 0 ... 00 012 0 ... 00 001 2 ... 00 
.. ..
...
.......... 

.......... 

A = 

.. 

.

.

... 

.

. 

. 

.

.. 

.. 

0  0  . . .  0  1  2  0  
0  0  . . .  0  0  1  2  
0  0  . . .  0  0  0  1  

(1) Find the solution x =(x1,...,xn) of the linear system Ax = b, 
for

.
. 

.... 

b1 b2 
. 
. 
. 
bn 
.... 

b = 

. 

(2) Prove that the matrix A is invertible and .nd its inverse A.1 . Given that the number of atoms in the universe is estimated to be ≤ 1082, compare the size of the coe.cients the inverse of A to 1082, if n ≥ 300. 
(3) Assume b is perturbed by a small amount Δb (note that Δb is a vector). Find the 
new solution of the system A(x +Δx)= b +Δb, 
where Δx is also a vector. In the case where b = (0,..., 0, 1), and Δb = (0,..., 0,ε), show that 
|(Δx)1| =2n.1|ε|. 
(where (Δx)1 is the .rst component of Δx). 
(4) Prove that (A . I)n = 0. 
Problem 2.13. An n × n matrix N is nilpotent if there is some integer r ≥ 1 such that 
Nr 
= 0. 
(1) 
Prove that if N is a nilpotent matrix, then the matrix I . N is invertible and 

(I . N).1 = I + N + N2 + ··· + Nr.1 . 

(2) 
Compute the inverse of the following matrix A using (1): 


.
. 
A = 

..... 

1  2  3  4  5  
0  1  2  3  4  
0  0  1  2  3  

00012 
00001 

..... 

. 

Problem 2.14. (1) Let A be an n × n matrix. If A is invertible, prove that for any x ∈ Rn , if Ax = 0, then x = 0. 
(2) Let A be an m × n matrix and let B be an n × m matrix. Prove that Im . AB is invertible i. In . BA is invertible. Hint. If for all x ∈ Rn , Mx = 0 implies that x = 0, then M is invertible. 
Problem 2.15. Consider the following n × n matrix, for n ≥ 3: 
.
. 
B = 

.......... 

1 .1 .1 .1 ··· .1 .1 1 .11 1 ··· 11 11 .11 ··· 11 11 1 .1 ··· 11 
.  .  .  .  .  .  .  
.  .  .  .  .  .  .  
.  .  .  .  .  .  .  
1  1  1  1  · · ·  .1  1  
1  1  1  1  · · ·  1  .1  

.......... 

(1) If we denote the columns of B by b1,...,bn, prove that 
(n . 3)b1 . (b2 + ··· + bn) = 2(n . 2)e1 b1 . b2 = 2(e1 + e2) b1 . b3 = 2(e1 + e3) 
.. 
.. 
.. b1 . bn = 2(e1 + en), 
2.10. PROBLEMS 
where e1,...,en are the canonical basis vectors of Rn . 
(2) Prove that B is invertible and that its inverse A =(aij) is given by 
(n . 3) 1 
a11 = ,ai1 = . 2 ≤ i ≤ n 
2(n . 2)2(n . 2) 
and  
aii  = . (n . 3) 2(n . 2),  2 ≤ i ≤ n  
1  
aji  =  2(n . 2),  2 ≤ i ≤ n, j  = i.  

(3) Show that the n diagonal n × n matrices Di de.ned such that the diagonal entries of Di are equal the entries (from top down) of the ith column of B form a basis of the space of n × n diagonal matrices (matrices with zeros everywhere except possibly on the diagonal). For example, when n = 4, we have 
.
.
.
. 
1000 .1 0 00 0
...
... 

... 

...

.1

010 

0 

00 

D1 
= 

D2 
= 

,

0 0001 0 001 
001 

0 010 

.
.
.
. 
.10 0 0 .100 0 

...
0001 000 .1 
Problem 2.16. Given any m×n matrix A and any n×p matrix B, if we denote the columns of A by A1,...,An and the rows of B by B1,...,Bn, prove that 
AB = A1B1 + ··· + AnBn. 
Problem 2.17. Let f : E → F be a linear map which is also a bijection (it is injective and surjective). Prove that the inverse function f.1 : F → E is linear. 
Problem 2.18. Given two vectors spaces E and F , let (ui)i∈I be any basis of E and let (vi)i∈I be any family of vectors in F . Prove that the unique linear map f : E → F such that f(ui)= vi for all i ∈ I is surjective i. (vi)i∈I spans F . 
Problem 2.19. Let f : E → F be a linear map with dim(E)= n and dim(F )= m. Prove that f has rank 1 i. f is represented by an m × n matrix of the form 
T
A = uv 
... 

... 

...

010 

0 00 .10 
0 10 0 

D3 D4 
= 

=

, 

. 

0 01 0 

with u a nonzero column vector of dimension m and v a nonzero column vector of dimension 
n. 
Problem 2.20. Find a nontrivial linear dependence among the linear forms 
.1(x, y, z)=2x . y +3z, .2(x, y, z)=3x . 5y + z, .3(x, y, z)=4x . 7y + z. Problem 2.21. Prove that the linear forms 
.1(x, y, z)= x +2y + z, .2(x, y, z)=2x +3y +3z, .3(x, y, z)=3x +7y + z are linearly independent. Express the linear form .(x, y, z)= x+y+z as a linear combination of .1,.2,.3. 




