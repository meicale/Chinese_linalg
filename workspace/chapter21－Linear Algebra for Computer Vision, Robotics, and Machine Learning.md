Chapter 21 
Applications of SVD and Pseudo-Inverses 
De tous les principes qu’on peut proposer pour cet objet, je pense qu’il n’en est pas de plus g▽en▽eral, de plus exact, ni d’une application plus facile, que celui dont nous avons fait usage dans les recherches pr▽ec▽edentes, et qui consiste `a rendre minimum la somme des carr▽es des erreurs. Par ce moyen il s’▽etablit entre les erreurs une sorte d’▽equilibre qui, emp.echant les extr.evaloir, est tr``etat du 
emes de pr▽es propre as faire connaitre l’▽
syst`eme le plus proche de la v▽erit▽e. 

―Legendre, 1805, Nouvelles M▽ethodes pour la d▽etermination des Orbites des Com`etes 
21.1 Least Squares Problems and the Pseudo-Inverse 
This chapter presents several applications of SVD. The .rst one is the pseudo-inverse, which plays a crucial role in solving linear systems by the method of least squares. The second ap-plication is data compression. The third application is principal component analysis (PCA), whose purpose is to identify patterns in data and understand the variance쭯covariance struc-ture of the data. The fourth application is the best a.ne approximation of a set of data, a problem closely related to PCA. 
The method of least squares is a way of “solving” an overdetermined system of linear equations Ax = b, 
i.e., a system in which A is a rectangular m ≠ n matrix with more equations than unknowns (when m>n). Historically, the method of least squares was used by Gauss and Legendre to solve problems in astronomy and geodesy. The method was .rst published by Legendre in 1805 in a paper on methods for determining the orbits of comets. However, Gauss had already used the method of least squares as early as 1801 to determine the orbit of the asteroid 
659 

Ceres, and he published a paper about it in 1810 after the discovery of the asteroid Pallas. Incidentally, it is in that same paper that Gaussian elimination using pivots is introduced. 
The reason why more equations than unknowns arise in such problems is that repeated measurements are taken to minimize errors. This produces an overdetermined and often in-consistent system of linear equations. For example, Gauss solved a system of eleven equations in six unknowns to determine the orbit of the asteroid Pallas. 
Example 21.1. As a concrete illustration, suppose that we observe the motion of a small object, assimilated to a point, in the plane. From our observations, we suspect that this point moves along a straight line, say of equation y = dx + c. Suppose that we observed the moving point at three di.erent locations (x1,y1), (x2,y2), and (x3,y3). Then we should have 
c + dx1 = y1, c + dx2 = y2, c + dx3 = y3. 
If there were no errors in our measurements, these equations would be compatible, and c and d would be determined by only two of the equations. However, in the presence of errors, the system may be inconsistent. Yet we would like to .nd c and d! 
The idea of the method of least squares is to determine (c, d) such that it minimizes the sum of the squares of the errors, namely, 
(c + dx1 . y1)2 +(c + dx2 . y2)2 +(c + dx3 . y3)2 . 
See Figure 21.1. 

21.1. LEAST SQUARES PROBLEMS AND THE PSEUDO-INVERSE 
In general, for an overdetermined m ≠ n system Ax = b, what Gauss and Legendre discovered is that there are solutions x minimizing 
lAx . bl22 
(where lul22 = u21 +···+u2 n, the square of the Euclidean norm of the vector u =(u1,...,un)), and that these solutions are given by the square n ≠ n system 
ATAx = ATb, 
called the normal equations. Furthermore, when the columns of A are linearly independent, it turns out that ATA is invertible, and so x is unique and given by 
x =(ATA).1ATb. 
Note that ATA is a symmetric matrix, one of the nice features of the normal equations of a least squares problem. For instance, since the above problem in matrix form is represented 
as .. .. 1 x1   y1 ..c ..
1 x2 = y2 ,
d
1 x3 y3 
the normal equations are
      
3 x1 + x2 + x3 cy1 + y2 + y3 
=.
222
x1 + x2 + x3 x1 + x2 + x3dx1y1 + x2y2 + x3y3
In fact, given any real m ≠ n matrix A, there is always a unique x+ of minimum norm that minimizes lAx . bl22, even when the columns of A are linearly dependent. How do we prove this, and how do we .nd x+? 
Theorem 21.1. Every linear system Ax = b, where A is an m ≠ n matrix, has a unique least squares solution x+ of smallest norm. 
Proof. Geometry o.ers a nice proof of the existence and uniqueness of x+ . Indeed, we can interpret b as a point in the Euclidean (a.ne) space Rm, and the image subspace of A (also called the column space of A) as a subspace U of Rm (passing through the origin). Then it is clear that 
inf lAx . bl22 = inf ly . bl22 , 
xÅRn yÅU 
with U = Im A, and we claim that x minimizes lAx . bl22 i. Ax = p, where p the orthogonal projection of b onto the subspace U. Recall from Section 12.1 that the orthogonal projection pU : U 쮵 U￥ ∪ U is the linear map given by pU (u + v)= u, 
with u Å U and v Å U￥ . If we let p = pU (b) Å U, then for any point y Å U, the vectors 
.
∪.
∪ 
py = y . p Å U and bp = p . b Å U￥ are orthogonal, which implies that .∪ .
∪ 
∪ 
2 + l.2,
lbyl22 = lbpl2 pyl2 
.
∪ 
where by = y . b. Thus, p is indeed the unique point in U that minimizes the distance from b to any point in U. See Figure 21.2. 

Thus the problem has been reduced to proving that there is a unique x+ of minimum norm such that Ax+ = p, with p = pU (b) Å U, the orthogonal projection of b onto U. We use the fact that 
Rn = Ker A 쮵 (Ker A)￥ . 
Consequently, every x Å Rn can be written uniquely as x = u + v, where u Å Ker A and v Å (Ker A)￥, and since u and v are orthogonal, 
lxl22 = lul22 + lvl22. 
Furthermore, since u Å Ker A, we have Au = 0, and thus Ax = p i. Av = p, which shows that the solutions of Ax = p for which x has minimum norm must belong to (Ker A)￥ . However, the restriction of A to (Ker A)￥ is injective. This is because if Av1 = Av2, where v1,v2 Å (Ker A)￥, then A(v2 . v2) = 0, which implies v2 . v1 Å Ker A, and since v1,v2 Å (Ker A)￥, we also have v2 . v1 Å (Ker A)￥, and consequently, v2 . v1 = 0. This shows that there is a unique x+ of minimum norm such that Ax+ = p, and that x+ must belong to (Ker A)￥ . By our previous reasoning, x+ is the unique vector of minimum norm minimizing lAx . bl22. 

21.1. LEAST SQUARES PROBLEMS AND THE PSEUDO-INVERSE 
.
∪ 
The proof also shows that x minimizes lAx . bl22 i. pb = b . Ax is orthogonal to U, which can be expressed by saying that b . Ax is orthogonal to every column of A. However, this is equivalent to 
AT(b . Ax)=0, i.e., ATAx = ATb. 
Finally, it turns out that the minimum norm least squares solution x+ can be found in terms of the pseudo-inverse A+ of A, which is itself obtained from any SVD of A. 
De.nition 21.1. Given any nonzero m ≠ n matrix A of rank r, if A = V DUT is an SVD 
of A such that  
D =  ┼ 0m.r,r  0r,n.r 0m.r,n.r  ,  
with  

┼ = diag(┡1,...,┡r) 
an r ≠ r diagonal matrix consisting of the nonzero singular values of A, then if we let D+ be the n ≠ m matrix 
┼.1 
0r,m.r
D+ 
=0n.r,r 0n.r,m.r , with ┼.1 = diag(1/┡1,..., 1/┡r), 
the pseudo-inverse of A is de.ned by 
A+ = UD+V T . 
If A =0m,n is the zero matrix, we set A+ =0n,m. Observe that D+ is obtained from D by inverting the nonzero diagonal entries of D, leaving all zeros in place, and then transposing the matrix. For example, given the matrix 
.
. 
10000 
02000 

00300 
00000 

...

...

D = 

, 

its pseudo-inverse is 

.
. 
1000 
1
0 00
2 
1 
..... 

..... 

D+ 
00 

0

= 

.

3 
0000 
0000 

The pseudo-inverse of a matrix is also known as the Moore쭯Penrose pseudo-inverse. Actually, it seems that A+ depends on the speci.c choice of U and V in an SVD (U, D, V ) for A, but the next theorem shows that this is not so. 
Theorem 21.2. The least squares solution of smallest norm of the linear system Ax = b, where A is an m ≠ n matrix, is given by 
x + = A+b = UD+V Tb. Proof. First assume that A is a (rectangular) diagonal matrix D, as above. Then since x minimizes lDx . bl22 i. Dx is the projection of b onto the image subspace F of D, it is fairly obvious that x+ = D+b. Otherwise, we can write 
A = V DUT , where U and V are orthogonal. However, since V is an isometry, lAx . bl2 = lV DUT x . bl2 = lDUT x . V Tbl2. Letting y = UTx, we have lxl2 = lyl2, since U is an isometry, and since U is surjective, lAx . bl2 is minimized i. lDy . V Tbl2 is minimized, and we have shown that the least solution is y + = D+V Tb. Since y = UTx, with lxl2 = lyl2, we get x + = UD+V Tb = A+b. Thus, the pseudo-inverse provides the optimal solution to the least squares problem. By Theorem 21.2 and Theorem 21.1, A+b is uniquely de.ned by every b, and thus A+ depends only on A. The Matlab command for computing the pseudo-inverse B of the matrix A is B = pinv(A). Example 21.2. If A is the rank 2 matrix 
.
. 
A = 

... 

1234 
2345 

3456 
4567 

... 

whose eigenvalues are .1.1652, 0, 0, 17.1652, using Matlab we obtain the SVD A = V DUT with 
.
. 
... 

.0.3147 0.7752 0.2630 .0.4805 
.0.4275 0.3424 0.0075 0.8366 

.0.5402 .0.0903 .0.8039 .0.2319 
.0.6530 .0.5231 0.5334 .0.1243 

...

,

U = 

.
.
.
. 
.0.3147 .0.7752 0.5452 0.0520 17.1652 0 00 

... 

.0.4275 .0.3424 .0.7658 0.3371 

.0.5402 0.0903 .0.1042 .0.8301 

...

, 

D = 

... 

01.1652 0 0 

0 000 

...

= 

. 

.0.6530 0.5231 0.3247 0.4411 0 000 


21.1. LEAST SQUARES PROBLEMS AND THE PSEUDO-INVERSE 
Then

.
. 
0.0583 0 00 00.8583 0 0 
D+ 
= 
...

0 0 000 
...
0 00 

, 

and

.
. 
A+ 
= UD+V T = 
... 

.0.5100 .0.2200 0.0700 0.3600 
.0.2200 .0.0900 0.0400 0.1700 

0.0700 0.0400 0.0100 .0.0200 
0.3600 0.1700 .0.0200 .0.2100 

...

, 

which is also the result obtained by calling pinv(A). If A is an m ≠ n matrix of rank n (and so m ◎ n), it is immediately shown that the QR-decomposition in terms of Householder transformations applies as follows: There are nm ≠ m matrices H1,...,Hn, Householder matrices or the identity, and an upper triangular m ≠ n matrix R of rank n such that A = H1 ··· HnR. Then because each Hi is an isometry, lAx . bl2 = lRx . Hn ··· H1bl2, and the least squares problem Ax = b is equivalent to the system Rx = Hn ··· H1b. Now the system Rx = Hn ··· H1b 
is of the form  
R1  c  
0m.n  x =  d  ,  

where R1 is an invertible n ≠ n matrix (since A has rank n), c Å Rn, and d Å Rm.n, and the least squares solution of smallest norm is 
+ = R.1 
x c.
1 
Since R1 is a triangular matrix, it is very easy to invert R1. 
The method of least squares is one of the most e.ective tools of the mathematical sciences. There are entire books devoted to it. Readers are advised to consult Strang [64], Golub and Van Loan [30], Demmel [16], and Trefethen and Bau [68], where extensions and applications of least squares (such as weighted least squares and recursive least squares) are described. Golub and Van Loan [30] also contains a very extensive bibliography, including a list of books on least squares. 


21.2 Properties of the Pseudo-Inverse 
We begin this section with a proposition which provides a way to calculate the pseudo-inverse of an m ≠ n matrix A without .rst determining an SVD factorization. 
Proposition 21.3. When A has full rank, the pseudo-inverse A+ can be expressed as A+ = (ATA).1AT when m ◎ n, and as A+ = AT(AAT).1 when n ◎ m. In the .rst case (m ◎ n), observe that A+A = I, so A+ is a left inverse of A; in the second case (n ◎ m), we have AA+ = I, so A+ is a right inverse of A. 
Proof. If m ◎ n and A has full rank n, we have 
┼ 
UT
A = V 0m.n,n 
with ┼ an n ≠ n diagonal invertible matrix (with positive entries), so 
  
A+ ┼.1 V T 
= U0n,m.n. 
We .nd that 
  
ATA = U┼0n,m.nV TV ┼ UT = U┼2UT ,0m.n,n 
which yields 
    
(ATA).1AT = U┼.2UTU┼0n,m.nV T = U┼.1 0n,m.nV T = A+ . 
Therefore, if m ◎ n and A has full rank n, then 
A+ =(ATA).1AT . 
If n ◎ m and A has full rank m, then 
  
A = V┼0m,n.mUT 
with ┼ an m ≠ m diagonal invertible matrix (with positive entries), so 
┼.1 
A+ V T 
= U. 0n.m,m 
We .nd that 
  
AAT = V┼0m,n.mUTU ┼ V T = V ┼2V T ,0n.m,m 
which yields 
AT(AAT).1 = U ┼ V TV ┼.2V T = U ┼.1 V T = A+ . 0n.m,m 0n.m,m 
Therefore, if n ◎ m and A has full rank m, then A+ = AT(AAT).1 . 
21.2. PROPERTIES OF THE PSEUDO-INVERSE 
.
. 

12 

For example, if A =

.

23., then A has rank 2 and since m ◎ n, A+ =(ATA).1AT 
0  1  
where  
A+ =  5 8  8 14  .1 AT =  7/3 4/3  .4/3 5/6  1 2  2 3  0 1  =  .1/3 1/3  2/3 .1/6  .4/3 5/6  .  

123 0 
If A = , since A has rank 2 and n ◎ m, then A+ = AT(AAT).1 where 
011 .1 
.
.
.. 
10 3/17 .5/17 

.1 
= 

... 

21 

31 

... 

= 

... 

1/17 4/17 

4/17 

.1/17 

...

. 

14 5 

3/17 

.5/17 

A+ = AT 
53 

.5/17 

14/17 

0 .15/17 .14/17 Let A = V ┣UT be an SVD for any m ≠ n matrix A. It is easy to check that both AA+ and A+A are symmetric matrices. In fact, Ir 0 
V TAA+ = V ┣UTU┣+V T = V ┣┣+V T = V 
00m.r 
and 
Ir 0 
A+A = U┣+V TV ┣UT = U┣+┣UT = UUT . 
00n.r From the above expressions we immediately deduce that AA+A = A, A+AA+ = A+ , 
and that (AA+)2 = AA+ , (A+A)2 = A+A, 
so both AA+ and A+A are orthogonal projections (since they are both symmetric). 
Proposition 21.4. The matrix AA+ is the orthogonal projection onto the range of A and A+A is the orthogonal projection onto Ker(A)￥ = Im(AT), the range of AT . Proof. Obviously, we have range(AA+) . range(A), and for any y = Ax Å range(A), since 
AA+A = A, we have 
AA+ y = AA+Ax = Ax = y, so the image of AA+ is indeed the range of A. It is also clear that Ker(A) . Ker(A+A), and since AA+A = A, we also have Ker(A+A) . Ker(A), and so 
Ker(A+A) = Ker(A). Since A+A is symmetric, range(A+A) = range((A+A)T) = Ker(A+A)￥ = Ker(A)￥ , as claimed. 
Proposition 21.5. The set range(A) = range(AA+) consists of all vectors y Å Rm such that 
z 
V T y = ,
0 
with z Å Rr . 
Proof. Indeed, if y = Ax, then 

V T y = V TAx = V TV ┣UT x =┣UT x =┣r 0 UT x = z,
00m.r 0 
where ┣r is the r ≠ r diagonal matrix diag(┮1,...,┮r). Conversely, if V Ty =( z 0 ), then y = V ( z 0 ), and 
Ir 0 
AA+ y = VV T y
00m.r Ir 0 z 
= VV TV 
00m.r 0 Ir 0 z = V 00m.r 0 
z 
= V = y, 
0 which shows that y belongs to the range of A. 
Similarly, we have the following result. Proposition 21.6. The set range(A+A) = Ker(A)￥ consists of all vectors y Å Rn such that 
z 
UT y = ,
0 
with z Å Rr . Proof. If y = A+Au, then Ir 0 z 
y = A+Au = UUT u = U,
00n.r 0 
for some z Å Rr . Conversely, if UTy =( z ), then y = U ( z ), and so 
0  0  
A+AU  z 0  = U  Ir 0  0 0n.r  UTU  z 0  
= U  Ir 0  0 0n.r  z 0  
= U  z 0  = y,  

which shows that y Å range(A+A). 

21.2. PROPERTIES OF THE PSEUDO-INVERSE 
Analogous results hold for complex matrices, but in this case, V and U are unitary matrices and AA+ and A+A are Hermitian orthogonal projections. 
If A is a normal matrix, which means that AAT = ATA, then there is an intimate relationship between SVD’s of A and block diagonalizations of A. As a consequence, the pseudo-inverse of a normal matrix A can be obtained directly from a block diagonalization of A. 
If A is a (real) normal matrix, then we know from Theorem 16.18 that A can be block diagonalized with respect to an orthogonal matrix U as 
A = U┼UT , 
where ┼ is the (real) block diagonal matrix 
┼ = diag(B1,...,Bn), 
consisting either of 2 ≠ 2 blocks of the form 
┡j .┢j
Bj = 
┢j ┡j 
with ┢j = 0, or of one-dimensional blocks Bk =(┡k). Then we have the following proposition: 
Proposition 21.7. For any (real) normal matrix A and any block diagonalization A = U┼UT of A as above, the pseudo-inverse of A is given by 
A+ = U┼+UT , 
where ┼+ is the pseudo-inverse of ┼. Furthermore, if 
┼r 0 
┼= ,
00 
where ┼r has rank r, then ┼.1 
┼+ r 0 
= . 
00 
Proof. Assume that B1,...,Bp are 2 ≠ 2 blocks and that ┡2p+1,...,┡n are the scalar entries. We know that the numbers ┡j ÷ i┢j, and the ┡2p+k are the eigenvalues of A. Let ┭2j.1 =
o  
┭2j = ┡j 2 + ┢j 2 =det(Bi) for j =1,...,p, ┭j = |┡j| for j =2p +1,...,r. Multiplying U by a suitable permutation matrix, we may assume that the blocks of ┼ are ordered so that ┭1 ◎ ┭2 ◎· ·· ◎ ┭r > 0. Then it is easy to see that 
AAT = ATA = U┼UTU┼TUT = U┼┼TUT , 
with 
┼┼T = diag(┭21,...,┭2 r, 0,..., 0), 
so ┭1 ◎ ┭2 ◎ · ·· ◎ ┭r > 0 are the singular values ┮1 ◎ ┮2 ◎ · ·· ◎ ┮r > 0 of A. De.ne the diagonal matrix 
┣ = diag(┮1,...,┮r, 0,..., 0), 
where r = rank(A), ┮1 ◎ · ·· ◎ ┮r > 0 and the block diagonal matrix ┬ de.ned such that the block Bi in ┼ is replaced by the block ┮.1Bi where ┮ = det(Bi), the nonzero scalar ┡j is replaced ┡j/|┡j|, and a diagonal zero is replaced by 1. Observe that ┬ is an orthogonal matrix and 
┼ = ┬┣. 
But then we can write A = U┼UT = U┬┣UT , 
and we if let V = U┬, since U is orthogonal and ┬ is also orthogonal, V is also orthogonal and A = V ┣UT is an SVD for A. Now we get 
A+ = U┣+┬TUT 
= U┣+V T . 
However, since ┬ is an orthogonal matrix, ┬T =┬.1, and a simple calculation shows that 
┣+┬T =┣+┬.1 =┼+ , 
which yields the formula 
A+ = U┼+UT 
. 
Also observe that ┼r is invertible and 
┼.1 
0 ┼+ r
= . 
00 
Therefore, the pseudo-inverse of a normal matrix can be computed directly from any block diagonalization of A, as claimed. 
Example 21.3. Consider the following real diagonal form of the normal matrix 
.
. 
A = 

... 

.2.7500 2.1651 .0.8660 0.5000 
2.1651 .0.2500 .1.5000 0.8660 

0.8660 1.5000 0.7500 .0.4330 .0.5000 .0.8660 .0.4330 0.2500 
... 

= U┼UT 
, 

with 

.
.
.
. 
cos(┪/3) 0 sin(┪/3)0 1 .200 

U 

= 

... 

sin(┪/3) 0 . cos(┪/3) 0 

0 cos(┪/6) 0 sin(┪/6) 

...

, 

┼= 

... 

21 
0 
00 

.4 

0 
0 

. 

...
0 . cos(┪/6) 0 sin(┪/6) 0000 

21.3. DATA COMPRESSION AND SVD 
We obtain 

.
. 
1/52/500 
.2/51/500

... 

...

,

┼+ 
= 
00 .1/40 
0 0 00 

and the pseudo-inverse of A is 

.
. 
.0.1375 0.1949 0.1732 .0.1000 .0.1732
...
.0.0866 0.1000 0.1732 .0.0866 0.0500 
which agrees with pinv(A). 
The following properties, due to Penrose, characterize the pseudo-inverse of a matrix. We have already proved that the pseudo-inverse satis.es these equations. For a proof of the converse, see Kincaid and Cheney [39]. 
Proposition 21.8. Given any m ≠ n matrix A (real or complex), the pseudo-inverse A+ of A is the unique n ≠ m matrix satisfying the following properties: AA+A = A, A+AA+ = A+ , (AA+)T = AA+ , (A+A)T = A+A. 


21.3 Data Compression and SVD 
Among the many applications of SVD, a very useful one is data compression, notably for images. In order to make precise the notion of closeness of matrices, we use the notion of matrix norm. This concept is de.ned in Chapter 8, and the reader may want to review it before reading any further. 
Given an m ≠ n matrix of rank r, we would like to .nd a best approximation of A by a matrix B of rank k ● r (actually, k<r) such that lA . Bl2 (or lA . BlF ) is minimized. The following proposition is known as the Eckart쭯Young theorem. 
Proposition 21.9. Let A be an m ≠ n matrix of rank r and let V DUT = A be an SVD for 
A. Write ui for the columns of U, vi for the columns of V , and ┮1 ◎ ┮2 ◎ · ·· ◎ ┮p for the singular values of A (p = min(m, n)). Then a matrix of rank k<r closest to A (in the ll2 norm) is given by 
... 

0.1949 0.0875 0.3000 

A+ = U┼+UT 
= 
.0.1732 .0.3000 

0.1500 

, 

k
Ak = ┮iviu T i = V diag(┮1,...,┮k, 0,..., 0)UT 
i=1 
and lA . Akl2 = ┮k+1. 
Proof. By construction, Ak has rank k, and we have 
  
  
  
p 
lA . Akl
2 
=

┮iviu 
T 
i
2 
i=k+1 
=

  

V diag(0,..., 0,┮k+1,...,┮p)UT
  

2 
= ┮k+1. 
It remains to show that lA . Bl2 ◎ ┮k+1 for all rank k matrices B. Let B be any rank k matrix, so its kernel has dimension n . k. The subspace Uk+1 spanned by (u1,...,uk+1) has dimension k + 1, and because the sum of the dimensions of the kernel of B and of Uk+1 is (n . k)+ k +1 = n + 1, these two subspaces must intersect in a subspace of dimension at least 1. Pick any unit vector h in Ker(B) ℃ Uk+1. Then since Bh = 0, and since U and V are isometries, we have 
lA . Bl2 ◎l(A . B)hl2 = lAhl2 
2 22 
=

  

V DUTh
  

2 
2 
=

  

DUTh
  

2 
2 ◎ ┮k2+1
  

UTh
  

2 
2 = ┮k2+1, 
which proves our claim. 
Note that Ak can be stored using (m + n)k entries, as opposed to mn entries. When k《 m, this is a substantial gain. Example 21.4. Consider the badly conditioned symmetric matrix 
.
. 
A = 

... 

1078 7 
756 5 

8 610 9 

... 

7 5 9 10 from Section 8.5. Since A is SPD, we have the SVD A = UDUT , with 
.
. 
... 

.0.5286 .0.6149 0.3017 .0.5016 
.0.3803 .0.3963 .0.0933 0.8304 

.0.5520 0.2716 .0.7603 .0.2086 
.0.5209 0.6254 0.5676 0.1237 

...

,

U = 

.
. 
... 

30.28870 0 0 03.85810 0 
0 00.8431 0 
0 000.0102 

...

D = 

. 

21.4. PRINCIPAL COMPONENTS ANALYSIS (PCA) 
If we set ┮3 = ┮4 = 0, we obtain the best rank 2 approximation 
.
. 
A2 = U(:, 1 : 2) . D(:, 1 : 2) . U(:, 1 : 2)I = 
... 

9.9207 7.0280 8.1923 6.8563 
7.0280 4.9857 5.9419 5.0436 

8.1923 5.9419 9.5122 9.3641 
6.8563 5.0436 9.3641 9.7282 

...

. 

A nice example of the use of Proposition 21.9 in image compression is given in Demmel [16], Chapter 3, Section 3.2.3, pages 113쭯115; see the Matlab demo. 
Proposition 21.9 also holds for the Frobenius norm; see Problem 21.4. 
An interesting topic that we have not addressed is the actual computation of an SVD. This is a very interesting but tricky subject. Most methods reduce the computation of an SVD to the diagonalization of a well-chosen symmetric matrix which is not ATA; see Problem 20.1 and Problem 20.3. Interested readers should read Section 5.4 of Demmel’s excellent book [16], which contains an overview of most known methods and an extensive list of references. 


21.4 Principal Components Analysis (PCA) 
Suppose we have a set of data consisting of n points X1,...,Xn, with each Xi Å Rd viewed as a row vector. Think of the Xi’s as persons, and if Xi =(xi 1,...,xid), each xij is the value of some feature (or attribute) of that person. 
Example 21.5. For example, the Xi’s could be mathematicians, d = 2, and the .rst com-ponent, xi 1, of Xi could be the year that Xi was born, and the second component, xi 2, the length of the beard of Xi in centimeters. Here is a small data set. 
Name  year  length  
Carl Friedrich Gauss  1777  0  
Camille Jordan  1838  12  
Adrien-Marie Legendre  1752  0  
Bernhard Riemann  1826  15  
David Hilbert  1862  2  
Henri Poincar▽e  1854  5  
Emmy Noether  1882  0  
Karl Weierstrass  1815  0  
Eugenio Beltrami  1835  2  
Hermann Schwarz  1843  20  

We usually form the n ≠ d matrix X whose ith row is Xi, with 1 ● i ● n. Then the jth column is denoted by Cj (1 ● j ● d). It is sometimes called a feature vector, but this terminology is far from being universally accepted. In fact, many people in computer vision call the data points Xi feature vectors! 
The purpose of principal components analysis, for short PCA, is to identify patterns in data and understand the variance쭯covariance structure of the data. This is useful for the following tasks: 
1. 
Data reduction: Often much of the variabi lity of the data can be accounted for by a smaller number of principal components. 

2. 
Interpretation: PCA can show relationships that were not previously suspected. 


Given a vector (a sample of measurements) x =(x1,...,xn) Å Rn, recall that the mean (or average) x of x is given by 
 
n 
xi
i=1 
x =. 
n We let x . x denote the centered data point 
x . x =(x1 . x, . . . , xn . x). 
In order to measure the spread of the xi’s around the mean, we de.ne the sample variance (for short, variance) var(x) (or s2) of the sample x by 
 n 
(xi . x)2 
var(x)=i=1. 
n . 1 
1+3.1
Example 21.6. If x = (1, 3, .1), x = 3 = 1, x . x = (0, 2, .2), and var(x)= 
02+22+(.2)2 1+2+3 (.1)2+02+12 
= 4. If y = (1, 2, 3), y = =2, y .y =(.1, 0, 1), and var(y)= = 
23 2 
2. 
There is a reason for using n . 1 instead of n. The above de.nition makes var(x) an unbiased estimator of the variance of the random variable being sampled. However, we don’t need to worry about this. Curious readers will .nd an explanation of these peculiar de.nitions in Epstein [20] (Chapter 14, Section 14.5) or in any decent statistics book. 
Given two vectors x =(x1,...,xn) and y =(y1,...,yn), the sample covariance (for short, covariance) of x and y is given by 
 n 
i=1(xi . x)(yi . y)
cov(x, y)=. 
n . 1 
Example 21.7. If we take x = (1, 3, .1) and y = (0, 2, .2), we know from Example 21.6 
0(.1)+2(0)+(.2)(1)
that x . x = (0, 2, .2) and y . y =(.1, 0, 1). Thus, cov(x, y)= 2 = .1. 
The covariance of x and y measures how x and y vary from the mean with respect to each other. Obviously, cov(x, y) = cov(y, x) and cov(x, x) = var(x). 
21.4. PRINCIPAL COMPONENTS ANALYSIS (PCA) 
Note that (x . x)T(y . y)
cov(x, y)= . 
n . 1 
We say that x and y are uncorrelated i. cov(x, y) = 0. 
Finally, given an n ≠ d matrix X of n points Xi, for PCA to be meaningful, it will be necessary to translate the origin to the centroid (or center of gravity) ┢ of the Xi’s, de.ned by 
1 
┢ =(X1 + ··· + Xn). 
n 
Observe that if ┢ =(┢1,...,┢d), then ┢j is the mean of the vector Cj (the jth column of X). 
We let X . ┢ denote the matrix whose ith row is the centered data point Xi . ┢ (1 ● i ● n). Then the sample covariance matrix (for short, covariance matrix) of X is the d ≠ d symmetric matrix 
┣= 1(X . ┢)T(X . ┢) = (cov(Ci,Cj)). 
n . 1
.. 
11 Example 21.8. Let X = . 32., the 3 ≠ 2 matrix whose columns are the vector x and .13 
y of Example 21.6. Then 
1 
┢ = [(1, 1) + (3, 2) + (.1, 3)] = (1, 2),
3 
.. 
0 .1 X . ┢ = . 20 . , .21 
and  .  .  
┣ =  1 2  0 .1  2 0  .2 1  . 0 2 .2  .1 0 1 . =  4 .1  .1 1  .  

1
Remark: The factor n.1 is irrelevant for our purposes and can be ignored. 
Example 21.9. Here is the matrix X . ┢ in the case of our bearded mathematicians: since ┢1 = 1828.4,┢2 =5.6, 
we get the following centered data set. 

See Figure 21.3. 

Name  year  length  
Carl Friedrich Gauss  .51.4  .5.6  
Camille Jordan  9.6  6.4  
Adrien-Marie Legendre  .76.4  .5.6  
Bernhard Riemann  .2.4  9.4  
David Hilbert  33.6  .3.6  
Henri Poincar▽e  25.6  .0.6  
Emmy Noether  53.6  .5.6  
Karl Weierstrass  13.4  .5.6  
Eugenio Beltrami  6.6  .3.6  
Hermann Schwarz  14.6  14.4  


We can think of the vector Cj as representing the features of X in the direction ej (the jth canonical basis vector in Rd, namely ej = (0,..., 1,... 0), with a 1 in the jth position). 
If v Å Rd is a unit vector, we wish to consider the projection of the data points X1,...,Xn onto the line spanned by v. Recall from Euclidean geometry that if x Å Rd is any vector and v Å Rd is a unit vector, the projection of x onto the line spanned by v is
（x, v령v. 

21.4. PRINCIPAL COMPONENTS ANALYSIS (PCA) 
Thus, with respect to the basis v, the projection of x has coordinate（x, v령. If x is represented by a row vector and v by a column vector, then
（x, v령 = xv. Therefore, the vector Y Å Rn consisting of the coordinates of the projections of X1,...,Xn onto the line spanned by v is given by Y = Xv, and this is the linear combination 
Xv = v1C1 + ··· + vdCd of the columns of X (with v =(v1,...,vd)). Observe that because ┢j is the mean of the vector Cj (the jth column of X), we get Y = Xv = v1┢1 + ··· + vd┢d, and so the centered point Y . Y is given by Y . Y = v1(C1 . ┢1)+ ··· + vd(Cd . ┢d)=(X . ┢)v. Furthermore, if Y = Xv and Z = Xw, then ((X . ┢)v)T(X . ┢)w 
cov(Y, Z)  =  n . 1  
=  v T 1 n . 1(X . ┢)T(X . ┢)w  
=  v T┣w,  

where ┣ is the covariance matrix of X. Since Y . Y has zero mean, we have 
var(Y ) = var(Y . Y )= v T 1(X . ┢)T(X . ┢)v. 
n . 1
The above suggests that we should move the origin to the centroid ┢ of the Xi’s and consider the matrix X . ┢ of the centered data points Xi . ┢. 
From now on beware that we denote the columns of X . ┢ by C1,...,Cd and that Y denotes the centered point Y =(X . ┢)v = jd =1 vjCj, where v is a unit vector. 
Basic idea of PCA: The principal components of X are uncorrelated projections Y of the data points X1, ..., Xn onto some directions v (where the v’s are unit vectors) such that var(Y ) is maximal. 
This suggests the following de.nition: 
De.nition 21.2. Given an n ≠ d matrix X of data points X1,...,Xn, if ┢ is the centroid of the Xi’s, then a .rst principal component of X (.rst PC) is a centered point Y1 =(X .┢)v1, the projection of X1,...,Xn onto a direction v1 such that var(Y1) is maximized, where v1 is a unit vector (recall that Y1 =(X . ┢)v1 is a linear combination of the Cj’s, the columns of X . ┢). 
More generally, if Y1,...,Yk are k principal components of X along some unit vectors v1,...,vk, where 1 ● k<d,a(k +1)th principal component of X ((k +1)th PC) is a centered point Yk+1 =(X . ┢)vk+1, the projection of X1,...,Xn onto some direction vk+1 such that var(Yk+1) is maximized, subject to cov(Yh,Yk+1) = 0 for all h with 1 ● h ● k, and where vk+1 is a unit vector (recall that Yh =(X . ┢)vh is a linear combination of the Cj’s). The vh are called principal directions. 
The following proposition is the key to the main result about PCA. This result was already proven in Proposition 16.23 except that the eigenvalues were listed in increasing order. For the reader’s convenience we prove it again. 
Proposition 21.10. If A is a symmetric d ≠ d matrix with eigenvalues ┡1 ◎ ┡2 ◎ ··· ◎ ┡d and if (u1,...,ud) is any orthonormal basis of eigenvectors of A, where ui is a unit eigenvector associated with ┡i, then 
xTAx 
max = ┡1
T
x x
=0 x 
(with the maximum attained for x = u1) and 
xTAx 
max = ┡k+1
T
x x
=0,xÅ{u1,...,uk}￥ x 
(with the maximum attained for x = uk+1), where 1 ● k ● d . 1. 
Proof. First observe that 
xTAx 
max = max {x TAx | x T x =1},
T
x xx
=0 x 
and similarly, 
xTAx  T  
max = max x TAx | (x Å{u1,...,uk}￥) ∞ (xx = 1). 
x xTx
=0,xÅ{u1,...,uk}￥ x
Since A is a symmetric matrix, its eigenvalues are real and it can be diagonalized with respect to an orthonormal basis of eigenvectors, so let (u1,...,ud) be such a basis. If we write 
d 
x = xiui, 
i=1 

a simple computation shows that 
d x TAx = ┡ixi 2 . i=1 

21.4. PRINCIPAL COMPONENTS ANALYSIS (PCA) 
If xTx = 1, then d x2 = 1, and since we assumed that ┡1 ◎ ┡2 ◎ · ·· ◎ ┡d, we get 
i=1 i 
dd x TAx = ┡ix 2 i ● ┡1 x 2 i = ┡1. i=1 i=1 
Thus, max x TAx | x T x =1 ● ┡1, 
x 
and since this maximum is achieved for e1 = (1, 0,..., 0), we conclude that 
max x TAx | x T x =1 = ┡1. 
x 
Next observe that x Å{u1,...,uk}￥ and xTx = 1 i. x1 = ··· = xk = 0 and d xi = 1.
i=1 
Consequently, for such an x, we have 
dd x TAx = ┡ixi 2 ● ┡k+1 xi 2 = ┡k+1. i=k+1 i=k+1 
Thus, max x TAx | (x Å{u1,...,uk}￥) ∞ (x T x = 1) ● ┡k+1, 
x 
and since this maximum is achieved for ek+1 = (0,..., 0, 1, 0,..., 0) with a 1 in position k +1, we conclude that 
max x TAx | (x Å{u1,...,uk}￥) ∞ (x T x =1) = ┡k+1, 
x 
as claimed. 
The quantity 
xTAx T
xx 
is known as the Rayleigh ratio or Rayleigh쭯Ritz ratio (see Section 16.6 ) and Proposition 
21.10 is often known as part of the Rayleigh쭯Ritz theorem. 
Proposition 21.10 also holds if A is a Hermitian matrix and if we replace xTAx by x .Ax and xTx by x . x. The proof is unchanged, since a Hermitian matrix has real eigenvalues and is diagonalized with respect to an orthonormal basis of eigenvectors (with respect to the Hermitian inner product). 
We then have the following fundamental result showing how the SVD of X yields the PCs: 
Theorem 21.11. (SVD yields PCA) Let X be an n ≠ d matrix of data points X1,...,Xn, and let ┢ be the centroid of the Xi’s. If X . ┢ = V DUT is an SVD decomposition of X . ┢ and if the main diagonal of D consists of the singular values ┮1 ◎ ┮2 ◎ · ·· ◎ ┮d, then the centered points Y1,...,Yd, where 
Yk =(X . ┢)uk = kth column of VD 
and uk is the kth column of U, are d principal components of X. Furthermore, 
┮2 var(Yk)= k 
n . 1 and cov(Yh,Yk)=0, whenever h = k and 1 ● k, h ● d. Proof. Recall that for any unit vector v, the centered projection of the points X1,...,Xn onto the line of direction v is Y =(X . ┢)v and that the variance of Y is given by 
var(Y )= v T 1(X . ┢)T(X . ┢)v. 
n . 1
Since X . ┢ = V DUT, we get  
var(Y )  =  v T 1 (n . 1)(X . ┢)T(X . ┢)v  
=  v T  1  UDV TV DUT v  
(n . 1)  
=  v TU  1  D2UT v.  
(n . 1)  

Similarly, if Y =(X . ┢)v and Z =(X . ┢)w, then the covariance of Y and Z is given by 1
TUD2UT
cov(Y, Z)= v w. 
(n . 1) 
2
21
┮
1 D2UT is a symmetric matrix whose eigenvalues are ┮
◎ · ·· ◎ 

d
Obviously, U 

, and 

(n.1) 
n.1 n.1 
the columns of U form an orthonormal basis of unit eigenvectors. We proceed by induction on k. For the base case, k = 1, maximizing var(Y ) is equivalent to maximizing 
1
TUD2UT 
v v, 
(n . 1) where v is a unit vector. By Proposition 21.10, the maximum of the above quantity is the 
21
largest eigenvalue of U 1 D2UT, namely ┮
(n.1) 
, and it is achieved for u1, the .rst columnn 
n.1 
of U. Now we get Y1 =(X . ┢)u1 = V DUT u1, 
and since the columns of U form an orthonormal basis, UTu1 = e1 = (1, 0,..., 0), and so Y1 is indeed the .rst column of VD. 

21.4. PRINCIPAL COMPONENTS ANALYSIS (PCA) 
By the induction hypothesis, the centered points Y1,...,Yk, where Yh =(X . ┢)uh and u1,...,uk are the .rst k columns of U, are k principal components of X. Because 1
TUD2UT
cov(Y, Z)= v w, 
(n . 1) 
where Y =(X . ┢)v and Z =(X . ┢)w, the condition cov(Yh,Z) = 0 for h =1,...,k is equivalent to the fact that w belongs to the orthogonal complement of the subspace spanned by {u1,...,uk}, and maximizing var(Z) subject to cov(Yh,Z)=0 for h =1,...,k is equivalent to maximizing 
1
TUD2UT 
w w, 
(n . 1) where w is a unit vector orthogonal to the subspace spanned by {u1,...,uk}. By Proposition 21.10, the maximum of the above quantity is the (k+1)th eigenvalue of U 1 D2UT, namely 
(n.1) 
┮2 
k+1 
n.1 , and it is achieved for uk+1, the (k + 1)th columnn of U. Now we get 
Yk+1 =(X . ┢)uk+1 = V DUT uk+1, 
and since the columns of U form an orthonormal basis, UTuk+1 = ek+1, and Yk+1 is indeed the (k + 1)th column of VD, which completes the proof of the induction step. 
The d columns u1,...,ud of U are usually called the principal directions of X . ┢ (and X). We note that not only do we have cov(Yh,Yk) = 0 whenever h = k, but the directions u1,...,ud along which the data are projected are mutually orthogonal. 
Example 21.10. For the centered data set of our bearded mathematicians (Example 21.9) we have X . ┢ = V ┣UT, where ┣ has two nonzero singular values, ┮1 = 116.9803,┮2 = 21.7812, and with 
0.9995 0.0325 
U = ,
0.0325 .0.9995 so the principal directions are u1 = (0.9995, 0.0325) and u2 = (0.0325, .0.9995). Observe that u1 is almost the direction of the x-axis, and u2 is almost the opposite direction of the y-axis. We also .nd that the projections Y1 and Y2 along the principal directions are 
.
.
.
. 
.51.5550 3.9249 .51.4000 .5.6000 

VD = 

............... 

9.8031 .6.0843 

.76.5417 

3.1116 

.2.0929 .9.4731 

33.4651 4.6912 

25.5669 1.4325 

53.3894 7.3408 

13.2107 6.0330 

6.4794 3.8128 

............... 

, with 

X . ┢ = 

............... 

............... 

. 

9.6000 6.4000 

.76.4000 .5.6000 

.2.4000 9.4000 

33.6000 

.3.6000 

25.6000 

.0.6000 

53.6000 

.5.6000 

13.4000 

.5.6000 

6.6000 

.3.6000 

15.0607 .13.9174 14.6000 14.4000 See Figures 21.4, 21.5, and 21.6. 


21.5. BEST AFFINE APPROXIMATION 

We know from our study of SVD that ┮12,...,┮d 2 are the eigenvalues of the symmetric positive semide.nite matrix (X . ┢)T(X . ┢) and that u1,...,ud are corresponding eigen-vectors. Numerically, it is preferable to use SVD on X . ┢ rather than to compute explicitly (X . ┢)T(X . ┢) and then diagonalize it. Indeed, the explicit computation of ATA from a matrix A can be numerically quite unstable, and good SVD algorithms avoid computing ATA explicitly. 
In general, since an SVD of X is not unique, the principal directions u1,...,ud are not unique. This can happen when a data set has some rotational symmetries, and in such a case, PCA is not a very good method for analyzing the data set. 
21.5 Best A.ne Approximation 
A problem very close to PCA (and based on least squares) is to best approximate a data set of n points X1,...,Xn, with Xi Å Rd, by a p-dimensional a.ne subspace A of Rd, with 1 ● p ● d . 1 (the terminology rank d . p is also used). 
First consider p = d . 1. Then A = A1 is an a.ne hyperplane (in Rd), and it is given by an equation of the form 
a1x1 + ··· + adxd + c =0. 
By best approximation, we mean that (a1,...,ad,c) solves the homogeneous linear system 
....
.
.
a1 0 
··· 1
x11 x1 d 
.. 

.... 

. 

. 

. 

ad 
.... 

= 

.... 

. 

. 

. 

0 

.... 

..

. . .. 

. . .. 

. . .. 

··· 1
xn 1 xnd 
c 0 

in the least squares sense, subject to the condition that a =(a1,...,ad) is a unit vector, that is, aTa = 1, where Xi =(xi 1, ··· ,xid). 
If we form the symmetric matrix 
.
.
T .
. 

··· 1 ··· 1
x11 x1 d x11 x1 d 
..

. . .. 

. . .. 

. . .. 

.. 
..

. . .. 

. . .. 

. . .. 

.. 

xn 1 ··· xnd 1 xn 1 ··· xnd 1 
involved in the normal equations, we see that the bottom row (and last column) of that matrix is 
n┢1 ··· n┢d n, 
where n┢j = in =1 xij is n times the mean of the column Cj of X. Therefore, if (a1,...,ad,c) is a least squares solution, that is, a solution of the normal equations, we must have n┢1a1 + ··· + n┢dad + nc =0, 
that is, a1┢1 + ··· + ad┢d + c =0, 
which means that the hyperplane A1 must pass through the centroid ┢ of the data points X1,...,Xn. Then we can rewrite the original system with respect to the centered data Xi . ┢, .nd that the variable c drops out, get the system 
(X . ┢)a =0, 
where a =(a1,...,ad). 
Thus, we are looking for a unit vector a solving (X . ┢)a = 0 in the least squares sense, that is, some a such that aTa = 1 minimizing 
a T(X . ┢)T(X . ┢)a. 
Compute some SVD V DUT of X . ┢, where the main diagonal of D consists of the singular values ┮1 ◎ ┮2 ◎ · ·· ◎ ┮d of X . ┢ arranged in descending order. Then 
a T(X . ┢)T(X . ┢)a = a TUD2UT a, 
where D2 = diag(┮12,...,┮d2) is a diagonal matrix, so pick a to be the last column in U (corresponding to the smallest eigenvalue ┮d 2 of (X . ┢)T(X . ┢)). This is a solution to our best .t problem. 

21.5. BEST AFFINE APPROXIMATION 
Therefore, if Ud.1 is the linear hyperplane de.ned by a, that is, 
Ud.1 = {u Å Rd|（ u, a령 =0}, 
where a is the last column in U for some SVD V DUT of X . ┢, we have shown that the a.ne hyperplane A1 = ┢ + Ud.1 is a best approximation of the data set X1,...,Xn in the least squares sense. 
It is easy to show that this hyperplane A1 = ┢ + Ud.1 minimizes the sum of the square distances of each Xi to its orthogonal projection onto A1. Also, since Ud.1 is the orthogonal complement of a, the last column of U, we see that Ud.1 is spanned by the .rst d.1 columns of U, that is, the .rst d . 1 principal directions of X . ┢. 
All this can be generalized to a best (d.k)-dimensional a.ne subspace Ak approximating X1,...,Xn in the least squares sense (1 ● k ● d . 1). Such an a.ne subspace Ak is cut out by k independent hyperplanes Hi (with 1 ● i ● k), each given by some equation 
ai 1x1 + ··· + aidxd + ci =0. 
If we write ai =(ai 1, ··· ,aid), to say that the Hi are independent means that a1,...,ak are linearly independent. In fact, we may assume that a1,...,ak form an orthonormal system. 
Then .nding a best (d . k)-dimensional a.ne subspace Ak amounts to solving the ho-mogeneous linear system 
.
. 

a1
.
.
.
.

...... 

...... 

= 

X 1 0 ··· 000 

0

c1 
. 

. 

. 

ak ck 
.. 

.. 

....
. 
. 
. ... ...
.
. .. 

.... 

,

.

. .. 

... 

0 00 ··· 0 X 1 

0 

T
in the least squares sense, subject to the conditions ai aj = ┙ij, for all i, j with 1 ● i, j ● k, where the matrix of the system is a block diagonal matrix consisting of k diagonal blocks (X, 1), where 1 denotes the column vector (1,..., 1) Å Rn . 
Again it is easy to see that each hyperplane Hi must pass through the centroid ┢ of X1,...,Xn, and by switching to the centered data Xi . ┢ we get the system 
..
.
..
. 
X . ┢ 0 ··· 0 a1 0 
.. 

.. 
..

. 

. 

. 

.... 
..
. 
=. 
. 
, 

.. .
.
.. 

..

.

.. 

. 

00 ··· X . ┢ak 0 
with aT i aj = ┙ij for all i, j with 1 ● i, j ● k. 
If V DUT = X .┢ is an SVD decomposition, it is easy to see that a least squares solution of this system is given by the last k columns of U, assuming that the main diagonal of D consists of the singular values ┮1 ◎ ┮2 ◎· ··◎ ┮d of X .┢ arranged in descending order. But now the (d . k)-dimensional subspace Ud.k cut out by the hyperplanes de.ned by a1,...,ak 
is simply the orthogonal complement of (a1,...,ak), which is the subspace spanned by the .rst d . k columns of U. 
So the best (d .k)-dimensional a.ne subpsace Ak approximating X1,...,Xn in the least squares sense is 
Ak = ┢ + Ud.k, 
where Ud.k is the linear subspace spanned by the .rst d.k principal directions of X .┢, that is, the .rst d . k columns of U. Consequently, we get the following interesting interpretation of PCA (actually, principal directions): 
Theorem 21.12. Let X be an n ≠ d matrix of data points X1,...,Xn, and let ┢ be the centroid of the Xi’s. If X . ┢ = V DUT is an SVD decomposition of X . ┢ and if the main diagonal of D consists of the singular values ┮1 ◎ ┮2 ◎ · ·· ◎ ┮d, then a best (d . k)-dimensional a.ne approximation Ak of X1,...,Xn in the least squares sense is given by 
Ak = ┢ + Ud.k, 
where Ud.k is the linear subspace spanned by the .rst d . k columns of U, the .rst d . k principal directions of X . ┢ (1 ● k ● d . 1). 
Example 21.11. Going back to Example 21.10, a best 1-dimensional a.ne approximation A1 is the a.ne line passing through (┢1,┢2) = (1824.4, 5.6) of direction u1 = (0.9995, 0.0325). 
Example 21.12. Suppose in the data set of Example 21.5 that we add the month of birth of every mathematician as a feature. We obtain the following data set. 
Name  month  year  length  
Carl Friedrich Gauss  4  1777  0  
Camille Jordan  1  1838  12  
Adrien-Marie Legendre  9  1752  0  
Bernhard Riemann  9  1826  15  
David Hilbert  1  1862  2  
Henri Poincar▽e  4  1854  5  
Emmy Noether  3  1882  0  
Karl Weierstrass  10  1815  0  
Eugenio Beltrami  10  1835  2  
Hermann Schwarz  1  1843  20  

The mean of the .rst column is 5.2, and the centered data set is given below. 

21.5. BEST AFFINE APPROXIMATION 
Name  month  year  length  
Carl Friedrich Gauss  .1.2  .51.4  .5.6  
Camille Jordan  .4.2  9.6  6.4  
Adrien-Marie Legendre  3.8  .76.4  .5.6  
Bernhard Riemann  3.8  .2.4  9.4  
David Hilbert  .4.2  33.6  .3.6  
Henri Poincar▽e  .1.2  25.6  .0.6  
Emmy Noether  .2.2  53.6  .5.6  
Karl Weierstrass  4.8  13.4  .5.6  
Eugenio Beltrami  4.8  6.6  .3.6  
Hermann Schwarz  .4.2  14.6  14.4  

Running SVD on this data set we get 

.
. 
, 

D = 

............... 

117.07060 0 0 22.0390 0 0 010.1571 00 
0 
0  0  0  
0  0  0  
0  0  0  
0  0  0  
0  0  0  
0  0  0  

............... 

.
. 
0.0394 0.1717 0.9844 

.

.

U 

= 

.0.9987 

0.0390 0.0332 

, 

.0.0327 .0.9844 0.1730 

and 
.
. 
............... 

51.4683 3.3013 .3.8569 .9.9623 .6.6467 .2.7082 76.6327 3.1845 0.2348 2.2393 .8.6943 5.2872 .33.6038 4.1334 .3.6415 
.25.5941 1.3833 .0.4350 .53.4333 7.2258 .1.3547 .13.0100 6.8594 4.2010 .6.2843 4.6254 4.3212 .15.2173 .14.3266 .1.1581 
............... 

VD = 

, 

.
. 
............... 

.1.2000 .51.4000 .5.6000 .4.2000 9.6000 6.4000 3.8000 .76.4000 .5.6000 3.8000 .2.4000 9.4000 .4.2000 33.6000 .3.6000 
.1.2000 25.6000 .0.6000 .2.2000 53.6000 .5.6000 4.8000 13.4000 .5.6000 4.8000 6.6000 .3.6000 
.4.2000 14.6000 14.4000 
............... 

X . ┢ = 

. 

The .rst principal direction u1 = (0.0394, .0.9987, .0.0327) is basically the opposite of the y-axis, and the most signi.cant feature is the year of birth. The second principal direction u2 = (0.1717, 0.0390, .0.9844) is close to the opposite of the z-axis, and the second most signi.cant feature is the lenght of beards. A best a.ne plane is spanned by the vectors u1 and u2. 
There are many applications of PCA to data compression, dimension reduction, and pattern analysis. The basic idea is that in many cases, given a data set X1,...,Xn, with Xi Å Rd, only a “small” subset of m<d of the features is needed to describe the data set accurately. 
If u1,...,ud are the principal directions of X . ┢, then the .rst m projections of the data (the .rst m principal components, i.e., the .rst m columns of VD) onto the .rst m principal directions represent the data without much loss of information. Thus, instead of using the original data points X1,...,Xn, with Xi Å Rd, we can use their projections onto the .rst m principal directions Y1,...,Ym, where Yi Å Rm and m<d, obtaining a compressed version of the original data set. 
For example, PCA is used in computer vision for face recognition. Sirovitch and Kirby (1987) seem to be the .rst to have had the idea of using PCA to compress facial images. They introduced the term eigenpicture to refer to the principal directions, ui. However, an explicit face recognition algorithm was given only later by Turk and Pentland (1991). They renamed eigenpictures as eigenfaces. 

21.6. SUMMARY 
For details on the topic of eigenfaces, see Forsyth and Ponce [21] (Chapter 22, Section 22.3.2), where you will also .nd exact references to Turk and Pentland’s papers. 
Another interesting application of PCA is to the recognition of handwritten digits. Such an application is described in Hastie, Tibshirani, and Friedman, [34] (Chapter 14, Section 14.5.1). 


21.6 Summary 
The main concepts and results of this chapter are listed below: 
. 	
Least squares problems. 

. 	
Existence of a least squares solution of smallest norm (Theorem 21.1). 

. 	
The pseudo-inverse A+ of a matrix A. 

. 	
The least squares solution of smallest norm is given by the pseudo-inverse (Theorem 21.2) 

. 	
Projection properties of the pseudo-inverse. 

. 	
The pseudo-inverse of a normal matrix. 

. 	
The Penrose characterization of the pseudo-inverse. 

. 	
Data compression and SVD. 

. 	
Best approximation of rank <r of a matrix. 

. 	
Principal component analysis. 

. 	
Review of basic statistical concepts: mean, variance, covariance, covariance matrix. 

. 	
Centered data, centroid. 

. 	
The principal components (PCA). 

. 	
The Rayleigh쭯Ritz theorem (Theorem 21.10). 

. 	
The main theorem: SVD yields PCA (Theorem 21.11). 

. 	
Best a.ne approximation. 

. 	
SVD yields a best a.ne approximation (Theorem 21.12). 

. 	
Face recognition, eigenfaces. 



21.7 Problems 
Problem 21.1. Consider the overdetermined system in the single variable x: a1x = b1,...,amx = bm. Prove that the least squares solution of smallest norm is given by a1b1 + ··· + ambm
+ 
x = . 
a12 + ··· + a2 
m 
Problem 21.2. Let X be an m ≠ n real matrix. For any strictly positive constant K> 0, the matrix XTX + KIn is invertible. Prove that the limit of the matrix (XTX + KIn).1XT when K goes to zero is equal to the pseudo-inverse X+ of X. 
Problem 21.3. Use Matlab to .nd the pseudo-inverse of the 8 ≠ 6 matrix 
.
. 
A = 

........... 

64  2  3  61  60  6  
9  55  54  12  13  51  
17  47  46  20  21  43  
40  26  27  37  36  30  

32  34  35  29  28  38  
41  23  22  44  45  19  
49  15  14  52  53  11  
8  58  59  5  4  62  

........... 

. 

Observe that the sums of the columns are all equal to to 256. Let b be the vector of dimension 6 whose coordinates are all equal to 256. Find the solution x+ of the system Ax = b. 
Problem 21.4. The purpose of this problem is to show that Proposition 21.9 (the Eckart쭯 Young theorem) also holds for the Frobenius norm. This problem is adapted from Strang [65], Section I.9. 
Suppose the m ≠ n matrix B of rank at most k minimizes lA . BlF . Start with an SVD 
of B,  
B = V  D 0  0 0  UT ,  
where D is a diagonal k ≠ k matrix. We can write  
A = V  L + E + R G  F H  UT ,  

where L is strictly lower triangular in the .rst k rows, E is diagonal, and R is strictly upper triangular, and let 
L + D + RF UT
C = V,
00 
21.7. PROBLEMS 
which clearly has rank ● k. 
(1) Prove that 
lA . Bl2 = lA . Cl2 + lLl2 + lRl2 + lF l2 F .
F FFF 
Since lA . BlF is minimal, show that L = R = F = 0. Similarly, show that G = 0. 
(2) We have E 0 D 0 
V TAU = ,V TBU = ,
0 H 00 where E is diagonal, so deduce that 
1. 
D = diag(┮1,...,┮k). 

2. 
The singular values of H must be the smallest n . k singular values of A. 

3. 
The minimum of lA . Blmust be lHl=(┮2 ·· + ┮2)1/2 .


F Fk+1 + · r 
Problem 21.5. Prove that the closest rank 1 approximation (in ll2) of the matrix 
30 
A = 
45 
is 
3 11 
A1 = . 
2 33 
Show that the Eckart쭯Young theorem fails for the operator norm ll◇ by .nding a rank 1 matrix B such that lA . Bl< lA . A1l◇.
◇ 
Problem 21.6. Find a closest rank 1 approximation (in ll2) for the matrices 
.  .  
A = . 3 0 0  0 2 0  0 0 1 . ,  A =  0 2  3 0  ,  A =  2 1  1 2  .  

Problem 21.7. Find a closest rank 1 approximation (in ll2) for the matrix cos ┍ . sin ┍ 
A = . 
sin ┍ cos ┍ Problem 21.8. Let S be a real symmetric positive de.nite matrix and let S = U┣UT be a 
T
diagonalization of S. Prove that the closest rank 1 matrix (in the L2-norm) to S is u1┮1u1 , where u1 is the .rst column of U. 



