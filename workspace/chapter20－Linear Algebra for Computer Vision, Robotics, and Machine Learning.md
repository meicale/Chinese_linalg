Chapter 20 
Singular Value Decomposition and Polar Form 
20.1 Properties of f.. f 
In this section we assume that we are dealing with real Euclidean spaces. Let f : E ↙ E be any linear map. In general, it may not be possible to diagonalize f. We show that every linear map can be diagonalized if we are willing to use two orthonormal bases. This is the celebrated singular value decomposition (SVD). A close cousin of the SVD is the polar form of a linear map, which shows how a linear map can be decomposed into its purely rotational component (perhaps with a .ip) and its purely stretching part. 
The key observation is that f. . f is self-adjoint since 
((f . . f)(u),vㄘ = (f(u),f(v)ㄘ = (u, (f . . f)(v)ㄘ. 
Similarly, f . f. is self-adjoint. 
The fact that f. . f and f . f. are self-adjoint is very important, because by Theorem 16.8, it implies that f. .f and f .f. can be diagonalized and that they have real eigenvalues. In fact, these eigenvalues are all nonnegative as shown in the following proposition. 
Proposition 20.1. The eigenvalues of f. . f and f . f. are nonnegative. 
Proof. If u is an eigenvector of f. . f for the eigenvalue 竹, then 
((f . . f)(u),uㄘ = (f(u),f(u)ㄘ 
and ((f . . f)(u),uㄘ = 竹(u, uㄘ, 
and thus 竹(u, uㄘ = (f(u),f(u)ㄘ, 
which implies that 竹 ≡ 0, since (.,.ㄘ is positive de.nite. A similar proof applies to f . f. . 
641 

Thus, the eigenvalues of f. .f are of the form 考12,...,考r 2 or 0, where 考i > 0, and similarly for f . f. . 
The above considerations also apply to any linear map f : E ↙ F between two Euclidean spaces (E, (.,.ㄘ1) and (F, (.,.ㄘ2). Recall that the adjoint f. : F ↙ E of f is the unique linear map f. such that 
(f(u),vㄘ2 = (u, f . (v)ㄘ1, for all u ﹋ E and all v ﹋ F. 
Then f. . f and f . f. are self-adjoint (the proof is the same as in the previous case), and the eigenvalues of f. . f and f . f. are nonnegative. 
Proof. If 竹 is an eigenvalue of f. . f and u ( = 0) is a corresponding eigenvector, we have 
((f . . f)(u),uㄘ1 = (f(u),f(u)ㄘ2, 
and also ((f . . f)(u),uㄘ1 = 竹(u, uㄘ1, 
so 竹(u, uㄘ1, = (f(u),f(u)ㄘ2, 
which implies that 竹 ≡ 0. A similar proof applies to f . f. . 
The situation is even better, since we will show shortly that f. . f and f . f. have the same nonzero eigenvalues. 
Remark: Given any two linear maps f : E ↙ F and g : F ↙ E, where dim(E)= n and dim(F )= m, it can be shown that 
竹m det(竹In . g . f)= 竹n det(竹Im . f . g), 
and thus g . f and f . g always have the same nonzero eigenvalues; see Problem 14.14. 
De.nition 20.1. Given any linear map f : E ↙ F , the square roots 考i > 0 of the positive eigenvalues of f. . f (and f . f.) are called the singular values of f. 
De.nition 20.2. A self-adjoint linear map f : E ↙ E whose eigenvalues are nonnegative is called positive semide.nite (or positive), and if f is also invertible, f is said to be positive de.nite. In the latter case, every eigenvalue of f is strictly positive. 
If f : E ↙ F is any linear map, we just showed that f. . f and f . f. are positive semide.nite self-adjoint linear maps. This fact has the remarkable consequence that every linear map has two important decompositions: 
1. The polar form. 

2. The singular value decomposition (SVD). 


20.1. PROPERTIES OF F . . F 
The wonderful thing about the singular value decomposition is that there exist two or-thonormal bases (u1,...,un) and (v1,...,vm) such that, with respect to these bases, f is a diagonal matrix consisting of the singular values of f or 0. Thus, in some sense, f can always be diagonalized with respect to two orthonormal bases. The SVD is also a useful tool for solving overdetermined linear systems in the least squares sense and for data analysis, as we show later on. 
First we show some useful relationships between the kernels and the images of f, f. , f. . f, and f . f. . Recall that if f : E ↙ F is a linear map, the image Im f of f is the subspace f(E) of F , and the rank of f is the dimension dim(Im f) of its image. Also recall that (Theorem 5.8) 
dim (Ker f) + dim (Im f) = dim(E), and that (Propositions 11.11 and 13.13) for every subspace W of E, 
dim (W ) + dim(W ﹠) = dim(E). 
Proposition 20.2. Given any two Euclidean spaces E and F , where E has dimension n and F has dimension m, for any linear map f : E ↙ F , we have 
Ker f = Ker (f . . f), Ker f . = Ker (f . f . ), Ker f = (Im f . )﹠ , Ker f . = (Im f)﹠ , dim(Im f) = dim(Im f . ), 
and f, f. , f. . f, and f . f. have the same rank. 
Proof. To simplify the notation, we will denote the inner products on E and F by the same symbol (.,.ㄘ (to avoid subscripts). If f(u)=0, then (f. . f)(u)= f.(f(u)) = f.(0) = 0, and so Ker f . Ker (f. . f). By de.nition of f., we have 
(f(u),f(u)ㄘ = ((f . . f)(u),uㄘ 
for all u ﹋ E. If(f. . f)(u) = 0, since (.,.ㄘ is positive de.nite, we must have f(u) = 0, and so Ker(f. . f) . Ker f. Therefore, 
Ker f = Ker (f . . f). The proof that Ker f. = Ker (f . f.) is similar. By de.nition of f., we have (f(u),vㄘ = (u, f . (v)ㄘ for all u ﹋ E and all v ﹋ F. (.) This immediately implies that 
Ker f = (Im f . )﹠ and Ker f . = (Im f)﹠ . 
Let us explain why Ker f = (Im f.)﹠, the proof of the other equation being similar. Because the inner product is positive de.nite, for every u ﹋ E, we have 
. u ﹋ Ker f 

. i. f(u)=0 


. i. (f(u),vㄘ = 0 for all v, 

. by (.) i. (u, f.(v)ㄘ = 0 for all v, 

. i. u ﹋ (Im f.)﹠ . 
Since 



dim(Im f)= n . dim(Ker f) and 
dim(Im f . )= n . dim((Im f . )﹠), from 
Ker f = (Im f . )﹠ we also have 
dim(Ker f) = dim((Im f . )﹠), from which we obtain 
dim(Im f) = dim(Im f . ). Since 
dim(Ker (f . . f)) + dim(Im (f . . f)) = dim(E), Ker (f. . f) = Ker f and Ker f = (Im f.)﹠, we get dim((Im f . )﹠) + dim(Im (f . . f)) = dim(E). Since 
dim((Im f . )﹠) + dim(Im f . ) = dim(E), we deduce that 
dim(Im f . ) = dim(Im (f . . f)). A similar proof shows that 
dim(Im f) = dim(Im (f . f . )). Consequently, f, f. , f. . f, and f . f. have the same rank. 

20.2. SINGULAR VALUE DECOMPOSITION FOR SQUARE MATRICES 


20.2 	Singular Value Decomposition for Square Matrices 
We will now prove that every square matrix has an SVD. Stronger results can be obtained if we .rst consider the polar form and then derive the SVD from it (there are uniqueness properties of the polar decomposition). For our purposes, uniqueness results are not as important so we content ourselves with existence results, whose proofs are simpler. Readers interested in a more general treatment are referred to Gallier [25]. 
The early history of the singular value decomposition is described in a fascinating paper by Stewart [61]. The SVD is due to Beltrami and Camille Jordan independently (1873, 1874). Gauss is the grandfather of all this, for his work on least squares (1809, 1823) (but Legendre also published a paper on least squares!). Then come Sylvester, Schmidt, and Hermann Weyl. Sylvester＊s work was apparently ※opaque.§ He gave a computational method to .nd an SVD. Schmidt＊s work really has to do with integral equations and symmetric and asymmetric kernels (1907). Weyl＊s work has to do with perturbation theory (1912). Autonne came up with the polar decomposition (1902, 1915). Eckart and Young extended SVD to rectangular matrices (1936, 1939). 
Theorem 20.3. (Singular value decomposition) For every real n ℅ n matrix A there are two orthogonal matrices U and V and a diagonal matrix D such that A = V DUT, where D is of 
the form  ..  
考1  . . .  
D = ..... . .  考2 . . .  . . . ...  . . . ....,  
. . .  考n  

where 考1,...,考r are the singular values of f, i.e., the (positive) square roots of the nonzero eigenvalues of ATA and AAT, and 考r+1 = ﹞﹞﹞ = 考n =0. The columns of U are eigenvectors of ATA, and the columns of V are eigenvectors of AAT. 
Proof. Since ATA is a symmetric matrix, in fact, a positive semide.nite matrix, there exists an orthogonal matrix U such that 
ATA = UD2UT, with D = diag(考1,...,考r, 0,..., 0), where 考12,...,考r 2 are the nonzero eigenvalues of ATA, and where r is the rank of A; that is, 考1,...,考r are the singular values of A. It follows that 
UTATAU =(AU)TAU = D2 , and if we let fj be the jth column of AU for j =1,...,n, then we have (fi,fjㄘ = 考i 2汛ij, 1 ≒ i, j ≒ r 
and fj =0,r +1 ≒ j ≒ n. 
If we de.ne (v1,...,vr) by vj = 考j .1fj, 1 ≒ j ≒ r, 
then we have 
(vi,vjㄘ = 汛ij, 1 ≒ i, j ≒ r, so complete (v1,...,vr) into an orthonormal basis (v1,...,vr,vr+1,...,vn) (for example, using Gram每Schmidt). Now since fj = 考jvj for j =1 ...,r, we have (vi,fjㄘ = 考j(vi,vjㄘ = 考j汛i,j, 1 ≒ i ≒ n, 1 ≒ j ≒ r and since fj = 0 for j = r +1,...,n, (vi,fjㄘ =0 1 ≒ i ≒ n, r +1 ≒ j ≒ n. If V is the matrix whose columns are v1,...,vn, then V is orthogonal and the above equations prove that V TAU = D, which yields A = V DUT, as required. The equation A = V DUT implies that ATA = UD2UT , AAT = VD2V T , which shows that ATA and AAT have the same eigenvalues, that the columns of U are eigenvectors of ATA, and that the columns of V are eigenvectors of AAT . Example 20.1. Here is a simple example of how to use the proof of Theorem 20.3 to obtain 
      
11 1011 
an SVD decomposition. Let A =. Then AT =, ATA =, and 
001011
  
AAT =2 0 . A simple calculation shows that the eigenvalues of ATA are 2 and 0, and 
00
 ﹟  
1/ 2 
for the eigenvalue 2, a unit eigenvector is﹟ , while a unit eigenvector for the eigenvalue 
1/ 2
 ﹟  
﹟
1/ 2 
0 is﹟ . Observe that the singular values are 考1 = 2 and 考2 = 0. Furthermore, 
.1/ 2
 ﹟﹟  
1/ 21/ 2 
U =﹟﹟ = UT . To determine V , the proof of Theorem 20.3 tells us to .rst 
1/ 2 .1/ 2calculate 
 ﹟  
20 
AU =,
00and then set 
 ﹟    
﹟ 
21 
v1 = (1/ 2)=. 
00
20.2. SINGULAR VALUE DECOMPOSITION FOR SQUARE MATRICES 
Once v1 is determined, since 考2 = 0, we have the freedom to choose v2 such that (v1,v2) 0 10 
forms an orthonormal basis for R2 . Naturally, we chose v2 = and set V = .Of 
1 01 course we could have found V by directly computing the eigenvalues and eigenvectors for AAT . We leave it to the reader to check that 
﹟ 
20 
UT
A = V. 
00 
Theorem 20.3 suggests the following de.nition. 
De.nition 20.3. A triple (U, D, V ) such that A = VDUT, where U and V are orthogonal and D is a diagonal matrix whose entries are nonnegative (it is positive semide.nite) is called a singular value decomposition (SVD) of A. 
The Matlab command for computing an SVD A = V DUT of a matrix A is [V, D, U] = svd(A). 
The proof of Theorem 20.3 shows that there are two orthonormal bases (u1,...,un) and (v1,...,vn), where (u1,...,un) are eigenvectors of ATA and (v1,...,vn) are eigenvectors of AAT . Furthermore, (u1,...,ur) is an orthonormal basis of Im AT,(ur+1,...,un) is an orthonormal basis of Ker A,(v1,...,vr) is an orthonormal basis of Im A, and (vr+1,...,vn) is an orthonormal basis of Ker AT . 
Using a remark made in Chapter 3, if we denote the columns of U by u1,...,un and the columns of V by v1,...,vn, then we can write 
TT
A = VDUT = 考1v1u1 + ﹞﹞﹞ + 考rvrur . 
As a consequence, if r is a lot smaller than n (we write r▲ n), we see that A can be reconstructed from U and V using a much smaller number of elements. This idea will be used to provide ※low-rank§ approximations of a matrix. The idea is to keep only the k top singular values for some suitable k▲ r for which 考k+1,...考r are very small. 
Remarks: 
(1) 
In Strang [64] the matrices U, V, D are denoted by U = Q2, V = Q1, and D = 曳, and an SVD is written as A = Q1曳QT 2 . This has the advantage that Q1 comes before Q2 in A = Q1曳QT 2 . This has the disadvantage that A maps the columns of Q2 (eigenvectors of ATA) to multiples of the columns of Q1 (eigenvectors of AAT). 

(2) 
Algorithms for actually computing the SVD of a matrix are presented in Golub and Van Loan [30], Demmel [16], and Trefethen and Bau [68], where the SVD and its applications are also discussed quite extensively. 

(3) 
If A is a symmetric matrix, then in general, there is no SVD V 曳UT of A with V = U. However, if A is positive semide.nite, then the eigenvalues of A are nonnegative, and so the nonzero eigenvalues of A are equal to the singular values of A and SVDs of A are of the form 

A = V 曳V T . 

(4) 
The SVD also applies to complex matrices. In this case, for every complex n℅n matrix A, there are two unitary matrices U and V and a diagonal matrix D such that 


A = VDU . , 
where D is a diagonal matrix consisting of real entries 考1,...,考n, where 考1,...,考r are the singular values of A, i.e., the positive square roots of the nonzero eigenvalues of A.A and AA., and 考r+1 = ... = 考n = 0. 


20.3 Polar Form for Square Matrices 
A notion closely related to the SVD is the polar form of a matrix. De.nition 20.4. A pair (R, S) such that A = RS with R orthogonal and S symmetric positive semide.nite is called a polar decomposition of A. Theorem 20.3 implies that for every real n ℅ n matrix A, there is some orthogonal matrix R and some positive semide.nite symmetric matrix S such that A = RS. This is easy to show and we will prove it below. Furthermore, R, S are unique if A is invertible, but this is harder to prove; see Problem 20.9. 
For example, the matrix 

.
. 
11 1 1 

A = 

1 

2 

... 

11 .1 .1 

1 .11 .1 

... 

1 .1 .11 
is both orthogonal and symmetric, and A = RS with R = A and S = I, which implies that some of the eigenvalues of A are negative. 
Remark: In the complex case, the polar decomposition states that for every complex n ℅ n matrix A, there is some unitary matrix U and some positive semide.nite Hermitian matrix H such that 
A = UH. 
It is easy to go from the polar form to the SVD, and conversely. 
Given an SVD decomposition A = VDUT, let R = VUT and S = UDUT . It is clear that R is orthogonal and that S is positive semide.nite symmetric, and 
RS = VUTUDUT = VDUT = A. 
20.3. POLAR FORM FOR SQUARE MATRICES 
Example 20.2. Recall from Example 20.1 that A = V DUT where V = I2 and 
 
 
﹟ 

11 

1122 
﹟
﹟
20 

A = ,U 

D = 

. 

00 

=

,

1. 1
22
﹟
﹟
00 

Set R = VUT = U and
 
 
11
S = UDUT =12 12 . 
﹟﹟ ﹟﹟ 
22
﹟ 
A = RS. 
1﹟ 2 
ATA, S has eigenvalues 
2 and 0. We leave it to the reader to check that 

Since S 

= 

T 
Going the other way, given a polar decomposition A = R1S, where R1 is orthogonal and S is positive semide.nite symmetric, there is an orthogonal matrix R2 and a positive 2 , and thus 
semide.nite diagonal matrix D such that S = R2DR
A = R1R2DR
T 
2 
= VDUT 
, 

where V = R1R2 and U = R2 are orthogonal. 
﹟﹟ 
11 1/ 21/ 2 
Example 20.3. Let A = and A = R1S, where R1 = ﹟﹟ and S = 
00 1/ 2 .1/ 2
﹟﹟ 
1/ 21/ 2
﹟﹟ . This is the polar decomposition of Example 20.2. Observe that 
1/ 21/ 2 
 
  
 
﹟ 

﹟
﹟
1122 
﹟
﹟
1122
20

T
S 

= R2DR
=

2 . 
﹟
﹟
﹟
﹟
1. 11. 1
2222
10 
Set U = R2 and V = R1R2 = to obtain the SVD decomposition of Example 20.1. 
01 
The eigenvalues and the singular values of a matrix are typically not related in any obvious way. For example, the n ℅ n matrix 
00

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

.. 

.. 

0  0  . . .  0  1  2  0  
0  0  . . .  0  0  1  2  
0  0  . . .  0  0  0  1  

has the eigenvalue 1 with multiplicity n, but its singular values, 考1 ≡ ﹞﹞﹞ ≡ 考n, which are the positive square roots of the eigenvalues of the matrix B = ATA with 
.
. 
120 0 ... 00 252 0 ... 00 025 2 ... 00 
.. 	..
...
B = 

.......... 

.......... 

.. 

.

.

...

. 

. 

.

.. 

.. 

0  0  . . .  2  5  2  0  
0  0  . . .  0  2  5  2  
0  0  . . .  0  0  2  5  

have a wide spread, since 
考1 
= cond2(A) ≡ 2n.1 . 
考n If A is a complex n ℅ n matrix, the eigenvalues 竹1,...,竹n and the singular values 考1 ≡ 考2 ≡﹞ ﹞﹞≡ 考n of A are not unrelated, since 
考2 考2 
1 ﹞﹞﹞ = det(A . A)= | det(A)|2 
n 
and |竹1|﹞﹞﹞|竹n| = | det(A)|, 
so we have 
|竹1|﹞﹞﹞|竹n| = 考1 ﹞﹞﹞ 考n. 
More generally, Hermann Weyl proved the following remarkable theorem: 
Theorem 20.4. (Weyl＊s inequalities, 1949 ) For any complex n℅n matrix, A, if 竹1,...,竹n ﹋ C are the eigenvalues of A and 考1,...,考n ﹋ R+ are the singular values of A, listed so that |竹1|≡ ﹞ ﹞﹞ ≡ |竹n| and 考1 ≡ ﹞ ﹞﹞ ≡ 考n ≡ 0, then 
|竹1|﹞﹞﹞|竹n| = 考1 ﹞﹞﹞ 考n and 
|竹1|﹞﹞﹞|竹k|≒ 考1 ﹞﹞﹞ 考k, for k =1,...,n . 1. 

A proof of Theorem 20.4 can be found in Horn and Johnson [37], Chapter 3, Section 3.3, where more inequalities relating the eigenvalues and the singular values of a matrix are given. 
Theorem 20.3 can be easily extended to rectangular m ℅ n matrices, as we show in the next section. For various versions of the SVD for rectangular matrices, see Strang [64] Golub and Van Loan [30], Demmel [16], and Trefethen and Bau [68]. 


20.4 	Singular Value Decomposition for Rectangular Matrices 
Here is the generalization of Theorem 20.3 to rectangular matrices. 
20.4. SINGULAR VALUE DECOMPOSITION FOR RECTANGULAR MATRICES 651 
Theorem 20.5. (Singular value decomposition) For every real m ℅ n matrix A, there are two orthogonal matrices U (n ℅ n) and V (m ℅ m) and a diagonal m ℅ n matrix D such that A = VDUT, where D is of the form 
.
. 
考1 ... 
考2 ... 
.. .
.
.. ..
.
.. . 
... 考n 
. 
.
0 . ... 0 .. .
.
.. ..
.
.. . 
............ 
............ 
. .  
考1  . . .  0  . . .  0  
..... . .  考2 . . .  . . . ...  . . .  0 0  . . . . . .  0 0 .... 

D = 

or D = 

, 

... 考m 0 ... 0 
. 
.
0 . ... 0 
where 考1,...,考r are the singular values of f, i.e. the (positive) square roots of the nonzero eigenvalues of ATA and AAT, and 考r+1 = ... = 考p =0, where p = min(m, n). The columns of U are eigenvectors of ATA, and the columns of V are eigenvectors of AAT . 
Proof. As in the proof of Theorem 20.3, since ATA is symmetric positive semide.nite, there exists an n ℅ n orthogonal matrix U such that 
ATA = U曳2UT , 
with 曳 = diag(考1,...,考r, 0,..., 0), where 考12,...,考r 2 are the nonzero eigenvalues of ATA, and where r is the rank of A. Observe that r ≒ min{m, n}, and AU is an m ℅ n matrix. It follows that 
UTATAU =(AU)TAU =曳2 , 
and if we let fj ﹋ Rm be the jth column of AU for j =1,...,n, then we have 
(fi,fjㄘ = 考i 2汛ij, 1 ≒ i, j ≒ r 
and fj =0,r +1 ≒ j ≒ n. 
If we de.ne (v1,...,vr) by vj = 考j .1fj, 1 ≒ j ≒ r, 
then we have 
(vi,vjㄘ = 汛ij, 1 ≒ i, j ≒ r, 
so complete (v1,...,vr) into an orthonormal basis (v1,...,vr,vr+1,...,vm) (for example, using Gram每Schmidt). 
Now since fj = 考jvj for j =1 ...,r, we have 
(vi,fjㄘ = 考j(vi,vjㄘ = 考j汛i,j, 1 ≒ i ≒ m, 1 ≒ j ≒ r 
and since fj = 0 for j = r +1,...,n, we have (vi,fjㄘ =0 1 ≒ i ≒ m, r +1 ≒ j ≒ n. If V is the matrix whose columns are v1,...,vm, then V is an m ℅ m orthogonal matrix and 
if m ≡ n, we let 

. 
. 
D =  曳 0m.n  = . . . . . . . . . . . . 考1 . . . 0 . . .  考2 . . . . . . . . .  . . . . . . ... . . . . . . ...  . . . 考n 0 . . . . . . . . . . . . . . . ,  
0  . . .  . . .  0  
else if n ≡ m, then we let  .  .  

D = . . . . 考1 . . .  考2 . . .  . . . . . . ...  . . .  0 0 0  . . . . . . . . .  0 0 0 . . . . .  
. . .  考m  0  . . .  0  
In either case, the above equations prove that  
V TAU = D,  
which yields A = V DUT, as required. The equation A = V DUT implies that  

 爛 .nr . 爛 .mr 
, 0,..., 0
.  
and AAT = V DDTV T = V diag(考12,...,考r 2 , 0,..., 0)V T , 
which shows that ATA and AAT have the same nonzero eigenvalues, that the columns of U are eigenvectors of ATA, and that the columns of V are eigenvectors of AAT . 
A triple (U, D, V ) such that A = VDUT is called a singular value decomposition (SVD) of A. 
ATA = UDTDUT = Udiag(考12,...,考r 2 
)UT 
. 
. 

. 
. 

11 100 11 
Example 20.4. Let A = 0 0 .Then AT = ATA = , and AAT = 
100 11 
00

. 
. 
. 

. 
200 
0 0 0 . The reader should verify that ATA = U曳2UT where 曳2 = 20 and 
00 
000 

20.4. SINGULAR VALUE DECOMPOSITION FOR RECTANGULAR MATRICES 653 

.﹟ ..﹟ . ..
﹟﹟ 
20 21 1/ 21/ 2 
U = UT = ﹟﹟ . Since AU = . 00. , set v1 = ﹟1. 0 . = .0. ,
1/ 2 .1/ 2 2 
00 00
.. .. 
00 and complete an orthonormal basis for R3 by assigning v2 = .1., and v3 = .0.. Thus 
01
.﹟ . 
20 V = I3, and the reader should verify that A = V DUT, where D = . 00.. 00 
Even though the matrix D is an m ℅ n rectangular matrix, since its only nonzero entries are on the descending diagonal, we still say that D is a diagonal matrix. 
The Matlab command for computing an SVD A = V DUT of a matrix A is also [V, D, U] = svd(A). 
If we view A as the representation of a linear map f : E ↙ F , where dim(E)= n and dim(F )= m, the proof of Theorem 20.5 shows that there are two orthonormal bases (u1,..., un) and (v1,...,vm) for E and F , respectively, where (u1,...,un) are eigenvectors of f. . f and (v1,...,vm) are eigenvectors of f . f. . Furthermore, (u1,...,ur) is an orthonormal basis of Im f.,(ur+1,...,un) is an orthonormal basis of Ker f,(v1,...,vr) is an orthonormal basis of Im f, and (vr+1,...,vm) is an orthonormal basis of Ker f. . 
The SVD of matrices can be used to de.ne the pseudo-inverse of a rectangular matrix; we will do so in Chapter 21. The reader may also consult Strang [64], Demmel [16], Trefethen and Bau [68], and Golub and Van Loan [30]. 
One of the spectral theorems states that a symmetric matrix can be diagonalized by an orthogonal matrix. There are several numerical methods to compute the eigenvalues of a symmetric matrix A. One method consists in tridiagonalizing A, which means that there exists some orthogonal matrix P and some symmetric tridiagonal matrix T such that A = PTP T . In fact, this can be done using Householder transformations; see Theorem 
17.2. It is then possible to compute the eigenvalues of T using a bisection method based on Sturm sequences. One can also use Jacobi＊s method. For details, see Golub and Van Loan [30], Chapter 8, Demmel [16], Trefethen and Bau [68], Lecture 26, Ciarlet [14], and Chapter 
17. Computing the SVD of a matrix A is more involved. Most methods begin by .nding orthogonal matrices U and V and a bidiagonal matrix B such that A = V BUT; see Problem 
12.8 and Problem 20.3. This can also be done using Householder transformations. Observe that BTB is symmetric tridiagonal. Thus, in principle, the previous method to diagonalize a symmetric tridiagonal matrix can be applied. However, it is unwise to compute BTB explicitly, and more subtle methods are used for this last step; the matrix of Problem 20.1 can be used, and see Problem 20.3. Again, see Golub and Van Loan [30], Chapter 8, Demmel [16], and Trefethen and Bau [68], Lecture 31. 
The polar form has applications in continuum mechanics. Indeed, in any deformation it is important to separate stretching from rotation. This is exactly what QS achieves. The orthogonal part Q corresponds to rotation (perhaps with an additional re.ection), and the symmetric matrix S to stretching (or compression). The real eigenvalues 考1,...,考r of S are the stretch factors (or compression factors) (see Marsden and Hughes [47]). The fact that S can be diagonalized by an orthogonal matrix corresponds to a natural choice of axes, the principal axes. 
The SVD has applications to data compression, for instance in image processing. The idea is to retain only singular values whose magnitudes are signi.cant enough. The SVD can also be used to determine the rank of a matrix when other methods such as Gaussian elimination produce very small pivots. One of the main applications of the SVD is the computation of the pseudo-inverse. Pseudo-inverses are the key to the solution of various optimization problems, in particular the method of least squares. This topic is discussed in the next chapter (Chapter 21). Applications of the material of this chapter can be found in Strang [64, 63]; Ciarlet [14]; Golub and Van Loan [30], which contains many other references; Demmel [16]; and Trefethen and Bau [68]. 


20.5 Ky Fan Norms and Schatten Norms 
The singular values of a matrix can be used to de.ne various norms on matrices which have found recent applications in quantum information theory and in spectral graph theory. Following Horn and Johnson [37] (Section 3.4) we can make the following de.nitions: 
De.nition 20.5. For any matrix A ﹋ Mm,n(C), let q = min{m, n}, and if 考1 ≡ ﹞ ﹞﹞ ≡ 考q are the singular values of A, for any k with 1 ≒ k ≒ q, let 
Nk(A)= 考1 + ﹞﹞﹞ + 考k, 
called the Ky Fan k-norm of A. More generally, for any p ≡ 1 and any k with 1 ≒ k ≒ q, let 
)1/p
Nk;p(A)=(考1 p + ﹞﹞﹞ + 考kp , 
called the Ky Fan p-k-norm of A. When k = q, Nq;p is also called the Schatten p-norm. 
Observe that when k = 1, N1(A)= 考1, and the Ky Fan norm N1 is simply the spectral norm from Chapter 8, which is the subordinate matrix norm associated with the Euclidean norm. When k = q, the Ky Fan norm Nq is given by 
Nq(A)= 考1 + ﹞﹞﹞ + 考q = tr((A . A)1/2) 
and is called the trace norm or nuclear norm. When p = 2 and k = q, the Ky Fan Nq;2 norm is given by 
 
Nk;2(A)=(考12 + ﹞﹞﹞ + 考q 2)1/2 =tr(A.A)= lAlF , 
which is the Frobenius norm of A. 
It can be shown that Nk and Nk;p are unitarily invariant norms, and that when m = n, they are matrix norms; see Horn and Johnson [37] (Section 3.4, Corollary 3.4.4 and Problem 3). 
20.6. SUMMARY 


20.6 Summary 
The main concepts and results of this chapter are listed below: 
. 	
For any linear map f : E ↙ E on a Euclidean space E, the maps f. . f and f . f. are self-adjoint and positive semide.nite. 

. 	
The singular values of a linear map. 

. 	
Positive semide.nite and positive de.nite self-adjoint maps. 

. 	
Relationships between Im f, Ker f, Im f., and Ker f. . 

. 	
The singular value decomposition theorem for square matrices (Theorem 20.3). 

. 	
The SVD of matrix. 

. 	
The polar decomposition of a matrix. 

. 	
The Weyl inequalities. 

. 	
The singular value decomposition theorem for m ℅ n matrices (Theorem 20.5). 

. 	
Ky Fan k-norms, Ky Fan p-k-norms, Schatten p-norms. 



20.7 Problems 
Problem 20.1. (1) Let A be a real n℅n matrix and consider the (2n)℅(2n) real symmetric 
matrix  
S =  0 AT  A 0  .  

Suppose that A has rank r. If A = V 曳UT is an SVD for A, with 曳 = diag(考1,...,考n) and 考1 ≡ ﹞ ﹞﹞ ≡ 考r > 0, denoting the columns of U by uk and the columns of V by vk, prove that 
考k is an eigenvalue of S with corresponding eigenvector vk for k =1,...,n, and that .考k 
uk 
is an eigenvalue of S with corresponding eigenvector vk for k =1,...,n.
.uk Hint. We have Auk = 考kvk for k =1,...,n. Show that ATvk = 考kuk for k =1,...,r, and that ATvk = 0 for k = r +1,...,n. Recall that Ker (AT) = Ker(AAT). 
(2) 
Prove that the 2n eigenvectors of S in (1) are pairwise orthogonal. Check that if A has rank r, then S has rank 2r. 

(3) 
Now assume that A is a real m ℅ n matrix and consider the (m + n) ℅ (m + n) real symmetric matrix 


0 A 
S = . 
AT 0 
Suppose that A has rank r. If A = V 曳UT is an SVD for A, prove that 考k is an eigenvalue of S with corresponding eigenvector vk for k =1,...,r, and that .考k is an eigenvalue of 
uk 
S with corresponding eigenvector vk for k =1,...,r.
.uk Find the remaining m + n . 2r eigenvectors of S associated with the eigenvalue 0. 
(4) Prove that these m + n eigenvectors of S are pairwise orthogonal. Problem 20.2. Let A be a real m ℅ n matrix of rank r and let q = min(m, n). 
(1) Consider the (m + n) ℅ (m + n) real symmetric matrix 0 A 
S = 
AT  0  
and prove that  
Im  z.1A  zIm  .A  zIm . z.1AAT  0  
0  In  .AT  zIn  =  .AT  zIn  .  
Use the above equation to prove that  

det(zIm+n . S)= tn.m det(t2Im . AAT). 
(2) Prove that the eigenvalues of S are ㊣考1,..., ㊣考q, with |m . n| additional zeros. Problem 20.3. Let B be a real bidiagonal matrix of the form 
.
. 
a1 b1 0 ﹞﹞﹞ 0 .
.
0 a2 b2 .0 
B 

= 

...... 

...... 

.

..
...
. 

.

.

..

. 

. 

.

. 

. 

0 ﹞﹞﹞ 0 an.1 bn.1 00 ﹞﹞﹞ 0 an 
Let A be the (2n) ℅ (2n) symmetric matrix 0 BT 
A = ,
B 0 and let P be the permutation matrix given by P =[e1,en+1,e2,en+2, ﹞﹞﹞ ,en,e2n]. 
(1) Prove that T = P TAP is a symmetric tridiagonal (2n) ℅ (2n) matrix with zero main 
diagonal of the form 

.
. 
............ 

0 a1 000 0 ﹞﹞﹞ 0 a1 0 b1 00 0 ﹞﹞﹞ 0 0 b1 0 a2 00 ﹞﹞﹞ 0 00 a2 0 b2 0 ﹞﹞﹞ 0 
... .
....
........

.. ..
... . 
0  0  0  ﹞ ﹞ ﹞  an.1  0  bn.1  0  
0  0  0  ﹞ ﹞ ﹞  0  bn.1  0  an  
0  0  0  ﹞ ﹞ ﹞  0  0  an  0  

............ 

T = 

. 

20.7. PROBLEMS 
(2) Prove that if xi is a unit eigenvector for an eigenvalue 竹i of T , then 竹i = ㊣考i where 考i is a singular value of B, and that 
1 ui
Pxi = ﹟ , 2 ㊣vi where the ui are unit eigenvectors of BTB and the vi are unit eigenvectors of BBT . Problem 20.4. Find the SVD of the matrix 
.  .  
0  2  0  
A = .0  0  3..  
0  0  0  

Problem 20.5. Let u, v ﹋ Rn be two nonzero vectors, and let A = uvT be the corresponding rank 1 matrix. Prove that the nonzero singular value of A is lul2 lvl2. 
Problem 20.6. Let A be a n℅n real matrix. Prove that if 考1,...,考n are the singular values of A, then 考13,...,考n 3 are the singular values of AATA. 
Problem 20.7. Let A be a real n ℅ n matrix. 
(1) Prove that the largest singular value 考1 of A is given by 
lAxl2
考1 = sup ,奀x=0 lxl2 
and that this supremum is achieved at x = u1, the .rst column in U in an SVD A = V 曳UT . 
(2) Extend the above result to real m ℅ n matrices. 
Problem 20.8. Let A be a real m ℅ n matrix. Prove that if B is any submatrix of A (by keeping M ≒ m rows and N ≒ n columns of A), then (考1)B ≒ (考1)A (where (考1)A is the largest singular value of A and similarly for (考1)B). 
Problem 20.9. Let A be a real n ℅ n matrix. 
(1) 
Assume A is invertible. Prove that if A = Q1S1 = Q2S2 are two polar decompositions of A, then Q1 = Q2 and S1 = S2. Hint. ATA = S12 = S22, with S1 and S2 symmetric positive de.nite. Then use Problem 16.7. 

(2) 
Now assume that A is singular. Prove that if A = Q1S1 = Q2S2 are two polar decompositions of A, then S1 = S2, but Q1 may not be equal to Q2. 


Problem 20.10. (1) Let A be any invertible (real) n ℅ n matrix. Prove that for every SVD, A = V DUT of A, the product VUT is the same (i.e., if V1DU1 T = V2DU2 T , then V1UT = V2U2 T). What does VUT have to do with the polar form of A? 
(2) Given any invertible (real) n ℅ n matrix, A, prove that there is a unique orthogonal matrix, Q ﹋ O(n), such that lA . QlF is minimal (under the Frobenius norm). In fact, prove that Q = VUT, where A = V DUT is anSVDof A. Moreover, if det(A) > 0, show 
that Q ﹋ SO(n). What can you say if A is singular (i.e., non-invertible)? 
Problem 20.11. (1) Prove that for any n ℅ n matrix A and any orthogonal matrix Q, we have 
max{tr(QA) | Q ﹋ O(n)} = 考1 + ﹞﹞﹞ + 考n, 
where 考1 ≡ ﹞ ﹞﹞ ≡ 考n are the singular values of A. Furthermore, this maximum is achieved by Q = UV T, where A = V 曳UT is any SVD for A. 
(2) By applying the above result with A = ZTX and Q = RT, deduce the following result : For any two .xed n ℅ k matrices X and Z, the minimum of the set 
{lX . ZRl| R ﹋ O(k)}
F 
is achieved by R = VUT for any SVD decomposition V 曳UT = ZTX of ZTX. 
Remark: The problem of .nding an orthogonal matrix R such that ZR comes as close as possible to X is called the orthogonal Procrustes problem; see Strang [65] (Section IV.9) for the history of this problem. 



