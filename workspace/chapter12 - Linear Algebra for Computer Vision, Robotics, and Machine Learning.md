Chapter 12 
QR-Decomposition for Arbitrary Matrices 
12.1 Orthogonal Re.ections 
Hyperplane re.ections are represented by matrices called Householder matrices. These ma-trices play an important role in numerical methods, for instance for solving systems of linear equations, solving least squares problems, for computing eigenvalues, and for transforming a symmetric matrix into a tridiagonal matrix. We prove a simple geometric lemma that imme-diately yields the QR-decomposition of arbitrary matrices in terms of Householder matrices. 
Orthogonal symmetries are a very important example of isometries. First let us review the de.nition of projections, introduced in Section 5.2, just after Proposition 5.5. Given a vector space E, let F and G be subspaces of E that form a direct sum E = F �� G. Since every u ﹋ E can be written uniquely as u = v + w, where v ﹋ F and w ﹋ G, we can de.ne the two projections pF : E ↙ F and pG : E ↙ G such that pF (u)= v and pG(u)= w. In Section 5.2 we used the notation 羽1 and 羽2, but in this section it is more convenient to use pF and pG. 
It is immediately veri.ed that pG and pF are linear maps, and that 
p 2 = pF ,p 2 = pG,pF . pG =0, and pF + pG = id.
FG = pG . pF 
. 
De.nition 12.1. Given a vector space E, for any two subspaces F and G that form a direct sum E = F �� G, the symmetry (or re.ection) with respect to F and parallel to G is the linear map s: E ↙ E de.ned such that 
s(u)=2pF (u) . u, 
for every u ﹋ E. 
429 

Because pF + pG = id, note that we also have 
s(u)= pF (u) . pG(u) 
and s(u)= u . 2pG(u), 
s2 = id, s is the identity on F , and s = .id on G. We now assume that E is a Euclidean space of .nite dimension. 
De.nition 12.2. Let E be a Euclidean space of .nite dimension n. For any two subspaces F and G, if F and G form a direct sum E = F �� G and F and G are orthogonal, i.e., F = G﹠, the orthogonal symmetry (or re.ection) with respect to F and parallel to G is the linear map s: E ↙ E de.ned such that 
s(u)=2pF (u) . u = pF (u) . pG(u), 
for every u ﹋ E. When F is a hyperplane, we call s a hyperplane symmetry with respect to F (or re.ection about F ), and when G is a plane (and thus dim(F )= n . 2), we call s a .ip about F . 
A re.ection about a hyperplane F is shown in Figure 12.1. 

Figure 12.1: A re.ection about the peach hyperplane F . Note that u is purple, pF (u) is blue and pG(u) is red. 
For any two vectors u, v ﹋ E, it is easily veri.ed using the bilinearity of the inner product that請u + v請2.請 u . v請2 = 4(u ﹞ v). (.) 
In particular, if u ﹞ v = 0, then請u + v請 =請u . v請. Then since 
u = pF (u)+ pG(u) 
12.1. ORTHOGONAL REFLECTIONS 
and s(u)= pF (u) . pG(u), 
and since F and G are orthogonal, it follows that 
pF (u) ﹞ pG(v)=0, 
and thus by (.)
請s(u)請 =請pF (u) . pG(u)請 =請pF (u)+ pG(u)請 =請u請, 
so that s is an isometry. 
Using Proposition 11.10, it is possible to .nd an orthonormal basis (e1,...,en) of E consisting of an orthonormal basis of F and an orthonormal basis of G. Assume that F has dimension p, so that G has dimension n . p. With respect to the orthonormal basis (e1,...,en), the symmetry s has a matrix of the form
    
Ip  0  
. 

.In.p
0 
Thus, det(s)=(.1)n.p, and s is a rotation i. n . p is even. In particular, when F is a hyperplane H, we have p = n . 1 and n . p = 1, so that s is an improper orthogonal transformation. When F = {0}, we have s = .id, which is called the symmetry with respect to the origin. The symmetry with respect to the origin is a rotation i. n is even, and an improper orthogonal transformation i. n is odd. When n is odd, since s . s = id and det(s)=(.1)n = .1, we observe that every improper orthogonal transformation f is the composition f =(f . s) . s of the rotation f . s with s, the symmetry with respect to the origin. When G is a plane, p = n . 2, and det(s)=(.1)2 = 1, so that a .ip about F is a rotation. In particular, when n = 3, F is a line, and a .ip about the line F is indeed a rotation of measure 羽 as illustrated by Figure 12.2. 
Remark: Given any two orthogonal subspaces F, G forming a direct sum E = F �� G, let f be the symmetry with respect to F and parallel to G, and let g be the symmetry with respect to G and parallel to F . We leave as an exercise to show that 
f . g = g . f = .id. 
When F = H is a hyperplane, we can give an explicit formula for s(u) in terms of any nonnull vector w orthogonal to H. Indeed, from 
u = pH (u)+ pG(u), 
since pG(u) ﹋ G and G is spanned by w, which is orthogonal to H, we have 
pG(u)= 竹w 

for some 竹 ﹋ R, and we get u ﹞ w = 竹請w請2 , 
and thus 
(u ﹞ w)pG(u)= 請w請2 w. Since s(u)= u . 2pG(u), 
we get 
(u ﹞ w)s(u)= u . 2 請w請2 w. Since the above formula is important, we record it in the following proposition. 
Proposition 12.1. Let E be a .nite-dimensional Euclidean space and let H be a hyperplane in E. For any nonzero vector w orthogonal to H, the hyperplane re.ection s about H is given by 
(u ﹞ w)
s(u)= u . 2 請w請2 w, u ﹋ E. 

Such re.ections are represented by matrices called Householder matrices, which play an important role in numerical matrix analysis (see Kincaid and Cheney [39] or Ciarlet [14]). 
De.nition 12.3. A Householder matrix if a matrix of the form 
WW TWW T H = In . 2 請W請2 = In . 2 W TW, 
where W ﹋ Rn is a nonzero vector. 

12.1. ORTHOGONAL REFLECTIONS 
Householder matrices are symmetric and orthogonal. It is easily checked that over an orthonormal basis (e1,...,en), a hyperplane re.ection about a hyperplane H orthogonal to a nonzero vector w is represented by the matrix 
WW TH = In . 2 請W請2 , where W is the column vector of the coordinates of w over the basis (e1,...,en). Since (u ﹞ w)pG(u)= 請w請2 w, the matrix representing pG is 
WW T 
W TW, 
and since pH + pG = id, the matrix representing pH is WW T 
In . . 
W TW 
These formulae can be used to derive a formula for a rotation of R3, given the direction w of its axis of rotation and given the angle 牟 of rotation. 
The following fact is the key to the proof that every isometry can be decomposed as a product of re.ections. 
Proposition 12.2. Let E be any nontrivial Euclidean space. For any two vectors u, v ﹋ E, if請u請 =請v請, then there is a hyperplane H such that the re.ection s about H maps u to v, and if u 
= v, then this re.ection is unique. See Figure 12.3. 
Proof. If u = v, then any hyperplane containing u does the job. Otherwise, we must have H = {v . u}﹠, and by the above formula, s(u)= u . 2(u ﹞ (v . u))(v . u)= u +2請u請2 . 2u ﹞ v(v . u),請(v . u)請2 請(v . u)請2 and since請(v . u)請2 =請u請2 +請v請2 . 2u ﹞ v and請u請 =請v請, we have請(v . u)請2 =2請u請2 . 2u ﹞ v, and thus, s(u)= v.
 If E is a complex vector space and the inner product is Hermitian, Proposition 12.2 
is false. The problem is that the vector v . u does not work unless the inner product u ﹞v is real! The proposition can be salvaged enough to yield the QR-decomposition in terms of Householder transformations; see Section 13.5. 
We now show that hyperplane re.ections can be used to obtain another proof of the QR-decomposition. 


12.2 QR-Decomposition Using Householder Matrices 
First we state the result geometrically. When translated in terms of Householder matrices, we obtain the fact advertised earlier that every matrix (not necessarily invertible) has a QR-decomposition. 
Proposition 12.3. Let E be a nontrivial Euclidean space of dimension n. For any orthonor-mal basis (e1,..., en) and for any n-tuple of vectors (v1, ..., vn), there is a sequence of n isometries h1,...,hn such that hi is a hyperplane re.ection or the identity, and if (r1,...,rn) are the vectors given by 
rj = hn .﹞ ﹞﹞. h2 . h1(vj), 
then every rj is a linear combination of the vectors (e1,...,ej), 1 ≒ j ≒ n. Equivalently, the matrix R whose columns are the components of the rj over the basis (e1,...,en) is an upper triangular matrix. Furthermore, the hi can be chosen so that the diagonal entries of R are nonnegative. 
Proof. We proceed by induction on n. For n = 1, we have v1 = 竹e1 for some 竹 ﹋ R. If 竹 ≡ 0, we let h1 = id, else if 竹< 0, we let h1 = .id, the re.ection about the origin. 
For n ≡ 2, we .rst have to .nd h1. Let 
r1,1 =請v1請. 

12.2. QR-DECOMPOSITION USING HOUSEHOLDER MATRICES 
If v1 = r1,1e1, we let h1 = id. Otherwise, there is a unique hyperplane re.ection h1 such that 
h1(v1)= r1,1 e1, 
de.ned such that 
h1(u)= u . 2(請 uw﹞ 1w請 12 )w1 for all u ﹋ E, where w1 = r1,1 e1 . v1. The map h1 is the re.ection about the hyperplane H1 orthogonal to the vector w1 = r1,1 e1 . v1. See Figure 12.4. Letting 

r1 = h1(v1)= r1,1 e1, it is obvious that r1 belongs to the subspace spanned by e1, and r1,1 =請v1請 is nonnegative. Next assume that we have found k linear maps h1,...,hk, hyperplane re.ections or the identity, where 1 ≒ k ≒ n . 1, such that if (r1,...,rk) are the vectors given by rj = hk .﹞ ﹞﹞. h2 . h1(vj), then every rj is a linear combination of the vectors (e1,...,ej), 1 ≒ j ≒ k. See Figure 
12.5. The vectors (e1,...,ek) form a basis for the subspace denoted by Uk㏒, the vectors (ek+1,...,en) form a basis for the subspace denoted by U㏒㏒, the subspaces U㏒and U㏒㏒are
k kk 
= U㏒�� U㏒㏒
orthogonal, and E kk . Let uk+1 = hk .﹞ ﹞﹞. h2 . h1(vk+1). 
We can write 
uk+1 = uk㏒+1 + uk㏒㏒ +1, 

㏒ ㏒㏒ ﹋ U㏒㏒
where u﹋ U㏒ and uk . See Figure 12.6. Let 
k+1 kk+1 
rk+1,k+1 =請uk㏒㏒ +1請. 
If u㏒㏒ k+1 = rk+1,k+1 ek+1, we let hk+1 = id. Otherwise, there is a unique hyperplane re.ection hk+1 such that 
㏒㏒
hk+1(uk+1)= rk+1,k+1 ek+1, 
de.ned such that 
(u ﹞ wk+1)hk+1(u)= u . 2 請wk+1請2 wk+1 for all u ﹋ E, where 
㏒㏒ 
wk+1 = rk+1,k+1 ek+1 . uk+1. The map hk+1 is the re.ection about the hyperplane Hk+1 orthogonal to the vector wk+1 = 
㏒㏒ ㏒㏒
rk+1,k+1 ek+1 .uk+1. However, since uk+1,ek+1 ﹋ U㏒㏒ and Uk ㏒ is orthogonal to Uk ㏒㏒, the subspace Uk ㏒ iscontainedinHk+1,andthus,thevectors(r1,.k..,rk) and uk㏒ +1, which belong to Uk㏒ , are invariant under hk+1. This proves that 
㏒ ㏒㏒㏒
hk+1(uk+1)= hk+1(u )+ hk+1(u )= uk+1 + rk+1,k+1 ek+1
k+1k+1
is a linear combination of (e1,...,ek+1). Letting 
rk+1 = hk+1(uk+1)= uk㏒+1 + rk+1,k+1 ek+1, 
since uk+1 = hk .﹞ ﹞﹞. h2 . h1(vk+1), the vector 
rk+1 = hk+1 .﹞ ﹞﹞. h2 . h1(vk+1) 
is a linear combination of (e1,...,ek+1). See Figure 12.7. The coe.cient of rk+1 over ek+1 is rk+1,k+1 =請uk㏒㏒ +1請, which is nonnegative. This concludes the induction step, and thus the proof. 
12.2. QR-DECOMPOSITION USING HOUSEHOLDER MATRICES 

Remarks: 
(1) Since every hi is a hyperplane re.ection or the identity, 
老 = hn .﹞ ﹞﹞. h2 . h1 
is an isometry. 
(2) If we allow negative diagonal entries in R, the last isometry hn may be omitted. 

(3) Instead of picking rk,k =請uk㏒㏒請, which means that 


㏒㏒ 
wk = rk,k ek . uk, 
where 1 ≒ k ≒ n, it might be preferable to pick rk,k =.請u㏒㏒ k請 if this makes請wk請2 larger, in which case wk = rk,k ek + uk㏒㏒ . 
Indeed, since the de.nition of hk involves division by請wk請2, it is desirable to avoid division by very small numbers. 
(4) The method also applies to any m-tuple of vectors (v1,...,vm), with m ≒ n. Then R is an upper triangular m ℅ m matrix and Q is an n ℅ m matrix with orthogonal columns (QTQ = Im). We leave the minor adjustments to the method as an exercise to the reader 
Proposition 12.3 directly yields the QR-decomposition in terms of Householder transfor-mations (see Strang [63, 64], Golub and Van Loan [30], Trefethen and Bau [68], Kincaid and Cheney [39], or Ciarlet [14]). 

Figure 12.7: The construction of h2 and r2 = h2 . h1(v2) in Proposition 12.3. 
R = Hn ﹞﹞﹞ H2H1A. 
As a corollary, there is a pair of matrices Q, R, where Q is orthogonal and R is upper triangular, such that A = QR (a QR-decomposition of A). Furthermore, R can be chosen so that its diagonal entries are nonnegative. 
Proof. The jth column of A can be viewed as a vector vj over the canonical basis (e1,...,en) of En (where (ej)i = 1 if i = j, and 0 otherwise, 1 ≒ i, j ≒ n). Applying Proposition 12.3 to (v1,...,vn), there is a sequence of n isometries h1,...,hn such that hi is a hyperplane re.ection or the identity, and if (r1,...,rn) are the vectors given by 
rj = hn .﹞ ﹞﹞. h2 . h1(vj), 
then every rj is a linear combination of the vectors (e1,...,ej), 1 ≒ j ≒ n. Letting R be the matrix whose columns are the vectors rj, and Hi the matrix associated with hi, it is clear that 
R = Hn ﹞﹞﹞ H2H1A, 

12.2. QR-DECOMPOSITION USING HOUSEHOLDER MATRICES 
where R is upper triangular and every Hi is either a Householder matrix or the identity. However, hi . hi = id for all i,1 ≒ i ≒ n, and so 
vj = h1 . h2 .﹞ ﹞﹞. hn(rj) 
for all j,1 ≒ j ≒ n. But 老 = h1 . h2 .﹞ ﹞﹞. hn is an isometry represented by the orthogonal matrix Q = H1H2 ﹞﹞﹞ Hn. It is clear that A = QR, where R is upper triangular. As we noted in Proposition 12.3, the diagonal entries of R can be chosen to be nonnegative. 
Remarks: 
(1) Letting 
Ak+1 = Hk ﹞﹞﹞ H2H1A, 
with A1 = A,1 ≒ k ≒ n, the proof of Proposition 12.3 can be interpreted in terms of the computation of the sequence of matrices A1,...,An+1 = R. The matrix Ak+1 has 
the shape 

.
. 
............. 
k+1
℅℅℅ u ℅℅℅℅
1 
. . .... 
. . ....
0 ℅ . . .... 
k+1
00 ℅ u ℅℅℅℅
k 
k+1 
......000 ℅℅℅℅u....... 
k+1
Ak+1 = k+1 ,
000 u ℅℅℅℅
k+2 
... . .... 
... . .... 
... . .... 
k+1
000 u ℅℅℅℅
n.1 k+1
000 un ℅℅℅℅ where the (k + 1)th column of the matrix is the vector uk+1 = hk .﹞ ﹞﹞. h2 . h1(vk+1), 
and thus 

 
 
 
㏒ k+1 k+1 
u =u ,...,u 
k+1 1 k
and 
㏒㏒ k+1 k+1 k+1 
 

u =uk+1,u k+2,...,u .
k+1 n
If the last n . k . 1 entries in column k + 1 are all zero, there is nothing to do, and we let Hk+1 = I. Otherwise, we kill these n . k . 1 entries by multiplying Ak+1 on the left by the Householder matrix Hk+1 sending
 
 
k+1 k+1
0,..., 0,u k+1,...,u nto (0,..., 0,rk+1,k+1, 0,..., 0), 
k+1 k+1
where rk+1,k+1 =請(uk+1,...,un )請. 
(2) If A is invertible and the diagonal entries of R are positive, it can be shown that Q and R are unique. 

(3) If we allow negative diagonal entries in R, the matrix Hn may be omitted (Hn = I). 

(4) The method allows the computation of the determinant of A. We have 


det(A)=(.1)m r1,1 ﹞﹞﹞ rn,n, 
where m is the number of Householder matrices (not the identity) among the Hi. 
(5) The ※condition number§ of the matrix A is preserved (see Strang [64], Golub and Van Loan [30], Trefethen and Bau [68], Kincaid and Cheney [39], or Ciarlet [14]). This is very good for numerical stability. 

(6) The method also applies to a rectangular m ℅ n matrix. If m ≡ n, then R is an n ℅ n upper triangular matrix and Q is an m ℅ n matrix such that QTQ = In. 


The following Matlab functions implement the QR-factorization method of a real square (possibly singular) matrix A using Householder re.ections 
The main function houseqr computes the upper triangular matrix R obtained by applying Householder re.ections to A. It makes use of the function house, which computes a unit vector u such that given a vector x ﹋ Rp, the Householder transformation P = I . 2uuT sets to zero all entries in x but the .rst entry x1. It only applies if請x(2 : p)請 = |x2|+﹞﹞﹞+|xp| >
1 
0. Since computations are done in .oating point, we use a tolerance factor tol, and if請x(2 : p)請1 ≒ tol, then we return u = 0, which indicates that the corresponding Householder transformation is the identity. To make sure that請Px請 is as large as possible, we pick uu = x + sign(x1)請x請2 e1, where sign(z)=1 if z ≡ 0 and sign(z)= .1 if z< 0. Note that as a result, diagonal entries in R may be negative. We will take care of this issue later. 
function s = signe(x) % if x >= 0, then signe(x) = 1 % else if x < 0 then signe(x) = -1 % 
if x < 0 s = -1; else 
s =1; end end 
function [uu, u] = house(x) 
% This constructs the unnormalized vector uu 
% defining the Householder reflection that 
% zeros all but the first entries in x. 
% u is the normalized vector uu/||uu|| 


12.2. QR-DECOMPOSITION USING HOUSEHOLDER MATRICES 
% 

tol = 2*10^(-15); % tolerance 
uu = x; 
p = size(x,1); 
% computes l^1-norm of x(2:p,1) 
n1 = sum(abs(x(2:p,1))); 
ifn1 <=tol 

u = zeros(p,1); uu = u; 
else l = sqrt(x＊*x); % l^2 norm of x uu(1) = x(1) + signe(x(1))*l; u = uu/sqrt(uu＊*uu); 
end end 
The Householder transformations are recorded in an array u of n . 1 vectors. There are more e.cient implementations, but for the sake of clarity we present the following version. 
function [R, u] = houseqr(A) % This function computes the upper triangular R in the QR % factorization of A using Householder reflections, and an % implicit representation of Q as a sequence of n -1 % vectors u_i representing Householder reflections 
n = size(A, 1); R =A; u = zeros(n,n-1); for i = 1:n-1 
[~, u(i:n,i)] = house(R(i:n,i)); 
if u(i:n,i) == zeros(n -i + 1,1) 
R(i+1:n,i) = zeros(n -i,1); 
else 
R(i:n,i:n) = R(i:n,i:n) 
-2*u(i:n,i)*(u(i:n,i)＊*R(i:n,i:n)); 

end end end 
If only R is desired, then houseqr does the job. In order to obtain R, we need to compose the Householder transformations. We present a simple method which is not the most e.cient (there is a way to avoid multiplying explicity the Householder matrices). 
The function buildhouse creates a Householder re.ection from a vector v. 
function P = buildhouse(v,i) 
% This function builds a Householder reflection 
% [I0] 
% [0 PP] 
% from a Householder reflection 
% PP = I -2uu*uu＊ 
% where uu = v(i:n) 
% Ifuu=0thenP-I 
% 

n = size(v,1); if v(i:n) == zeros(n -i + 1,1) P = eye(n); 
else PP = eye(n -i + 1) -2*v(i:n)*v(i:n)＊; P = [eye(i-1) zeros(i-1, n -i + 1); zeros(n-i+ 1,i-1)PP]; 
end end 
The function buildQ builds the matrix Q in the QR-decomposition of A. 
function Q = buildQ(u) % Builds the matrix Q in the QR decomposition % of an nxn matrix A using Householder matrices, % where u is a representation of the n -1 % Householder reflection by a list u of vectors produced by % houseqr 
n = size(u,1); Q = buildhouse(u(:,1),1); for i = 2:n-1 
Q = Q*buildhouse(u(:,i),i); end end 
The function buildhouseQR computes a QR-factorization of A. At the end, if some entries on the diagonal of R are negative, it creates a diagonal orthogonal matrix P such that PR has nonnegative diagonal entries, so that A =(QP )(PR) is the desired QR-factorization of A. 
function [Q,R] = buildhouseQR(A) % % Computes the QR decomposition of a square 

12.3. SUMMARY 
% matrix A (possibly singular) using Householder reflections 
n = size(A,1); 
[R,u] = houseqr(A); 
Q = buildQ(u); 
% Produces a matrix R whose diagonal entries are 
% nonnegative 
P = eye(n); 
fori =1:n if R(i,i) < 0 P(i,i) = -1; end 
end 
Q = Q*P; R = P*R; 
end 
Example 12.1. Consider the matrix 
.
. 
...
1234 5
... 

234 

A = 

. 

6 4567 
Running the function buildhouseQR, we get 
345 

.
. 
... 

0.1826 0.8165 0.4001 0.3741 
0.3651 0.4082 .0.2546 .0.7970 

0.5477 .0.0000 .0.6910 0.4717 
0.7303 .0.4082 0.5455 .0.0488 

...

Q = 

and

.
. 
... 

5.4772 7.3030 9.1287 10.9545 00.8165 1.6330 2.4495 
0 .0.0000 0.0000 0.0000 
0 .0.00000 0.0000 

...

R = 

. 

Observe that A has rank 2. The reader should check that A = QR. 
Remark: Curiously, running Matlab built-in function qr, the same R is obtained (up to column signs) but a di.erent Q is obtained (the last two columns are di.erent). 

12.3 Summary 
The main concepts and results of this chapter are listed below: 
. 	
Symmetry (or re.ection) with respect to F and parallel to G. 

. 	
Orthogonal symmetry (or re.ection) with respect to F and parallel to G; re.ections, .ips. 

. 	
Hyperplane re.ections and Householder matrices. 

. 	
A key fact about re.ections (Proposition 12.2). 

. 	
QR-decomposition in terms of Householder transformations (Theorem 12.4). 



12.4 Problems 
Problem 12.1. (1) Given a unit vector (. sin 牟, cos 牟), prove that the Householder matrix determined by the vector (. sin 牟, cos 牟) is 
cos 2牟 sin 2牟 
. 
sin 2牟 . cos 2牟 
Give a geometric interpretation (i.e., why the choice (. sin 牟, cos 牟)?). 

(2) Given any matrix 
ab 
A = , 
cd Prove that there is a Householder matrix H such that AH is lower triangular, i.e., 
㏒
a0 
AH = 
㏒ 	d㏒
cfor some a㏒,c㏒,d㏒ ﹋ R. Problem 12.2. Given a Euclidean space E of dimension n, if h is a re.ection about some 
hyperplane orthogonal to a nonzero vector u and f is any isometry, prove that f . h . f.1 is the re.ection about the hyperplane orthogonal to f(u). Problem 12.3. (1) Given a matrix 
ab 
A = , 
cd prove that there are Householder matrices G, H such that cos 牟 sin 牟 ab cos . sin . 
GAH =	= D, 
sin 牟 . cos 牟 cd sin . . cos . where D is a diagonal matrix, i. the following equations hold: (b + c) cos(牟 + .)=(a . d) sin(牟 + .), (c . b) cos(牟 . .)=(a + d) sin(牟 . .). 

12.4. PROBLEMS 
(2) Discuss the solvability of the system. Consider the following cases: 
Case 1: a . d = a + d = 0. 
Case 2a: a . d = b + c = 0, a + d = 0. 
Case 2b: a . d = 0, b + c = 0, a + d = 0. 
Case 3a: a + d = c . b = 0, a . d = 0. 
Case 3b: a + d = 0, c . b = 0, a . d = 0. 
Case 4: a + d = 0, a . d = 0. Show that the solution in this case is 

1 b + cc . b 
牟 = arctan + arctan ,
2 a . da + d
1 b + cc . b 
. = arctan . arctan . 
2 a . da + d
If b = 0, show that the discussion is simpler: basically, consider c =0 or c = 0. 
(3) Expressing everything in terms of u = cot 牟 and v = cot ., show that the equations in (2) become 
(b + c)(uv . 1) = (u + v)(a . d), 
(c . b)(uv +1) = (.u + v)(a + d). 
Problem 12.4. Let A be an n ℅ n real invertible matrix. 
(1) Prove that ATA is symmetric positive de.nite. 
(2) Use the Cholesky factorization ATA = RTR with R upper triangular with positive di-agonal entries to prove that Q = AR.1 is orthogonal, so that A = QR is the QR-factorization of A. 
Problem 12.5. Modify the function houseqr so that it applies to an m ℅ n matrix with m ≡ n, to produce an m ℅ n upper-triangular matrix whose last m . n rows are zeros. 
Problem 12.6. The purpose of this problem is to prove that given any self-adjoint linear map f : E ↙ E (i.e., such that f. = f), where E is a Euclidean space of dimension n ≡ 3, given an orthonormal basis (e1,...,en), there are n . 2 isometries hi, hyperplane re.ections or the identity, such that the matrix of 
hn.2 .﹞ ﹞﹞. h1 . f . h1 .﹞ ﹞﹞. hn.2 
is a symmetric tridiagonal matrix. 
(1) Prove that for any isometry f : E ↙ E, we have f = f. = f.1 i. f . f = id. 
Prove that if f and h are self-adjoint linear maps (f. = f and h. = h), then h . f . h is a self-adjoint linear map. 
(2) Let Vk be the subspace spanned by (ek+1,...,en). Proceed by induction. For the base case, proceed as follows. 
Let 
f(e1)= a10 e1 + ﹞﹞﹞ + an0 en, 

and let 

r1, 2 =
  

a 

en
  

. 

00 
2e2 + ﹞﹞﹞ + a
n
Find an isometry h1 (re.ection or id) such that h1(f(e1) . a10 e1)= r1, 2 e2. Observe that 
w1 = r1, 2 e2 + a10 e1 . f(e1) ﹋ V1, and prove that h1(e1)= e1, so that h1 . f . h1(e1)= a 01e1 + r1, 2 e2. 
Let f1 = h1 . f . h1. Assuming by induction that 
fk = hk .﹞ ﹞﹞. h1 . f . h1 .﹞ ﹞﹞. hk has a tridiagonal matrix up to the kth row and column, 1 ≒ k ≒ n . 3, let 
kk k
fk(ek+1)= akek + ak+1ek+1 + ﹞﹞﹞ + anen, 
and let 

rk+1,k+2 =
  

a 

kk 
﹞﹞ + a
k+2ek+2 + ﹞ 
en
  

.

n
Find an isometry hk+1 (re.ection or id) such that 
kk
hk+1(fk(ek+1) . akek . ak+1ek+1)= rk+1,k+2 ek+2. 
Observe that 
kk 
wk+1 = rk+1,k+2 ek+2 + akek + ak+1ek+1 . fk(ek+1) ﹋ Vk+1, and prove that hk+1(ek)= ek and hk+1(ek+1)= ek+1, so that 
kk
hk+1 . fk . hk+1(ek+1)= akek + ak+1ek+1 + rk+1,k+2 ek+2. Let fk+1 = hk+1 . fk . hk+1, and .nish the proof. 
(3) Prove that given any symmetric n℅n-matrix A, there are n.2 matrices H1,...,Hn.2, Householder matrices or the identity, such that 
B = Hn.2 ﹞﹞﹞ H1AH1 ﹞﹞﹞ Hn.2 
is a symmetric tridiagonal matrix. 
(4) Write a computer program implementing the above method. 

12.4. PROBLEMS 
Problem 12.7. Recall from Problem 5.6 that an n ℅ n matrix H is upper Hessenberg if hjk =0 forall (j, k) such that j . k ≡ 0. Adapt the proof of Problem 12.6 to prove that given any n ℅ n-matrix A, there are n . 2 ≡ 1 matrices H1,...,Hn.2, Householder matrices or the identity, such that 
B = Hn.2 ﹞﹞﹞ H1AH1 ﹞﹞﹞ Hn.2 
is upper Hessenberg. 
Problem 12.8. The purpose of this problem is to prove that given any linear map f : E ↙ E, where E is a Euclidean space of dimension n ≡ 2, given an orthonormal basis (e1,...,en), there are isometries gi,hi, hyperplane re.ections or the identity, such that the matrix of 
gn .﹞ ﹞﹞. g1 . f . h1 .﹞ ﹞﹞. hn 
is a lower bidiagonal matrix, which means that the nonzero entries (if any) are on the main descending diagonal and on the diagonal below it. 
(1) Let U㏒ be the subspace spanned by (e1,...,ek) and Uk ㏒㏒ be the subspace spanned by 
(ek+1,...,en),k1 ≒ k ≒ n . 1. Proceed by induction For the base case, proceed as follows. Let v1 = f.(e1) and r1, 1 =請v1請. Find an isometry h1 (re.ection or id) such that 
h1(f . (e1)) = r1, 1e1. 
Observe that h1(f.(e1)) ﹋ U1㏒ , so that 
(h1(f . (e1)),ejㄘ =0 
for all j, 2 ≒ j ≒ n, and conclude that 
(e1,f . h1(ej)ㄘ =0 
for all j, 2 ≒ j ≒ n. Next let u1 = f . h1(e1)= u1㏒+ u1㏒㏒ , 
where u1 ㏒ ﹋ U1 ㏒ and u1 ㏒㏒ ﹋ U1 ㏒㏒, and let r2, 1 =請u1㏒㏒請. Find an isometry g1 (re.ection or id) such that 
g1(u1㏒㏒
)= r2, 1e2. 
Show that g1(e1)= e1, g1 . f . h1(e1)= u1 ㏒ + r2, 1e2, 
and that (e1,g1 . f . h1(ej)ㄘ =0 
for all j, 2 ≒ j ≒ n. At the end of this stage, show that g1 . f . h1 has a matrix such that all entries on its .rst row except perhaps the .rst are zero, and that all entries on the .rst column, except perhaps the .rst two, are zero. 
Assume by induction that some isometries g1,...,gk and h1,...,hk have been found, either re.ections or the identity, and such that 
fk = gk .﹞ ﹞﹞. g1 . f . h1 ﹞﹞﹞. hk 
has a matrix which is lower bidiagonal up to and including row and column k, where 1 ≒ k ≒ n . 2. Let vk+1 = fk . (ek+1)= vk㏒+1 + vk㏒㏒ +1, ㏒ ﹋ U㏒㏒
where v﹋ U㏒ and v㏒㏒ , and let rk+1,k+1 = v㏒㏒ . Find an isometry hk+1 (re.ection
k+1 kk+1 kk+1 
or id) such that 
hk+1(vk㏒㏒ +1)= rk+1,k+1ek+1. Show that if hk+1 is a re.ection, then U㏒ . Hk+1, where Hk+1 is the hyperplane de.ning the re.ection hk+1. Deduce that hk+1(v㏒ )k= v㏒ , and that 
k+1k+1
hk+1(fk . (ek+1)) = vk㏒ +1 + rk+1,k+1ek+1. 
Observe that hk+1(f. , so that 
k (ek+1)) ﹋ Uk㏒ +1
(hk+1(fk . (ek+1)),ejㄘ =0 
for all j, k +2 ≒ j ≒ n, and thus, 
(ek+1,fk . hk+1(ej)ㄘ =0 
for all j, k +2 ≒ j ≒ n. Next let 
㏒㏒ 
uk+1 = fk . hk+1(ek+1)= u +1 + u
k㏒k+1, 
㏒ ﹋ U㏒ ㏒㏒ ﹋ U㏒㏒ ㏒㏒

where uk+1 k+1 and uk+1 k+1, and let rk+2,k+1 = uk+1 . Find an isometry gk+1 (re.ection or id) such that 
㏒㏒ 
gk+1(uk+1)= rk+2,k+1ek+2. Show that if gk+1 is a re.ection, then Uk㏒ +1 . Gk+1, where Gk+1 is the hyperplane de.ning the re.ection gk+1. Deduce that gk+1(ei)= ei for all i,1 ≒ i ≒ k + 1, and that 
㏒ 
gk+1 . fk . hk+1(ek+1)= uk+1 + rk+2,k+1ek+2. 
Since by induction hypothesis, (ei,fk . hk+1(ej)ㄘ =0 
for all i, j,1 ≒ i ≒ k + 1, k +2 ≒ j ≒ n, and since gk+1(ei)= ei for all i,1 ≒ i ≒ k + 1, conclude that 
(ei,gk+1 . fk . hk+1(ej)ㄘ =0 
for all i, j,1 ≒ i ≒ k + 1, k +2 ≒ j ≒ n. Finish the proof. 


