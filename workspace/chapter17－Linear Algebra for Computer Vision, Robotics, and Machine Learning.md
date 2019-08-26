Chapter 17 
Computing Eigenvalues and Eigenvectors 
After the problem of solving a linear system, the problem of computing the eigenvalues and the eigenvectors of a real or complex matrix is one of most important problems of numerical linear algebra. Several methods exist, among which we mention Jacobi, Givens¨CHouseholder, divide-and-conquer, QR iteration, and Rayleigh¨CRitz; see Demmel [16], Trefethen and Bau [68], Meyer [48], Serre [57], Golub and Van Loan [30], and Ciarlet [14]. Typically, better performing methods exist for special kinds of matrices, such as symmetric matrices. 
In theory, given an n¡Án complex matrix A, if we could compute a Schur form A = UTU. , where T is upper triangular and U is unitary, we would obtain the eigenvalues of A, since they are the diagonal entries in T . However, this would require .nding the roots of a polynomial, but methods for doing this are known to be numerically very unstable, so this is not a practical method. 
= P .1
A common paradigm is to construct a sequence (Pk) of matrices such that Akk APk converges, in some sense, to a matrix whose eigenvalues are easily determined. For example, 
P .1
Ak = k APk could become upper triangular in the limit. Furthermore, Pk is typically a product of ¡°nice¡± matrices, for example, orthogonal matrices. 
For general matrices, that is, matrices that are not symmetric, the QR iteration algo-rithm, due to Rutishauser, Francis, and Kublanovskaya in the early 1960s, is one of the most e.cient algorithms for computing eigenvalues. A fascinating account of the history of the QR algorithm is given in Golub and Uhlig [29]. The QR algorithm constructs a se-quence of matrices (Ak), where Ak+1 is obtained from Ak by performing a QR-decomposition Ak = QkRk of Ak and then setting Ak+1 = RkQk, the result of swapping Qk and Rk. It is immediately veri.ed that Ak+1 = Q. kAkQk, so Ak and Ak+1 have the same eigenvalues, which are the eigenvalues of A. 
The basic version of this algorithm runs into di.culties with matrices that have several eigenvalues with the same modulus (it may loop or not ¡°converge¡± to an upper triangular matrix). There are ways of dealing with some of these problems, but for ease of exposition, we .rst present a simpli.ed version of the QR algorithm which we call basic QR algorithm. 
577 

We prove a convergence theorem for the basic QR algorithm, under the rather restrictive hypothesis that the input matrix A is diagonalizable and that its eigenvalues are nonzero and have distinct moduli. The proof shows that the part of Ak strictly below the diagonal converges to zero and that the diagonal entries of Ak converge to the eigenvalues of A. 
Since the convergence of the QR method depends crucially only on the fact that the part of Ak below the diagonal goes to zero, it would be highly desirable if we could replace A by a similar matrix U.AU easily computable from A and having lots of zero strictly below the diagonal. It turns out that there is a way to construct a matrix H = U.AU which is almost triangular, except that it may have an extra nonzero diagonal below the main diagonal. Such matrices called, Hessenberg matrices, are discussed in Section 17.2. An n ¡Á n diagonalizable Hessenberg matrix H having the property that hi+1i 
= 0 for i =1,...,n . 1 (such a matrix is called unreduced) has the nice property that its eigenvalues are all distinct. Since every Hessenberg matrix is a block diagonal matrix of unreduced Hessenberg blocks, it su.ces to compute the eigenvalues of unreduced Hessenberg matrices. There is a special case of particular interest: symmetric (or Hermitian) positive de.nite tridiagonal matrices. Such matrices must have real positive distinct eigenvalues, so the QR algorithm converges to a diagonal matrix. 
In Section 17.3, we consider techniques for making the basic QR method practical and more e.cient. The .rst step is to convert the original input matrix A to a similar matrix H in Hessenberg form, and to apply the QR algorithm to H (actually, to the unreduced blocks of H). The second and crucial ingredient to speed up convergence is to add shifts. 
A shift is the following step: pick some ¦Òk, hopefully close to some eigenvalue of A (in general, ¦Ën), QR-factor Ak . ¦ÒkI as 
Ak . ¦ÒkI = QkRk, 
and then form Ak+1 = RkQk + ¦ÒkI. 
It is easy to see that we still have Ak+1 = Q. kAkQk. A judicious choice of ¦Òk can speed up convergence considerably. If H is real and has pairs of complex conjugate eigenvalues, we can perform a double shift, and it can be arranged that we work in real arithmetic. 
The last step for improving e.ciency is to compute Ak+1 = Qk.AkQk without even per-forming a QR-factorization of Ak .¦ÒkI. This can be done when Ak is unreduced Hessenberg. Such a method is called QR iteration with implicit shifts. There is also a version of QR iteration with implicit double shifts. 
If the dimension of the matrix A is very large, we can .nd approximations of some of the eigenvalues of A by using a truncated version of the reduction to Hessenberg form due to Arnoldi in general and to Lanczos in the symmetric (or Hermitian) tridiagonal case. Arnoldi iteration is discussed in Section 17.4. If A is an m ¡Á m matrix, for n¡¶ m (n much smaller than m) the idea is to generate the n ¡Á n Hessenberg submatrix Hn of the full Hessenberg matrix H (such that A = UHU.) consisting of its .rst n rows and n columns; the matrix Un consisting of the .rst n columns of U is also produced. The Rayleigh¨CRitz method consists 
17.1. THE BASIC QR ALGORITHM 
in computing the eigenvalues of Hn using the QR-method with shifts. These eigenvalues, called Ritz values, are approximations of the eigenvalues of A. Typically, extreme eigenvalues are found .rst. 
Arnoldi iteration can also be viewed as a way of computing an orthonormal basis of a Krylov subspace, namely the subspace Kn(A, b) spanned by (b, Ab, . . . , Anb). We can also use Arnoldi iteration to .nd an approximate solution of a linear equation Ax = b by minimizing lb . Axnl2 for all xn is the Krylov space Kn(A, b). This method named GMRES is discussed in Section 17.5. 
The special case where H is a symmetric (or Hermitian) tridiagonal matrix is discussed in Section 17.6. In this case, Arnoldi¡¯s algorithm becomes Lanczos¡¯ algorithm. It is much more e.cient than Arnoldi iteration. 
We close this chapter by discussing two classical methods for computing a single eigen-vector and a single eigenvalue: power iteration and inverse (power) iteration; see Section 
17.7. 

17.1 The Basic QR Algorithm 
Let A be an n ¡Á n matrix which is assumed to be diagonalizable and invertible. The basic QR algorithm makes use of two very simple steps. Starting with A1 = A, we construct sequences of matrices (Ak), (Qk)(Rk) and (Pk) as follows: 
Factor A1 = Q1R1 Set A2 = R1Q1 Factor A2 = Q2R2 Set A3 = R2Q2 
. 
. 
. Factor Ak = QkRk Set Ak+1 = RkQk 
. 
. 
. 
Thus, Ak+1 is obtained from a QR-factorization Ak = QkRk of Ak by swapping Qk and Rk. De.ne Pk by Pk = Q1Q2 ¡¤¡¤¡¤ Qk. 
Since Ak = QkRk, we have Rk = Qk.Ak, and since Ak+1 = RkQk, we obtain 
= Q . 
Ak+1 kAkQk. (. 1) 
An obvious induction shows that 
= Q . ¡¤¡¤ Q . A1Q1 ¡¤¡¤¡¤ Qk = P . APk,
Ak+1 k ¡¤ 1k 
that is Ak+1 = Pk . APk. (. 2) 
Therefore, Ak+1 and A are similar, so they have the same eigenvalues. 
The basic QR iteration method consists in computing the sequence of matrices Ak, and in the ideal situation, to expect that Ak ¡°converges¡± to an upper triangular matrix, more precisely that the part of Ak below the main diagonal goes to zero, and the diagonal entries converge to the eigenvalues of A. 
This ideal situation is only achieved in rather special cases. For one thing, if A is unitary (or orthogonal in the real case), since in the QR decomposition we have R = I, we get A2 = IQ = Q = A1, so the method does not make any progress. Also, if A is a real matrix, since the Ak are also real, if A has complex eigenvalues, then the part of Ak below the main diagonal can¡¯t go to zero. Generally, the method runs into troubles whenever A has distinct eigenvalues with the same modulus. 
The convergence of the sequence (Ak) is only known under some fairly restrictive hy-potheses. Even under such hypotheses, this is not really genuine convergence. Indeed, it can be shown that the part of Ak below the main diagonal goes to zero, and the diagonal entries converge to the eigenvalues of A, but the part of Ak above the diagonal may not converge. However, for the purpose of .nding the eigenvalues of A, this does not matter. 
The following convergence result is proven in Ciarlet [14] (Chapter 6, Theorem 6.3.10 and Serre [57] (Chapter 13, Theorem 13.2). It is rarely applicable in practice, except for symmetric (or Hermitian) positive de.nite matrices, as we will see shortly. 
Theorem 17.1. Suppose the (complex) n¡Án matrix A is invertible, diagonalizable, and that its eigenvalues ¦Ë1,...,¦Ën have di.erent moduli, so that 
|¦Ë1| > |¦Ë2| > ¡¤¡¤¡¤ > |¦Ën| > 0. 
If A = P ¦«P .1, where ¦« = diag(¦Ë1,...,¦Ën), and if P .1 has an LU-factorization, then the strictly lower-triangular part of Ak converges to zero, and the diagonal of Ak converges to ¦«. 
Proof. We reproduce the proof in Ciarlet [14] (Chapter 6, Theorem 6.3.10). The strategy is to study the asymptotic behavior of the matrices Pk = Q1Q2 ¡¤¡¤¡¤ Qk. For this, it turns out that we need to consider the powers Ak . 
Step 1 . Let Rk = Rk ¡¤¡¤¡¤ R2R1. We claim that 
Ak =(Q1Q2 ¡¤¡¤¡¤ Qk)(Rk ¡¤¡¤¡¤ R2R1)= PkRk. (. 3) 
We proceed by induction. The base case k = 1 is trivial. For the induction step, from (. 2), we have PkAk+1 = APk. 
Since Ak+1 = RkQk = Qk+1Rk+1, we have 
= Ak+1
Pk+1Rk+1 = PkQk+1Rk+1Rk = PkAk+1Rk = APkRk = AAk 

17.1. THE BASIC QR ALGORITHM 
establishing the induction step. 
Step 2 . We will express the matrix Pk as Pk = QQ
kDk, in terms of a diagonal matrix Dk with unit entries, with Q and Q
k unitary. 
Let P = QR,a QR-factorization of P (with R an upper triangular matrix with positive diagonal entries), and P .1 = LU, an LU-factorization of P .1 . Since A = P ¦«P .1, we have 
Ak = P ¦«kP .1 = QR¦«kLU = QR(¦«kL¦«.k)¦«kU. (. 4) 
Here, ¦«.k is the diagonal matrix with entries ¦Ë.ik . The reason for introducing the matrix ¦«kL¦«.k is that its asymptotic behavior is easy to determine. Indeed, we have 
. .
. 

0 if i<j 1 
if i = j

(¦«kL¦«.k)ij 
=

 
.
. 
 k
¦Ëi Lij if i > j. ¦Ëj
QQ
The hypothesis that |¦Ë1| > |¦Ë2| > ¡¤¡¤¡¤ > |¦Ën| > 0 implies that 
lim ¦«kL¦«.k = I. (.) 
kh¡ú¡Þ 
Note that it is to obtain this limit that we made the hypothesis on the moduli of the eigenvalues. We can write 
¦«kL¦«.k = I + Fk, with lim Fk =0, 
kh¡ú¡Þ 
and consequently, since R(¦«kL¦«.k)= R(I + Fk)= R + RFkR.1R =(I + RFkR.1)R, we have 
R(¦«kL¦«.k)=(I + RFkR.1)R. (. 5) 
By Proposition 8.11(1), since limkh¡ú¡Þ Fk = 0, and thus limkh¡ú¡Þ RFkR.1 = 0, the matrices I + RFkR.1 are invertible for k large enough. Consequently for k large enough, we have a QR-factorization 
Rk, 
Qk 
I + RFkR.1 
= Q
k(. 6)
   

   

with (RQk)ii 
Q
Qk
> 0 for i =1,...,n. Since the matrices

are unitary, we have

1,

= 
2 
Q
so the sequence (Q
k) is bounded. It follows that it has a convergent subsequence (Q
´¨) that converges to some matrix Q
, which is also unitary. Since 
R´¨ =(Q
´¨) . (I + RF´¨R.1), 
we deduce that the subsequence (RQ´¨) also converges to some matrix RQ, which is also upper triangular with positive diagonal entries. By passing to the limit (using the subsequences), we get RQ=(Q
)., that is, 
I R.
= Q
Q
By the uniqueness of a QR-decomposition (when the diagonal entries of R are positive), we get 
QQ
Q = R = I. 
Since the above reasoning applies to any subsequences of ( QQk) and ( RQk), by the uniqueness of limits, we conclude that the ¡°full¡± sequences ( QQk) and ( RQk) converge: 
QQ
lim Qk = I, lim Rk = I. 
kh¡ú¡Þ kh¡ú¡Þ 
Since by (. 4), 
Ak = QR(¦«kL¦«.k)¦«kU, 

by (. 5), R(¦«kL¦«.k)=(I + RFkR.1)R, 
and by (. 6) 
QQ
I + RFkR.1 = QkRk, 
we proved that Ak =(QQQk)(RQkR¦«kU). (. 7) 
Observe that QQQk is a unitary matrix, and RQkR¦«kU is an upper triangular matrix, as a product of upper triangular matrices. However, some entries in ¦« may be negative, so 
Q
we can¡¯t claim that RkR¦«kU has positive diagonal entries. Nevertheless, we have another QR-decomposition of Ak , 
Ak =(QQQk)(RQkR¦«kU)= PkRk. 
It is easy to prove that there is diagonal matrix Dk with |(Dk)ii| = 1 for i =1,...,n, such that 
Pk = QQQkDk. (. 8) 
The existence of Dk is consequence of the following fact: If an invertible matrix B has two QR factorizations B = Q1R1 = Q2R2, then there is a diagonal matrix D with unit entries such that Q2 = DQ1. 
The expression for Pk in (. 8) is that which we were seeking. 
Step 3 . Asymptotic behavior of the matrices Ak+1 = Pk .APk. 
Since A = P ¦«P .1 = QR¦«R.1Q.1 and by (. 8), Pk = QQQkDk, we get 
Ak+1 = Dk. (QQk) . Q . QR¦«R.1Q.1QQQkDk = Dk. (QQk) . R¦«R.1QQkDk. (. 9) 
Q
Since limkh¡ú¡Þ Qk = I, we deduce that 
.. 
¦Ë1 . ¡¤¡¤¡¤ . 
..
0 ¦Ë2 ¡¤¡¤¡¤ .
..
lim (QQk) . R¦«R.1QQk = R¦«R.1 = . .. . ,.
kh¡ú¡Þ. . .. .
. 00 ¡¤¡¤¡¤ ¦Ën .. 

17.1. THE BASIC QR ALGORITHM 
an upper triangular matrix with the eigenvalues of A on the diagonal. Since R is upper triangular, the order of the eigenvalues is preserved. If we let 
Dk =( QQk) . R¦«R.1QQk, (. 10) 
then by (. 9) we have Ak+1 = Dk.DkDk, and since the matrices Dk are diagonal matrices, we have 
(Ak+1)jj =(D .DkDk)ij =(Dk)ii(Dk)jj(Dk)ij,
k
which implies that 
(Ak+1)ii =(Dk)ii,i =1, . . . , n, (. 11) since |(Dk)ii| = 1 for i =1,...,n. Since limkh¡ú¡Þ Dk = R¦«R.1, we conclude that the strictly lower-triangular part of Ak+1 converges to zero, and the diagonal of Ak+1 converges to ¦«. 
Observe that if the matrix A is real, then the hypothesis that the eigenvalues have distinct moduli implies that the eigenvalues are all real and simple. 
The following Matlab program implements the basic QR-method using the function qrv4 from Section 11.8. 
function T = qreigen(A,m) 
T =A; 
fork =1:m 

[Q R] = qrv4(T); 
T = R*Q; end end 
Example 17.1. If we run the function qreigen with 100 iterations on the 8 ¡Á 8 symmetric 
matrix

.
. 
........... 

4  1  0  0  0  0  0  0  
1  4  1  0  0  0  0  0  
0  1  4  1  0  0  0  0  
0  0  1  4  1  0  0  0  

0  0  0  1  4  1  0  0  
0  0  0  0  1  4  1  0  
0  0  0  0  0  1  4  1  
0  0  0  0  0  0  1  4  

........... 

A = 

, 

we .nd the matrix 

.
. 
........... 

5.8794 0.0015 0.0000 .0.0000 0.0000 .0.0000 0.0000 .0.0000 0.0015 5.5321 0.0001 0.0000 .0.0000 0.0000 0.0000 0.0000 0 
0.0001 5.0000 0.0000 .0.0000 0.0000 0.0000 0.0000 0 00.0000 4.3473 0.0000 0.0000 0.0000 0.0000 
0  0  0  0.0000  3.6527  0.0000  0.0000  .0.0000  
0  0  0  0  0.0000  3.0000  0.0000  .0.0000  
0  0  0  0  0  0.0000  2.4679  0.0000  
0  0  0  0  0  0  0.0000  2.1.206  

........... 

T = 

. 

The diagonal entries match the eigenvalues found by running the Matlab function eig(A). 
If several eigenvalues have the same modulus, then the proof breaks down, we can no longer claim (.), namely that 
¦«kL¦«.k
lim = I. 
kh¡ú¡Þ 
If we assume that P .1 has a suitable ¡°block LU-factorization,¡± it can be shown that the matrices Ak+1 converge to a block upper-triangular matrix, where each block corresponds to eigenvalues having the same modulus. For example, if A is a 9 ¡Á 9 matrix with eigenvalues ¦Ëi such that |¦Ë1| = |¦Ë2| = |¦Ë3| > |¦Ë4| > |¦Ë5| = |¦Ë6| = |¦Ë7| = |¦Ë8| = |¦Ë9|, then Ak converges to a block diagonal matrix (with three blocks, a 3 ¡Á 3 block, a 1 ¡Á 1 block, and a 5 ¡Á 5 block) 
of the form

.
. 
*** ...... *** ...... *** ...... 
............. 
............. 

000 * ..... 

0000 ***** 

. 

0000 ***** 0000 ***** 0000 ***** 0000 ***** 
See Ciarlet [14] (Chapter 6 Section 6.3) for more details. 
Under the conditions of Theorem 17.1, in particular, if A is a symmetric (or Hermitian) positive de.nite matrix, the eigenvectors of A can be approximated. However, when A is not a symmetric matrix, since the upper triangular part of Ak does not necessarily converge, one has to be cautious that a rigorous justi.cation is lacking. 
Suppose we apply the QR algorithm to a matrix A satisfying the hypotheses of Theorem Theorem 17.1. For k large enough, Ak+1 = Pk .APk is nearly upper triangular and the diagonal entries of Ak+1 are all distinct, so we can consider that they are the eigenvalues of Ak+1, and thus of A. To avoid too many subscripts, write T for the upper triangular matrix obtained by settting the entries of the part of Ak+1 below the diagonal to 0. Then we can .nd the corresponding eigenvectors by solving the linear system 
Tv = tiiv, 
and since T is upper triangular, this can be done by bottom-up elimination. We leave it as an exercise to show that the following vectors vi =(v1i ,...,vni ) are eigenvectors: 
1 
v = e1, and if i =2,...,n, then 
. .
. .
. 

0 if i +1 ¡Ü j ¡Ü n 

i 1 if j = i 
vj = i ¡¤ i
tjj+1vj+1 + ¡¤¡¤ + tjivi
. if i . 1 ¡Ý j ¡Ý 1. tjj . tii 

17.2. HESSENBERG MATRICES 
Then the vectors (Pkv1,...,Pkvn) are a basis of (approximate) eigenvectors for A. In the special case where T is a diagonal matrix, then vi = ei for i =1,...,n and the columns of Pk are an orthonormal basis of (approximate) eigenvectors for A. 
If A is a real matrix whose eigenvalues are not all real, then there is some complex pair of eigenvalues ¦Ë + i¦Ì (with ¦Ì = 0), and the QR-algorithm cannot converge to a matrix whose strictly lower-triangular part is zero. There is a way to deal with this situation using upper Hessenberg matrices which will be discussed in the next section. 
Since the convergence of the QR method depends crucially only on the fact that the part of Ak below the diagonal goes to zero, it would be highly desirable if we could replace A by a similar matrix U.AU easily computable from A having lots of zero strictly below the diagonal. We can¡¯t expect U.AU to be a diagonal matrix (since this would mean that A was easily diagonalized), but it turns out that there is a way to construct a matrix H = U.AU which is almost triangular, except that it may have an extra nonzero diagonal below the main diagonal. Such matrices called Hessenberg matrices are discussed in the next section. 

17.2 Hessenberg Matrices 
De.nition 17.1. An n ¡Á n matrix (real or complex) H is an (upper) Hessenberg matrix if it is almost triangular, except that it may have an extra nonzero diagonal below the main diagonal. Technically, hjk =0 for all (j, k) such that j . k ¡Ý 2. 
The 5 ¡Á 5 matrix below is an example of a Hessenberg matrix. 
.
. 
H = 

..... 

. . . .. h21 . . .. 
0  h32  .  .  .  
0  0  h43  .  .  
0  0  0  h54  .  

..... 

. 

The following result can be shown. 
Theorem 17.2. Every n ¡Á n complex or real matrix A is similar to an upper Hessenberg matrix H, that is, A = UHU. for some unitary matrix U. Furthermore, H can be constructed as a product of Householder matrices (the de.nition is the same as in Section 12.1, except that W is a complex vector, and that the inner product is the Hermitian inner product on Cn). If A is a real matrix, then H is an orthogonal matrix (and H is a real matrix). 
Theorem 17.2 and algorithms for converting a matrix to Hessenberg form are discussed in Trefethen and Bau [68] (Lecture 26), Demmel [16] (Section 4.4.6, in the real case), Serre [57] (Theorem 13.1), and Meyer [48] (Example 5.7.4, in the real case). The proof of correctness is not di.cult and will be the object of a homework problem. 
The following functions written in Matlab implement a function to compute a Hessenberg form of a matrix. 
The function house constructs the normalized vector u de.ning the Householder re.ection that zeros all but the .rst entries in a vector x. 
function [uu, u] = house(x) 
tol = 2*10^(-15); % tolerance 
uu = x; 
p = size(x,1); 
% computes l^1-norm of x(2:p,1) 
n1 = sum(abs(x(2:p,1))); 
ifn1 <=tol 

u = zeros(p,1); uu = u; 
else l = sqrt(x¡¯*x); % l^2 norm of x uu(1) = x(1) + signe(x(1))*l; u = uu/sqrt(uu¡¯*uu); 
end end 
The function signe(z) returms .1 if z< 0, else +1. 
The function buildhouse builds a Householder re.ection from a vector uu. 

function P = buildhouse(v,i) 
% This function builds a Householder reflection 
% [I0] 
% [0 PP] 
% from a Householder reflection 
% PP = I -2uu*uu¡¯ 
% where uu = v(i:n) 
% Ifuu=0thenP-I 
% 

n = size(v,1); if v(i:n) == zeros(n -i + 1,1) P = eye(n); 
else PP = eye(n -i + 1) -2*v(i:n)*v(i:n)¡¯; P = [eye(i-1) zeros(i-1, n -i + 1); 
zeros(n -i + 1, i -1) PP]; end end 
The function Hessenberg1 computes an upper Hessenberg matrix H and an orthogonal matrix Q such that A = QTHQ. 
function [H, Q] = Hessenberg1(A) % 
17.2. HESSENBERG MATRICES 
%  This  function  constructs  an  upper  Hessenberg  
%  matrix  H  and  an  orthogonal  matrix  Q  such  that  
%  A  =  Q¡¯  H  Q  
n  =  size(A,1);  
H  =  A;  
Q  =  eye(n);  
for  i  =  1:n-2  
%  H(i+1:n,i)  
[~,u]  =  house(H(i+1:n,i));  
%  u  
P  =  buildhouse(u,1);  
Q(i+1:n,i:n)  =  P*Q(i+1:n,i:n);  
H(i+1:n,i:n)  =  H(i+1:n,i:n)  -2*u*(u¡¯)*H(i+1:n,i:n);  
H(1:n,i+1:n)  =  H(1:n,i+1:n)  -2*H(1:n,i+1:n)*u*(u¡¯);  
end  
end  
Example 17.2. If



.
. 
...
1234 5
... 

234 

A = 

,

6 4567 
running Hessenberg1 we .nd 
345 

.
. 
... 

1.0000 .5.38520 0 .5.3852 15.2069 .1.6893 .0.0000 
.0.0000 .1.6893 .0.2069 .0.0000 0 .0.0000 0.0000 0.0000 
...

H = 

.
. 
... 

1.00000 0 0 0 .0.3714 .0.5571 .0.7428 
00.8339 0.1516 .0.5307 
00.4082 .0.8165 0.4082 

...

Q = 

. 

An important property of (upper) Hessenberg matrices is that if some subdiagonal entry Hp+1p = 0, then H is of the form 
H11 H12
H = ,
0 H22 
where both H11 and H22 are upper Hessenberg matrices (with H11 a p ¡Á p matrix and H22 a (n . p) ¡Á (n . p) matrix), and the eigenvalues of H are the eigenvalues of H11 and H22. For 
example, in the matrix 

.
. 
. . . .. h21 . . .. 
H = 

..... 

0  h32  .  .  .  
0  0  h43  .  .  
0  0  0  h54  .  

..... 

, 

if h43 = 0, then we have the block matrix 
.
. 
. .... h21 .... 0 ... 
..... 

..... 

H 

= 

h32 
. 

0 00 .. 0 00 h54 . 
Then the list of eigenvalues of H is the concatenation of the list of eigenvalues of H11 and the list of the eigenvalues of H22. This is easily seen by induction on the dimension of the block H11. 
More generally, every upper Hessenberg matrix can be written in such a way that it has diagonal blocks that are Hessenberg blocks whose subdiagonal is not zero. 
De.nition 17.2. An upper Hessenberg n ¡Á n matrix H is unreduced if hi+1i = 0 for i = 1,...,n . 1. A Hessenberg matrix which is not unreduced is said to be reduced. 
The following is an example of an 8 ¡Á 8 matrix consisting of three diagonal unreduced 
Hessenberg blocks: 

.
. 
........... 

* ** . .... h21 ** . .... 0h32 * . .... 0 00 * ** .. 
0 00 h54 ** .. 0 00 0h65 * .. 0 000 00 ** 0 000 00 h87 * 
........... 

H = 

. 

An interesting and important property of unreduced Hessenberg matrices is the following. 
Proposition 17.3. Let H be an n ¡Á n complex or real unreduced Hessenberg matrix. Then every eigenvalue of H is geometrically simple, that is, dim(E¦Ë)=1 for every eigenvalue ¦Ë, where E¦Ë is the eigenspace associated with ¦Ë. Furthermore, if H is diagonalizable, then every eigenvalue is simple, that is, H has n distinct eigenvalues. 
Proof. We follow Serre¡¯s proof [57] (Proposition 3.26). Let ¦Ë be any eigenvalue of H, let M = ¦ËIn . H, and let N be the (n . 1) ¡Á (n . 1) matrix obtained from M by deleting its .rst row and its last column. Since H is upper Hessenberg, N is a diagonal matrix with entries .hi+1i = 0, i =1,...,n . 1. Thus N is invertible and has rank n . 1. But a matrix 

17.2. HESSENBERG MATRICES 
has rank greater than or equal to the rank of any of its submatrices, so rank(M)= n . 1, since M is singular. By the rank-nullity theorem, rank(Ker N) = 1, that is, dim(E¦Ë) = 1, as claimed. 
If H is diagonalizable, then the sum of the dimensions of the eigenspaces is equal to n, which implies that the eigenvalues of H are distinct. 
As we said earlier, a case where Theorem 17.1 applies is the case where A is a symmetric (or Hermitian) positive de.nite matrix. This follows from two facts. 
The .rst fact is that if A is Hermitian (or symmetric in the real case), then it is easy to show that the Hessenberg matrix similar to A is a Hermitian (or symmetric in real case) tridiagonal matrix. The conversion method is also more e.cient. Here is an example of a symmetric tridiagonal matrix consisting of three unreduced blocks: 
.
. 
H = 

........... 

¦Á1 ¦Â1 0 00000 ¦Â1 ¦Á2 ¦Â2 00000 0 ¦Â2 ¦Á3 00000 000 ¦Á4 ¦Â4 0 00 
0  0  0  ¦Â4  ¦Á5  ¦Â5  0  0  
0  0  0  0  ¦Â5  ¦Á6  0  0  
0  0  0  0  0  0  ¦Á7  ¦Â7  
0  0  0  0  0  0  ¦Â7  ¦Á8  

........... 

. 

Thus the problem of .nding the eigenvalues of a symmetric (or Hermitian) matrix reduces to the problem of .nding the eigenvalues of a symmetric (resp. Hermitian) tridiagonal matrix, and this can be done much more e.ciently. 
The second fact is that if H is an upper Hessenberg matrix and if it is diagonalizable, then there is an invertible matrix P such that H = P ¦«P .1 with ¦« a diagonal matrix consisting of the eigenvalues of H, such that P .1 has an LU-decomposition; see Serre [57] (Theorem 13.3). 
As a consequence, since any symmetric (or Hermitian) tridiagonal matrix is a block diag-onal matrix of unreduced symmetric (resp. Hermitian) tridiagonal matrices, by Proposition 17.3, we see that the QR algorithm applied to a tridiagonal matrix which is symmetric (or Hermitian) positive de.nite converges to a diagonal matrix consisting of its eigenvalues. Let us record this important fact. 
Theorem 17.4. Let H be a symmetric (or Hermitian) positive de.nite tridiagonal matrix. If H is unreduced, then the QR algorithm converges to a diagonal matrix consisting of the eigenvalues of H. 
Since every symmetric (or Hermitian) positive de.nite matrix is similar to tridiagonal symmetric (resp. Hermitian) positive de.nite matrix, we deduce that we have a method for .nding the eigenvalues of a symmetric (resp. Hermitian) positive de.nite matrix (more accurately, to .nd approximations as good as we want for these eigenvalues). 
If A is a symmetric (or Hermitian) matrix, since its eigenvalues are real, for some ¦Ì> 0 large enough (pick ¦Ì>¦Ñ(A)), A + ¦ÌI is symmetric (resp. Hermitan) positive de.nite, so we can apply the QR algorithm to an upper Hessenberg matrix similar to A + ¦ÌI to .nd its eigenvalues, and then the eigenvalues of A are obtained by subtracting ¦Ì. 
The problem of .nding the eigenvalues of a symmetric matrix is discussed extensively in Parlett [50], one of the best references on this topic. 
The upper Hessenberg form also yields a way to handle singular matrices. First, checking the proof of Proposition 13.21 that an n ¡Á n complex matrix A (possibly singular) can be factored as A = QR where Q is a unitary matrix which is a product of Householder re.ections and R is upper triangular, it is easy to see that if A is upper Hessenberg, then Q is also upper Hessenberg. If H is an unreduced upper Hessenberg matrix, since Q is upper Hessenberg and R is upper triangular, we have hi+1i = qi+1irii for i =1 ...,n . 1, and since H is unreduced, rii = 0 for i =1,...,n . 1. Consequently H is singular i. rnn = 0. Then the matrix RQ is a matrix whose last row consists of zero¡¯s thus we can de.ate the problem by considering the (n . 1) ¡Á (n . 1) unreduced Hessenberg matrix obtained by deleting the last row and the last column. After .nitely many steps (not larger that the multiplicity of the eigenvalue 0), there remains an invertible unreduced Hessenberg matrix. As an alternative, see Serre [57] (Chapter 13, Section 13.3.2). 
As is, the QR algorithm, although very simple, is quite ine.cient for several reasons. In the next section, we indicate how to make the method more e.cient. This involves a lot of work and we only discuss the main ideas at a high level. 
17.3 	Making the QR Method More E.cient Using Shifts 
To improve e.ciency and cope with pairs of complex conjugate eigenvalues in the case of real matrices, the following steps are taken: 
(1) 
Initially reduce the matrix A to upper Hessenberg form, as A = UHU. . Then apply the QR-algorithm to H (actually, to its unreduced Hessenberg blocks). It is easy to see that the matrices Hk produced by the QR algorithm remain upper Hessenberg. 

(2) 
To 	accelerate convergence, use shifts, and to deal with pairs of complex conjugate eigenvalues, use double shifts. 

(3) 
Instead of computing 	a QR-factorization explicitly while doing a shift, perform an implicit shift which computes Ak+1 = Q. kAkQk without having to compute a QR-factorization (of Ak . ¦ÒkI), and similarly in the case of a double shift. This is the most intricate modi.cation of the basic QR algorithm and we will not discuss it here. This method is usually referred as bulge chasing. Details about this technique for real matrices can be found in Demmel [16] (Section 4.4.8) and Golub and Van Loan 


[30] (Section 7.5). Watkins discusses the QR algorithm with shifts as a bulge chasing method in the more general case of complex matrices [73, 74]. 

17.3. MAKING THE QR METHOD MORE EFFICIENT USING SHIFTS 
Let us repeat an important remark made in the previous section. If we start with a matrix H in upper Hessenberg form, if at any stage of the QR algorithm we .nd that some subdiagonal entry (Hk)p+1p =0 or is very small, then Hk is of the form 
H11 H12
Hk = ,
0 H22 
where both H11 and H22 are upper Hessenberg matrices (with H11 a p ¡Á p matrix and H22 a(n . p) ¡Á (n . p) matrix), and the eigenvalues of Hk are the eigenvalues of H11 and H22. For example, in the matrix 
.
. 
H = 

..... 

. . . .. h21 . . .. 
0  h32  .  .  .  
0  0  h43  .  .  
0  0  0  h54  .  

..... 

, 

if h43 = 0, then we have the block matrix 
.
. 
..... 

. .... h21 .... 
0  h32  .  .  .  
0  0  0  .  .  
0  0  0  h54  .  

..... 

H = 

. 

Then we can recursively apply the QR algorithm to H11 and H22. 
In particular, if (Hk)nn.1 = 0 or is very small, then (Hk)nn is a good approximation of an eigenvalue, so we can delete the last row and the last column of Hk and apply the QR algorithm to this submatrix. This process is called de.ation. If(Hk)n.1n.2 = 0 or is very small, then the 2 ¡Á 2 ¡°corner block¡± 
(Hk)n.1n.1 (Hk)n.1n (Hk)nn.1 (Hk)nn 
appears, and its eigenvalues can be computed immediately by solving a quadratic equation. Then we de.ate Hk by deleting its last two rows and its last two columns and apply the QR algorithm to this submatrix. 
Thus it would seem desirable to modify the basic QR algorithm so that the above situ-ations arises, and this is what shifts are designed for. More precisely, under the hypotheses of Theorem 17.1, it can be shown (see Ciarlet [14], Section 6.3) that the entry (Ak)ij with 
        
i>j converges to 0 as |¦Ëi/¦Ëj|k ¦Ë2 converges to 0. Also, if we let ri be de.ned by 
 

    

    

    

¦Ëi+1 
    

 

  

   
 ,
r1 =¦Ë1
, 

ri = max
¦Ëi ¦Ëi.1
,

¦Ëi
, 

2 ¡Ü i ¡Ü n . 1, 

rn =
  

¦Ën ¦Ën.1
then there is a constant C (independent of k) such that |(Ak)ii . ¦Ëi|¡Ü Crik , 1 ¡Ü i ¡Ü n. 
In particular, if H is upper Hessenberg, then the entry (Hk)i+1i converges to 0 as |¦Ëi+1/¦Ëi|k converges to 0. Thus if we pick ¦Òk close to ¦Ëi, we expect that (Hk . ¦ÒkI)i+1i converges to 0 as |(¦Ëi+1 .¦Òk)/(¦Ëi .¦Òk)|k converges to 0, and this ratio is much smaller than 1 as ¦Òk is closer to ¦Ëi. Typically, we apply a shift to accelerate convergence to ¦Ën (so i = n . 1). In this case, both (Hk . ¦ÒkI)nn.1 and |(Hk . ¦ÒkI)nn . ¦Ën| converge to 0 as |(¦Ën . ¦Òk)/(¦Ën.1 . ¦Òk)|k converges to 0. 
A shift is the following modi.ed QR-steps (switching back to an arbitrary matrix A, since the shift technique applies in general). Pick some ¦Òk, hopefully close to some eigenvalue of A (in general, ¦Ën), and QR-factor Ak . ¦ÒkI as 
Ak . ¦ÒkI = QkRk, 
and then form Ak+1 = RkQk + ¦ÒkI. 
Since 
Ak+1 = RkQk + ¦ÒkI = Q . kQkRkQk + Q . kQk¦Òk 
= Q . 
k(QkRk + ¦ÒkI)Qk = Q . 
kAkQk, 
Ak+1 is similar to Ak, as before. If Ak is upper Hessenberg, then it is easy to see that Ak+1 is also upper Hessenberg. 
If A is upper Hessenberg and if ¦Òi is exactly equal to an eigenvalue, then Ak . ¦ÒkI is singular, and forming the QR-factorization will detect that Rk has some diagonal entry equal to 0. Assuming that the QR-algorithm returns (Rk)nn = 0 (if not, the argument is easily adapted), then the last row of RkQk is 0, so the last row of Ak+1 = RkQk + ¦ÒkI ends with ¦Òk (all other entries being zero), so we are in the case where we can de.ate Ak (and ¦Òk is indeed an eigenvalue). 
The question remains, what is a good choice for the shift ¦Òk? 
Assuming again that H is in upper Hessenberg form, it turns out that when (Hk)nn.1 is small enough, then a good choice for ¦Òk is (Hk)nn. In fact, the rate of convergence is quadratic, which means roughly that the number of correct digits doubles at every iteration. The reason is that shifts are related to another method known as inverse iteration, and such a method converges very fast. For further explanations about this connection, see Demmel 
[16] (Section 4.4.4) and Trefethen and Bau [68] (Lecture 29). 
One should still be cautious that the QR method with shifts does not necessarily converge, and that our convergence proof no longer applies, because instead of having the identity Ak = PkRk, we have 
(A . ¦ÒkI) ¡¤¡¤¡¤ (A . ¦Ò2I)(A . ¦Ò1I)= PkRk. 

17.3. MAKING THE QR METHOD MORE EFFICIENT USING SHIFTS 
Of course, the QR algorithm loops immediately when applied to an orthogonal matrix 
A. This is also the case when A is symmetric but not positive de.nite. For example, both the QR algorithm and the QR algorithm with shifts loop on the matrix 
0  1  
A =  .  
1  0  

In the case of symmetric matrices, Wilkinson invented a shift which helps the QR algo-rithm with shifts to make progress. Again, looking at the lower corner of Ak, say 
an.1 bn.1
B = ,
bn.1 an 
the Wilkinson shift picks the eigenvalue of B closer to an. If we let 
an.1 . an
¦Ä = ,
2 
it is easy to see that the eigenvalues of B are given by 
 
an + an.1
¦Ë = ¡À¦Ä2 + b2 n.1. 
2 
It follows that 
 
¦Ë . an = ¦Ä ¡À¦Ä2 + bn2 .1, 
and from this it is easy to see that the eigenvalue closer to an is given by 
sign(¦Ä)b2 n.1 
¦Ì = an .  . (|¦Ä| +¦Ä2 + b2 )
n.1
If ¦Ä = 0, then we pick arbitrarily one of the two eigenvalues. Observe that the Wilkinson shift applied to the matrix 
01 
A = 
10 
is either +1 or .1, and in one step, de.ation occurs and the algorithm terminates successfully. 
We now discuss double shifts, which are intended to deal with pairs of complex conjugate eigenvalues. 
Let us assume that A is a real matrix. For any complex number ¦Òk with nonzero imaginary part, a double shift consists of the following steps: 
Ak . ¦ÒkI = QkRk 
Ak+1 = RkQk + ¦ÒkI 
Ak+1 . ¦ÒkI = Qk+1Rk+1 
Ak+2 = Rk+1Qk+1 + ¦ÒkI. 
From the computation made for a single shift, we have Ak+1 = Q. kAkQk and Ak+2 = Q. k+1Ak+1Qk+1, so we obtain 
= Q . Q . 
Ak+2 k+1kAkQkQk+1. 
The matrices Qk are complex, so we would expect that the Ak are also complex, but remarkably we can keep the products QkQk+1 real, and so the Ak also real. This is highly desirable to avoid complex arithmetic, which is more expensive. 
Observe that since 
Qk+1Rk+1 = Ak+1 . ¦ÒkI = RkQk +(¦Òk . ¦Òk)I, 
we have 
QkQk+1Rk+1Rk = Qk(RkQk +(¦Òk . ¦Òk)I)Rk 
= QkRkQkRk +(¦Òk . ¦Òk)QkRk 
=(Ak . ¦ÒkI)2 +(¦Òk . ¦Òk)(Ak . ¦ÒkI) 
= A2 . 2(Ìá¦Òk)Ak + |¦Òk|2I. 
k 
If we assume by induction that matrix Ak is real (with k =2M+1,M ¡Ý 0), then the matrix S = A2 . 2(Ìá¦Òk)Ak + |¦Òk|2I is also real, and since QkQk+1 is unitary and Rk+1Rk is upper 
k 
triangular, we see that 
S = QkQk+1Rk+1Rk 
is a QR-factorization of the real matrix S, thus QkQk+1 and Rk+1Rk can be chosen to be real matrices, in which case (QkQk+1). is also real, and thus 
= Q . Q . 
Ak+2 k+1kAkQkQk+1 =(QkQk+1) . AkQkQk+1 
is real. Consequently, if A1 = A is real, then A2´¨+1 is real for all M ¡Ý 0. The strategy that consists in picking ¦Òk and ¦Òk as the complex conjugate eigenvalues of the corner block (Hk)n.1n.1 (Hk)n.1n 
(Hk)nn.1 (Hk)nn 
is called the Francis shift (here we are assuming that A has be reduced to upper Hessenberg form). 
It should be noted that there are matrices for which neither a shift by (Hk)nn nor the Francis shift works. For instance, the permutation matrix 
.  .  
0  0  1  
A = .1  0  0.  

010 
i2¦Ð/3i4¦Ð/3
has eigenvalues e,e, +1, and neither of the above shifts apply to the matrix 00 
. 
10 


17.4. KRYLOV SUBSPACES; ARNOLDI ITERATION 
However, a shift by 1 does work. There are other kinds of matrices for which the QR algorithm does not converge. Demmel gives the example of matrices of the form 
.
. 
... 

0  1  0  0  
1  0  h  0  
0  .h  0  1  
0  0  1  0  

... 

where h is small. 
Algorithms implementing the QR algorithm with shifts and double shifts perform ¡°ex-ceptional¡± shifts every 10 shifts. Despite the fact that the QR algorithm has been perfected since the 1960¡¯s, it is still an open problem to .nd a shift strategy that ensures convergence of all matrices. 
Implicit shifting is based on a result known as the implicit Q theorem. This theorem says that if A is reduced to upper Hessenberg form as A = UHU. and if H is unreduced (hi+1i = 0 for i =1,...,n.1), then the columns of index 2,...,n of U are determined by the .rst column of U up to sign; see Demmel [16] (Theorem 4.9) and Golub and Van Loan [30] (Theorem 7.4.2) for the proof in the case of real matrices. Actually, the proof is not di.cult and will be the object of a homework exercise. In the case of a single shift, an implicit shift generates Ak+1 = Qk.AkQk without having to compute a QR-factorization of Ak . ¦ÒkI. For real matrices, this is done by applying a sequence of Givens rotations which perform a bulge chasing process (a Givens rotation is an orthogonal block diagonal matrix consisting of a single block which is a 2D rotation, the other diagonal entries being equal to 1). Similarly, in the case of a double shift, Ak+2 =(QkQk+1).AkQkQk+1 is generated without having to compute the QR-factorizations of Ak . ¦ÒkI and Ak+1 . ¦ÒkI. Again, (QkQk+1).AkQkQk+1 is generated by applying some simple orthogonal matrices which perform a bulge chasing process. See Demmel [16] (Section 4.4.8) and Golub and Van Loan [30] (Section 7.5) for further explanations regarding implicit shifting involving bulge chasing in the case of real matrices. Watkins [73, 74] discusses bulge chasing in the more general case of complex matrices. 
The Matlab function for .nding the eigenvalues and the eigenvectors of a matrix A is eig and is called as [U, D] = eig(A). It is implemented using an optimized version of the QR-algorithm with implicit shifts. 
If the dimension of the matrix A is very large, we can .nd approximations of some of the eigenvalues of A by using a truncated version of the reduction to Hessenberg form due to Arnoldi in general and to Lanczos in the symmetric (or Hermitian) tridiagonal case. 

17.4 Krylov Subspaces; Arnoldi Iteration 
In this section, we denote the dimension of the square real or complex matrix A by m rather than n, to make it easier for the reader to follow Trefethen and Bau exposition [68], which is particularly lucid. 
Suppose that the m ¡Á m matrix A has been reduced to the upper Hessenberg form H, as A = UHU. . For any n ¡Ü m (typically much smaller than m), consider the (n + 1) ¡Á n 
upper left block 

.
. 
h11 h12 h13 ¡¤¡¤¡¤ h1n h21 h22 h23 ¡¤¡¤¡¤ h2n 0 h32 h33 ¡¤¡¤¡¤ h3n 
........ 

........ 

Hn 
. 
of H, and the n ¡Á n upper Hessenberg matrix Hn obtained by deleting the last row of HQn, 
Q
= 

..
...
.....
.. .
.. 0 ¡¤¡¤¡¤ 0 hnn.1 hnn 0 ¡¤¡¤¡¤ 00 hn+1n 
. 

h11 h12 h13 ¡¤¡¤¡¤ h1n h21 h22 h23 ¡¤¡¤¡¤ h2n 0 h32 h33 ¡¤¡¤¡¤ h3n
Hn = 
...... 

...... 

. 

..
...
.....
.. .
.. 
0 ¡¤¡¤¡¤ 0 hnn.1 hnn 
Q
If we denote by Un the m¡Án matrix consisting of the .rst n columns of U, denoted u1,...,un, then matrix consisting of the .rst n columns of the matrix UH = AU can be expressed as Hn. It follows that the nth column of this matrix can be expressed as 
AUn = Un+1
(. 1) 
Aun = h1nu1 + ¡¤ ¡¤ ¡¤ + hnnun + hn+1nun+1.  (. 2)  
Since (u1, . . . , un) form an orthonormal basis, we deduce from (. 2) that  
(uj, Aun£© = u . j Aun = hjn,  j = 1, . . . , n.  (. 3)  

Q
Hn 
following algorithm due to Arnoldi, known as Arnoldi iteration: Given an arbitrary nonzero vector b ¡Ê Cm, let u1 = b/ lbl; for n =1, 2, 3,... do z := Aun; for j =1 to n do hjn := uj . z; 
z := z . hjnuj 
endfor 
hn+1n := lzl; 
if hn+1n =0 quit 
Equations (. 2) and (. 3) show that Un+1 and
can be computed iteratively using the 


17.4. KRYLOV SUBSPACES; ARNOLDI ITERATION 
un+1 = z/hn+1n When hn+1n = 0, we say that we have a breakdown of the Arnoldi iteration. Arnoldi iteration is an algorithm for producing the n¡Án Hessenberg submatrix Hn of the full Hessenberg matrix H consisting of its .rst n rows and n columns (the .rst n columns of U are also produced), not using Householder matrices. As long as hj+1j = 0 for j =1,...,n, Equation (. 2) shows by an easy induction that un+1 belong to the span of (b, Ab, . . . , Anb), and obviously Aun belongs to the span of (u1,...,un+1), and thus the following spaces are identical: 
Span(b, Ab, . . . , Anb) = Span(u1,...,un+1). 
The space Kn(A, b) = Span(b, Ab, . . . , An.1b) is called a Krylov subspace. We can view Arnoldi¡¯s algorithm as the construction of an orthonormal basis for Kn(A, b). Itisasortof Gram¨CSchmidt procedure. 
Equation (. 2) shows that if Kn is the m ¡Á n matrix whose columns are the vectors (b, Ab, . . . , An.1b), then there is a n ¡Á n upper triangular matrix Rn such that 
Kn = UnRn. (. 4) 
The above is called a reduced QR factorization of Kn. Since (u1,...,un) is an orthonormal system, the matrix Un.Un+1 is the n ¡Á (n + 1) matrix 
Q
consisting of the identity matrix In plus an extra column of 0¡¯s, so Un.Un+1Hn = Un.AUn is obtained by deleting the last row of HQn, namely Hn, and so 
U . nAUn = Hn. (. 5) 
We summarize the above facts in the following proposition. 
Proposition 17.5. If Arnoldi iteration run on an m ¡Á m matrix A starting with a nonzero vector b ¡Ê Cm does not have a breakdown at stage n ¡Ü m, then the following properties hold: 
(1) If Kn is the m ¡Á n Krylov matrix associated with the vectors (b, Ab, . . . , An.1b) and if Un is the m ¡Á n matrix of orthogonal vectors produced by Arnoldi iteration, then there is a QR-factorization 
Kn = UnRn, 
for some n ¡Á n upper triangular matrix Rn. 
(2) 
The m ¡Á n upper Hessenberg matrices Hn produced by Arnoldi iteration are the projec-tion of A onto the Krylov space Kn(A, b), that is, 

Hn = Un. AUn. 

(3) 
The successive iterates are related by the formula 


Q
AUn = Un+1Hn. 
Remark: If Arnoldi iteration has a breakdown at stage n, that is, hn+1 = 0, then we found the .rst unreduced block of the Hessenberg matrix H. It can be shown that the eigenvalues of Hn are eigenvalues of A. So a breakdown is actually a good thing. In this case, we can pick some new nonzero vector un+1 orthogonal to the vectors (u1,...,un) as a new starting vector and run Arnoldi iteration again. Such a vector exists since the (n + 1)th column of U works. So repeated application of Arnoldi yields a full Hessenberg reduction of A. However, this is not what we are after, since m is very large an we are only interested in a ¡°small¡± number of eigenvalues of A. 
There is another aspect of Arnoldi iteration, which is that it solves an optimization problem involving polynomials of degree n. Let Pn denote the set of (complex) monic polynomials of degree n, that is, polynomials of the form 
nn.1 
p(z)= z + cn.1z + ¡¤¡¤¡¤ + c1z + c0 (ci ¡Ê C). 
For any m ¡Á m matrix A, we write 
p(A)= An + cn.1An.1 + ¡¤¡¤¡¤ + c1A + c0I. 
The following result is proven in Trefethen and Bau [68] (Lecture 34, Theorem 34.1). 
Theorem 17.6. If Arnoldi iteration run on an m ¡Á m matrix A starting with a nonzero vector b does not have a breakdown at stage n ¡Ü m, then there is a unique polynomial p ¡ÊPn such that lp(A)bl2 is minimum, namely the characteristic polynomial det(zI . Hn) of Hn. 
Theorem 17.6 can be viewed as the ¡°justi.cation¡± for a method to .nd some of the eigenvalues of A (say n¡¶ m of them). Intuitively, the closer the roots of the character-istic polynomials of Hn are to the eigenvalues of A, the smaller lp(A)bl2 should be, and conversely. In the extreme case where m = n, by the Cayley¨CHamilton theorem, p(A)=0 (where p is the characteristic polynomial of A), so this idea is plausible, but this is far from constituting a proof (also, b should have nonzero coordinates in all directions associated with the eigenvalues). 
The method known as the Rayleigh¨CRitz method is to run Arnoldi iteration on A and some b = 0 chosen at random for n¡¶ m steps before or until a breakdown occurs. Then run the QR algorithm with shifts on Hn. The eigenvalues of the Hessenberg matrix Hn may then be considered as approximations of the eigenvalues of A. The eigenvalues of Hn are called Arnoldi estimates or Ritz values. One has to be cautious because Hn is a truncated version of the full Hessenberg matrix H, so not all of the Ritz values are necessary close to eigenvalues of A. It has been observed that the eigenvalues that are found .rst are the extreme eigenvalues of A, namely those close to the boundary of the spectrum of A plotted in 
C. So if A has real eigenvalues, the largest and the smallest eigenvalues appear .rst as Ritz values. In many problems where eigenvalues occur, the extreme eigenvalues are the one that need to be computed. Similarly, the eigenvectors of Hn may be considered as approximations of eigenvectors of A. 
17.5. GMRES 
The Matlab function eigs is based on the computation of Ritz values. It computes the six eigenvalues of largest magnitude of a matrix A, and the call is [V, D] = eigs(A). More generally, to get the top k eigenvalues, use [V, D] = eigs(A, k). 
In the absence of rigorous theorems about error estimates, it is hard to make the above statements more precise; see Trefethen and Bau [68] (Lecture 34) for more on this subject. 
However, if A is a symmetric (or Hermitian) matrix, then Hn is a symmetric (resp. Hermitian) tridiagonal matrix and more precise results can be shown; see Demmel [16] (Chapter 7, especially Section 7.2). We will consider the symmetric (and Hermitan) case in the next section, but .rst we show how Arnoldi iteration can be used to .nd approximations for the solution of a linear system Ax = b where A is invertible but of very large dimension 
m. 
17.5 GMRES 
Suppose A is an invertible m¡Ám matrix and let b be a nonzero vector in Cm . Let x0 = A.1b, the unique solution of Ax = b. It is not hard to show that x0 ¡ÊKn(A, b) for some n ¡Ü m. In fact, there is a unique monic polynomial p(z) of minimal degree s ¡Ü m such that p(A)b = 0, so x0 ¡ÊKs(A, b). Thus it makes sense to search for a solution of Ax = b in Krylov spaces of dimension m ¡Ü s. The idea is to .nd an approximation xn ¡ÊKn(A, b) of x0 such that rn = b . Axn is minimized, that is, lrnl2 = lb . Axnl2 is minimized over xn ¡ÊKn(A, b). This minimization problem can be stated as 
minimize lrnl2 = lAxn . bl2 ,xn ¡ÊKn(A, b). 
This is a least-squares problem, and we know how to solve it (see Section 21.1). The quantity rn is known as the residual and the method which consists in minimizing lrnl2 is known as GMRES, for generalized minimal residuals. 
Now since (u1,...,un) is a basis of Kn(A, b) (since n ¡Ü s, no breakdown occurs, except for n = s), we may write xn = Uny, so our minimization problem is 
minimize lAUny . bl2 ,y ¡Ê Cn . 
Q
Since by (. 1) of Section 17.4, we have AUn = Un+1Hn, minimizing lAUny . bl2 is equiv-
QQ
alent to minimizing lUn+1Hny . bl2 over Cm . Since Un+1Hny and b belong to the column space of Un+1, minimizing lUn+1HQny . bl2 is equivalent to minimizing lHQny . Un. +1bl2. However, by construction, 
b = lbl2e1 ¡Ê Cn+1
Un. +1, 
so our minimization problem can be stated as 
minimize lHQny .lbl2e1l2,y ¡Ê Cn . 
The approximate solution of Ax = b is then 
xn = Uny. 
Starting with u1 = b/ lbl2 and with n = 1, the GMRES method runs n ¡Ü s Arnoldi Hn
minimize
Q
iterations to .nd Un and
, and then runs a method to solve the least squares problem 

lHQny .lbl2e1l2,y ¡Ê Cn 
. 

When lrnl= lHQny .lbl2e1l2 is considered small enough, we stop and the approximate 
2 
solution of Ax = b is then 
xn = Uny. 
There are ways of improving e.ciency of the ¡°naive¡± version of GMRES that we just presented; see Trefethen and Bau [68] (Lecture 35). We now consider the case where A is a Hermitian (or symmetric) matrix. 
17.6 The Hermitian Case; Lanczos Iteration 
If A is an m¡Ám symmetric or Hermitian matrix, then Arnoldi¡¯s method is simpler and much more e.cient. Indeed, in this case, it is easy to see that the upper Hessenberg matrices Hn are also symmetric (Hermitian respectively), and thus tridiagonal. Also, the eigenvalues of A and Hn are real. It is convenient to write 
.
. 
¦Á1 ¦Â1 ¦Â1 ¦Á2 ¦Â2 
..
Hn 
= 

...... 

...... 

.

¦Â2 ¦Á3 
. 

..
..
.. 
¦Ân.1 ¦Ân.1 ¦Án 
The recurrence (. 2) of Section 17.4 becomes the three-term recurrence 
Aun = ¦Ân.1un.1 + ¦Ánun + ¦Ânun+1. (. 6) 
We also have ¦Án = un. AUn, so Arnoldi¡¯s algorithm become the following algorithm known as Lanczos¡¯ algorithm (or Lanczos iteration). The inner loop on j from 1 to n has been eliminated and replaced by a single assignment. 
Given an arbitrary nonzero vector b ¡Ê Cm, let u1 = b/ lbl; for n =1, 2, 3,... do 
z := Aun; 
¦Án := u . nz; 
z := z . ¦Ân.1un.1 . ¦Ánun 

¦Ân := lzl; 
if ¦Ân =0 quit 
un+1 = z/¦Ân 
17.7. POWER METHODS 
When ¦Ân = 0, we say that we have a breakdown of the Lanczos iteration. 
Versions of Proposition 17.5 and Theorem 17.6 apply to Lanczos iteration. 
Besides being much more e.cient than Arnoldi iteration, Lanczos iteration has the advan-tage that the Rayleigh¨CRitz method for .nding some of the eigenvalues of A as the eigenvalues of the symmetric (respectively Hermitian) tridiagonal matrix Hn applies, but there are more methods for .nding the eigenvalues of symmetric (respectively Hermitian) tridiagonal matri-ces. Also theorems about error estimates exist. The version of Lanczos iteration given above may run into problems in .oating point arithmetic. What happens is that the vectors uj may lose the property of being orthogonal, so it may be necessary to reorthogonalize them. For more on all this, see Demmel [16] (Chapter 7, in particular Section 7.2-7.4). The version of GMRES using Lanczos iteration is called MINRES. 
We close our brief survey of methods for computing the eigenvalues and the eigenvectors of a matrix with a quick discussion of two methods known as power methods. 
17.7 Power Methods 
Let A be an m ¡Á m complex or real matrix. There are two power methods, both of which yield one eigenvalue and one eigenvector associated with this vector: 
(1) 
Power iteration. 

(2) 
Inverse (power) iteration. 


Power iteration only works if the matrix A has an eigenvalue ¦Ë of largest modulus, which means that if ¦Ë1,...,¦Ëm are the eigenvalues of A, then 
|¦Ë1| > |¦Ë2|¡Ý ¡¤ ¡¤¡¤¡Ý |¦Ëm|¡Ý 0. 
In particular, if A is a real matrix, then ¦Ë1 must be real (since otherwise there are two complex conjugate eigenvalues of the same largest modulus). If the above condition is satis.ed, then power iteration yields ¦Ë1 and some eigenvector associated with it. The method is simple enough: 
Pick some initial unit vector x0 and compute the following sequence (xk), where 
k+1 Axk 
x = ,k ¡Ý 0. 
lAxkl
We would expect that (xk) converges to an eigenvector associated with ¦Ë1, but this is not quite correct. The following results are proven in Serre [57] (Section 13.5.1). First assume that ¦Ë1 = 0. 
We have lim Axk = |¦Ë1|. 
kh¡ú¡Þ 
If A is a complex matrix which has a unique complex eigenvalue ¦Ë1 of largest modulus, then 
k ¦Ë1 k 
v = lim x 
kh¡ú¡Þ |¦Ë1| 
is a unit eigenvector of A associated with ¦Ë1. If ¦Ë1 is real, then 
k 
v = lim x 
kh¡ú¡Þ 
is a unit eigenvector of A associated with ¦Ë1. Actually some condition on x0 is needed: x0 must have a nonzero component in the eigenspace E associated with ¦Ë1 (in any direct sum of Cm in which E is a summand). 
The eigenvalue ¦Ë1 is found as follows. If ¦Ë1 is complex, and if vj = 0 is any nonzero coordinate of v, then (Axk)j
¦Ë1 = lim k . 
kh¡ú¡Þ x
j 
If ¦Ë1 is real, then we can de.ne the sequence (¦Ë(k)) by 
¦Ë(k+1) =(x k+1) . Axk+1 ,k ¡Ý 0, 
and we have 
¦Ë(k)
¦Ë1 = lim . 
kh¡ú¡Þ 
Indeed, in this case, since v = limkh¡ú¡Þ xk and v is a unit eigenvector for ¦Ë1, we have 
¦Ë(k) k+1) . Axk+1 . 
lim = lim(x = v . Av = ¦Ë1vv = ¦Ë1. 
kh¡ú¡Þ kh¡ú¡Þ
Note that since xk+1 is a unit vector, (xk+1).Axk+1 is a Rayleigh ratio. 
If A is a Hermitian matrix, then the eigenvalues are real and we can say more about the rate of convergence, which is not great (only linear). For details, see Trefethen and Bau [68] (Lecture 27). 
If ¦Ë1 = 0, then there is some power M<m such that Ax´¨ 
= 0. The inverse iteration method is designed to .nd an eigenvector associated with an eigen-value ¦Ë of A for which we know a good approximation ¦Ì. Pick some initial unit vector x0 and compute the following sequences (wk) and (xk), where wk+1 is the solution of the system 
k+1 kk
(A . ¦ÌI)w = x equivalently w k+1 =(A . ¦ÌI).1 x, k ¡Ý 0, 
and 
k+1
w
x k+1 = ,k ¡Ý 0. 
k+1l
lw
The following result is proven in Ciarlet [14] (Theorem 6.4.1). 
17.8. SUMMARY 
Proposition 17.7. Let A be an m ¡Á m diagonalizable (complex or real) matrix with eigen-values ¦Ë1,...,¦Ëm, and let ¦Ë = ¦Ë´¨ be an arbitrary eigenvalue of A (not necessary simple). For any ¦Ì such that 
¦Ì = ¦Ë and |¦Ì . ¦Ë| < |¦Ì . ¦Ëj| for all j = M, 
if x0 does not belong to the subspace spanned by the eigenvectors associated with the eigen-values ¦Ëj with j = M, then 
(¦Ë . ¦Ì)kk
lim x = v, 
kh¡ú¡Þ |¦Ë . ¦Ì|k 
where v is an eigenvector associated with ¦Ë. Furthermore, if both ¦Ë and ¦Ì are real, we have 
lim x k = v if ¦Ì < ¦Ë, 
kh¡ú¡Þ 
lim (.1)k x k = v if ¦Ì > ¦Ë. 
kh¡ú¡Þ
Also, if we de.ne the sequence (¦Ë(k)) by 
¦Ë(k+1) =(x k+1) . Axk+1 
, 
then 
¦Ë(k+1)
lim = ¦Ë. 
kh¡ú¡Þ 
The condition of x0 may seem quite stringent, but in practice, a vector x0 chosen at random usually satis.es it. 
If A is a Hermitian matrix, then we can say more. In particular, the inverse iteration algorithm can be modi.ed to make use of the newly computed ¦Ë(k+1) instead of ¦Ì, and an even faster convergence is achieved. Such a method is called the Rayleigh quotient iteration. When it converges (which is for almost all x0), this method eventually achieves cubic convergence, which is remarkable. Essentially, this means that the number of correct digits is tripled at every iteration. For more details, see Trefethen and Bau [68] (Lecture 27) and Demmel [16] (Section 5.3.2). 
17.8 Summary 
The main concepts and results of this chapter are listed below: 
. QR iteration, QR algorithm. 

. Upper Hessenberg matrices. 

. Householder matrix. 

. Unreduced and reduced Hessenberg matrices. 

. 
De.ation. 

. 
Shift. 

. 
Wilkinson shift. 

. 
Double shift. 

. 
Francis shift. 

. 
Implicit shifting. 

. 
Implicit Q-theorem. 

. 
Arnoldi iteration. 

. 
Breakdown of Arnoldi iteration. 

. 
Krylov subspace. 

. 
Rayleigh¨CRitz method. 

. 
Ritz values, Arnoldi estimates. 

. 
Residual. 

. 
GMRES 

. 
Lanczos iteration. 

. 
Power iteration. 

. 
Inverse power iteration. 

. 
Rayleigh ratio. 


17.9 Problems 
Problem 17.1. Prove Theorem 17.2; see Problem 12.7. 
Problem 17.2. Prove that if a matrix A is Hermitian (or real symmetric), then any Hes-senberg matrix H similar to A is Hermitian tridiagonal (real symmetric tridiagonal). 
Problem 17.3. For any matrix (real or complex) A, if A = QR is a QR-decomposition of A using Householder re.ections, prove that if A is upper Hessenberg then so is Q. 
Problem 17.4. Prove that if A is upper Hessenberg, then the matrices Ak obtained by applying the QR-algorithm are also upper Hessenberg. 
17.9. PROBLEMS 
Problem 17.5. Prove the implicit Q theorem. This theorem says that if A is reduced to upper Hessenberg form as A = UHU. and if H is unreduced (hi+1i = 0 for i =1,...,n . 1), then the columns of index 2,...,n of U are determined by the .rst column of U up to sign; 
Problem 17.6. Read Section 7.5 of Golub and Van Loan [30] and implement their version of the QR-algorithm with shifts. 
Problem 17.7. If an Arnoldi iteration has a breakdown at stage n, that is, hn+1 = 0, then we found the .rst unreduced block of the Hessenberg matrix H. Prove that the eigenvalues of Hn are eigenvalues of A. 
Problem 17.8. Prove Theorem 17.6. 
Problem 17.9. Implement GRMES and test it on some linear systems. 
Problem 17.10. State and prove versions of Proposition 17.5 and Theorem 17.6 for the Lanczos iteration. 
Problem 17.11. Prove the results about the power iteration method stated in Section 17.7. 
Problem 17.12. Prove the results about the inverse power iteration method stated in Section 17.7. 
Problem 17.13. Implement and test the power iteration method and the inverse power iteration method. 
Problem 17.14. Read Lecture 27 in Trefethen and Bau [68] and implement and test the Rayleigh quotient iteration method. 


