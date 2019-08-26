Chapter 14 
Eigenvectors and Eigenvalues 
In this chapter all vector spaces are de.ned over an arbitrary .eld K. For the sake of concreteness, the reader may safely assume that K = R or K = C. 
14.1 Eigenvectors and Eigenvalues of a Linear Map 
Given a .nite-dimensional vector space E, let f : E ¡ú E be any linear map. If by luck there is a basis (e1,...,en) of E with respect to which f is represented by a diagonal matrix 
.
. 
D = 

.... 

¦Ë1 0 	... 0 .
. 
..
0 ¦Ë2 .. 
. 
..
..	.
. ..0 
0 ... 0 ¦Ën 
....

, 

then the action of f on E is very simple; in every ¡°direction¡± ei, we have 
f(ei)= ¦Ëiei. 
We can think of f as a transformation that stretches or shrinks space along the direction e1,...,en (at least if E is a real vector space). In terms of matrices, the above property translates into the fact that there is an invertible matrix P and a diagonal matrix D such that a matrix A can be factored as 
A = P DP .1 . 
When this happens, we say that f (or A) is diagonalizable, the ¦Ëi¡¯s are called the eigenvalues of f, and the ei¡¯s are eigenvectors of f. For example, we will see that every symmetric matrix can be diagonalized. Unfortunately, not every matrix can be diagonalized. For example, the 
matrix

 
 

A1 = 1 0  1 1 
489  

can¡¯t be diagonalized. Sometimes a matrix fails to be diagonalizable because its eigenvalues do not belong to the .eld of coe.cients, such as 
0 .1 
A2 = ,
10 
whose eigenvalues are ¡Ài. This is not a serious problem because A2 can be diagonalized over the complex numbers. However, A1 is a ¡°fatal¡± case! Indeed, its eigenvalues are both 1 and the problem is that A1 does not have enough eigenvectors to span E. 
The next best thing is that there is a basis with respect to which f is represented by an upper triangular matrix. In this case we say that f can be triangularized, or that f is triangulable. As we will see in Section 14.2, if all the eigenvalues of f belong to the .eld of coe.cients K, then f can be triangularized. In particular, this is the case if K = C. 
Now an alternative to triangularization is to consider the representation of f with respect to two bases (e1,...,en) and (f1,...,fn), rather than a single basis. In this case, if K = R or K = C, it turns out that we can even pick these bases to be orthonormal, and we get a diagonal matrix ¦² with nonnegative entries, such that 
f(ei)= ¦Òifi, 1 ¡Ü i ¡Ü n. 
The nonzero ¦Òi¡¯s are the singular values of f, and the corresponding representation is the singular value decomposition, or SVD. The SVD plays a very important role in applications, and will be considered in detail in Chapter 20. 
In this section we focus on the possibility of diagonalizing a linear map, and we introduce the relevant concepts to do so. Given a vector space E over a .eld K, let id denote the identity map on E. 
The notion of eigenvalue of a linear map f : E ¡ú E de.ned on an in.nite-dimensional space E is quite subtle because it cannot be de.ned in terms of eigenvectors as in the .nite-dimensional case. The problem is that the map ¦Ë id . f (with ¦Ë ¡Ê C) could be noninvertible (because it is not surjective) and yet injective. In .nite dimension this cannot happen, so until further notice we assume that E is of .nite dimension n. 
De.nition 14.1. Given any vector space E of .nite dimension n and any linear map f : E ¡ú E, a scalar ¦Ë ¡Ê K is called an eigenvalue, or proper value, or characteristic value of f if there is some nonzero vector u ¡Ê E such that 
f(u)= ¦Ëu. 
Equivalently, ¦Ë is an eigenvalue of f if Ker(¦Ë id . f) is nontrivial (i.e., Ker (¦Ë id . f) 
= {0}) i. ¦Ë id.f is not invertible (this is where the fact that E is .nite-dimensional is used; a linear map from E to itself is injective i. it is invertible). A vector u ¡Ê E is called an eigenvector, or proper vector, or characteristic vector of f if u there is some ¦Ë ¡Ê K such that 
= 0 and if 
f(u)= ¦Ëu; 

14.1. EIGENVECTORS AND EIGENVALUES OF A LINEAR MAP 
the scalar ¦Ë is then an eigenvalue, and we say that u is an eigenvector associated with ¦Ë. Given any eigenvalue ¦Ë ¡Ê K, the nontrivial subspace Ker (¦Ë id . f) consists of all the eigenvectors associated with ¦Ë together with the zero vector; this subspace is denoted by E¦Ë(f), or E(¦Ë, f), or even by E¦Ë, and is called the eigenspace associated with ¦Ë, or proper subspace associated with ¦Ë. 
Note that distinct eigenvectors may correspond to the same eigenvalue, but distinct eigenvalues correspond to disjoint sets of eigenvectors. 
Remark: As we emphasized in the remark following De.nition 8.4, we require an eigenvector to be nonzero. This requirement seems to have more bene.ts than inconveniences, even though it may considered somewhat inelegant because the set of all eigenvectors associated with an eigenvalue is not a subspace since the zero vector is excluded. 
The next proposition shows that the eigenvalues of a linear map f : E ¡ú E are the roots of a polynomial associated with f. 
Proposition 14.1. Let E be any vector space of .nite dimension n and let f be any linear map f : E ¡ú E. The eigenvalues of f are the roots (in K) of the polynomial 
det(¦Ë id . f). 
Proof. A scalar ¦Ë ¡Ê K is an eigenvalue of f i. there is some vector u =0 in E such that f(u)= ¦Ëu i. (¦Ë id . f)(u)=0 i. (¦Ë id . f) is not invertible i., by Proposition 6.14, det(¦Ë id . f)=0. 
In view of the importance of the polynomial det(¦Ë id.f), we have the following de.nition. 
De.nition 14.2. Given any vector space E of dimension n, for any linear map f : E ¡ú E, the polynomial Pf (X)= ¦Öf (X) = det(X id . f) is called the characteristic polynomial of 
f. For any square matrix A, the polynomial PA(X)= ¦ÖA(X) = det(XI . A) is called the characteristic polynomial of A. 
Note that we already encountered the characteristic polynomial in Section 6.7; see De.-nition 6.11. 
Given any basis (e1,...,en), if A = M(f) is the matrix of f w.r.t. (e1,...,en), we can compute the characteristic polynomial ¦Öf (X) = det(X id . f) of f by expanding the following determinant: 
det(XI . A)=

         

X . a11 .a12 ... .a1 n .a21 X . a22 ... .a2 n 
.. .
.
. ...
.
.. . 
         

. 

.an 1 .an 2 ... X . ann
If we expand this determinant, we .nd that 
¦ÖA(X) = det(XI . A)= Xn . (a11 + ¡¤¡¤¡¤ + ann)Xn.1 + ¡¤¡¤¡¤ +(.1)n det(A). 
The sum tr(A)= a11 + ¡¤¡¤¡¤ + ann of the diagonal elements of A is called the trace of A. Since we proved in Section 6.7 that the characteristic polynomial only depends on the linear map f, the above shows that tr(A) has the same value for all matrices A representing f. Thus, the trace of a linear map is well-de.ned; we have tr(f) = tr(A) for any matrix A representing 
f. 
Remark: The characteristic polynomial of a linear map is sometimes de.ned as det(f . X id). Since 
det(f . X id) = (.1)n det(X id . f), 
this makes essentially no di.erence but the version det(X id . f) has the small advantage that the coe.cient of Xn is +1. 
If we write 
¦ÖA(X) = det(XI . A) 
= Xn . ¦Ó1(A)Xn.1 + ¡¤¡¤¡¤ +(.1)k¦Ók(A)Xn.k + ¡¤¡¤¡¤ +(.1)n¦Ón(A), 
then we just proved that 
¦Ó1(A) = tr(A) and ¦Ón(A) = det(A). 
It is also possible to express ¦Ók(A) in terms of determinants of certain submatrices of A. For any nonempty subset, I .{1,...,n}, say I = {i1 < ... < ik}, let AI,I be the k¡Ák submatrix of A whose jth column consists of the elements aih ij , where h =1,...,k. Equivalently, AI,I is the matrix obtained from A by .rst selecting the columns whose indices belong to I, and then the rows whose indices also belong to I. Then it can be shown that 
 
¦Ók(A) =det(AI,I ). 
I.{1,...,n}|I|=k 
If all the roots, ¦Ë1,...,¦Ën, of the polynomial det(XI . A) belong to the .eld K, then we can write ¦ÖA(X) = det(XI . A)=(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ën), 
where some of the ¦Ëi¡¯s may appear more than once. Consequently, 
¦ÖA(X) = det(XI . A) 
= Xn . ¦Ò1(¦Ë)Xn.1 + ¡¤¡¤¡¤ +(.1)k¦Òk(¦Ë)Xn.k + ¡¤¡¤¡¤ +(.1)n¦Òn(¦Ë), 

14.1. EIGENVECTORS AND EIGENVALUES OF A LINEAR MAP 
where 
 
¦Òk(¦Ë)= ¦Ëi, 
I.{1,...,n} i¡ÊI |I|=k
the kth elementary symmetric polynomial (or function) of the ¦Ëi¡¯s, where ¦Ë =(¦Ë1,...,¦Ën). The elementary symmetric polynomial ¦Òk(¦Ë) is often denoted Ek(¦Ë), but this notation may be confusing in the context of linear algebra. For n = 5, the elementary symmetric polynomials are listed below: 
¦Ò0(¦Ë)=1 
¦Ò1(¦Ë)= ¦Ë1 + ¦Ë2 + ¦Ë3 + ¦Ë4 + ¦Ë5 
¦Ò2(¦Ë)= ¦Ë1¦Ë2 + ¦Ë1¦Ë3 + ¦Ë1¦Ë4 + ¦Ë1¦Ë5 + ¦Ë2¦Ë3 + ¦Ë2¦Ë4 + ¦Ë2¦Ë5 
+ 
¦Ë3¦Ë4 + ¦Ë3¦Ë5 + ¦Ë4¦Ë5 ¦Ò3(¦Ë)= ¦Ë3¦Ë4¦Ë5 + ¦Ë2¦Ë4¦Ë5 + ¦Ë2¦Ë3¦Ë5 + ¦Ë2¦Ë3¦Ë4 + ¦Ë1¦Ë4¦Ë5 

+ 
¦Ë1¦Ë3¦Ë5 + ¦Ë1¦Ë3¦Ë4 + ¦Ë1¦Ë2¦Ë5 + ¦Ë1¦Ë2¦Ë4 + ¦Ë1¦Ë2¦Ë3 ¦Ò4(¦Ë)= ¦Ë1¦Ë2¦Ë3¦Ë4 + ¦Ë1¦Ë2¦Ë3¦Ë5 + ¦Ë1¦Ë2¦Ë4¦Ë5 + ¦Ë1¦Ë3¦Ë4¦Ë5 + ¦Ë2¦Ë3¦Ë4¦Ë5 ¦Ò5(¦Ë)= ¦Ë1¦Ë2¦Ë3¦Ë4¦Ë5. 


Since 
¦ÖA(X)= Xn . ¦Ó1(A)Xn.1 + ¡¤¡¤¡¤ +(.1)k¦Ók(A)Xn.k + ¡¤¡¤¡¤ +(.1)n¦Ón(A) 
= Xn . ¦Ò1(¦Ë)Xn.1 + ¡¤¡¤¡¤ +(.1)k¦Òk(¦Ë)Xn.k + ¡¤¡¤¡¤ +(.1)n¦Òn(¦Ë), 
we have ¦Òk(¦Ë)= ¦Ók(A),k =1, . . . , n, 
and in particular, the product of the eigenvalues of f is equal to det(A) = det(f), and the sum of the eigenvalues of f is equal to the trace tr(A) = tr(f), of f; for the record, 
tr(f)= ¦Ë1 + ¡¤¡¤¡¤ + ¦Ën 
det(f)= ¦Ë1 ¡¤¡¤¡¤ ¦Ën, 
where ¦Ë1,...,¦Ën are the eigenvalues of f (and A), where some of the ¦Ëi¡¯s may appear more than once. In particular, f is not invertible i. it admits 0 has an eigenvalue (since f is singular i. ¦Ë1 ¡¤¡¤¡¤ ¦Ën = det(f) = 0). 
Remark: Depending on the .eld K, the characteristic polynomial ¦ÖA(X) = det(XI . A) may or may not have roots in K. This motivates considering algebraically closed .elds, which are .elds K such that every polynomial with coe.cients in K has all its root in K. For example, over K = R, not every polynomial has real roots. If we consider the matrix 
cos ¦È . sin ¦È 
A = ,
sin ¦È cos ¦È 
then the characteristic polynomial det(XI . A) has no real roots unless ¦È = k¦Ð. However, over the .eld C of complex numbers, every polynomial has roots. For example, the matrix above has the roots cos ¦È ¡À i sin ¦È = e¡Ài¦È . 
Remark: It is possible to show that every linear map f over a complex vector space E must have some (complex) eigenvalue without having recourse to determinants (and the characteristic polynomial). Let n = dim(E), pick any nonzero vector u ¡Ê E, and consider the sequence 
u, f(u),f2(u),...,fn(u). 
Since the above sequence has n + 1 vectors and E has dimension n, these vectors must be linearly dependent, so there are some complex numbers c0,...,cm, not all zero, such that 
c0fm(u)+ c1fm.1(u)+ ¡¤¡¤¡¤ + cmu =0, 
where m ¡Ü n is the largest integer such that the coe.cient of fm(u) is nonzero (m must exits since we have a nontrivial linear dependency). Now because the .eld C is algebraically closed, the polynomial 
c0Xm + c1Xm.1 + ¡¤¡¤¡¤ + cm 
can be written as a product of linear factors as 
c0Xm + c1Xm.1 + ¡¤¡¤¡¤ + cm = c0(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ëm) 
for some complex numbers ¦Ë1,...,¦Ëm ¡Ê C, not necessarily distinct. But then since c0 = 0, 
c0fm(u)+ c1fm.1(u)+ ¡¤¡¤¡¤ + cmu =0 
is equivalent to (f . ¦Ë1 id) .¡¤ ¡¤¡¤. (f . ¦Ëm id)(u)=0. 
If all the linear maps f . ¦Ëi id were injective, then (f . ¦Ë1 id) .¡¤ ¡¤¡¤. (f . ¦Ëm id) would be injective, contradicting the fact that u = 0. Therefore, some linear map f . ¦Ëi id must have a nontrivial kernel, which means that there is some v = 0 so that 
f(v)= ¦Ëiv; 
that is, ¦Ëi is some eigenvalue of f and v is some eigenvector of f. 
As nice as the above argument is, it does not provide a method for .nding the eigenvalues of f, and even if we prefer avoiding determinants as a much as possible, we are forced to deal with the characteristic polynomial det(X id . f). 
De.nition 14.3. Let A be an n ¡Á n matrix over a .eld K. Assume that all the roots of the characteristic polynomial ¦ÖA(X) = det(XI . A) of A belong to K, which means that we can write 
det(XI . A)=(X . ¦Ë1)k1 ¡¤¡¤¡¤ (X . ¦Ëm)km , 

14.1. EIGENVECTORS AND EIGENVALUES OF A LINEAR MAP 
where ¦Ë1,...,¦Ëm ¡Ê K are the distinct roots of det(XI . A) and k1 + ¡¤¡¤¡¤ + km = n. The integer ki is called the algebraic multiplicity of the eigenvalue ¦Ëi, and the dimension of the eigenspace E¦Ëi = Ker(¦ËiI . A) is called the geometric multiplicity of ¦Ëi. We denote the algebraic multiplicity of ¦Ëi by alg(¦Ëi), and its geometric multiplicity by geo(¦Ëi). 
By de.nition, the sum of the algebraic multiplicities is equal to n, but the sum of the geometric multiplicities can be strictly smaller. 
Proposition 14.2. Let A be an n ¡Á n matrix over a .eld K and assume that all the roots of the characteristic polynomial ¦ÖA(X) = det(XI.A) of A belong to K. For every eigenvalue ¦Ëi of A, the geometric multiplicity of ¦Ëi is always less than or equal to its algebraic multiplicity, that is, 
geo(¦Ëi) ¡Ü alg(¦Ëi). 
Proof. To see this, if ni is the dimension of the eigenspace E¦Ëi associated with the eigenvalue ¦Ëi, we can form a basis of Kn obtained by picking a basis of E¦Ëi and completing this linearly independent family to a basis of Kn . With respect to this new basis, our matrix is of the form 
½Ð¦ËiIni B 
A= ,
0 D and a simple determinant calculation shows that det(XI . A) = det(XI . A½Ð)=(X . ¦Ëi)ni det(XIn.ni . D). 
Therefore, (X .¦Ëi)ni divides the characteristic polynomial of A½Ð, and thus, the characteristic polynomial of A. It follows that ni is less than or equal to the algebraic multiplicity of ¦Ëi. 
The following proposition shows an interesting property of eigenspaces. 
Proposition 14.3. Let E be any vector space of .nite dimension n and let f be any linear map. If u1,...,um are eigenvectors associated with pairwise distinct eigenvalues ¦Ë1,...,¦Ëm, then the family (u1,...,um) is linearly independent. 
Proof. Assume that (u1,...,um) is linearly dependent. Then there exists ¦Ì1,...,¦Ìk ¡Ê K such that 
¦Ì1ui1 + ¡¤¡¤¡¤ + ¦Ìkuik =0, where 1 ¡Ü k ¡Ü m, ¦Ìi = 0 for all i,1 ¡Ü i ¡Ü k, {i1,...,ik}.{1,...,m}, and no proper subfamily of (ui1 ,...,uik ) is linearly dependent (in other words, we consider a dependency relation with k minimal). Applying f to this dependency relation, we get 
¦Ì1¦Ëi1 ui1 + ¡¤¡¤¡¤ + ¦Ìk¦Ëik uik =0, 
and if we multiply the original dependency relation by ¦Ëi1 and subtract it from the above, we get 
¦Ì2(¦Ëi2 . ¦Ëi1 )ui2 + ¡¤¡¤¡¤ + ¦Ìk(¦Ëik . ¦Ëi1 )uik =0, which is a nontrivial linear dependency among a proper subfamily of (ui1 ,...,uik ) since the ¦Ëj are all distinct and the ¦Ìi are nonzero, a contradiction. 
As a corollary of Proposition 14.3 we have the following result. 
Corollary 14.4. If ¦Ë1,...,¦Ëm are all the pairwise distinct eigenvalues of f (where m ¡Ü n), we have a direct sum 
E¦Ë1 ¨’¡¤ ¡¤¡¤¨’ E¦Ëm 
of the eigenspaces E¦Ëi . 
Unfortunately, it is not always the case that 
E = E¦Ë1 ¨’¡¤ ¡¤¡¤¨’ E¦Ëm . 
De.nition 14.4. When 
E = E¦Ë1 ¨’¡¤ ¡¤¡¤¨’ E¦Ëm , 
we say that f is diagonalizable (and similarly for any matrix associated with f). 
Indeed, picking a basis in each E¦Ëi , we obtain a matrix which is a diagonal matrix consisting of the eigenvalues, each ¦Ëi occurring a number of times equal to the dimension of E¦Ëi . This happens if the algebraic multiplicity and the geometric multiplicity of every eigenvalue are equal. In particular, when the characteristic polynomial has n distinct roots, then f is diagonalizable. It can also be shown that symmetric matrices have real eigenvalues and can be diagonalized. 
For a negative example, we leave it as exercise to show that the matrix 
1  1  
M =  
0  1  

cannot be diagonalized, even though 1 is an eigenvalue. The problem is that the eigenspace of 1 only has dimension 1. The matrix 
cos ¦È . sin ¦È 
A = 
sin ¦È cos ¦È 
cannot be diagonalized either, because it has no real eigenvalues, unless ¦È = k¦Ð. However, over the .eld of complex numbers, it can be diagonalized. 


14.2 Reduction to Upper Triangular Form 
Unfortunately, not every linear map on a complex vector space can be diagonalized. The next best thing is to ¡°triangularize,¡± which means to .nd a basis over which the matrix has zero entries below the main diagonal. Fortunately, such a basis always exist. 
14.2. REDUCTION TO UPPER TRIANGULAR FORM 
We say that a square matrix A is an upper triangular matrix if it has the following shape, 
.
. 
........ 

a11 a12 a13 ... a1 n.1 a1 n 
0 ... 

a22 a23 a2 n.1 a2 n 

00 ... 


a33 a3 n.1 a3 n 
. . .  . . .  . . .  ...  . . .  . . .  
0  0  0  . . .  an.1 n.1  an.1 n  
0  0  0  . . .  0  an n  

........ 

, 

i.e., aij = 0 whenever j<i,1 ¡Ü i, j ¡Ü n. 
Theorem 14.5. Given any .nite dimensional vector space over a .eld K, for any linear map f : E ¡ú E, there is a basis (u1,...,un) with respect to which f is represented by an upper triangular matrix (in Mn(K)) i. all the eigenvalues of f belong to K. Equivalently, for every n ¡Á n matrix A ¡Ê Mn(K), there is an invertible matrix P and an upper triangular matrix T (both in Mn(K)) such that 
A = PTP .1 
i. all the eigenvalues of A belong to K. 
Proof. If there is a basis (u1,...,un) with respect to which f is represented by an upper triangular matrix T in Mn(K), then since the eigenvalues of f are the diagonal entries of T , all the eigenvalues of f belong to K. 
For the converse, we proceed by induction on the dimension n of E. For n = 1 the result is obvious. If n> 1, since by assumption f has all its eigenvalue in K, pick some eigenvalue ¦Ë1 ¡Ê K of f, and let u1 be some corresponding (nonzero) eigenvector. We can .nd n . 1 vectors (v2,...,vn) such that (u1,v2,...,vn) is a basis of E, and let F be the subspace of dimension n . 1 spanned by (v2,...,vn). In the basis (u1,v2 ...,vn), the matrix of f is of 
the form 

.
. 
U = 

.... 

¦Ë1 a12 ... a1 n 
0 ... 
a22 a2 n 
.. .
.
....
. 0 an 2 ... ann 
.. . 
....

, 

since its .rst column contains the coordinates of ¦Ë1u1 over the basis (u1,v2, ...,vn). Ifwe let p: E ¡ú F be the projection de.ned such that p(u1)=0 and p(vi)= vi when 2 ¡Ü i ¡Ü n, the linear map g : F ¡ú F de.ned as the restriction of p . f to F is represented by the (n . 1) ¡Á (n . 1) matrix V =(aij)2¡Üi,j¡Ün over the basis (v2,...,vn). We need to prove that all the eigenvalues of g belong to K. However, since the .rst column of U has a single nonzero entry, we get 
¦ÖU (X) = det(XI . U)=(X . ¦Ë1) det(XI . V )=(X . ¦Ë1)¦ÖV (X), 
where ¦ÖU (X) is the characteristic polynomial of U and ¦ÖV (X) is the characteristic polynomial of V . It follows that ¦ÖV (X) divides ¦ÖU (X), and since all the roots of ¦ÖU (X) are in K, all the roots of ¦ÖV (X) are also in K. Consequently, we can apply the induction hypothesis, and there is a basis (u2,...,un) of F such that g is represented by an upper triangular matrix (bij)1¡Üi,j¡Ün.1. However, 
E = Ku1 ¨’ F, 
and thus (u1,...,un) is a basis for E. Since p is the projection from E = Ku1 ¨’ F onto F and g : F ¡ú F is the restriction of p . f to F , we have 
f(u1)= ¦Ë1u1 
and 
i 
f(ui+1)= a1 iu1 + bijuj+1 j=1 
for some a1 i ¡Ê K, when 1 ¡Ü i ¡Ü n . 1. But then the matrix of f with respect to (u1,...,un) is upper triangular. 
For the matrix version, we assume that A is the matrix of f with respect to some basis, Then we just proved that there is a change of basis matrix P such that A = PTP .1 where T is upper triangular. 
If A = PTP .1 where T is upper triangular, note that the diagonal entries of T are the eigenvalues ¦Ë1,...,¦Ën of A. Indeed, A and T have the same characteristic polynomial. Also, if A is a real matrix whose eigenvalues are all real, then P can be chosen to real, and if A is a rational matrix whose eigenvalues are all rational, then P can be chosen rational. Since any polynomial over C has all its roots in C, Theorem 14.5 implies that every complex n ¡Á n matrix can be triangularized. 
If E is a Hermitian space (see Chapter 13), the proof of Theorem 14.5 can be easily adapted to prove that there is an orthonormal basis (u1,...,un) with respect to which the matrix of f is upper triangular. This is usually known as Schur¡¯s lemma. 
Theorem 14.6. (Schur decomposition) Given any linear map f : E ¡ú E over a complex Hermitian space E, there is an orthonormal basis (u1,...,un) with respect to which f is represented by an upper triangular matrix. Equivalently, for every n ¡Á n matrix A ¡Ê Mn(C), there is a unitary matrix U and an upper triangular matrix T such that 
A = UTU . . 
If A is real and if all its eigenvalues are real, then there is an orthogonal matrix Q and a real upper triangular matrix T such that 
A = QT QT. 

14.2. REDUCTION TO UPPER TRIANGULAR FORM 
Proof. During the induction, we choose F to be the orthogonal complement of Cu1 and we pick orthonormal bases (use Propositions 13.13 and 13.12). If E is a real Euclidean space and if the eigenvalues of f are all real, the proof also goes through with real matrices (use Propositions 11.11 and 11.10). 
If ¦Ë is an eigenvalue of the matrix A and if u is an eigenvector associated with ¦Ë, from 
Au = ¦Ëu, 
we obtain A2 u = A(Au)= A(¦Ëu)= ¦ËAu = ¦Ë2 u, 
which shows that ¦Ë2 is an eigenvalue of A2 for the eigenvector u. An obvious induction shows that ¦Ëk is an eigenvalue of Ak for the eigenvector u, for all k ¡Ý 1. Now, if all eigenvalues ¦Ë1,...,¦Ën of A are in K, it follows that ¦Ëk 1,...,¦Ëkn are eigenvalues of Ak . However, it is not obvious that Ak does not have other eigenvalues. In fact, this can¡¯t happen, and this can be proven using Theorem 14.5. 
Proposition 14.7. Given any n ¡Á n matrix A ¡Ê Mn(K) with coe.cients in a .eld K, if all eigenvalues ¦Ë1,...,¦Ën of A are in K, then for every polynomial q(X) ¡Ê K[X], the eigenvalues of q(A) are exactly (q(¦Ë1),...,q(¦Ën)). 
Proof. By Theorem 14.5, there is an upper triangular matrix T and an invertible matrix P (both in Mn(K)) such that 
A = PTP .1 . 
Since A and T are similar, they have the same eigenvalues (with the same multiplicities), so the diagonal entries of T are the eigenvalues of A. Since 
Ak = PT kP .1 ,k ¡Ý 1, 
for any polynomial q(X)= c0Xm + ¡¤¡¤¡¤ + cm.1X + cm, we have 
q(A)= c0Am + ¡¤¡¤¡¤ + cm.1A + cmI 
= c0PT mP .1 + ¡¤¡¤¡¤ + cm.1PTP .1 + cmPIP .1 
= P (c0T m + ¡¤¡¤¡¤ + cm.1T + cmI)P .1 
= Pq(T )P .1 . 
Furthermore, it is easy to check that q(T ) is upper triangular and that its diagonal entries are q(¦Ë1),...,q(¦Ën), where ¦Ë1,...,¦Ën are the diagonal entries of T , namely the eigenvalues of A. It follows that q(¦Ë1),...,q(¦Ën) are the eigenvalues of q(A). 
Remark: There is another way to prove Proposition 14.7 that does not use Theorem 14.5, but instead uses the fact that given any .eld K, there is .eld extension K of K (K . K) such that every polynomial q(X)= c0Xm + ¡¤¡¤¡¤ + cm.1X + cm (of degree m ¡Ý 1) with coe.cients ci ¡Ê K factors as 
q(X)= c0(X . ¦Á1) ¡¤¡¤¡¤ (X . ¦Án),¦Ái ¡Ê K,i =1, . . . , n. 
The .eld K is called an algebraically closed .eld (and an algebraic closure of K). 
Assume that all eigenvalues ¦Ë1,...,¦Ën of A belong to K. Let q(X) be any polynomial (in K[X]) and let ¦Ì ¡Ê K be any eigenvalue of q(A) (this means that ¦Ì is a zero of the characteristic polynomial ¦Öq(A)(X) ¡Ê K[X] of q(A). Since K is algebraically closed, ¦Öq(A)(X) has all its roots in K). We claim that ¦Ì = q(¦Ëi) for some eigenvalue ¦Ëi of A. 
Proof. (After Lax [44], Chapter 6). Since K is algebraically closed, the polynomial ¦Ì . q(X) factors as 
¦Ì . q(X)= c0(X . ¦Á1) ¡¤¡¤¡¤ (X . ¦Án), 
for some ¦Ái ¡Ê K. Now ¦ÌI . q(A) is a matrix in Mn(K), and since ¦Ì is an eigenvalue of q(A), it must be singular. We have 
¦ÌI . q(A)= c0(A . ¦Á1I) ¡¤¡¤¡¤ (A . ¦ÁnI), 
and since the left-hand side is singular, so is the right-hand side, which implies that some factor A . ¦ÁiI is singular. This means that ¦Ái is an eigenvalue of A, say ¦Ái = ¦Ëi. As ¦Ái = ¦Ëi is a zero of ¦Ì . q(X), we get 
¦Ì = q(¦Ëi), 
which proves that ¦Ì is indeed of the form q(¦Ëi) for some eigenvalue ¦Ëi of A. 
Using Theorem 14.6, we can derive two very important results. 
Proposition 14.8. If A is a Hermitian matrix (i.e. A. = A), then its eigenvalues are real and A can be diagonalized with respect to an orthonormal basis of eigenvectors. In matrix terms, there is a unitary matrix U and a real diagonal matrix D such that A = UDU. . If A is a real symmetric matrix (i.e. AT = A), then its eigenvalues are real and A can be diagonalized with respect to an orthonormal basis of eigenvectors. In matrix terms, there is an orthogonal matrix Q and a real diagonal matrix D such that A = QDQT . 
Proof. By Theorem 14.6, we can write A = UTU. where T =(tij) is upper triangular and U is a unitary matrix. If A. = A, we get 
. U . 
UTU . = UT , 
and this implies that T = T . . Since T is an upper triangular matrix, T . is a lower triangular matrix, which implies that T is a diagonal matrix. Furthermore, since T = T ., we have tii = tii for i =1,...,n, which means that the tii are real, so T is indeed a real diagonal matrix, say D. 
If we apply this result to a (real) symmetric matrix A, we obtain the fact that all the eigenvalues of a symmetric matrix are real, and by applying Theorem 14.6 again, we conclude that A = QDQT, where Q is orthogonal and D is a real diagonal matrix. 

14.3. LOCATION OF EIGENVALUES 
More general versions of Proposition 14.8 are proven in Chapter 16. 
When a real matrix A has complex eigenvalues, there is a version of Theorem 14.6 involving only real matrices provided that we allow T to be block upper-triangular (the diagonal entries may be 2 ¡Á 2 matrices or real entries). 
Theorem 14.6 is not a very practical result but it is a useful theoretical result to cope with matrices that cannot be diagonalized. For example, it can be used to prove that every complex matrix is the limit of a sequence of diagonalizable matrices that have distinct eigenvalues! 


14.3 Location of Eigenvalues 
If A is an n ¡Á n complex (or real) matrix A, it would be useful to know, even roughly, where the eigenvalues of A are located in the complex plane C. The Gershgorin discs provide some precise information about this. 
De.nition 14.5. For any complex n ¡Á n matrix A, for i =1,...,n, let 
n 
R½Ð(A)= |aij|
ij=1 
Óõj
=i 
and let 
G(A)= 
n{z ¡Ê C ||z . aii|¡Ü R½Ð (A)}.
i
i=1
Each disc {z ¡Ê C ||z . aii|¡Ü R½Ð(A)} is called a Gershgorin disc and their union G(A) is 
i.. 
12 3 
called the Gershgorin domain.  An example of Gershgorin domain for A = .4  i  6  .  
7  8  1 + i  
is illustrated in Figure 14.1.  

Although easy to prove, the following theorem is very useful: 
Theorem 14.9. (Gershgorin¡¯s disc theorem) For any complex n ¡Án matrix A, all the eigen-values of A belong to the Gershgorin domain G(A). Furthermore the following properties hold: 
(1) If A is strictly row diagonally dominant, that is 
n |aii| > |aij|, for i =1, . . . , n, j=1,j=i
Óõ 
then A is invertible. 

Figure14.1: Let A bethe3 3matrixspeci.edattheendofDe.nition14.5. Forthis ¡Á ½Ð½Ð½Ðparticular A,we.ndthat R(A)5, R(A)10,and R(A)15.Theblue/purpledisk === 3	Óõ 
is|z.1|¡Ü5,thepinkdiski1s|z.i|¡Ü120, the peach disk is |z . 1 . i|¡Ü 15, and G(A) is the union of these three disks. 

(2) If A is strictly row diagonally dominant, and if 	aii > 0 for i =1,...,n, then every eigenvalue of A has a strictly positive real part. 
Proof. Let ¦Ë be any eigenvalue of A and let u be a corresponding eigenvector (recall that we must have u = 0). Let k be an index such that 
|uk| = max |ui|. 
1¡Üi¡Ün 
Since Au = ¦Ëu, we have 
n 
(¦Ë . akk)uk = akjuj, 
j=1 j=k 
which implies that 
n	n 
|¦Ë . akk||uk|¡Ü |akj||uj|¡Ü|uk||akj|. j=1 j=1 
Óõ
Óõ
j=k	j=k 
Since u = 0 and |uk| = max1¡Üi¡Ün |ui|, we must have |uk| = 0, and it follows that 
n 
|¦Ë . akk|¡Ü |akj| = Rk½Ð (A), 
j=1 
Óõ
j=k 
and thus 
¦Ë ¡Ê{z ¡Ê C ||z . akk|¡Ü R½Ð (A)}. G(A),
k

14.3. LOCATION OF EIGENVALUES 
as claimed. 
positive real part. 
version of Theorem 14.9 for the discs of radius 
n Cj½Ð (A)= |aij|, 
i=1 i=j 
(1)Strictrowdiagonaldominanceimpliesthat0doesnotbelongtoanyoftheGershgorin discs,soalleigenvaluesof A arenonzero,and A isinvertible. (2)IfAisstrictlyrowdiagonallydominantand 0for i =1,theneachofthe >a,...,nii Gershgorindiscsliesstrictlyintherighthalf-plane,soeveryeigenvalueof A hasastrictly Inparticular,Theorem14.9impliesthatifasymmetricmatrixisstrictlyrowdiagonally dominantandhasstrictlypositivediagonalentries,thenitispositivede.nite.Theorem14.9 issometimescalledthe Gershgorin¨CHadamardtheorem. TSince A and Ahavethesameeigenvalues(evenforcomplexmatrices)wealsohavea ÓõTwhosedomain G(A)isgivenby Theorem14.10. Foranycomplex matrix A,alltheeigenvaluesof A belongtothe ¡Ánn T¡ÉintersectionoftheGershgorindomains G(A) G(A)SeeFigure14.3. Furthermorethe . (1)IfAisstrictlycolumndiagonallydominant,thatis (2)IfAisstrictlycolumndiagonallydominant,andif 0 for i =1,thenevery >a,...,nii eigenvalueof A hasastrictlypositiverealpart. Therearere.nementsofGershgorin¡¯stheoremandeigenvaluelocationresultsinvolving otherdomainsbesidesdiscs;formoreonthissubject,seeHornandJohnson[36],Sections Remark: Neitherstrictrowdiagonaldominancenorstrictcolumndiagonaldominanceare necessaryforinvertibility. Also,ifwerelaxallstrictinequalitiestoinequalities,thenrow diagonaldominance(orcolumndiagonaldominance)isnotasu.cientconditionforinvert-

nT½ÐC{¡Ê||.|¡Ü }G(A)= C(A)zza.iiiÓõ 
i=1
.. 
12 3 Figure 14.2 shows G(AT) for A = .4 i 6 .. 7 8 1+ i Thus we get the following: 
following properties hold: 
n 
|aii| > |aij|, for j =1, . . . , n, i=1,i=j 
then A is invertible. 
6.1 and 6.2. 
ibility. 

Figure 14.2: Let A be the 3 ¡Á 3 matrix speci.ed at the end of De.nition 14.5. For this particular A, we .nd that C½Ð (A) = 11, C½Ð (A) = 10, and C3½Ð (A) = 9. The pale blue disk is |z.1|¡Ü1,thepinkdiskis|z1.i|¡Ü10,th2e ocher disk is |z . 1 . i|¡Ü 9, and G(AT) is the union of these three disks. 


14.4 Conditioning of Eigenvalue Problems 
The following n ¡Á n matrix 
. 

0  
A = ........ 1  0 1  0 ...  .. . 1  0  
1  0  

. ........ 

has the eigenvalue 0 with multiplicity n. However, if we perturb the top rightmost entry of A by ¦Å, it is easy to see that the characteristic polynomial of the matrix 
.

. ........ 

........ 
0 ¦Å 10 10 
A(¦Å)= 
..
..
.. 10 10 
It follows that if n = 40 and ¦Å = 10.40is Xn . ¦Å. 
k2¦Ði/40 
, A(10.40) has the eigenvalues 10.1e
with k =1,..., 40. Thus, we see that a very small change (¦Å = 10.40) to the matrix A causes 
k2¦Ði/40
a signi.cant change to the eigenvalues of A (from 0 to 10.1e). Indeed, the relative 
14.4. CONDITIONING OF EIGENVALUE PROBLEMS 

Figure 14.3: Let A be the 3 ¡Á3 matrix speci.ed at the end of De.nition 14.5. The dusty rose region is G(A) ¡ÉG(AT). 
error is 10.39 . Worse, due to machine precision, since very small numbers are treated as 0, the error on the computation of eigenvalues (for example, of the matrix A(10.40)) can be very large. 
This phenomenon is similar to the phenomenon discussed in Section 8.5 where we studied the e.ect of a small perturbation of the coe.cients of a linear system Ax = b on its solution. In Section 8.5, we saw that the behavior of a linear system under small perturbations is governed by the condition number cond(A) of the matrix A. In the case of the eigenvalue problem (.nding the eigenvalues of a matrix), we will see that the conditioning of the problem depends on the condition number of the change of basis matrix P used in reducing the matrix A to its diagonal form D = P .1AP , rather than on the condition number of A itself. The following proposition in which we assume that A is diagonalizable and that the matrix norm IIsatis.es a special condition (satis.ed by the operator norms IIp for p =1, 2, ¡Þ), is due to Bauer and Fike (1960). 
Proposition 14.11. Let A ¡ÊMn(C) be a diagonalizable matrix, P be an invertible matrix, and D be a diagonal matrix D = diag(¦Ë1,...,¦Ën) such that 
A = P DP .1 , 
and let IIbe a matrix norm such that 
Idiag(¦Á1,...,¦Án)I= max |¦Ái|, 
1¡Üi¡Ün 
for every diagonal matrix. Then for every perturbation matrix ¦¤A, if we write 
Bi = {z ¡ÊC ||z .¦Ëi|¡Ücond(P ) I¦¤AI}, 
for every eigenvalue ¦Ë of A +¦¤A, we have 

n
¦Ë ¡Ê Bk. 
k=1 

Proof. Let ¦Ë be any eigenvalue of the matrix A +¦¤A. If ¦Ë = ¦Ëj for some j, then the result is trivial. Thus assume that ¦Ë = ¦Ëj for j =1,...,n. In this case the matrix D . ¦ËI is invertible (since its eigenvalues are ¦Ë . ¦Ëj for j =1,...,n), and we have 
P .1(A +¦¤A . ¦ËI)P = D . ¦ËI + P .1(¦¤A)P =(D . ¦ËI)(I +(D . ¦ËI).1P .1(¦¤A)P ). Since ¦Ë is an eigenvalue of A +¦¤A, the matrix A +¦¤A . ¦ËI is singular, so the matrix I +(D . ¦ËI).1P .1(¦¤A)P must also be singular. By Proposition 8.11(2), we have 
  

(D . ¦ËI).1P .1(¦¤A)P
  

,

1 ¡Ü

and since II is a matrix norm,

  

  

¡Ü

  

(D . ¦ËI).1
  
  
P .1
  

(D . ¦ËI).1P .1(¦¤A)P
I¦¤AIIP I , 

so we have 

1 ¡Ü

  

(D . ¦ËI).1
  
  
P .1
  

I¦¤AIIP I . 

Now (D . ¦ËI).1 is a diagonal matrix with entries 1/(¦Ëi . ¦Ë), so by our assumption on the 
norm,

  

(D . ¦ËI).1
  

1 

= . 
mini(|¦Ëi . ¦Ë|) 
As a consequence, since there is some index k for which mini(|¦Ëi . ¦Ë|)= |¦Ëk . ¦Ë|, we have
  

  

1 

= ,
|¦Ëk . ¦Ë|
(D . ¦ËI).1
and we obtain 

  

P .1
  

|¦Ë . ¦Ëk|¡Ü
I¦¤AIIP I = cond(P ) I¦¤AI , 

which proves our result. Proposition 14.11 implies that for any diagonalizable matrix A, if we de.ne ¦£(A) by ¦£(A) = inf{cond(P ) | P .1AP = D}, then for every eigenvalue ¦Ë of A +¦¤A, we have 

n
¦Ë ¡Ê{z ¡Ê Cn ||z . ¦Ëk|¡Ü ¦£(A) I¦¤AI}. 
k=1
14.5. EIGENVALUES OF THE MATRIX EXPONENTIAL 
De.nition 14.6. The number ¦£(A) = inf{cond(P ) | P .1AP = D} is called the conditioning of A relative to the eigenvalue problem. 
If A is a normal matrix, since by Theorem 16.22, A can be diagonalized with respect to a unitary matrix U, and since for the spectral norm IUI2 = 1, we see that ¦£(A) = 1. Therefore, normal matrices are very well conditionned w.r.t. the eigenvalue problem. In fact, for every eigenvalue ¦Ë of A +¦¤A (with A normal), we have 

n¦Ë ¡Ê{z ¡Ê Cn ||z . ¦Ëk|¡ÜI¦¤AI2}. k=1
If A and A+¦¤A are both symmetric (or Hermitian), there are sharper results; see Proposition 
16.28. Note that the matrix A(¦Å) from the beginning of the section is not normal. 


14.5 Eigenvalues of the Matrix Exponential 
The Schur decomposition yields a characterization of the eigenvalues of the matrix exponen-tial eA in terms of the eigenvalues of the matrix A. First we have the following proposition. 
Proposition 14.12. Let A and U be (real or complex) matrices and assume that U is invertible. Then 
UAU.1 
e = UeAU.1 . 
Proof. A trivial induction shows that 
UApU.1 =(UAU.1)p, 
and thus 
(UAU.1)p UApU.1 
UAU.1 
e == 
p! p! 
p¡Ý0 p¡Ý0 
U.1 
= U Ap! p = UeAU.1 , p¡Ý0 
as claimed. 
Proposition 14.13. Given any complex n ¡Á n matrix A, if ¦Ë1,...,¦Ën are the eigenvalues of A, then e¦Ë1 ,...,e¦Ën are the eigenvalues of eA . Furthermore, if u is an eigenvector of A for ¦Ëi, then u is an eigenvector of eA for e¦Ëi . 
Proof. By Theorem 14.5, there is an invertible matrix P and an upper triangular matrix T such that 
A = PTP .1 . 
By Proposition 14.12, 
PTP .1 
e = PeT P .1 . 
 
Note that eT =T p is upper triangular since T p is upper triangular for all p ¡Ý 0. If 
p¡Ý0 p! ¦Ë1,¦Ë2,...,¦Ën are the diagonal entries of T , the properties of matrix multiplication, when combined with an induction on p, imply that the diagonal entries of T p are ¦Ë1p,¦Ëp 2,...,¦Ëp . 
 ¦Ëp nT i ¦Ëi
This in turn implies that the diagonal entries of eare= efor 1 ¡Ü i ¡Ü n. Since 
p¡Ý0 p! 
A and T are similar matrices, we know that they have the same eigenvalues, namely the 
A PTP .1 T
diagonal entries ¦Ë1,...,¦Ën of T . Since e= e= PeT P .1, and eis upper triangular, we use the same argument to conclude that both eA and eT have the same eigenvalues, which are the diagonal entries of eT , where the diagonal entries of eT are of the form e¦Ë1 ,...,e¦Ën . Now, if u is an eigenvector of A for the eigenvalue ¦Ë, a simple induction shows that u is an eigenvector of An for the eigenvalue ¦Ën, from which is follows that 
  
A AA2 A3 A2 A3 
eu =I ++ + + ...u = u + Au + u + u + ... 
1!2! 3! 2! 3! 
  
¦Ë2 ¦Ë3 ¦Ë2 ¦Ë3 = u + ¦Ëu + u + u + ¡¤¡¤¡¤ =1+ ¦Ë +++ ...u = e ¦Ë u,
2! 3! 2!3! 
which shows that u is an eigenvector of eA for e¦Ë . 
As a consequence, we obtain the following result. 
Proposition 14.14. For every complex (or real) square matrix A, we have 
tr(A)
det(e A)= e, 
where tr(A) is the trace of A, i.e., the sum a11 + ¡¤¡¤¡¤ + ann of its diagonal entries. 
Proof. The trace of a matrix A is equal to the sum of the eigenvalues of A. The determinant of a matrix is equal to the product of its eigenvalues, and if ¦Ë1,...,¦Ën are the eigenvalues of A, then by Proposition 14.13, e¦Ë1 ,...,e¦Ën are the eigenvalues of eA, and thus 
  
A¦Ë1 ¦Ën ¦Ë1+¡¤¡¤¡¤+¦Ën tr(A)
dete = e ¡¤¡¤¡¤ e = e = e, 
as desired. 
If B is a skew symmetric matrix, since tr(B) = 0, we deduce that det(eB)= e0 = 1. This allows us to obtain the following result. Recall that the (real) vector space of skew symmetric matrices is denoted by so(n). 
Proposition 14.15. For every skew symmetric matrix B ¡Ê so(n), we have eB ¡Ê SO(n), that is, eB is a rotation. 
Proof. By Proposition 8.23, eB is an orthogonal matrix. Since tr(B) = 0, we deduce that det(eB)= e0 = 1. Therefore, eB ¡Ê SO(n). 
14.6. SUMMARY 
Proposition 14.15 shows that the map B ¡ú eB is a map exp: so(n) ¡ú SO(n). Itisnot injective, but it can be shown (using one of the spectral theorems) that it is surjective. 
If B is a (real) symmetric matrix, then 
BT B
(e B)T = e = e, 
so eB is also symmetric. Since the eigenvalues ¦Ë1,...,¦Ën of B are real, by Proposition 
B¦Ën 	¦Ëi
14.13, since the eigenvalues of eare e¦Ë1 ,...,eand the ¦Ëi are real, we have e> 0 for i =1,...,n, which implies that eB is symmetric positive de.nite. In fact, it can be shown that for every symmetric positive de.nite matrix A, there is a unique symmetric matrix B such that A = eB; see Gallier [25]. 


14.6 Summary 
The main concepts and results of this chapter are listed below: 
. 	
Diagonal matrix. 

. 	
Eigenvalues, eigenvectors; the eigenspace associated with an eigenvalue. 

. 	
Characteristic polynomial. 

. 	
Trace. 

. 	
Algebraic and geometric multiplicity. 

. 	
Eigenspaces associated with distinct eigenvalues form a direct sum (Proposition 14.3). 

. 	
Reduction of a matrix to an upper-triangular matrix. 

. 	
Schur decomposition. 

. 	
The Gershgorin¡¯s discs can be used to locate the eigenvalues of a complex matrix; see Theorems 14.9 and 14.10. 

. 	
The conditioning of eigenvalue problems. 

tr(A)

. 	
Eigenvalues of the matrix exponential. The formula det(eA)= e. 



14.7 Problems 
Problem 14.1. Let A be the following 2 ¡Á 2 matrix 1 .1 
A = . 
1 .1 
(1) 
Prove that A has the eigenvalue 0 with multiplicity 2 and that A2 = 0. 

(2) 
Let A be any real 2 ¡Á 2 matrix 


ab 
A = . 
cd 
Prove that if bc > 0, then A has two distinct real eigenvalues. Prove that if a,b,c,d > 0, then there is a positive eigenvector u associated with the largest of the two eigenvalues of A, which means that if u =(u1,u2), then u1 > 0 and u2 > 0. 
(3) Suppose now that A is any complex 2 ¡Á 2 matrix as in (2). Prove that if A has the eigenvalue 0 with multiplicity 2, then A2 = 0. Prove that if A is real symmetric, then A = 0. 
Problem 14.2. Let A be any complex n ¡Á n matrix. Prove that if A has the eigenvalue 0 with multiplicity n, then An = 0. Give an example of a matrix A such that An = 0 but A = 0. 
Problem 14.3. Let A be a complex 2 ¡Á 2 matrix, and let ¦Ë1 and ¦Ë2 be the eigenvalues of 
A. Prove that if ¦Ë1 = ¦Ë2, then 
¦Ë2 ¦Ë1 ¦Ë1 ¦Ë2 
A ¦Ë1e. ¦Ë2ee. e
e = I + A. 
¦Ë1 . ¦Ë2 ¦Ë1 . ¦Ë2 
Problem 14.4. Let A be the real symmetric 2 ¡Á 2 matrix 
ab 
A = . 
bc 
(1) Prove that the eigenvalues of A are real and given by 
a + c +4b2 +(a . c)2 a + c . 4b2 +(a . c)2 
¦Ë1 = ,¦Ë2 = . 
22 
(2) Prove that A has a double eigenvalue (¦Ë1 = ¦Ë2 = a) if and only if b = 0 and a = c; that is, A is a diagonal matrix. 
(3) 
Prove that the eigenvalues of A are nonnegative i. b2 ¡Ü ac and a + c ¡Ý 0. 

(4) 
Prove that the eigenvalues of A are positive i. b2 < ac, a> 0 and c> 0. 


Problem 14.5. Find the eigenvalues of the matrices 
3011 41 
A = ,B = ,C = A + B = . 
1103 14 
Check that the eigenvalues of A + B are not equal to the sums of eigenvalues of A plus eigenvalues of B. 
14.7. PROBLEMS 
Problem 14.6. Let A be a real symmetric n ¡Á n matrix and B be a real symmetric positive de.nite n ¡Á n matrix. We would like to solve the generalized eigenvalue problem: .nd ¦Ë ¡Ê R and u = 0 such that 
Au = ¦ËBu. (.) 
(1) Use the Cholseky decomposition B = CCT to show that ¦Ë and u are solutions of the generalized eigenvalue problem (.) i. ¦Ë and v are solutions the (ordinary) eigenvalue problem 
C.1A(CT).1 v = ¦Ëv, with v = CTu. 
Check that C.1A(CT).1 is symmetric. 
(2) Prove that if Au1 = ¦Ë1Bu1, Au2 = ¦Ë2Bu2, with u1 = 0, u2 = 0 and ¦Ë1 = ¦Ë2, then 
T
u1 Bu2 = 0. 
(3) Prove that B.1A and C.1A(CT).1 have the same eigenvalues. 
Problem 14.7. The sequence of Fibonacci numbers,0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55,..., is given by the recurrence 
Fn+2 = Fn+1 + Fn, 
with F0 = 0 and F1 = 1. In matrix form, we can write 
Fn+1 11 Fn F1 1 
= ,n ¡Ý 1, = . 
Fn 10 Fn.1 F0 0 
(1) Show that 
Fn+1 11 n 1 

= . 
Fn 10 0 
(2) Prove that the eigenvalues of the matrix 
11 
A = 
10 
are ¡Ì 1 ¡À 5 
¦Ë = . 
2 The number 
¡Ì 
1+ 5 
. = 
2 is called the golden ratio. Show that the eigenvalues of A are . and ...1 . 
(3) Prove that A is diagonalized as 
11 1 . ...1 . 01 ..1 
A == ¡Ì . 
10 110 ...1 .1 .
5 
Prove that 
1 ...1 .n
Fn+1 . 
= ¡Ì ,
Fn 51 1 .(...1)n 
and thus 
 

¡Ì n ¡Ì n
 

1 11+ 5 1 . 5 
Fn = ¡Ì (.n . (...1)n)= ¡Ì. ,n ¡Ý 0. 
5522 
Problem 14.8. Let A be an n ¡Á n matrix. For any subset I of {1,...,n}, let AI,I be the matrix obtained from A by .rst selecting the columns whose indices belong to I, and then the rows whose indices also belong to I. Prove that 
¦Ók(A) = det(AI,I ). 
I.{1,...,n}|I|=k 
Problem 14.9. (1) Consider the matrix 
.
. 
00 .a3 
A =

.

10 

.a2
.

. 

01 .a1 Prove that the characteristic polynomial ¦ÖA(z) = det(zI . A) of A is given by ¦ÖA(z)= z 3 + a1z 2 + a2z + a3. 
(2) Consider the matrix 

.
. 
... 

000 .a4 100 .a3 
010 .a2 001 .a1 
...

A = 

. 

Prove that the characteristic polynomial ¦ÖA(z) = det(zI . A) of A is given by 
432
¦ÖA(z)= z + a1z + a2z + a3z + a4. 
(3) Consider the n ¡Á n matrix (called a companion matrix) 
.
. 
00 0 ¡¤¡¤¡¤ 0 .an 10 0 ¡¤¡¤¡¤ 0 .an.1 01 0 ¡¤¡¤¡¤ 0 .an.2 
A = 

........ 

........ 

.

. ..
...
. 

.

.

.. .

. 

. 

.

. 

.. 

.
.
000 .0 .a2 00 0 ¡¤¡¤¡¤ 1 .a1 
Prove that the characteristic polynomial ¦ÖA(z) = det(zI . A) of A is given by 
nn.1 n.2
¦ÖA(z)= z + a1z + a2z + ¡¤¡¤¡¤ + an.1z + an. 

14.7. PROBLEMS 
Hint. Use induction. 
Explain why .nding the roots of a polynomial (with real or complex coe.cients) and .nding the eigenvalues of a (real or complex) matrix are equivalent problems, in the sense that if we have a method for solving one of these problems, then we have a method to solve the other. 
Problem 14.10. Let A be a complex n ¡Á n matrix. Prove that if A is invertible and if the eigenvalues of A are (¦Ë1,...,¦Ën), then the eigenvalues of A.1 are (¦Ë.11,...,¦Ë.n 1). Prove that if u is an eigenvector of A for ¦Ëi, then u is an eigenvector of A.1 for ¦Ë.i 1 . 
Problem 14.11. Prove that every complex matrix is the limit of a sequence of diagonalizable matrices that have distinct eigenvalues 
Problem 14.12. Consider the following tridiagonal n ¡Á n matrices 
.
.
.
. 
2 .10 010 

A = 

...... 

...... 

, 

S 

= 

...... 

...... 

. 

.12 .1 

10 1 

...

.

. 

.

. 

.

...

. 

.

. 

.

. 

.

. 

.12 .1 

1 01 

0 .12 010 
Observe that A =2I . S and show that the eigenvalues of A are ¦Ëk =2 . ¦Ìk, where the ¦Ìk are the eigenvalues of S. 
(2) Using Problem 9.6, prove that the eigenvalues of the matrix A are given by 
¦Ëk = 4 sin2 k¦Ð ,k =1, . . . , n. 
2(n + 1) Show that A is symmetric positive de.nite. 
(3) 
Find the condition number of A with respect to the 2-norm. 

(k)(k)

(4) 
Show that an eigenvector (y1 ,...,yn ) associated wih the eigenvalue ¦Ëk is given by 


(k) kj¦Ð 
y= sin ,j =1, . . . , n. 
j 
n +1 Problem 14.13. Consider the following real tridiagonal symmetric n ¡Á n matrix 
.
. 
c 10 1 c 1 ...
A = 

...... 

...... 

.

.

. 

.

. 

.

. 

1 c 1 
01 c 

(1) Using Problem 9.6, prove that the eigenvalues of the matrix A are given by k¦Ð 
¦Ëk = c + 2 cos ,k =1, . . . , n. 
n +1 
(2) Find a condition on c so that A is positive de.nite. It is satis.ed by c = 4? 
Problem 14.14. Let A be an m ¡Á n matrix and B be an n ¡Á m matrix (over C). 
(1) Prove that det(Im . AB) = det(In . BA). 
Hint. Consider the matrices  
X =  Im B  A In  and  Y  =  Im .B  0 In  .  
(2) Prove that  

¦Ën det(¦ËIm . AB)= ¦Ëm det(¦ËIn . BA). Hint. Consider the matrices ¦ËIm AIm 0 
X = and Y = . 
BIn .B ¦ËIn Deduce that AB and BA have the same nonzero eigenvalues with the same multiplicity. Problem 14.15. The purpose of this problem is to prove that the characteristic polynomial 
of the matrix 

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

is  PA(¦Ë) = ¦Ën.2  ¦Ë2 . n 2¦Ë . 1 12 n 2(n 2 . 1)  .  
(1) Prove that the characteristic polynomial PA(¦Ë) is given by PA(¦Ë) = ¦Ën.2P (¦Ë),  
with  ¦Ë . 1  .2  .3  .4  ¡¤ ¡¤ ¡¤  .n + 3  .n + 2  .n + 1  .n  
.¦Ë . 1  ¦Ë . 1  .1  .1  ¡¤ ¡¤ ¡¤  .1  .1  .1  .1  
1  .2  1  0  ¡¤ ¡¤ ¡¤  0  0  0  0  
0  1  .2  1  ¡¤ ¡¤ ¡¤  0  0  0  0  

.. ....
...
. ..... . ..
...
.. ....
P (¦Ë)= . 
0  0  0  0  ...  1  0  0  0  
0  0  0  0  ...  .2  1  0  0  
0  0  0  0  ¡¤ ¡¤ ¡¤  1  .2  1  0  
0  0  0  0  ¡¤ ¡¤ ¡¤  0  1  .2  1  


14.7. PROBLEMS 
(2) Prove that the sum of the roots ¦Ë1,¦Ë2 of the (degree two) polynomial P (¦Ë) is 
¦Ë1 + ¦Ë2 = n 2 . 
The problem is thus to compute the product ¦Ë1¦Ë2 of these roots. Prove that 
¦Ë1¦Ë2 = P (0). 
(3) The problem is now to evaluate dn = P (0), where 
.1 .2 .3 .4 ¡¤¡¤¡¤ .n +3 .n +2 .n +1 .n .1 .1 .1 .1 ¡¤¡¤¡¤ .1 .1 .1 .1 1 .21 0 ¡¤¡¤¡¤ 0 0 00 01 .21 ¡¤¡¤¡¤ 0 0 00 
.. ....
...
...... . ..
dn =...... . . . .
.
0000.1 0 00 .
.
0000 . .21 00 0000 ¡¤¡¤¡¤ 1 .2 10 0000 ¡¤¡¤¡¤ 01 .21 
I suggest the following strategy: cancel out the .rst entry in row 1 and row 2 by adding 
a suitable multiple of row 3 to row 1 and row 2, and then subtract row 2 from row 1. Do this twice. You will notice that the .rst two entries on row 1 and the .rst two entries on row 2 
change, but the rest of the matrix looks the same, except that the dimension is reduced. This suggests setting up a recurrence involving the entries uk,vk,xk,yk in the determinant 
uk xk .3 .4 ¡¤¡¤¡¤ .n + k . 3 .n + k . 2 .n + k . 1 .n + k 
vk yk .1 .1 ¡¤¡¤¡¤ .1 .1 .1 .1 1 .21 0 ¡¤¡¤¡¤ 0 0 00 01 .21 ¡¤¡¤¡¤ 0 0 00 
.. ....
...
...... . ..
Dk =...... . . . , .
.
0000.1 0 00 .
.
0000 . .21 00 
00 0 0 ¡¤¡¤¡¤ 1 .21 
0 00 0 0 ¡¤¡¤¡¤ 01 .2 
1 
starting with k = 0, with u0 = .1,v0 = .1,x0 = .2,y0 = .1, and ending with k = n . 2, so that un.3 xn.3 .3 un.2 xn.2
dn = Dn.2 = vn.3 yn.3 .1= . vn.2 yn.2
1 .21 
Prove that we have the recurrence relations 

...
..
.
.
. 
uk+1 2 .21 .1 uk 0 
... 

vk+1 
xk+1 
... 

= 

... 

0 201 

.1 

100

... 
... 

vk 
xk 
...

+ 

... 

0 

.2 

...

. 

yk+1 0 .10 0 yk .1 
These appear to be nasty a.ne recurrence relations, so we will use the trick to convert this a.ne map to a linear map. 
(4) Consider the linear map given by 
.
..
... 
uk+1 2 .21 .10 uk 
.....
.1 1 00001 1 
and show that its action on uk,vk,xk,yk is the same as the a.ne action of Part (3). Use Matlab to .nd the eigenvalues of the matrix 
..... 

..... 

..... 

..... 

..... 

0 201 0 .2 
vk+1 
vk 
.1 

100

xk+1 
= 

xk 
, 

.1

0 

00

yk+1 
yk 
.
. 
..... 

2  .2  1  .1  0  
0  2  0  1  0  
.1  1  0  0  .2  

0 .10 0 .1 
0 000 1 

..... 

T = 

. 

You will be stunned! Let N be the matrix given by  N = T . I.  
Prove that  N4 = 0.  
Use this to prove that  

T k = I + kN +1 k(k . 1)N2 +1 k(k . 1)(k . 2)N3 ,
26 

14.7. PROBLEMS 
for all k ¡Ý 0. 
(5) Prove that 
.
...
..
k .
. 

.12 .21 .1 
..... 
.1

0

uk 
..... 

..... 

= T k 
.1 .2 .1 
..... 
..... 

..... 
..... 

..... 

.1 

.2 

.1 

0 201 0

vk 
xk 
yk 
.1100 .2

= 

, 

.10 0 .1

0 

1 1 00001 1 
for k ¡Ý 0. Prove that 
.

. 

k + 1  .k(k + 1)  k  .k2  1 6 (k . 1)k(2k . 7)  
0  k + 1  0  k  .1 2 (k . 1)k  
.k  k2  1 . k  (k . 1)k  .1 3 k((k . 6)k + 11)  
0  .k  0  1 . k  1 2 (k . 3)k  
0  0  0  0  1  

........ 

........ 

T k 
= 

, 

and thus that

.
.
.

.

1 
6 (2k3 +3k2 . 5k . 6)
uk 
..... 

..... 

= 

..... 

..... 

.12 (k2 +3k + 2) 
1 
3 (.k3 + k . 6) 
vk 
xk 
, 

yk 21 (k2 + k . 2) and that 
uk xk 7232 1 
= .1 . k . k2 . k3 . k4 . 
vk yk 312 3 12 As a consequence, prove that amazingly 
dn = Dn.2 = . 1 n 2(n 2 . 1). 
12 
(6) Prove that the characteristic polynomial of A is indeed PA(¦Ë)= ¦Ën.2 ¦Ë2 . n 2¦Ë . 1 n 2(n 2 . 1) . 
12 Use the above to show that the two nonzero eigenvalues of A are 
¡Ì ¡Ì
n 3 
¦Ë = n ¡À 4n2 . 1 . 
23 

The negative eigenvalue ¦Ë1 can also be expressed as 
¡Ì 

2 (3 . 23) 1 
¦Ë1 = n 1 . . 
2
64n
Use this expression to explain the following phenomenon: if we add any number greater than or equal to (2/25)n2 to every diagonal entry of A we get an invertible matrix. What about 0.077351n2? Try it! 
Problem 14.16. Let A be a symmetric tridiagonal n ¡Á n-matrix 
.
. 
A = 

........ 

b1 c1 c1 b2 c2 c2 b3 c3 
........ 

,

...  ...  ...  
cn.2  bn.1  cn.1  
cn.1  bn  

where it is assumed that ci = 0 for all i,1 ¡Ü i ¡Ü n . 1, and let Ak be the k ¡Á k-submatrix consisting of the .rst k rows and columns of A,1 ¡Ü k ¡Ü n. We de.ne the polynomials Pk(x) as follows: (0 ¡Ü k ¡Ü n). 
P0(x)=1, 
P1(x)= b1 . x, 
Pk(x)=(bk . x)Pk.1(x) . c 2 k.1Pk.2(x), 

where 2 ¡Ü k ¡Ü n. 
(1) 
Prove the following properties: 

(i) 
Pk(x) is the characteristic polynomial of Ak, where 1 ¡Ü k ¡Ü n. 

(ii) 
limx¡ú.¡Þ Pk(x)=+¡Þ, where 1 ¡Ü k ¡Ü n. 


(iii) If Pk(x) = 0, then Pk.1(x)Pk+1(x) < 0, where 1 ¡Ü k ¡Ü n . 1. 
(iv) 
Pk(x) has k distinct real roots that separate the k + 1 roots of Pk+1(x), where 1 ¡Ü k ¡Ü n . 1. 

(2) 
Given any real number ¦Ì> 0, for every k,1 ¡Ü k ¡Ü n, de.ne the function sgk(¦Ì) as follows: 


 

sign of Pk(¦Ì) if Pk(¦Ì) = 0, 
sgk(¦Ì)=
sign of Pk.1(¦Ì) if Pk(¦Ì) = 0. We encode the sign of a positive number as +, and the sign of a negative number as .. Then let E(k, ¦Ì) be the ordered list E(k, ¦Ì)=£¨+, sg1(¦Ì), sg2(¦Ì), ..., sgk(¦Ì)) , and let N(k, ¦Ì) be the number changes of sign between consecutive signs in E(k, ¦Ì). 

14.7. PROBLEMS 
Prove that sgk(¦Ì) is well de.ned and that N(k, ¦Ì) is the number of roots ¦Ë of Pk(x) such that ¦Ë<¦Ì. 
Remark: The above can be used to compute the eigenvalues of a (tridiagonal) symmetric matrix (the method of Givens-Householder). 



