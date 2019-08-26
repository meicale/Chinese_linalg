Chapter 6 
Determinants 
In this chapter all vector spaces are de.ned over an arbitrary .eld K. For the sake of concreteness, the reader may safely assume that K = R. 
6.1 Permutations, Signature of a Permutation 
This chapter contains a review of determinants and their use in linear algebra. We begin with permutations and the signature of a permutation. Next we de.ne multilinear maps and alternating multilinear maps. Determinants are introduced as alternating multilinear maps taking the value 1 on the unit matrix (following Emil Artin). It is then shown how to compute a determinant using the Laplace expansion formula, and the connection with the usual de.nition is made. It is shown how determinants can be used to invert matrices and to solve (at least in theory!) systems of linear equations (the Cramer formulae). The determinant of a linear map is de.ned. We conclude by de.ning the characteristic polynomial of a matrix (and of a linear map) and by proving the celebrated Cayley�CHamilton theorem which states that every matrix is a ��zero�� of its characteristic polynomial (we give two proofs; one computational, the other one more conceptual). 
Determinants can be de.ned in several ways. For example, determinants can be de.ned in a fancy way in terms of the exterior algebra (or alternating algebra) of a vector space. We will follow a more algorithmic approach due to Emil Artin. No matter which approach is followed, we need a few preliminaries about permutations on a .nite set. We need to show that every permutation on n elements is a product of transpositions and that the parity of the number of transpositions involved is an invariant of the permutation. Let [n]= {1, 2 ...,n}, where n �� N, and n> 0. 
De.nition 6.1. A permutation on n elements is a bijection �� :[n] �� [n]. When n = 1, the only function from [1] to [1] is the constant map: 1 �� 1. Thus, we will assume that n �� 2. A transposition is a permutation �� :[n] �� [n] such that, for some i<j (with 1 �� i<j �� n), ��(i)= j, ��(j)= i, and ��(k)= k, for all k �� [n] .{i, j}. In other words, a transposition exchanges two distinct elements i, j �� [n]. 
153 

If �� is a transposition, clearly, �� . �� = id. We will also use the terminology product of permutations (or transpositions) as a synonym for composition of permutations. 
A permutation �� on n elements, say ��(i)= ki for i =1,...,n, can be represented in functional notation by the 2 �� n array 
1 ������ i ������ n k1 ������ ki ������ kn 
known as Cauchy two-line notation. For example, we have the permutation �� denoted by 
.  .  
1  2  3  4  5  6  
2  4  3  6  5  1  .  

A more concise notation often used in computer science and in combinatorics is to rep-resent a permutation by its image, namely by the sequence 
��(1) ��(2) ������ ��(n) 
written as a row vector without commas separating the entries. The above is known as the one-line notation. For example, in the one-line notation, our previous permutation �� is represented by 
243651. 
The reason for not enclosing the above sequence within parentheses is avoid confusion with the notation for cycles, for which is it customary to include parentheses. 
Clearly, the composition of two permutations is a permutation and every permutation has an inverse which is also a permutation. Therefore, the set of permutations on [n] isa group often denoted Sn and called the symmetric group on n elements. 
It is easy to show by induction that the group Sn has n! elements. The following propo-sition shows the importance of transpositions. 
Proposition 6.1. For every n �� 2, every permutation �� :[n] �� [n] can be written as a nonempty composition of transpositions. 
Proof. We proceed by induction on n. If n = 2, there are exactly two permutations on [2], the transposition �� exchanging 1 and 2, and the identity. However, id2 = �� 2 . Now let n �� 3. If ��(n)= n, since by the induction hypothesis, the restriction of �� to [n . 1] can be written as a product of transpositions, �� itself can be written as a product of transpositions. If ��(n)= k n, letting �� be the transposition such that ��(n)= k and ��(k)= n, it is clear 
= that �� . �� leaves n invariant, and by the induction hypothesis, we have �� . �� = ��m . ... . ��1 for some transpositions, and thus 
�� = �� . ��m . ... . ��1, 
a product of transpositions (since �� . �� = idn). 
6.1. PERMUTATIONS, SIGNATURE OF A PERMUTATION 
Remark: When �� = idn is the identity permutation, we can agree that the composition of 0 transpositions is the identity. Proposition 6.1 shows that the transpositions generate the group of permutations Sn. 
A transposition �� that exchanges two consecutive elements k and k +1 of [n] (1 �� k �� n.1) may be called a basic transposition. We leave it as a simple exercise to prove that every transposition can be written as a product of basic transpositions. In fact, the transposition that exchanges k and k + p (1 �� p �� n . k) can be realized using 2p . 1 basic transpositions. Therefore, the group of permutations Sn is also generated by the basic transpositions. 
Given a permutation written as a product of transpositions, we now show that the parity of the number of transpositions is an invariant. For this, we introduce the following function. 
De.nition 6.2. For every n �� 2, let ��: Zn �� Z be the function given by 
 
��(x1,...,xn)=(xi . xj). 1��i<j��n 
More generally, for any permutation �� �� Sn, de.ne ��(x��(1),...,x��(n)) by 
 
��(x��(1),...,x��(n))=(x��(i) . x��(j)). 
1��i<j��n 
The expression ��(x1,...,xn) is often called the discriminant of (x1,...,xn). 
  
��(x1,...,xn) = 0. The discriminant consists ofn 2factors. When n = 3, 
��(x1,x2,x3)=(x1 . x2)(x1 . x3)(x2 . x3). 
If �� is the permutation  
1  2  3  
2  3  2  ,  
then  

��(x��(1),x��(2),x��(3))=(x��(1) . x��(2))(x��(1) . x��(3))(x��(2) . x��(3)) =(x2 . x3)(x2 . x1)(x3 . x1). 
Observe that ��(x��(1),x��(2),x��(3))=(.1)2��(x1,x2,x3), 
since two transpositions applied to the identity permutation 123 (written in one-line notation) give rise to 2 3 1. This result regarding the parity of ��(x��(1),...,x��(n)) is generalized by the following proposition. 
Proposition 6.2. For every basic transposition �� of [n] (n �� 2), we have 
��(x��(1),...,x��(n))= .��(x1,...,xn). 
The above also holds for every transposition, and more generally, for every composition of transpositions �� = ��p .�� ����. ��1, we have 
��(x��(1),...,x��(n))=(.1)p��(x1,...,xn). 
Consequently, for every permutation �� of [n], the parity of the number p of transpositions involved in any decomposition of �� as �� = ��p .�� ����. ��1 is an invariant (only depends on ��). 
Proof. Suppose �� exchanges xk and xk+1. The terms xi . xj that are a.ected correspond to i = k, or i = k + 1, or j = k, or j = k + 1. The contribution of these terms in ��(x1,...,xn) is 
(xk . xk+1)[(xk . xk+2) ������ (xk . xn)][(xk+1 . xk+2) ������ (xk+1 . xn)] [(x1 . xk) ������ (xk.1 . xk)][(x1 . xk+1) ������ (xk.1 . xk+1)]. 
When we exchange xk and xk+1, the .rst factor is multiplied by .1, the second and the third factor are exchanged, and the fourth and the .fth factor are exchanged, so the whole product ��(x1,...,xn) is is indeed multiplied by .1, that is, 
��(x��(1),...,x��(n))= .��(x1,...,xn). 
For the second statement, .rst we observe that since every transposition �� can be written as the composition of an odd number of basic transpositions (see the the remark following Proposition 6.1), we also have 
��(x��(1),...,x��(n))= .��(x1,...,xn). 
Next we proceed by induction on the number p of transpositions involved in the decompo-sition of a permutation ��. 
The base case p = 1 has just been proven. If p �� 2, if we write �� = ��p.1 .�� ����. ��1, then �� = ��p . �� and 
��(x��(1),...,x��(n)) = ��(x��p(��(1)),...,x��p(��(n))) 
= .��(x��(1),...,x��(n)) 
= .(.1)p.1��(x1,...,xn) 
=(.1)p��(x1,...,xn), 
where we used the induction hypothesis from the second to the third line, establishing the induction hypothesis. Since ��(x��(1),...,x��(n)) only depends on ��, the equation 
��(x��(1),...,x��(n))=(.1)p��(x1,...,xn). 
shows that the parity (.1)p of the number of transpositions in any decomposition of �� is an invariant. 

6.1. PERMUTATIONS, SIGNATURE OF A PERMUTATION 
In view of Proposition 6.2, the following de.nition makes sense: 
De.nition 6.3. For every permutation �� of [n], the parity��(��) (or sgn(��)) of the number of transpositions involved in any decomposition of �� is called the signature (or sign) of ��. 
Obviously��(�� )= .1 for every transposition �� (since (.1)1 = .1). 
A simple way to compute the signature of a permutation is to count its number of inversions. 
De.nition 6.4. Given any permutation �� on n elements, we say that a pair (i, j) of indices i, j ��{1,...,n} such that i<j and ��(i) >��(j) is an inversion of the permutation ��. 
For example, the permutation �� given by 
123456 243651 
has seven inversions 
(1, 6), (2, 3), (2, 6), (3, 6), (4, 5), (4, 6), (5, 6). 
Proposition 6.3. The signature��(��) of any permutation �� is equal to the parity (.1)I(��) of the number I(��) of inversions in ��. 
Proof. In the product 
��(x��(1),...,x��(n))= (x��(i) . x��(j)), 1��i<j��n 
the terms x��(i) . x��(j) for which ��(i) <��(j) occur in ��(x1,...,xn), whereas the terms x��(i) . x��(j) for which ��(i) >��(j) occur in ��(x1,...,xn) with a minus sign. Therefore, the number �� of terms in ��(x��(1),...,x��(n)) whose sign is the opposite of a term in ��(x1,...,xn), is equal to the number I(��) of inversions in ��, which implies that 
��(x��(1),...,x��(n))=(.1)I(��)��(x1,...,xn). 
By Proposition 6.2, the sign of (.1)I(��) is equal to the signature of ��. 
For example, the permutation 
123456 243651 
has odd signature since it has seven inversions and (.1)7 = .1. 
Remark: When �� = idn is the identity permutation, since we agreed that the composition of 0 transpositions is the identity, it it still correct that (.1)0 =��(id) = +1. From Proposition 
6.2, it is immediate that��(��f . ��)=f)��(��). In particular, since ��.1 . �� idn, we get
��(��= ��(��.1)=��(��). We can now proceed with the de.nition of determinants. 


6.2 Alternating Multilinear Maps 
First we de.ne multilinear maps, symmetric multilinear maps, and alternating multilinear maps. 
Remark: Most of the de.nitions and results presented in this section also hold when K is a commutative ring and when we consider modules over K (free modules, when bases are needed). 
Let E1,...,En, and F , be vector spaces over a .eld K, where n �� 1. 
De.nition 6.5. A function f : E1 �� ... �� En �� F is a multilinear map (or an n-linear map) if it is linear in each argument, holding the others .xed. More explicitly, for every i, 1 �� i �� n, for all x1 �� E1,..., xi.1 �� Ei.1, xi+1 �� Ei+1,..., xn �� En, for all x, y �� Ei, for all �� �� K, 
f(x1,...,xi.1,x + y, xi+1,...,xn)= f(x1,...,xi.1, x, xi+1,...,xn) 
+ f(x1,...,xi.1, y, xi+1,...,xn), f(x1,...,xi.1, ��x, xi+1,...,xn)= ��f(x1,...,xi.1, x, xi+1,...,xn). 
When F = K, we call f an n-linear form (or multilinear form). If n �� 2 and E1 = E2 = ... = En, an n-linear map f : E �� ... �� E �� F is called symmetric, if f(x1,...,xn)= f(x��(1),...,x��(n)) for every permutation �� on {1,...,n}. An n-linear map f : E��...��E �� F is called alternating, if f(x1,...,xn) = 0 whenever xi = xi+1 for some i,1 �� i �� n . 1 (in other words, when two adjacent arguments are equal). It does no harm to agree that when n = 1, a linear map is considered to be both symmetric and alternating, and we will do so. 
When n = 2, a 2-linear map f : E1 �� E2 �� F is called a bilinear map. We have already seen several examples of bilinear maps. Multiplication ��: K �� K �� K is a bilinear map, treating K as a vector space over itself. 
The operation (.,.��: E. ��E �� K applying a linear form to a vector is a bilinear map. 
Symmetric bilinear maps (and multilinear maps) play an important role in geometry (inner products, quadratic forms) and in di.erential calculus (partial derivatives). 
A bilinear map is symmetric if f(u, v)= f(v, u), for all u, v �� E. 
Alternating multilinear maps satisfy the following simple but crucial properties. 
Proposition 6.4. Let f : E ��... ��E �� F be an n-linear alternating map, with n �� 2. The following properties hold: 
(1) 


f(...,xi,xi+1,...)= .f(...,xi+1,xi,...) 
(2) 

f(...,xi,...,xj,...)=0, 
where xi = xj, and 1 �� i<j �� n. 

6.2. ALTERNATING MULTILINEAR MAPS 
(3) 


f(...,xi,...,xj,...)= .f(...,xj,...,xi,...), 
where 1 �� i<j �� n. 

(4) 

f(...,xi,...)= f(...,xi + ��xj,...), for any �� �� K, and where i = j. Proof. (1) By multilinearity applied twice, we have 
f(...,xi + xi+1,xi + xi+1,...)= f(...,xi,xi,...)+ f(...,xi,xi+1,...) 
+ f(...,xi+1,xi,...)+ f(...,xi+1,xi+1,...), and since f is alternating, this yields 
0= f(...,xi,xi+1,...)+ f(...,xi+1,xi,...), that is, f(...,xi,xi+1,...)= .f(...,xi+1,xi,...). 
(2) If xi = xj and i and j are not adjacent, we can interchange xi and xi+1, and then xi and xi+2, etc, until xi and xj become adjacent. By (1), 

f(...,xi,...,xj,...)=��f(...,xi,xj,...), where�� = +1 or .1, but f(...,xi,xj,...) = 0, since xi = xj, and (2) holds. 

(3) follows from (2) as in (1). (4) is an immediate consequence of (2). Proposition 6.4 will now be used to show a fundamental property of alternating multilin-ear maps. First we need to extend the matrix notation a little bit. Let E be a vector space over K. Given an n �� n matrix A =(aij) over K, we can de.ne a map L(A): En �� En as follows: L(A)1(u)= a11u1 + ������ + a1 nun, ... L(A)n(u)= an 1u1 + ������ + annun, for all u1,...,un �� E and with u =(u1,...,un). It is immediately veri.ed that L(A) is linear. Then given two n��n matrices A =(aij) and B =(bij), by repeating the calculations establishing the product of matrices (just before De.nition 2.14), we can show that 


L(AB)= L(A) . L(B). It is then convenient to use the matrix notation to describe the e.ect of the linear map L(A), 
as

.
..
..
. 
L(A)1(u) a11 a12 ... a1 n u1 
.... 

L(A)2(u) 
... 

.... .... 
= 
.... 
.... 

u2 
... 

.... 

a21 a22 ... a2 n 
... 
... ... 
.

.
.
. 
L(A)n(u) an 1 an 2 ... ann un 
Lemma 6.5. Let f : E �� ... �� E �� F be an n-linear alternating map. Let (u1,...,un) and (v1,...,vn) be two families of n vectors, such that, 
v1 = a11u1 + ������ + an 1un, 
... 
vn = a1 nu1 + ������ + annun. 

Equivalently, letting 
.
. 
a1 1  a1 2  . . .  a1 n  
A = .... a2 1 . . .  a2 2 . . .  . . . ...  a2 n . . .  
an 1  an 2  . . .  an n  
assume that we have  
. . . .  

...,. 

.... 

.... .... 
v1 u1 
.vu.22 .. 
= AT .
.. 
.. 
.. 
vn un 
Then, 
  
 

f(v1,...,vn)=��(��)a��(1) 1 ������ a��(n) nf(u1,...,un), 
�С�Sn
where the sum ranges over all permutations �� on {1,...,n}. 
Proof. Expanding f(v1,...,vn) by multilinearity, we get a sum of terms of the form 
a��(1) 1 ������ a��(n) nf(u��(1),...,u��(n)), 
for all possible functions �� : {1,...,n}��{1,...,n}. However, because f is alternating, only the terms for which �� is a permutation are nonzero. By Proposition 6.1, every permutation �� is a product of transpositions, and by Proposition 6.2, the parity��(��) of the number of transpositions only depends on ��. Then applying Proposition 6.4 (3) to each transposition in ��, we get 
a��(1) 1 ������ a��(n) nf(u��(1),...,u��(n))=��(��)a��(1) 1 ������ a��(n) nf(u1,...,un). 
Thus, we get the expression of the lemma. 

6.2. ALTERNATING MULTILINEAR MAPS 
For the case of n = 2, the proof details of Lemma 6.5 become 
f(v1,v2)= f(a11u1 + a21u2,a12u1 + a22u2) = f(a11u1 + a21u2,a12u1)+ f(a11u1 + a21u2,a22u2) = f(a11u1,a12u1)+ f(a21u2,a12u1) 
+ f(a11ua,a22u2)+ f(a21u2,a22u2) = a11a12f(u1,u1)+ a21a12f(u2,u1)+ a11a22f(u1,u2) 
+ a21a22f(u2,u2) = a21a12f(u2,u1)a11a22f(u1,u2) =(a11a22 . a12a22) f(u1,u2). 
Hopefully the reader will recognize the quantity a11a22 . a12a22. It is the determinant of the 
2 �� 2 matrix  
A =  a11  a12  .  
a21  a22  
This is no accident. The quantity  
det(A) =  ��(��)a��(1) 1 �� �� �� a��(n) n  

�С�Sn
is in fact the value of the determinant of A (which, as we shall see shortly, is also equal to the determinant of AT). However, working directly with the above de.nition is quite awkward, and we will proceed via a slightly indirect route 
Remark: The reader might have been puzzled by the fact that it is the transpose matrix AT rather than A itself that appears in Lemma 6.5. The reason is that if we want the generic term in the determinant to be
��(��)a��(1) 1 ������ a��(n) n, 
where the permutation applies to the .rst index, then we have to express the vjs in terms of the uis in terms of AT as we did. Furthermore, since 
vj = a1 ju1 + ������ + aijui + ������ + anjun, 
we see that vj corresponds to the jth column of the matrix A, and so the determinant is viewed as a function of the columns of A. 
The literature is split on this point. Some authors prefer to de.ne a determinant as we did. Others use A itself, which amounts to viewing det as a function of the rows, in which case we get the expression 
��(��)a1 ��(1) ������ an��(n). �ҡ�Sn
Corollary 6.8 show that these two expressions are equal, so it doesn��t matter which is chosen. This is a matter of taste. 
6.3 De.nition of a Determinant 
Recall that the set of all square n �� n-matrices with coe.cients in a .eld K is denoted by Mn(K). 
De.nition 6.6. A determinant is de.ned as any map 
D :Mn(K) �� K, 
which, when viewed as a map on (Kn)n, i.e., a map of the n columns of a matrix, is n-linear alternating and such that D(In) = 1 for the identity matrix In. Equivalently, we can consider a vector space E of dimension n, some .xed basis (e1,...,en), and de.ne 
D : En �� K 
as an n-linear alternating map such that D(e1,...,en) = 1. 
First we will show that such maps D exist, using an inductive de.nition that also gives a recursive method for computing determinants. Actually, we will de.ne a family (Dn)n��1 of (.nite) sets of maps D :Mn(K) �� K. Second we will show that determinants are in fact uniquely de.ned, that is, we will show that each Dn consists of a single map. This will show the equivalence of the direct de.nition det(A) of Lemma 6.5 with the inductive de.nition D(A). Finally, we will prove some basic properties of determinants, using the uniqueness theorem. 
Given a matrix A �� Mn(K), we denote its n columns by A1,...,An . In order to describe the recursive process to de.ne a determinant we need the notion of a minor. 
De.nition 6.7. Given any n��n matrix with n �� 2, for any two indices i, j with 1 �� i, j �� n, let Aij be the (n . 1) �� (n . 1) matrix obtained by deleting Row i and Column j from A 
and called a minor:

.
. 
......... 

�� �� �������������� 
�� 

�� 
�� 
�� 

......... 

Aij 
= 

. 

For example, if 

.
. 
..... 

2 .10 0 0 
.12 .10 0 
0 .12 .10 

00 .12 .1 
000 .12 

..... 

A = 

6.3. DEFINITION OF A DETERMINANT 
then

.
. 
2 .10 0 
0 .1 .10 

A23 = 
... 

...

. 

00 2 .1 
00 .12 

De.nition 6.8. For every n �� 1, we de.ne a .nite set Dn of maps D :Mn(K) �� K inductively as follows: When n = 1, D1 consists of the single map D such that, D(A)= a, where A =(a), with a �� K. Assume that Dn.1 has been de.ned, where n �� 2. Then Dn consists of all the maps D such that, for some i,1 �� i �� n, 
D(A)=(.1)i+1 ai 1D(Ai 1)+ ������ +(.1)i+n ainD(Ain), 
where for every j,1 �� j �� n, D(Aij) is the result of applying any D in Dn.1 to the minor Aij. We confess that the use of the same letter D for the member of Dn being de.ned, and for members of Dn.1, may be slightly confusing. We considered using subscripts to distinguish, but this seems to complicate things unnecessarily. One should not worry too much anyway, since it will turn out that each Dn contains just one map. 
Each (.1)i+jD(Aij) is called the cofactor of aij, and the inductive expression for D(A) is called a Laplace expansion of D according to the i-th Row. Given a matrix A �� Mn(K), each D(A) is called a determinant of A. 
We can think of each member of Dn as an algorithm to evaluate ��the�� determinant of A. The main point is that these algorithms, which recursively evaluate a determinant using all possible Laplace row expansions, all yield the same result, det(A). 
We will prove shortly that D(A) is uniquely de.ned (at the moment, it is not clear that Dn consists of a single map). Assuming this fact, given a n �� n-matrix A =(aij), 
.
. 
A = 

.... 

a11 a12 ... a1 n a21 a22 ... a2 n 
.. .
.
. ...
.
.. . 
an 1 an 2 ... ann 
....

, 

its determinant is denoted by D(A) or det(A), or more explicitly by 

det(A)=

         

a11 a12 ... a1 n a21 a22 ... a2 n 
.. .
.
. ...
.
.. . 
an 1 an 2 ... ann
         

. 

Let us .rst consider some examples. 

Example 6.1. 
1. When n = 2, if 
ab 

A = , 
cd 
then by expanding according to any row, we have 
D(A)= ad . bc. 
2. When n = 3, if 
.. 
a11 a12 a13 
..
A = ,
a21 a22 a23 
a31 a32 a33 
then by expanding according to the .rst row, we have 
a22 a23 a21 a23 a21 a22 
D(A)= a11 . a12 + a13 , a32 a33 a31 a33 a31 a32 
that is, 
D(A)= a11(a22a33 . a32a23) . a12(a21a33 . a31a23) 
+ a13(a21a32 . a31a22), 
which gives the explicit formula 
D(A)= a11a22a33 + a21a32a13 + a31a12a23 
. a11a32a23 . a21a12a33 . a31a22a13. 
We now show that each D ��Dn is a determinant (map). 
Lemma 6.6. For every n �� 1, for every D ��Dn as de.ned in De.nition 6.8, D is an alternating multilinear map such that D(In)=1. 
Proof. By induction on n, it is obvious that D(In) = 1. Let us now prove that D is multilinear. Let us show that D is linear in each column. Consider any Column k. Since 
D(A)=(.1)i+1 ai 1D(Ai 1)+ ������ +(.1)i+jaijD(Aij)+ ������ 
+(.1)i+n ainD(Ain), 
if j = k, then by induction, D(Aij) is linear in Column k, and aij does not belong to Column k, so (.1)i+jaijD(Aij) is linear in Column k. If j = k, then D(Aij) does not depend on Column k = j, since Aij is obtained from A by deleting Row i and Column j = k, and aij 
6.3. DEFINITION OF A DETERMINANT 
belongs to Column j = k. Thus, (.1)i+jaijD(Aij) is linear in Column k. Consequently, in all cases, (.1)i+jaijD(Aij) is linear in Column k, and thus, D(A) is linear in Column k. 
Let us now prove that D is alternating. Assume that two adjacent columns of A are 
Ak+1
equal, say Ak = . Assume that j = k and j = k + 1. Then the matrix Aij has two identical adjacent columns, and by the induction hypothesis, D(Aij) = 0. The remaining terms of D(A) are 
(.1)i+k aikD(Aik)+(.1)i+k+1 aik+1D(Aik+1). 
However, the two matrices Aik and Aik+1 are equal, since we are assuming that Columns k and k +1 of A are identical and Aik is obtained from A by deleting Row i and Column k while Aik+1 is obtained from A by deleting Row i and Column k + 1. Similarly, aik = aik+1, since Columns k and k +1 of A are equal. But then, 
(.1)i+k aikD(Aik)+(.1)i+k+1 aik+1D(Aik+1) 
=(.1)i+k aikD(Aik) . (.1)i+k aikD(Aik)=0. 
This shows that D is alternating and completes the proof. 
Lemma 6.6 shows the existence of determinants. We now prove their uniqueness. 
Theorem 6.7. For every n �� 1, for every D ��Dn, for every matrix A �� Mn(K), we have 
D(A)= ��(��)a��(1) 1 ������ a��(n) n, 
�С�Sn
where the sum ranges over all permutations �� on {1,...,n}. As a consequence, Dn consists of a single map for every n �� 1, and this map is given by the above explicit formula. 
Proof. Consider the standard basis (e1,...,en) of Kn, where (ei)i = 1 and (ei)j = 0, for j = i. Then each column Aj of A corresponds to a vector vj whose coordinates over the basis (e1,...,en) are the components of Aj, that is, we can write 
v1 = a11e1 + ������ + an 1en, 
... 
vn = a1 ne1 + ������ + annen. 

Since by Lemma 6.6, each D is a multilinear alternating map, by applying Lemma 6.5, we get 
D(A)= D(v1,...,vn)= ��(��)a��(1) 1 ������ a��(n) n D(e1,...,en), 
�С�Sn
where the sum ranges over all permutations �� on {1,...,n}. But D(e1,...,en)= D(In), and by Lemma 6.6, we have D(In) = 1. Thus, 
D(A)= ��(��)a��(1) 1 ������ a��(n) n, 
�С�Sn
where the sum ranges over all permutations �� on {1,...,n}. 

From now on we will favor the notation det(A) over D(A) for the determinant of a square matrix. 
Remark: There is a geometric interpretation of determinants which we .nd quite illumi-nating. Given n linearly independent vectors (u1,...,un) in Rn, the set 
Pn = {��1u1 + ������ + ��nun | 0 �� ��i �� 1, 1 �� i �� n} 
is called a parallelotope. If n = 2, then P2 is a parallelogram and if n = 3, then P3 is a parallelepiped, a skew box having u1,u2,u3 as three of its corner sides. See Figures 6.1 and 
6.2. 

Figure 6.1: The parallelogram in Rw spanned by the vectors u1 = (1, 0) and u2 = (1, 1). 
Then it turns out that det(u1,...,un) is the signed volume of the parallelotope Pn (where volume means n-dimensional volume). The sign of this volume accounts for the orientation of Pn in Rn . 
We can now prove some properties of determinants. 
Corollary 6.8. For every matrix A �� Mn(K), we have det(A) = det(AT). 
Proof. By Theorem 6.7, we have 
det(A)= ��(��)a��(1) 1 ������ a��(n) n, 
�С�Sn
where the sum ranges over all permutations �� on {1,...,n}. Since a permutation is invertible, every product 
a��(1) 1 ������ a��(n) n 
can be rewritten as 
a1 ��.1(1) ������ an��.1(n), 
6.3. DEFINITION OF A DETERMINANT 

Figure 6.2: The parallelepiped in R3 spanned by the vectors u1 = (1, 1, 0), u2 = (0, 1, 0), and u3 = (0, 0, 1). 
and since��(��.1)=��(��) and the sum is taken over all permutations on {1,...,n}, we have 
��(��)a��(1) 1 ������ a��(n) n = ��(��)a1 ��(1) ������ an��(n), �С�Sn�ҡ�Sn
where �� and �� range over all permutations. But it is immediately veri.ed that 
det(AT)= ��(��)a1 ��(1) ������ an��(n). 
�ҡ�Sn
A useful consequence of Corollary 6.8 is that the determinant of a matrix is also a multi-linear alternating map of its rows. This fact, combined with the fact that the determinant of a matrix is a multilinear alternating map of its columns, is often useful for .nding short-cuts in computing determinants. We illustrate this point on the following example which shows up in polynomial interpolation. 
Example 6.2. Consider the so-called Vandermonde determinant 
1  1  . . .  1  
x1  x2  . . .  xn  
V (x1, . . . , xn) =  x2 1  x2 2  . . .  x2 n  .  

.. .
.
. ...
.
.. . 
n.1 n.1 n.1
x x ... x
12 n 
We claim that V (x1,...,xn)= (xj . xi), 1��i<j��n 
with V (x1,...,xn) = 1, when n = 1. We prove it by induction on n �� 1. The case n =1 is obvious. Assume n �� 2. We proceed as follows: multiply Row n . 1 by x1 and subtract it from Row n (the last row), then multiply Row n . 2 by x1 and subtract it from Row n . 1, etc, multiply Row i . 1 by x1 and subtract it from row i, until we reach Row 1. We obtain the following determinant: 
1  1  . . .  1  
0  x2 . x1  . . .  xn . x1  
V (x1, . . . , xn) =  0  x2(x2 . x1)  . . .  xn(xn . x1)  .  
. . .  . . .  ...  . . .  
0  x n.2 2 (x2 . x1)  . . .  xn.2 n (xn . x1)  

Now expanding this determinant according to the .rst column and using multilinearity, we can factor (xi . x1) from the column of index i . 1 of the matrix obtained by deleting the .rst row and the .rst column, and thus 
V (x1,...,xn)=(x2 . x1)(x3 . x1) ������ (xn . x1)V (x2,...,xn), 
which establishes the induction step. 
Remark: Observe that ��(x1,...,xn)= V (xn,...,x1)=(.1)(n)V (x1,...xn),
2
where ��(x1,...,xn) is the discriminant of (x1,...,xn) introduced in De.nition 6.2. Lemma 6.5 can be reformulated nicely as follows. 
Proposition 6.9. Let f : E �� ... �� E �� F be an n-linear alternating map. Let (u1,...,un) and (v1,...,vn) be two families of n vectors, such that 
v1 = a11u1 + ������ + a1 nun, 
... 
vn = an 1u1 + ������ + annun. 

Equivalently, letting 

.
. 
.... 

a11 a12 ... a1 n a21 a22 ... a2 n 
.. .
.
. ...
.
.. . 
....

A = 

, 

an 1 an 2 ... ann 
6.3. DEFINITION OF A DETERMINANT 
assume that we have

..
.. 
v1 u1 
v2 u2 
....

. 

. 

. 

.... 

= A 

....

. 

. 

. 

.... 

. 

vn un 
Then, 
f(v1,...,vn) = det(A)f(u1,...,un). 
Proof. The only di.erence with Lemma 6.5 is that here we are using AT instead of A. Thus, by Lemma 6.5 and Corollary 6.8, we get the desired result. 
As a consequence, we get the very useful property that the determinant of a product of matrices is the product of the determinants of these matrices. Proposition 6.10. For any two n��n-matrices A and B, we have det(AB) = det(A) det(B). Proof. We use Proposition 6.9 as follows: let (e1,...,en) be the standard basis of Kn, and 
let

..
.. 
w1 e1 
.... 

w2 
. 

. 

. 

.... 

= AB 

.... 

e2 
. 

. 

. 

.... 

. 

wn en 
Then we get 
det(w1,...,wn) = det(AB) det(e1,...,en) = det(AB), since det(e1,...,en) = 1. Now letting 
..
.. 
v1 e1 
.... 

v2 
. 

. 

. 

.... 

= B 

.... 

e2 
. 

. 

. 

....

, 

vn en 
we get 
det(v1,...,vn) = det(B), 
and since

..
.. 
w1 v1 
.... 

w2 
. 

. 

. 

.... 

= A 

.... 

v2 
. 

. 

. 

....

, 

wn vn 
we get det(w1,...,wn) = det(A) det(v1,...,vn) = det(A) det(B). 
It should be noted that all the results of this section, up to now, also hold when K is a commutative ring and not necessarily a .eld. We can now characterize when an n��n-matrix A is invertible in terms of its determinant det(A). 
6.4 Inverse Matrices and Determinants 
In the next two sections, K is a commutative ring and when needed a .eld. De.nition 6.9. Let K be a commutative ring. Given a matrix A �� Mn(K), let A 
=(bij) be the matrix de.ned such that bij =(.1)i+j det(Aji), 
the cofactor of aji. The matrix A 
is called the adjugate of A, and each matrix Aji is called a minor of the matrix A. 
For example, if 
.. 
11 1 
..
A =2 .2 .2 , 33 .3 we have .2 .2 11 
b11 = det(A11)= =12 b12 = . det(A21)= . =6 
3 .33 .3 11 2 .2 
b13 = det(A31)= =0 b21 = . det(A12)= . =0 
.2 .23 .3 11 11 
b22 = det(A22)= = .6 b23 = . det(A32)= . =4 
3 .32 .2 2 .2 11 
b31 = det(A13)= =12 b32 = . det(A23)= . =0 
33 33 11 
b33 = det(A33)= = .4,
2 .2 
we .nd that .. 126 0 
..
A 
=0 .64 . 
12 0 .4 
Note the reversal of the indices in 
bij =(.1)i+j det(Aji). 

Thus, A 
is the transpose of the matrix of cofactors of elements of A. We have the following proposition. Proposition 6.11. Let K be a commutative ring. For every matrix A �� Mn(K), we have AA 
= AA 
= det(A)In. As a consequence, A is invertible i. det(A) is invertible, and if so, A.1 = (det(A)).1A
. 
6.4. INVERSE MATRICES AND DETERMINANTS 
Proof. If AA=(bij) and AAA=(cij), we know that the entry cij in row i and column j of AAAis 
cij = ai 1b1 j + ������ + aikbkj + ������ + ainbnj, 
which is equal to 
ai 1(.1)j+1 det(Aj 1)+ ������ + ain(.1)j+n det(Ajn). 
If j = i, then we recognize the expression of the expansion of det(A) according to the i-th row: 
cii = det(A)= ai 1(.1)i+1 det(Ai 1)+ ������ + ain(.1)i+n det(Ain). 
If j = i, we can form the matrix Af by replacing the j-th row of A by the i-th row of A. Now the matrix Ajk obtained by deleting row j and column k from A is equal to the matrix Af 
jk obtained by deleting row j and column k from Af, since A and Af only di.er by the j-th row. Thus, det(Ajk) = det(Af ),
jk
and we have cij = ai 1(.1)j+1 det(Af )+ ������ + ain(.1)j+n det(Af ).
j 1jn
However, this is the expansion of det(Af) according to the j-th row, since the j-th row of Af is equal to the i-th row of A. Furthermore, since Af has two identical rows i and j, because det is an alternating map of the rows (see an earlier remark), we have det(Af) = 0. Thus, we have shown that cii = det(A), and cij = 0, when j = i, and so 
AAA= det(A)In. 
A
It is also obvious from the de.nition of A, that 
A
AAT = AT . 
Then applying the .rst part of the argument to AT, we have 
ATAAT = det(AT)In, 
A
and since det(AT) = det(A), AAT AT, and ( A= ATAAT, we get 
= AA)T 
AT AT =( A
det(A)In = AT A= AT AAA)T , 
that is, ( A
AA)T = det(A)In, 
which yields 
A
AA = det(A)In, since In T = In. This proves that 
AA = det(A)In. 
A
AAA=

As a consequence, if det(A) is invertible, we have A.1 = (det(A)).1AA. Conversely, if A is invertible, from AA.1 = In, by Proposition 6.10, we have det(A) det(A.1) = 1, and det(A) is invertible. 
For example, we saw earlier that 
.
.
.
. 
111 1260 

A = 33 .3 120 .4 
A
A =

.

.

.

.

.2 .2 

.6

2 

and

0 

4 

, 

and we have

.
..
..
. 
111 1260 100 

.

2 

.2 .2

.
.

0 

.6 

4 

. = 24 

.

010

. 

33 .3 120 .4 001 
with det(A) = 24. 
When K is a .eld, an element a �� K is invertible i. a = 0. In this case, the second part of the proposition can be stated as A is invertible i. det(A) = 0. Note in passing that this method of computing the inverse of a matrix is usually not practical. 
6.5 Systems of Linear Equations and Determinants 
We now consider some applications of determinants to linear independence and to solving systems of linear equations. Although these results hold for matrices over certain rings, their proofs require more sophisticated methods. Therefore, we assume again that K is a .eld (usually, K = R or K = C). 
Let A be an n �� n-matrix, x a column vectors of variables, and b another column vector, and let A1,...,An denote the columns of A. Observe that the system of equations Ax = b, 
.
...
.
. 
a11 a12 ... a1 n x1 b1 
.... 

.... 
.... 

x2 
. 

. 

. 

.... 

= 

.... 

b2 
. 

. 

. 

.... 

a21 a22 ... a2 n 
.. .
.
.. 

..

.

.. 

. 

an 1 an 2 ... ann xn bn 
is equivalent to x1A1 + ������ + xjAj + ������ + xnAn = b, since the equation corresponding to the i-th row is in both cases ai 1x1 + ������ + aijxj + ������ + ainxn = bi. First we characterize linear independence of the column vectors of a matrix A in terms of its determinant. 
6.5. SYSTEMS OF LINEAR EQUATIONS AND DETERMINANTS 
Proposition 6.12. Given an n �� n-matrix A over a .eld K, the columns A1,...,An of A are linearly dependent i. det(A) = det(A1,...,An)=0. Equivalently, A has rank n i. det(A)=0. 
Proof. First assume that the columns A1,...,An of A are linearly dependent. Then there are x1,...,xn �� K, such that 
x1A1 + ������ + xjAj + ������ + xnAn =0, 
where xj = 0 for some j. If we compute 
det(A1 ,...,x1A1 + ������ + xjAj + ������ + xnAn,...,An) = det(A1 ,..., 0,...,An)=0, 
where 0 occurs in the j-th position. By multilinearity, all terms containing two identical columns Ak for k = j vanish, and we get 
det(A1 ,...,x1A1 + ������ + xjAj + ������ + xnAn,...,An)= xj det(A1,...,An)=0. 
Since xj = 0 and K is a .eld, we must have det(A1,...,An) = 0. 
Conversely, we show that if the columns A1,...,An of A are linearly independent, then det(A1,...,An) = 0. If the columns A1,...,An of A are linearly independent, then they form a basis of Kn, and we can express the standard basis (e1,...,en) of Kn in terms of 
A1,...,An . 
.
..
.. 
Thus,wehave . 
A1
e1 b11 b12 ... b1 n 
.... 

e2 
. 

. 

. 

.... 

= 

.... 

.... 
.... 

....

A2 
. 

. 

. 
An b21 b22 ... b2 n 
,

.. .
.
.. 

..

.

.. 

. 

en bn 1 bn 2 ... bnn for some matrix B =(bij), and by Proposition 6.9, we get 
det(e1,...,en) = det(B) det(A1,...,An), 
and since det(e1,...,en) = 1, this implies that det(A1,...,An) = 0 (and det(B) = 0). For the second assertion, recall that the rank of a matrix is equal to the maximum number of linearly independent columns, and the conclusion is clear. 
We now characterize when a system of linear equations of the form Ax = b has a unique solution. 
Proposition 6.13. Given an n �� n-matrix A over a .eld K, the following properties hold: 
(1) 
For every column vector b, there is a unique column vector x such that Ax = b i. the only solution to Ax =0 is the trivial vector x =0, i. det(A)=0. 

(2) 
If det(A)=0, the unique solution of Ax = b is given by the expressions 


det(A1,...,Aj.1, b, Aj+1,...,An) 
xj = ,
det(A1,...,Aj.1,Aj,Aj+1,...,An)
known as Cramer��s rules. 
(3) The system of linear equations Ax =0 has a nonzero solution i. det(A)=0. 
Proof. (1) Assume that Ax = b has a single solution x0, and assume that Ay = 0 with y = 0. Then, 
A(x0 + y)= Ax0 + Ay = Ax0 +0= b, 
and x0 + y = x0 is another solution of Ax = b, contradicting the hypothesis that Ax = b has a single solution x0. Thus, Ax = 0 only has the trivial solution. Now assume that Ax =0 only has the trivial solution. This means that the columns A1,...,An of A are linearly independent, and by Proposition 6.12, we have det(A) = 0. Finally, if det(A)= 0,by Proposition 6.11, this means that A is invertible, and then for every b, Ax = b is equivalent to x = A.1b, which shows that Ax = b has a single solution. 
(2) Assume that Ax = b. If we compute 
det(A1 ,...,x1A1 + ������ + xjAj + ������ + xnAn,...,An) = det(A1, . . . , b, . . . , An), 
where b occurs in the j-th position, by multilinearity, all terms containing two identical columns Ak for k = j vanish, and we get 
xj det(A1,...,An) = det(A1,...,Aj.1, b, Aj+1,...,An), 
for every j,1 �� j �� n. Since we assumed that det(A) = det(A1,...,An) = 0, we get the desired expression. 
(3) Note that Ax = 0 has a nonzero solution i. A1,...,An are linearly dependent (as observed in the proof of Proposition 6.12), which, by Proposition 6.12, is equivalent to det(A) = 0. 
As pleasing as Cramer��s rules are, it is usually impractical to solve systems of linear equations using the above expressions. However, these formula imply an interesting fact, which is that the solution of the system Ax = b are continuous in A and b. If we assume that the entries in A are continuous functions aij(t) and the entries in b are are also continuous functions bj(t) of a real parameter t, since determinants are polynomial functions of their entries, the expressions 
det(A1,...,Aj.1, b, Aj+1,...,An) 
xj(t)= 
det(A1,...,Aj.1,Aj,Aj+1,...,An) 
are ratios of polynomials, and thus are also continuous as long as det(A(t)) is nonzero. Similarly, if the functions aij(t) and bj(t) are di.erentiable, so are the xj(t). 
6.6. DETERMINANT OF A LINEAR MAP 
6.6 Determinant of a Linear Map 
Given a vector space E of .nite dimension n, given a basis (u1,...,un) of E, for every linear map f : E �� E, if M(f) is the matrix of f w.r.t. the basis (u1,...,un), we can de.ne det(f) = det(M(f)). If(v1,...,vn) is any other basis of E, and if P is the change of basis matrix, by Corollary 3.5, the matrix of f with respect to the basis (v1,...,vn) is P .1M(f)P . By Proposition 6.10, we have 
det(P .1M(f)P ) = det(P .1) det(M(f)) det(P )= det(P .1) det(P ) det(M(f)) = det(M(f)). 
Thus, det(f) is indeed independent of the basis of E. 
De.nition 6.10. Given a vector space E of .nite dimension, for any linear map f : E �� E, we de.ne the determinant det(f) of f as the determinant det(M(f)) of the matrix of f in any basis (since, from the discussion just before this de.nition, this determinant does not depend on the basis). 
Then we have the following proposition. 
Proposition 6.14. Given any vector space E of .nite dimension n, a linear map f : E �� E is invertible i. det(f)=0. 
Proof. The linear map f : E �� E is invertible i. its matrix M(f) in any basis is invertible (by Proposition 3.2), i. det(M(f)) = 0, by Proposition 6.11. 
Given a vector space of .nite dimension n, it is easily seen that the set of bijective linear maps f : E �� E such that det(f) = 1 is a group under composition. This group is a subgroup of the general linear group GL(E). It is called the special linear group (of E), and it is denoted by SL(E), or when E = Kn, by SL(n, K), or even by SL(n). 
6.7 The Cayley�CHamilton Theorem 
We next discuss an interesting and important application of Proposition 6.11, the Cayley�C Hamilton theorem. The results of this section apply to matrices over any commutative ring 
K. First we need the concept of the characteristic polynomial of a matrix. 
De.nition 6.11. If K is any commutative ring, for every n �� n matrix A �� Mn(K), the characteristic polynomial PA(X) of A is the determinant 
PA(X) = det(XI . A). 
The characteristic polynomial PA(X) is a polynomial in K[X], the ring of polynomials in the indeterminate X with coe.cients in the ring K. For example, when n = 2, if 
ab 
A = , 
cd 
then 
X . a .b 
PA(X)= = X2 . (a + d)X + ad . bc. 
.cX . d 
We can substitute the matrix A for the variable X in the polynomial PA(X), obtaining a matrix PA. If we write 
PA(X)= Xn + c1Xn.1 + ������ + cn, 
then PA = An + c1An.1 + ������ + cnI. 
We have the following remarkable theorem. 
Theorem 6.15. (Cayley�CHamilton) If K is any commutative ring, for every n �� n matrix A �� Mn(K), if we let 
PA(X)= Xn + c1Xn.1 + ������ + cn be the characteristic polynomial of A, then 
PA = An + c1An.1 + ������ + cnI =0. Proof. We can view the matrix B = XI . A as a matrix with coe.cients in the polynomial 
A
ring K[X], and then we can form the matrix B which is the transpose of the matrix of cofactors of elements of B. Each entry in BAisan (n . 1) �� (n . 1) determinant, and thus a polynomial of degree a most n . 1, so we can write BAas 
BA= Xn.1B0 + Xn.2B1 + ������ + Bn.1, 
for some n �� n matrices B0,...,Bn.1 with coe.cients in K. For example, when n = 2, we have 
X . a .bX . db 10 .db
A
B = ,B == X + .
.cX . d cX . a 01 c .a By Proposition 6.11, we have BBA= det(B)I = PA(X)I. On the other hand, we have BBA=(XI . A)(Xn.1B0 + Xn.2B1 + ������ + Xn.j.1Bj + ������ + Bn.1), and by multiplying out the right-hand side, we get BBA= XnD0 + Xn.1D1 + ������ + Xn.jDj + ������ + Dn, 
6.7. THE CAYLEY�CHAMILTON THEOREM 
with 

D0 = B0 D1 = B1 . AB0 . 
. 
. Dj = Bj . ABj.1 . 
. 
. 
Dn.1 = Bn.1 . ABn.2 Dn = .ABn.1. 
Since PA(X)I =(Xn + c1Xn.1 + ������ + cn)I, 
the equality XnD0 + Xn.1D1 + ������ + Dn =(Xn + c1Xn.1 + ������ + cn)I 
is an equality between two matrices, so it requires that all corresponding entries are equal, and since these are polynomials, the coe.cients of these polynomials must be identical, which is equivalent to the set of equations 
I = B0 c1I = B1 . AB0 . 
. 
. cjI = Bj . ABj.1 . 
. 
. 
cn.1I = Bn.1 . ABn.2 cnI = .ABn.1, 
for all j, with 1 �� j �� n . 1. If, as in the table below, 
An = AnB0 
c1An.1 = An.1(B1 . AB0) . 
. 
. 
cjAn.j = An.j(Bj . ABj.1) . 
. 
. cn.1A = A(Bn.1 . ABn.2) cnI = .ABn.1, 
we multiply the .rst equation by An, the last by I, and generally the (j + 1)th by An.j, when we add up all these new equations, we see that the right-hand side adds up to 0, and we get our desired equation 
An + c1An.1 + ������ + cnI =0, 
as claimed. 
As a concrete example, when n = 2, the matrix 
a  b  
A =  
c  d  

satis.es the equation A2 . (a + d)A +(ad . bc)I =0. 
Most readers will probably .nd the proof of Theorem 6.15 rather clever but very myste-rious and unmotivated. The conceptual di.culty is that we really need to understand how polynomials in one variable ��act�� on vectors in terms of the matrix A. This can be done and yields a more ��natural�� proof. Actually, the reasoning is simpler and more general if we free ourselves from matrices and instead consider a .nite-dimensional vector space E and some given linear map f : E �� E. Given any polynomial p(X)= a0Xn + a1Xn.1 + ������ + an with coe.cients in the .eld K, we de.ne the linear map p(f): E �� E by 
p(f)= a0fn + a1fn.1 + ������ + anid, 
where fk = f .�� ����. f, the k-fold composition of f with itself. Note that 
p(f)(u)= a0fn(u)+ a1fn.1(u)+ ������ + anu, 
for every vector u �� E. Then we de.ne a new kind of scalar multiplication ��: K[X] ��E �� E by polynomials as follows: for every polynomial p(X) �� K[X], for every u �� E, 
p(X) �� u = p(f)(u). 
It is easy to verify that this is a ��good action,�� which means that 
p �� (u + v)= p �� u + p �� v (p + q) �� u = p �� u + q �� u (pq) �� u = p �� (q �� u) 
1 �� u = u, 
for all p, q �� K[X] and all u, v �� E. With this new scalar multiplication, E is a K[X]-module. If p = �� is just a scalar in K (a polynomial of degree 0), then 
�� �� u =(��id)(u)= ��u, 
6.7. THE CAYLEY�CHAMILTON THEOREM 
which means that K acts on E by scalar multiplication as before. If p(X)= X (the monomial X), then 
X �� u = f(u). 
Now if we pick a basis (e1,...,en) of E, if a polynomial p(X) �� K[X] has the property that p(X) �� ei =0,i =1, . . . , n, 
then this means that p(f)(ei)=0 for i =1,...,n, which means that the linear map p(f) vanishes on E. We can also check, as we did in Section 6.2, that if A and B are two n �� n matrices and if (u1,...,un) are any n vectors, then 
..
..
.
. 
u1 u1 
A �� 

..

B �� 

..

. 

. 

. 

..
..

=(AB) �� 

..

. 

. 

. 

..

. 

un un 
This suggests the plan of attack for our second proof of the Cayley�CHamilton theorem. For simplicity, we prove the theorem for vector spaces over a .eld. The proof goes through for a free module over a commutative ring. 
Theorem 6.16. (Cayley�CHamilton) For every .nite-dimensional vector space over a .eld K, for every linear map f : E �� E, for every basis (e1,...,en), if A is the matrix over f over the basis (e1,...,en) and if 
PA(X)= Xn + c1Xn.1 + ������ + cn 
is the characteristic polynomial of A, then 
PA(f)= fn + c1fn.1 + ������ + cnid = 0. 
Proof. Since the columns of A consist of the vector f(ej) expressed over the basis (e1,...,en), we have 
n f(ej)= aijei, 1 �� j �� n. i=1 
Using our action of K[X] on E, the above equations can be expressed as 
n X �� ej = aij �� ei, 1 �� j �� n, i=1 
which yields 
j.1 n .aij �� ei +(X . ajj) �� ej + .aij �� ei =0, 1 �� j �� n. i=1 i=j+1 
Observe that the transpose of the characteristic polynomial shows up, so the above system can be written as 
..
...
. 
X . a11 .a21 ������ .an 1 e1 0 
.... 

.a12 X . a22 ������ .an 2 
. ... 

. ... 

. ... 

.... 

�� 

.... 

e2 
. 

. 

. 

.... 

= 

.... 

0 

. 

. 

. 

.... 

. 

.a1 n .a2 n ������ X . ann en 0 If we let B = XI . AT, then as in the previous proof, if BAis the transpose of the matrix of cofactors of B, we have BB = det(B)I = det(XI . AT)I = det(XI . A)I
A
= PAI. 
But since 

.... 
e1 0 
.... 

.... 

= 

.... 

0 

. 

. 

. 
0 

....

e2 
. 

. 

. 

B �� 

, 

en 
and since BAis matrix whose entries are polynomials in K[X], it makes sense to multiply on 
the left by BAand we get 

..
..
..
..
.. 
e1 e1 e1 00 
AB �� B �� = ( ABB) �� = P. . . . . .  AB �� AI �� ;== . . . . . . . . .  
en en  en 0 0  
that is, PA �� ej = 0, which proves that PA(f) = 0, as claimed.  j  = 1, . . . , n,  

.... 

....

.... 

.... 

.... 

.... 

.... 

.... 

.... 

....

0 

0

e2 e2 e2 
If K is a .eld, then the characteristic polynomial of a linear map f : E �� E is independent of the basis (e1,...,en) chosen in E. To prove this, observe that the matrix of f over another basis will be of the form P .1AP , for some inverible matrix P , and then 
det(XI . P .1AP ) = det(XP .1IP . P .1AP ) 
= det(P .1(XI . A)P ) 
= det(P .1) det(XI . A) det(P ) 
= det(XI . A). 
Therefore, the characteristic polynomial of a linear map is intrinsic to f, and it is denoted by Pf . 
The zeros (roots) of the characteristic polynomial of a linear map f are called the eigen-values of f. They play an important role in theory and applications. We will come back to this topic later on. 
6.8. PERMANENTS 
6.8 Permanents 
Recall that the explicit formula for the determinant of an n �� n matrix is 
det(A)= ��(��)a��(1) 1 ������ a��(n) n. 
�С�Sn
If we drop the sign��(��) of every permutation from the above formula, we obtain a quantity known as the permanent: 
per(A)= a��(1) 1 ������ a��(n) n. 
�С�Sn 
Permanents and determinants were investigated as early as 1812 by Cauchy. It is clear from the above de.nition that the permanent is a multilinear symmetric form. We also have 
per(A) = per(AT), 
and the following unsigned version of the Laplace expansion formula: 
per(A)= ai 1per(Ai 1)+ ������ + aijper(Aij)+ ������ + ainper(Ain), 
for i =1,...,n. However, unlike determinants which have a clear geometric interpretation as signed volumes, permanents do not have any natural geometric interpretation. Furthermore, determinants can be evaluated e.ciently, for example using the conversion to row reduced echelon form, but computing the permanent is hard. 
Permanents turn out to have various combinatorial interpretations. One of these is in terms of perfect matchings of bipartite graphs which we now discuss. 
See De.nition 18.5 for the de.nition of an undirected graph. A bipartite (undirected) graph G =(V, E) is a graph whose set of nodes V can be partitioned into two nonempty disjoint subsets V1 and V2, such that every edge e �� E has one endpoint in V1 and one endpoint in V2. 
An example of a bipartite graph with 14 nodes is shown in Figure 6.3; its nodes are partitioned into the two sets {x1,x2,x3,x4,x5,x6,x7} and {y1,y2,y3,y4,y5,y6,y7}. 
A matching in a graph G =(V, E) (bipartite or not) is a set M of pairwise non-adjacent edges, which means that no two edges in M share a common vertex. A perfect matching is a matching such that every node in V is incident to some edge in the matching M (every node in V is an endpoint of some edge in M). Figure 6.4 shows a perfect matching (in red) in the bipartite graph G. 
Obviously, a perfect matching in a bipartite graph can exist only if its set of nodes has a partition in two blocks of equal size, say {x1,...,xm} and {y1,...,ym}. Then there is a bijection between perfect matchings and bijections �� : {x1,...,xm}��{y1,...,ym} such that ��(xi)= yj i. there is an edge between xi and yj. 
Now every bipartite graph G with a partition of its nodes into two sets of equal size as above is represented by an m �� m matrix A =(aij) such that aij = 1 i. there is an edge between xi and yj, and aij = 0 otherwise. Using the interpretation of perfect matchings as 

Figure 6.3: A bipartite graph G. 


Figure 6.4: A perfect matching in the bipartite graph G. 

6.9. SUMMARY 
bijections �� : {x1,...,xm}��{y1,...,ym}, we see that the permanent per(A) of the (0, 1)-matrix A representing the bipartite graph G counts the number of perfect matchings in G. 
In a famous paper published in 1979, Leslie Valiant proves that computing the permanent is a #P-complete problem. Such problems are suspected to be intractable. It is known that if a polynomial-time algorithm existed to solve a #P-complete problem, then we would have P = NP , which is believed to be very unlikely. 
Another combinatorial interpretation of the permanent can be given in terms of systems of distinct representatives. Given a .nite set S, let (A1,...,An) be any sequence of nonempty subsets of S (not necessarily distinct). A system of distinct representatives (for short SDR) of the sets A1,...,An is a sequence of n distinct elements (a1,...,an), with ai �� Ai for i = 1,...,n. The number of SDR��s of a sequence of sets plays an important role in combinatorics. Now, if S = {1, 2,...,n} and if we associate to any sequence (A1,...,An) of nonempty subsets of S the matrix A =(aij) de.ned such that aij =1 if j �� Ai and aij = 0 otherwise, then the permanent per(A) counts the number of SDR��s of the sets A1,...,An. 
This interpretation of permanents in terms of SDR��s can be used to prove bounds for the permanents of various classes of matrices. Interested readers are referred to van Lint and Wilson [71] (Chapters 11 and 12). In particular, a proof of a theorem known as Van der Waerden conjecture is given in Chapter 12. This theorem states that for any n �� n matrix A with nonnegative entries in which all row-sums and column-sums are 1 (doubly stochastic matrices), we have 
n! 
per(A) �� , 
nn with equality for the matrix in which all entries are equal to 1/n. 
6.9 Summary 
The main concepts and results of this chapter are listed below: 
. 	Permutations, transpositions, basics transpositions. 

. 	Every permutation can be written as a composition of permutations. 


. 	The parity of the number of transpositions involved in any decomposition of a permu-tation �� is an invariant; it is the signature��(��) of the permutation ��. 

. 	Multilinear maps (also called n-linear maps); bilinear maps. 

. 	Symmetric and alternating multilinear maps. 

. 	A basic property of alternating multilinear maps (Lemma 6.5) and the introduction of the formula expressing a determinant. 

. 	De.nition of a determinant as a multlinear alternating map D :Mn(K) �� K such that D(I) = 1. 

. 
We de.ne the set of algorithms Dn, to compute the determinant of an n �� n matrix. 

. 
Laplace expansion according to the ith row; cofactors. 

. 
We prove that the algorithms in Dn compute determinants (Lemma 6.6). 

. 
We prove that all algorithms in Dn compute the same determinant (Theorem 6.7). 

. 
We give an interpretation of determinants as signed volumes. 

. 
We prove that det(A) = det(AT). 

. 
We prove that det(AB) = det(A) det(B). 

A

. 
The adjugate A of a matrix A. 

. 
Formula for the inverse in terms of the adjugate. 

. 
A matrix A is invertible i. det(A) = 0. 

. 
Solving linear equations using Cramer��s rules. 

. 
Determinant of a linear map. 

. 
The characteristic polynomial of a matrix. 

. 
The Cayley�CHamilton theorem. 

. 
The action of the polynomial ring induced by a linear map on a vector space. 

. 
Permanents. 

. 
Permanents count the number of perfect matchings in bipartite graphs. 

. 
Computing the permanent is a #P-perfect problem (L. Valiant). 

. 
Permanents count the number of SDRs of sequences of subsets of a given set. 


6.10 Further Readings 
Thorough expositions of the material covered in Chapter 2�C5 and 6 can be found in Strang [64, 63], Lax [44], Lang [41], Artin [3], Mac Lane and Birkho. [46], Ho.man and Kunze [38], Dummit and Foote [19], Bourbaki [8, 9], Van Der Waerden [70], Serre [57], Horn and Johnson [36], and Bertin [7]. These notions of linear algebra are nicely put to use in classical geometry, see Berger [5, 6], Tisseron [67] and Dieudonn��e [17]. 
6.11. PROBLEMS 
6.11 Problems 
Problem 6.1. Prove that every transposition can be written as a product of basic transpo-sitions. 
Problem 6.2. (1) Given two vectors in R2 of coordinates (c1.a1,c2.a2) and (b1.a1,b2.a2), prove that they are linearly dependent i. 
a1 b1 c1 a2 b2 c2 =0. 
111 
(2) Given three vectors in R3 of coordinates (d1.a1,d2.a2,d3.a3), (c1.a1,c2.a2,c3.a3), and (b1 . a1,b2 . a2,b3 . a3), prove that they are linearly dependent i. 
a1  b1  c1  d1  
a2 a3  b2 b3  c2 c3  d2 d3  = 0.  
1  1  1  1  

Problem 6.3. Let A be the (m + n) �� (m + n) block matrix (over any .eld K) given by A1 A2
A = ,
0 A4 
where A1 is an m �� m matrix, A2 is an m �� n matrix, and A4 is an n �� n matrix. Prove that det(A) = det(A1) det(A4). 
Use the above result to prove that if A is an upper triangular n��n matrix, then det(A)= a11a22 ������ ann. 
Problem 6.4. Prove that if n �� 3, then 
.
. 
det 

.... 

1+ x1y1 1+ x1y2 ������ 1+ x1yn 1+ x2y1 1+ x2y2 ������ 1+ x2yn 
. ... 
. ... 
. ... 
1+ xny1 1+ xny2 ������ 1+ xnyn 

.... 

=0. 

Problem 6.5. Prove that  
1  4  9  16  
4 9  9 16  16 25  25 36  = 0.  
16  25  36  49  

Problem 6.6. Consider the n �� n symmetric matrix 
.
. 
120 0 ... 00 252 0 ... 00 025 2 ... 00 
.. ..
...
A = 

.......... 

.......... 

.

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

(1) 
Find an upper-triangular matrix R such that A = RTR. 

(2) 
Prove that det(A) = 1. 

(3) 
Consider the sequence 


p0(��)=1 p1(��)=1 . �� pk(��) = (5 . ��)pk.1(��) . 4pk.2(��)2 �� k �� n. 
Prove that det(A . ��I)= pn(��). 
Remark: It can be shown that pn(��) has n distinct (real) roots and that the roots of pk(��) separate the roots of pk+1(��). 
Problem 6.7. Let B be the n �� n matrix (n �� 3) given by 
.
. 
B = 

.......... 

1 .1 .1 .1 ������ .1 .1 1 .11 1 ������ 11 11 .11 ������ 11 11 1 .1 ������ 11 
.  .  .  .  .  .  .  
.  .  .  .  .  .  .  
.  .  .  .  .  .  .  
1  1  1  1  �� �� ��  .1  1  
1  1  1  1  �� �� ��  1  .1  

.......... 

. 

Prove that det(B)=(.1)n(n . 2)2n.1 . 
Problem 6.8. Given a .eld K (say K = R or K = C), given any two polynomials p(X),q(X) �� K[X], we says that q(X) divides p(X) (and that p(X) is a multiple of q(X)) i. there is some polynomial s(X) �� K[X] such that 
p(X)= q(X)s(X). 
6.11. PROBLEMS 
In this case we say that q(X) is a factor of p(X), and if q(X) has degree at least one, we say that q(X) is a nontrivial factor of p(X). 
Let f(X) and g(X) be two polynomials in K[X] with 
f(X)= a0Xm + a1Xm.1 + ������ + am 
of degree m �� 1 and g(X)= b0Xn + b1Xn.1 + ������ + bn 
of degree n �� 1 (with a0,b0 = 0). 
You will need the following result which you need not prove: 
Two polynomials f(X) and g(X) with deg(f)= m �� 1 and deg(g)= n �� 1 have some common nontrivial factor i. there exist two nonzero polynomials p(X) and q(X) such that 
fp = gq, 
with deg(p) �� n . 1 and deg(q) �� m . 1. 
(1) Let Pm denote the vector space of all polynomials in K[X] of degree at most m . 1, and let T : Pn ��Pm ��Pm+n be the map given by 
T (p, q)= fp + gq, p ��Pn,q ��Pm, 
where f and g are some .xed polynomials of degree m �� 1 and n �� 1. Prove that the map T is linear. 
(2) 
Prove that T is not injective i. f and g have a common nontrivial factor. 

(3) 
Prove that f and g have a nontrivial common factor i. R(f, g) = 0, where R(f, g) is 


the determinant given by  
a0 0 �� �� ��  a1 a0 �� �� ��  �� �� �� a1 �� �� ��  �� �� �� �� �� �� �� �� ��  am �� �� �� �� �� ��  0 am �� �� ��  �� �� �� 0 �� �� ��  �� �� �� �� �� �� �� �� ��  �� �� �� �� �� �� �� �� ��  �� �� �� �� �� �� �� �� ��  0 0 �� �� ��  
�� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  
�� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  
R(f, g) =  �� �� �� 0 b0 0 �� �� ��  �� �� �� �� �� �� b1 b0 �� �� ��  �� �� �� �� �� �� �� �� �� b1 �� �� ��  �� �� �� �� �� �� �� �� �� �� �� �� �� �� ��  �� �� �� �� �� �� �� �� �� �� �� �� �� �� ��  �� �� �� 0 �� �� �� �� �� �� �� �� ��  �� �� �� a0 �� �� �� �� �� �� �� �� ��  �� �� �� a1 bn �� �� �� �� �� ��  �� �� �� �� �� �� 0 bn �� �� ��  �� �� �� �� �� �� �� �� �� 0 �� �� ��  �� �� �� am 0 �� �� �� �� �� �� .  
0  �� �� ��  0  b0  b1  �� �� ��  �� �� ��  �� �� ��  �� �� ��  �� �� ��  bn  

The above determinant is called the resultant of f and g. 
Note that the matrix of the resultant is an (n + m) �� (n + m) matrix, with the .rst row (involving the ais) occurring n times, each time shifted over to the right by one column, and the (n + 1)th row (involving the bjs) occurring m times, each time shifted over to the right by one column. Hint. Express the matrix of T over some suitable basis. 
(4) 
Compute the resultant in the following three cases: 

(a) 
m = n = 1, and write f(X)= aX + b and g(X)= cX + d. 

(b) 
m = 1 and n �� 2 arbitrary. 

(c) 
f(X)= aX2 + bX + c and g(X)=2aX + b. 

(5) 
Compute the resultant of f(X)= X3 + pX + q and g(X)=3X2 + p, and 


f(X)= a0X2 + a1X + a2 g(X)= b0X2 + b1X + b2. 
In the second case, you should get 
4R(f, g) = (2a0b2 . a1b1 +2a2b0)2 . (4a0a2 . a12)(4b0b2 . b21). 
Problem 6.9. Let A, B, C, D be n �� n real or complex matrices. 
(1) Prove that if A is invertible and if AC = CA, then 
AB 

det = det(AD . CB). 
CD 
(2) 
Prove that if H is an n �� n Hadamard matrix (n �� 2), then | det(H)| = nn/2 . 

(3) 
Prove that if H is an n �� n Hadamard matrix (n �� 2), then 
HH 



det =(2n)n . 
H .H 
Problem 6.10. Compute the product of the following determinants 
a .b .c .dx .y .z .t ba .dcyx .tz cd a .bzt x .y d .cb at .zy x 
to prove the following identity (due to Euler): 
2 2 222
(a + b2 + c + d2)(x 	+ y + z + t2) 
=(ax + by + cz + dt)2 +(ay . bx + ct . dz)2 

+(az . bt . cx + dy)2 +(at + bz . cy + dx)2 . 
Problem 6.11. Let A be an n �� n matrix with integer entries. Prove that A.1 exists and has integer entries if and only if det(A)= ��1. 
Problem 6.12. Let A be an n �� n real or complex matrix. 
(1) 
Prove that if AT = .A (A is skew-symmetric) and if n is odd, then det(A) = 0. 

(2) 
Prove that 


0 a bc .a 0 de =(af . be + dc)2 . 
.b .d 0 f 
.c .e .f 0 

6.11. PROBLEMS 
Problem 6.13. A Cauchy matrix is a matrix of the form 
.

. 

........ 

11 1 
������ 
��1 . ��1 ��1 . ��2 ��1 . ��n 
11 1 
������ 
��2 . ��1 ��2 . ��2 ��2 . ��n 
. ... 
. ... 
. ... 
11 1 
������ ��n . ��1 ��n . ��2 ��n . ��n 
........ 

where ��i = ��j, for all i, j, with 1 �� i, j �� n. Prove that the determinant Cn of a Cauchy matrix as above is given by 
i=2
 
 
n  i.1 
Cn =
j=1(��i . ��j)(��j . ��i)nn (��i . ��j)i=1j=1
 
. 

Problem 6.14. Let (��1,...,��m+1) be a sequence of pairwise distinct scalars in R and let (��1,...,��m+1) be any sequence of scalars in R, not necessarily distinct. 
(1) Prove that there is a unique polynomial P of degree at most m such that 
P (��i)= ��i, 1 �� i �� m +1. Hint. Remember Vandermonde! 
(2) Let Li(X) be the polynomial of degree m given by 
(X . ��1) ������ (X . ��i.1)(X . ��i+1) ������ (X . ��m+1)
Li(X)= ,
(��i . ��1) ������ (��i . ��i.1)(��i . ��i+1) ������ (��i . ��m+1)1 �� i �� m +1. The polynomials Li(X) are known as Lagrange polynomial interpolants. Prove that Li(��j)= ��ij 1 �� i, j �� m +1. Prove that P (X)= ��1L1(X)+ ������ + ��m+1Lm+1(X) is the unique polynomial of degree at most m such that P (��i)= ��i, 1 �� i �� m +1. 
(3) Prove that L1(X),...,Lm+1(X) are linearly independent, and that they form a basis 
of all polynomials of degree at most m. How is 1 (the constant polynomial 1) expressed over the basis (L1(X),...,Lm+1(X))? Give the expression of every polynomial P (X) of degree at most m over the basis 
(L1(X),...,Lm+1(X)). 
(4) Prove that the dual basis (L.,...,L. ) of the basis (L1(X),...,Lm+1(X)) consists 
1m+1of the linear forms L. i given by 
L . 
i (P )= P (��i), for every polynomial P of degree at most m; this is simply evaluation at ��i. 


