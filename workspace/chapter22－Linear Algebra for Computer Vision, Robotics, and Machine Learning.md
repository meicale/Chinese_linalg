Chapter 22 
Annihilating Polynomials and the Primary Decomposition 
In this chapter all vector spaces are de.ned over an arbitrary .eld K. 
In Section 6.7 we explained that if f : E ¡ú E is a linear map on a K-vector space E, then for any polynomial p(X)= a0Xd + a1Xd.1 + ¡¤¡¤¡¤ + ad with coe.cients in the .eld K, we can de.ne the linear map p(f): E ¡ú E by 
p(f)= a0fd + a1fd.1 + ¡¤¡¤¡¤ + adid, 
where fk = f .¡¤ ¡¤¡¤. f, the k-fold composition of f with itself. Note that 
p(f)(u)= a0fd(u)+ a1fd.1(u)+ ¡¤¡¤¡¤ + adu, 
for every vector u ¡Ê E. Then we showed that if E is .nite-dimensional and if ¦Öf (X)= det(XI .f) is the characteristic polynomial of f, by the Cayley¨CHamilton theorem, we have 
¦Öf (f)=0. 
This fact suggests looking at the set of all polynomials p(X) such that 
p(f)=0. 
Such polynomials are called annihilating polynomials of f, the set of all these polynomials, denoted Ann(f), is called the annihilator of f, and the Cayley-Hamilton theorem shows that it is nontrivial since it contains a polynomial of positive degree. It turns out that Ann(f) contains a polynomial mf of smallest degree that generates Ann(f), and this polynomial divides the characteristic polynomial. Furthermore, the polynomial mf encapsulates a lot of information about f, in particular whether f can be diagonalized. One of the main reasons for this is that a scalar ¦Ë ¡Ê K is a zero of the minimal polynomial mf if and only if ¦Ë is an eigenvalue of f. 
693 

The .rst main result is Theorem 22.12 which states that if f : E ¡ú E is a linear map on a .nite-dimensional space E, then f is diagonalizable i. its minimal polynomial m is of the form 
m =(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ëk), 
where ¦Ë1,...,¦Ëk are distinct elements of K. 
One of the technical tools used to prove this result is the notion of f-conductor; see De.nition 22.7. As a corollary of Theorem 22.12 we obtain results about .nite commuting families of diagonalizable or triangulable linear maps. 
If f : E ¡ú E is a linear map and ¦Ë ¡Ê K is an eigenvalue of f, recall that the eigenspace E¦Ë associated with ¦Ë is the kernel of the linear map ¦Ëid . f. If all the eigenvalues ¦Ë1 ...,¦Ëk of f are in K and if f is diagonalizable, then 
E = E¦Ë1 ¨’¡¤ ¡¤¡¤¨’ E¦Ëk , 
but in general there are not enough eigenvectors to span E. A remedy is to generalize the notion of eigenvector and look for (nonzero) vectors u (called generalized eigenvectors) such that 
(¦Ëid . f)r(u)=0, for some r ¡Ý 1. 
Then it turns out that if the minimal polynomial of f is of the form 
m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk , 
then r = ri does the job for ¦Ëi; that is, if we let 
Wi = Ker (¦Ëiid . f)ri , 
then 
E = W1 ¨’¡¤ ¡¤¡¤¨’ Wk. 
The above facts are parts of the primary decomposition theorem (Theorem 22.17). It is a special case of a more general result involving the factorization of the minimal polynomial m into its irreducible monic factors; see Theorem 22.16. 
Theorem 22.17 implies that every linear map f that has all its eigenvalues in K can be written as f = D + N, where D is diagonalizable and N is nilpotent (which means that Nr = 0 for some positive integer r). Furthermore D and N commute and are unique. This is the Jordan decomposition, Theorem 22.18. 
The Jordan decomposition suggests taking a closer look at nilpotent maps. We prove that for any nilpotent linear map f : E ¡ú E on a .nite-dimensional vector space E of dimension n over a .eld K, there is a basis of E such that the matrix N of f is of the form 
.
. 
N = 

...... 

0 ¦Í1 0 ¡¤¡¤¡¤ 00 00 ¦Í2 ¡¤¡¤¡¤ 00 ... ... 
... ... 

.  .  .  .  .  .  
0  0  0  ¡¤ ¡¤ ¡¤  0  ¦Ín  
0  0  0  ¡¤ ¡¤ ¡¤  0  0  

...... 

, 

22.1. BASIC PROPERTIES OF POLYNOMIALS; IDEALS, GCD¡¯S 
where ¦Íi = 1 or ¦Íi = 0; see Theorem 22.22. As a corollary we obtain the Jordan form; which involves matrices of the form 
.
. 
¦Ë 10 ¡¤¡¤¡¤ 0 0 ¦Ë 1 ¡¤¡¤¡¤ 0 .. .
..
.....Jr(¦Ë)= 
...... 

...... 

,

. 

.

.. 

. 

..
000 .1 00 0 ¡¤¡¤¡¤ ¦Ë 
called Jordan blocks; see Theorem 22.23. 

22.1 Basic Properties of Polynomials; Ideals, GCD¡¯s 
In order to understand the structure of Ann(f), we need to review three basic properties of polynomials. We refer the reader to Ho.man and Kunze, [38], Artin [3], Dummit and Foote [19], and Godement [26] for comprehensive discussions of polynomials and their properties. 
We begin by recalling some basic nomenclature. Given a .eld K, any nonzero polynomial p(X) ¡Ê K[X] has some monomial of highest degree a0Xn with a0 = 0, and the integer n = deg(p) ¡Ý 0 is called the degree of p. It is convenient to set the degree of the zero polynomial (denoted by 0) to be 
deg(0) = .¡Þ. 
A polynomial p(X) such that the coe.cient a0 of its monomial of highest degree is 1 is called a monic polynomial. For example, let K = R. The polynomial p(X)=4X7+2X5 is of degree 7 but is not monic since a0 = 4. On the other hand, the polynomial p(X)= X3 . 3X +1 is a monic polynomial of degree 3. 
We now discuss three key concepts of polynomial algebra: 
1. 
Ideals 

2. 
Greatest common divisors and the Bezout identity. 

3. 
Irreducible polynomials and prime factorization. 


Recall the de.nition a of ring (see De.nition 2.2). 
De.nition 22.1. A ring is a set A equipped with two operations +: A ¡Á A ¡ú A (called addition) and .: A ¡Á A ¡ú A (called multiplication) having the following properties: 
(R1) A is an abelian group w.r.t. +; 
(R2) . is associative and has an identity element 1 ¡Ê A; 
(R3) . is distributive w.r.t. +. 
The identity element for addition is denoted 0, and the additive inverse of a ¡Ê A is denoted by .a. More explicitly, the axioms of a ring are the following equations which hold for all a, b, c ¡Ê A: 
a +(b + c)=(a + b)+ c (associativity of +) (22.1) a + b = b + a (commutativity of +) (22.2) a +0=0+ a = a (zero) (22.3) a +(.a)=(.a)+ a = 0 (additive inverse) (22.4) a . (b . c)=(a . b) . c (associativity of .) (22.5) a . 1=1 . a = a (identity for .) (22.6) (a + b) . c =(a . c)+(b . c) (distributivity) (22.7) a . (b + c)=(a . b)+(a . c) (distributivity) (22.8) 
The ring A is commutative if a . b = b . a for all a, b ¡Ê A. From (22.7) and (22.8), we easily obtain a . 0=0 . a = 0 (22.9) a . (.b)=(.a) . b = .(a . b). (22.10) The .rst crucial notion is that of an ideal. De.nition 22.2. Given a commutative ring A with unit 1, an ideal of A is a nonempty 
subset I of A satisfying the following properties: (ID1) If a, b ¡Ê I, then b . a ¡Ê I. (ID2) If a ¡Ê I, then ax ¡Ê I for all x ¡Ê A. 
An ideal I is a principal ideal if there is some a ¡Ê I, called a generator, such that I = {ax | x ¡Ê A}. In this case we usually write I = aA or I =(a). The ideal I = (0) = {0} is called the null ideal (or zero ideal). The following proposition is a fundamental result about polynomials over a .eld. 
Proposition 22.1. If K is a .eld, then every polynomial ideal I . K[X] is a principal ideal. As a consequence, if I is not the zero ideal, then there is a unique monic polynomial 
p(X)= Xn + a1Xn.1 + ¡¤¡¤¡¤ + an.1X + an 
in I such that I =(p). 


22.1. BASIC PROPERTIES OF POLYNOMIALS; IDEALS, GCD¡¯S 
Proof. This result is not hard to prove if we recall that polynomials can divided. Given any two nonzero polynomials f, g ¡Ê K[X], there are unique polynomials q, r such that 
f = qg + r, and deg(r) < deg(g). (*) 
If I is not the zero ideal, there is some polynomial of smallest degree in I, and since K is a .eld, by suitable multiplication by a scalar, we can make sure that this polynomial is monic. Thus, let f be a monic polynomial of smallest degree in I. By (ID2), it is clear that (f) . I. Now let g ¡Ê I. Using (.), there exist unique q, r ¡Ê K[X] such that 
g = qf + r and deg(r) < deg(f). 
If r = 0, there is some ¦Ë = 0 in K such that ¦Ër is a monic polynomial, and since ¦Ër = ¦Ëg . ¦Ëqf, with f, g ¡Ê I, by (ID1) and (ID2), we have ¦Ër ¡Ê I, where deg(¦Ër) < deg(f) and ¦Ër is a monic polynomial, contradicting the minimality of the degree of f. Thus, r = 0, and g ¡Ê (f). The uniqueness of the monic polynomial f is left as an exercise. 
We will also need to know that the greatest common divisor of polynomials exist. Given any two nonzero polynomials f, g ¡Ê K[X], recall that f divides g if g = qf for some q ¡Ê K[X]. 
De.nition 22.3. Given any two nonzero polynomials f, g ¡Ê K[X], a polynomial d ¡Ê K[X] is a greatest common divisor of f and g (for short, a gcd of f and g) if d divides f and g and whenever h ¡Ê K[X] divides f and g, then h divides d. We say that f and g are relatively prime if 1 is a gcd of f and g. 
Note that f and g are relatively prime i. all of their gcd¡¯s are constants (scalars in K), or equivalently, if f, g have no common divisor q of degree deg(q) ¡Ý 1. For example, over R, gcd(X2 . 1,X3 + X2 . X . 1) = (X . 1)(X + 1) since X3 + X2 . X . 1=(X . 1)(X + 1)2 , while gcd(X3 +1,X . 1) = 1. 
We can characterize gcd¡¯s of polynomials as follows. 
Proposition 22.2. Let K be a .eld and let f, g ¡Ê K[X] be any two nonzero polynomials. For every polynomial d ¡Ê K[X], the following properties are equivalent: 
(1) 
The polynomial d is a gcd of f and g. 

(2) 
The polynomial d divides f and g and there exist u, v ¡Ê K[X] such that 

d = uf + vg. 

(3) 
The ideals (f), (g), and (d) satisfy the equation 


(d)=(f)+(g). 
In addition, d =0, and d is unique up to multiplication by a nonzero scalar in K. 

As a consequence of Proposition 22.2, two nonzero polynomials f, g ¡Ê K[X] are relatively prime i. there exist u, v ¡Ê K[X] such that 
uf + vg =1. 
The identity 
d = uf + vg 

of Part (2) of Lemma 22.2 is often called the Bezout identity. For an example of Bezout¡¯s identity, take K = R. Since X3 + 1 and X . 1 are relatively prime, we have 1=1/2(X3 + 1) . 1/2(X2 + X + 1)(X . 1). 
An important consequence of the Bezout identity is the following result. 
Proposition 22.3. (Euclid¡¯s proposition) Let K be a .eld and let f, g, h ¡Ê K[X] be any nonzero polynomials. If f divides gh and f is relatively prime to g, then f divides h. 
Proposition 22.3 can be generalized to any number of polynomials. 
Proposition 22.4. Let K be a .eld and let f, g1,...,gm ¡Ê K[X] be some nonzero polyno-mials. If f and gi are relatively prime for all i, 1 ¡Ü i ¡Ü m, then f and g1 ¡¤¡¤¡¤ gm are relatively prime. 
De.nition 22.3 is generalized to any .nite number of polynomials as follows. 
De.nition 22.4. Given any nonzero polynomials f1,...,fn ¡Ê K[X], where n ¡Ý 2, a poly-nomial d ¡Ê K[X] is a greatest common divisor of f1,...,fn (for short, a gcd of f1,...,fn) if d divides each fi and whenever h ¡Ê K[X] divides each fi, then h divides d. We say that f1,...,fn are relatively prime if 1 is a gcd of f1,...,fn. 
It is easily shown that Proposition 22.2 can be generalized to any .nite number of poly-nomials. 
Proposition 22.5. Let K be a .eld and let f1,...,fn ¡Ê K[X] be any n ¡Ý 2 nonzero polynomials. For every polynomial d ¡Ê K[X], the following properties are equivalent: 
(1) 
The polynomial d is a gcd of f1,...,fn. 

(2) 
The polynomial d divides each fi and there exist u1,...,un ¡Ê K[X] such that 

d = u1f1 + ¡¤¡¤¡¤ + unfn. 

(3) 
The ideals (fi), and (d) satisfy the equation 


(d)=(f1)+ ¡¤¡¤¡¤ +(fn). 
In addition, d =0, and d is unique up to multiplication by a nonzero scalar in K. 


22.1. BASIC PROPERTIES OF POLYNOMIALS; IDEALS, GCD¡¯S 
As a consequence of Proposition 22.5, any n ¡Ý 2 nonzero polynomials f1,...,fn ¡Ê K[X] are relatively prime i. there exist u1,...,un ¡Ê K[X] such that 
u1f1 + ¡¤¡¤¡¤+ unfn =1, 
the Bezout identity. 
We will also need to know that every nonzero polynomial (over a .eld) can be factored into irreducible polynomials, which is the generalization of the prime numbers to polynomials. 
De.nition 22.5. Given a .eld K, a polynomial p ¡Ê K[X] is irreducible or indecomposable or prime if deg(p) ¡Ý 1 and if p is not divisible by any polynomial q ¡Ê K[X] such that 1 ¡Ü deg(q) < deg(p). Equivalently, p is irreducible if deg(p) ¡Ý 1 and if p = q1q2, then either q1 ¡Ê K or q2 ¡Ê K (and of course, q1 = 0, q2 = 0). 
Every polynomial aX + b of degree 1 is irreducible. Over the .eld R, the polynomial X2 + 1 is irreducible (why?), but X3 + 1 is not irreducible, since 
X3 +1=(X + 1)(X2 . X + 1). 
The polynomial X2 . X + 1 is irreducible over R (why?). It would seem that X4 +1 is irreducible over R, but in fact, 
¡Ì¡Ì 
X4 +1=(X2 . 2X + 1)(X2 +2X + 1). 
However, in view of the above factorization, X4 + 1 is irreducible over Q. 
It can be shown that the irreducible polynomials over R are the polynomials of degree 1 or the polynomials of degree 2 of the form aX2 + bX + c, for which b2 .4ac < 0 (i.e., those having no real roots). This is not easy to prove! Over the complex numbers C, the only irreducible polynomials are those of degree 1. This is a version of a fact often referred to as the ¡°Fundamental Theorem of Algebra.¡± 
Observe that the de.nition of irreducibilty implies that any .nite number of distinct irreducible polynomials are relatively prime. 
The following fundamental result can be shown 
Theorem 22.6. Given any .eld K, for every nonzero polynomial 
f = adXd + ad.1Xd.1 + ¡¤¡¤¡¤+ a0 
of degree d = deg(f) ¡Ý 1 in K[X], there exists a unique set {(p1,k1£©, ..., (pm,km£©} such that 
k1 km
f = adp1 ¡¤¡¤¡¤pm , 
where the pi ¡Ê K[X] are distinct irreducible monic polynomials, the ki are (not necessarily distinct) integers, and with m ¡Ý 1, ki ¡Ý 1. 
We can now return to minimal polynomials. 


22.2 	Annihilating Polynomials and the Minimal Poly-nomial 
Given a linear map f : E ¡ú E, it is easy to check that the set Ann(f) of polynomials that annihilate f is an ideal. Furthermore, when E is .nite-dimensional, the Cayley¨CHamilton theorem implies that Ann(f) is not the zero ideal. Therefore, by Proposition 22.1, there is a unique monic polynomial mf that generates Ann(f). 
De.nition 22.6. If f : E ¡ú E is a linear map on a .nite-dimensional vector space E, the unique monic polynomial mf (X) that generates the ideal Ann(f) of polynomials which annihilate f (the annihilator of f) is called the minimal polynomial of f. 
The minimal polynomial mf of f is the monic polynomial of smallest degree that an-nihilates f. Thus, the minimal polynomial divides the characteristic polynomial ¦Öf , and deg(mf ) ¡Ý 1. For simplicity of notation, we often write m instead of mf . 
If A is any n ¡Á n matrix, the set Ann(A) of polynomials that annihilate A is the set of polynomials p(X)= a0Xd + a1Xd.1 + ¡¤¡¤¡¤ + ad.1X + ad 
such that a0Ad + a1Ad.1 + ¡¤¡¤¡¤ + ad.1A + adI =0. 
It is clear that Ann(A) is a nonzero ideal and its unique monic generator is called the minimal polynomial of A. We check immediately that if Q is an invertible matrix, then A and Q.1AQ have the same minimal polynomial. Also, if A is the matrix of f with respect to some basis, then f and A have the same minimal polynomial. 
The zeros (in K) of the minimal polynomial of f and the eigenvalues of f (in K) are intimately related. 
Proposition 22.7. Let f : E ¡ú E be a linear map on some .nite-dimensional vector space 
E. Then ¦Ë ¡Ê K is a zero of the minimal polynomial mf (X) of f i. ¦Ë is an eigenvalue of f i. ¦Ë is a zero of ¦Öf (X). Therefore, the minimal and the characteristic polynomials have the same zeros (in K), except for multiplicities. 
Proof. First assume that m(¦Ë) = 0 (with ¦Ë ¡Ê K, and writing m instead of mf ). If so, using polynomial division, m can be factored as 
m =(X . ¦Ë)q, 
with deg(q) < deg(m). Since m is the minimal polynomial, q(f) = 0, so there is some nonzero vector v ¡Ê E such that u = q(f)(v) = 0. But then, because m is the minimal polynomial, 
0= m(f)(v) 
=(f . ¦Ëid)(q(f)(v)) 
=(f . ¦Ëid)(u), 

22.3. MINIMAL POLYNOMIALS OF DIAGONALIZABLE LINEAR MAPS 
which shows that ¦Ë is an eigenvalue of f. 
Conversely, assume that ¦Ë ¡Ê K is an eigenvalue of f. This means that for some u = 0, we have f(u)= ¦Ëu. Now it is easy to show that 
m(f)(u)= m(¦Ë)u, 
and since m is the minimal polynomial of f, we have m(f)(u)=0, so m(¦Ë)u = 0, and since u = 0, we must have m(¦Ë) = 0. 
Proposition 22.8. Let f : E ¡ú E be a linear map on some .nite-dimensional vector space 
E. If f diagonalizable, then its minimal polynomial is a product of distinct factors of degree 
1. 
Proof. If we assume that f is diagonalizable, then its eigenvalues are all in K, and if ¦Ë1,...,¦Ëk are the distinct eigenvalues of f, and then by Proposition 22.7, the minimal polynomial m of f must be a product of powers of the polynomials (X . ¦Ëi). Actually, we claim that 
m =(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ëk). 
For this we just have to show that m annihilates f. However, for any eigenvector u of f, one of the linear maps f . ¦Ëiid sends u to 0, so 
m(f)(u)=(f . ¦Ë1id) .¡¤ ¡¤¡¤. (f . ¦Ëkid)(u)=0. 
Since E is spanned by the eigenvectors of f, we conclude that 
m(f)=0. 
It turns out that the converse of Proposition 22.8 is true, but this will take a little work to establish it. 

22.3 	Minimal Polynomials of Diagonalizable Linear Maps 
In this section we prove that if the minimal polynomial mf of a linear map f is of the form 
mf =(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ëk) 
for distinct scalars ¦Ë1,...,¦Ëk ¡Ê K, then f is diagonalizable. This is a powerful result that has a number of implications. But .rst we need of few properties of invariant subspaces. 
Given a linear map f : E ¡ú E, recall that a subspace W of E is invariant under f if f(u) ¡Ê W for all u ¡Ê W . For example, if f : R2 ¡ú R2 is f(x, y)=(.x, y), the y-axis is invariant under f. 
Proposition 22.9. Let W be a subspace of E invariant under the linear map f : E ¡ú E (where E is .nite-dimensional). Then the minimal polynomial of the restriction f | W of f to W divides the minimal polynomial of f, and the characteristic polynomial of f | W divides the characteristic polynomial of f. 
Sketch of proof. The key ingredient is that we can pick a basis (e1,...,en) of E in which (e1,...,ek) is a basis of W . The matrix of f over this basis is a block matrix of the form 
  
BC 
A =,
0 D
where B is a k ¡Á k matrix, D is an (n . k) ¡Á (n . k) matrix, and C is a k ¡Á (n . k) matrix. Then 
det(XI . A) = det(XI . B) det(XI . D), 
which implies the statement about the characteristic polynomials. Furthermore, 
  
Bi 
CiAi 
=,
0 Di
for some k ¡Á (n . k) matrix Ci. It follows that any polynomial which annihilates A also annihilates B and D. So the minimal polynomial of B divides the minimal polynomial of 
A. 
For the next step, there are at least two ways to proceed. We can use an old-fashion argument using Lagrange interpolants, or we can use a slight generalization of the notion of annihilator. We pick the second method because it illustrates nicely the power of principal ideals. 
What we need is the notion of conductor (also called transporter). 
De.nition 22.7. Let f : E ¡ú E be a linear map on a .nite-dimensional vector space E, let W be an invariant subspace of f, and let u be any vector in E. The set Sf (u, W ) consisting of all polynomials q ¡Ê K[X] such that q(f)(u) ¡Ê W is called the f-conductor of u into W . 
Observe that the minimal polynomial mf of f always belongs to Sf (u, W ), sothisisa nontrivial set. Also, if W = (0), then Sf (u, (0)) is just the annihilator of f. The crucial property of Sf (u, W ) is that it is an ideal. 
Proposition 22.10. If W is an invariant subspace for f, then for each u ¡Ê E, the f-conductor Sf (u, W ) is an ideal in K[X]. 
We leave the proof as a simple exercise, using the fact that if W invariant under f, then W is invariant under every polynomial q(f) in Sf (u, W ). 
Since Sf (u, W ) is an ideal, it is generated by a unique monic polynomial q of smallest degree, and because the minimal polynomial mf of f is in Sf (u, W ), the polynomial q divides 
m. 

22.3. MINIMAL POLYNOMIALS OF DIAGONALIZABLE LINEAR MAPS 
De.nition 22.8. The unique monic polynomial which generates Sf (u, W ) is called the conductor of u into W . 
Example 22.1. For example, suppose f : R2 ¡ú R2 where f(x, y)=(x, 0). Observe that 10 
W = {(x, 0) ¡Ê R2} is invariant under f. By representing f as , we see that mf (X)= 
00 ¦Öf (X)= X2 . X. Let u = (0,y). Then Sf (u, W )=(X), and we say X is the conductor of u into W . 
Proposition 22.11. Let f : E ¡ú E be a linear map on a .nite-dimensional space E and assume that the minimal polynomial m of f is of the form 
m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk , 
where the eigenvalues ¦Ë1,...,¦Ëk of f belong to K. If W is a proper subspace of E which is invariant under f, then there is a vector u ¡Ê E with the following properties: 
(a) 
u/¡Ê W ; 

(b) 
(f . ¦Ëid)(u) ¡Ê W , for some eigenvalue ¦Ë of f. 


Proof. Observe that (a) and (b) together assert that the conductor of u into W is a polyno-mial of the form X . ¦Ëi. Pick any vector v ¡Ê E not in W , and let g be the conductor of v into W , i.e. g(f)(v) ¡Ê W . Since g divides m and v/¡Ê W , the polynomial g is not a constant, and thus it is of the form 
g =(X . ¦Ë1)s1 ¡¤¡¤¡¤ (X . ¦Ëk)sk , 
with at least some si > 0. Choose some index j such that sj > 0. Then X . ¦Ëj is a factor of g, so we can write 
g =(X . ¦Ëj)q. (*) 
By de.nition of g, the vector u = q(f)(v) cannot be in W , since otherwise g would not be of minimal degree. However, (.) implies that 
(f . ¦Ëjid)(u)=(f . ¦Ëjid)(q(f)(v)) = g(f)(v) 
is in W , which concludes the proof. 
We can now prove the main result of this section. 
Theorem 22.12. Let f : E ¡ú E be a linear map on a .nite-dimensional space E. Then f is diagonalizable i. its minimal polynomial m is of the form 
m =(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ëk), 
where ¦Ë1,...,¦Ëk are distinct elements of K. 
Proof. We already showed in Proposition 22.8 that if f is diagonalizable, then its minimal polynomial is of the above form (where ¦Ë1,...,¦Ëk are the distinct eigenvalues of f). 
For the converse, let W be the subspace spanned by all the eigenvectors of f. If W = E, since W is invariant under f, by Proposition 22.11, there is some vector u ¡Ê/W such that for some ¦Ëj, we have 
(f . ¦Ëjid)(u) ¡Ê W. 
Let v =(f . ¦Ëjid)(u) ¡Ê W . Since v ¡Ê W , we can write 
v = w1 + ¡¤¡¤¡¤ + wk 
where f(wi)= ¦Ëiwi (either wi = 0 or wi is an eigenvector for ¦Ëi), and so for every polynomial h, we have 
h(f)(v)= h(¦Ë1)w1 + ¡¤¡¤¡¤ + h(¦Ëk)wk, 
which shows that h(f)(v) ¡Ê W for every polynomial h. We can write 
m =(X . ¦Ëj)q 
for some polynomial q, and also 
q . q(¦Ëj)= p(X . ¦Ëj) 
for some polynomial p. We know that p(f)(v) ¡Ê W , and since m is the minimal polynomial of f, we have 
0= m(f)(u)=(f . ¦Ëjid)(q(f)(u)), 
which implies that q(f)(u) ¡Ê W (either q(f)(u) = 0, or it is an eigenvector associated with ¦Ëj). However, 
q(f)(u) . q(¦Ëj)u = p(f)((f . ¦Ëjid)(u)) = p(f)(v), 
and since p(f)(v) ¡Ê W and q(f)(u) ¡Ê W , we conclude that q(¦Ëj)u ¡Ê W . But, u/¡Ê W , which implies that q(¦Ëj) = 0, so ¦Ëj is a double root of m, a contradiction. Therefore, we must have W = E. 
Remark: 	Proposition 22.11 can be used to give a quick proof of Theorem 14.5. 

22.4 	Commuting Families of Diagonalizable and Trian-gulable Maps 
Using Theorem 22.12, we can give a short proof about commuting diagonalizable linear maps. 
De.nition 22.9. If F is a family of linear maps on a vector space E, we say that F is a commuting family i. f . g = g . f for all f, g ¡ÊF. 

22.4. COMMUTING FAMILIES OF LINEAR MAPS 
Proposition 22.13. Let F be a commuting family of diagonalizable linear maps on a vector space E. There exists a basis of E such that every linear map in F is represented in that basis by a diagonal matrix. 
Proof. We proceed by induction on n = dim(E). If n = 1, there is nothing to prove. If n> 1, there are two cases. If all linear maps in F are of the form ¦Ëid for some ¦Ë ¡Ê K, then the proposition holds trivially. In the second case, let f ¡ÊF be some linear map in F which is not a scalar multiple of the identity. In this case, f has at least two distinct eigenvalues ¦Ë1,...,¦Ëk, and because f is diagonalizable, E is the direct sum of the corresponding eigenspaces E¦Ë1 ,...,E¦Ëk . For every index i, the eigenspace E¦Ëi is invariant under f and under every other linear map g in F, since for any g ¡ÊF and any u ¡Ê E¦Ëi , because f and g commute, we have 
f(g(u)) = g(f(u)) = g(¦Ëiu)= ¦Ëig(u) 
so g(u) ¡Ê E¦Ëi . Let Fi be the family obtained by restricting each f ¡ÊF to E¦Ëi . By Proposition 22.9, the minimal polynomial of every linear map f | E¦Ëi in Fi divides the minimal polynomial mf of f, and since f is diagonalizable, mf is a product of distinct linear factors, so the minimal polynomial of f | E¦Ëi is also a product of distinct linear factors. By Theorem 22.12, the linear map f | E¦Ëi is diagonalizable. Since k> 1, we have dim(E¦Ëi ) < dim(E) for i =1,...,k, and by the induction hypothesis, for each i there is a basis of E¦Ëi over which f | E¦Ëi is represented by a diagonal matrix. Since the above argument holds for all i, by combining the bases of the E¦Ëi , we obtain a basis of E such that the matrix of every linear map f ¡ÊF is represented by a diagonal matrix. 
There is also an analogous result for commuting families of linear maps represented by upper triangular matrices. To prove this we need the following proposition. 
Proposition 22.14. Let F be a nonempty commuting family of triangulable linear maps on a .nite-dimensional vector space E. Let W be a proper subspace of E which is invariant under F. Then there exists a vector u ¡Ê E such that: 
1. 
u/¡Ê W . 

2. 
For every f ¡ÊF, the vector f(u) belongs to the subspace W ¨’ Ku spanned by W and 


u. 
Proof. By renaming the elements of F if necessary, we may assume that (f1,...,fr) is a basis of the subspace of End(E) spanned by F. We prove by induction on r that there exists some vector u ¡Ê E such that 
1. 
u/¡Ê W . 

2. 
(fi . ¦Áiid)(u) ¡Ê W for i =1,...,r, for some scalars ¦Ái ¡Ê K. 


Consider the base case r = 1. Since f1 is triangulable, its eigenvalues all belong to K since they are the diagonal entries of the triangular matrix associated with f1 (this is the easy direction of Theorem 14.5), so the minimal polynomial of f1 is of the form 
m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk , 
where the eigenvalues ¦Ë1,...,¦Ëk of f1 belong to K. We conclude by applying Proposition 
22.11. 
Next assume that r ¡Ý 2 and that the induction hypothesis holds for f1,...,fr.1. Thus, there is a vector ur.1 ¡Ê E such that 
1. 
ur.1 ¡Ê/W . 

2. 
(fi . ¦Áiid)(ur.1) ¡Ê W for i =1,...,r . 1, for some scalars ¦Ái ¡Ê K. 


Let 
Vr.1 = {w ¡Ê E | (fi . ¦Áiid)(w) ¡Ê W, i =1,...,r . 1}. 

Clearly, W . Vr.1 and ur.1 ¡Ê Vr.1. We claim that Vr.1 is invariant under F. This is because, for any v ¡Ê Vr.1 and any f ¡ÊF, since f and fi commute, we have 
(fi . ¦Áiid)(f(v)) = f((fi . ¦Áiid)(v)), 1 ¡Ü i ¡Ü r . 1. 
Now (fi . ¦Áiid)(v) ¡Ê W because v ¡Ê Vr.1, and W is invariant under F, so f(fi . ¦Áiid)(v)) ¡Ê W , that is, (fi . ¦Áiid)(f(v)) ¡Ê W . 
Consider the restriction gr of fr to Vr.1. The minimal polynomial of gr divides the minimal polynomial of fr, and since fr is triangulable, just as we saw for f1, the minimal polynomial of fr is of the form 
m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk , 
where the eigenvalues ¦Ë1,...,¦Ëk of fr belong to K, so the minimal polynomial of gr is of the same form. By Proposition 22.11, there is some vector ur ¡Ê Vr.1 such that 
1. 
ur ¡Ê/W . 

2. 
(gr . ¦Árid)(ur) ¡Ê W for some scalars ¦Ár ¡Ê K. 


Now since ur ¡Ê Vr.1, we have (fi . ¦Áiid)(ur) ¡Ê W for i =1,...,r . 1, so (fi . ¦Áiid)(ur) ¡Ê W for i =1,...,r (since gr is the restriction of fr), which concludes the proof of the induction step. Finally, since every f ¡ÊF is the linear combination of (f1,...,fr), Condition (2) of the inductive claim implies Condition (2) of the proposition. 
We can now prove the following result. 
Proposition 22.15. Let F be a nonempty commuting family of triangulable linear maps on a .nite-dimensional vector space E. There exists a basis of E such that every linear map in F is represented in that basis by an upper triangular matrix. 

22.5. THE PRIMARY DECOMPOSITION THEOREM 
Proof. Let n = dim(E). We construct inductively a basis (u1,...,un) of E such that if Wi is the subspace spanned by (u1 ...,ui), then for every f ¡ÊF, 
f(ui)= a f ¡¤¡¤ + a f 
1iu1 + ¡¤ iiui, 
for some aijf ¡Ê K; that is, f(ui) belongs to the subspace Wi. We begin by applying Proposition 22.14 to the subspace W0 = (0) to get u1 so that for all f ¡ÊF, f(u1)= ¦Á1 f u1. 
For the induction step, since Wi invariant under F, we apply Proposition 22.14 to the subspace Wi, to get ui+1 ¡Ê E such that 
1. 
ui+1 ¡Ê/Wi. 

2. 
For every f ¡ÊF, the vector f(ui+1) belong to the subspace spanned by Wi and ui+1. 


Condition (1) implies that (u1,...,ui,ui+1) is linearly independent, and Condition (2) means that for every f ¡ÊF, 
f(ui+1)= a f ¡¤¡¤ + a f 
1i+1u1 + ¡¤ i+1i+1ui+1, for some aif +1j ¡Ê K, establishing the induction step. After n steps, each f ¡ÊF is represented by an upper triangular matrix. 
Observe that if F consists of a single linear map f and if the minimal polynomial of f is of the form m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk , 
with all ¦Ëi ¡Ê K, using Proposition 22.11 instead of Proposition 22.14, the proof of Proposition 
22.15 yields another proof of Theorem 14.5. 

22.5 The Primary Decomposition Theorem 
If f : E ¡ú E is a linear map and ¦Ë ¡Ê K is an eigenvalue of f, recall that the eigenspace E¦Ë associated with ¦Ë is the kernel of the linear map ¦Ëid . f. If all the eigenvalues ¦Ë1 ...,¦Ëk of f are in K, it may happen that 
E = E¦Ë1 ¨’¡¤ ¡¤¡¤¨’ E¦Ëk , 
but in general there are not enough eigenvectors to span E. What if we generalize the notion of eigenvector and look for (nonzero) vectors u such that 
(¦Ëid . f)r(u)=0, for some r ¡Ý 1? 
It turns out that if the minimal polynomial of f is of the form 
m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk , 
then r = ri does the job for ¦Ëi; that is, if we let 
Wi = Ker (¦Ëiid . f)ri , 
then 
E = W1 ¨’¡¤ ¡¤¡¤¨’ Wk. 
This result is very nice but seems to require that the eigenvalues of f all belong to K. Actually, it is a special case of a more general result involving the factorization of the minimal polynomial m into its irreducible monic factors (see Theorem 22.6), 
r1 rk
m = p ¡¤¡¤¡¤ p,
1 k 
where the pi are distinct irreducible monic polynomials over K. 
Theorem 22.16. (Primary Decomposition Theorem) Let f : E ¡ú E be a linear map on the .nite-dimensional vector space E over the .eld K. Write the minimal polynomial m of f as 
r1 rk
m = p ¡¤¡¤¡¤ p,
1 k 
where the pi are distinct irreducible monic polynomials over K, and the ri are positive inte-gers. Let 
Wi = Ker (p ri i (f)),i =1, . . . , k. 
Then 
(a) 
E = W1 ¨’¡¤ ¡¤¡¤¨’ Wk. 

(b) 
Each Wi is invariant under f. 

(c) 
The minimal polynomial of the restriction f | Wi of f to Wi is p ri i . 


Proof. The trick is to construct projections ¦Ði using the polynomials p rj so that the range 
j of ¦Ði is equal to Wi. Let 
 
gi = m/piri =pjrj . 
j=i
M
Note that p ri i gi = m. 
Since p1,...,pk are irreducible and distinct, they are relatively prime. Then using Proposi-tion 22.4, it is easy to show that g1,...,gk are relatively prime. Otherwise, some irreducible polynomial p would divide all of g1,...,gk, so by Proposition 22.4 it would be equal to one of the irreducible factors pi. But that pi is missing from gi, a contradiction. Therefore, by Proposition 22.5, there exist some polynomials h1,...,hk such that 
g1h1 + ¡¤¡¤¡¤ + gkhk =1. 

22.5. THE PRIMARY DECOMPOSITION THEOREM 
Let qi = gihi and let ¦Ði = qi(f)= gi(f)hi(f). We have q1 + ¡¤¡¤¡¤ + qk =1, and since m divides qiqj for i = j, we get ¦Ð1 + ¡¤¡¤¡¤ + ¦Ðk = id ¦Ði¦Ðj =0,i = j. (We implicitly used the fact that if p, q are two polynomials, the linear maps p(f) . q(f) and q(f) . p(f) are the same since p(f) and q(f) are polynomials in the powers of f, which commute.) Composing the .rst equation with ¦Ði and using the second equation, we get ¦Ð2 
i = ¦Ði. Therefore, the ¦Ði are projections, and E is the direct sum of the images of the ¦Ði. Indeed, every u ¡Ê E can be expressed as u = ¦Ð1(u)+ ¡¤¡¤¡¤ + ¦Ðk(u). Also, if ¦Ð1(u)+ ¡¤¡¤¡¤ + ¦Ðk(u)=0, then by applying ¦Ði we get 0= ¦Ði 2(u)= ¦Ði(u),i =1, . . . k. To .nish proving (a), we need to show that Wi = Ker (piri (f)) = ¦Ði(E). If v ¡Ê ¦Ði(E), then v = ¦Ði(u) for some u ¡Ê E, so 
ri ri
pi (f)(v)= pi (f)(¦Ði(u)) = piri (f)gi(f)hi(f)(u) = hi(f)piri (f)gi(f)(u) = hi(f)m(f)(u)=0, 
because m is the minimal polynomial of f. Therefore, v ¡Ê Wi. 
ri ri
Conversely, assume that v ¡Ê Wi = Ker (pi (f)). If j = i, then gjhj is divisible by pi , so gj(f)hj(f)(v)= ¦Ðj(v)=0,j = i. Then since ¦Ð1 + ¡¤¡¤¡¤ + ¦Ðk = id, we have v = ¦Ðiv, which shows that v is in the range of ¦Ði. Therefore, Wi = Im(¦Ði), and this .nishes the proof of (a). 
ri ri ri
If pi (f)(u) = 0, then pi (f)(f(u)) = f(pi (f)(u)) = 0, so (b) holds. 
ri ri
If we write fi = f | Wi, then p (fi) = 0, because p (f) = 0 on Wi (its kernel). Therefore, 
ii the minimal polynomial of fi divides piri . Conversely, let q be any polynomial such that q(fi)=0 (on Wi). Since m = p ri i gi, the fact that m(f)(u)=0 for all u ¡Ê E shows that 
piri (f)(gi(f)(u)) = 0,u ¡Ê E, 
and thus Im(gi(f)) . Ker (piri (f)) = Wi. Consequently, since q(f) is zero on Wi, 
q(f)gi(f)=0 forall u ¡Ê E. 
ri ri
But then qgi is divisible by the minimal polynomial m = pi gi of f, and since pi and gi are relatively prime, by Euclid¡¯s proposition, p ri i must divide q. This .nishes the proof that the minimal polynomial of fi is piri , which is (c). 
To best understand the projection constructions of Theorem 22.16, we provide the fol-lowing two explicit examples of the primary decomposition theorem. 
Example 22.2. First let f : R3 ¡ú R3 be de.ned as f(x, y, z)=(y, .x, z). In terms of the 
.. 
0 .10 
..
standard basis f is represented by the 3 ¡Á 3 matrix Xf := 1 0 0 . Thenasimple 
001 calculation shows that mf (x)= ¦Öf (x)=(x2 + 1)(x . 1). Using the notation of the preceding proof set 
m = p1p2,p1 = x 2 +1,p2 = x . 1. 
Then 
mm 
g1 == x . 1,g2 == x 2 +1. p1 p2 
We must .nd h1,h2 ¡Ê R[x] such that g1h1 + g2h2 = 1. In general this is the hard part of the projection construction. But since we are only working with two relatively prime polynomials g1,g2, we may apply the Euclidean algorithm to discover that 
x +1 1 
. (x . 1)+ (x 2 +1) = 1,
22
= .x+1 1

where h12 while h2 = 2 . By de.nition 
.. 
100 ¦Ð1 = g1(f)h1(f)= . 1(Xf . id)(Xf + id) = . 1(Xf 2 . id) = .010. ,
22
000 
and 
.. 
000 ¦Ð2 = g2(f)h2(f)= 1(Xf 2 + id) = .000. . 
2001 

22.5. THE PRIMARY DECOMPOSITION THEOREM 
Then R3 = W1 ¨’ W2, where 
W1 = ¦Ð1(R3) = Ker(p1(Xf )) = Ker(Xf 2 + id) 
.	. 
000 
.	.
=Ker 	000 = {(x, y, 0) ¡Ê R3}, 001 
W2 = ¦Ð2(R3) = Ker(p2(Xf )) = Ker(Xf . id) 
.. 
.1 .1	0 
..
= Ker 	1 .10 = {(0, 0,z) ¡Ê R3}. 0 00 
Example 22.3. For our second example of the primary decomposition theorem let f : R3 ¡ú 
R3 be de.ned as f(x, y, z)=(y, .x + z, .y), with standard matrix representation Xf = 

.. 
0 .1	0 
.	10 .1.. A simple calculation shows that mf (x)= ¦Öf (x)= x(x2 + 2). Set 01 0 
mf 	mf 
p1 = x 2 +2,p2 = x, g1 == x, g2 == x 2 +2. p1 p2 
Since gcd(g1,g2) = 1, we use the Euclidean algorithm to .nd 
11 
h1 = . x, h2 = ,
22
such that g1h1 + g2h2 = 1. Then 
.. 
1 .1
0
2	2
1 
¦Ð1 = g1(f)h1(f)= . Xf 2 = . 010 . ,
2 	1
.1 
0
22 
while 
.	. 
11
0
22 ¦Ð2 = g2(f)h2(f)= 1(Xf 2 + 2id) = .000. . 211
0
22 
Although it is not entirely obvious, ¦Ð1 and ¦Ð2 are indeed projections since 
. .. .. . 
1 .11 .11 .1
00 	0
2222 22 
¦Ð12 = . 010 .. 010 . = . 010 . = ¦Ð1, 
.111 	1
.1 	.1
00 	0
2222 22 
and 
. .. .. . 
1111 11
00 0
2222 22 
¦Ð22 = .000..000. = .000. = ¦Ð2. 
1111 11
00 0
2222 22 Furthermore observe that ¦Ð1 + ¦Ð2 = id. The primary decomposition theorem implies that R3 = W1 ¨’ W2 where 
W1 = ¦Ð1(R3) = Ker(p1(f)) = Ker(X2 + 2) 
.	. 
101 
.
=Ker 	000. = span{(0, 1, 0), (1, 0, .1)}, 101 
W2 = ¦Ð2(R3) = Ker(p2(f)) = Ker(X) = span{(1, 0, 1)}. See Figure 22.1. 

Figure 22.1: The direct sum decomposition of R3 = W1 ¨’W2 where W1 is the plane x+z =0 and W2 is line t(1, 0, 1). The spanning vectors of W1 are in blue. 
If all the eigenvalues of f belong to the .eld K, we obtain the following result. 
Theorem 22.17. (Primary Decomposition Theorem, Version 2) Let f : E ¡ú E be a lin-ear map on the .nite-dimensional vector space E over the .eld K. If all the eigenvalues ¦Ë1,...,¦Ëk of f belong to K, write 
m =(X . ¦Ë1)r1 ¡¤¡¤¡¤ (X . ¦Ëk)rk 

22.5. THE PRIMARY DECOMPOSITION THEOREM 
for the minimal polynomial of f, 
¦Öf =(X . ¦Ë1)n1 ¡¤¡¤¡¤ (X . ¦Ëk)nk 
for the characteristic polynomial of f, with 1 ¡Ü ri ¡Ü ni, and let 
Wi = Ker (¦Ëiid . f)ri ,i =1, . . . , k. 
Then 
(a) 
E = W1 ¨’¡¤ ¡¤¡¤¨’ Wk. 

(b) 
Each Wi is invariant under f. 

(c) 
dim(Wi)= ni. 

(d) 
The minimal polynomial of the restriction f | Wi of f to Wi is (X . ¦Ëi)ri . 


Proof. Parts (a), (b) and (d) have already been proven in Theorem 22.16, so it remains to prove (c). Since Wi is invariant under f, let fi be the restriction of f to Wi. The characteristic polynomial ¦Öfi of fi divides ¦Ö(f), and since ¦Ö(f) has all its roots in K, so does ¦Öi(f). By Theorem 14.5, there is a basis of Wi in which fi is represented by an upper triangular matrix, and since (¦Ëiid . f)ri = 0, the diagonal entries of this matrix are equal to ¦Ëi. Consequently, 
=(X . ¦Ëi)dim(Wi)
¦Öfi , 
and since ¦Öfi divides ¦Ö(f), we conclude hat 
dim(Wi) ¡Ü ni,i =1, . . . , k. 
Because E is the direct sum of the Wi, we have dim(W1)+ ¡¤¡¤¡¤ + dim(Wk)= n, and since n1 + ¡¤¡¤¡¤ + nk = n, we must have 
dim(Wi)= ni,i =1, . . . , k, proving (c). De.nition 22.10. If ¦Ë ¡Ê K is an eigenvalue of f, we de.ne a generalized eigenvector of f as a nonzero vector u ¡Ê E such that (¦Ëid . f)r(u)=0, for some r ¡Ý 1. The index of ¦Ë is de.ned as the smallest r ¡Ý 1 such that 
Ker (¦Ëid . f)r = Ker (¦Ëid . f)r+1 . 
It is clear that Ker (¦Ëid . f)i . Ker (¦Ëid . f)i+1 for all i ¡Ý 1. By Theorem 22.17(d), if ¦Ë = ¦Ëi, the index of ¦Ëi is equal to ri. 

22.6 Jordan Decomposition 
Recall that a linear map g : E ¡ú E is said to be nilpotent if there is some positive integer r such that gr = 0. Another important consequence of Theorem 22.17 is that f can be written as the sum of a diagonalizable and a nilpotent linear map (which commute). For example 
f : R2 ¡ú R2 be the R-linear map f(x, y)=(x, x + y) with standard matrix representation Xf = 1 1 . A basic calculation shows that mf (x)= ¦Öf (x)=(x . 1)2 . By Theorem 
01 
22.12 we know that f is not diagonalizable over R. But since the eigenvalue ¦Ë1 = 1 of f does belong to R, we may use the projection construction inherent within Theorem 22.17 to write f = D + N, where D is a diagonalizable linear map and N is a nilpotent linear map. The proof of Theorem 22.16 implies that 
r1
p1 =(x . 1)2 ,g1 =1= h1,¦Ð1 = g1(f)h1(f) = id. Then D = ¦Ë1¦Ð1 = id, N = f . D = f(x, y) . id(x, y)=(x, x + y) . (x, y) = (0,y), which is equivalent to the matrix decomposition 11 1001 
Xf ==+ . 
01 0100 This example suggests that the diagonal summand of f is related to the projection 
constructions associated with the proof of the primary decomposition theorem. If we write 
D = ¦Ë1¦Ð1 + ¡¤¡¤¡¤ + ¦Ëk¦Ðk, 
where ¦Ði is the projection from E onto the subspace Wi de.ned in the proof of Theorem 

22.16, since ¦Ð1 + ¡¤¡¤¡¤ + ¦Ðk = id, we have f = f¦Ð1 + ¡¤¡¤¡¤ + f¦Ðk, and so we get N = f . D =(f . ¦Ë1id)¦Ð1 + ¡¤¡¤¡¤ +(f . ¦Ëkid)¦Ðk. We claim that N = f . D is a nilpotent operator. Since by construction the ¦Ði are polyno-mials in f, they commute with f, using the properties of the ¦Ði, we get Nr =(f . ¦Ë1id)r¦Ð1 + ¡¤¡¤¡¤ +(f . ¦Ëkid)r¦Ðk. Therefore, if r = max{ri}, we have (f . ¦Ëkid)r = 0 for i =1,...,k, which implies that Nr 
=0. 

22.6. JORDAN DECOMPOSITION 
It remains to show that D is diagonalizable. Since N is a polynomial in f, it commutes with f, and thus with D. From 
D = ¦Ë1¦Ð1 + ¡¤¡¤¡¤ + ¦Ëk¦Ðk, 
and ¦Ð1 + ¡¤¡¤¡¤ + ¦Ðk = id, 
we see that 
D . ¦Ëiid = ¦Ë1¦Ð1 + ¡¤¡¤¡¤ + ¦Ëk¦Ðk . ¦Ëi(¦Ð1 + ¡¤¡¤¡¤ + ¦Ðk) =(¦Ë1 . ¦Ëi)¦Ð1 + ¡¤¡¤¡¤ +(¦Ëi.1 . ¦Ëi)¦Ði.1 +(¦Ëi+1 . ¦Ëi)¦Ði+1 
+ ¡¤¡¤¡¤ +(¦Ëk . ¦Ëi)¦Ðk. 
Since the projections ¦Ðj with j = i vanish on Wi, the above equation implies that D . ¦Ëiid vanishes on Wi and that (D . ¦Ëjid)(Wi) . Wi, and thus that the minimal polynomial of D is 
(X . ¦Ë1) ¡¤¡¤¡¤ (X . ¦Ëk). 
Since the ¦Ëi are distinct, by Theorem 22.12, the linear map D is diagonalizable. 
In summary we have shown that when all the eigenvalues of f belong to K, there exist a diagonalizable linear map D and a nilpotent linear map N such that 
f = D + N DN = ND, 
and N and D are polynomials in f. 
De.nition 22.11. A decomposition of f as f = D + N as above is called a Jordan decom-position. 
In fact, we can prove more: the maps D and N are uniquely determined by f. 
Theorem 22.18. (Jordan Decomposition) Let f : E ¡ú E be a linear map on the .nite-dimensional vector space E over the .eld K. If all the eigenvalues ¦Ë1,...,¦Ëk of f belong to K, then there exist a diagonalizable linear map D and a nilpotent linear map N such that 
f = D + N DN = ND. 
Furthermore, D and N are uniquely determined by the above equations and they are polyno-mials in f. 
Proof. We already proved the existence part. Suppose we also have f = D' + N', with D'N' = N'D', where D' is diagonalizable, N' is nilpotent, and both are polynomials in f. We need to prove that D = D' and N = N'. 
Since D ' and N ' commute with one another and f = D ' + N ' , we see that D ' and N ' 
commute with f. Then D ' and N ' commute with any polynomial in f; hence they commute with D and N. From 
D + N = D ' + N ' , 
we get 
D . D ' = N ' . N, 
and D, D ' , N, N ' commute with one another. Since D and D ' are both diagonalizable and commute, by Proposition 22.13, they are simultaneousy diagonalizable, so D . D ' is diago-nalizable. Since N and N ' commute, by the binomial formula, for any r ¡Ý 1, 
r
(N ' . N)r = Ù² (.1)j 	r (N ' )r.jNj. j
j=0 
'	' )r2
Since both N and N are nilpotent, we have Nr1 =0 and (N = 0, for some r1,r2 > 0, so for r ¡Ý r1 + r2, the right-hand side of the above expression is zero, which shows that N ' . N is nilpotent. (In fact, it is easy that r1 = r2 = n works). It follows that D . D ' = N ' . N is both diagonalizable and nilpotent. Clearly, the minimal polynomial of a nilpotent linear map is of the form Xr for some r> 0 (and r ¡Ü dim(E)). But D . D ' is diagonalizable, so its minimal polynomial has simple roots, which means that r = 1. Therefore, the minimal polynomial of D . D ' is X, which says that D . D ' = 0, and then N = N ' . 
If K is an algebraically closed .eld, then Theorem 22.18 holds. This is the case when K = C. This theorem reduces the study of linear maps (from E to itself) to the study of nilpotent operators. There is a special normal form for such operators which is discussed in the next section. 

22.7 Nilpotent Linear Maps and Jordan Form 
This section is devoted to a normal form for nilpotent maps. We follow Godement¡¯s expo-sition [27]. Let f : E ¡ú E be a nilpotent linear map on a .nite-dimensional vector space over a .eld K, and assume that f is not the zero map. There is a smallest positive integer r ¡Ý 1 such fr = 0 and fr+1 = 0. Clearly, the polynomial Xr+1 annihilates f, and it is the minimal polynomial of f since fr = 0. It follows that r +1 ¡Ü n = dim(E). Let us de.ne the subspaces Ni by 
Ni = Ker (fi),i ¡Ý 0. 
Note that N0 = (0), N1 = Ker (f), and Nr+1 = E. Also, it is obvious that 
Ni . Ni+1,i ¡Ý 0. 
Proposition 22.19. Given a nilpotent linear map f with fr =0 and fr+1 =0 as above, the inclusions in the following sequence are strict: 
(0) = N0 . N1 . ¡¤ ¡¤¡¤ . Nr . Nr+1 = E. 

22.7. NILPOTENT LINEAR MAPS AND JORDAN FORM 
Proof. We proceed by contradiction. Assume that Ni = Ni+1 for some i with 0 ¡Ü i ¡Ü r. Since fr+1 = 0, for every u ¡Ê E, we have 
0= fr+1(u)= fi+1(fr.i(u)), 
which shows that fr.i(u) ¡Ê Ni+1. Since Ni = Ni+1, we get fr.i(u) ¡Ê Ni, and thus fr(u) = 0. Since this holds for all u ¡Ê E, we see that fr = 0, a contradiction. 
Proposition 22.20. Given a nilpotent linear map f with fr =0 and fr+1 =0, for any integer i with 1 ¡Ü i ¡Ü r, for any subspace U of E, if U ¡É Ni = (0), then f(U) ¡É Ni.1 = (0), and the restriction of f to U is an isomorphism onto f(U). 
Proof. Pick v ¡Ê f(U) ¡É Ni.1. We have v = f(u) for some u ¡Ê U and fi.1(v) = 0, which means that fi(u) = 0. Then u ¡Ê U ¡É Ni, so u = 0 since U ¡É Ni = (0), and v = f(u) = 0. Therefore, f(U) ¡É Ni.1 = (0). The restriction of f to U is obviously surjective on f(U). Suppose that f(u) = 0 for some u ¡Ê U. Then u ¡Ê U ¡É N1 . U ¡É Ni = (0) (since i ¡Ý 1), so u = 0, which proves that f is also injective on U. 
Proposition 22.21. Given a nilpotent linear map f with fr =0 and fr+1 =0, there exists a sequence of subspace U1,...,Ur+1 of E with the following properties: 
(1) Ni = Ni.1 ¨’ Ui, for i =1,...,r +1. 
(2) We have f(Ui) . Ui.1, and the restriction of f to Ui is an injection, for i =2,...,r+1. See Figure 22.2. 
Proof. We proceed inductively, by de.ning the sequence Ur+1,Ur,...,U1. We pick Ur+1 to be any supplement of Nr in Nr+1 = E, so that 
E = Nr+1 = Nr ¨’ Ur+1. 
Since fr+1 = 0 and Nr = Ker (fr), we have f(Ur+1) . Nr, and by Proposition 22.20, as Ur+1 ¡É Nr = (0), we have f(Ur+1) ¡É Nr.1 = (0). As a consequence, we can pick a supplement Ur of Nr.1 in Nr so that f(Ur+1) . Ur. We have 
Nr = Nr.1 ¨’ Ur and f(Ur+1) . Ur. 
By Proposition 22.20, f is an injection from Ur+1 to Ur. Assume inductively that Ur+1,...,Ui have been de.ned for i ¡Ý 2 and that they satisfy (1) and (2). Since 
Ni = Ni.1 ¨’ Ui, 
we have Ui . Ni, so fi.1(f(Ui)) = fi(Ui) = (0), which implies that f(Ui) . Ni.1. Also, since Ui ¡É Ni.1 = (0), by Proposition 22.20, we have f(Ui) ¡ÉNi.2 = (0). It follows that there is a supplement Ui.1 of Ni.2 in Ni.1 that contains f(Ui). We have 
Ni.1 = Ni.2 ¨’ Ui.1 and f(Ui) . Ui.1. 
The fact that f is an injection from Ui into Ui.1 follows from Proposition 22.20. Therefore, the induction step is proven. The construction stops when i = 1. 

Figure 22.2: A schematic illustration of Ni = Ni.1¨’Ui with f(Ui) . Ui.1 for i = r+1, r, r.1. 
Because N0 = (0) and Nr+1 = E, we see that E is the direct sum of the Ui: E = U1 ¨’¡¤ ¡¤¡¤¨’ Ur+1, with f(Ui) . Ui.1, and f an injection from Ui to Ui.1, for i = r +1,..., 2. By a clever choice of bases in the Ui, we obtain the following nice theorem. Theorem 22.22. For any nilpotent linear map f : E ¡ú E on a .nite-dimensional vector space E of dimension n over a .eld K, there is a basis of E such that the matrix N of f is 
of the form 

.
. 
N = 

...... 

0 ¦Í1 0 ¡¤¡¤¡¤ 00 00 ¦Í2 ¡¤¡¤¡¤ 00 
... ... 
... ... 

.  .  .  .  .  .  
0  0  0  ¡¤ ¡¤ ¡¤  0  ¦Ín  
0  0  0  ¡¤ ¡¤ ¡¤  0  0  

...... 

, 


22.7. NILPOTENT LINEAR MAPS AND JORDAN FORM 
where ¦Íi =1 or ¦Íi =0. 
 

r+1
Proof. First apply Proposition 22.21 to obtain a direct sum E 

=
i=1 
Ui. Then we de.ne 
a basis of E inductively as follows. First we choose a basis 
r+1 r+1 
e ,...,e 
1 nr+1 of Ur+1. Next, for i = r +1,..., 2, given the basis 
ii 
e1,...,e ni 
of Ui, since f is injective on Ui and f(Ui) . Ui.1, the vectors f(e1i ),...,f(eini ) are linearly independent, so we de.ne a basis of Ui.1 by completing f(e1i ),...,f(eni i ) to a basis in Ui.1: 
i.1 i.1 i.1 i.1 
e ,...,e ,e ni+1,...,e 
1 ni ni.1 
with 
i.1 i 
e = f(ej),j =1 ...,ni.
j 
Since U1 = N1 = Ker (f), we have f(ej 1)=0,j =1,...,n1. These basis vectors can be arranged as the rows of the following matrix: 
. 

r+1 r+1
e ¡¤¡¤¡¤ e
1 
nr+1 
.. 
.. 
.. 
rrr r
e¡¤¡¤¡¤ ee¡¤¡¤¡¤ e
1 nr+1 nr+1+1 nr 
... . 
... . 
... . 
r.1 r.1 r.1 r.1 r.1 r.1
e ¡¤¡¤¡¤ ee ¡¤¡¤¡¤ ee ¡¤¡¤¡¤ e
1 nr+1 nr+1+1 nr nr+1 nr.1 
... ... 
... ... 
... ... 
... ... 
... ... 
... ... 
111 111 1
e¡¤¡¤¡¤ ee¡¤¡¤¡¤ ee¡¤¡¤¡¤ e¡¤¡¤¡¤ ¡¤¡¤¡¤ e
1 nr+1 nr+1+1 nr nr+1 nr.1 n1 
. 

.............. 

.............. 

Finally, we de.ne the basis (e1,...,en) by listing each column of the above matrix from the bottom-up, starting with column one, then column two, etc. This means that we list the vectors ei in the following order: 
j 
1 r+1
For j =1,...,nr+1, list ej ,...,e ;
j 
In general, for i = r, . . . , 1, 
for j = ni+1 +1,...,ni, list ej 1,...,eji . 
1 i.1 i

Then because f(ej )=0 and ej = f(ej) for i ¡Ý 2, either 
f(ei)=0 or f(ei)= ei.1, 
which proves the theorem. 

As an application of Theorem 22.22, we obtain the Jordan form of a linear map. De.nition 22.12. A Jordan block is an r ¡Á r matrix Jr(¦Ë), of the form 
.
. 
¦Ë 10 ¡¤¡¤¡¤ 0 0 ¦Ë 1 ¡¤¡¤¡¤ 0 .. .
..
..... 
...... 

...... 

Jr(¦Ë)= 
. 

.

.. 

. 

, 

..
000 .1 00 0 ¡¤¡¤¡¤ ¦Ë 
where ¦Ë ¡Ê K, with J1(¦Ë)=(¦Ë) if r =1. A Jordan matrix, J, is an n ¡Á n block diagonal 
matrix of the form 

.
. 
Jr1 (¦Ë1) ¡¤¡¤¡¤ 0 ..
.
.. 

..

J 

= 

. 

.. 

,

.

. 

. 

0 ¡¤¡¤¡¤ Jrm (¦Ëm) where each Jrk (¦Ëk) is a Jordan block associated with some ¦Ëk ¡Ê K, and with r1 +¡¤¡¤¡¤+rm = n. To simplify notation, we often write J(¦Ë) for Jr(¦Ë). Here is an example of a Jordan matrix with four blocks: 
.
. 
J = 

........... 

¦Ë  1  0  0  0  0  0  0  
0  ¦Ë  1  0  0  0  0  0  
0  0  ¦Ë  0  0  0  0  0  
0  0  0  ¦Ë  1  0  0  0  

0  0  0  0  ¦Ë  0  0  0  
0  0  0  0  0  ¦Ë  0  0  
0  0  0  0  0  0  ¦Ì  1  
0  0  0  0  0  0  0  ¦Ì  
Theorem 22.23. (Jordan form) Let E be a vector space of dimension n over a .eld K and let f : E ¡ú E be a linear map. The following properties are equivalent: 


........... 

. 

(1) The eigenvalues of f all belong to K (i.e. the roots of the characteristic polynomial ¦Öf all belong to K). 
(2) There is a basis of E in which the matrix of f is a Jordan matrix. Proof. Assume (1). First we apply Theorem 22.17, and we get a direct sum E = k Wk,
j=1 
such that the restriction of gi = f . ¦Ëjid to Wi is nilpotent. By Theorem 22.22, there is a basis of Wi such that the matrix of the restriction of gi is of the form 
.
. 
Gi = 
...... 

0 ¦Í1 0 ¡¤¡¤¡¤ 00 00 ¦Í2 ¡¤¡¤¡¤ 00 ... ... 
... ... 

.  .  .  .  .  .  
0  0  0  ¡¤ ¡¤ ¡¤  0  ¦Íni  
0  0  0  ¡¤ ¡¤ ¡¤  0  0  

...... 

, 


22.7. NILPOTENT LINEAR MAPS AND JORDAN FORM 
where ¦Íi = 1 or ¦Íi = 0. Furthermore, over any basis, ¦Ëiid is represented by the diagonal matrix Di with ¦Ëi on the diagonal. Then it is clear that we can split Di + Gi into Jordan blocks by forming a Jordan block for every uninterrupted chain of 1s. By putting the bases of the Wi together, we obtain a matrix in Jordan form for f. 
Now assume (2). If f can be represented by a Jordan matrix, it is obvious that the diagonal entries are the eigenvalues of f, so they all belong to K. 
Observe that Theorem 22.23 applies if K = C. It turns out that there are uniqueness properties of the Jordan blocks but more machinery is needed to prove this result. If a complex n ¡Á n matrix A is expressed in terms of its Jordan decomposition as A = D + N, since D and N commute, by Proposition 8.21, the exponential of A is given by 
A DN 
e = ee, 
and since N is an n ¡Á n nilpotent matrix, Nn.1 = 0, so we obtain 
NN2 Nn.1 
AD 
e = eI ++ + ¡¤¡¤¡¤ + . 
1!2! (n . 1)! 
In particular, the above applies if A is a Jordan matrix. This fact can be used to solve (at least in theory) systems of .rst-order linear di.erential equations. Such systems are of the form 
dX 
= AX, (.)
dt where A is an n ¡Á n matrix and X is an n-dimensional vector of functions of the parameter 
t. 
It can be shown that the columns of the matrix etA form a basis of the vector space of solutions of the system of linear di.erential equations (.); see Artin [3] (Chapter 4). Furthermore, for any matrix B and any invertible matrix P , if A = P BP .1, then the system (.) is equivalent to 
P .1 
dX = BP .1X, 
dt so if we make the change of variable Y = P .1X, we obtain the system dY 
= BY. (..)
dt Consequently, if B is such that the exponential etB can be easily computed, we obtain an explicit solution Y of (..) , and X = PY is an explicit solution of (.). This is the case when B is a Jordan form of A. In this case, it su.ces to consider the Jordan blocks of B. Then 
we have 

.
. 
01 0 ¡¤¡¤¡¤ 0 00 1 ¡¤¡¤¡¤ 0 .. .
..
..... 
...... 

...... 

Jr(¦Ë)= ¦ËIr 
+ 

= ¦ËIr + N, 
. 

.

.. 

. 

..
000 .1 00 0 ¡¤¡¤¡¤ 0 
and the powers Nk are easily computed. For example, if 
.. .. 
310 010 B = .031. =3I3 + .001. 003 000 
we obtain .. .. 310 0 t 0 tB = t .031. =3tI3 + .00 t. 003 000 
and so ..... . 
3t 3t 3t
e00 1 t (1/2)t2 ete3t (1/2)t2etB 3t 3t
e = . 0 e0 ..01 t . = . 0 ete3t .. 
3t 3t
00 e001 00 e
The columns of etB form a basis of the space of solutions of the system of linear di.erential equations 
dY1 
=3Y1 + Y2
dt dY2 
=3Y2 + Y3
dt dY3 
=3Y3,
dt 
in matrix form, 
.. 
dY1 . ... dt 310 Y1
.. 
... ...
dY2 = 031 Y2 .
. dt . dY3 003 Y3 dt 
Explicitly, the general solution of the above system is 
....... . 
3t 3t
Y1 ete3t (1/2)t2e.Y2. = c1 . 0 . + c2 . e3t . + c3 . te3t ., 
Y3 00 e3t 
with c1,c2,c3 ¡Ê R. Solving systems of .rst-order linear di.erential equations is discussed in Artin [3] and more extensively in Hirsh and Smale [35]. 

22.8 Summary 
The main concepts and results of this chapter are listed below: 
. 
Ideals, principal ideals, greatest common divisors. 

22.9. PROBLEMS 

. 
Monic polynomial, irreducible polynomial, relatively prime polynomials. 

. 
Annihilator of a linear map. 

. 
Minimal polynomial of a linear map. 

. 
Invariant subspace. 

. 
f-conductor of u into W ; conductor of u into W . 

. 
Diagonalizable linear maps. 

. 
Commuting families of linear maps. 

. 
Primary decomposition. 

. 
Generalized eigenvectors. 

. 
Nilpotent linear map. 

. 
Normal form of a nilpotent linear map. 

. 
Jordan decomposition. 

. 
Jordan block. 

. 
Jordan matrix. 

. 
Jordan normal form. 

. 
Systems of .rst-order linear di.erential equations. 


22.9 Problems 
Problem 22.1. Prove that the minimal monic polynomial of Proposition 22.1 is unique. 
Problem 22.2. Given a linear map f : E ¡ú E, prove that the set Ann(f) of polynomials 
that annihilate f is an ideal. 
Problem 22.3. Provide the details of Proposition 22.9. 
Problem 22.4. Prove that the f-conductor Sf (u, W ) is an ideal in K[X] (Proposition 

22.10). 
Problem 22.5. Prove that the polynomials g1,...,gk used in the proof of Theorem 22.16 
are relatively prime. 

.
. 
A =

. 

6 .3 .2 
4 .1 .2

.

. 

10 .5 .3 Problem 22.7. Find the Jordan decomposition of the matrix 
.
. 
A =

. 

31 .1 
22 .1

.

. 

22 0 
Problem 22.8. Let f : E ¡ú E be a linear map on a .nite-dimensional vector space. Prove that if f has rank 1, then either f is diagonalizable or f is nilpotent but not both. Problem 22.9. Find the Jordan form of the matrix 
.
. 
...
3 0000 
Problem 22.10. Let N be a 3 ¡Á 3 nilpotent matrix over C. Prove that the matrix A = I + (1/2)N . (1/8)N2 satis.es the equation 
A2 
= I + N. 
In other words, A is a square root of I + N. 
Generalize the above fact to any n ¡Á n nilpotent matrix N over C using the binomial series for (1 + t)1/2 . 
Problem 22.11. Let K be an algebraically closed .eld (for example, K = C). Prove that every 4 ¡Á 4 matrix is similar to a Jordan matrix of the following form: 
... 

0100 0
002 

A = 

. 

000 

.
.
.
.
.
. 
¦Ë1 000 ¦Ë 10 0 ¦Ë 10 0 
... 

0 ¦Ë2 00 
00 ¦Ë3 0 
...

, 

... 

0 ¦Ë 00 

00 ¦Ë3 0 
...

, 

... 

0 ¦Ë 10 

00 ¦Ë 0 

...

, 

000 ¦Ë4 00 0 ¦Ë4 000 ¦Ë4 
.
.
.
. 
¦Ë 100 ¦Ë 100 

... 

0 ¦Ë 00 

00 ¦Ì 1 

...

, 

... 

0 ¦Ë 10 

00 ¦Ë 1 

...

. 

000 ¦Ì 000 ¦Ë 

22.9. PROBLEMS 
Problem 22.12. In this problem the .eld K is of characteristic 0. Consider an (r ¡Á r) 
Jordan block 

.
. 
¦Ë 10 ¡¤¡¤¡¤ 0 0 ¦Ë 1 ¡¤¡¤¡¤ 0 .. .
..
.....Jr(¦Ë)= 
...... 

...... 

.

. 

.

.. 

. 

..
000 .1 00 0 ¡¤¡¤¡¤ ¦Ë 
Prove that for any polynomial f(X), we have 

.
. 
f(¦Ë) f1(¦Ë) f2(¦Ë) ¡¤¡¤¡¤ fr.1(¦Ë) 
0 f(¦Ë) f1(¦Ë) ¡¤¡¤¡¤ fr.2(¦Ë) 
.. .
..
.. .. . 
...... 

...... 

f(Jr(¦Ë)) = 
. 

.

.. 

. 

, 

..
000 . f1(¦Ë) 00 0 ¡¤¡¤¡¤ f(¦Ë) 
where f(k)(X)fk(X)= ,
k! 
and f(k)(X) is the kth derivative of f(X). 


