Chapter 8 Vector Norms and Matrix Norms 
8.1 Normed Vector Spaces 
In order to de.ne how close two vectors or two matrices are, and in order to de.ne the convergence of sequences of vectors or matrices, we can use the notion of a norm. Recall that R+ = {x ¡ÊR |x ¡Ý0}. Also recall that if z = a + ib ¡ÊC is a complex number, with 
¡Ì¡Ì 
a, b ¡ÊR, then z = a .ib and |z|= zz = a2 + b2 (|z|is the modulus of z). 
De.nition 8.1. Let E be a vector space over a .eld K, where K is either the .eld R of 
reals, or the .eld C of complex numbers. A norm on E is a function ll: E ¡úR+, assigning 

a nonnegative real number lulto any vector u ¡ÊE, and satisfying the following conditions 
for all x, y, z ¡ÊE and ¦Ë ¡ÊK: 
(N1) lxl¡Ý0, and lxl= 0 i. x =0. (positivity) 
(N2) l¦Ëxl= |¦Ë|lxl. (homogeneity (or scaling)) 
(N3) lx + yl¡Ülxl+ lyl. (triangle inequality) 

A vector space E together with a norm llis called a normed vector space. By (N2), setting ¦Ë = .1, we obtain l.xl= l(.1)xl= |.1|lxl= lxl; that is, l.xl= lxl. From (N3), we have lxl= lx .y + yl¡Ülx .yl+ lyl, which implies that 
lxl.lyl¡Ülx .yl. By exchanging x and y and using the fact that by (N2), ly .xl= l.(x .y)l= lx .yl, 
271 

we also have 
lyl.lxl¡Ülx . yl . 
Therefore, |lxl.lyl|¡Ü lx . yl, for all x, y ¡Ê E. (.) 
Observe that setting ¦Ë = 0 in (N2), we deduce that l0l = 0 without assuming (N1). Then by setting y =0 in (.), we obtain 
|lxl| ¡Ü lxl , for all x ¡Ê E. 
Therefore, the condition lxl¡Ý 0 in (N1) follows from (N2) and (N3), and (N1) can be replaced by the weaker condition 
(N1¡¯) For all x ¡Ê E, if lxl = 0, then x = 0, 
A function ll : E ¡ú R satisfying Axioms (N2) and (N3) is called a seminorm. From the 
above discussion, a seminorm also has the properties 
lxl¡Ý 0 for all x ¡Ê E, and l0l = 0. 
However, there may be nonzero vectors x ¡Ê E such that lxl = 0. 
Let us give some examples of normed vector spaces. 

Example 8.1. 
1. 
Let E = R, and lxl = |x|, the absolute value of x. 

2. 
Let E = C, and lzl = |z|, the modulus of z. 

3. 
Let E = Rn (or E = Cn). There are three standard norms. For every (x1,...,xn) ¡Ê E, we have the norm lxl1, de.ned such that, 


lxl1 = |x1| + ¡¤¡¤¡¤ + |xn|, 
we have the Euclidean norm lxl2, de.ned such that, 
  1 
2
lxl2 =|x1|2 + ¡¤¡¤¡¤ + |xn|2, 
and the sup-norm lxl¡Þ, de.ned such that, 
lxl¡Þ = max{|xi|| 1 ¡Ü i ¡Ü n}. 
More generally, we de.ne the fp-norm (for p ¡Ý 1) by 
lxlp =(|x1|p + ¡¤¡¤¡¤ + |xn|p)1/p . 
See Figures 8.1 through 8.4. There are other norms besides the fp-norms. Here are some examples. 
8.1. NORMED VECTOR SPACES 

1. 
For E = R2 , 
l(u1,u2)l = |u1| +2|u2|. 


See Figure 8.5. 

2. 
For E = R2 , 

l(u1,u2)l =(u1 + u2)2 + u 21 1/2 . 
See Figure 8.6. 


3. 
For E = C2 , 
l(u1,u2)l = |u1 + iu2| + |u1 . iu2|. 



The reader should check that they satisfy all the axioms of a norm. Some work is required to show the triangle inequality for the fp-norm. 
Proposition 8.1. If E = Cn or E = Rn, for every real number p ¡Ý 1, the fp-norm is indeed a norm. 

Proof. The cases p = 1 and p = ¡Þ are easy and left to the reader. If p> 1, then let q> 1 such that 
11 
+ =1. 
pq We will make use of the following fact: for all ¦Á, ¦Â ¡Ê R, if ¦Á, ¦Â ¡Ý 0, then 
¦Áp  ¦Âq  
¦Á¦Â ¡Ü  +  .  (.)  
p  q  

To prove the above inequality, we use the fact that the exponential function t ¡ú et satis.es the following convexity inequality: 
y
e ¦Èx+(1.¦È)y ¡Ü ¦Èex + (1 . ¦È)e, 
for all x, y ¡Ê R and all ¦È with 0 ¡Ü ¦È ¡Ü 1. 

8.1. NORMED VECTOR SPACES 

Since the case ¦Á¦Â = 0 is trivial, let us assume that ¦Á> 0 and ¦Â> 0. If we replace ¦È by 1/p, x by p log ¦Á and y by q log ¦Â, then we get 
11 
11 
p p log ¦Á+ q q log ¦Âp log ¦Áq log ¦Â 
e ¡Ü e+ e, 
pq 
which simpli.es to 
¦Áp ¦Âq 
¦Á¦Â ¡Ü + , 
pq 
as claimed. 
We will now prove that for any two vectors u, v ¡Ê E, (where E is of dimension n), we have 
n
n 
|uivi|¡Ülullvl . (..)
pq i=1 
Since the above is trivial if u = 0 or v = 0, let us assume that u 0 and v= 0. Then
=  

Inequality (.) with ¦Á = |ui|/ lul and ¦Â = |vi|/ lvl yields
pq |uivi||ui|p |vi|q 
¡Ü + ,
lullvl p lulp q lulq 
pq pq 
for i =1,...,n, and by summing up these inequalities, we get 
n
n|uivi|¡Ülullvl ,
pq i=1 
as claimed. To .nish the proof, we simply have to prove that property (N3) holds, since (N1) and (N2) are clear. For i =1,...,n, we can write (|ui| + |vi|)p = |ui|(|ui| + |vi|)p.1 + |vi|(|ui| + |vi|)p.1 , so that by summing up these equations we get 
n
n
n
nnn
(|ui| + |vi|)p = |ui|(|ui| + |vi|)p.1 + |vi|(|ui| + |vi|)p.1 , 
i=1 i=1 i=1 
8.1. NORMED VECTOR SPACES 

and using Inequality (..), with V ¡Ê E where Vi =(|ui| + |vi|)p.1, we get 
n
n 
(|ui| + |vi|)p ¡ÜlullV l + lvllV l 
pq pq i=1  n 1/q
n 
(|ui| + |vi|)(p.1)q
=(lul p + lvl p). 
i=1 
However, 1/p +1/q = 1 implies pq = p + q, that is, (p . 1)q = p, so we have 
n n 1/q
nn 
(|ui| + |vi|)p ¡Ü (lul + lvl )(|ui| + |vi|)p,
ppi=1 i=1 
which yields
 n 1.1/q  n 1/p
nn 
(|ui| + |vi|)p=(|ui| + |vi|)p¡Ülul + lvl . 
pp i=1 i=1 
Since |ui + vi|¡Ü|ui| + |vi|, the above implies the triangle inequality lu + vl ¡Ülul + lvl ,
p pp
as claimed. 
For p> 1 and 1/p +1/q = 1, the inequality 
n n 1/p n 1/q
nn n 
|uivi|¡Ü|ui|p|vi|qi=1 i=1 i=1 
is known as H¡§older¡¯s inequality. For p = 2, it is the Cauchy¨CSchwarz inequality. 


Actually, if we de.ne the Hermitian inner product£¨.,.£© on Cn by
n
n
£¨u, v£© = uivi, 
i=1 
where u =(u1,...,un) and v =(v1,...,vn), then
n
n
nn|£¨u, v£©|¡Ü |uivi| = |uivi|, 
i=1 i=1 
so H¡§older¡¯s inequality implies the following inequalities. 
Corollary 8.2. (H¡§older¡¯s inequalities) For any real numbers p, q, such that p, q ¡Ý 1 and 
11 
+ =1, 
pq (with q =+¡Þ if p =1 and p =+¡Þ if q =1), we have the inequalities 
n
n
n
nnn|uivi|¡Ü |ui|p |vi|q i=1 i=1 i=1 
and
|£¨u, v£©| ¡Ül ullvl , u,v ¡Ê Cn . 
pq 
For p = 2, this is the standard Cauchy¨CSchwarz inequality. The triangle inequality for the fp-norm, 
1/p 1/q 
n
n
n
nnn
(|ui + vi|)p ¡Ü|ui|p + |vi|q , 
i=1 i=1 i=1 
1/p 1/p 1/q 
8.1. NORMED VECTOR SPACES 
is known as Minkowski¡¯s inequality. When we restrict the Hermitian inner product to real vectors, u, v ¡Ê Rn , we get the 
Euclidean inner product
n
n 
£¨u, v£© = uivi. i=1 
It is very useful to observe that if we represent (as usual) u =(u1,...,un) and v =(v1,...,vn) (in Rn) by column vectors, then their Euclidean inner product is given by
£¨u, v£© = uTv = vTu, 
and when u, v ¡Ê Cn, their Hermitian inner product is given by
.
£¨u, v£© = v . u = u v. 
In particular, when u = v, in the complex case we get 
. 
lul22 = u u, 
and in the real case this becomes 
T
lul22 = uu. 
As convenient as these notations are, we still recommend that you do not abuse them; the notation£¨u, v£© is more intrinsic and still ¡°works¡± when our vector space is in.nite dimen-sional. 
Remark: If 0 <p< 1, then x ¡úlxl p is not a norm because the triangle inequality 
fails. For example, consider x = (2, 0) and y = (0, 2). Then x + y = (2, 2), and we have 2(p+1)/p
lxl = (2p +0p)1/p = 2, lyl = (0p +2p)1/p = 2, and lx + yl = (2p +2p)1/p =. 
pp p 
Thus lx + yl =2(p+1)/p, lxl + lyl =4=22 . 
p pp 
Since 0 <p< 1, we have 2p<p + 1, that is, (p + 1)/p > 2, so 2(p+1)/p > 22 = 4, and the triangle inequality lx + yl ¡Ülxl + lyl fails. 
p pp 
Observe that 
l(1/2)xl = (1/2) lxl = l(1/2)yl = (1/2) lyl =1,
p ppp 
l(1/2)(x + y)l p =21/p, 
and since p< 1, we have 21/p > 2, so 
l(1/2)(x + y)l p =21/p > 2 = (1/2) lxl p + (1/2) lyl p , 
and the map x ¡úlxl p is not convex. 
For p = 0, for any x ¡Ê Rn, we have lxl= |{i ¡Ê{1,...,n}| xi =0}|,
0 
the number of nonzero components of x. The map x ¡úlxl0 is not a norm this time because Axiom (N2) fails. For example, 
l(1, 0)l0 = l(10, 0)l0 = 1 = 10 = 10 l(1, 0)l0 . 
The map x ¡úlxl0 is also not convex. For example, 
l(1/2)(2, 2)l0 = l(1, 1)l0 =2, 
and l(2, 0)l0 = l(0, 2)l0 =1, but l(1/2)(2, 2)l0 =2 > 1 = (1/2) l(2, 0)l0 + (1/2) l(0, 2)l0 . 
Nevertheless, the ¡°zero-norm¡± x ¡úlxl0 is used in machine learning as a regularizing term which encourages sparsity, namely increases the number of zero components of the vector x. 
The following proposition is easy to show. 
Proposition 8.3. The following inequalities hold for all x ¡Ê Rn (or x ¡Ê Cn): 
lxl¡Þ ¡Ülxl1 ¡Ü nlxl¡Þ, 
¡Ì 
lxl¡Þ ¡Ülxl2 ¡Ü nlxl¡Þ, 
¡Ì 
lxl2 ¡Ülxl1 ¡Ü nlxl2. 
Proposition 8.3 is actually a special case of a very important result: in a .nite-dimensional vector space, any two norms are equivalent. 
De.nition 8.2. Given any (real or complex) vector space E, two norms ll a and llb are equivalent i. there exists some positive reals C1,C2 > 0, such that 
lul¡Ü C1 luland lul¡Ü C2 lul , for all u ¡Ê E. 
ab ba 
Given any norm ll on a vector space of dimension n, for any basis (e1,...,en) of E, observe that for any vector x = x1e1 + ¡¤¡¤¡¤ + xnen, we have 
lxl = lx1e1 + ¡¤¡¤¡¤ + xnenl¡Ü|x1|le1l + ¡¤¡¤¡¤ + |xn|lenl ¡Ü C(|x1| + ¡¤¡¤¡¤ + |xn|)= C lxl1 , 
with C = max1¡Üi¡Ün leil and with the norm lxl1 de.ned as 
lxl1 = lx1e1 + ¡¤¡¤¡¤ + xnenl = |x1| + ¡¤¡¤¡¤ + |xn|. The above implies that 
|lul.lvl| ¡Ü lu . vl¡Ü C lu . vl1 , 
and this implies the following corollary. 

8.1. NORMED VECTOR SPACES 
Corollary 8.4. For any norm u ¡úlul on a .nite-dimensional (complex or real) vector space E, the map u ¡úlul is continuous with respect to the norm ll1. 
Let S1 n.1 be the unit sphere with respect to the norm ll1, namely 
Sn.1 = {x ¡Ê E |lxl=1}.
11 
Now S1 n.1 is a closed and bounded subset of a .nite-dimensional vector space, so by Heine¨C Borel (or equivalently, by Bolzano¨CWeiertrass), S1 n.1 is compact. On the other hand, it is a well known result of analysis that any continuous real-valued function on a nonempty compact set has a minimum and a maximum, and that they are achieved. Using these facts, we can prove the following important theorem: 
Theorem 8.5. If E is any real or complex vector space of .nite dimension, then any two norms on E are equivalent. 
Proof. It is enough to prove that any norm ll is equivalent to the 1-norm. We already proved that the function x ¡úlxl is continuous with respect to the norm ll1, and we observed that the unit sphere S1 n.1 is compact. Now we just recalled that because the function f : x ¡úlxl is continuous and because S1 n.1 is compact, the function f has a minimum m and a maximum M, and because lxl is never zero on S1 n.1, we must have m> 0. Consequently, we just proved that if lxl1 = 1, then 
0 <m ¡Ülxl¡Ü M, 
so for any x ¡Ê E with x = 0, we get 
m ¡Ülx/ lxl1l¡Ü M, 
which implies 
m lxl1 ¡Ülxl¡Ü M lxl1 . 
Since the above inequality holds trivially if x = 0, we just proved that ll and ll1 are equivalent, as claimed. 
Remark: Let P be a n ¡Á n symmetric positive de.nite matrix. It is immediately veri.ed that the map x ¡úlxlP given by 
lxlP =(x TPx)1/2 
is a norm on Rn called a quadratic norm. Using some convex analysis (the L¡§owner¨CJohn ellipsoid), it can be shown that any norm ll on Rn can be approximated by a quadratic norm in the sense that there is a quadratic norm llP such that 
¡Ì 
lxlP ¡Ülxl¡Ü n lxlP for all x ¡Ê Rn; 
see Boyd and Vandenberghe [11], Section 8.4.1. Next we will consider norms on matrices. 


8.2 Matrix Norms 
For simplicity of exposition, we will consider the vector spaces Mn(R) and Mn(C) of square n ¡Á n matrices. Most results also hold for the spaces Mm,n(R) and Mm,n(C) of rectangular m ¡Á n matrices. Since n ¡Á n matrices can be multiplied, the idea behind matrix norms is that they should behave ¡°well¡± with respect to matrix multiplication. 
De.nition 8.3. A matrix norm ll on the space of square n ¡Á n matrices in Mn(K), with K = R or K = C, is a norm on the vector space Mn(K), with the additional property called submultiplicativity that 
lABl¡ÜlAllBl , 
for all A, B ¡Ê Mn(K). A norm on matrices satisfying the above property is often called a submultiplicative matrix norm. 
Since I2 = I, from lIl = lI2l¡ÜlIl2 , we get lIl¡Ý 1, for every matrix norm. 
Before giving examples of matrix norms, we need to review some basic de.nitions about matrices. Given any matrix A =(aij) ¡Ê Mm,n(C), the conjugate A of A is the matrix such that 
Aij = aij, 1 ¡Ü i ¡Ü m, 1 ¡Ü j ¡Ü n. The transpose of A is the n ¡Á m matrix AT such that 
AT = aji, 1 ¡Ü i ¡Ü m, 1 ¡Ü j ¡Ü n.
ij 
The adjoint of A is the n ¡Á m matrix A. such that 
A . =(AT)=(A)T . 
When A is a real matrix, A. = AT . A matrix A ¡Ê Mn(C) is Hermitian if 
A . 
= A. 
If A is a real matrix (A ¡Ê Mn(R)), we say that A is symmetric if 
AT 
= A. 
A matrix A ¡Ê Mn(C) is normal if AA . = A . A, 
and if A is a real matrix, it is normal if 
AAT = ATA. 
A matrix U ¡Ê Mn(C) is unitary if 
UU . = U . U = I. 
8.2. MATRIX NORMS 
A real matrix Q ¡Ê Mn(R) is orthogonal if 
QQT = QTQ = I. 
Given any matrix A =(aij) ¡Ê Mn(C), the trace tr(A) of A is the sum of its diagonal elements tr(A)= a11 + ¡¤¡¤¡¤ + ann. 
It is easy to show that the trace is a linear map, so that 
tr(¦ËA)= ¦Ëtr(A) 
and tr(A + B) = tr(A) + tr(B). 
Moreover, if A is an m ¡Á n matrix and B is an n ¡Á m matrix, it is not hard to show that 
tr(AB) = tr(BA). 
We also review eigenvalues and eigenvectors. We content ourselves with de.nition in-volving matrices. A more general treatment will be given later on (see Chapter 14). 
De.nition 8.4. Given any square matrix A ¡Ê Mn(C), a complex number ¦Ë ¡Ê C is an eigenvalue of A if there is some nonzero vector u ¡Ê Cn, such that 
Au = ¦Ëu. 
If ¦Ë is an eigenvalue of A, then the nonzero vectors u ¡Ê Cn such that Au = ¦Ëu are called eigenvectors of A associated with ¦Ë; together with the zero vector, these eigenvectors form a subspace of Cn denoted by E¦Ë(A), and called the eigenspace associated with ¦Ë. 
Remark: Note that De.nition 8.4 requires an eigenvector to be nonzero. A somewhat unfortunate consequence of this requirement is that the set of eigenvectors is not a subspace, since the zero vector is missing! On the positive side, whenever eigenvectors are involved, there is no need to say that they are nonzero. The fact that eigenvectors are nonzero is implicitly used in all the arguments involving them, so it seems safer (but perhaps not as elegant) to stipulate that eigenvectors should be nonzero. 
If A is a square real matrix A ¡Ê Mn(R), then we restrict De.nition 8.4 to real eigenvalues ¦Ë ¡Ê R and real eigenvectors. However, it should be noted that although every complex matrix always has at least some complex eigenvalue, a real matrix may not have any real 
eigenvalues. For example, the matrix  
A =  0 1  .1 0  

has the complex eigenvalues i and .i, but no real eigenvalues. Thus, typically even for real matrices, we consider complex eigenvalues. 
Observe that ¦Ë ¡Ê C is an eigenvalue of A 
. 
i. Au = ¦Ëu for some nonzero vector u ¡Ê Cn 

. 
i. (¦ËI . A)u =0 

. 
i. the matrix ¦ËI . A de.nes a linear map which has a nonzero kernel, that is, 


. i. ¦ËI . A not invertible. However, from Proposition 6.11, ¦ËI . A is not invertible i. det(¦ËI . A)=0. Now det(¦ËI . A) is a polynomial of degree n in the indeterminate ¦Ë, in fact, of the form ¦Ën . tr(A)¦Ën.1 + ¡¤¡¤¡¤ +(.1)n det(A). Thus we see that the eigenvalues of A are the zeros (also called roots) of the above polyno-
mial. Since every complex polynomial of degree n has exactly n roots, counted with their multiplicity, we have the following de.nition: De.nition 8.5. Given any square n ¡Á n matrix A ¡Ê Mn(C), the polynomial 
det(¦ËI . A)= ¦Ën . tr(A)¦Ën.1 + ¡¤¡¤¡¤ +(.1)n det(A) is called the characteristic polynomial of A. The n (not necessarily distinct) roots ¦Ë1,...,¦Ën of the characteristic polynomial are all the eigenvalues of A and constitute the spectrum of 
A. We let ¦Ñ(A) = max |¦Ëi|
1¡Üi¡Ün 
be the largest modulus of the eigenvalues of A, called the spectral radius of A. Since the eigenvalue ¦Ë1,...,¦Ën of A are the zeros of the polynomial det(¦ËI . A)= ¦Ën . tr(A)¦Ën.1 + ¡¤¡¤¡¤ +(.1)n det(A), we deduce (see Section 14.1 for details) that tr(A)= ¦Ë1 + ¡¤¡¤¡¤ + ¦Ën det(A)= ¦Ë1 ¡¤¡¤¡¤ ¦Ën. Proposition 8.6. For any matrix norm ll on Mn(C) and for any square n ¡Á n matrix A ¡Ê Mn(C), we have ¦Ñ(A) ¡ÜlAl . 

8.2. MATRIX NORMS 
Proof. Let ¦Ë be some eigenvalue of A for which |¦Ë| is maximum, that is, such that |¦Ë| = ¦Ñ(A). If u (= 0) is any eigenvector associated with ¦Ë and if U is the n ¡Á n matrix whose columns are all u, then Au = ¦Ëu implies 
AU = ¦ËU, 
and since |¦Ë|lUl = l¦ËUl = lAUl¡ÜlAllUl 
and U = 0, we have lUl = 0, and get 
¦Ñ(A)= |¦Ë|¡ÜlAl , 
as claimed. 
Proposition 8.6 also holds for any real matrix norm ll on Mn(R) but the proof is more subtle and requires the notion of induced norm. We prove it after giving De.nition 8.7. 
It turns out that if A is a real n ¡Á n symmetric matrix, then the eigenvalues of A are all real and there is some orthogonal matrix Q such that 
A = Qdiag(¦Ë1,...,¦Ën)QT , 
where diag(¦Ë1,...,¦Ën) denotes the matrix whose only nonzero entries (if any) are its diagonal entries, which are the (real) eigenvalues of A. Similarly, if A is a complex n ¡Á n Hermitian matrix, then the eigenvalues of A are all real and there is some unitary matrix U such that 
A = Udiag(¦Ë1,...,¦Ën)U . , 
where diag(¦Ë1,...,¦Ën) denotes the matrix whose only nonzero entries (if any) are its diagonal entries, which are the (real) eigenvalues of A. See Chapter 16 for the proof of these results. 
We now return to matrix norms. We begin with the so-called Frobenius norm, which is just the norm ll2 on Cn2 , where the n ¡Á n matrix A is viewed as the vector obtained by concatenating together the rows (or the columns) of A. The reader should check that for any n ¡Á n complex matrix A =(aij), 
n1/2  
n  
|aij|2 =tr(A.A) =tr(AA.). 
i,j=1 
De.nition 8.6. The Frobenius norm llF is de.ned so that for every square n ¡Á n matrix A ¡Ê Mn(C), 
n1/2
n   
lAlF = |aij|2 =tr(AA.) =tr(A.A). 
i,j=1 
The following proposition show that the Frobenius norm is a matrix norm satisfying other nice properties. 
Proposition 8.7. The Frobenius norm llF on Mn(C) satis.es the following properties: 
(1) 
It is a matrix norm; that is, lABl¡ÜlAllBl, for all A, B ¡Ê Mn(C).

F FF 

(2) 
It is unitarily invariant, which means that for all unitary matrices U, V , we have 
lAl= lUAl= lAV l= lUAV lF .



FFF 
¡Ì 
(3) ¦Ñ(A.A) ¡ÜlAlF ¡Ü n¦Ñ(A.A), for all A ¡Ê Mn(C). 
Proof. (1) The only property that requires a proof is the fact lABlF ¡ÜlAlF lBlF . This follows from the Cauchy¨CSchwarz inequality: 
½Ð½Ð½Ð½Ð 
½Ð½Ð
i,j=1k=1 
nnnn
2 
nn 
½Ð
½Ð
lABl2 F 
aikbkj
nn¡Ü|aih|2 |bkj|2 i,j=1 h=1 k=1 
n
n 
= 

nn nn  
=  i,h=1 |aih|2  k,j=1 |bkj|2  = lAl2 F  lBl2 F .  
(2) We have  

lAl2 F = tr(A . A) = tr(VV . A . A) = tr(V . A . AV )= lAV l2 F , 
and lAl2 F = tr(A . A) = tr(A . U . UA)= lUAl2 F . 
The identity 
lAlF = lUAV lF 
follows from the previous two. 
(3) It is well known that the trace of a matrix is equal to the sum of its eigenvalues. Furthermore, A.A is symmetric positive semide.nite (which means that its eigenvalues are nonnegative), so ¦Ñ(A.A) is the largest eigenvalue of A.A and 
¦Ñ(A . A) ¡Ü tr(A . A) ¡Ü n¦Ñ(A . A), 
which yields (3) by taking square roots. 
Remark: The Frobenius norm is also known as the Hilbert-Schmidt norm or the Schur norm. So many famous names associated with such a simple thing! 
8.3. SUBORDINATE NORMS 
8.3 Subordinate Norms 
We now give another method for obtaining matrix norms using subordinate norms. First we need a proposition that shows that in a .nite-dimensional space, the linear map induced by a matrix is bounded, and thus continuous. 
Proposition 8.8. For every norm ll on Cn (or Rn), for every matrix A ¡Ê Mn(C) (or A ¡Ê Mn(R)), there is a real constant CA ¡Ý 0, such that 
lAul¡Ü CA lul , 
for every vector u ¡Ê Cn (or u ¡Ê Rn if A is real). 
Proof. For every basis (e1,...,en) of Cn (or Rn), for every vector u = u1e1 + ¡¤¡¤¡¤ + unen, we have 
lAul = lu1A(e1)+ ¡¤¡¤¡¤ + unA(en)l 
¡Ü|u1|lA(e1)l + ¡¤¡¤¡¤ + |un|lA(en)l 
¡Ü C1(|u1| + ¡¤¡¤¡¤ + |un|)= C1 lul1 , 
where C1 = max1¡Üi¡Ün lA(ei)l. By Theorem 8.5, the norms ll and ll1 are equivalent, so there is some constant C2 > 0 so that lul1 ¡Ü C2 lul for all u, which implies that 
lAul¡Ü CA lul , 
where CA = C1C2. 
Proposition 8.8 says that every linear map on a .nite-dimensional space is bounded. This implies that every linear map on a .nite-dimensional space is continuous. Actually, it is not hard to show that a linear map on a normed vector space E is bounded i. it is continuous, regardless of the dimension of E. 
Proposition 8.8 implies that for every matrix A ¡Ê Mn(C) (or A ¡Ê Mn(R)), 
lAxl 
sup ¡Ü CA. 
x¡ÊCn lxl 
x 
=0 
Since l¦Ëul = |¦Ë|lul, for every nonzero vector x, we have 
lAxllxllA(x/ lxl)l 
== lA(x/ lxl)l ,
lxllxl 
which implies that 
lAxl 
sup = sup lAxl . 
x¡ÊCn lxl x¡ÊCn x lxl=1
=0 
Similarly 
lAxl 
sup = sup lAxl . 
x¡ÊRn lxl x¡ÊRn x=0 lxl=1 
The above considerations justify the following de.nition. 
De.nition 8.7. If ll is any norm on Cn, we de.ne the function ll op on Mn(C) by 
lAxl 
lAl op = sup = sup lAxl . 
x¡ÊCn lxl x¡ÊCn x=0 lxl=1 
The function A ¡úlAl op is called the subordinate matrix norm or operator norm induced by the norm ll. 
Another notation for the operator norm of a matrix A (in particular, used by Horn and Johnson [36]), is |
A|
. It is easy to check that the function A ¡úlAl op is indeed a norm, and by de.nition, it satis.es the property lAxl¡ÜlAl op lxl , for all x ¡Ê Cn . 
A norm ll op on Mn(C) satisfying the above property is said to be subordinate to the vector norm ll on Cn . As a consequence of the above inequality, we have 
lABxl¡ÜlAllBxl¡ÜlAllBllxl ,
op opop 
for all x ¡Ê Cn, which implies that 
lABl ¡ÜlAllBl for all A, B ¡Ê Mn(C),
op opop 
showing that A ¡úlAl op is a matrix norm (it is submultiplicative). Observe that the operator norm is also de.ned by 
lAl = inf{¦Ë ¡Ê R |lAxl¡Ü ¦Ë lxl , for all x ¡Ê Cn}. 
op 
Since the function x ¡úlAxl is continuous (because |lAyl.lAxl| ¡Ü lAy . Axl¡Ü CA lx . yl) and the unit sphere Sn.1 = {x ¡Ê Cn |lxl =1} is compact, there is some x ¡Ê Cn such that lxl = 1 and 
lAxl = lAl . 
op 
Equivalently, there is some x ¡Ê Cn such that x = 0 and 
lAxl = lAllxl . 
op 
The de.nition of an operator norm also implies that 
lIl =1. 
op 
8.3. SUBORDINATE NORMS 
The above shows that the Frobenius norm is not a subordinate matrix norm (why?). 
If ll is a vector norm on Cn, the operator norm ll op that it induces applies to matrices in Mn(C). If we are careful to denote vectors and matrices so that no confusion arises, for example, by using lower case letters for vectors and upper case letters for matrices, it should be clear that lAl op is the operator norm of the matrix A and that lxl is the vector norm of 
x. Consequently, following common practice to alleviate notation, we will drop the subscript 
¡°op¡± and simply write lAl instead of lAl op . The notion of subordinate norm can be slightly generalized. 
De.nition 8.8. If K = R or K = C, for any norm ll on Mm,n(K), and for any two norms ll on Kn and llon Km, we say that the norm ll is subordinate to the norms ll and 
a  b  a  
l lb if  
lAxlb ¡Ü lAl lxl a  for all A ¡Ê Mm,n(K) and all x ¡Ê Kn .  

Remark: For any norm ll on Cn, we can de.ne the function llR on Mn(R) by 
lAxl 
lAlR = sup = sup lAxl . 
x¡ÊRn lxl x¡ÊRn x=0 lxl=1 
The function A ¡úlAlR is a matrix norm on Mn(R), and 
lAlR ¡ÜlAl , 
for all real matrices A ¡Ê Mn(R). However, it is possible to construct vector norms ll on Cn and real matrices A such that 
lAlR < lAl . 
In order to avoid this kind of di.culties, we de.ne subordinate matrix norms over Mn(C). Luckily, it turns out that lAlR = lAl for the vector norms, ll1 , ll2, and ll¡Þ. 
We now prove Proposition 8.6 for real matrix norms. 
Proposition 8.9. For any matrix norm ll on Mn(R) and for any square n ¡Á n matrix A ¡Ê Mn(R), we have 
¦Ñ(A) ¡ÜlAl . 
Proof. We follow the proof in Denis Serre¡¯s book [57]. If A is a real matrix, the problem is that the eigenvectors associated with the eigenvalue of maximum modulus may be complex. We use a trick based on the fact that for every matrix A (real or complex), 
¦Ñ(Ak)=(¦Ñ(A))k , 
which is left as an exercise (use Proposition 14.7 which shows that if (¦Ë1,...,¦Ën) are the (not necessarily distinct) eigenvalues of A, then (¦Ëk 1,...,¦Ëkn) are the eigenvalues of Ak, for k ¡Ý 1). 
Pick any complex matrix norm ll c on Cn (for example, the Frobenius norm, or any subordinate matrix norm induced by a norm on Cn). The restriction of ll c to real matrices is a real norm that we also denote by ll c . Now by Theorem 8.5, since Mn(R) has .nite dimension n2, there is some constant C> 0 so that 
lBl c ¡Ü C lBl , for all B ¡Ê Mn(R). 
Furthermore, for every k ¡Ý 1 and for every real n ¡Á n matrix A, by Proposition 8.6, ¦Ñ(Ak) ¡Ü Ak Ak ¡ÜlAlk
Ak 
llllll
llcll
, and because ll is a matrix norm, 

, so we have 

ll

ll

ll

¡Ü C lAlk
(¦Ñ(A))k = ¦Ñ(Ak) ¡Ü 
Ak
¡Ü C 

,

c 
for all k ¡Ý 1. It follows that ¦Ñ(A) ¡Ü C1/k lAl , 
for all k ¡Ý 1. However because C> 0, we have limk ¡ú¡Þ C1/k = 1 (we have limk ¡ú¡Þ k 1 log(C) = 0). There-fore, we conclude that ¦Ñ(A) ¡ÜlAl , as desired. 
We now determine explicitly what are the subordinate matrix norms associated with the vector norms ll1 , ll2, and ll¡Þ. Proposition 8.10. For every square matrix A =(aij) ¡Ê Mn(C), we have 
n
n
lAl1 = sup lAxl1 = max |aij|
j lxl1=1 i=1
x¡ÊCn 
n
n
lAl¡Þ = sup lAxl¡Þ = max |aij|
i
x¡ÊCn 
j=1
lxl¡Þ=1 
lAl2 = sup lAxl2 = ¦Ñ(A.A)= ¦Ñ(AA.). 
x¡ÊCn lxl2=1 
Note that lAl1 is the maximum of the f1-norms of the columns of A and lAl¡Þ is the maximum of the f1-norms of the rows of A. Furthermore, lA.l2 = lAl2, the norm ll2 is unitarily invariant, which means that 
lAl2 = lUAV l2 
for all unitary matrices U, V , and if A is a normal matrix, then lAl2 = ¦Ñ(A). 
8.3. SUBORDINATE NORMS 
Proof. For every vector u, we have 
n
nn 
= 

½Ð½Ðn½Ð½Ð n 
aijuj
½Ð½Ð½Ð½Ð
lAul

¡Ü|uj||aij|¡Ü max |aij|lul
j
1 1 , 
ij j i i 
which implies that 

n
lAl¡Ü max
1 j |aij|. i=1 
It remains to show that equality can be achieved. For this let j0 be some index such that 
n 
n
n 
max |aij| = |aij0 |, 
j ii 
and let ui = 0 for all i = j0 and uj0 = 1. In a similar way, we have 
½Ð½Ð½Ð½Ð
lAul¡Þ = max aijuj¡Ü 
i
½Ð½Ðn½Ð½Ð
max 

i 
n 

jj 
which implies that 
|aij|lul¡Þ , 
n
n
lAl¡Þ ¡Ü max |aij|. 
i 
j=1 
To achieve equality, let i0 be some index such that 
n
n 
max |aij| = |ai0j|. 
i 
jj 
The reader should check that the vector given by 
 

ai0j if ai0j =0 
|ai0j |
uj =
1 if ai0j =0 
works. We have lAl22 = sup lAxl22 = sup x . A . Ax. 
x¡ÊCn x¡ÊCn 
.. 
xx=1 xx=1 
Since the matrix A.A is symmetric, it has real eigenvalues and it can be diagonalized with respect to a unitary matrix. These facts can be used to prove that the function x ¡ú x .A.Ax has a maximum on the sphere x . x = 1 equal to the largest eigenvalue of A.A, namely, ¦Ñ(A.A). We postpone the proof until we discuss optimizing quadratic functions. Therefore, 
lAl2 = ¦Ñ(A.A). 
Let use now prove that ¦Ñ(A.A)= ¦Ñ(AA.). First assume that ¦Ñ(A.A) > 0. In this case, there is some eigenvector u (= 0) such that A . Au = ¦Ñ(A . A)u, and since ¦Ñ(A.A) > 0, we must have Au = 0. Since Au = 0, AA . (Au)= A(A . Au)= ¦Ñ(A . A)Au which means that ¦Ñ(A.A) is an eigenvalue of AA., and thus ¦Ñ(A . A) ¡Ü ¦Ñ(AA . ). Because (A.). = A, by replacing A by A. , we get ¦Ñ(AA . ) ¡Ü ¦Ñ(A . A), and so ¦Ñ(A.A)= ¦Ñ(AA.). If ¦Ñ(A.A) = 0, then we must have ¦Ñ(AA.) = 0, since otherwise by the previous reasoning we would have ¦Ñ(A.A)= ¦Ñ(AA.) > 0. Hence, in all case 
lAl22 = ¦Ñ(A . A)= ¦Ñ(AA . )= lA .l22 . For any unitary matrices U and V , it is an easy exercise to prove that V .A.AV and A.A have the same eigenvalues, so 
lAl22 = ¦Ñ(A . A)= ¦Ñ(V . A . AV )= lAV l22 , and also lAl22 = ¦Ñ(A . A)= ¦Ñ(A . U . UA)= lUAl22 . Finally, if A is a normal matrix (AA. = A.A), it can be shown that there is some unitary matrix U so that 
A = UDU . , 
where D = diag(¦Ë1,...,¦Ën) is a diagonal matrix consisting of the eigenvalues of A, and thus 
A . A =(UDU . ) . UDU . = UD . U . UDU . = UD . DU . . 

However, D.D = diag(|¦Ë1|2 ,..., |¦Ën|2), which proves that ¦Ñ(A . A)= ¦Ñ(D . D) = max|¦Ëi|2 =(¦Ñ(A))2 , 
i 
so that lAl2 = ¦Ñ(A). 
De.nition 8.9. For A =(aij) ¡Ê Mn(C), the norm lAl2 = is often called the spectral norm. 

8.3. SUBORDINATE NORMS 
Observe that Property (3) of Proposition 8.7 says that 
¡Ì 
lAl2 ¡ÜlAlF ¡Ü n lAl2 , 
which shows that the Frobenius norm is an upper bound on the spectral norm. The Frobenius norm is much easier to compute than the spectral norm. 
The reader will check that the above proof still holds if the matrix A is real (change unitary to orthogonal), con.rming the fact that lAlR = lAl for the vector norms ll1 , ll2, and ll¡Þ. It is also easy to verify that the proof goes through for rectangular m¡Án matrices, with the same formulae. Similarly, the Frobenius norm given by 
mn1/2
nn 
lAlF = |aij|2 = tr(A.A) = tr(AA.) 
i=1 j=1 
is also a norm on rectangular matrices. For these norms, whenever AB makes sense, we have 
lABl¡ÜlAllBl . 
11 
Remark: It can be shown that for any two real numbers p, q ¡Ý 1 such that + = 1, we 
pq 
have 
lA .l q = lAl p = sup{Ìá(y . Ax) |lxl p =1, lyl q =1} 
= sup{|£¨Ax, y£©| |l xl =1, lyl =1},
pq 
where lA.l q and lAl p are the operator norms. 
Remark: Let (E, ll) and (F, ll) be two normed vector spaces (for simplicity of notation, we use the same symbol ll for the norms on E and F ; this should not cause any confusion). Recall that a function f : E ¡ú F is continuous if for every a ¡Ê E, for every ¦É> 0, there is some ¦Ç> 0 such that for all x ¡Ê E, 
if lx . al¡Ü ¦Ç then lf(x) . f(a)l¡Ü ¦É. 
It is not hard to show that a linear map f : E ¡ú F is continuous i. there is some constant C ¡Ý 0 such that 
lf(x)l¡Ü C lxl for all x ¡Ê E. 
If so, we say that f is bounded (or a linear bounded operator). We let L(E; F ) denote the set of all continuous (equivalently, bounded) linear maps from E to F . Then we can de.ne the operator norm (or subordinate norm) ll on L(E; F ) as follows: for every f ¡ÊL(E; F ), 
lf(x)l 
lfl = sup = sup lf(x)l , x¡ÊE lxl x¡ÊE x=0 lxl=1 or equivalently by 
lfl = inf{¦Ë ¡Ê R |lf(x)l¡Ü ¦Ë lxl , for all x ¡Ê E}. 
It is not hard to show that the map f ¡úlfl is a norm on L(E; F ) satisfying the property lf(x)l¡Ülfllxl for all x ¡Ê E, and that if f ¡ÊL(E; F ) and g ¡ÊL(F ; G), then lg . fl¡Ülgllfl . Operator norms play an important role in functional analysis, especially when the spaces E and F are complete. 
8.4 Inequalities Involving Subordinate Norms 
In this section we discuss two technical inequalities which will be needed for certain proofs in the last three sections of this chapter. First we prove a proposition which will be needed when we deal with the condition number of a matrix. 
Proposition 8.11. Let ll be any matrix norm, and let B ¡Ê Mn(C) such that lBl < 1. 
(1) If ll is a subordinate matrix norm, then the matrix I + B is invertible and 
ll

(I + B).1 
ll

1 

¡Ü . 
1 .lBl 
(2) If a matrix of the form I + B is singular, then lBl¡Ý 1 for every matrix norm (not necessarily subordinate). Proof. (1) Observe that (I + B)u = 0 implies Bu = .u, so 
lul = lBul . Recall that 
lBul¡ÜlBllul for every subordinate norm. Since lBl < 1, if u = 0, then lBul < lul , which contradicts lul = lBul. Therefore, we must have u = 0, which proves that I + B is injective, and thus bijective, i.e., invertible. Then we have (I + B).1 + B(I + B).1 =(I + B)(I + B).1 = I, 
8.4. INEQUALITIES INVOLVING SUBORDINATE NORMS 
so we get 
(I + B).1 = I . B(I + B).1 , 
which yields

ll

(I + B).1 
ll

ll

(I + B).1 
ll

¡Ü 1+ lBl 

, 

and .nally,

ll

(I + B).1 
ll

1 

¡Ü . 
1 .lBl 
(2) If I + B is singular, then .1 is an eigenvalue of B, and by Proposition 8.6, we get ¦Ñ(B) ¡ÜlBl, which implies 1 ¡Ü ¦Ñ(B) ¡ÜlBl. 
The second inequality is a result is that is needed to deal with the convergence of se-quences of powers of matrices. 
Proposition 8.12. For every matrix A ¡Ê Mn(C) and for every ¦É> 0, there is some subor-dinate matrix norm ll such that 
lAl¡Ü ¦Ñ(A)+ ¦É. 
Proof. By Theorem 14.5, there exists some invertible matrix U and some upper triangular matrix T such that 
A = UTU.1 , 
and say that 

.
. 
¦Ë1 t12 t13 ¡¤¡¤¡¤ t1n 0 ¦Ë2 t23 ¡¤¡¤¡¤ t2n 
.. ..
.
...... 

...... 

T 

= 

.. 

.. . 

,

.

.. 

.. 

00 ¡¤¡¤¡¤ ¦Ën.1 tn.1 n 00 ¡¤¡¤¡¤ 0 ¦Ën 
where ¦Ë1,...,¦Ën are the eigenvalues of A. For every ¦Ä = 0, de.ne the diagonal matrix D¦Ä = diag(1, ¦Ä, ¦Ä2,...,¦Än.1), 
and consider the matrix 
.
. 
¦Ë1 ¦Ät12 ¦Ä2t13 ¡¤¡¤¡¤ ¦Än.1t1n 0 ¦Ë2 ¦Ät23 ¡¤¡¤¡¤ ¦Än.2t2n 
(UD¦Ä).1A(UD¦Ä)= D¦Ä .1TD¦Ä = 
...... 

...... 

. 

.. ..
.
.. 

.. .

.

.. 

.. 

00 ¡¤¡¤¡¤ ¦Ën.1 ¦Ätn.1 n 00 ¡¤¡¤¡¤ 0 ¦Ën 
Now de.ne the function ll:Mn(C) ¡ú R by 
ll

ll

lBl = 

(UD¦Ä).1B(UD¦Ä) 
,
¡Þ 
for every B ¡Ê Mn(C). Then it is easy to verify that the above function is the matrix norm subordinate to the vector norm 
v ¡ú 

ll

(UD¦Ä).1 
v 

ll

.
¡Þ 
Furthermore, for every ¦É> 0, we can pick ¦Ä so that 

n|¦Äj.itij|¡Ü ¦É, 1 ¡Ü i ¡Ü n . 1, j=i+1 
and by de.nition of the norm ll¡Þ, we get lAl¡Ü ¦Ñ(A)+ ¦É, which shows that the norm that we have constructed satis.es the required properties. Note that equality is generally not possible; consider the matrix 01 
A = ,
00 for which ¦Ñ(A)=0 < lAl, since A = 0. 
n 
8.5 Condition Numbers of Matrices 
Unfortunately, there exist linear systems Ax = b whose solutions are not stable under small perturbations of either b or A. For example, consider the system 
.
...
.
. 
1078 7 x1 32 
... 

756 5 

8 610 9 

... 
... 

x2 
x3 
... 

= 

... 

23 

33 

...

. 

7 5 9 10 x4 31 
The reader should check that it has the solution x = (1, 1, 1, 1). If we perturb slightly the 

right-hand side as b +¦¤b, where 

.
. 

¦¤b = 

... 

0.1 .0.1 0.1 .0.1 
...

, 

we obtain the new system 

...
..
. 
1078 7 x1 +¦¤x1 32.1 
... 

756 5 

8 610 9 

... 
... 

x2 +¦¤x2 
x3 +¦¤x3 
... 

= 

... 

22.9 

33.1 

...

. 

7 5 9 10 x4 +¦¤x4 30.9 
8.5. CONDITION NUMBERS OF MATRICES 
The new solution turns out to be x +¦¤x = (9.2, .12.6, 4.5, .1.1), where ¦¤x = (9.2, .12.6, 4.5, .1.1) . (1, 1, 1, 1) = (8.2, .13.6, 3.5, .2.1). Then a relative error of the data in terms of the one-norm, l¦¤bl0.44 1
1 == ¡Ö ,
lbl1 119 1190 300produces a relative error in the input l¦¤xl1 27.4 
= ¡Ö 7. 
lxl1 4 
So a relative order of the order 1/300 in the data produces a relative error of the order 7/1 in the solution, which represents an ampli.cation of the relative error of the order 2100. 
Now let us perturb the matrix slightly, obtaining the new system 
...
..
. 
10 78.17.2 x1 +¦¤x1 32 
... 

7.08 5.046 5 

85.98 9.98 9 

... 
... 

x2 +¦¤x2 
x3 +¦¤x3 
... 

= 

... 

23 

33 

...

. 

6.99 4.99 9 9.98 x4 +¦¤x4 31 
This time the solution is x +¦¤x =(.81, 137, .34, 22). Again a small change in the data alters the result rather drastically. Yet the original system is symmetric, has determinant 1, and has integer entries. The problem is that the matrix of the system is badly conditioned, a concept that we will now explain. 
Given an invertible matrix A, .rst assume that we perturb b to b+¦¤b, and let us analyze the change between the two exact solutions x and x +¦¤x of the two systems 
Ax = b A(x +¦¤x)= b +¦¤b. 
We also assume that we have some norm ll and we use the subordinate matrix norm on matrices. From 
Ax = b 
Ax + A¦¤x = b +¦¤b, 

we get ¦¤x = A.1¦¤b, 
and we conclude that 
l¦¤xl¡Ü 

ll

A.1 
ll

l¦¤bl 

lbl¡ÜlAllxl . 

Consequently, the relative error in the result l¦¤xl / lxl is bounded in terms of the relative error l¦¤bl / lbl in the data as follows: 
l¦¤xl 

¡ÜlAl 

ll

A.1 
ll

l¦¤bl 

. 

lxllbl 
Now let us assume that A is perturbed to A+¦¤A, and let us analyze the change between the exact solutions of the two systems 
Ax = b (A +¦¤A)(x +¦¤x)= b. 
The second equation yields Ax + A¦¤x +¦¤A(x +¦¤x)= b, and by subtracting the .rst equation we get 
¦¤x = .A.1¦¤A(x +¦¤x). 
It follows that 

l¦¤xl¡Ü 

ll

A.1 
ll

l¦¤Allx +¦¤xl , 

which can be rewritten as 

l¦¤xl 

¡ÜlAl 

ll

A.1 
ll

l¦¤Al 

. 

lx +¦¤xllAl 
Observe that the above reasoning is valid even if the matrix A +¦¤A is singular, as long as x +¦¤x is a solution of the second system. Furthermore, if l¦¤Al is small enough, it is not unreasonable to expect that the ratio l¦¤xl / lx +¦¤xl is close to l¦¤xl / lxl. This will be made more precise later. 
In summary, for each of the two perturbations, we see that the relative error in the result is bounded by the relative error in the data, multiplied the number lAllA.1l. In fact, this factor turns out to be optimal and this suggests the following de.nition: 
De.nition 8.10. For any subordinate matrix norm ll, for any invertible matrix A, the 
number 

cond(A)= lAl 

ll

A.1 
ll 

is called the condition number of A relative to ll. 
The condition number cond(A) measures the sensitivity of the linear system Ax = b to variations in the data b and A; a feature referred to as the condition of the system. Thus, when we says that a linear system is ill-conditioned, we mean that the condition number of its matrix is large. We can sharpen the preceding analysis as follows: 
Proposition 8.13. Let A be an invertible matrix and let x and x +¦¤x be the solutions of the linear systems 
Ax = b A(x +¦¤x)= b +¦¤b. 
8.5. CONDITION NUMBERS OF MATRICES 
If b =0, then the inequality 
l¦¤xll¦¤bl 
¡Ü cond(A)
lxllbl holds and is the best possible. This means that for a given matrix A, there exist some vectors b =0 and ¦¤b =0 for which equality holds. 
Proof. We already proved the inequality. Now, because ll is a subordinate matrix norm, there exist some vectors x =0 and ¦¤b = 0 for which 
ll

A.1¦¤b 
ll

= 

ll

A.1 
ll

l¦¤bl and lAxl = lAllxl . 

Proposition 8.14. Let A be an invertible matrix and let x and x +¦¤x be the solutions of the two systems 
Ax = b (A +¦¤A)(x +¦¤x)= b. 
If b =0, then the inequality 
l¦¤xll¦¤Al 
¡Ü cond(A)
lx +¦¤xllAl 
holds and is the best possible. This means that given a matrix A, there exist a vector b =0 and a matrix ¦¤A =0 for which equality holds. Furthermore, if l¦¤Al is small enough (for instance, if l¦¤Al < 1/ lA.1l), we have 
l¦¤xll¦¤Al 
¡Ü cond(A) (1+ O(l¦¤Al));
lxllAl 
in fact, we have 
l¦¤xll¦¤Al 1 
¡Ü cond(A) . 
lxllAl 1 .lA.1ll¦¤Al 
Proof. The .rst inequality has already been proven. To show that equality can be achieved, let w be any vector such that w = 0 and 
ll

A.1 
w 

ll

= 

ll

A.1 
ll

lwl , 

and let ¦Â = 0 be any real number. Now the vectors 
¦¤x = .¦ÂA.1 w x +¦¤x = w b =(A + ¦ÂI)w 
300  CHAPTER 8.  VECTOR NORMS AND MATRIX NORMS  
and the matrix sastisfy the equations  Ax = b  ¦¤A = ¦ÂI  

(A +¦¤A)(x +¦¤x)= b 

l¦¤xl = |¦Â| 

ll

A.1 
w 

ll

= l¦¤Al 

ll

A.1 
ll

lx +¦¤xl . 

Finally we can pick ¦Â so that .¦Â is not equal to any of the eigenvalues of A, so that A +¦¤A = A + ¦ÂI is invertible and b is is nonzero. 
If l¦¤Al < 1/ lA.1l, then 
ll

A.1¦¤A 
ll

¡Ü 

ll

A.1 
ll

l¦¤Al < 1, 

so by Proposition 8.11, the matrix I + A.1¦¤A is invertible and 
ll

(I + A.1¦¤A).1 
ll

11 

¡Ü¡Ü . 
1 .lA.1¦¤Al 1 .lA.1ll¦¤Al 
Recall that we proved earlier that 
¦¤x = .A.1¦¤A(x +¦¤x), 
and by adding x to both sides and moving the right-hand side to the left-hand side yields 
(I + A.1¦¤A)(x +¦¤x)= x, 

and thus 
x +¦¤x =(I + A.1¦¤A).1 x, which yields ¦¤x = ((I + A.1¦¤A).1 . I)x =(I + A.1¦¤A).1(I . (I + A.1¦¤A))x = .(I + A.1¦¤A).1A.1(¦¤A)x. 
From this and

ll

(I + A.1¦¤A).1 
ll

1 

¡Ü ,
1 .lA.1ll¦¤Al
we get 
lA.1ll¦¤Al 
l¦¤xl¡Ü lxl ,
1 .lA.1ll¦¤Al which can be written as l¦¤xll¦¤Al 1 
¡Ü cond(A) ,
lxllAl 1 .lA.1ll¦¤Al 
which is the kind of inequality that we were seeking. 

8.5. CONDITION NUMBERS OF MATRICES 
Remark: If A and b are perturbed simultaneously, so that we get the ¡°perturbed¡± system 

(A +¦¤A)(x +¦¤x)= b +¦¤b, 
it can be shown that if l¦¤Al < 1/ lA.1l (and b = 0), then 
l¦¤xl cond(A) l¦¤All¦¤bl 
¡Ü +;
lxl 1 .lA.1ll¦¤AllAllbl 
see Demmel [16], Section 2.2 and Horn and Johnson [36], Section 5.8. 
We now list some properties of condition numbers and .gure out what cond(A) is in the case of the spectral norm (the matrix norm induced by ll2). First, we need to introduce a very important factorization of matrices, the singular value decomposition, for short, SVD. 
It can be shown (see Section 20.2) that given any n ¡Á n matrix A ¡Ê Mn(C), there exist two unitary matrices U and V , and a real diagonal matrix ¦² = diag(¦Ò1,...,¦Òn), with ¦Ò1 ¡Ý ¦Ò2 ¡Ý¡¤ ¡¤¡¤¡Ý ¦Òn ¡Ý 0, such that 
A = V ¦²U . . 
De.nition 8.11. Given a complex n ¡Á n matrix A, a triple (U, V, ¦²) such that A = V ¦²UT , where U and V are n ¡Á n unitary matrices and ¦² = diag(¦Ò1,...,¦Òn) is a diagonal matrix of real numbers ¦Ò1 ¡Ý ¦Ò2 ¡Ý ¡¤ ¡¤¡¤ ¡Ý ¦Òn ¡Ý 0, is called a singular decomposition (for short SVD) of 
A. If A is a real matrix, then U and V are orthogonal matrices The nonnegative numbers ¦Ò1,...,¦Òn are called the singular values of A. 
The factorization A = V ¦²U. implies that 
A . A = U¦²2U . and AA . = V ¦²2V . , 
which shows that ¦Ò12,...,¦Òn 2 are the eigenvalues of both A.A and AA., that the columns of U are corresponding eivenvectors for A.A, and that the columns of V are corresponding eivenvectors for AA. . 
Since ¦Ò12 is the largest eigenvalue of A.A (and AA.), note that ¦Ñ(A.A)= ¦Ñ(AA.)= ¦Ò1. 
Corollary 8.15. The spectral norm lAl2 of a matrix A is equal to the largest singular value of A. Equivalently, the spectral norm lAl2 of a matrix A is equal to the f¡Þ-norm of its vector of singular values, 
lAl2 = max ¦Òi = l(¦Ò1,...,¦Òn)l¡Þ . 
1¡Üi¡Ün 
Since the Frobenius norm of a matrix A is de.ned by lAlF = tr(A.A) and since 

tr(A . A)= ¦Ò12 + ¡¤¡¤¡¤ + ¦Òn 2 
where ¦Ò12,...,¦Òn 2 are the eigenvalues of A.A, we see that 
lAlF =(¦Ò12 + ¡¤¡¤¡¤ + ¦Òn2)1/2 = l(¦Ò1,...,¦Òn)l2 . 
Corollary 8.16. The Frobenius norm of a matrix is given by the f2-norm of its vector of singular values; lAlF = l(¦Ò1,...,¦Òn)l2. In the case of a normal matrix if ¦Ë1,...,¦Ën are the (complex) eigenvalues of A, then ¦Òi = |¦Ëi|, 1 ¡Ü i ¡Ü n. Proposition 8.17. For every invertible matrix A ¡Ê Mn(C), the following properties hold: 
(1) 


cond(A) ¡Ý 1, cond(A) = cond(A.1) cond(¦ÁA) = cond(A) for all ¦Á ¡Ê C.{0}. 
(2) 
If cond2(A) denotes the condition number of A with respect to the spectral norm, then ¦Ò1

cond2(A)= ,¦Òn 
where ¦Ò1 ¡Ý ¡¤ ¡¤¡¤ ¡Ý ¦Òn are the singular values of A. 
(3) If the matrix A is normal, then 
|¦Ë1|

cond2(A)= ,
|¦Ën|
where ¦Ë1,...,¦Ën are the eigenvalues of A sorted so that |¦Ë1|¡Ý¡¤ ¡¤¡¤¡Ý |¦Ën|. 

(4) 
If A is a unitary or an orthogonal matrix, then 
cond2(A)=1. 


(5) 
The 	condition number cond2(A) is invariant under unitary transformations, which means that 


cond2(A) = cond2(UA) = cond2(AV ), for all unitary matrices U and V . Proof. The properties in (1) are immediate consequences of the properties of subordinate matrix norms. In particular, AA.1 = I implies 
1= lIl¡ÜlAl 

ll

A.1 
ll

= cond(A). 

(2) We showed earlier that lAl22 = ¦Ñ(A.A), which is the square of the modulus of the largest eigenvalue of A.A. Since we just saw that the eigenvalues of A.A are ¦Ò12 ¡Ý ¡¤ ¡¤¡¤ ¡Ý ¦Òn2, where ¦Ò1,...,¦Òn are the singular values of A, we have 
lAl2 = ¦Ò1. 
8.5. CONDITION NUMBERS OF MATRICES 
Now if A is invertible, then ¦Ò1 ¡Ý¡¤ ¡¤¡¤ ¡Ý ¦Òn > 0, and it is easy to show that the eigenvalues of (A.A).1 are ¦Òn .2 ¡Ý¡¤ ¡¤¡¤ ¡Ý ¦Ò1 .2, which shows that 
ll

A.1 
ll 

2 = ¦Òn .1 , 
and thus 
¦Ò1
cond2(A)= . ¦Òn 
(3) This follows from the fact that lAl2 = ¦Ñ(A) for a normal matrix. 
(4) 
If A is a unitary matrix, then A.A = AA. = I, so ¦Ñ(A.A) = 1, and lAl2 = ¦Ñ(A.A) = 1. We also have lA.1l2 = lA.l2 = ¦Ñ(AA.) = 1, and thus cond(A) = 1. 

(5) 
This follows immediately from the unitary invariance of the spectral norm. 


Proposition 8.17 (4) shows that unitary and orthogonal transformations are very well-conditioned, and Part (5) shows that unitary transformations preserve the condition number. 
In order to compute cond2(A), we need to compute the top and bottom singular values of A, which may be hard. The inequality 
¡Ì 
lAl2 ¡ÜlAlF ¡Ü n lAl2 , 
may be useful in getting an approximation of cond2(A)= lAl2 lA.1l2, if A.1 can be determined. 
Remark: There is an interesting geometric characterization of cond2(A). If ¦È(A) denotes the least angle between the vectors Au and Av as u and v range over all pairs of orthonormal vectors, then it can be shown that 
cond2(A) = cot(¦È(A)/2)). 
Thus if A is nearly singular, then there will be some orthonormal pair u, v such that Au and Av are nearly parallel; the angle ¦È(A) will the be small and cot(¦È(A)/2)) will be large. For more details, see Horn and Johnson [36] (Section 5.8 and Section 7.4). 
It should be noted that in general (if A is not a normal matrix) a matrix could have a very large condition number even if all its eigenvalues are identical! For example, if we consider the n ¡Á n matrix 
.
. 
120 0 ... 00 012 0 ... 00 001 2 ... 00 
.. ..
...
A = 

.......... 

.......... 

,

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

it turns out that cond2(A) ¡Ý 2n.1 . 
A classical example of matrix with a very large condition number is the Hilbert matrix H(n), the n ¡Á n matrix with 
1 
Hij = i + j . 1 . 
(n) 
For example, when n = 5, 

.
. 
1

1 
5 
1 
4 
1 
3 
1 
2 
....... 
1 
61 
71 
81 
9 
1 
51 
61 
71 
8 
1 
41 
51 
61 
7 
1 
31 
41 
51 
6 
1 
21 
31 
41 
5 
....... 
H(5) = 
. 

It can be shown that cond2(H(5)) ¡Ö 4.77 ¡Á 105 . 
Hilbert introduced these matrices in 1894 while studying a problem in approximation theory. The Hilbert matrix H(n) is symmetric positive de.nite. A closed-form formula can be given for its determinant (it is a special form of the so-called Cauchy determinant); see Problem 8.15. The inverse of H(n) can also be computed explicitly; see Problem 8.15. It can be shown that 
¡Ì¡Ì 
cond2(H(n))= O((1+ 2)4n/n). 
Going back to our matrix 

.
. 
... 

1078 7 
756 5 

8 610 9 
7 5 9 10 

...

A = 

, 

which is a symmetric positive de.nite matrix, it can be shown that its eigenvalues, which in this case are also its singular values because A is SPD, are ¦Ë1 ¡Ö 30.2887 >¦Ë2 ¡Ö 3.858 >¦Ë3 ¡Ö 0.8431 >¦Ë4 ¡Ö 0.01015, so that 
¦Ë1
cond2(A)= ¡Ö 2984. 
¦Ë4 The reader should check that for the perturbation of the right-hand side b used earlier, the relative errors l¦¤xl /lxl and l¦¤xl /lxl satisfy the inequality 
l¦¤xll¦¤bl 
¡Ü cond(A)
lxllbl 
and comes close to equality. 
8.6. AN APPLICATION OF NORMS: INCONSISTENT LINEAR SYSTEMS 
8.6 	An Application of Norms: Solving Inconsistent Linear Systems 
The problem of solving an inconsistent linear system Ax = b often arises in practice. This is a system where b does not belong to the column space of A, usually with more equations than variables. Thus, such a system has no solution. Yet we would still like to ¡°solve¡± such a system, at least approximately. 
Such systems often arise when trying to .t some data. For example, we may have a set of 3D data points 
{p1,...,pn}, 
and we have reason to believe that these points are nearly coplanar. We would like to .nd a plane that best .ts our data points. Recall that the equation of a plane is 
¦Áx + ¦Ây + ¦Ãz + ¦Ä =0, 
with (¦Á, ¦Â, ¦Ã) = (0, 0, 0). Thus, every plane is either not parallel to the x-axis (¦Á = 0) or not parallel to the y-axis (¦Â = 0) or not parallel to the z-axis (¦Ã = 0). 
Say we have reasons to believe that the plane we are looking for is not parallel to the z-axis. If we are wrong, in the least squares solution, one of the coe.cients, ¦Á, ¦Â, will be very large. If ¦Ã = 0, then we may assume that our plane is given by an equation of the form 
z = ax + by + d, 
and we would like this equation to be satis.ed for all the pi¡¯s, which leads to a system of n equations in 3 unknowns a, b, d, with pi =(xi,yi,zi); 
ax1 + by1 + d = z1 .. 
.. 
.. axn + byn + d = zn. 
However, if n is larger than 3, such a system generally has no solution. Since the above system can¡¯t be solved exactly, we can try to .nd a solution (a, b, d) that minimizes the least-squares error 
n
n 
(axi + byi + d . zi)2 . i=1 
This is what Legendre and Gauss .gured out in the early 1800¡¯s! In general, given a linear system Ax = b, 
we solve the least squares problem: minimize lAx . bl22. Fortunately, every n ¡Á m-matrix A can be written as 
A = V DUT 
where U and V are orthogonal and D is a rectangular diagonal matrix with non-negative entries (singular value decomposition, or SVD); see Chapter 20. 
The SVD can be used to solve an inconsistent system. It is shown in Chapter 21 that there is a vector x of smallest norm minimizing lAx . bl2. It is given by the (Penrose) pseudo-inverse of A (itself given by the SVD). 
It has been observed that solving in the least-squares sense may give too much weight to ¡°outliers,¡± that is, points clearly outside the best-.t plane. In this case, it is preferable to minimize (the f1-norm) 
n
n 
|axi + byi + d . zi|. i=1 
This does not appear to be a linear problem, but we can use a trick to convert this minimization problem into a linear program (which means a problem involving linear con-straints). 
Note that |x| = max{x, .x}. So by introducing new variables e1,...,en, our minimiza-tion problem is equivalent to the linear program (LP): 
minimize e1 + ¡¤¡¤¡¤ + en 
subject to axi + byi + d . zi ¡Ü ei 
.(axi + byi + d . zi) ¡Ü ei 
1 ¡Ü i ¡Ü n. 
Observe that the constraints are equivalent to 
ei ¡Ý|axi + byi + d . zi|, 1 ¡Ü i ¡Ü n. 
For an optimal solution, we must have equality, since otherwise we could decrease some ei and get an even better solution. Of course, we are no longer dealing with ¡°pure¡± linear algebra, since our constraints are inequalities. 
We prefer not getting into linear programming right now, but the above example provides a good reason to learn more about linear programming! 
8.7 Limits of Sequences and Series 
 
n
If x ¡Ê R or x ¡Ê Cand if |x| < 1, it is well known that the sumsn xk = 1+x+x2 +¡¤¡¤¡¤+x
k=0 
converge to the limit 1/(1 . x) when n goes to in.nity, and we write 
¡Þ
n 
1
k 
x = . 
1 . x 
k=0 
For example, 
¡Þ
n 
1 
=2. 
2k 
k=0 
8.7. LIMITS OF SEQUENCES AND SERIES 
Similarly, the sums 
nk
n 
x
Sn = k! 
k=0 
x

converge to ewhen n goes to in.nity, for every x (in R or C). What if we replace x by a real of complex n ¡Á n matrix A? 
n nAk 
The partial sums Ak and still make sense, but we have to de.ne what is 
k=0 k=0 k! 
the limit of a sequence of matrices. This can be done in any normed vector space. 
De.nition 8.12. Let (E, ll) be a normed vector space. A sequence (un)n¡ÊN in E is any function u: N ¡ú E. For any v ¡Ê E, the sequence (un) converges to v (and v is the limit of the sequence (un)) if for every ¦É> 0, there is some integer N> 0 such that 
lun . vl <¦É for all n ¡Ý N. 
Often we assume that a sequence is indexed by N.{0}, that is, its .rst term is u1 rather than u0. 
If the sequence (un) converges to v, then since by the triangle inequality 
lum . unl¡Ülum . vl + lv . unl , 
we see that for every ¦É> 0, we can .nd N> 0 such that lum . vl < ¦É/2 and lun . vl < ¦É/2, and so 
lum . unl <¦É for all m, n ¡Ý N. 
The above property is necessary for a convergent sequence, but not necessarily su.cient. For example, if E = Q, there are sequences of rationals satisfying the above condition, but whose limit is not a rational number. For example, the sequence n 1 converges to e, and 
k=1 k! 
the sequence n (.1)k 1 converges to ¦Ð/4, but e and ¦Ð/4 are not rational (in fact, they 
k=02k+1 
are transcendental). However, R is constructed from Q to guarantee that sequences with the above property converge, and so is C. 
De.nition 8.13. Given a normed vector space (E, ll), a sequence (un) is a Cauchy sequence if for every ¦É> 0, there is some N> 0 such that 
lum . unl <¦É for all m, n ¡Ý N. 
If every Cauchy sequence converges, then we say that E is complete. A complete normed vector spaces is also called a Banach space. 
A fundamental property of R is that it is complete. It follows immediately that C is also complete. If E is a .nite-dimensional real or complex vector space, since any two norms are equivalent, we can pick the f¡Þ norm, and then by picking a basis in E, a sequence (un) of vectors in E converges i. the n sequences of coordinates (uni ) (1 ¡Ü i ¡Ü n) converge, so any .nite-dimensional real or complex vector space is a Banach space. 
Let us now consider the convergence of series. 
¡Þ 
uk
k=0
De.nition 8.14. Given a normed vector space (E, ll), a series is an in.nite sum

of elements uk ¡Ê E. We denote by Sn the partial sum of the .rst n + 1 elements, 
n
Sn = uk. 
k=0 
n 
¡Þ 
k=0 
uk converges to the limit v ¡Ê E if the sequence 
De.nition 8.15. We say that the series

(Sn) converges to v, i.e., given any ¦É> 0, there exists a positive integer N such that for all n ¡Ý N, 
lSn . vl < ¦É. 
¡Þ 
uk converges 
k=0
In this case, we say that the series is convergent. We say that the series

¡Þ 
k=0 
lukl is convergent. 
absolutely if the series of norms

If the series

¡Þ 
k=0 
uk converges to v, since for all m, n with m>n we have 
n
mmm
uk . Sn = uk . uk = uk, 
k=0 k=0 k=0 k=n+1 
n
n
nn 
if we let m go to in.nity (with n .xed), we see that the series

¡Þ 
k=n+1 
uk converges and that 
¡Þn 
v . Sn = uk. 
k=n+1 
There are series that are convergent but not absolutely convergent; for example, the series 
¡Þn
(.1)k.1 
k 
k=1 
1
¡Þ
converges to ln 2, but

does not converge (this sum is in.nite). 

k=1 k 
If E is complete, the converse is an enormously useful result. 

¡Þ 
uk
k=0
Proposition 8.18. Assume (E, ll) is a complete normed vector space. If a series

is absolutely convergent, then it is convergent. 

Proof. If

¡Þ 
k=0 
uk is absolutely convergent, then we prove that the sequence (Sm) is a Cauchy 
sequence; that is, for every ¦É> 0, there is some p> 0 such that for all n ¡Ý m ¡Ý p, 
lSn . Sml¡Ü ¦É. Observe that 
lSn . Sml = lum+1 + ¡¤¡¤¡¤ + unl¡Ülum+1l + ¡¤¡¤¡¤ + lunl , 
and since the sequence

¡Þ 
k=0 
lukl converges, it satis.es Cauchy¡¯s criterion. 
Thus, the se-
quence (Sm) also satis.es Cauchy¡¯s criterion, and since E is a complete vector space, the sequence (Sm) converges. 
8.8. THE MATRIX EXPONENTIAL 
Remark: It can be shown that if (E, ll) is a normed vector space such that every absolutely convergent series is also convergent, then E must be complete (see Schwartz [54]). 
¡Þ
An important corollary of absolute convergence is that if the terms in series

uk
k=0 
are rearranged, then the resulting series is still absolutely convergent and has the same 
sum. More precisely, let ¦Ò be any permutation (bijection) of the natural numbers. The 

series

¡Þ 
k=0 
u¦Ò(k) is called a rearrangement of the original series. The following result can be 
shown (see Schwartz [54]). 

Proposition 8.19. Assume (E, ll) is a normed vector space. If a series

¡Þ 
k=0 
uk is conver-
¡Þ 
k=0 u¦Ò(k)
gent as well as absolutely convergent, then for every permutation ¦Ò of N, the series

is convergent and absolutely convergent, and its sum is equal to the sum of the original series: 

¡Þn
¡Þn 
u¦Ò(k) = uk. k=0 k=0 
In particular, if (E, ll) is a complete normed vector space, then Proposition 8.19 holds. We now apply Proposition 8.18 to the matrix exponential. 
8.8 The Matrix Exponential 
Proposition 8.20. For any n ¡Á n real or complex matrix A, the series 
¡Þn
Ak 
k! 

k=0 
converges absolutely for any operator norm on Mn(C) (or Mn(R)). 
Proof. Pick any norm on Cn (or Rn) and let ll be the corresponding operator norm on 

lll
lll

2, it is complete. By Proposition 8.18, it su.ces to 
Ak 
Mn(C). Since Mn(C) has dimension n
n
show that the series of nonnegative reals 

Since ll is an operator 

converges. 

k=0 
k! 
n 
norm, this a matrix norm, so we have 
n
llll 

llll

lAlk
n
nk=0 k=0 
Ak 
¡Ü e lAl 
.

¡Ü 

k! 

k! 

lll

Ak 
k! 
lll

is bounded by e

lAl 
,

Thus, the nondecreasing sequence of positive real numbers n 
k=0 
and by a fundamental property of R, it has a least upper bound which is its limit. 
De.nition 8.16. Let E be a .nite-dimensional real of complex normed vector space. For any n ¡Á n matrix A, the limit of the series 
¡Þn
Ak 
k! 

k=0 
is the exponential of A and is denoted eA . 
x
A basic property of the exponential x ¡ú ewith x ¡Ê C is 
x+y xy
e = ee, for all x, y ¡Ê C. 
xx).1 .x
As a consequence, eis always invertible and (e= e. For matrices, because matrix multiplication is not commutative, in general, 
A+B AB 
e = ee 
fails! This result is salvaged as follows. 
Proposition 8.21. For any two n ¡Á n complex matrices A and B, if A and B commute, that is, AB = BA, then 
A+B AB 
e = ee. 
A proof of Proposition 8.21 can be found in Gallier [25]. 
Since A and .A commute, as a corollary of Proposition 8.21, we see that eA is always invertible and that 
A).1 .A
(e = e. 
It is also easy to see that (e A)T = e AT . 
In general, there is no closed-form formula for the exponential eA of a matrix A, but for skew symmetric matrices of dimension 2 and 3, there are explicit formulae. Everyone should enjoy computing the exponential eA where 
0  .¦È  
A =  .  
¦È  0  

If we write 
0 .1 
J = ,
10 
then 
A = ¦ÈJ 
The key property is that 
J2 = .I. 
Proposition 8.22. If A = ¦ÈJ, then 
A cos ¦È . sin ¦È 
e = cos ¦ÈI + sin ¦ÈJ = . 
sin ¦È cos ¦È 
8.8. THE MATRIX EXPONENTIAL 
Proof. We have 
A4n = ¦È4nI2, A4n+1 = ¦È4n+1J, A4n+2 = .¦È4n+2I2, A4n+3 = .¦È4n+3J, 
and so 
A ¦È¦È2 ¦È3 ¦È4 ¦È5 ¦È6 ¦È7 
e = I2 + J . I2 . J + I2 + J . I2 . J + ¡¤¡¤¡¤ . 
1!2! 3!4! 5!6! 7! Rearranging the order of the terms, we have 
¦È2 ¦È4 ¦È6 ¦È¦È3 ¦È5 ¦È7 
e A =1 . + . + ¡¤¡¤¡¤ I2 + . + . + ¡¤¡¤¡¤ J. 
2!4!6! 1!3!5!7! We recognize the power series for cos ¦È and sin ¦È, and thus 
e A = cos ¦ÈI2 + sin ¦ÈJ, 
that is 
A cos ¦È . sin ¦È 
e = ,
sin ¦È cos ¦È 
as claimed. 
Thus, we see that the exponential of a 2 ¡Á 2 skew-symmetric matrix is a rotation matrix. This property generalizes to any dimension. An explicit formula when n = 3 (the Rodrigues¡¯ formula) is given in Section 11.7. 
Proposition 8.23. If B is an n ¡Á n (real) skew symmetric matrix, that is, BT = .B, then Q = eB is an orthogonal matrix, that is 
QTQ = QQT = I. 
Proof. Since BT = .B, we have B)T 
BT .B
QT =(e = e = e. 
Since B and .B commute, we have 
.BB .B+B 0
QTQ = ee = e = e = I. 
Similarly, 
B .BB.B 0
QQT = ee = e = e = I, 
which concludes the proof. 

It can also be shown that det(Q) = det(eB) = 1, but this requires a better understanding of the eigenvalues of eB (see Section 14.5). Furthermore, for every n ¡Á n rotation matrix Q (an orthogonal matrix Q such that det(Q) = 1), there is a skew symmetric matrix B such that Q = eB . This is a fundamental property which has applications in robotics for n = 3. 
All familiar series have matrix analogs. For example, if lAl < 1 (where ll is an operator 
¡Þ
norm), then the series k=0 Ak converges absolutely, and it can be shown that its limit is (I . A).1 . 
Another interesting series is the logarithm. For any n ¡Á n complex matrix A, if lAl < 1 (where ll is an operator norm), then the series 
¡Þ
n Ak log(I + A)= (.1)k+1 
k 
k=1 
converges absolutely. 
8.9 Summary 
The main concepts and results of this chapter are listed below: 
. 
Norms and normed vector spaces. 

. 
The triangle inequality. 

. 
The Euclidean norm; the fp-norms. 

. 
H¡§older¡¯s inequality; the Cauchy¨CSchwarz inequality; Minkowski¡¯s inequality. 

. 
Hermitian inner product and Euclidean inner product. 

. 
Equivalent norms. 

. 
All norms on a .nite-dimensional vector space are equivalent (Theorem 8.5). 

. 
Matrix norms. 

. 
Hermitian, symmetric and normal matrices. Orthogonal and unitary matrices. 

. 
The trace of a matrix. 

. 
Eigenvalues and eigenvectors of a matrix. 

. 
The characteristic polynomial of a matrix. 

. 
The spectral radius ¦Ñ(A) of a matrix A. 

. The Frobenius norm. 
8.10. PROBLEMS 

. 	
The Frobenius norm is a unitarily invariant matrix norm. 

. 	
Bounded linear maps. 

. 	
Subordinate matrix norms. 

. 	
Characterization of the subordinate matrix norms for the vector norms ll1 , ll2, and ll¡Þ. 

. 	
The spectral norm. 

. 	
For every matrix A ¡Ê Mn(C) and for every ¦É> 0, there is some subordinate matrix norm ll such that lAl¡Ü ¦Ñ(A)+ ¦É. 

. 	
Condition numbers of matrices. 

. 	
Perturbation analysis of linear systems. 

. 	
The singular value decomposition (SVD). 

. 	
Properties of conditions numbers. Characterization of cond2(A) in terms of the largest and smallest singular values of A. 

. 	
The Hilbert matrix: a very badly conditioned matrix. 

. 	
Solving inconsistent linear systems by the method of least-squares; linear programming. 

. 	
Convergence of sequences of vectors in a normed vector space. 

. 	
Cauchy sequences, complex normed vector spaces, Banach spaces. 

. 	
Convergence of series. Absolute convergence. 

. 	
The matrix exponential. 

. 	
Skew symmetric matrices and orthogonal matrices. 


8.10 Problems 
Problem 8.1. Let A be the following matrix: 
¡Ì 
11/ 2 
B = ¡Ì . 
1/ 23/2 
Compute the operator 2-norm lAl2 of A. 
Problem 8.2. Prove Proposition 8.3, namely that the following inequalities hold for all x ¡Ê Rn (or x ¡Ê Cn): 
lxl¡Þ ¡Ülxl1 ¡Ü nlxl¡Þ, 
¡Ì 
lxl¡Þ ¡Ülxl2 ¡Ü nlxl¡Þ, 
¡Ì 
lxl2 ¡Ülxl1 ¡Ü nlxl2. Problem 8.3. For any p ¡Ý 1, prove that for all x ¡Ê Rn , lim lxl p = lxl¡Þ . 
p¡ú¡Þ 
Problem 8.4. Let A be an n ¡Á n matrix which is strictly row diagonally dominant, which 
means that 

n 

|aii| > |aij|, 
j=i 
for i =1,...,n, and let 

 n 
 
¦Ä = min |aii|. |aij|. 
ij=i 
The fact that A is strictly row diagonally dominant is equivalent to the condition ¦Ä> 0. 
(1) For any nonzero vector v, prove that 
lAvl¡Þ ¡Ýlvl¡Þ ¦Ä. Use the above to prove that A is invertible. 
(2) Prove that

ll

A.1 
ll

.
¡Þ ¡Ü ¦Ä.1 
Hint. Prove that 
lA.1vl¡Þ lwl¡Þ
sup = sup . v=0 lvl¡Þ w=0 lAwl¡Þ 
Problem 8.5. Let A be any invertible complex n ¡Á n matrix. 
(1) For any vector norm ll on Cn, prove that the function llA : Cn ¡ú R given by 
lxlA = lAxl for all x ¡Ê Cn , is a vector norm. 
(2) Prove that the operator norm induced by llA, also denoted by llA, is given by 
lBl

= 
A 
ll

ABA.1 
ll

for every n ¡Á n matrix B, 

where lABA.1l uses the operator norm induced by ll. 
8.10. PROBLEMS 
Problem 8.6. Give an example of a norm on Cn and of a real matrix A such that 
lAlR < lAl , 
where l.lR and l.l are the operator norms associated with the vector norm l.l. Hint. This can already be done for n = 2. 
Problem 8.7. Let ll be any operator norm. Given an invertible n ¡Á n matrix A, if c =1/(2 lA.1l), then for every n ¡Á n matrix H, if lHl¡Ü c, then A + H is invertible. Furthermore, show that if lHl¡Ü c, then l(A + H).1l¡Ü 1/c. 
Problem 8.8. Let A be any m ¡Á n matrix and let ¦Ë ¡Ê R be any positive real number ¦Ë> 0. 
(1) 
Prove that ATA + ¦ËIn and AAT + ¦ËIm are invertible. 

(2) 
Prove that 
AT(AAT + ¦ËIm).1 =(ATA + ¦ËIn).1AT . 



Remark: The expressions above correspond to the matrix for which the function 
¦µ(x)=(Ax . b)T(Ax . b)+ ¦ËxT x 
achieves a minimum. It shows up in machine learning (kernel methods). 
Problem 8.9. Let Z be a q ¡Á p real matrix. Prove that if Ip . ZTZ is positive de.nite, then the (p + q) ¡Á (p + q) matrix 
ZTIp
S = 
ZIq 
is symmetric positive de.nite. 
Problem 8.10. Prove that for any real or complex square matrix A, we have 
lAl22 ¡ÜlAl1 lAl¡Þ , 
where the above norms are operator norms. 

Hint. Use Proposition 8.10 (among other things, it shows that lAl

= 
1 
ll

AT 
ll 

¡Þ
). 

Problem 8.11. Show that the map A ¡ú ¦Ñ(A) (where ¦Ñ(A) is the spectral radius of A) is neither a norm nor a matrix norm. In particular, .nd two 2 ¡Á 2 matrices A and B such that 
¦Ñ(A + B) >¦Ñ(A)+ ¦Ñ(B)=0 and ¦Ñ(AB) >¦Ñ(A)¦Ñ(B)=0. 
Problem 8.12. De.ne the map A ¡ú M(A) (de.ned on n¡Án real or complex n¡Án matrices) by 
M(A) = max{|aij|| 1 ¡Ü i, j ¡Ü n}. 
(1) Prove that 
M(AB) ¡Ü nM(A)M(B) 

for all n ¡Á n matrices A and B. 
(2) 
Give a counter-example of the inequality 

M(AB) ¡Ü M(A)M(B). 

(3) 
Prove that the map A ¡úlAlM given by 


lAl= nM(A)= n max{|aij|| 1 ¡Ü i, j ¡Ü n}
M 
is a matrix norm. 
Problem 8.13. Let S be a real symmetric positive de.nite matrix. 
(1) Use the Cholesky factorization to prove that there is some upper-triangular matrix C, unique if its diagonal elements are strictly positive, such that S = CTC. 
(2) For any x ¡Ê Rn, de.ne 
lxlS =(x TSx)1/2 . 

Prove that 
lxlS = lCxl2 , 
and that the map x ¡úlxlS is a norm. 
Problem 8.14. Let A be a real 2 ¡Á 2 matrix 
a11 a12 
A = . a21 a22 
(1) Prove that the squares of the singular values ¦Ò1 ¡Ý ¦Ò2 of A are the roots of the 
quadratic equation X2 . tr(ATA)X + | det(A)|2 =0. 
(2) If we let 
2222
a11 + a12 + a21 + a22 
¦Ì(A)= ,2|a11a22 . a12a21| 
prove that 
¦Ò1
cond2(A)= = ¦Ì(A)+(¦Ì(A)2 . 1)1/2 . 
¦Ò2 
(3) Consider the subset S of 2 ¡Á 2 invertible matrices whose entries aij are integers such that 0 ¡Ü aij ¡Ü 100. Prove that the functions cond2(A) and ¦Ì(A) reach a maximum on the set S for the same values of A. 
Check that for the matrix 
100 99 

Am = 
99 98 
we have ¦Ì(Am) = 19, 603 det(Am)= .1 
8.10. PROBLEMS 
and 
cond2(Am) ¡Ö 39, 206. 
(4) Prove that for all A ¡ÊS, if | det(A)|¡Ý 2 then ¦Ì(A) ¡Ü 10, 000. Conclude that the maximum of ¦Ì(A) on S is achieved for matrices such that det(A)= ¡À1. Prove that .nding matrices that maximize ¦Ì on S is equivalent to .nding some integers n1,n2,n3,n4 such that 
0 ¡Ü n4 ¡Ü n3 ¡Ü n2 ¡Ü n1 ¡Ü 100 
2222 

n1 + n2 + n3 + n4 ¡Ý 1002 + 992 + 992 + 982 = 39, 206 |n1n4 . n2n3| =1. 
You may use without proof that the fact that the only solution to the above constraints is the multiset {100, 99, 99, 98}. 
(5) Deduce from part (4) that the matrices in S for which ¦Ì has a maximum value are 10099 98 99 99100 99 98 
Am = 
99 98 99100 98 99 10099 and check that ¦Ì has the same value for these matrices. Conclude that max cond2(A) = cond2(Am). 
A¡ÊS 
(6) Solve the system 
100 99 x1 199 

= . 
99 98 x2 197 
Perturb the right-hand side b by 
.0.0097 

¦¤b = 
0.0106 and solve the new system Amy = b +¦¤b where y =(y1,y2). Check that 2 
¦¤x = y . x = . 
.2.0203 
Compute lxl2, l¦¤xl2, l¦¤bl, and estimate 

2, lbl2l¦¤xl2 l¦¤bl2 .1 
c = . 
lbl
lxl22 
Check that 
c ¡Ö cond2(Am) = 39, 206. 
Problem 8.15. Consider a real 2 ¡Á 2 matrix with zero trace of the form 
ab 
A = . 
c .a 
(1) Prove that 
A2 =(a 2 + bc)I2 = . det(A)I2. 

If a2 + bc = 0, prove that e A = I2 + A. 
(2) If a2 + bc < 0, let ¦Ø> 0 be such that ¦Ø2 = .(a2 + bc). Prove that 
sin ¦Ø 

e A = cos ¦ØI2 + A. 
¦Ø 
(3) If a2 + bc > 0, let ¦Ø> 0 be such that ¦Ø2 = a2 + bc. Prove that 
sinh ¦Ø 

e A = cosh ¦ØI2 + A. 
¦Ø 
(3) 
Prove that in all cases 


det e A = 1 and tr(A) ¡Ý.2. 
(4) 
Prove that there exist some real 2 ¡Á 2 matrix B with det(B) = 1 such that there is 

A
no real 2 ¡Á 2 matrix A with zero trace such that e= B. Problem 8.16. Recall that the Hilbert matrix is given by 
1
(n)
H= 
ij . 
i + j . 1 
(1) Prove that 
(1!2! ¡¤¡¤¡¤ (n . 1)!)4 

det(H(n))= ,
1!2! ¡¤¡¤¡¤ (2n . 1)! thus the reciprocal of an integer. Hint. Use Problem 6.13. 
(2) Amazingly, the entries of the inverse of H(n) are integers. Prove that (H(n)).1 =(¦Áij), with 
2 
n + i . 1 n + j . 1 i + j . 2 
¦Áij =(.1)i+j(i + j . 1) . 
n . jn . ii . 1 




