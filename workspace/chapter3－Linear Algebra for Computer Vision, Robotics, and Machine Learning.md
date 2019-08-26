Chapter 3 
Matrices and Linear Maps 
In this chapter, all vector spaces are de.ned over an arbitrary .eld K. For the sake of concreteness, the reader may safely assume that K = R. 
3.1 Representation of Linear Maps by Matrices 
Proposition 2.17 shows that given two vector spaces E and F and a basis (uj)j¡ÊJ of E, every linear map f : E ¡ú F is uniquely determined by the family (f(uj))j¡ÊJ of the images under f of the vectors in the basis (uj)j¡ÊJ . 
If we also have a basis (vi)i¡ÊI of F , then every vector f(uj) can be written in a unique way as 
 
f(uj)=aijvi, 
i¡ÊI 
where j ¡Ê J, for a family of scalars (aij)i¡ÊI . Thus, with respect to the two bases (uj)j¡ÊJ of E and (vi)i¡ÊI of F , the linear map f is completely determined by a ¡°I ¡Á J-matrix¡± M(f)=(aij)i¡ÊI, j¡ÊJ . 
Remark: Note that we intentionally assigned the index set J to the basis (uj)j¡ÊJ of E, and the index set I to the basis (vi)i¡ÊI of F , so that the rows of the matrix M(f) associated with f : E ¡ú F are indexed by I, and the columns of the matrix M(f) are indexed by J. Obviously, this causes a mildly unpleasant reversal. If we had considered the bases (ui)i¡ÊI of E and (vj)j¡ÊJ of F , we would obtain a J ¡Á I-matrix M(f)=(aji)j¡ÊJ, i¡ÊI . No matter what we do, there will be a reversal! We decided to stick to the bases (uj)j¡ÊJ of E and (vi)i¡ÊI of F , so that we get an I ¡Á J-matrix M(f), knowing that we may occasionally su.er from this decision! 
When I and J are .nite, and say, when |I| = m and |J| = n, the linear map f is determined by the matrix M(f) whose entries in the j-th column are the components of the 
77 

vector f(uj) over the basis (v1,...,vm), that is, the matrix 
.
. 
M(f)= 

.... 

.. .
.
. ...
.
.. . 
am 1 am 2 ... amn 
.... 

whose entry on Row i and Column j is aij (1 ¡Ü i ¡Ü m,1 ¡Ü j ¡Ü n). 
We will now show that when E and F have .nite dimension, linear maps can be very conveniently represented by matrices, and that composition of linear maps corresponds to matrix multiplication. We will follow rather closely an elegant presentation method due to 
. 
Emil Artin. 
Let E and F be two vector spaces, and assume that E has a .nite basis (u1,...,un) and that F has a .nite basis (v1,...,vm). Recall that we have shown that every vector x ¡Ê E can be written in a unique way as 
x = x1u1 + ¡¤¡¤¡¤ + xnun, 
and similarly every vector y ¡Ê F can be written in a unique way as 
y = y1v1 + ¡¤¡¤¡¤ + ymvm. 
Let f : E ¡ú F be a linear map between E and F . Then for every x = x1u1 + ¡¤¡¤¡¤ + xnun in E, by linearity, we have 
f(x)= x1f(u1)+ ¡¤¡¤¡¤ + xnf(un). 
Let f(uj)= a1 jv1 + ¡¤¡¤¡¤ + amjvm, 
or more concisely, 
m f(uj)= aijvi, 
i=1 for every j,1 ¡Ü j ¡Ü n. This can be expressed by writing the coe.cients a1j,a2j,...,amj of f(uj) over the basis (v1,...,vm), as the jth column of a matrix, as shown below: 
)
. 
f(u1) f(unv1 a11 a12 ... a1n 
f(u2) 
... 

v2 
. 

. 

. 

.... 

a21 a22 ... a2n 
.. .
.
.... 

. 

.. 

.. 

..
. 
. 

vm am1 am2 ... amn Then substituting the right-hand side of each f(uj) into the expression for f(x), we get 
mm 
f(x)= x1( ai 1vi)+ ¡¤¡¤¡¤ + xn( ainvi), 
i=1 i=1 
3.1. REPRESENTATION OF LINEAR MAPS BY MATRICES 
which, by regrouping terms to obtain a linear combination of the vi, yields 
nn f(x)=( a1 jxj)v1 + ¡¤¡¤¡¤ +( amjxj)vm. j=1 j=1 
Thus, letting f(x)= y = y1v1 + ¡¤¡¤¡¤ + ymvm, we have 
n yi = aijxj (1) j=1 
for all i,1 ¡Ü i ¡Ü m. To make things more concrete, let us treat the case where n = 3 and m = 2. In this case, 
f(u1)= a11v1 + a21v2 f(u2)= a12v1 + a22v2 f(u3)= a13v1 + a23v2, 
which in matrix form is expressed by 
f(u1) f(u2) f(u3)
  
v1 a11 a12 a13 
, 
v2a21 a22 a23
and for any x = x1u1 + x2u2 + x3u3, we have 
f(x)= f(x1u1 + x2u2 + x3u3) = x1f(u1)+ x2f(u2)+ x3f(u3) = x1(a11v1 + a21v2)+ x2(a12v1 + a22v2)+ x3(a13v1 + a23v2) =(a11x1 + a12x2 + a13x3)v1 +(a21x1 + a22x2 + a23x3)v2. 
Consequently, since y = y1v1 + y2v2, we have 
y1 = a11x1 + a12x2 + a13x3 y2 = a21x1 + a22x2 + a23x3. 
This agrees with the matrix equation
       . x1 .  
y1  = a11  a12  a13  .x2..  
y2 a21  a22  a23 x3  

We now formalize the representation of linear maps by matrices. 

De.nition 3.1. Let E and F be two vector spaces, and let (u1,...,un) be a basis for E, and (v1,...,vm) be a basis for F . Each vector x ¡Ê E expressed in the basis (u1,...,un) as x = x1u1 + ¡¤¡¤¡¤ + xnun is represented by the column matrix 
.
. 

M(x)= 

.. 

x1 
. 
. 
. 
xn 
.. 

and similarly for each vector y ¡Ê F expressed in the basis (v1,...,vm). Every linear map f : E ¡ú F is represented by the matrix M(f)=(aij), where aij is the i-th component of the vector f(uj) over the basis (v1,...,vm), i.e., where 
m f(uj)= aijvi, for every j,1 ¡Ü j ¡Ü n. i=1 
The coe.cients a1j,a2j,...,amj of f(uj) over the basis (v1,...,vm) form the jth column of the matrix M(f) shown below: 
f() f(uu1n. vaaa... 111121n 
f(u2) 
... 

v2 
. 

. 

. 

.... 

a21 a22 ... a2n 
.. .
.
.. 

..
. 
. .... 
) 
. 

.. 

. 

vm am1 am2 ... 
amn 
The matrix M(f) associated with the linear map f : E ¡ú F is called the matrix of f with respect to the bases (u1,...,un) and (v1,...,vm). When E = F and the basis (v1,...,vm) is identical to the basis (u1,...,un) of E, the matrix M(f) associated with f : E ¡ú E (as above) is called the matrix of f with respect to the basis (u1,...,un). 
Remark: As in the remark after De.nition 2.14, there is no reason to assume that the vectors in the bases (u1,...,un) and (v1,...,vm) are ordered in any particular way. However, it is often convenient to assume the natural ordering. When this is so, authors sometimes refer to the matrix M(f) as the matrix of f with respect to the ordered bases (u1,...,un) and (v1,...,vm). 
Let us illustrate the representation of a linear map by a matrix in a concrete situation. Let E be the vector space R[X]4 of polynomials of degree at most 4, let F be the vector space R[X]3 of polynomials of degree at most 3, and let the linear map be the derivative map d: that is, 
d(P + Q)= dP + dQ d(¦ËP )= ¦ËdP, 

3.1. REPRESENTATION OF LINEAR MAPS BY MATRICES 
with ¦Ë ¡Ê R. We choose (1, x, x2,x3,x4) as a basis of E and (1, x, x2,x3) as a basis of F . Then the 4 ¡Á 5 matrix D associated with d is obtained by expressing the derivative dxi of each basis vector xi for i =0, 1, 2, 3, 4 over the basis (1, x, x2,x3). We .nd 
.
. 
D = 

... 

01000 
00200 

00030 
00004 

...

. 

If P denotes the polynomial 
P =3x 4 . 5x 3 + x 2 . 7x +5, 
we have dP = 12x 3 . 15x 2 +2x . 7. 
The polynomial P is represented by the vector (5, .7, 1, .5, 3), the polynomial dP is repre-sented by the vector (.7, 2, .15, 12), and we have 
.
. 

.
.
.
.

5 

.7

01000

..... 

..... 

= 

.7 

1 

.5 

... 

... 

... 

...

00200 

00030 

2 

.15 

, 

00004 12 
3 as expected! The kernel (nullspace) of d consists of the polynomials of degree 0, that is, the constant polynomials. Therefore dim(Ker d) = 1, and from dim(E) = dim(Ker d) + dim(Im d) (see Theorem 5.8), we get dim(Im d) = 4 (since dim(E) = 5). For fun, let us .gure out the linear map from the vector space R[X]3 to the vector space 
i, for i =0, 1, 2, 3).
¶Ö
R[X]4 given by integration (.nding the primitive, or anti-derivative) of xThe 5 ¡Á 4 matrix S representingwith respect to the same bases as before is 
.
. 
..... 

00 0 0 
10 0 0 

01/20 0 

001/30 
00 01/4 

..... 

S = 

. 

We verify that DS = I4, 
.
. .
..
.
0  0  0  0  
1  0  0  0  
0  1/2  0  0  

... 
01000 0 
1000

..... 

.....

001/30 

... 

... 

...

0020 

0100 

0010 

= 

. 

0003 

0 00004 
0001 

00 01/4 

This is to be expected by the fundamental theorem of calculus since the derivative of an integral returns the function. As we will shortly see, the above matrix product corresponds to this functional composition. The equation DS = I4 shows that S is injective and has D as a left inverse. However, SD 
= I5, and instead 
.
.
.
. .
.
0  0  0  0  
1  0  0  0  
0  1/2  0  0  

00000 

01000

..... 

..... 

= 

..... 

..... 

, 

001/30 

01000 

00100 

00010 

... 

... 

00200 

00030 

00004 

00 	01/4 00001 
because constant polynomials (polynomials of degree 0) belong to the kernel of D. 


3.2 	Composition of Linear Maps and Matrix Multipli-cation 
Let us now consider how the composition of linear maps is expressed in terms of bases. 
Let E, F , and G, be three vectors spaces with respective bases (u1,...,up) for E, (v1,...,vn) for F , and (w1,...,wm) for G. Let g : E ¡ú F and f : F ¡ú G be linear maps. As explained earlier, g : E ¡ú F is determined by the images of the basis vectors uj, and 
f : F ¡ú G is determined by the images of the basis vectors vk. We would like to understand how f . g : E ¡ú G is determined by the images of the basis vectors uj. 
Remark: Note that we are considering linear maps g : E ¡ú F and f : F ¡ú G, instead of f : E ¡ú F and g : F ¡ú G, which yields the composition f . g : E ¡ú G instead of g . f : E ¡ú G. Our perhaps unusual choice is motivated by the fact that if f is represented by a matrix M(f)=(aik) and g is represented by a matrix M(g)=(bkj), then f .g : E ¡ú G is represented by the product AB of the matrices A and B. If we had adopted the other choice where f : E ¡ú F and g : F ¡ú G, then g . f : E ¡ú G would be represented by the product BA. Personally, we .nd it easier to remember the formula for the entry in Row i and Column j of the product of two matrices when this product is written by AB, rather than BA. Obviously, this is a matter of taste! We will have to live with our perhaps unorthodox 
choice.  
Thus, let  
m  
f(vk) =  ai kwi,  
i=1  
for every k, 1 ¡Ü k ¡Ü n, and let  
n  
g(uj) =  bk jvk,  
k=1  

3.2. COMPOSITION OF LINEAR MAPS AND MATRIX MULTIPLICATION 
for every j,1 ¡Ü j ¡Ü p; in matrix form, we have 
f(v2) ... f(vn
f(v1)) w1 a11 a12 ... a1n 
.
. .... 

w2 
. 

. 

. 

.... 

a21 a22 ... a2n 
.. .
.
.. 

..
.
.. 

. 

wm am1 am2 ... amn 
and 
)
.
.
g(u1) g(upb11 b12 ... b1p g(u2) 
... 

v1 
.... 

.... 

b21 b22 ... b2pv2 
. 

. 

. 

.

.. .
.
.. 

..
.
.. 

. 

vn bn1 bn2 ... bnp 
By previous considerations, for every x = x1u1 + ¡¤¡¤¡¤ + xpup, letting g(x)= y = y1v1 + ¡¤¡¤¡¤ + ynvn, we have 
p yk = bkjxj (2) j=1 
for all k,1 ¡Ü k ¡Ü n, and for every y = y1v1 + ¡¤¡¤¡¤ + ynvn, letting f(y)= z = z1w1 + ¡¤¡¤¡¤ + zmwm, we have 
n zi = aikyk (3) k=1 
for all i,1 ¡Ü i ¡Ü m. Then if y = g(x) and z = f(y), we have z = f(g(x)), and in view of (2) and (3), we have 
np zi = aik( bkjxj) k=1 j=1 np 
= aikbkjxj k=1 j=1 pn 
= aikbkjxj j=1 k=1 pn =( aikbkj)xj. j=1 k=1 
Thus, de.ning cij such that 
n 
cij = aikbkj, k=1 
for 1 ¡Ü i ¡Ü m, and 1 ¡Ü j ¡Ü p, we have 
p zi = cijxj (4) j=1 
Identity (4) shows that the composition of linear maps corresponds to the product of matrices. Then given a linear map f : E ¡ú F represented by the matrix M(f)=(aij) w.r.t. the bases (u1,...,un) and (v1,...,vm), by Equation (1), namely n yi = aijxj 1 ¡Ü i ¡Ü m, j=1 
and the de.nition of matrix multiplication, the equation y = f(x) corresponds to the matrix equation M(y)= M(f)M(x), that is, 
.
..
.
.. 
y1 a11 ... a1 n x1 
..

. 

. 

. 

..

= 

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
. 

..

.

. 

. 

ym am 1 ... amn xn 
Recall that 
..
..
..
.
..
. 
a11 a12 ... a1 n x1 a11 a12 a1 n 
.... 

.... 
.... 

x2 
. 

. 

. 

.... 

= x1 
.... 

a21 
. 

. 

. 

.... 

+ x2 
.... 

a22 
. 

. 

. 

.... 

+ ¡¤¡¤¡¤ + xn 
.... 

a2 n 
. 

. 

. 

.... 

a21 a22 ... a2 n 
.

.. .
.
.. 

..

.

.. 

. 

am 1 am 2 ... amn xn am 1 am 2 amn 
Sometimes, it is necessary to incorporate the bases (u1,...,un) and (v1,...,vm) in the notation for the matrix M(f) expressing f with respect to these bases. This turns out to be a messy enterprise! 
We propose the following course of action: 
De.nition 3.2. Write U =(u1,...,un) and V =(v1,...,vm) for the bases of E and F , and denote by MU,V (f) the matrix of f with respect to the bases U and V. Furthermore, write xU for the coordinates M(x)=(x1,...,xn) of x ¡Ê E w.r.t. the basis U and write yV for the coordinates M(y)=(y1,...,ym) of y ¡Ê F w.r.t. the basis V . Then 
y = f(x) 
is expressed in matrix form by yV = MU,V (f) xU . When U = V, we abbreviate MU,V (f) as MU (f). 

3.2. COMPOSITION OF LINEAR MAPS AND MATRIX MULTIPLICATION 
The above notation seems reasonable, but it has the slight disadvantage that in the expression MU,V (f)xU , the input argument xU which is fed to the matrix MU,V (f) does not appear next to the subscript U in MU,V (f). We could have used the notation MV,U (f), and some people do that. But then, we .nd a bit confusing that V comes before U when f maps from the space E with the basis U to the space F with the basis V. So, we prefer to use the notation MU,V (f). 
Be aware that other authors such as Meyer [48] use the notation [f]U,V , and others such as Dummit and Foote [19] use the notation MUV (f), instead of MU,V (f). This gets worse! You may .nd the notation MVU (f) (as in Lang [41]), or U [f]V , or other strange notations. 
De.nition 3.2 shows that the function which associates to a linear map f : E ¡ú F the matrix M(f) w.r.t. the bases (u1,...,un) and (v1,...,vm) has the property that matrix mul-tiplication corresponds to composition of linear maps. This allows us to transfer properties of linear maps to matrices. Here is an illustration of this technique: 
Proposition 3.1. (1) Given any matrices A ¡Ê Mm,n(K), B ¡Ê Mn,p(K), and C ¡Ê Mp,q(K), we have 
(AB)C = A(BC); 
that is, matrix multiplication is associative. 
(2) Given any matrices A, B ¡Ê Mm,n(K), and C, D ¡Ê Mn,p(K), for all ¦Ë ¡Ê K, we have 
(A + B)C = AC + BC A(C + D)= AC + AD (¦ËA)C = ¦Ë(AC) A(¦ËC)= ¦Ë(AC), 
so that matrix multiplication ¡¤:Mm,n(K) ¡Á Mn,p(K) ¡ú Mm,p(K) is bilinear. 
Proof. (1) Every m ¡Á n matrix A =(aij) de.nes the function fA : Kn ¡ú Km given by 
fA(x)= Ax, 
for all x ¡Ê Kn . It is immediately veri.ed that fA is linear and that the matrix M(fA) representing fA over the canonical bases in Kn and Km is equal to A. Then Formula (4) proves that 
M(fA . fB)= M(fA)M(fB)= AB, 
so we get M((fA . fB) . fC )= M(fA . fB)M(fC )=(AB)C 
and M(fA . (fB . fC)) = M(fA)M(fB . fC )= A(BC), 
and since composition of functions is associative, we have (fA . fB) . fC = fA . (fB . fC ), which implies that 
(AB)C = A(BC). 
(2) It is immediately veri.ed that if f1,f2 ¡Ê HomK (E, F ), A, B ¡Ê Mm,n(K), (u1,...,un) is any basis of E, and (v1,...,vm) is any basis of F , then 
M(f1 + f2)= M(f1)+ M(f2) fA+B = fA + fB. 
Then we have 
(A + B)C = M(fA+B)M(fC ) 
= M(fA+B . fC ) 
= M((fA + fB) . fC )) 
= M((fA . fC )+(fB . fC )) 
= M(fA . fC )+ M(fB . fC ) 
= M(fA)M(fC)+ M(fB)M(fC ) 
= AC + BC. 
The equation A(C + D)= AC + AD is proven in a similar fashion, and the last two equations are easily veri.ed. We could also have veri.ed all the identities by making matrix computations. 
Note that Proposition 3.1 implies that the vector space Mn(K) of square matrices is a (noncommutative) ring with unit In. (It even shows that Mn(K) is an associative algebra.) 
The following proposition states the main properties of the mapping f ¡ú M(f) between Hom(E, F ) and Mm,n. In short, it is an isomorphism of vector spaces. 
Proposition 3.2. Given three vector spaces E, F , G, with respective bases (u1,...,up), (v1,...,vn), and (w1,...,wm), the mapping M : Hom(E, F ) ¡ú Mn,p that associates the ma-trix M(g) to a linear map g : E ¡ú F satis.es the following properties for all x ¡Ê E, all g, h: E ¡ú F , and all f : F ¡ú G: 
M(g(x)) = M(g)M(x) M(g + h)= M(g)+ M(h) M(¦Ëg)= ¦ËM(g) M(f . g)= M(f)M(g), 
where M(x) is the column vector associated with the vector x and M(g(x)) is the column vector associated with g(x), as explained in De.nition 3.1. 
Thus, M : Hom(E, F ) ¡ú Mn,p is an isomorphism of vector spaces, and when p = n and the basis (v1,...,vn) is identical to the basis (u1,...,up), M : Hom(E, E) ¡ú Mn is an isomorphism of rings. 
Proof. That M(g(x)) = M(g)M(x) was shown by De.nition 3.2 or equivalently by Formula (1). The identities M(g + h)= M(g)+ M(h) and M(¦Ëg)= ¦ËM(g) are straightforward, and 

3.3. CHANGE OF BASIS MATRIX 
M(f . g)= M(f)M(g) follows from Identity (4) and the de.nition of matrix multiplication. The mapping M : Hom(E, F ) ¡ú Mn,p is clearly injective, and since every matrix de.nes a linear map (see Proposition 3.1), it is also surjective, and thus bijective. In view of the above identities, it is an isomorphism (and similarly for M : Hom(E, E) ¡ú Mn, where Proposition 
3.1 is used to show that Mn is a ring). 
In view of Proposition 3.2, it seems preferable to represent vectors from a vector space of .nite dimension as column vectors rather than row vectors. Thus, from now on, we will denote vectors of Rn (or more generally, of Kn) as column vectors. 


3.3 Change of Basis Matrix 
It is important to observe that the isomorphism M : Hom(E, F ) ¡ú Mn,p given by Proposition 
3.2 depends on the choice of the bases (u1,...,up) and (v1,...,vn), and similarly for the isomorphism M : Hom(E, E) ¡ú Mn, which depends on the choice of the basis (u1,...,un). Thus, it would be useful to know how a change of basis a.ects the representation of a linear map f : E ¡ú F as a matrix. The following simple proposition is needed. 
Proposition 3.3. Let E be a vector space, and let (u1,...,un) be a basis of E. For every 
 
family (v1,...,vn), let P =(aij) be the matrix de.ned such that vj =n The matrix 
i=1 aijui. 
P is invertible i. (v1,...,vn) is a basis of E. 
Proof. Note that we have P = M(f), the matrix associated with the unique linear map 
f : E ¡ú E such that f(ui)= vi. By Proposition 2.17, f is bijective i. (v1,...,vn) is a basis of E. Furthermore, it is obvious that the identity matrix In is the matrix associated with the identity id: E ¡ú E w.r.t. any basis. If f is an isomorphism, then f .f.1 = f.1 .f = id, and by Proposition 3.2, we get M(f)M(f.1)= M(f.1)M(f)= In, showing that P is invertible and that M(f.1)= P .1 . 
Proposition 3.3 suggests the following de.nition. 
De.nition 3.3. Given a vector space E of dimension n, for any two bases (u1,...,un) and (v1,...,vn) of E, let P =(aij) be the invertible matrix de.ned such that 
n vj = aijui, i=1 
which is also the matrix of the identity id: E ¡ú E with respect to the bases (v1,...,vn) and (u1,...,un), in that order. Indeed, we express each id(vj)= vj over the basis (u1,...,un). The coe.cients a1j,a2j,...,anj of vj over the basis (u1,...,un) form the jth column of the matrix P shown below: 
.

. 

u1 a11 a12 ... a1n 
u2 
. 

. 

. 

.... 

.... 

a21 a22 ... a2n 
.

.. .
.
.. 

..
.
.. 

. 

un an1 an2 ... ann 
The matrix P is called the change of basis matrix from (u1,...,un) to (v1,...,vn). 
Clearly, the change of basis matrix from (v1,...,vn) to (u1,...,un) is P .1 . Since P = (aij) is the matrix of the identity id: E ¡ú E with respect to the bases (v1,...,vn) and (u1,...,un), given any vector x ¡Ê E, if x = x1u1 + ¡¤¡¤¡¤ + xnun over the basis (u1,...,un) and 
x = x1Ôýv1 + ¡¤¡¤¡¤ + xnÔý vn over the basis (v1,...,vn), from Proposition 3.2, we have 
.
..
.
.. 
x1 a11 ... a1 n x1 Ôý 
..

. 

. 

. 

..

= 

.. 

.. 
..

. 

. 

. 

..

, 

..
.
. 

..

.

. 

. 

xn an 1 ... ann xÔý
n 
showing that the old coordinates (xi) of x (over (u1,...,un)) are expressed in terms of the new coordinates (xiÔý ) of x (over (v1,...,vn)). 
Now we face the painful task of assigning a ¡°good¡± notation incorporating the bases U =(u1,...,un) and V =(v1,...,vn) into the notation for the change of basis matrix from U to V. Because the change of basis matrix from U to V is the matrix of the identity map idE with respect to the bases V and U in that order, we could denote it by MV,U (id) (Meyer 
[48] uses the notation [I]V,U ). We prefer to use an abbreviation for MV,U (id). 
De.nition 3.4. The change of basis matrix from U to V is denoted 
PV,U . 
Note that = P .1
PU,VV,U . 
Then, if we write xU =(x1,...,xn) for the old coordinates of x with respect to the basis U 
ÔýÔý
and xV =(x1,...,xn) for the new coordinates of x with respect to the basis V, we have 
= P .1 
xU = PV,U xV ,xVV,U xU . 
The above may look backward, but remember that the matrix MU,V (f) takes input expressed over the basis U to output expressed over the basis V. Consequently, PV,U takes input expressed over the basis V to output expressed over the basis U, and xU = PV,U xV matches this point of view! 
3.3. CHANGE OF BASIS MATRIX 
 Beware that some authors (such as Artin [3]) de.ne the change of basis matrix from U to V as PU,V = PV.,1 U . Under this point of view, the old basis U is expressed in terms of the new basis V. We .nd this a bit unnatural. Also, in practice, it seems that the new basis is often expressed in terms of the old basis, rather than the other way around. 
Since the matrix P = PV,U expresses the new basis (v1,...,vn) in terms of the old basis (u1,..., un), we observe that the coordinates (xi) of a vector x vary in the opposite direction of the change of basis. For this reason, vectors are sometimes said to be contravariant. However, this expression does not make sense! Indeed, a vector in an intrinsic quantity that does not depend on a speci.c basis. What makes sense is that the coordinates of a vector vary in a contravariant fashion. 
Let us consider some concrete examples of change of bases. 
...
Example 3.1. Let E = F = R2, with u1 = (1, 0), u2 = (0, 1), v1 = (1, 1) and v2 =(.1, 1). The change of basis matrix P from the basis U =(u1,u2) to the basis V =(v1,v2) is 
1 .1 
P = 
11 
and its inverse is 
... 
1/21/2 
P .1 
= .
.1/21/2 
The old coordinates (x1,x2) with respect to (u1,u2) are expressed in terms of the new coordinates (x1Ôý ,x2Ôý ) with respect to (v1,v2) by 
x1 1 .1 xÔý 
1 
x2 =1 1 x2 Ôý , 
and the new coordinates (x1Ôý ,x2Ôý ) with respect to (v1,v2) are expressed in terms of the old 
coordinates (x1,x2) with respect to (u1,u2) by 
xÔý 1/21/2 x1
1 
x2 Ôý = .1/21/2 x2 . 
Example 3.2. Let E = F = R[X]3 be the set of polynomials of degree at most 3, and consider the bases U = (1, x, x2,x3) and V =(B3(x),B3(x),B23(x),B33(x)), where B03(x),B13(x),B23(x),B33(x)aretheBernsteinpolynomialso0fdegre1e 3, given by 
B3 B3 23
(x) = (1 . x)3 (x) = 3(1 . x)2 xB3(x) = 3(1 . x)xB3(x)= x.
01 2 3 
By expanding the Bernstein polynomials, we .nd that the change of basis matrix PV,U is given by 1 0 00 .33 00 
PV,U = . 
3 .630 .13 .31 
.
. 
We also .nd that the inverse of PV,U is 
.
. 
10 00 11/300 
P .1 
= 	.
V,U 	12/31/30 11 11 
Therefore, the coordinates of the polynomial 2x3 . x + 1 over the basis V are 
...
... .
..
... 
1 1000 1 

... 

2/3 

1/3 

= 
... ... 
11/300 

12/31/30 

... 
... 

..1..
0 
, 

and so 

2 1111 2 

2x 3 . x +1= B03(x)+ 2 B13(x)+ 1 B23(x)+2B33(x). 
33 
3.4 	The E.ect of a Change of Bases on Matrices 
The e.ect of a change of bases on the representation of a linear map is described in the following proposition. 
Proposition 3.4. Let E and F be vector spaces, let U =(u1,...,un) and UÔý =(u1Ôý ,...,uÔý )
n
be two bases of E, and let V =(v1,...,vm) and VÔý =(v1Ôý ,...,vmÔý ) be two bases of F . Let P = PUy,U be the change of basis matrix from U to UÔý, and let Q = PVy,V be the change of basis matrix from V to VÔý . For any linear map f : E ¡ú F , let M(f)= MU,V (f) be the matrix associated to f w.r.t. the bases U and V, and let MÔý(f)= MUy,Vy (f) be the matrix associated 
to f w.r.t. the bases UÔý  and VÔý .  We have  
MÔý(f) = Q.1M(f)P,  
or more explicitly  

MUyy (f)= PV.1 ,V MU,V (f)PUyy MU,V (f)PUy
,Vy,U = PV,V,U . 
Proof. Since f : E ¡ú F can be written as f = idF . f . idE, since P is the matrix of idE 
w.r.t. the bases (uÔý 1,...,uÔý n) and (u1,...,un), and Q.1 is the matrix of idF w.r.t. the bases (v1,...,vm) and (v1Ôý ,...,vÔý ), by Proposition 3.2, we have MÔý(f)= Q.1M(f)P .
m
As a corollary, we get the following result. 
Corollary 3.5. Let E be a vector space, and let U =(u1,...,un) and UÔý =(uÔý 1,...,uÔý ) be 
n
two bases of E. Let P = PUy,U be the change of basis matrix from U to UÔý . For any linear 

3.4. THE EFFECT OF A CHANGE OF BASES ON MATRICES 
map f : E ¡ú E, let M(f)= MU (f) be the matrix associated to f w.r.t. the basis U, and let MÔý(f)= MUy (f) be the matrix associated to f w.r.t. the basis UÔý . We have 
MÔý(f)= P .1M(f)P, 
or more explicitly, 
MUy (f)= PU.y1 MU (f)PUy,U = PU,Uy MU (f)PUy,U . 
,U 
Example 3.3. Let E = R2 , U =(e1,e2) where e1 = (1, 0) and e2 = (0, 1) are the canonical basis vectors, let V =(v1,v2)=(e1,e1 . e2), and let 
21 
A = . 
01 
The change of basis matrix P = PV,U from U to V is 
11 
P = ,
0 .1 
and we check that 
P .1 = P. 
Therefore, in the basis V, the matrix representing the linear map f de.ned by A is 
AÔý = P .1AP = P AP =11 21 11 =20= D, 0 .1 01 0 .1 01 
a diagonal matrix. In the basis V, it is clear what the action of f is: it is a stretch by a factor of 2 in the v1 direction and it is the identity in the v2 direction. Observe that v1 and v2 are not orthogonal. 
What happened is that we diagonalized the matrix A. The diagonal entries 2 and 1 are the eigenvalues of A (and f), and v1 and v2 are corresponding eigenvectors. We will come back to eigenvalues and eigenvectors later on. 
The above example showed that the same linear map can be represented by di.erent matrices. This suggests making the following de.nition: 
De.nition 3.5. Two n¡Án matrices A and B are said to be similar i. there is some invertible matrix P such that 
B = P .1AP. 
It is easily checked that similarity is an equivalence relation. From our previous consid-erations, two n ¡Á n matrices A and B are similar i. they represent the same linear map with respect to two di.erent bases. The following surprising fact can be shown: Every square matrix A is similar to its transpose AT. The proof requires advanced concepts (the Jordan form or similarity invariants). 
If U =(u1,...,un) and V =(v1,...,vn) are two bases of E, the change of basis matrix 
.
. 
P = PV,U = 
.... 

.. .
.
....
.
.. . 
an1 an2 ¡¤¡¤¡¤ ann 
.... 

from (u1,...,un) to (v1,...,vn) is the matrix whose jth column consists of the coordinates of vj over the basis (u1,...,un), which means that 
n vj = aijui. i=1 
.
. 

.. 

v1 
. 
. 
. 
vn 
..

in En 
as the 

It is natural to extend the matrix notation and to express the vector 

.
. 

product of a matrix times the vector 

.. 

u1 
. 
. 
. 
un 
..

in En, namely as 
.
..
... 
v1 a11 a21 ¡¤¡¤¡¤ an1 u1 
.... 

v2 
. 

. 

. 

.... 

= 

.... 

.... 
.... 

u2 
. 

. 

. 

....

a12 a22 
¡¤¡¤¡¤ 

an2 
,

.. .
.
.. 

..

.

.. 

. 

vn a1n a2n ¡¤¡¤¡¤ ann un 
but notice that the matrix involved is not P , but its transpose P T . This observation has the following consequence: if U =(u1,...,un) and V =(v1,...,vn) 
are two bases of E and if

..
.. 
v1 u1 
. 

. 

. 

..

= A

. 

. 

. 

, 

vn  un  
that is,  
n  
vi =  aijuj,  
j=1  
for any vector w ¡Ê E, if  
n  n  
w =  xiui =  ykvk,  
i=1  k=1  

then

..
.. 
x1 y1 
..

. 

. 

. 

..

= AT
..

. 

. 

. 

..

, 

xn yn 

3.4. THE EFFECT OF A CHANGE OF BASES ON MATRICES 
and so

..
.. 
y1 x1 
..

. 

. 

..

..

. 

. 

..

.

. 

=(AT).1 
. 

yn xn 
It is easy to see that (AT).1 =(A.1)T . Also, if U =(u1,...,un), V =(v1,...,vn), and W =(w1,...,wn) are three bases of E, and if the change of basis matrix from U to V is P = PV,U and the change of basis matrix from V to W is Q = PW,V , then 
..
.
...
.. 
v1 u1 w1 v1 
..

. 

. 

. 

..

= P T
..

. 

. 

. 

..

, 

..

. 

. 

. 

..

= QT
..

. 

. 

. 

..

, 

vn un wn vn 
so

..
..
.. 
w1 u1 u1 
..

. 

. 

. 

..

= QTP T
..

. 

. 

. 

..

=(PQ)T
..

. 

. 

. 

..

, 

wn un un 
which means that the change of basis matrix PW,U from U to W is PQ. This proves that 
PW,U = PV,U PW,V . 
Even though matrices are indispensable since they are the major tool in applications of linear algebra, one should not lose track of the fact that 

linear maps are more fundamental because they are intrinsic objects that do not depend on the choice of bases. Consequently, we advise the reader to try to think in terms of linear maps rather than reduce everything to matrices. 
In our experience, this is particularly e.ective when it comes to proving results about linear maps and matrices, where proofs involving linear maps are often more ¡°conceptual.¡± These proofs are usually more general because they do not depend on the fact that the dimension is .nite. Also, instead of thinking of a matrix decomposition as a purely algebraic operation, it is often illuminating to view it as a geometric decomposition. This is the case of the SVD, which in geometric terms says that every linear map can be factored as a rotation, followed by a rescaling along orthogonal axes and then another rotation. 
After all, 

a matrix is a representation of a linear map, 
and most decompositions of a matrix re.ect the fact that with a suitable choice of a basis (or bases), the linear map is a represented by a matrix having a special shape. The problem is then to .nd such bases. 
Still, for the beginner, matrices have a certain irresistible appeal, and we confess that it takes a certain amount of practice to reach the point where it becomes more natural to deal with linear maps. We still recommend it! For example, try to translate a result stated in terms of matrices into a result stated in terms of linear maps. Whenever we tried this exercise, we learned something. 
Also, always try to keep in mind that 
linear maps are geometric in nature; they act on space. 


3.5 Summary 
The main concepts and results of this chapter are listed below: 
. 	
The representation of linear maps by matrices. 

. 	
The matrix representation mapping M : Hom(E, F ) ¡ú Mn,p and the representation isomorphism (Proposition 3.2). 

. 	
Change of basis matrix and Proposition 3.4. 



3.6 Problems 
Problem 3.1. Prove that the column vectors of the matrix A1 given by 
.. 
123 
..
A1 = 	237 131 
are linearly independent. Prove that the coordinates of the column vectors of the matrix B1 over the basis consisting 
of the column vectors of A1 given by B1 =  . . 3 1 4  5 2 3  1 1 .6  . .  
are the columns of the matrix P1 given by .  .  

.27 .61 .41 
.	.
P1 =9 18 9 . 4 10 8 
Give a nontrivial linear dependence of the columns of P1. Check that B1 = A1P1. Is the matrix B1 invertible? 
3.6. PROBLEMS 
Problem 3.2. Prove that the column vectors of the matrix A2 given by 
.
. 
A2 = 
... 

1111 
1213 

1122 
1113 

... 

are linearly independent. 
Prove that the column vectors of the matrix B2 given by 

.
. 
B2 = 
... 

1 .22 .2 
0 .32 .3 

3 .55 .4 
3 .44 .4 

... 

are linearly independent. 
Prove that the coordinates of the column vectors of the matrix B2 over the basis consisting of the column vectors of A2 are the columns of the matrix P2 given by 
.
. 
... 

201 .1 .31 .21 
1 .22 .1 
1 .11 .1 

...

.

P2 = 
Check that A2P2 = B2. Prove that 
.
. 
... 

.1 .1 .11 211 .2 
212 .3 
.1 .10 .1 

...

P .1 
= 
2 
. 

What are the coordinates over the basis consisting of the column vectors of B2 of the vector whose coordinates over the basis consisting of the column vectors of A1 are (2, .3, 0, 0)? 
Problem 3.3. Consider the polynomials 
B2 B2 B22 
0 (t) = (1 . t)21 (t) = 2(1 . t)t 2 (t)= tB3 B3 B32 B33
(t) = (1 . t)3 (t) = 3(1 . t)2t (t) = 3(1 . t)t(t)= t,
01 2 3 
known as the Bernstein polynomials of degree 2 and 3. 
(1) Show that the Bernstein polynomials B2(t),B2(t),B2(t) are expressed as linear com-
012 binations of the basis (1, t, t2) of the vector space of polynomials of degree at most 2 as 
follows:

.
..
..
. 
B02(t)1 .21 1 
.

B2 
1 (t)
.

=

.

02 

.2

.
.

t

.

. 

B22(t) 001 t2 
Prove that 
20
(t)+ B

21
(t)+ B

22
(t)=1. 

(2) Show that the Bernstein polynomials B

30
(t),B

31
(t),B

32
(t),B

33
(t) are expressed as linear 

combinations of the basis (1, t, t2,t3) of the vector space of polynomials of degree at most 3 
as follows:

.
..
..
. 
30
B

(t)1 

.33 .1 

1 

... 

(t) 


(t) 



... 

= 

... 

... 
... 

...

313233
B

B

03 

.6 

3 

t 

t2 
. 

00 3 

.3 

t3
B

(t) 0001 

Prove that 
B
30
(t)+ B

31
(t)+ B

32
(t)+ B

33
(t)=1. 

(3) Prove that the Bernstein polynomials of degree 2 are linearly independent, and that the Bernstein polynomials of degree 3 are linearly independent. 
 
 

Problem 3.4. Recall that the binomial coe.cientm is given by 
k
mm! 
= ,
kk!(m . k)!with 0 ¡Ü k ¡Ü m. For any m ¡Ý 1, we have the m +1 Bernstein polynomials of degree m given by m 
Bmk 
k (t)= (1 . t)m.kt, 0 ¡Ü k ¡Ü m. 
k 
(1) Prove that 
m mj
Bkm(t)= (.1)j.k tj. (.)
jk 
j=k 
0
Use the above to prove that Bm 
(2) Prove that (t),...,Bm(t) are linearly independent. 
m 
0
Bm 
(3) What can you say about the symmetries of the (m + 1) ¡Á (m + 1) matrix expressing (t)+ ¡¤¡¤¡¤ + Bm(t)=1.
m 
00
Bm 
the (k+1)th row of the (j +1)th column, since 0 ¡Ü k, j ¡Ü m. Make appropriate modi.cations to the indices). 
What can you say about the sum of the entries on each row of the above matrix? What about the sum of the entries on each column? 
(4) The purpose of this question is to express the ti in terms of the Bernstein polynomials 
Bm 
,...,Bm in terms of the basis 1, t, . . . , tm?
m 
Prove your claim (beware that in equation (.) the coe.cient of tj in Bkm is the entry on 
(t),...,Bm(t), with 0 ¡Ü i ¡Ü m.
m 

3.6. PROBLEMS 
First, prove that 
m.i iiBm.i
t= tj (t), 0 ¡Ü i ¡Ü m. 
j=0 
Then prove that 
mm . i mi + j 
= . 
ij i + ji Use the above facts to prove that 
m.ii+j ti = i Bm (t). 
mi+jj=0 i 
Conclude that the Bernstein polynomials Bm(t),...,Bm(t) form a basis of the vector 
0 m 
space of polynomials of degree ¡Ü m. 
Compute the matrix expressing 1, t, t2 in terms of B02(t),B12(t),B22(t), and the matrix expressing 1, t, t2,t3 in terms of B3(t),B3(t),B3(t),B3(t).
0123 
You should .nd

.
. 
1 .
0 001 
11 

1/21

. 

and

.
. 
... 

11 11 
01/32/31 

001/31 
00 01 

...

. 

(5) A polynomial curve C(t) of degree m in the plane is the set of points 
x(t)

C(t) = given by two polynomials of degree ¡Ü m, 
y(t) 
x(t)= ¦Á0tm1 + ¦Á1tm1.1 + ¡¤¡¤¡¤ + ¦Ám1 y(t)= ¦Â0tm2 + ¦Â1tm2.1 + ¡¤¡¤¡¤ + ¦Âm2 , 
with 1 ¡Ü m1,m2 ¡Ü m and ¦Á0,¦Â0 = 0. Prove that there exist m + 1 points b0,...,bm ¡Ê R2 so that 
x(t)
C(t)= = Bm(t)b0 + Bm(t)b1 + ¡¤¡¤¡¤ + Bm(t)bm
01 m
y(t) 
for all t ¡Ê R, with C(0) = b0 and C(1) = bm. Are the points b1,...,bm.1 generally on the curve? We say that the curve C is a B¡äezier curve and (b0,...,bm) is the list of control points of the curve (control points need not be distinct). 
Remark: Because Bm(t)+ ¡¤¡¤¡¤ + Bm(t) = 1 and Bm(t) ¡Ý 0 when t ¡Ê [0, 1], the curve 
0 mi 
segment C[0, 1] corresponding to t ¡Ê [0, 1] belongs to the convex hull of the control points. This is an important property of B¡äezier curves which is used in geometric modeling to .nd the intersection of curve segments. B¡äezier curves play an important role in computer graphics and geometric modeling, but also in robotics because they can be used to model the trajectories of moving objects. 
Problem 3.5. Consider the n ¡Á n matrix 
.
. 
00 0 ¡¤¡¤¡¤ 0 .an 10 0 ¡¤¡¤¡¤ 0 .an.1 01 0 ¡¤¡¤¡¤ 0 .an.2 
A = 

........ 

........ 

,

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
with an = 0. 
(1) Find a matrix P such that 
AT 
= P .1AP. 
What happens when an = 0? 
Hint. First, try n =3, 4, 5. Such a matrix must have zeros above the ¡°antidiagonal,¡± and 
identical entries pij for all i, j ¡Ý 0 such that i + j = n + k, where k =1,...,n. 

(2) Prove that if an = 1 and if a1,...,an.1 are integers, then P can be chosen so that the entries in P .1 are also integers. 
Problem 3.6. For any matrix A ¡Ê Mn(C), let RA and LA be the maps from Mn(C) to itself de.ned so that 
LA(B)= AB, RA(B)= BA, for all B ¡Ê Mn(C). 
(1) Check that LA and RA are linear, and that LA and RB commute for all A, B. Let adA :Mn(C) ¡ú Mn(C) be the linear map given by 
adA(B)= LA(B) . RA(B)= AB . BA =[A, B], for all B ¡Ê Mn(C). 
Note that [A, B] is the Lie bracket. 
(2) 
Prove that if A is invertible, then LA and RA are invertible; in fact, (LA).1 = LA.1 and (RA).1 = RA.1 . Prove that if A = P BP .1 for some invertible matrix P , then 

LA = LP . LB . LP .1 ,RA = RP .1 . RB . RP . 

(3) 
Recall that the n2 matrices Eij de.ned such that all entries in Eij are zero except the (i, j)th entry, which is equal to 1, form a basis of the vector space Mn(C). Consider the partial ordering of the Eij de.ned such that for i =1,...,n, if n ¡Ý j>k ¡Ý 1, then then Eij precedes Eik, and for j =1,...,n, if 1 ¡Ü i<h ¡Ü n, then Eij precedes Ehj. 



3.6. PROBLEMS 
Draw the Hasse diagram of the partial order de.ned above when n = 3. 
There are total orderings extending this partial ordering. How would you .nd them algorithmically? Check that the following is such a total order: 
(1, 3), (1, 2), (1, 1), (2, 3), (2, 2), (2, 1), (3, 3), (3, 2), (3, 1). 
(4) Let the total order of the basis (Eij) extending the partial ordering de.ned in (2) be 
given by 

 

i = h and j>k 
(i, j) < (h, k) i.
or i<h. Let R be the n ¡Á n permutation matrix given by 
.
. 
00 ... 01 

00 ... 10 .. ..


.
R = 

...... 

...... 

.

.. 

... 

.. 

. 

.. 

01 ... 00 

10 ... 00 


Observe that R.1 = R. Prove that for any n ¡Ý 1, the matrix of LA is given by A.In, and the matrix of RA is given by In . RATR (over the basis (Eij) ordered as speci.ed above), where . is the Kronecker product (also called tensor product) of matrices de.ned in De.nition 4.4. Hint. Figure out what are RB(Eij)= EijB and LB(Eij)= BEij. 
(5) Prove that if A is upper triangular, then the matrices representing LA and RA are 
also upper triangular. Note that if instead of the ordering 
E1n,E1n.1,...,E11,E2n,...,E21,...,Enn,...,En1, 
that I proposed you use the standard lexicographic ordering 
E11,E12,...,E1n,E21,...,E2n,...,En1,...,Enn, 
then the matrix representing LA is still A . In, but the matrix representing RA is In . AT . In this case, if A is upper-triangular, then the matrix of RA is lower triangular. This is the motivation for using the .rst basis (avoid upper becoming lower). 



