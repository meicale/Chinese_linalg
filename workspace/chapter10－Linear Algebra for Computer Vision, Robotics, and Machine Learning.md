Chapter 10 
The Dual Space and Duality 
In this chapter all vector spaces are de.ned over an arbitrary .eld K. For the sake of concreteness, the reader may safely assume that K = R. 
10.1 The Dual Space E. and Linear Forms 
In Section 2.8 we de.ned linear forms, the dual space E. = Hom(E, K) of a vector space E, and showed the existence of dual bases for vector spaces of .nite dimension. 
In this chapter we take a deeper look at the connection between a space E and its dual space E. . As we will see shortly, every linear map f : E ¡ú F gives rise to a linear map 
.
fT : F ¡ú E., and it turns out that in a suitable basis, the matrix of fT is the transpose of the matrix of f. Thus, the notion of dual space provides a conceptual explanation of the phenomena associated with transposition. 
But it does more, because it allows us to view a linear equation as an element of the dual space E., and thus to view subspaces of E as solutions of sets of linear equations and vice-versa. The relationship between subspaces and sets of linear forms is the essence of duality, a term which is often used loosely, but can be made precise as a bijection between the set of subspaces of a given vector space E and the set of subspaces of its dual E. . In this correspondence, a subspace V of E yields the subspace V 0 of E. consisting of all linear forms that vanish on V (that is, have the value zero for all input in V ). 
Consider the following set of two ¡°linear equations¡± in R3 , 
x . y + z =0 x . y . z =0, 
and let us .nd out what is their set V of common solutions (x, y, z) ¡Ê R3 . By subtracting the second equation from the .rst, we get 2z = 0, and by adding the two equations, we .nd that 2(x . y) = 0, so the set V of solutions is given by 
y = x z =0. 
343 

This is a one dimensional subspace of R3 . Geometrically, this is the line of equation y = x in the plane z = 0 as illustrated by Figure 10.1. 

Figure 10.1: The intersection of the magenta plane x . y + z = 0 with the blue-gray plane x . y . z = 0 is the pink line y = x. 
Now why did we say that the above equations are linear? Because as functions of (x, y, z), both maps f1 :(x, y, z) ¡ú x . y + z and f2 :(x, y, z) ¡ú x . y . z are linear. The set of all such linear functions from R3 to R is a vector space; we used this fact to form linear combinations of the ¡°equations¡± f1 and f2. Observe that the dimension of the subspace V is 1. The ambient space has dimension n = 3 and there are two ¡°independent¡± equations f1,f2, so it appears that the dimension dim(V ) of the subspace V de.ned by m independent equations is 
dim(V )= n . m, 
which is indeed a general fact (proven in Theorem 10.4). 
More generally, in Rn, a linear equation is determined by an n-tuple (a1,...,an) ¡Ê Rn , and the solutions of this linear equation are given by the n-tuples (x1,...,xn) ¡Ê Rn such that 
a1x1 + ¡¤¡¤¡¤ + anxn = 0; 
these solutions constitute the kernel of the linear map (x1,...,xn) ¡ú a1x1 + ¡¤¡¤¡¤ + anxn. The above considerations assume that we are working in the canonical basis (e1,...,en) of Rn, but we can de.ne ¡°linear equations¡± independently of bases and in any dimension, by viewing them as elements of the vector space Hom(E, K) of linear maps from E to the .eld 
K. 
De.nition 10.1. Given a vector space E, the vector space Hom(E, K) of linear maps from E to the .eld K is called the dual space (or dual) of E. The space Hom(E, K) is also denoted 
10.1. THE DUAL SPACE E. AND LINEAR FORMS 345 
by E., and the linear maps in E. are called the linear forms, or covectors. The dual space E.. 
of the space E. is called the bidual of E. 
As a matter of notation, linear forms f : E ¡ú K will also be denoted by starred symbol, such as u . , x . , etc. 
Given a vector space E and any basis (ui)i¡ÊI for E, we can associate to each ui a linear form u . i ¡Ê E., and the u . i have some remarkable properties. 
De.nition 10.2. Given a vector space E and any basis (ui)i¡ÊI for E, by Proposition 2.17, for every i ¡Ê I, there is a unique linear form ui . such that 
 

1 if i = ju . (uj)=
i 0 if i= j, 
for every j ¡Ê I. The linear form ui . is called the coordinate form of index i w.r.t. the basis (ui)i¡ÊI . 
The reason for the terminology coordinate form was explained in Section 2.8. 
We proved in Theorem 2.20 that if (u1,...,un) is a basis of E, then (u1. ,...,u n. ) is a basis of E. called the dual basis. 
If (u1,...,un) is a basis of Rn (more generally Kn), it is possible to .nd explicitly the 
.. .
dual basis (u1,...,u n), where each ui is represented by a row vector. 
Example 10.1. For example, consider the columns of the B¡äezier matrix 
.
. 
B4 = 
... 

1 .33 .1 
03 .63 

00 3 .3 
00 0 1 

...

. 

In other words, we have the basis 

..
..
..
.. 
1 .33 .1 

u1 = 
... 

0 

0 

...

u2 = 
... 

3 

0 

...

u3 = 
... 

.6 

3 

...

u4 = 
... 

3 

.3 

...

. 

00 0 1 
. ....
Since the form u1 is de.ned by the conditions u1(u1)=1,u 1(u2)=0,u 1(u3)=0,u 1(u4) = 0, it is represented by a row vector (¦Ë1 ¦Ë2 ¦Ë3 ¦Ë4) such that
.
.  
 
¦Ë1 ¦Ë2 ¦Ë3 ¦Ë4
... 

1 .33 .1 
03 .63 

00 3 .3 
00 0 1 

... 

 
 
=1 0 0 0. 

This implies that u1 . is the .rst row of the inverse of B4. Since 
.
. 
B.1 
= 
4 
... 

11 11 
01/32/31 

001/31 
00 01 

...

, 

.... .
the linear forms (u1,u 2,u 3,u 4) correspond to the rows of B4 .1 . In particular, u1 is represented by(1 1 1 1). 
The above method works for any n. Given any basis (u1,...,un) of Rn, if P is the n ¡Á n matrix whose jth column is uj, then the dual form u . i is given by the ith row of the matrix P .1 
. When E is of .nite dimension n and (u1,...,un) is a basis of E, by Theorem 10.4 (1), the family (u . 1,...,u . n) is a basis of the dual space E. . Let us see how the coordinates of a ¡Ê E. ..
linear form .. over the dual basis (u1,...,u n) vary under a change of basis. 
Let (u1,...,un) and (v1,...,vn) be two bases of E, and let P =(aij) be the change of basis matrix from (u1,...,un) to (v1,...,vn), so that 
n
nvj = aijui, i=1 
and let P .1 =(bij) be the inverse of P , so that 
n
nui = bjivj. j=1 
For .xed j, where 1 ¡Ü j ¡Ü n, we want to .nd scalars (ci)ni=1 such that 
... . 
v = c1u1 + c2u2 + ¡¤¡¤¡¤ + cnu.
jn
To .nd each ci, we evaluate the above expression at ui. Since ui .(uj)= ¦Äij and vi .(vj)= ¦Äij, we get 
... . 
v (ui)=(c1u ¡¤ )(ui)= ci
j 1 + c2u2 + ¡¤¡¤ + cnun
n
nvj . (ui)= vj . ( bkivk)= bji, k=1 
and thus 
n
nv . = bjiui . .
j 
i=1 
Similar calculations show that 
n
nu . = aijvj . .
i 
j=1 

10.1. THE DUAL SPACE E. AND LINEAR FORMS 347 
This means that the change of basis from the dual basis (u1. ,...,u . n) to the dual basis (v1. ,...,v n.) is (P .1)T . Since 
 
 
nnnnnn. . = .iui . = .i aijvj . = aij.ivj = .iÒ÷ vi . , 
i=1 i=1 j=1 j=1i=1 i=1 
we get 
n
nn
nn
n 
n
.Ò÷j = aij.i, 

i=1 so the new coordinates .Ò÷j are expressed in terms of the old coordinates .i using the matrix P T . If we use the row vectors (.1,...,.n) and (.Ò÷,...,.Ò÷), we have 
1n
(.Ò÷1,...,.nÒ÷)=(.1,...,.n)P. 
These facts are summarized in the following proposition. 
n 
Proposition 10.1. Let (u1,...,un) and (v1,...,vn) be two bases of E, and let P =(aij) be the change of basis matrix from (u1,...,un) to (v1,...,vn), so that 
n
nvj = aijui. i=1 
.. ..
Then the change of basis from the dual basis (u1,...,u n) to the dual basis (v1,...,v n) is (P .1)T, and for any linear form ., the new coordinates .Ò÷j of . are expressed in terms of 
the old coordinates .i of . using the matrix P T; that is, 
(.Ò÷1,...,.nÒ÷)=(.1,...,.n)P. 
To best understand the preceding paragraph, recall Example 3.1, in which E = R2 , u1 = (1, 0), u2 = (0, 1), and v1 = (1, 1), v2 =(.1, 1). Then P , the change of basis matrix from (u1,u2) to (v1,v2), is given by 
 
 

1 .1 
P =,
11with (v1,v2)=(u1,u2)P , and (u1,u2)=(v1,v2)P .1, where 
 
 
1/21/2 
P .1 
=.
.1/21/2
.. ..
Let (u1,u 2) be the dual basis for (u1,u2) and (v1,v 2) be the dual basis for (v1,v2). We claim 
that

 
 
1/2 .1/2
.... .. 
(v1,v 2)=(u1,u 2)=(u1,u 2)(P .1)T ,
1/21/2
.... ..
Indeed, since v = c1u1 + c2u and v = C1u1 + C2u we .nd that 
122 2 
c1 = v1. (u1)= v1. (1/2v1 . 1/2v2)=1/2 
c2 = v1.(u2)= v1.(1/2v1 +1/2v2)=1/2 
C1 = v2. (u1)= v2. (1/2v1 . 1/2v2)= .1/2 
C2 = v2. (u2)= v1. (1/2v1 +1/2v2)=1/2. 

.... ....
Furthermore, since (u1,u 2)=(v1,v 2)P T (since (v1,v 2)=(u1,u 2)(P T).1), we .nd that 
.. .... 
. . = .1u1 + .2u = .1(v . v2)+ .(v1 + v =(.1 + .2)v . 2 +(..11 + .2)v . = .Ò÷ + 2.) Ò÷ . 
1 21v1.2v2 
Hence  
1 .1  1 1  .1 .2  =  .Ò÷ 1 .Ò÷ 2  ,  
where  
P T =  1 .1  1 1  .  
Comparing with the change of basis  

n
nvj = aijui, i=1 
we note that this time, the coordinates (.i) of the linear form .. change in the same direction as the change of basis. For this reason, we say that the coordinates of linear forms are covariant. By abuse of language, it is often said that linear forms are covariant, which explains why the term covector is also used for a linear form. 
Observe that if (e1,...,en) is a basis of the vector space E, then, as a linear map from E to K, every linear form f ¡Ê E. is represented by a 1 ¡Á n matrix, that is, by a row vector 
(¦Ë1 ¡¤¡¤¡¤ ¦Ën), 
with respect to the basis (e1,...,en) of E, and 1 of K, where f(ei)= ¦Ëi. A vector u =
wn 
i=1 uiei ¡Ê E is represented by a n ¡Á 1 matrix, that is, by a column vector 
.
. 

.. 

..
u1 
. 
. ,
. 
un 
and the action of f on u, namely f(u), is represented by the matrix product 

.
. 

.... 
u1 
.
¦Ë1 ¡¤¡¤¡¤ ¦Ën .= ¦Ë1u1 + ¡¤¡¤¡¤ + ¦Ënun.
. 
un 
10.1. THE DUAL SPACE E. AND LINEAR FORMS 349 On the other hand, with respect to the dual basis (e1. ,...,e n. ) of E., the linear form f is 
represented by the column vector

.
. 

.. 

¦Ë1 
. 
. 
.
. 
¦Ën 
..
Remark: In many texts using tensors, vectors are often indexed with lower indices. If so, it is more convenient to write the coordinates of a vector x over the basis (u1,...,un) as (xi), using an upper index, so that 
n
i 
x = xui, 
i=1 
and in a change of basis, we have 
n 
ni 
vj = ajui 
n i=1 
and 
n
n
x i = aji xÒ÷j. 

j=1 
Dually, linear forms are indexed with upper indices. Then it is more convenient to write the coordinates of a covector .. over the dual basis (u .1 ,...,u .n) as (.i), using a lower index, so that 
n
n
. . = .iu .i 
i=1 

and in a change of basis, we have 
n
n.ii .j
u = ajv j=1 
and 
n
n
.Ò÷ = a i .i.
jji=1 
With these conventions, the index of summation appears once in upper position and once in lower position, and the summation sign can be safely omitted, a trick due to Einstein. For example, we can write .Ò÷ j = a ij.i as an abbreviation for 
n
n
.Ò÷ = a i .i.
jji=1 
For another example of the use of Einstein¡¯s notation, if the vectors (v1,...,vn) are linear combinations of the vectors (u1,...,un), with 
n
n 
vi = aijuj, 1 ¡Ü i ¡Ü n, j=1 
then the above equations are written as vi = aijuj, 1 ¡Ü i ¡Ü n. Thus, in Einstein¡¯s notation, the n ¡Á n matrix (aij) is denoted by (aij), a (1, 1)-tensor. Beware that some authors view a matrix as a mapping between coordinates, in which case the matrix (aij) is denoted by (aji ). 


10.2 	Pairing and Duality Between E and E. ¡Ê E..
Given a linear form u . and a vector v ¡Ê E, the result u .(v) of applying u to v is also denoted by£¨u . ,v£©. This de.nes a binary operation£¨.,.£©: E. ¡Á E ¡ú K satisfying the following properties:
£¨u1 . + u2. ,v£© =£¨u1. ,v£© +£¨u2. ,v£©
£¨u . ,v1 + v2£© =£¨u . ,v1£© +£¨u . ,v2£©£¨¦Ëu . ,v£© = ¦Ë£¨u . ,v£©£¨u . , ¦Ëv£© = ¦Ë£¨u . ,v£©. 
The above identities mean that£¨.,.£© is a bilinear map, since it is linear in each argument. It is often called the canonical pairing between E. and E. In view of the above identities, given any .xed vector v ¡Ê E, the map evalv : E. ¡ú K (evaluation at v) de.ned such that 
evalv(u . )=£¨u . ,v£© = u . (v) for every u . ¡Ê E. 
is a linear map from E. to K, that is, evalv is a linear form in E.. . Again, from the above identities, the map evalE : E ¡ú E.., de.ned such that 
evalE(v) = evalv for every v ¡Ê E, 
is a linear map. Observe that 
evalE(v)(u . ) = evalv(u . )=£¨u . ,v£© = u . (v), for all v ¡Ê E and all u . ¡Ê E. . 
We shall see that the map evalE is injective, and that it is an isomorphism when E has .nite dimension. 
10.2. PAIRING AND DUALITY BETWEEN E AND E. 
We now formalize the notion of the set V 0 of linear equations vanishing on all vectors in a given subspace V . E, and the notion of the set U0 of common solutions of a given set U . E. of linear equations. The duality theorem (Theorem 10.4) shows that the dimensions of V and V 0, and the dimensions of U and U0, are related in a crucial way. It also shows that, in .nite dimension, the maps V ¡ú V 0 and U ¡ú U0 are inverse bijections from subspaces of E to subspaces of E. . 
De.nition 10.3. Given a vector space E and its dual E., we say that a vector v ¡Ê E and a linear form u . ¡Ê E. are orthogonal i.£¨u . ,v£© = 0. Given a subspace V of E and a subspace U of E., we say that V and U are orthogonal i.£¨u . ,v£© = 0 for every u . ¡Ê U and every v ¡Ê V . Given a subset V of E (resp. a subset U of E.), the orthogonal V 0 of V is the subspace V 0 of E. de.ned such that 
V 0 = {u . ¡Ê E .|£¨ u . ,v£© =0, for every v ¡Ê V } 
(resp. the orthogonal U0 of U is the subspace U0 of E de.ned such that 
U0 = {v ¡Ê E|£¨ u . ,v£© =0, for every u . ¡Ê U}). 
The subspace V 0 . E. is also called the annihilator of V . The subspace U0 . E annihilated by U . E. does not have a special name. It seems reasonable to call it the linear subspace (or linear variety) de.ned by U. 
Informally, V 0 is the set of linear equations that vanish on V , and U0 is the set of common zeros of all linear equations in U. We can also de.ne V 0 by 
V 0 = {u . ¡Ê E . | V . Ker u .} 
and U0 by 
 
U0 . 
=Ker u. 
u .¡ÊU 
Observe that E0 = {0} = (0), and {0}0 = E. . 
Proposition 10.2. If V1 . V2 . E, then V20 . V10 . E., and if U1 . U2 . E., then U20 . U10 . E. See Figure 10.2. 
Proof. Indeed, if V1 . V2 . E, then for any f. ¡Ê V 0 we have f.(v)=0 for all v ¡Ê V2, and thus f.(v)=0 for all v ¡Ê V1, so f. ¡Ê V10 . Similarly,2if U1 . U2 . E., then for any v ¡Ê U20 , we have f.(v) = 0 for all f. ¡Ê U2, so f.(v) = 0 for all f. ¡Ê U1, which means that v ¡Ê U10 . 
Here are some examples. 
Example 10.2. Let E =M2(R), the space of real 2 ¡Á 2 matrices, and let V be the subspace of M2(R) spanned by the matrices 
01 10 00 
,,. 
10 00 01 

Figure 10.2: The top pair of .gures schematically illustrates the relation if V1 . V2 . E, then V20 . V10 . E., while the bottom pair of .gures illustrates the relationship if U1 . U2 . E. , then U20 . U10 . E. 
We check immediately that the subspace V consists of all matrices of the form 
b  a  
,  
a  c  

that is, all symmetric matrices. The matrices 
a11 a12 a21 a22 
in V satisfy the equation a12 . a21 =0, 
and all scalar multiples of these equations, so V 0 is the subspace of E. spanned by the linear form given by u .(a11,a12,a21,a22)= a12 . a21. By the duality theorem (Theorem 10.4) we have 
dim(V 0) = dim(E) . dim(V )=4 . 3=1. 
Example 10.3. The above example generalizes to E =Mn(R) for any n ¡Ý 1, but this time, consider the space U of linear forms asserting that a matrix A is symmetric; these are the linear forms spanned by the n(n . 1)/2 equations 
aij . aji =0, 1 ¡Ü i<j ¡Ü n; 
Note there are no constraints on diagonal entries, and half of the equations 
aij . aji =0, 1 ¡Ü i = j ¡Ü n 

10.2. PAIRING AND DUALITY BETWEEN E AND E. 
are redundant. It is easy to check that the equations (linear forms) for which i<j are linearly independent. To be more precise, let U be the space of linear forms in E. spanned by the linear forms 
uij. (a11,...,a1n,a21,...,a2n,...,an1,...,ann)= aij . aji, 1 ¡Ü i<j ¡Ü n. 
The dimension of U is n(n . 1)/2. Then the set U0 of common solutions of these equations is the space S(n) of symmetric matrices. By the duality theorem (Theorem 10.4), this space has dimension 
n(n + 1) n(n . 1) 
= n 2 . . 
22 We leave it as an exercise to .nd a basis of S(n). 
Example 10.4. If E =Mn(R), consider the subspace U of linear forms in E. spanned by the linear forms 
u . (a11,...,a1n,a21,...,a2n,...,an1,...,ann)= aij + aji, 1 ¡Ü i<j ¡Ü n
iju . ii(a11,...,a1n,a21,...,a2n,...,an1,...,ann)= aii, 1 ¡Ü i ¡Ü n. 
It is easy to see that these linear forms are linearly independent, so dim(U)= n(n + 1)/2. The space U0 of matrices A ¡Ê Mn(R) satifying all of the above equations is clearly the space Skew(n) of skew-symmetric matrices. By the duality theorem (Theorem 10.4), the dimension of U0 is 
n(n . 1) n(n + 1) 
= n 2 . . 
22 We leave it as an exercise to .nd a basis of Skew(n). 
Example 10.5. For yet another example with E =Mn(R), for any A ¡Ê Mn(R), consider the linear form in E. given by 
tr(A)= a11 + a22 + ¡¤¡¤¡¤ + ann, 
called the trace of A. The subspace U0 of E consisting of all matrices A such that tr(A)=0 is a space of dimension n2 . 1. We leave it as an exercise to .nd a basis of this space. 
The dimension equations 
dim(V ) + dim(V 0) = dim(E) dim(U) + dim(U0) = dim(E) 
are always true (if E is .nite-dimensional). This is part of the duality theorem (Theorem 10.4). 
Remark: In contrast with the previous examples, given a matrix A ¡Ê Mn(R), the equations asserting that ATA = I are not linear constraints. For example, for n = 2, we have 
22 
a11 + a21 =1 a 2 222 =1 
21 + a 
a11a12 + a21a22 =0. 

Remarks: 
(1) 
The notation V 0 (resp. U0) for the orthogonal of a subspace V of E (resp. a subspace U of E.) is not universal. Other authors use the notation V ¡Í (resp. U¡Í). However, the notation V ¡Í is also used to denote the orthogonal complement of a subspace V with respect to an inner product on a space E, in which case V ¡Í is a subspace of E and not a subspace of E. (see Chapter 11). To avoid confusion, we prefer using the notation V 0 . 

(2) 
Since linear forms can be viewed as linear equations (at least in .nite dimension), given a subspace (or even a subset) U of E., we can de.ne the set Z(U) of common zeros of the equations in U by 


Z(U)= {v ¡Ê E | u . (v)=0, for all u . ¡Ê U}. 
Of course Z(U)= U0, but the notion Z(U) can be generalized to more general kinds of equations, namely polynomial equations. In this more general setting, U is a set of polynomials in n variables with coe.cients in a .eld K (where n = dim(E)). Sets of the form Z(U) are called algebraic varieties. Linear forms correspond to the special case where homogeneous polynomials of degree 1 are considered. 
If V is a subset of E, it is natural to associate with V the set of polynomials in K[X1,...,Xn] that vanish on V . This set, usually denoted I(V ), has some special properties that make it an ideal. If V is a linear subspace of E, it is natural to restrict our attention to the space V 0 of linear forms that vanish on V , and in this case we identify I(V ) and V 0 (although technically, I(V ) is no longer an ideal). 
For any arbitrary set of polynomials U . K[X1,...,Xn] (resp. subset V . E), the relationship between I(Z(U)) and U (resp. Z(I(V )) and V ) is generally not simple, even though we always have 
U .I(Z(U)) (resp. V .Z(I(V ))). 
However, when the .eld K is algebraically closed, then I(Z(U)) is equal to the radical of the ideal U, a famous result due to Hilbert known as the Nullstellensatz (see Lang 
[41] or Dummit and Foote [19]). The study of algebraic varieties is the main subject of algebraic geometry, a beautiful but formidable subject. For a taste of algebraic geometry, see Lang [41] or Dummit and Foote [19]. 

10.3. THE DUALITY THEOREM AND SOME CONSEQUENCES 
The duality theorem (Theorem 10.4) shows that the situation is much simpler if we restrict our attention to linear subspaces; in this case 
U = I(Z(U)) and V = Z(I(V )). 
. V 00
Proposition 10.3. We have V for every subspace V of E, and U . U00 for every subspace U of E. . 
Proof. Indeed, for any v ¡Ê V , to show that v ¡Ê V 00 we need to prove that u .(v)=0 for all u . ¡Ê V 0 . However, V 0 consists of all linear forms u . such that u .(y) = 0 for all y ¡Ê V ; in particular, for a .xed v ¡Ê V , we have u .(v)=0 for all u . ¡Ê V 0, as required. 
.. ¡Ê U00
Similarly, for any u ¡Ê U, to show that u we need to prove that u .(v) =0 for all v ¡Ê U0 . However, U0 consists of all vectors v such that f.(v) = 0 for all f. ¡Ê U; in particular, for a .xed u . ¡Ê U, we have u .(v)=0 for all v ¡Ê U0, as required. 
= U00
We will see shortly that in .nite dimension, we have V = V 00 and U . 


10.3 The Duality Theorem and Some Consequences 
Given a vector space E of dimension n ¡Ý 1 and a subspace U of E, by Theorem 2.13, every basis (u1,...,um) of U can be extended to a basis (u1,...,un) of E. We have the following important theorem adapted from E. Artin [2] (Chapter 1). 
Theorem 10.4. (Duality theorem) Let E be a vector space of dimension n. The following properties hold: 
(a) 
For every basis (u1,...,un) of E, the family of coordinate forms (u . 1,...,u . n) is a basis of E. (called the dual basis of (u1,...,un)). 

(b) 
For every subspace V of E, we have V 00 = V . 

(c) 
For every pair of subspaces V and W of E such that E = V ¨’ W , with V of dimen-sion m, for every basis (u1,...,un) of E such that (u1,...,um) is a basis of V and (um+1,...,un) is a basis of W , the family (u1. ,...,u . ) is a basis of the orthogonal W 0 


m
of W in E., so that 
dim(W ) + dim(W 0) = dim(E). 

Furthermore, we have W 00 = W . 
(d) For every subspace U of E., we have 
dim(U) + dim(U0) = dim(E), 
where U0 is the orthogonal of U in E, and U00 = U. 
Proof. (a) This part was proven in Theorem 2.20. 
(b) By Proposition 10.3 we have V . V 00 . If V = V 00, then let (u1,...,up) be a basis of V 00 such that (u1,...,um) is a basis of V , with m<p. Since um+1 ¡Ê V 00 , um+1 is orthogonal to every linear form in V 0 . By de.nition we have um. +1(ui) = 0 for all i =1,...,m, and thus u . ¡Ê V 0 . However, um. +1(um+1) = 1, contradicting the fact that um+1 is orthogonal 
m+1 = V 00
to every linear form in V 0 . Thus, V . 
(c) Every linear form f. ¡Ê W 0 is orthogonal to every uj for j = m +1,...,n, and thus, f.(uj)=0 for j = m +1,...,n. For such a linear form f. ¡Ê W 0, let 
.. . 
g = f . (u1)u1 + ¡¤¡¤¡¤ + f . (um)um. 
We have g .(ui)= f.(ui), for every i,1 ¡Ü i ¡Ü m. Furthermore, by de.nition, g . vanishes on 
all uj with j = m+1,...,n. Thus, f. and g . agree on the basis (u1,...,un) of E, and so g . = f. .. 
. This shows that (u1,...,u m) generates W 0, and since it is also a linearly independent family, (u1. ,...,u m. ) is a basis of W 0 . It is then obvious that dim(W ) + dim(W 0) = dim(E), and by Part (b), we have W 00 = W . 
(d) The only remaining fact to prove is that U00 = U. 1 ,...,f. ) be a basis of U.
Let (f. mNote that the map h: E ¡ú Km de.ned such that 
h(v)=(f1 . m
(v),...,f . (v)) 
for every v ¡Ê E is a linear map, and that its kernel Ker h is precisely U0 . Then by Proposition 5.8, 
n = dim(E) = dim(Ker h) + dim(Im h) ¡Ü dim(U0)+ m, 
since dim(Im h) ¡Ü m. Thus, n . dim(U0) ¡Ü m. By (c), we have dim(U0) + dim(U00)= dim(E)= n, so we get dim(U00) ¡Ü m. However, by Proposition 10.3 it is clear that U . U00 , which implies m = dim(U) ¡Ü dim(U00), so dim(U) = dim(U00)= m, and we must have 
= U00
U . 
Part (a) of Theorem 10.4 shows that 
dim(E) = dim(E . ), 
and if (u1,...,un) is a basis of E, then (u . 1,...,u . n) is a basis of the dual space E. called the dual basis of (u1,...,un). 
De.ne the function E (E for equations) from subspaces of E to subspaces of E. and the function Z (Z for zeros) from subspaces of E. to subspaces of E by 
E(V )= V 0 ,V . E 
Z(U)= U0 ,U . E . . 
By Parts (c) and (d) of Theorem 10.4, 
(Z.E)(V )= V 00 = V 
(E.Z)(U)= U00 = U, 
10.3. THE DUALITY THEOREM AND SOME CONSEQUENCES 
so Z.E = id and E.Z = id, and the maps E and Z are inverse bijections. These maps set up a duality between subspaces of E and subspaces of E. . In particular, every subspace V . E of dimension m is the set of common zeros of the space of linear forms (equations) V 0, which has dimension n . m. This con.rms the claim we made about the dimension of the subpsace de.ned by a set of linear equations. 
One should be careful that this bijection does not hold if E has in.nite dimension. Some restrictions on the dimensions of U and V are needed. 
Remark: However, even if E is in.nite-dimensional, the identity V = V 00 holds for every subspace V of E. The proof is basically the same but uses an in.nite basis of V 00 extending a basis of V . 
We now discuss some applications of the duality theorem. 
Problem 1 . Suppose that V is a subspace of Rn of dimension m and that (v1,...,vm) is a basis of V . The problem is to .nd a basis of V 0 . 
We .rst extend (v1,...,vm) to a basis (v1,...,vn) of Rn, and then by part (c) of Theorem 
..
10.4, we know that (vm+1,...,v n) is a basis of V 0 . 
Example 10.6. For example, suppose that V is the subspace of R4 spanned by the two linearly independent vectors 
..
.. 
11 

v1 = 
... 

1 

1 

...

v2 = 
... 

1 

.1 

...

, 

1 .1 the .rst two vectors of the Haar basis in R4 . The four columns of the Haar matrix 
.
. 
W = 

... 

11 1 0 
11 .10 

1 .10 1 
1 .10 .1 

... 

form a basis of R4, and the inverse of W is given by 
.
..
..
. 
1/4000 1111 1/41/41/41/4 

W .1 
= 
... 

01/40 0 

0 01/20 

... 
... 

... 

= 

... 

...

. 

11 

.1 .1 

1/41/4 

.1/4 .1/4 

1 

.1 

00 

1/2 

.1/2 

00 

0 0 01/2 001 .1 001/2 .1/2 

....
Since the dual basis (v1,v 2,v 3,v 4) is given by the rows of W .1, the last two rows of W .1 , 1/2 .1/20 0 
,
0 01/2 .1/2 
form a basis of V 0 . We also obtain a basis by rescaling by the factor 1/2, so the linear forms 
given by the row vectors 
1 .10 0 
001 .1 
form a basis of V 0, the space of linear forms (linear equations) that vanish on the subspace 
V . 
The method that we described to .nd V 0 requires .rst extending a basis of V and then inverting a matrix, but there is a more direct method. Indeed, let A be the n ¡Á m matrix whose columns are the basis vectors (v1,...,vm) of V . Then a linear form u represented by a row vector belongs to V 0 i. uvi = 0 for i =1,...,m i. 
uA =0 
i. 
AT 
u T =0. 
Therefore, all we need to do is to .nd a basis of the nullspace of AT . This can be done quite e.ectively using the reduction of a matrix to reduced row echelon form (rref); see Section 
7.10. 
Example 10.7. For example, if we reconsider the previous example, ATuT = 0 becomes 
.
. 

u1 
111 1 

11 .1 .1 

... 

u2 
u3 u4 
... 

0 

= . 

0 

Since the rref of AT is  1 0  1 0  0 1  0 1  ,  
the above system is equivalent to  

.
. 

u1 
1100 

0011 

... 

u2 
u3 u4 
... 

u1 + u2 0 
== , 

u3 + u4 0 
where the free variables are associated with u2 and u4. Thus to determine a basis for the kernel of AT, we set u2 =1,u4 = 0 and u2 =0,u4 = 1 and obtain a basis for V 0 as 
1 .100 , 001 .1 . 

10.3. THE DUALITY THEOREM AND SOME CONSEQUENCES 
Problem 2 . Let us now consider the problem of .nding a basis of the hyperplane H in Rn de.ned by the equation c1x1 + ¡¤¡¤¡¤ + cnxn =0. 
More precisely, if u .(x1,...,xn) is the linear form in (Rn). given by u .(x1,...,xn)= c1x1 + ¡¤¡¤¡¤ + cnxn, then the hyperplane H is the kernel of u . . Of course we assume that some cj is nonzero, in which case the linear form u . spans a one-dimensional subspace U of (Rn)., and U0 = H has dimension n . 1. 
Since u . is not the linear form which is identically zero, there is a smallest positive index j ¡Ü n such that cj = 0, so our linear form is really u .(x1,...,xn)= cjxj + ¡¤¡¤¡¤ + cnxn. We claim that the following n . 1 vectors (in Rn) form a basis of H: 
12 ... j . 1 jj +1 ... n . 1 
.
. 
1 10 ... 00 0 ... 0 

2 

. 

. 

. 

j . 1 

j 

j +1 

j +2 

. 

. 

.............. 

.............. 

. 

01 ... 00 0 ... 0 

..... .
..
.. 

... . 

..

. 

.

.. 

.. . 

. 

00 ... 10 0 ... 0 


00 ... 0 



.cj+1/cj .cj+2/cj ... .cn/cj 
00 ... 01 0 ... 0 


00 ... 00 1 ... 0 



..... .
..
.. 

... . 

..

. 

.

. .. 

.. . 

. 

n 00 ... 00 0 ... 1 
Observe that the (n.1)¡Á(n.1) matrix obtained by deleting row j is the identity matrix, so the columns of the above matrix are linearly independent. A simple calculation also shows that the linear form u .(x1,...,xn)= cjxj + ¡¤¡¤¡¤ + cnxn vanishes on every column of the above matrix. For a concrete example in R6, if u .(x1,...,x6)= x3 +2x4 +3x5 +4x6, we obtain the basis for the hyperplane H of equation 
x3 +2x4 +3x5 +4x6 =0 
given by the following matrix: 
.
. 
....... 

1  0  0  0  0  
0  1  0  0  0  
0  0  .2  .3  .4  

0  0  1  0  0  
0  0  0  1  0  
0  0  0  0  1  

....... 

. 

Problem 3 . Conversely, given a hyperplane H in Rn given as the span of n . 1 linearly vectors (u1,...,un.1), it is possible using determinants to .nd a linear form (¦Ë1,...,¦Ën) that vanishes on H. 
In the case n = 3, we are looking for a row vector (¦Ë1,¦Ë2,¦Ë3) such that if 
.. .. 
u1 v1 u = .u2. and v = .v2. u3 v3 
are two linearly independent vectors, then 
.. 
¦Ë1 
u1 u2 u2 0
..
¦Ë2 = , 
v1 v2 v2 0 
¦Ë3 
and the cross-product u ¡Á v of u and v given by 
.. 
u2v3 . u3v2 
..
u ¡Á v = 	u3v1 . u1v3 u1v2 . u2v1 
is a solution. In other words, the equation of the plane spanned by u and v is 
(u2v3 . u3v2)x +(u3v1 . u1v3)y +(u1v2 . u2v1)z =0. 
Problem 4 . Here is another example illustrating the power of Theorem 10.4. Let E =Mn(R), and consider the equations asserting that the sum of the entries in every row of a matrix A ¡Ê Mn(R) is equal to the same number. We have n . 1 equations 
n
n 
(aij . ai+1j)=0, 1 ¡Ü i ¡Ü n . 1, 
j=1 
and it is easy to see that they are linearly independent. Therefore, the space U of linear forms in E. spanned by the above linear forms (equations) has dimension n . 1, and the space U0 of matrices satisfying all these equations has dimension n2 . n + 1. It is not so obvious to .nd a basis for this space. 
We will now pin down the relationship between a vector space E and its bidual E.. . 


10.4 	The Bidual and Canonical Pairings 
Proposition 10.5. Let E be a vector space. The following properties hold: 
(a) The linear map evalE : E ¡ú E.. de.ned such that 
evalE(v) = evalv for all v ¡Ê E, 
that is, evalE(v)(u .)=£¨u . ,v£© = u .(v) for every u . ¡Ê E., is injective. 
10.4. THE BIDUAL AND CANONICAL PAIRINGS 
(b) When E is of .nite dimension n, the linear map evalE : E ¡ú E.. is an isomorphism (called the canonical isomorphism). 
w 
Proof. (a) Let (ui)i¡ÊI be a basis of E, and let v = i¡ÊI viui. If evalE(v) = 0, then in particular evalE(v)(ui .)=0 for all ui ., and since 
evalE(v)(ui . )=£¨ui.,v£© = vi, 
we have vi = 0 for all i ¡Ê I, that is, v = 0, showing that evalE : E ¡ú E.. is injective. If E is of .nite dimension n, by Theorem 10.4, for every basis (u1,...,un), the family 
.. ....
(u1,...,u n) is a basis of the dual space E., and thus the family (u ,...,u n ) is a basis of the bidual E.. . Thisshowsthatdim(E)=dim(E..)=n,andsinceb1y Part (a), we know that evalE : E ¡ú E.. is injective, in fact, evalE : E ¡ú E.. is bijective (by Proposition 5.11). 
When E is of .nite dimension and (u1,...,un) is a basis of E, in view of the canon-
.. ..
ical isomorphism evalE : E ¡ú E.., the basis (u1 ,...,u n ) of the bidual is identi.ed with (u1,...,un). 
Proposition 10.5 can be reformulated very fruitfully in terms of pairings, a remarkably useful concept discovered by Pontrjagin in 1931 (adapted from E. Artin [2], Chapter 1). Given two vector spaces E and F over a .eld K, we say that a function .: E ¡Á F ¡ú K is bilinear if for every v ¡Ê V , the map u ¡ú .(u, v) (from E to K) is linear, and for every u ¡Ê E, the map v ¡ú .(u, v) (from F to K) is linear. 
De.nition 10.4. Given two vector spaces E and F over K,a pairing between E and F is a bilinear map .: E ¡Á F ¡ú K. Such a pairing is nondegenerate i. 
(1) for every u ¡Ê E, if .(u, v)=0 for all v ¡Ê F , then u = 0, and 

(2) for every v ¡Ê F , if .(u, v)=0 for all u ¡Ê E, then v = 0. 


A pairing .: E ¡Á F ¡ú K is often denoted by£¨.,.£©: E ¡Á F ¡ú K. For example, the map£¨.,.£©: E. ¡Á E ¡ú K de.ned earlier is a nondegenerate pairing (use the proof of (a) in Proposition 10.5). If E = F and K = R, any inner product on E is a nondegenerate pairing (because an inner product is positive de.nite); see Chapter 11. Other interesting nondegenerate pairings arise in exterior algebra and di.erential geometry. 
Given a pairing . : E ¡Á F ¡ú K, we can de.ne two maps l. : E ¡ú F . and r. : F ¡ú E. as follows: For every u ¡Ê E, we de.ne the linear form l.(u) in F . such that 
l.(u)(y)= .(u, y) for every y ¡Ê F, 
and for every v ¡Ê F , we de.ne the linear form r.(v) in E. such that 
r.(v)(x)= .(x, v) for every x ¡Ê E. 
We have the following useful proposition. 
Proposition 10.6. Given two vector spaces E and F over K, for every nondegenerate pairing .: E ¡Á F ¡ú K between E and F , the maps l. : E ¡ú F . and r. : F ¡ú E. are linear and injective. Furthermore, if E and F have .nite dimension, then this dimension is the same and l. : E ¡ú F . and r. : F ¡ú E. are bijections. 
Proof. The maps l. : E ¡ú F . and r. : F ¡ú E. are linear because a pairing is bilinear. If l.(u) = 0 (the null form), then 
l.(u)(v)= .(u, v) = 0 for every v ¡Ê F, 
and since . is nondegenerate, u = 0. Thus, l. : E ¡ú F . is injective. Similarly, r. : F ¡ú E. is injective. When F has .nite dimension n, we have seen that F and F . have the same dimension. Since l. : E ¡ú F . is injective, we have m = dim(E) ¡Ü dim(F )= n. The same argument applies to E, and thus n = dim(F ) ¡Ü dim(E)= m. But then, dim(E) = dim(F ), and l. : E ¡ú F . and r. : F ¡ú E. are bijections. 
When E has .nite dimension, the nondegenerate pairing£¨.,.£©: E. ¡Á E ¡ú K yields another proof of the existence of a natural isomorphism between E and E.. . When E = F , the nondegenerate pairing induced by an inner product on E yields a natural isomorphism between E and E. (see Section 11.2). 
We now show the relationship between hyperplanes and linear forms. 


10.5 Hyperplanes and Linear Forms 
Actually Proposition 10.7 below follows from Parts (c) and (d) of Theorem 10.4, but we feel that it is also interesting to give a more direct proof. 
Proposition 10.7. Let E be a vector space. The following properties hold: 
(a) Given any nonnull linear form f. ¡Ê E., its kernel H = Ker f. is a hyperplane. 

(b) For any hyperplane H in E, there is a (nonnull) linear form f. ¡Ê E. such that H = Ker f. . 

(c) Given any hyperplane H in E and any (nonnull) linear form f. ¡Ê E. such that H = 


. ..
Ker f., for every linear form g ¡Ê E. , H = Ker g i. g = ¦Ëf. for some ¦Ë =0 in K. 
Proof. (a) If f. ¡Ê E. is nonnull, there is some vector v0 ¡Ê E such that f.(v0) = 0. Let H = Ker f. . For every v ¡Ê E, we have 
f.(v) f.(v)
f . v . v0 = f . (v) . f . (v0)= f . (v) . f . (v)=0. 
f.(v0) f.(v0) 
Thus, 
f.(v) 
v . v0 = h ¡Ê H, 
f.(v0) 
10.6. TRANSPOSE OF A LINEAR MAP AND OF A MATRIX 
and 
f.(v) 
v = h + v0,
f.(v0) 
that is, E = H + Kv0. Also since f.(v0) = 0, we have v0 ¡Ê/H, that is, H ¡É Kv0 = 0. Thus, E = H ¨’ Kv0, and H is a hyperplane. 
(b) 
If H is a hyperplane, E = H ¨’ Kv0 for some v0 ¡Ê/H. Then every v ¡Ê E can be written in a unique way as v = h + ¦Ëv0. Thus there is a well-de.ned function f. : E ¡ú K, such that, f.(v)= ¦Ë, for every v = h + ¦Ëv0. We leave as a simple exercise the veri.cation that f. is a linear form. Since f.(v0) = 1, the linear form f. is nonnull. Also, by de.nition, it is clear that ¦Ë = 0 i. v ¡Ê H, that is, Ker f. = H. 

(c) 
Let H be a hyperplane in E, and let f. ¡Ê E. be any (nonnull) linear form such that H = Ker f. . Clearly, if g . = ¦Ëf. for some ¦Ë = 0, then H = Ker g . . Conversely, assume that H = Ker g . for some nonnull linear form g . . From (a), we have E = H ¨’ Kv0, for some v0 such that f.(v0)=0 and g .(v0) = 0. Then observe that 


g .(v0)
. . f . 
g 
f.(v0) 
is a linear form that vanishes on H, since both f. and g . vanish on H, but also vanishes on Kv0. Thus, g . = ¦Ëf., with 
g .(v0)
¦Ë = . 
f.(v0) 
We leave as an exercise the fact that every subspace V = E of a vector space E is the intersection of all hyperplanes that contain V . We now consider the notion of transpose of a linear map and of a matrix. 


10.6 Transpose of a Linear Map and of a Matrix 
Given a linear map f : E ¡ú F , it is possible to de.ne a map fT : F . ¡ú E. which has some interesting properties. 
De.nition 10.5. Given a linear map f : E ¡ú F , the transpose fT : F . ¡ú E. of f is the linear map de.ned such that 
fT(v . )= v . . f, for every v . ¡Ê F . , 
as shown in the diagram below: 
f 
E ÈÏÈÏ. 
F 
ÈÏ
. 
f (v .)ÈÏÈÏÈÏÈÏÈÏ  l v 
K. 
Equivalently, the linear map fT : F . ¡ú E. is de.ned such that£¨v . ,f(u)£© =£¨fT(v . ),u£©, (.) for all u ¡Ê E and all v . ¡Ê F . . It is easy to verify that the following properties hold: (f + g)T = fT + g T (g . f)T = fT . g T idT = idE. .
E 
Note the reversal of composition on the right-hand side of (g . f)T = fT . gT . 
The equation (g . f)T = fT . gT implies the following useful proposition. Proposition 10.8. If f : E ¡ú F is any linear map, then the following properties hold: 
(1) If f is injective, then fT is surjective. 

(2) If f is surjective, then fT is injective. 


Proof. If f : E ¡ú F is injective, then it has a retraction r : F ¡ú E such that r . f = idE, and if f : E ¡ú F is surjective, then it has a section s: F ¡ú E such that f . s = idF . Now if 
f : E ¡ú F is injective, then we have (r . f)T = fT . r T = idE. , which implies that fT is surjective, and if f is surjective, then we have (f . s)T = s T . fT = idF . , which implies that fT is injective. 
The following proposition shows the relationship between orthogonality and transposi-tion. 
Proposition 10.9. Given a linear map f : E ¡ú F , for any subspace V of E, we have 
f(V )0 =(fT).1(V 0)= {w . ¡Ê F . | fT(w . ) ¡Ê V 0}. As a consequence, 
Ker fT = (Im f)0 . 
We also have 
Ker f = (Im fT)0 . 
10.6. TRANSPOSE OF A LINEAR MAP AND OF A MATRIX 
Proof. We have£¨w . ,f(v)£© =£¨fT(w . ),v£©, 
for all v ¡Ê E and all w . ¡Ê F . , and thus, we have£¨w .,f(v)£© = 0 for every v ¡Ê V , i.e. w . ¡Ê f(V )0 i.£¨fT(w .),v£© = 0 for every v ¡Ê V i. fT(w .) ¡Ê V 0, i.e. w . ¡Ê (fT).1(V 0), proving that 
f(V )0 =(fT).1(V 0). 
Since we already observed that E0 = (0), letting V = E in the above identity we obtain that Ker fT = (Im f)0 . 
From the equation
£¨w . ,f(v)£© =£¨fT(w . ),v£©, 

we deduce that v ¡Ê (Im fT)0 i.£¨fT(w .),v£© = 0 for all w . ¡Ê F . i.£¨w .,f(v)£© = 0 for all w . ¡Ê F . . Assume that v ¡Ê (Im fT)0 . If we pick a basis (wi)i¡ÊI of F , then we have the linear forms wi . : F ¡ú K such that w .(wj)= ¦Äij, and since we must have£¨w .,f(v)£© = 0 for all 
ii i ¡Ê I and (wi)i¡ÊI is a basis of F , we conclude that f(v) = 0, and thus v ¡Ê Ker f (this is because£¨w .,f(v)£© is the coe.cient of f(v) associated with the basis vector wi). Conversely, 
i 
if v ¡Ê Ker f, then£¨w .,f(v)£© = 0 for all w . ¡Ê F . , so we conclude that v ¡Ê (Im fT)0 . Therefore, v ¡Ê (Im fT)0 i. v ¡Ê Ker f; that is, 
Ker f = (Im fT)0 , 
as claimed. 
The following theorem shows the relationship between the rank of f and the rank of fT . 
Theorem 10.10. Given a linear map f : E ¡ú F , the following properties hold. 
(a) The dual (Im f). of Im f is isomorphic to Im fT = fT(F .); that is, 
¡«
(Im f) . = Im fT . 
(b) If F is .nite dimensional, then rk(f) = rk(fT). 
Proof. (a) Consider the linear maps 
pj
E .¡ú Im f .¡ú F, 
pfj
where E .¡ú Im f is the surjective map induced by E .¡ú F , and Im f .¡ú F is the injective inclusion map of Im f into F . By de.nition, f = j . p. To simplify the notation, 
pp
let I = Im f. By Proposition 10.8, since E .¡ú I is surjective, I. .¡ú E. is injective, and 
jj
since Im f .¡ú F is injective, F . .¡ú I. is surjective. Since f = j . p, we also have 
fT =(j . p)T = p T . jT , 
jp
and since F . .¡ú I. is surjective, and I. .¡ú E. is injective, we have an isomorphism between (Im f). and fT(F .). 
(b) We already noted that Part (a) of Theorem 10.4 shows that dim(F ) = dim(F .), for every vector space F of .nite dimension. Consequently, dim(Im f) = dim((Im f).), and thus, by Part (a) we have rk(f) = rk(fT). 
Remark: When both E and F are .nite-dimensional, there is also a simple proof of (b) that doesn¡¯t use the result of Part (a). By Theorem 10.4(c) 
dim(Im f) + dim((Im f)0) = dim(F ), 
and by Theorem 5.8 dim(Ker fT) + dim(Im fT) = dim(F . ). 
Furthermore, by Proposition 10.9, we have 
Ker fT = (Im f)0 , 
and since F is .nite-dimensional dim(F ) = dim(F .), so we deduce 
dim(Im f) + dim((Im f)0) = dim((Im f)0) + dim(Im fT), 
which yields dim(Im f) = dim(Im fT); that is, rk(f) = rk(fT). 
The following proposition can be shown, but it requires a generalization of the duality theorem, so its proof is omitted. 
Proposition 10.11. If f : E ¡ú F is any linear map, then the following identities hold: 
Im fT = (Ker (f))0 Ker (fT) = (Im f)0 Im f = (Ker (fT)0 Ker (f) = (Im fT)0 . 
Observe that the second and the fourth equation have already be proven in Proposition 
10.9. Since for any subspace V . E, even in.nite-dimensional, we have V 00 = V , the third equation follows from the second equation by taking orthogonals. Actually, the fourth equation follows from the .rst also by taking orthogonals. Thus the only equation to be proven is the .rst equation. We will give a proof later in the case where E is .nite-dimensional (see Proposition 10.18). 
The following proposition shows the relationship between the matrix representing a linear map f : E ¡ú F and the matrix representing its transpose fT : F . ¡ú E. . 
10.6. 	TRANSPOSE OF A LINEAR MAP AND OF A MATRIX 
Proposition 10.12. Let E and F be two vector spaces, and let (u1,...,un) be a basis for E and (v1,...,vm) be a basis for F . Given any linear map f : E ¡ú F , if M(f) is the m ¡Á n-matrix representing f w.r.t. the bases (u1,...,un) and (v1,...,vm), then the n ¡Á m-matrix 
. 	.. ..
M(fT) representing fT : F ¡ú E. w.r.t. the dual bases (v1,...,v m) and (u1,...,u n) is the transpose M(f)T of M(f). 
Proof. Recall that the entry aij in row i and column j of M(f) is the i-th coordinate of f(uj) over the basis (v1,...,vm). By de.nition of v ., we have£¨v .,f(uj)£© = aij. The entry 
ii 
T
aji in row j and column i of M(fT) is the j-th coordinate of 
. T. T. T. 
fT(vi )= a1 iu1 + ¡¤¡¤¡¤ + ajiuj + ¡¤¡¤¡¤ + aniun 
over the basis (u1. ,...,u . ), which is just aT = fT(vi .)(uj)=£¨fT(vi .),uj£©. Since
n	ji 
£¨v . ,f(uj)£© =£¨fT(vi . ),uj£©,
i 
T
we have aij = aji, proving that M(fT)= M(f)T . 
We now can give a very short proof of the fact that the rank of a matrix is equal to the rank of its transpose. 
Proposition 10.13. Given an m ¡Á n matrix A over a .eld K, we have rk(A) = rk(AT). 
Proof. The matrix A corresponds to a linear map f : Kn ¡ú Km, and by Theorem 10.10, rk(f) = rk(fT). By Proposition 10.12, the linear map fT corresponds to AT . Since rk(A)= rk(f), and rk(AT) = rk(fT), we conclude that rk(A) = rk(AT). 
Thus, given an m¡Án-matrix A, the maximum number of linearly independent columns is equal to the maximum number of linearly independent rows. There are other ways of proving this fact that do not involve the dual space, but instead some elementary transformations on rows and columns. 
Proposition 10.13 immediately yields the following criterion for determining the rank of a matrix: 
Proposition 10.14. Given any m¡Án matrix A over a .eld K (typically K = R or K = C), the rank of A is the maximum natural number r such that there is an invertible r¡Ár submatrix of A obtained by selecting r rows and r columns of A. 
For example, the 3 ¡Á 2 matrix 
.. 
a11 a12 
..
A = 	a21 a22 a31 a32 
has rank 2 i. one of the three 2 ¡Á 2 matrices 
a11 a12 a11 a12 a21 a22 a21 a22 a31 a32 a31 a32 
is invertible. 
If we combine Proposition 6.12 with Proposition 10.14, we obtain the following criterion for .nding the rank of a matrix. 
Proposition 10.15. Given any m¡Án matrix A over a .eld K (typically K = R or K = C), the rank of A is the maximum natural number r such that there is an r ¡Á r submatrix B of A obtained by selecting r rows and r columns of A, such that det(B)=0. 
This is not a very e.cient way of .nding the rank of a matrix. We will see that there are better ways using various decompositions such as LU, QR, or SVD. 
10.7 Properties of the Double Transpose 
First we have the following property showing the naturality of the eval map. 
Proposition 10.16. For any linear map f : E ¡ú F , we have 
fTT . evalE = evalF . f, 
or equivalently the following diagram commutes: 
E.. f ..
F
    
evalE
evalF
E F. 
f 
Proof. For every u ¡Ê E and every . ¡Ê F ., we have 
(fTT . evalE)(u)(.)=£¨fTT(evalE(u)),.£© =£¨evalE(u),fT(.)£© =£¨fT(.),u£© =£¨., f(u)£© =£¨evalF (f(u)),.£© =£¨(evalF . f)(u),.£©= (evalF . f)(u)(.), 
which proves that fTT . evalE = evalF . f, as claimed. 
If E and F are .nite-dimensional, then evalE and evalF are isomorphisms, so Proposition 
10.16 shows that fTT = evalF . f . eval.1 . (.)
E 
The above equation is often interpreted as follows: if we identify E with its bidual E.. and F with its bidual F .., then fTT = f. This is an abuse of notation; the rigorous statement is (.). 
As a corollary of Proposition 10.16, we obtain the following result. 
10.7. PROPERTIES OF THE DOUBLE TRANSPOSE 
Proposition 10.17. If dim(E) is .nite, then we have 
Ker (fTT) = evalE(Ker (f)). 

Proof. Indeed, if E is .nite-dimensional, the map evalE : E ¡ú E.. is an isomorphism, so every . ¡Ê E.. is of the form . = evalE(u) for some u ¡Ê E, the map evalF : F ¡ú F .. is injective, and we have 
fTT(.)=0 	i. fTT(evalE(u)) = 0 i. evalF (f(u)) = 0 i. f(u)=0 i. u ¡Ê Ker (f) i. . ¡Ê evalE(Ker (f)), 
which proves that Ker (fTT) = evalE(Ker (f)). 
Remarks: If dim(E) is .nite, following an argument of Dan Guralnik, the fact that rk(f)= rk(fT) can be proven using Proposition 10.17. Proof. We know from Proposition 10.9 applied to fT : F . ¡ú E. that Ker (fTT) = (Im fT)0 , and we showed in Proposition 10.17 that Ker (fTT) = evalE(Ker (f)). It follows (since evalE is an isomorphism) that dim((Im fT)0) = dim(Ker (fTT)) = dim(Ker (f)) = dim(E) . dim(Im f), and since dim(Im fT) + dim((Im fT)0) = dim(E), we get dim(Im fT) = dim(Im f). As indicated by Dan Guralnik, if dim(E) is .nite, the above result can be used to prove the following result. Proposition 10.18. If dim(E) is .nite, then for any linear map f : E ¡ú F , we have 
Im fT = (Ker (f))0 . 
Proof. From
£¨fT(.),u£© =£¨., f(u)£© for all . ¡Ê F . and all u ¡Ê E, we see that if u ¡Ê Ker (f), then£¨fT(.),u£© =£¨., 0£© = 0, which means that fT(.) ¡Ê (Ker (f))0, and thus, Im fT . (Ker (f))0 . For the converse, since dim(E) is .nite, we have 
dim((Ker (f))0) = dim(E) . dim(Ker (f)) = dim(Im f), but we just proved that dim(Im fT) = dim(Im f), so we get dim((Ker (f))0) = dim(Im fT), and since Im fT . (Ker (f))0, we obtain Im fT = (Ker (f))0 , as claimed. 
Remarks: 
1. By the duality theorem, since (Ker (f))00 = Ker (f), the above equation yields another 
proof of the fact that 
Ker (f) = (Im fT)0 , 

when E is .nite-dimensional. 
2. The equation 
Im fT = (Ker (f))0 

is actually valid even if when E if in.nite-dimensional, but we will not prove this here. 
10.8 The Four Fundamental Subspaces 
Given a linear map f : E ¡ú F (where E and F are .nite-dimensional), Proposition 10.9 revealed that the four spaces 
Im f, Im fT , Ker f, Ker fT 
play a special role. They are often called the fundamental subspaces associated with f. These spaces are related in an intimate manner, since Proposition 10.9 shows that 
Ker f = (Im fT)0 Ker fT = (Im f)0 , 
10.8. THE FOUR FUNDAMENTAL SUBSPACES 
and Theorem 10.10 shows that 
rk(f) = rk(fT). 
It is instructive to translate these relations in terms of matrices (actually, certain linear algebra books make a big deal about this!). If dim(E)= n and dim(F )= m, given any basis (u1,...,un) of E and a basis (v1,...,vm) of F , we know that f is represented by an m ¡Á n matrix A =(aij), where the jth column of A is equal to f(uj) over the basis (v1,...,vm). Furthermore, the transpose map fT is represented by the n ¡Á m matrix AT (with respect to the dual bases). Consequently, the four fundamental spaces 
Im f, Im fT , Ker f, Ker fT 
correspond to 
(1) 
The column space of A, denoted by Im A or R(A); this is the subspace of Rm spanned by the columns of A, which corresponds to the image Im f of f. 

(2) 
The kernel or nullspace of A, denoted by Ker A or N (A); this is the subspace of Rn consisting of all vectors x ¡Ê Rn such that Ax = 0. 

(3) 
The row space of A, denoted by Im AT or R(AT); this is the subspace of Rn spanned by the rows of A, or equivalently, spanned by the columns of AT, which corresponds to the image Im fT of fT . 

(4) 
The left kernel or left nullspace of A denoted by Ker AT or N (AT); this is the kernel (nullspace) of AT , the subspace of Rm consisting of all vectors y ¡Ê Rm such that ATy = 0, or equivalently, yTA = 0. 


Recall that the dimension r of Im f, which is also equal to the dimension of the column space Im A = R(A), is the rank of A (and f). Then, some our previous results can be reformulated as follows: 
1. 
The column space R(A) of A has dimension r. 

2. 
The nullspace N (A) of A has dimension n . r. 

3. 
The row space R(AT) has dimension r. 

4. 
The left nullspace N (AT) of A has dimension m . r. 


The above statements constitute what Strang calls the Fundamental Theorem of Linear Algebra, Part I (see Strang [64]). 
The two statements 
Ker f = (Im fT)0 Ker fT = (Im f)0 
translate to (1) The nullspace of A is the orthogonal of the row space of A. 
(2) The left nullspace of A is the orthogonal of the column space of A. The above statements constitute what Strang calls the Fundamental Theorem of Linear Algebra, Part II (see Strang [64]). 
Since vectors are represented by column vectors and linear forms by row vectors (over a basis in E or F ), a vector x ¡Ê Rn is orthogonal to a linear form y i. 
yx =0. 
Then, a vector x ¡Ê Rn is orthogonal to the row space of A i. x is orthogonal to every row of A, namely Ax = 0, which is equivalent to the fact that x belong to the nullspace of A. Similarly, the column vector y ¡Ê Rm (representing a linear form over the dual basis of F .) belongs to the nullspace of AT i. ATy = 0, i. yTA = 0, which means that the linear form given by yT (over the basis in F ) is orthogonal to the column space of A. 
Since (2) is equivalent to the fact that the column space of A is equal to the orthogonal of the left nullspace of A, we get the following criterion for the solvability of an equation of the form Ax = b: 
The equation Ax = b has a solution i. for all y ¡Ê Rm, if ATy =0, then yTb = 0. Indeed, the condition on the right-hand side says that b is orthogonal to the left nullspace of A; that is, b belongs to the column space of A. This criterion can be cheaper to check that checking directly that b is spanned by the columns of A. For example, if we consider the system 
x1 . x2 = b1 x2 . x3 = b2 x3 . x1 = b3 
which, in matrix form, is written Ax = b as below: 
.  ..  .  .  .  
1  .1  0  x1  b1  
. 0  1  .1..x2. = .b2.,  
.1  0  1  x3  b3  

we see that the rows of the matrix A add up to 0. In fact, it is easy to convince ourselves that the left nullspace of A is spanned by y = (1, 1, 1), and so the system is solvable i. yTb = 0, namely 
b1 + b2 + b3 =0. Note that the above criterion can also be stated negatively as follows: The equation Ax = b has no solution i. there is some y ¡Ê Rm such that ATy =0 and yTb = 0. Since ATy = 0 i. yTA = 0, we can view yT as a row vector representing a linear form, and yTA = 0 asserts that the linear form yT vanishes on the columns A1,...,An of A but does not vanish on b. Since the linear form yT de.nes the hyperplane H of equation yTz =0 (with z ¡Ê Rm), geometrically the equation Ax = b has no solution i. there is a hyperplane H containing A1,...,An and not containing b. 
10.9. SUMMARY 
10.9 Summary 
The main concepts and results of this chapter are listed below: 
. 	The dual space E. and linear forms (covector). The bidual E.. . 
. 	The bilinear pairing£¨.,.£©: E. ¡Á E ¡ú K (the canonical pairing). 

. 	Evaluation at v: evalv : E. ¡ú K. 

. 	The map evalE : E ¡ú E.. . 

. 	Othogonality between a subspace V of E and a subspace U of E.; the orthogonal V 0 and the orthogonal U0 . 

. 	Coordinate forms. 

. 	The Duality theorem (Theorem 10.4). 

. 	The dual basis of a basis. 

. 	The isomorphism evalE : E ¡ú E.. when dim(E) is .nite. 

. 	Pairing between two vector spaces; nondegenerate pairing; Proposition 10.6. 

. 	Hyperplanes and linear forms. 
. ¡ú E.


. 	The transpose fT : F of a linear map f : E ¡ú F . 

. 	The fundamental identities: 


Ker fT = (Im f)0 and Ker f = (Im fT)0 
(Proposition 10.9). 
. 	If F is .nite-dimensional, then 
rk(f) = rk(fT). 
(Theorem 10.10). 
. 	The matrix of the transpose map fT is equal to the transpose of the matrix of the map f (Proposition 10.12). 

. 	For any m ¡Á n matrix A, 
rk(A) = rk(AT). 


. 	Characterization of the rank of a matrix in terms of a maximal invertible submatrix (Proposition 10.14). 

. 	
The four fundamental subspaces: 
Im f, Im fT , Ker f, Ker fT . 


. 	
The column space, the nullspace, the row space, and the left nullspace (of a matrix). 

. 	
Criterion for the solvability of an equation of the form Ax = b in terms of the left nullspace. 


10.10 Problems 
Problem 10.1. Prove the following properties of transposition: (f + g)T = fT + g T (g . f)T = fT . g T idT 
E = idE. . 
Problem 10.2. Let (u1,...,un.1) be n . 1 linearly independent vectors ui ¡Ê Cn . Prove that the hyperlane H spanned by (u1,...,un.1) is the nullspace of the linear form 
x ¡ú det(u1,...,un.1,x),x ¡Ê Cn . 
Prove that if A is the n ¡Á n matrix whose columns are (u1,...,un.1,x), and if ci = (.1)i+n det(Ain) is the cofactor of ain = xi for i =1,...,n, then H is de.ned by the equation 
c1x1 + ¡¤¡¤¡¤ + cnxn =0. 
Problem 10.3. (1) Let .: Rn ¡Á Rn ¡ú R be the map de.ned by 
.((x1,...,xn), (y1,...,yn)) = x1y1 + ¡¤¡¤¡¤ + xnyn. 
Prove that . is a bilinear nondegenerate pairing. Deduce that (Rn). is isomorphic to Rn . Prove that .(x, x)=0 i. x = 0. 
(2) Let .L : R4 ¡Á R4 ¡ú R be the map de.ned by .L((x1,x2,x3,x4), (y1,y2,y3,,y4)) = x1y1 . x2y2 . x3y3 . x4y4. Prove that . is a bilinear nondegenerate pairing. Show that there exist nonzero vectors x ¡Ê R4 such that .L(x, x) = 0. 
Remark: The vector space R4 equipped with the above bilinear form called the Lorentz form is called Minkowski space. 
10.10. PROBLEMS 
Problem 10.4. Given any two subspaces V1,V2 of a .nite-dimensional vector space E, prove that 
= V 0 ¡É V 0
(V1 + V2)012 (V1 ¡É V2)0 = V 0 + V 0 
12 . 
Beware that in the second equation, V1 and V2 are subspaces of E, not E. . 
Hint. To prove the second equation, prove the inclusions V10 + V20 . (V1 ¡É V2)0 and (V1 ¡É V2)0 . V10 + V20 . Proving the second inclusion is a little tricky. First, prove that we can pick a subspace W1 of V1 and a subspace W2 of V2 such that 
1. 
V1 is the direct sum V1 =(V1 ¡É V2) ¨’ W1. 

2. 
V2 is the direct sum V2 =(V1 ¡É V2) ¨’ W2. 

3. 
V1 + V2 is the direct sum V1 + V2 =(V1 ¡É V2) ¨’ W1 ¨’ W2. 


Problem 10.5. (1) Let A be any n ¡Á n matrix such that the sum of the entries of every row of A is the same (say c1), and the sum of entries of every column of A is the same (say c2). Prove that c1 = c2. 
(2) Prove that for any n ¡Ý 2, the 2n . 2 equations asserting that the sum of the entries of every row of A is the same, and the sum of entries of every column of A is the same are lineary independent. For example, when n = 4, we have the following 6 equations 
a11 + a12 + a13 + a14 . a21 . a22 . a23 . a24 =0 a21 + a22 + a23 + a24 . a31 . a32 . a33 . a34 =0 a31 + a32 + a33 + a34 . a41 . a42 . a43 . a44 =0 a11 + a21 + a31 + a41 . a12 . a22 . a32 . a42 =0 a12 + a22 + a32 + a42 . a13 . a23 . a33 . a43 =0 a13 + a23 + a33 + a43 . a14 . a24 . a34 . a44 =0. 
Hint. Group the equations as above; that is, .rst list the n . 1 equations relating the rows, and then list the n . 1 equations relating the columns. Prove that the .rst n . 1 equations are linearly independent, and that the last n . 1 equations are also linearly independent. Then, .nd a relationship between the two groups of equations that will allow you to prove that they span subspace V r and V c such that V r ¡É V c = (0). 
(3) Now consider magic squares. Such matrices satisfy the two conditions about the sum of the entries in each row and in each column to be the same number, and also the additional two constraints that the main descending and the main ascending diagonals add up to this common number. Traditionally, it is also required that the entries in a magic square are positive integers, but we will consider generalized magic square with arbitrary real entries. For example, in the case n = 4, we have the following system of 8 equations: 
a11 + a12 + a13 + a14 . a21 . a22 . a23 . a24 =0 a21 + a22 + a23 + a24 . a31 . a32 . a33 . a34 =0 a31 + a32 + a33 + a34 . a41 . a42 . a43 . a44 =0 a11 + a21 + a31 + a41 . a12 . a22 . a32 . a42 =0 a12 + a22 + a32 + a42 . a13 . a23 . a33 . a43 =0 a13 + a23 + a33 + a43 . a14 . a24 . a34 . a44 =0 
a22 + a33 + a44 . a12 . a13 . a14 =0 a41 + a32 + a23 . a11 . a12 . a13 =0. 
In general, the equation involving the descending diagonal is 
a22 + a33 + ¡¤¡¤¡¤ + ann . a12 . a13 .¡¤ ¡¤¡¤. a1n =0 (r) 
and the equation involving the ascending diagonal is 
an1 + an.12 + ¡¤¡¤¡¤ + a2n.1 . a11 . a12 .¡¤ ¡¤¡¤. a1n.1 =0. (c) 
Prove that if n ¡Ý 3, then the 2n equations asserting that a matrix is a generalized magic square are linearly independent. Hint. Equations are really linear forms, so .nd some matrix annihilated by all equations except equation r, and some matrix annihilated by all equations except equation c. 
Problem 10.6. Let U1,...,Up be some subspaces of a vector space E, and assume that they form a direct sum U = U1 ¨’¡¤ ¡¤¡¤¨’ Up. Let ji : Ui ¡ú U1 ¨’¡¤ ¡¤¡¤¨’ Up be the canonical injections, and let ¦Ði : U1 . ¡Á¡¤ ¡¤¡¤¡Á Up . ¡ú Ui . be the canonical projections. Prove that there is an isomorphism f from (U1 ¨’¡¤ ¡¤¡¤¨’ Up). to U1 . ¡Á¡¤ ¡¤¡¤¡Á Up . such that 
¦Ði . f = ji T , 1 ¡Ü i ¡Ü p. 
Problem 10.7. Let U and V be two subspaces of a vector space E such that E = U ¨’ V . Prove that 
E . = U0 ¨’ V 0 
. 


