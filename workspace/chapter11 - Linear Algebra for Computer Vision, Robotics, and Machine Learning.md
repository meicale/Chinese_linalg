Chapter 11 
Euclidean Spaces 
Rien n¡¯est beau que le vrai. 
¡ªHermann Minkowski 
11.1 Inner Products, Euclidean Spaces 
So far the framework of vector spaces allows us to deal with ratios of vectors and linear combinations, but there is no way to express the notion of angle or to talk about orthogonality of vectors. A Euclidean structure allows us to deal with metric notions such as angles, orthogonality, and length (or distance). 
This chapter covers the bare bones of Euclidean geometry. Deeper aspects of Euclidean geometry are investigated in Chapter 12. One of our main goals is to give the basic properties of the transformations that preserve the Euclidean structure, rotations and re.ections, since they play an important role in practice. Euclidean geometry is the study of properties invariant under certain a.ne maps called rigid motions. Rigid motions are the maps that preserve the distance between points. 
We begin by de.ning inner products and Euclidean spaces. The Cauchy¨CSchwarz in-equality and the Minkowski inequality are shown. We de.ne orthogonality of vectors and of subspaces, orthogonal bases, and orthonormal bases. We prove that every .nite-dimensional Euclidean space has orthonormal bases. The .rst proof uses duality and the second one the Gram¨CSchmidt orthogonalization procedure. The QR-decomposition for invertible matrices is shown as an application of the Gram¨CSchmidt procedure. Linear isometries (also called or-thogonal transformations) are de.ned and studied brie.y. We conclude with a short section in which some applications of Euclidean geometry are sketched. One of the most important applications, the method of least squares, is discussed in Chapter 21. 
For a more detailed treatment of Euclidean geometry see Berger [5, 6], Snapper and Troyer [59], or any other book on geometry, such as Pedoe [51], Coxeter [15], Fresnel [22], Tisseron [67], or Cagnac, Ramis, and Commeau [12]. Serious readers should consult Emil Artin¡¯s famous book [2], which contains an in-depth study of the orthogonal group, as well 
377 

as other groups arising in geometry. It is still worth consulting some of the older classics, such as Hadamard [31, 32] and Rouch¡äe and de Comberousse [52]. The .rst edition of [31] was published in 1898 and .nally reached its thirteenth edition in 1947! In this chapter it is assumed that all vector spaces are de.ned over the .eld R of real numbers unless speci.ed otherwise (in a few cases, over the complex numbers C). 
First we de.ne a Euclidean structure on a vector space. Technically, a Euclidean structure over a vector space E is provided by a symmetric bilinear form on the vector space satisfying some extra properties. Recall that a bilinear form .: E ¡Á E ¡ú R is de.nite if for every u ¡Ê E, u  
= 0 implies that .(u, u)= 0, and positive if for every u ¡Ê E, .(u, u) ¡Ý 0. 
De.nition 11.1. A Euclidean space is a real vector space E equipped with a symmetric bilinear form .: E ¡ÁE ¡ú R that is positive de.nite. More explicitly, . : E ¡ÁE ¡ú R satis.es the following axioms: 
.(u1 + u2,v)= .(u1,v)+ .(u2,v), 
.(u, v1 + v2)= .(u, v1)+ .(u, v2), 
.(¦Ëu, v)= ¦Ë.(u, v), 
.(u, ¦Ëv)= ¦Ë.(u, v), 
.(u, v)= .(v, u), u 0 implies that .(u, u) > 0.
= 
The real number .(u, v) is also called the inner product (or scalar product) of u and v. We also de.ne the quadratic form associated with . as the function ¦µ: E ¡ú R+ such that 
¦µ(u)= .(u, u), 
for all u ¡Ê E. 
Since . is bilinear, we have .(0, 0) = 0, and since it is positive de.nite, we have the stronger fact that .(u, u)=0 i. u =0, 
that is, ¦µ(u)=0 i. u = 0. Given an inner product .: E ¡Á E ¡ú R on a vector space E, we also denote .(u, v) by 
u ¡¤ v or Ju, v¡¢ or (u|v), 
 
and¦µ(u) by lul. 
Example 11.1. The standard example of a Euclidean space is Rn, under the inner product ¡¤ de.ned such that 
(x1,...,xn) ¡¤ (y1,...,yn)= x1y1 + x2y2 + ¡¤¡¤¡¤ + xnyn. 
This Euclidean space is denoted by En . 
11.1. INNER PRODUCTS, EUCLIDEAN SPACES 
There are other examples. Example 11.2. For instance, let E be a vector space of dimension 2, and let (e1,e2) be a basis of E. If a> 0 and b2 . ac < 0, the bilinear form de.ned such that .(x1e1 + y1e2,x2e1 + y2e2)= ax1x2 + b(x1y2 + x2y1)+ cy1y2 yields a Euclidean structure on E. In this case, ¦µ(xe1 + ye2)= ax 2 +2bxy + cy 2 . Example 11.3. Let C[a, b] denote the set of continuous functions f :[a, b] ¡ú R. It is easily checked that C[a, b] is a vector space of in.nite dimension. Given any two functions f, g ¡ÊC[a, b], let 
 b 
Jf, g¡¢ =f(t)g(t)dt. 
a 
We leave it as an easy exercise that J.,.¡¢ is indeed an inner product on C[a, b]. In the case where a = .¦Ð and b = ¦Ð (or a = 0 and b =2¦Ð, this makes basically no di.erence), one should compute 
Jsin px, sin qx¡¢, Jsin px, cos qx¡¢, and Jcos px, cos qx¡¢, 
for all natural numbers p, q ¡Ý 1. The outcome of these calculations is what makes Fourier analysis possible! 
Example 11.4. Let E =Mn(R) be the vector space of real n ¡Á n matrices. If we view a matrix A ¡Ê Mn(R) as a ¡°long¡± column vector obtained by concatenating together its columns, we can de.ne the inner product of two matrices A, B ¡Ê Mn(R) as 
n
n 
JA, B¡¢ = aijbij, i,j=1 
which can be conveniently written as 
JA, B¡¢ = tr(ATB) = tr(BTA). 
Since this can be viewed as the Euclidean product on Rn2 , it is an inner product on Mn(R). The corresponding norm 
lAlF = TA)
tr(Ais the Frobenius norm (see Section 8.2). Let us observe that . can be recovered from ¦µ. 
Proposition 11.1. We have 
1 
.(u, v) = [¦µ(u + v) . ¦µ(u) . ¦µ(v)]
2
for all u, v ¡Ê E. We say that . is the polar form of ¦µ. 
Proof. By bilinearity and symmetry, we have 
¦µ(u + v)= .(u + v, u + v) = .(u, u + v)+ .(v, u + v) = .(u, u)+2.(u, v)+ .(v, v) = ¦µ(u)+2.(u, v) + ¦µ(v). 
Note that we are committing an abuse of notation since x = in =1 xiei is a vector in E, but the column vector associated with (x1,...,xn) belongs to Rn . To avoid this minor abuse, we could denote the column vector associated with (x1,...,xn) by x (and similarly y for the column vector associated with (y1,...,yn)), in wich case the ¡°correct¡± expression for .(x, y) is 
.(x, y)= x TGy. 
However, in view of the isomorphism between E and Rn , to keep notation as simple as possible, we will use x and y instead of x and y. 
Also observe that . is symmetric i. G = GT, and . is positive de.nite i. the matrix G is positive de.nite, that is, 
x TGx > 0 for all x ¡Ê Rn,x =0. 
The matrix G associated with an inner product is called the Gram matrix of the inner product with respect to the basis (e1,...,en). Conversely, if A is a symmetric positive de.nite n ¡Á n matrix, it is easy to check that the bilinear form Jx, y¡¢ = x TAy 
a
If E is .nite-dimensional and if .: E ¡Á E ¡ú R is a bilinear form on E, given any basis (e1,...,en) of E, we can write x = n xiei and y = n yjej, and we have 
i=1 j=1 
a 
 

 

n
nni=1 j=1 i,j=1 
n .()= ..()x,yxeye= xye,e,.iijjijijIfwelet G bethematrix G =(.()),andif andarethecolumnvectorsassociated e,exyij
with (x1,...,xn) and (y1,...,yn), then we can write TGT
.(x, y)= x TGy = y x. 
a 

nn 

11.1. INNER PRODUCTS, EUCLIDEAN SPACES 
is an inner product. If we make a change of basis from the basis (e1,...,en) to the basis (f1,...,fn), and if the change of basis matrix is P (where the jth column of P consists of 
'
the coordinates of fj over the basis (e1,...,en)), then with respect to coordinates x' and yover the basis (f1,...,fn), we have 
x TGy = x'TP TGP y', 
so the matrix of our inner product over the basis (f1,...,fn) is P TGP . We summarize these facts in the following proposition. 
Proposition 11.2. Let E be a .nite-dimensional vector space, and let (e1,...,en) be a basis of E. 
1. For any inner product J.,.¡¢ on E, if G =(Jei,ej¡¢) is the Gram matrix of the inner product J.,.¡¢ w.r.t. the basis (e1,...,en), then G is symmetric positive de.nite. 

2. For any change of basis matrix P , the Gram matrix of J.,.¡¢ with respect to the new basis is P TGP . 


3. If A is any n ¡Á n symmetric positive de.nite matrix, then 
Jx, y¡¢ = x TAy 
is an inner product on E. 
We will see later that a symmetric matrix is positive de.nite i. its eigenvalues are all positive. 
One of the very important properties of an inner product . is that the map u ¡ú ¦µ(u) is a norm. 
Proposition 11.3. Let E be a Euclidean space with inner product ., and let ¦µ be the corresponding quadratic form. For all u, v ¡Ê E, we have the Cauchy¨CSchwarz inequality 
.(u, v)2 ¡Ü ¦µ(u)¦µ(v), 
the equality holding i. u and v are linearly dependent. We also have the Minkowski inequality 
¦µ(u + v) ¡Ü ¦µ(u)+ ¦µ(v), 
the equality holding i. u and v are linearly dependent, where in addition if u =0 and v =0, then u = ¦Ëv for some ¦Ë> 0. 
Proof. For any vectors u, v ¡Ê E, we de.ne the function T : R ¡ú R such that 
T (¦Ë) = ¦µ(u + ¦Ëv), 
for all ¦Ë ¡Ê R. Using bilinearity and symmetry, we have ¦µ(u + ¦Ëv)= .(u + ¦Ëv, u + ¦Ëv) = .(u, u + ¦Ëv)+ ¦Ë.(v, u + ¦Ëv) = .(u, u)+2¦Ë.(u, v)+ ¦Ë2.(v, v) = ¦µ(u)+2¦Ë.(u, v)+ ¦Ë2¦µ(v). Since . is positive de.nite, ¦µ is nonnegative, and thus T (¦Ë) ¡Ý 0 for all ¦Ë ¡Ê R. If ¦µ(v) = 0, then v = 0, and we also have .(u, v) = 0. In this case, the Cauchy¨CSchwarz inequality is trivial, and v = 0 and u are linearly dependent. Now assume ¦µ(v) > 0. Since T (¦Ë) ¡Ý 0, the quadratic equation ¦Ë2¦µ(v)+2¦Ë.(u, v) + ¦µ(u)=0 cannot have distinct real roots, which means that its discriminant ¦¤ = 4(.(u, v)2 . ¦µ(u)¦µ(v)) is null or negative, which is precisely the Cauchy¨CSchwarz inequality 
.(u, v)2 ¡Ü ¦µ(u)¦µ(v). Let us now consider the case where we have the equality 
.(u, v)2 = ¦µ(u)¦µ(v). There are two cases. If ¦µ(v) = 0, then v = 0 and u and v are linearly dependent. If ¦µ(v) = 0, then the above quadratic equation has a double root ¦Ë0, and we have ¦µ(u + ¦Ë0v) = 0. Since . is positive de.nite, ¦µ(u + ¦Ë0v) = 0 implies that u + ¦Ë0v = 0, which shows that u and v are linearly dependent. Conversely, it is easy to check that we have equality when u and v 
are linearly dependent. The Minkowski inequality  
¦µ(u + v) ¡Ü  ¦µ(u) +  ¦µ(v)  
is equivalent to  ¦µ(u + v) ¡Ü ¦µ(u) + ¦µ(v) + 2  ¦µ(u)¦µ(v).  
However, we have shown that  

2.(u, v) = ¦µ(u + v) . ¦µ(u) . ¦µ(v), and so the above inequality is equivalent to .(u, v) ¡Ü ¦µ(u)¦µ(v), 

11.1. INNER PRODUCTS, EUCLIDEAN SPACES 
which is trivial when .(u, v) ¡Ü 0, and follows from the Cauchy¨CSchwarz inequality when .(u, v) ¡Ý 0. Thus, the Minkowski inequality holds. Finally assume that u = 0 and v = 0, and that 
¦µ(u + v)= ¦µ(u)+ ¦µ(v). 
When this is the case, we have 
.(u, v)= ¦µ(u)¦µ(v), 
and we know from the discussion of the Cauchy¨CSchwarz inequality that the equality holds i. u and v are linearly dependent. The Minkowski inequality is an equality when u or v is null. Otherwise, if u = 0 and v = 0, then u = ¦Ëv for some ¦Ë = 0, and since 
.(u, v)= ¦Ë.(v, v)= ¦µ(u)¦µ(v), 
by positivity, we must have ¦Ë> 0. 
Note that the Cauchy¨CSchwarz inequality can also be written as 
|.(u, v)|¡Ü ¦µ(u) ¦µ(v). 
Remark: It is easy to prove that the Cauchy¨CSchwarz and the Minkowski inequalities still hold for a symmetric bilinear form that is positive, but not necessarily de.nite (i.e., .(u, v) ¡Ý 0 for all u, v ¡Ê E). However, u and v need not be linearly dependent when the equality holds. 
The Minkowski inequality 
¦µ(u + v) ¡Ü ¦µ(u)+ ¦µ(v) 
shows that the map u ¡ú ¦µ(u) satis.es the convexity inequality (also known as triangle inequality), condition (N3) of De.nition 8.1, and since . is bilinear and positive de.nite, it also satis.es conditions (N1) and (N2) of De.nition 8.1, and thus it is a norm on E. The norm induced by . is called the Euclidean norm induced by .. 
The Cauchy¨CSchwarz inequality can be written as 
|u ¡¤ v|¡Ülullvl, 
and the Minkowski inequality as 
lu + vl¡Ülul + lvl. 
If u and v are nonzero vectors then the Cauchy¨CSchwarz inequality implies that 
u ¡¤ v 

.1 ¡Ü¡Ü +1. 
lullvl 
Then there is a unique ¦È ¡Ê [0,¦Ð] such that u ¡¤ v 
cos ¦È = . 
lullvl 
We have u = v i. ¦È = 0 and u = .v i. ¦È = ¦Ð. For 0 <¦È<¦Ð, the vectors u and v are linearly independent and there is an orientation of the plane spanned by u and v such that ¦È is the angle between u and v. See Problem 11.8 for the precise notion of orientation. If u is a unit vector (which means that lul = 1), then the vector 
(lvl cos ¦È)u =(u ¡¤ v)u =(v ¡¤ u)u 
is called the orthogonal projection of v onto the space spanned by u. 
Remark: One might wonder if every norm on a vector space is induced by some Euclidean inner product. In general this is false, but remarkably, there is a simple necessary and su.cient condition, which is that the norm must satisfy the parallelogram law: 
lu + vl2 + lu . vl2 = 2(lul2 + lvl2). 
See Figure 11.1. 

Figure 11.1: The parallelogram law states that the sum of the lengths of the diagonals of the parallelogram determined by vectors u and v equals the sum of all the sides. 
If J.,.¡¢ is an inner product, then we have 
lu + vl2 = lul2 + lvl2 +2Ju, v¡¢ 
lu . vl2 = lul2 + lvl2 . 2Ju, v¡¢, and by adding and subtracting these identities, we get the parallelogram law and the equation 
Ju, v¡¢ = 1(lu + vl2 .lu . vl2),
4

11.1. INNER PRODUCTS, EUCLIDEAN SPACES 
which allows us to recover J.,.¡¢ from the norm. 
Conversely, if ll is a norm satisfying the parallelogram law, and if it comes from an inner product, then this inner product must be given by 
Ju, v¡¢ = 14(lu + vl2 .lu . vl2). 
We need to prove that the above form is indeed symmetric and bilinear. 
Symmetry holds because lu . vl = l.(u . v)l = lv . ul. Let us prove additivity in the variable u. By the parallelogram law, we have 
2(lx + zl2 + lyl2)= lx + y + zl2 + lx . y + zl2 
which yields 
lx + y + zl2 = 2(lx + zl2 + lyl2) .lx . y + zl2 lx + y + zl2 = 2(ly + zl2 + lxl2) .ly . x + zl2 , 
where the second formula is obtained by swapping x and y. Then by adding up these equations, we get 
lx + y + zl2 = lxl2 + lyl2 + lx + zl2 + ly + zl2 . 1 lx . y + zl2 . 1 ly . x + zl2 . 
22 Replacing z by .z in the above equation, we get 
lx + y . zl2 = lxl2 + lyl2 + lx . zl2 + ly . zl2 
. 1 lx . y . zl2 . 1 ly . x . zl2 ,
22 
Since lx . y + zl = l.(x . y + z)l = ly . x . zl and ly . x + zl = l.(y . x + z)l = lx . y . zl, by subtracting the last two equations, we get 
Jx + y, z¡¢ = 1(lx + y + zl2 .lx + y . zl2)
4= 1(lx + zl2 .lx . zl2)+ 1(ly + zl2 .ly . zl2)
44= Jx, z¡¢ + Jy, z¡¢, 
as desired. Proving that J¦Ëx, y¡¢ = ¦ËJx, y¡¢ for all ¦Ë ¡Ê R 
is a little tricky. The strategy is to prove the identity for ¦Ë ¡Ê Z, then to promote it to Q, and then to R by continuity. 
Since 
J.u, v¡¢ = 14(l.u + vl2 . l.u . vl2) 
= 1(lu . vl2 .lu + vl2)
4= .Ju, v¡¢, 
the property holds for ¦Ë = .1. By linearity and by induction, for any n ¡Ê N with n ¡Ý 1, writing n = n . 1+1, we get 
J¦Ëx, y¡¢ = ¦ËJx, y¡¢ for all ¦Ë ¡Ê N, 
and since the above also holds for ¦Ë = .1, it holds for all ¦Ë ¡Ê Z. For ¦Ë = p/q with p, q ¡Ê Z and q = 0, we have 
qJ(p/q)u, v¡¢ = Jpu, v¡¢ = pJu, v¡¢, 
which shows that J(p/q)u, v¡¢ =(p/q)Ju, v¡¢, 
and thus J¦Ëx, y¡¢ = ¦ËJx, y¡¢ for all ¦Ë ¡Ê Q. 
To .nish the proof, we use the fact that a norm is a continuous map x ¡úlxl. Then, the continuous function t ¡ú 1 Jtu, v¡¢ de.ned on R.{0} agrees with Ju, v¡¢ on Q.{0}, so it is 
t 
equal to Ju, v¡¢ on R.{0}. The case ¦Ë = 0 is trivial, so we are done. We now de.ne orthogonality. 


11.2 Orthogonality and Duality in Euclidean Spaces 
An inner product on a vector space gives the ability to de.ne the notion of orthogonality. Families of nonnull pairwise orthogonal vectors must be linearly independent. They are called orthogonal families. In a vector space of .nite dimension it is always possible to .nd orthogonal bases. This is very useful theoretically and practically. Indeed, in an orthogonal basis, .nding the coordinates of a vector is very cheap: It takes an inner product. Fourier series make crucial use of this fact. When E has .nite dimension, we prove that the inner product on E induces a natural isomorphism between E and its dual space E. . This allows us to de.ne the adjoint of a linear map in an intrinsic fashion (i.e., independently of bases). It is also possible to orthonormalize any basis (certainly when the dimension is .nite). We give two proofs, one using duality, the other more constructive using the Gram¨CSchmidt orthonormalization procedure. 
De.nition 11.2. Given a Euclidean space E, any two vectors u, v ¡Ê E are orthogonal, or perpendicular, if u ¡¤ v = 0. Given a family (ui)i¡ÊI of vectors in E, we say that (ui)i¡ÊI is orthogonal if ui ¡¤ uj = 0 for all i, j ¡Ê I, where i = j. We say that the family (ui)i¡ÊI is 
11.2. ORTHOGONALITY AND DUALITY IN EUCLIDEAN SPACES 
orthonormal if ui ¡¤ uj = 0 for all i, j ¡Ê I, where i = j, and luil = ui ¡¤ ui = 1, for all i ¡Ê I. For any subset F of E, the set F ¡Í = {v ¡Ê E | u ¡¤ v =0, for all u ¡Ê F }, of all vectors orthogonal to all vectors in F , is called the orthogonal complement of F . Since inner products are positive de.nite, observe that for any vector u ¡Ê E, we have u ¡¤ v = 0 for all v ¡Ê E i. u =0. It is immediately veri.ed that the orthogonal complement F ¡Í of F is a subspace of E. Example 11.5. Going back to Example 11.3 and to the inner product 
¦Ð 
Jf, g¡¢ = f(t)g(t)dt 
.¦Ð 
on the vector space C[.¦Ð, ¦Ð], it is easily checked that 
 
¦Ð if p = q, p, q ¡Ý 1,
Jsin px, sin qx¡¢ =0if p = q, p, q ¡Ý 1, 
 
¦Ð if p = q, p, q ¡Ý 1,
Jcos px, cos qx¡¢ =0if p = q, p, q ¡Ý 0, 
and Jsin px, cos qx¡¢ =0, 
for all p ¡Ý 1 and q ¡Ý 0, and of course, J1, 1¡¢ = ¦Ð dx =2¦Ð.
.¦Ð 
As a consequence, the family (sin px)p¡Ý1 ¡È(cos qx)q¡Ý0 is orthogonal. It is not orthonormal, 
¡Ì¡Ì 
but becomes so if we divide every trigonometric function by ¦Ð,and1by 2¦Ð. 
Proposition 11.4. Given a Euclidean space E, for any family (ui)i¡ÊI of nonnull vectors in E, if (ui)i¡ÊI is orthogonal, then it is linearly independent. 
Proof. Assume there is a linear dependence 
n 
¦Ëjuj =0 j¡ÊJ 
for some ¦Ëj ¡Ê R and some .nite subset J of I. By taking the inner product with ui for any i ¡Ê J, and using the the bilinearity of the inner product and the fact that ui ¡¤ uj =0 whenever i = j, we get 
  
n 
0= ui ¡¤ 0= ui ¡¤¦Ëjujj¡ÊJ 
n 
= ¦Ëj(ui ¡¤ uj)= ¦Ëi(ui ¡¤ ui), j¡ÊJ 
so 
¦Ëi(ui ¡¤ ui)=0, for all i ¡Ê J, and since ui = 0 and an inner product is positive de.nite, ui ¡¤ ui = 0, so we obtain ¦Ëi =0, for all i ¡Ê J, which shows that the family (ui)i¡ÊI is linearly independent. We leave the following simple result as an exercise. 
Proposition 11.5. Given a Euclidean space E, any two vectors u, v ¡Ê E are orthogonal i. lu + vl2 = lul2 + lvl2 . 
See Figure 11.2 for a geometrical interpretation. 

Figure 11.2: The sum of the lengths of the two sides of a right triangle is equal to the length of the hypotenuse; i.e. the Pythagorean theorem. 
One of the most useful features of orthonormal bases is that they a.ord a very simple method for computing the coordinates of a vector over any basis vector. Indeed, assume that (e1,...,em) is an orthonormal basis. For any vector 
x = x1e1 + ¡¤¡¤¡¤ + xmem, if we compute the inner product x ¡¤ ei, we get x ¡¤ ei = x1e1 ¡¤ ei + ¡¤¡¤¡¤ + xiei ¡¤ ei + ¡¤¡¤¡¤ + xmem ¡¤ ei = xi, 
since 1 if i = j, 
ei ¡¤ ej = 
0 if i = j 
11.2. ORTHOGONALITY AND DUALITY IN EUCLIDEAN SPACES 

Figure 11.3: The orthogonal projection of the red vector x onto the black basis vector ei is the maroon vector xiei. Observe that x ¡¤ ei = lxl cos ¦È. 
is the property characterizing an orthonormal family. Thus, 
xi = x ¡¤ ei, 
which means that xiei =(x ¡¤ ei)ei is the orthogonal projection of x onto the subspace generated by the basis vector ei. See Figure 11.3. If the basis is orthogonal but not necessarily orthonormal, then 
x ¡¤ ei x ¡¤ ei 
xi == . 
ei ¡¤ ei leil2 
All this is true even for an in.nite orthonormal (or orthogonal) basis (ei)i¡ÊI .
A However, remember that every vector x is expressed as a linear combination 
n 
x = xiei 
i¡ÊI 
where the family of scalars (xi)i¡ÊI has .nite support, which means that xi = 0 for all i ¡Ê I . J, where J is a .nite set. Thus, even though the family (sin px)p¡Ý1 ¡È (cos qx)q¡Ý0 is orthogonal (it is not orthonormal, but becomes so if we divide every trigonometric function by
¡Ì¡Ì 
¦Ð,and1by 2¦Ð; we won¡¯t because it looks messy!), the fact that a function f ¡ÊC0[.¦Ð, ¦Ð] can be written as a Fourier series as 
¡Þ
n 
f(x)= a0 +(ak cos kx + bk sin kx) k=1 
does not mean that (sin px)p¡Ý1 ¡È (cos qx)q¡Ý0 is a basis of this vector space of functions, because in general, the families (ak) and (bk) do not have .nite support! In order for this in.nite linear combination to make sense, it is necessary to prove that the partial sums 
n
n 
a0 +(ak cos kx + bk sin kx) k=1 
of the series converge to a limit when n goes to in.nity. This requires a topology on the space. 
A very important property of Euclidean spaces of .nite dimension is that the inner product induces a canonical bijection (i.e., independent of the choice of bases) between the vector space E and its dual E. . The reason is that an inner product ¡¤: E ¡Á E ¡ú R de.nes a nondegenerate pairing, as de.ned in De.nition 10.4. Indeed, if u ¡¤ v = 0 for all v ¡Ê E then u = 0, and similarly if u ¡¤ v = 0 for all u ¡Ê E then v = 0 (since an inner product is positive de.nite and symmetric). By Proposition 10.6, there is a canonical isomorphism between E and E. . We feel that the reader will appreciate if we exhibit this mapping explicitly and reprove that it is an isomorphism. 
The mapping from E to E. is de.ned as follows. 
De.nition 11.3. For any vector u ¡Ê E, let .u : E ¡ú R be the map de.ned such that 
.u(v)= u ¡¤ v, for all v ¡Ê E. 
Since the inner product is bilinear, the map .u is a linear form in E. . Thus, we have a map 
b: E ¡ú E., de.ned such that b(u)= .u. 
Theorem 11.6. Given a Euclidean space E, the map b: E ¡ú E. de.ned such that 
b(u)= .u 
is linear and injective. When E is also of .nite dimension, the map b: E ¡ú E. is a canonical isomorphism. 
Proof. That b: E ¡ú E. is a linear map follows immediately from the fact that the inner product is bilinear. If .u = .v, then .u(w)= .v(w) for all w ¡Ê E, which by de.nition of .u means that u ¡¤ w = v ¡¤ w for all w ¡Ê E, which by bilinearity is equivalent to 
(v . u) ¡¤ w =0 
for all w ¡Ê E, which implies that u = v, since the inner product is positive de.nite. Thus, 
b: E ¡ú E. is injective. Finally, when E is of .nite dimension n, we know that E. is also of dimension n, and then b: E ¡ú E. is bijective. 
The inverse of the isomorphism b: E ¡ú E. is denoted by . : E. ¡ú E. 
As a consequence of Theorem 11.6 we have the following corollary. 
Corollary 11.7. If E is a Euclidean space of .nite dimension, every linear form f ¡Ê E. corresponds to a unique u ¡Ê E such that 
f(v)= u ¡¤ v, for every v ¡Ê E. 
In particular, if f is not the zero form, the kernel of f, which is a hyperplane H, is precisely the set of vectors that are orthogonal to u. 

11.2. ORTHOGONALITY AND DUALITY IN EUCLIDEAN SPACES 
Remarks: 
(1) 
The ¡°musical map¡± b: E ¡ú E. is not surjective when E has in.nite dimension. The result can be salvaged by restricting our attention to continuous linear maps, and by assuming that the vector space E is a Hilbert space (i.e., E is a complete normed vector space w.r.t. the Euclidean norm). This is the famous ¡°little¡± Riesz theorem (or Riesz representation theorem). 

(2) 
Theorem 11.6 still holds if the inner product on 	E is replaced by a nondegenerate symmetric bilinear form .. We say that a symmetric bilinear form .: E ¡Á E ¡ú R is nondegenerate if for every u ¡Ê E, 


if 	.(u, v)=0 forall v ¡Ê E, then u =0. 
For example, the symmetric bilinear form on R4 (the Lorentz form) de.ned such that 
.((x1,x2,x3,x4), (y1,y2,y3,y4)) = x1y1 + x2y2 + x3y3 . x4y4 
is nondegenerate. However, there are nonnull vectors u ¡Ê R4 such that .(u, u) = 0, which is impossible in a Euclidean space. Such vectors are called isotropic. 
Example 11.6. Consider Rn with its usual Euclidean inner product. Given any di.eren-tiable function f : U ¡ú R, where U is some open subset of Rn, by de.nition, for any x ¡Ê U, the total derivative dfx of f at x is the linear form de.ned so that for all u =(u1,...,un) ¡Ê Rn , 
.
. 

u1
n
ni=1
un 
.f .f 

.f 

.. 

..

.

dfx(u)= .x1 
(x) ¡¤¡¤¡¤ (x)

.xn (x) ui.
. 

= 

. 

.xi 
The unique vector v ¡Ê Rn such that v ¡¤ u = dfx(u) for all u ¡Ê Rn is the transpose of the Jacobian matrix of f at x, the 1 ¡Á n matrix .f .f 
(x) ¡¤¡¤¡¤ (x) . 
.x1 .xn This is the gradient grad(f)x of f at x, given by 
.
. 
grad(f)x = 
..... 

.f 
(x)
.x1 
. 
. 

. 
.f 
(x)
.xn 
..... 

. 

Example 11.7. Given any two vectors u, v ¡Ê R3, let c(u, v) be the linear form given by c(u, v)(w) = det(u, v, w) for all w ¡Ê R3 . 
Since 
u1 u2 
u2 u3 
u3 = w1(u2v3 . u3v2)+ w2(u3v1 . u1v3)+ w3(u1v2 . u2v1), we see that the unique vector z ¡Ê R3 such that z ¡¤ w = c(u, v)(w) = det(u, v, w) for all w ¡Ê R3 
LLLL 
LLLLLL 
LLLLLL 

v1 w1 
LLLL LLLu1 L 
u1 v1 v1
. w2 L + w3
L u2 
LLLLLLLLLL
det(u, v, w)= 

v2 
v2 w2 
= w1 
v3 u3 v3 v2 
v3 w3 
is the vector 

.
. 
.
u2v3 . u3v2 
..uvuvz = .3113 
u1v2 . u2v1 
This is just the cross-product u ¡Á v of u and v. Since det(u, v, u) = det(u, v, v)=0, we see that u¡Áv is orthogonal to both u and v. The above allows us to generalize the cross-product to Rn . Given any n . 1 vectors u1,...,un.1 ¡Ê Rn, the cross-product u1 ¡Á¡¤ ¡¤¡¤¡Á un.1 is the unique vector in Rn such that 
(u1 ¡Á¡¤ ¡¤¡¤¡Á un.1) ¡¤ w = det(u1,...,un.1,w) for all w ¡Ê Rn . 
Example 11.8. Consider the vector space Mn(R) of real n ¡Á n matrices with the inner product 
JA, B¡¢ = tr(ATB). 
Let s:Mn(R) ¡ú R be the function given by 
n
ns(A)= aij, i,j=1 
where A =(aij). It is immediately veri.ed that s is a linear form. It is easy to check that the unique matrix Z such that JZ, A¡¢ = s(A) for all A ¡Ê Mn(R) is the matrix Z = ones(n, n) whose entries are all equal to 1. 

11.3. ADJOINT OF A LINEAR MAP 


11.3 Adjoint of a Linear Map 
The existence of the isomorphism b: E ¡ú E. is crucial to the existence of adjoint maps. The importance of adjoint maps stems from the fact that the linear maps arising in physical problems are often self-adjoint, which means that f = f. . Moreover, self-adjoint maps can be diagonalized over orthonormal bases of eigenvectors. This is the key to the solution of many problems in mechanics and engineering in general (see Strang [63]). 
Let E be a Euclidean space of .nite dimension n, and let f : E ¡ú E be a linear map. For every u ¡Ê E, the map v ¡ú u ¡¤ f(v) 
is clearly a linear form in E., and by Theorem 11.6, there is a unique vector in E denoted by f.(u) such that 
f . (u) ¡¤ v = u ¡¤ f(v), 
for every v ¡Ê E. The following simple proposition shows that the map f. is linear. 
Proposition 11.8. Given a Euclidean space E of .nite dimension, for every linear map 
f : E ¡ú E, there is a unique linear map f. : E ¡ú E such that f . (u) ¡¤ v = u ¡¤ f(v), for all u, v ¡Ê E. Proof. Given u1,u2 ¡Ê E, since the inner product is bilinear, we have (u1 + u2) ¡¤ f(v)= u1 ¡¤ f(v)+ u2 ¡¤ f(v), for all v ¡Ê E, and (f . (u1)+ f . (u2)) ¡¤ v = f . (u1) ¡¤ v + f . (u2) ¡¤ v, for all v ¡Ê E, and since by assumption, f . (u1) ¡¤ v = u1 ¡¤ f(v) and f . (u2) ¡¤ v = u2 ¡¤ f(v), for all v ¡Ê E. Thus we get (f . (u1)+ f . (u2)) ¡¤ v =(u1 + u2) ¡¤ f(v)= f . (u1 + u2) ¡¤ v, for all v ¡Ê E. Since b is bijective, this implies that f . (u1 + u2)= f . (u1)+ f . (u2). Similarly, (¦Ëu) ¡¤ f(v)= ¦Ë(u ¡¤ f(v)), 
for all v ¡Ê E, and (¦Ëf . (u)) ¡¤ v = ¦Ë(f . (u) ¡¤ v), for all v ¡Ê E, and since by assumption, 
f . (u) ¡¤ v = u ¡¤ f(v), 
for all v ¡Ê E, we get 
(¦Ëf . (u)) ¡¤ v = ¦Ë(u ¡¤ f(v)) = (¦Ëu) ¡¤ f(v)= f . (¦Ëu) ¡¤ v 
for all v ¡Ê E. Since b is bijective, this implies that 
f . (¦Ëu)= ¦Ëf . (u). 
Thus, f. is indeed a linear map, and it is unique since b is a bijection. 
De.nition 11.4. Given a Euclidean space E of .nite dimension, for every linear map 

f : E ¡ú E, the unique linear map f. : E ¡ú E such that 
f . (u) ¡¤ v = u ¡¤ f(v), for all u, v ¡Ê E 
given by Proposition 11.8 is called the adjoint of f (w.r.t. to the inner product). Linear maps f : E ¡ú E such that f = f. are called self-adjoint maps. 
Self-adjoint linear maps play a very important role because they have real eigenvalues, and because orthonormal bases arise from their eigenvectors. Furthermore, many physical problems lead to self-adjoint linear maps (in the form of symmetric matrices). 
Remark: Proposition 11.8 still holds if the inner product on E is replaced by a nondegen-erate symmetric bilinear form .. 
Linear maps such that f.1 = f., or equivalently 
f . . f = f . f . = id, 
also play an important role. They are linear isometries, or isometries. Rotations are special kinds of isometries. Another important class of linear maps are the linear maps satisfying the property 
f . . f = f . f . , 
called normal linear maps. We will see later on that normal maps can always be diagonalized over orthonormal bases of eigenvectors, but this will require using a Hermitian inner product (over C). 
Given two Euclidean spaces E and F , where the inner product on E is denoted by J.,.¡¢1 and the inner product on F is denoted by J.,.¡¢2, given any linear map f : E ¡ú F , it is immediately veri.ed that the proof of Proposition 11.8 can be adapted to show that there is a unique linear map f. : F ¡ú E such that 
Jf(u),v¡¢2 = Ju, f . (v)¡¢1 
for all u ¡Ê E and all v ¡Ê F . The linear map f. is also called the adjoint of f. The following properties immediately follow from the de.nition of the adjoint map: 
11.4. EXISTENCE AND CONSTRUCTION OF ORTHONORMAL BASES 
(1) For any linear map f : E ¡ú F , we have 
f .. 
= f. 
(2) For any two linear maps f, g : E ¡ú F and any scalar ¦Ë ¡Ê R: 

(f + g) . = f . + g . (¦Ëf) . = ¦Ëf . . 

(3) If 	E, F, G are Euclidean spaces with respective inner products J.,.¡¢1, J.,.¡¢2, and J.,.¡¢3, and if f : E ¡ú F and g : F ¡ú G are two linear maps, then 


(g . f) . = f . . g . . 
Remark: Given any basis for E and any basis for F , it is possible to characterize the matrix of the adjoint f. of f in terms of the matrix of f and the Gram matrices de.ning the inner products; see Problem 11.5. We will do so with respect to orthonormal bases in Proposition 11.14(2). Also, since inner products are symmetric, the adjoint f. of f is also characterized by 
f(u) ¡¤ v = u ¡¤ f . (v), 
for all u, v ¡Ê E. 


11.4 	Existence and Construction of Orthonormal Bases 
We can also use Theorem 11.6 to show that any Euclidean space of .nite dimension has an orthonormal basis. 
Proposition 11.9. Given any nontrivial Euclidean space E of .nite dimension n ¡Ý 1, there is an orthonormal basis (u1,...,un) for E. 
Proof. We proceed by induction on n. When n = 1, take any nonnull vector v ¡Ê E, which exists since we assumed E nontrivial, and let 
v 
u = . 
lvl 
If n ¡Ý 2, again take any nonnull vector v ¡Ê E, and let v 
u1 = . 
lvl 
Consider the linear form .u1 associated with u1. Since u1 = 0, by Theorem 11.6, the linear form .u1 is nonnull, and its kernel is a hyperplane H. Since .u1 (w) = 0i. u1 ¡¤ w = 0, the hyperplane H is the orthogonal complement of {u1}. Furthermore, since u1 = 0 and the inner product is positive de.nite, u1 ¡¤ u1 = 0, and thus, u1 ¡Ê/H, which implies that E = H ¨’ Ru1. However, since E is of .nite dimension n, the hyperplane H has dimension n . 1, and by the induction hypothesis, we can .nd an orthonormal basis (u2,...,un) for H. Now because H and the one dimensional space Ru1 are orthogonal and E = H ¨’ Ru1, it is clear that (u1,...,un) is an orthonormal basis for E. 
As a consequence of Proposition 11.9, given any Euclidean space of .nite dimension n, if (e1,...,en) is an orthonormal basis for E, then for any two vectors u = u1e1 + ¡¤¡¤¡¤ + unen and v = v1e1 + ¡¤¡¤¡¤ + vnen, the inner product u ¡¤ v is expressed as 
n
n 
u ¡¤ v =(u1e1 + ¡¤¡¤¡¤ + unen) ¡¤ (v1e1 + ¡¤¡¤¡¤ + vnen)= uivi, 
i=1 
and the norm lul as 
n1/2
n 
lul = lu1e1 + ¡¤¡¤¡¤ + unenl = u 2 i . 
i=1 
The fact that a Euclidean space always has an orthonormal basis implies that any Gram matrix G can be written as G = QTQ, 
for some invertible matrix Q. Indeed, we know that in a change of basis matrix, a Gram matrix G becomes G ' = P TGP . If the basis corresponding to G ' is orthonormal, then G ' = I, so G =(P .1)TP .1 . 
There is a more constructive way of proving Proposition 11.9, using a procedure known as the Gram¨CSchmidt orthonormalization procedure. Among other things, the Gram¨CSchmidt orthonormalization procedure yields the QR-decomposition for matrices, an important tool in numerical methods. 
Proposition 11.10. Given any nontrivial Euclidean space E of .nite dimension n ¡Ý 1, from any basis (e1,...,en) for E we can construct an orthonormal basis (u1,...,un) for E, with the property that for every k, 1 ¡Ü k ¡Ü n, the families (e1,...,ek) and (u1,...,uk) generate the same subspace. 
Proof. We proceed by induction on n. For n = 1, let 
e1 
u1 = . 
le1l 
For n ¡Ý 2, we also let e1 
u1 = ,
le1l
11.4. EXISTENCE AND CONSTRUCTION OF ORTHONORMAL BASES 
and assuming that (u1,...,uk) is an orthonormal system that generates the same subspace as (e1,...,ek), for every k with 1 ¡Ü k<n, we note that the vector 
k
n 
' 
uk+1 = ek+1 . (ek+1 ¡¤ ui) ui i=1 
is nonnull, since otherwise, because (u1,...,uk) and (e1,...,ek) generate the same subspace, (e1,...,ek+1) would be linearly dependent, which is absurd, since (e1,..., en) is a basis. Thus, the norm of the vector u ' k+1 being nonzero, we use the following construction of the vectors uk and u ' k: 
' u1 ' 
u = e1,u1 = ,
1 lu ' 1l
and for the inductive step 
k
n ' 
u
' k+1 
u = ek+1 . (ek+1 ¡¤ ui) ui,uk+1 = ,
k+1 ' 
luk+1l
i=1 
where 1 ¡Ü k ¡Ü n . 1. It is clear that luk+1l = 1, and since (u1,...,uk) is an orthonormal system, we have 
u ' k+1 ¡¤ ui = ek+1 ¡¤ ui . (ek+1 ¡¤ ui)ui ¡¤ ui = ek+1 ¡¤ ui . ek+1 ¡¤ ui =0, 
for all i with 1 ¡Ü i ¡Ü k. This shows that the family (u1,...,uk+1) is orthonormal, and since (u1,...,uk) and (e1,...,ek) generates the same subspace, it is clear from the de.nition of uk+1 that (u1,...,uk+1) and (e1,...,ek+1) generate the same subspace. This completes the induction step and the proof of the proposition. 
Note that uk' +1 is obtained by subtracting from ek+1 the projection of ek+1 itself onto the orthonormal vectors u1,...,uk that have already been computed. Then uk' +1 is normalized. 
Example 11.9. For a speci.c example of this procedure, let E = R3 with the standard Euclidean norm. Take the basis 
.  .  .  .  .  .  
1  1  1  
e1 = .1.  e2 = .0.  e3 = .1. .  
1  1  0  

Then .. 1
1 
..
u1 = ¡Ì 1 , 
3 
1 and .... .. 
11 1 u ' 2 = e2 . (e2 ¡¤ u1)u1 = .0. . 2 .1. =1 ..2. . 
33
11 1 This implies that 
.. 
1
1 
..
u2 = ¡Ì.2 , 
6 
1 and that 
...... .. 
111 1
21 1 
u ' 3 = e3 . (e3 ¡¤ u1)u1 . (e3 ¡¤ u2)u2 = .1. . .1. + ..2. = . 0 .. 
36 2
01 1 .1 To complete the orthonormal basis, normalize u3 ' to obtain 
.. 
1
1 
..
u3 = ¡Ì 0 . 
2 
.1 An illustration of this example is provided by Figure 11.4. 

Figure 11.4: The top .gure shows the construction of the blue u2 ' as perpendicular to the orthogonal projection of e2 onto u1, while the bottom .gure shows the construction of the green u3 ' as normal to the plane determined by u1 and u2. 
Remarks: 
(1) 
The 	QR-decomposition can now be obtained very easily, but we postpone this until Section 11.6. 

11.4. EXISTENCE AND CONSTRUCTION OF ORTHONORMAL BASES 

(2) 
The proof of Proposition 11.10 also works for a countably in.nite basis for E, producing a countably in.nite orthonormal basis. 


It should also be said that the Gram¨CSchmidt orthonormalization procedure that we have presented is not very stable numerically, and instead, one should use the modi.ed Gram¨C Schmidt method. To compute u ' k+1, instead of projecting ek+1 onto u1,...,uk in a single 
k+1 k+1 k+1
step, it is better to perform k projections. We compute u ,u ,...,u as follows: 
12 k 
u k+1 = ek+1 . (ek+1 ¡¤ u1) u1,
1 
k+1 k+1 k+1 
u = u . (u ¡¤ ui+1) ui+1,
i+1 ii 
' k+1
where 1 ¡Ü i ¡Ü k . 1. It is easily shown that u = u .
k+1 k 
Example 11.10. Let us apply the modi.ed Gram¨CSchmidt method to the (e1,e2,e3) basis of Example 11.9. The only change is the computation of u3' . For the modi.ed Gram¨CSchmidt procedure, we .rst calculate 
.... .. 
11 1 u13 = e3 . (e3 ¡¤ u1)u1 = .1. . 2 .1. =1 . 1 . . 
33
01 .2 
Then .... .. 
111
111
33 3 
u = u . (u ¡¤ u2)u2 = . 1 . + ..2. = . 0 . ,
21 1 
362
.21 .1 
and observe that u32 = u ' 3. See Figure 11.5. 
The following Matlab program implements the modi.ed Gram¨CSchmidt procedure. 
function q = gramschmidt4(e) 
n = size(e,1); 
fori =1:n 

q(:,i) = e(:,i); 
for j = 1:i-1 
r = q(:,j)¡¯*q(:,i); 
q(:,i) = q(:,i) -r*q(:,j); 
end 
r = sqrt(q(:,i)¡¯*q(:,i)); 
q(:,i) = q(:,i)/r; end end 
If we apply the above function to the matrix 
.. 
111 
..
101 , 
110 


Figure 11.5: The top .gure shows the construction of the blue u31 as perpendicular to the orthogonal projection of e3 onto u1, while the bottom .gure shows the construction of the sky blue u23 as perpendicular to the orthogonal projection of u13 onto u2. 
the ouput is the matrix  .  .  
0.5774  0.4082  0.7071  
.0.5774  .0.8165  .0.0000. ,  
0.5774  0.4082  .0.7071  

which matches the result of Example 11.9. 
Example 11.11. If we consider polynomials and the inner product 
1 
Jf, g¡¢ = f(t)g(t)dt, 
.1 
applying the Gram¨CSchmidt orthonormalization procedure to the polynomials 
1, x, x 2 ,...,x n ,..., 
which form a basis of the polynomials in one variable with real coe.cients, we get a family of orthonormal polynomials Qn(x) related to the Legendre polynomials. 
The Legendre polynomials Pn(x) have many nice properties. They are orthogonal, but their norm is not always 1. The Legendre polynomials Pn(x) can be de.ned as follows. Letting fn be the function 
fn(x)=(x 2 . 1)n , 
11.4. EXISTENCE AND CONSTRUCTION OF ORTHONORMAL BASES 
we de.ne Pn(x) as follows: 
1 
f(n)
P0(x)=1, and (x)= (x),
Pnn
2n
n! 
(n)
where fn is the nth derivative of fn. They can also be de.ned inductively as follows: 
P0(x)=1, 
P1(x)= x, 2n +1 n 
Pn+1(x)= xPn(x) . Pn.1(x). 
n +1 n +1 
Here is an explicit summation fo Pn(x): 
ln/2.
1 n n 2n . 2k 
Pn(x)= (.1)k x n.2k . 
2n kn 
k=0 
The polynomials Qn are related to the Legendre polynomials Pn as follows: 
2n +1 
Qn(x)= Pn(x). 
2 
Example 11.12. Consider polynomials over [.1, 1], with the symmetric bilinear form 
1 
1 
Jf, g¡¢ = ¡Ì f(t)g(t)dt. .1 1 . t2 
We leave it as an exercise to prove that the above de.nes an inner product. It can be shown that the polynomials Tn(x) given by 
Tn(x) = cos(n arccos x),n ¡Ý 0, 
(equivalently, with x = cos ¦È, we have Tn(cos ¦È) = cos(n¦È)) are orthogonal with respect to the above inner product. These polynomials are the Chebyshev polynomials. Their norm is not equal to 1. Instead, we have 
 
JTn,Tn¡¢ =¦Ð 2 if n> 0, 
¦Ð if n =0. 
Using the identity (cos ¦È + i sin ¦È)n = cos n¦È + i sin n¦È and the binomial formula, we obtain the following expression for Tn(x): 
ln/2.
n 
n 
Tn(x)= (x 2 . 1)k x n.2k . 
2k 
k=0 
The Chebyshev polynomials are de.ned inductively as follows: 
T0(x)=1 
T1(x)= x 
Tn+1(x)=2xTn(x) . Tn.1(x),n ¡Ý 1. 
Using these recurrence equations, we can show that 
¡Ì¡Ì (x . x2 . 1)n +(x + x2 . 1)n 
Tn(x)= . 
2 
The polynomial Tn has n distinct roots in the interval [.1, 1]. The Chebyshev polynomials play an important role in approximation theory. They are used as an approximation to a best polynomial approximation of a continuous function under the sup-norm (¡Þ-norm). 
The inner products of the last two examples are special cases of an inner product of the form 
1 
Jf, g¡¢ = W (t)f(t)g(t)dt, 
.1 
where W (t) is a weight function. If W is a nonzero continuous function such that W (x) ¡Ý 0 on (.1, 1), then the above bilinear form is indeed positive de.nite. Families of orthogonal polynomials used in approximation theory and in physics arise by a suitable choice of the weight function W . Besides the previous two examples, the Hermite polynomials correspond to W (x)= e.x2 , the Laguerre polynomials to W (x)= e.x , and the Jacobi polynomials to W (x) = (1 . x)¦Á(1 + x)¦Â, with ¦Á,¦Â > .1. Comprehensive treatments of orthogonal polynomials can be found in Lebedev [45], Sansone [53], and Andrews, Askey and Roy [1]. 
We can also prove the following proposition regarding orthogonal spaces. 
Proposition 11.11. Given any nontrivial Euclidean space E of .nite dimension n ¡Ý 1, for any subspace F of dimension k, the orthogonal complement F ¡Í of F has dimension n . k, and E = F ¨’ F ¡Í . Furthermore, we have F ¡Í¡Í = F . 
Proof. From Proposition 11.9, the subspace F has some orthonormal basis (u1,...,uk). This linearly independent family (u1,...,uk) can be extended to a basis (u1,...,uk, vk+1,...,vn), and by Proposition 11.10, it can be converted to an orthonormal basis (u1,...,un), which contains (u1,...,uk) as an orthonormal basis of F . Now any vector w = w1u1+¡¤¡¤¡¤+wnun ¡Ê E is orthogonal to F i. w ¡¤ ui = 0, for every i, where 1 ¡Ü i ¡Ü k, i. wi = 0 for every i, where 1 ¡Ü i ¡Ü k. Clearly, this shows that (uk+1,...,un) is a basis of F ¡Í, and thus E = F ¨’F ¡Í, and F ¡Í has dimension n . k. Similarly, any vector w = w1u1 + ¡¤¡¤¡¤ + wnun ¡Ê E is orthogonal to F ¡Í i. w ¡¤ ui = 0, for every i, where k +1 ¡Ü i ¡Ü n, i. wi = 0 for every i, where k +1 ¡Ü i ¡Ü n. Thus, (u1,...,uk) is a basis of F ¡Í¡Í, and F ¡Í¡Í = F . 
11.5. LINEAR ISOMETRIES (ORTHOGONAL TRANSFORMATIONS) 
11.5 Linear Isometries (Orthogonal Transformations) 
In this section we consider linear maps between Euclidean spaces that preserve the Euclidean norm. These transformations, sometimes called rigid motions, play an important role in geometry. 
De.nition 11.5. Given any two nontrivial Euclidean spaces E and F of the same .nite dimension n, a function f : E ¡ú F is an orthogonal transformation, or a linear isometry, if it is linear and 
lf(u)l = lul, for all u ¡Ê E. 
Remarks: 
(1) A linear isometry is often de.ned as a linear map such that 
lf(v) . f(u)l = lv . ul, 
for all u, v ¡Ê E. Since the map f is linear, the two de.nitions are equivalent. The second de.nition just focuses on preserving the distance between vectors. 
(2) Sometimes, a linear map satisfying the condition of De.nition 11.5 is called a metric map, and a linear isometry is de.ned as a bijective metric map. 
An isometry (without the word linear) is sometimes de.ned as a function f : E ¡ú F (not necessarily linear) such that lf(v) . f(u)l = lv . ul, 
for all u, v ¡Ê E, i.e., as a function that preserves the distance. This requirement turns out to be very strong. Indeed, the next proposition shows that all these de.nitions are equivalent when E and F are of .nite dimension, and for functions such that f(0) = 0. 
Proposition 11.12. Given any two nontrivial Euclidean spaces E and F of the same .nite dimension n, for every function f : E ¡ú F , the following properties are equivalent: 
(1) 
f is a linear map and lf(u)l = lul, for all u ¡Ê E; 

(2) 
lf(v) . f(u)l = lv . ul, for all u, v ¡Ê E, and f(0) = 0; 

(3) 
f(u) ¡¤ f(v)= u ¡¤ v, for all u, v ¡Ê E. 


Furthermore, such a map is bijective. 
Proof. Clearly, (1) implies (2), since in (1) it is assumed that f is linear. Assume that (2) holds. In fact, we shall prove a slightly stronger result. We prove that if lf(v) . f(u)l = lv . ul 
for all u, v ¡Ê E, then for any vector ¦Ó ¡Ê E, the function g : E ¡ú F de.ned such that 

g(u)= f(¦Ó + u) . f(¦Ó) 

for all u ¡Ê E is a linear map such that g(0) = 0 and (3) holds. Clearly, g(0) = f(¦Ó).f(¦Ó) = 0. Note that from the hypothesis 
lf(v) . f(u)l = lv . ul 
for all u, v ¡Ê E, we conclude that 
lg(v) . g(u)l = lf(¦Ó + v) . f(¦Ó) . (f(¦Ó + u) . f(¦Ó))l, = lf(¦Ó + v) . f(¦Ó + u)l, = l¦Ó + v . (¦Ó + u)l, = lv . ul, 
for all u, v ¡Ê E. Since g(0) = 0, by setting u =0 in 
lg(v) . g(u)l = lv . ul, 
we get lg(v)l = lvl 
for all v ¡Ê E. In other words, g preserves both the distance and the norm. To prove that g preserves the inner product, we use the simple fact that 
2u ¡¤ v = lul2 + lvl2 .lu . vl2 
for all u, v ¡Ê E. Then since g preserves distance and norm, we have 
2g(u) ¡¤ g(v)= lg(u)l2 + lg(v)l2 .lg(u) . g(v)l2 
= lul2 + lvl2 .lu . vl2 
=2u ¡¤ v, 
and thus g(u) ¡¤ g(v)= u ¡¤ v, for all u, v ¡Ê E, which is (3). In particular, if f(0) = 0, by letting ¦Ó = 0, we have g = f, and f preserves the scalar product, i.e., (3) holds. 
Now assume that (3) holds. Since E is of .nite dimension, we can pick an orthonormal basis (e1,...,en) for E. Since f preserves inner products, (f(e1),..., f(en)) is also orthonor-mal, and since F also has dimension n, it is a basis of F . Then note that since (e1,...,en) and (f(e1),...,f(en)) are orthonormal bases, for any u ¡Ê E we have 
nn
nn 
u =(u ¡¤ ei)ei = uiei i=1 i=1 
11.5. LINEAR ISOMETRIES (ORTHOGONAL TRANSFORMATIONS) 
and 
n
f(u)= (f(u) ¡¤ f(ei))f(ei), 
i=1 
and since f preserves inner products, this shows that 
n 
nnn
f(u)= (f(u) ¡¤ f(ei))f(ei)= (u ¡¤ ei)f(ei)= uif(ei), 
i=1 i=1 i=1 

n 
which proves that f is linear. Obviously, f preserves the Euclidean norm, and (3) implies (1). 
Finally, if f(u)= f(v), then by linearity f(v . u) = 0, so that lf(v . u)l = 0, and since f preserves norms, we must have lv . ul = 0, and thus u = v. Thus, f is injective, and since E and F have the same .nite dimension, f is bijective. 
n 
Remarks: 
(i) 
The dimension assumption is needed only to prove that (3) implies (1) when f is not known to be linear, and to prove that f is surjective, but the proof shows that (1) implies that f is injective. 

(ii) 
The implication that (3) implies (1) holds if we also assume that f is surjective, even if E has in.nite dimension. 


n 
In (2), when f does not satisfy the condition f(0) = 0, the proof shows that f is an a.ne map. Indeed, taking any vector ¦Ó as an origin, the map g is linear, and 
f(¦Ó + u)= f(¦Ó)+ g(u) for all u ¡Ê E. 
By Proposition 5.15, this shows that f is a.ne with associated linear map g. This fact is worth recording as the following proposition. 
Proposition 11.13. Given any two nontrivial Euclidean spaces E and F of the same .nite dimension n, for every function f : E ¡ú F , if 
lf(v) . f(u)l = lv . ul for all u, v ¡Ê E, 
then f is an a.ne map, and its associated linear map g is an isometry. 
In view of Proposition 11.12, we usually abbreviate ¡°linear isometry¡± as ¡°isometry,¡± unless we wish to emphasize that we are dealing with a map between vector spaces. 
We are now going to take a closer look at the isometries f : E ¡ú E of a Euclidean space of .nite dimension. 
11.6 The Orthogonal Group, Orthogonal Matrices 
In this section we explore some of the basic properties of the orthogonal group and of orthogonal matrices. 
Proposition 11.14. Let E be any Euclidean space of .nite dimension n, and let f : E ¡ú E be any linear map. The following properties hold: 
(1) 
The linear map f : E ¡ú E is an isometry i. 
f . f . = f . . f = id. 


(2) 
For every orthonormal basis (e1,...,en) of E, if the matrix of f is A, then the matrix of f. is the transpose AT of A, and f is an isometry i. A satis.es the identities 


AAT = ATA = In, 
where In denotes the identity matrix of order n, i. the columns of A form an orthonor-mal basis of Rn, i. the rows of A form an orthonormal basis of Rn . 
Proof. (1) The linear map f : E ¡ú E is an isometry i. f(u) ¡¤ f(v)= u ¡¤ v, for all u, v ¡Ê E, i. f . (f(u)) ¡¤ v = f(u) ¡¤ f(v)= u ¡¤ v for all u, v ¡Ê E, which implies (f . (f(u)) . u) ¡¤ v =0 for all u, v ¡Ê E. Since the inner product is positive de.nite, we must have f . (f(u)) . u =0 for all u ¡Ê E, that is, f . . f = id. But an endomorphism f of a .nite-dimensional vector space that has a left inverse is an isomorphism, so f . f. = id. The converse is established by doing the above steps backward. 
(2) If (e1,...,en) is an orthonormal basis for E, let A =(aij) be the matrix of f, and let B =(bij) be the matrix of f. . Since f. is characterized by 
f . (u) ¡¤ v = u ¡¤ f(v) 
for all u, v ¡Ê E, using the fact that if w = w1e1 + ¡¤¡¤¡¤ + wnen we have wk = w ¡¤ ek for all k, 1 ¡Ü k ¡Ü n, letting u = ei and v = ej, we get 
bji = f . (ei) ¡¤ ej = ei ¡¤ f(ej)= aij, 
11.6. THE ORTHOGONAL GROUP, ORTHOGONAL MATRICES 
for all i, j,1 ¡Ü i, j ¡Ü n. Thus, B = AT . Now if X and Y are arbitrary matrices over the basis (e1,...,en), denoting as usual the jth column of X by Xj, and similarly for Y , a simple calculation shows that 
XTY =(Xi ¡¤ Y j)1¡Üi,j¡Ün. Then it is immediately veri.ed that if X = Y = A, then 
ATA = AAT = In 
i. the column vectors (A1,...,An) form an orthonormal basis. Thus, from (1), we see that 
(2) is clear (also because the rows of A are the columns of AT). 
Proposition 11.14 shows that the inverse of an isometry f is its adjoint f. . Recall that the set of all real n ¡Á n matrices is denoted by Mn(R). Proposition 11.14 also motivates the following de.nition. 
De.nition 11.6. A real n ¡Á n matrix is an orthogonal matrix if 
AAT = ATA = In. 
Remark: It is easy to show that the conditions AAT = In, ATA = In, and A.1 = AT, are equivalent. Given any two orthonormal bases (u1,...,un) and (v1,...,vn), if P is the change of basis matrix from (u1,...,un) to (v1,...,vn), since the columns of P are the coordinates of the vectors vj with respect to the basis (u1,...,un), and since (v1,...,vn) is orthonormal, the columns of P are orthonormal, and by Proposition 11.14 (2), the matrix P is orthogonal. 
The proof of Proposition 11.12 (3) also shows that if f is an isometry, then the image of an orthonormal basis (u1,...,un) is an orthonormal basis. Students often ask why orthogonal matrices are not called orthonormal matrices, since their columns (and rows) are orthonormal bases! I have no good answer, but isometries do preserve orthogonality, and orthogonal matrices correspond to isometries. 
Recall that the determinant det(f) of a linear map f : E ¡ú E is independent of the choice of a basis in E. Also, for every matrix A ¡Ê Mn(R), we have det(A) = det(AT), and for any two n ¡Á n matrices A and B, we have det(AB) = det(A) det(B). Then if f is an isometry, and A is its matrix with respect to any orthonormal basis, AAT = ATA = In implies that det(A)2 = 1, that is, either det(A) = 1, or det(A)= .1. It is also clear that the isometries of a Euclidean space of dimension n form a group, and that the isometries of determinant +1 form a subgroup. This leads to the following de.nition. 
De.nition 11.7. Given a Euclidean space E of dimension n, the set of isometries f : E ¡ú E forms a subgroup of GL(E) denoted by O(E), or O(n) when E = Rn, called the orthogonal group (of E). For every isometry f, we have det(f)= ¡À1, where det(f) denotes the deter-minant of f. The isometries such that det(f) = 1 are called rotations, or proper isometries, or proper orthogonal transformations, and they form a subgroup of the special linear group SL(E) (and of O(E)), denoted by SO(E), or SO(n) when E = Rn, called the special or-thogonal group (of E). The isometries such that det(f)= .1 are called improper isometries, or improper orthogonal transformations, or .ip transformations. 
11.7 The Rodrigues Formula 
When n = 3 and A is a skew symmetric matrix, it is possible to work out an explicit formula for eA . For any 3 ¡Á 3 real skew symmetric matrix 
.. 
0 .cb A = . c 0 .a. , .ba 0 
¡Ì 
if we let ¦È = a2 + b2 + c2 and 
.. 
2
aab ac 
B = .ab b2 bc. , 2
acbc c
then we have the following result known as Rodrigues¡¯ formula (1840). The (real) vector space of n ¡Á n skew symmetric matrices is denoted by so(n). 
Proposition 11.15. The exponential map exp: so(3) ¡ú SO(3) is given by 
sin ¦È (1 . cos ¦È) 
e A = cos ¦ÈI3 + A + B, 
¦È¦È2 
or, equivalently, by 
sin ¦È (1 . cos ¦È) 
e A = I3 + A + A2 
¦È¦È2 if ¦È =0, with e03 = I3. 
Proof sketch. First observe that A2 = .¦È2I3 + B, 
since 
..... . 
0 .cb 0 .cb .c2 . b2 ba ca 
A2 = . c 0 .a.. c 0 .a. = . ab .c2 . a2 cb . .b2 . a2
.ba 0 .ba 0 ac cb 
. ... 
.a2 . b2 . c2 00 a2 ba ca = . 0 .a2 . b2 . c2 0 . + .ab b2 cb. 00 .a2 . b2 . c2 accb c2 
= .¦È2I3 + B, 
and that AB = BA =0. 
From the above, deduce that 
A3 = .¦È2A, 
11.7. THE RODRIGUES FORMULA 
and for any k ¡Ý 0, 
A4k+1 
= ¦È4kA, A4k+2 
= ¦È4kA2 , A4k+3 = .¦È4k+2A, 
A4k+4 = .¦È4k+2A2 
. 
Then prove the desired result by writing the power series for eA and regrouping terms so that the power series for cos ¦È and sin ¦È show up. In particular 
n Ap n A2p+1 n A2p e A = I3 += I3 ++ 
p! (2p + 1)! (2p)!
p¡Ý1 p¡Ý0 p¡Ý1 
nn (.1)p.1¦È2(p.1)
(.1)p¦È2p 
A2 
= I3 + A + 
(2p + 1)! (2p)!
p¡Ý0 p¡Ý1 n (.1)p¦È2p+1 A2 n
A (.1)p¦È2p 
= I3 + . 
¦È (2p + 1)! ¦È2 (2p)!
p¡Ý0 p¡Ý1 A2 n A2
sin ¦È (.1)p¦È2p 
= I3 + A . + 
¦È¦È2 p¡Ý0 (2p)! ¦È2 
sin ¦È (1 . cos ¦È) 
= I3 + A + A2 ,
¦È¦È2 
as claimed. 
The above formulae are the well-known formulae expressing a rotation of axis speci.ed by the vector (a, b, c) and angle ¦È. 
The Rodrigues formula can used to show that the exponential map exp: so(3) ¡ú SO(3) is surjective. 
Given any rotation matrix R ¡Ê SO(3), we have the following cases: 
(1) 
The case R = I is trivial. 

(2) 
If R = I and tr(R)= .1, then 


L  
Lexp .1(R)= (R . RT ) LL 1 + 2cos ¦È = tr(R). 
¦È 
2 sin ¦È 
(Recall that tr(R)= r11 + r22 + r33, the trace of the matrix R). 
Then there is a unique skew-symmetric B with corresponding ¦È satisfying 0 <¦È<¦Ð 
B
such that e= R. 
(3) If R = I and tr(R)= .1, then R is a rotation by the angle ¦Ð and things are more complicated, but a matrix B can be found. We leave this part as a good exercise: see Problem 16.8. 
The computation of a logarithm of a rotation in SO(3) as sketched above has applications in kinematics, robotics, and motion interpolation. 
As an immediate corollary of the Gram¨CSchmidt orthonormalization procedure, we obtain the QR-decomposition for invertible matrices. 
11.8 QR-Decomposition for Invertible Matrices 
Now that we have the de.nition of an orthogonal matrix, we can explain how the Gram¨C Schmidt orthonormalization procedure immediately yields the QR-decomposition for matri-ces. 
De.nition 11.8. Given any real n ¡Á n matrix A,a QR-decomposition of A is any pair of n ¡Á n matrices (Q, R), where Q is an orthogonal matrix and R is an upper triangular matrix such that A = QR. 
Note that if A is not invertible, then some diagonal entry in R must be zero. 
Proposition 11.16. Given any real n ¡Á n matrix A, if A is invertible, then there is an orthogonal matrix Q and an upper triangular matrix R with positive diagonal entries such that A = QR. 
Proof. We can view the columns of A as vectors A1,...,An in En . If A is invertible, then they are linearly independent, and we can apply Proposition 11.10 to produce an orthonor-mal basis using the Gram¨CSchmidt orthonormalization procedure. Recall that we construct 
'k
vectors Qk and Qas follows: 
'1 
'1 Q
Q= A1 ,Q1 = ,
lQ'1l
and for the inductive step 
k
n 'k+1 
'k+1 = Ak+1 . (Ak+1 Qk+1 Q
Q¡¤ Qi) Qi , = ,
lQ'k+1l
i=1 
where 1 ¡Ü k ¡Ü n . 1. If we express the vectors Ak in terms of the Qi and Q'i, we get the triangular system 
A1 = lQ'1lQ1 , 
. 
. 
. Aj =(Aj ¡¤ Q1) Q1 + ¡¤¡¤¡¤ +(Aj ¡¤ Qi) Qi + ¡¤¡¤¡¤ +(Aj ¡¤ Qj.1) Qj.1 + lQ'jlQj, . 
. 
. 
'
An =(An ¡¤ Q1) Q1 + ¡¤¡¤¡¤ +(An ¡¤ Qn.1) Qn.1 + lQnlQn . 
11.8. QR-DECOMPOSITION FOR INVERTIBLE MATRICES 
Letting rkk = lQ ' kl, and rij = Aj ¡¤ Qi (the reversal of i and j on the right-hand side is intentional!), where 1 ¡Ü k ¡Ü n,2 ¡Ü j ¡Ü n, and 1 ¡Ü i ¡Ü j . 1, and letting qij be the ith component of Qj, we note that aij, the ith component of Aj, is given by 
aij = r1 jqi 1 + ¡¤¡¤¡¤ + rijqii + ¡¤¡¤¡¤ + rjjqij = qi 1r1 j + ¡¤¡¤¡¤ + qiirij + ¡¤¡¤¡¤ + qijrjj. 
If we let Q =(qij), the matrix whose columns are the components of the Qj, and R =(rij), the above equations show that A = QR, where R is upper triangular. The diagonal entries rkk = lQ ' kl = Ak ¡¤ Qk are indeed positive. 
The reader should try the above procedure on some concrete examples for 2 ¡Á2 and 3 ¡Á3 matrices. 
Remarks: 
(1) 
Because the diagonal entries of R are positive, it can be shown that Q and R are unique. More generally, if A is invertible and if A = Q1R1 = Q2R2 are two QR-decompositions for A, then 

R1R.21 = QT 1 Q2. The matrix QT 1 Q2 is orthogonal and it is easy to see that R1R2 .1 is upper triangular. But an upper triangular matrix which is orthogonal must be a diagonal matrix D with diagonal entries ¡À1, so Q2 = Q1D and R2 = DR1. 

(2) 
The QR-decomposition holds even when A is not invertible. In this case, R has some zero on the diagonal. However, a di.erent proof is needed. We will give a nice proof using Householder matrices (see Proposition 12.4, and also Strang [63, 64], Golub and Van Loan [30], Trefethen and Bau [68], Demmel [16], Kincaid and Cheney [39], or Ciarlet [14]). 


For better numerical stability, it is preferable to use the modi.ed Gram¨CSchmidt method to implement the QR-factorization method. Here is a Matlab program implementing QR-factorization using modi.ed Gram¨CSchmidt. 
function [Q,R] = qrv4(A) n = size(A,1); fori =1:n 
Q(:,i) = A(:,i); 
for j = 1:i-1 
R(j,i) = Q(:,j)¡¯*Q(:,i); 
Q(:,i) = Q(:,i) -R(j,i)*Q(:,j); 
end 
R(i,i) = sqrt(Q(:,i)¡¯*Q(:,i)); 
Q(:,i) = Q(:,i)/R(i,i); end end 
Example 11.13. Consider the matrix 
.. 
005 
..
A = 	041 . 111 
To determine the QR-decomposition of A, we .rst use the Gram-Schmidt orthonormalization procedure to calculate Q =(Q1Q2Q3). By de.nition 
.. 
0 
A1 = Q '1 = Q1 = .	0. , 1 
.. 
0 
and since A2 = .	4., we discover that 1 
.... .	. 
000 Q '2 = A2 . (A2 ¡¤ Q1)Q1 = .4. . .0. = .4.. 110 
.. 
0 
Hence, Q2 = .	1.. Finally, 0 
.....	. .. 
5005 
Q '3 = A3 . (A3 ¡¤ Q1)Q1 . (A3 ¡¤ Q2)Q2 = .	1. . .0. . .1. = .0., 1100 
.. 
1 
which implies that Q3 = .	0.. According to Proposition 11.16, in order to determine R we 0 
need to calculate 
  
r11 = Q '1 =1 r12 = A2 ¡¤ Q1 =1 r13 = A3 ¡¤ Q1 =1 
  
r22 = Q '2 =4 	r23 = A3 ¡¤ Q2 =1 
  
r33 = Q '3 =5. 
.. 
00	5 
In summary, we have found that the QR-decomposition of A = .	041. is 111 
.. 	.. 
001 	111 
.. 	..
Q = 	010 and R = 041 . 100 005 
11.8. QR-DECOMPOSITION FOR INVERTIBLE MATRICES Example 11.14. Another example of QR-decomposition is 
..
.
¡Ì¡Ì

.
.

¡Ì ¡Ì¡Ì

. 

112 1/ 21/ 20 21/ 22

¡Ì 

¡Ì 

A =

.

001

.

=

.

0 01
¡Ì¡Ì 
.
.

01/ 22

.

. 

100 1/ 2 .1/ 20 001 Example 11.15. If we apply the above Matlab function to the matrix 
.
. 
A = 

..... 

4  1  0  0  0  
1  4  1  0  0  
0  1  4  1  0  

00141 
00014 

..... 

, 

we obtain 

.
. 
0.9701 .0.2339 0.0619 .0.0166 0.0046 0.2425 0.9354 .0.2477 0.0663 .0.0184 00.2650 0.9291 .0.2486 0.0691 
0 00.2677 0.9283 .0.2581 
00 00.2679 0.9634 

..... 

..... 

Q = 

and

.
. 
..... 

4.1231 1.9403 0.24250 0 03.7730 1.9956 0.2650 0 
0 03.7361 1.9997 0.2677 

0 0 073.7324 2.0000 
000 03.5956 

..... 

R = 

. 

Remark: The Matlab function qr, called by [Q, R] = qr(A), does not necessarily return an upper-triangular matrix whose diagonal entries are positive. 
The QR-decomposition yields a rather e.cient and numerically stable method for solving systems of linear equations. Indeed, given a system Ax = b, where A is an n ¡Á n invertible matrix, writing A = QR, since Q is orthogonal, we get 
Rx = QTb, 
and since R is upper triangular, we can solve it by Gaussian elimination, by solving for the last variable xn .rst, substituting its value into the system, then solving for xn.1, etc. The QR-decomposition is also very useful in solving least squares problems (we will come back to this in Chapter 21), and for .nding eigenvalues; see Chapter 17. It can be easily adapted to the case where A is a rectangular m ¡Á n matrix with independent columns (thus, n ¡Ü m). In this case, Q is not quite orthogonal. It is an m ¡Á n matrix whose columns are orthogonal, and R is an invertible n¡Án upper triangular matrix with positive diagonal entries. For more on QR, see Strang [63, 64], Golub and Van Loan [30], Demmel [16], Trefethen and Bau [68], or Serre [57]. 
A somewhat surprising consequence of the QR-decomposition is a famous determinantal inequality due to Hadamard. 
Proposition 11.17. (Hadamard) For any real n ¡Á n matrix A =(aij), we have 
nn1/2 nn1/2
nn nn 
| det(A)|¡Ü a 2 and | det(A)|¡Ü a 2 .
ij ij i=1 j=1 j=1 i=1 
Moreover, equality holds i. either A has a zero row in the left inequality or a zero column in the right inequality, or A is orthogonal. 
Proof. If det(A) = 0, then the inequality is trivial. In addition, if the righthand side is also 0, then either some column or some row is zero. If det(A) = 0, then we can factor A as A = QR, with Q is orthogonal and R =(rij) upper triangular with positive diagonal entries. Then since Q is orthogonal det(Q)= ¡À1, so 
n 
| det(A)| = | det(Q)|| det(R)| = rjj. 
j=1 
Now as Q is orthogonal, it preserves the Euclidean norm, so 
nn
nn
2 22
2 22 
a = Aj = QRj = Rj = r ¡Ý r
ij 2 22 ij jj, i=1 i=1 
which implies that 
nnnn1/2
nn nn 
| det(A)| = rjj ¡Ü Rj = aij 2 . 
2 j=1 j=1 j=1 i=1 
The other inequality is obtained by replacing A by AT . Finally, if det(A) = 0 and equality holds, then we must have 
rjj = Aj 2 , 1 ¡Ü j ¡Ü n, 
which can only occur if A is orthogonal. 
Another version of Hadamard¡¯s inequality applies to symmetric positive semide.nite matrices. 
Proposition 11.18. (Hadamard) For any real n ¡Á n matrix A =(aij), if A is symmetric positive semide.nite, then we have 
n
n 
det(A) ¡Ü aii. i=1 
Moreover, if A is positive de.nite, then equality holds i. A is a diagonal matrix. 

11.9. SOME APPLICATIONS OF EUCLIDEAN GEOMETRY 
Proof. If det(A) = 0, the inequality is trivial. Otherwise, A is positive de.nite, and by Theorem 7.10 (the Cholesky Factorization), there is a unique upper triangular matrix B with positive diagonal entries such that 
A = BTB. 
Thus, det(A) = det(BTB) = det(BT) det(B) = det(B)2 . If we apply the Hadamard inequal-ity (Proposition 11.17) to B, we obtain 
det(B) ¡Ü b2 . (.)
ij 
j=1 i=1 However, the diagonal entries ajj of A = BTB are precisely the square norms lBjl2 
nnnn 
a 
1/2 
= 
2 
n b2 
ij, so by squaring (.), we obtain 
i=1 
n
n
nnn
det(A) = det(B)2 ¡Ü b2 = 
n a.jjij 
j=1 i=1 j=1 
If det(A) = 0 and equality holds, then B must be orthogonal, which implies that B is a diagonal matrix, and so is A. 
We derived the second Hadamard inequality (Proposition 11.18) from the .rst (Proposi-tion 11.17). We leave it as an exercise to prove that the .rst Hadamard inequality can be deduced from the second Hadamard inequality. 
11.9 Some Applications of Euclidean Geometry 
Euclidean geometry has applications in computational geometry, in particular Voronoi dia-grams and Delaunay triangulations. In turn, Voronoi diagrams have applications in motion planning (see O¡¯Rourke [49]). 
Euclidean geometry also has applications to matrix analysis. Recall that a real n ¡Á n matrix A is symmetric if it is equal to its transpose AT . One of the most important properties of symmetric matrices is that they have real eigenvalues and that they can be diagonalized by an orthogonal matrix (see Chapter 16). This means that for every symmetric matrix A, there is a diagonal matrix D and an orthogonal matrix P such that 
A = P DP T . 
Even though it is not always possible to diagonalize an arbitrary matrix, there are various decompositions involving orthogonal matrices that are of great practical interest. For exam-ple, for every real matrix A, there is the QR-decomposition, which says that a real matrix A can be expressed as 
A = QR, 
where Q is orthogonal and R is an upper triangular matrix. This can be obtained from the Gram¨CSchmidt orthonormalization procedure, as we saw in Section 11.8, or better, using Householder matrices, as shown in Section 12.2. There is also the polar decomposition, which says that a real matrix A can be expressed as 
A = QS, 
where Q is orthogonal and S is symmetric positive semide.nite (which means that the eigen-values of S are nonnegative). Such a decomposition is important in continuum mechanics and in robotics, since it separates stretching from rotation. Finally, there is the wonderful singular value decomposition, abbreviated as SVD, which says that a real matrix A can be expressed as 
A = V DUT , 
where U and V are orthogonal and D is a diagonal matrix with nonnegative entries (see Chapter 20). This decomposition leads to the notion of pseudo-inverse, which has many applications in engineering (least squares solutions, etc). For an excellent presentation of all these notions, we highly recommend Strang [64, 63], Golub and Van Loan [30], Demmel [16], Serre [57], and Trefethen and Bau [68]. 
The method of least squares, invented by Gauss and Legendre around 1800, is another great application of Euclidean geometry. Roughly speaking, the method is used to solve inconsistent linear systems Ax = b, where the number of equations is greater than the number of variables. Since this is generally impossible, the method of least squares consists in .nding a solution x minimizing the Euclidean norm lAx . bl2, that is, the sum of the squares of the ¡°errors.¡± It turns out that there is always a unique solution x+ of smallest norm minimizing lAx . bl2, and that it is a solution of the square system 
ATAx = ATb, 
called the system of normal equations. The solution x+ can be found either by using the QR-decomposition in terms of Householder transformations, or by using the notion of pseudo-inverse of a matrix. The pseudo-inverse can be computed using the SVD decomposition. Least squares methods are used extensively in computer vision. More details on the method of least squares and pseudo-inverses can be found in Chapter 21. 
11.10 Summary 
The main concepts and results of this chapter are listed below: 
. 
Bilinear forms; positive de.nite bilinear forms. 

. 
Inner products, scalar products, Euclidean spaces. 

. 
Quadratic form associated with a bilinear form. 

11.10. SUMMARY 

. 	
The Euclidean space En . 

. 	
The polar form of a quadratic form. 

. 	
Gram matrix associated with an inner product. 

. 	
The Cauchy¨CSchwarz inequality; the Minkowski inequality. 

. 	
The parallelogram law. 

. 	
Orthogonality, orthogonal complement F ¡Í; orthonormal family. 

. 	
The musical isomorphisms b: E ¡ú E. and : E. ¡ú E (when E is .nite-dimensional); Theorem 11.6. 

. 	
The adjoint of a linear map (with respect to an inner product). 

. 	
Existence of an orthonormal basis in a .nite-dimensional Euclidean space (Proposition 11.9). 

. 	
The Gram¨CSchmidt orthonormalization procedure (Proposition 11.10). 

. 	
The Legendre and the Chebyshev polynomials. 

. 	
Linear isometries (orthogonal transformations, rigid motions). 

. 	
The orthogonal group, orthogonal matrices. 

. 	
The matrix representing the adjoint f. of a linear map f is the transpose of the matrix representing f. 

. 	
The orthogonal group O(n) and the special orthogonal group SO(n). 

. 	
QR-decomposition for invertible matrices. 

. 	
The Hadamard inequality for arbitrary real matrices. 

. 	
The Hadamard inequality for symmetric positive semide.nite matrices. 

. 	
The Rodrigues formula for rotations in SO(3). 


11.11 Problems 
Problem 11.1. E be a vector space of dimension 2, and let (e1,e2) be a basis of E. Prove that if a> 0 and b2 . ac < 0, then the bilinear form de.ned such that .(x1e1 + y1e2,x2e1 + y2e2)= ax1x2 + b(x1y2 + x2y1)+ cy1y2 
is a Euclidean inner product. Problem 11.2. Let C[a, b] denote the set of continuous functions f :[a, b] ¡ú R. Given any two functions f, g ¡ÊC[a, b], let 
b 
Jf, g¡¢ = f(t)g(t)dt. 
a 
Prove that the above bilinear form is indeed a Euclidean inner product. Problem 11.3. Consider the inner product 
¦Ð 
Jf, g¡¢ = f(t)g(t)dt 
.¦Ð 
of Problem 11.2 on the vector space C[.¦Ð, ¦Ð]. Prove that Jsin px, sin qx¡¢ = ¦Ð if p = q, p, q ¡Ý 1, 
0 if p = q, p, q ¡Ý 1, Jcos px, cos qx¡¢ = ¦Ð if p = q, p, q ¡Ý 1, 
0 if p = q, p, q ¡Ý 0, Jsin px, cos qx¡¢ =0, for all p ¡Ý 1 and q ¡Ý 0, and J1, 1¡¢ = ¦Ð dx =2¦Ð.
.¦Ð 
Problem 11.4. Prove that the following matrix is orthogonal and skew-symmetric: 
0  1  1  1  
M = 1 ¡Ì 3 ... .1 .1  0 1  .1 0  1 .1 ....  
.1  .1  1  0  

Problem 11.5. Let E and F be two .nite Euclidean spaces, let (u1,...,un) be a basis of E, and let (v1,...,vm) be a basis of F . For any linear map f : E ¡ú F , if A is the matrix of f w.r.t. the basis (u1,...,un) and B is the matrix of f. w.r.t. the basis (v1,...,vm), if G1 is the Gram matrix of the inner product on E (w.r.t. (u1,...,un)) and if G2 is the Gram matrix of the inner product on F (w.r.t. (v1,...,vm)), then 
B = G.11ATG2. 
.
. 
11.11. PROBLEMS 
Problem 11.6. Let A be an invertible matrix. Prove that if A = Q1R1 = Q2R2 are two QR-decompositions of A and if the diagonal entries of R1 and R2 are positive, then Q1 = Q2 and R1 = R2. 
Problem 11.7. Prove that the .rst Hadamard inequality can be deduced from the second Hadamard inequality. 
Problem 11.8. Let E be a real vector space of .nite dimension, n ¡Ý 1. Say that two bases, (u1,...,un) and (v1,...,vn), of E have the same orientation i. det(P ) > 0, where P the change of basis matrix from (u1,...,un) and (v1,...,vn), namely, the matrix whose jth columns consist of the coordinates of vj over the basis (u1,...,un). 
(1) Prove that having the same orientation is an equivalence relation with two equivalence 
classes. An orientation of a vector space, E, is the choice of any .xed basis, say (e1,...,en), of 
E. Any other basis, (v1,...,vn), has the same orientation as (e1,...,en) (and is said to be positive or direct) i. det(P ) > 0, else it is said to have the opposite orientation of (e1,...,en) (or to be negative or indirect), where P is the change of basis matrix from (e1,...,en) to (v1,...,vn). An oriented vector space is a vector space with some chosen orientation (a positive basis). 
(2) Let B1 =(u1,...,un) and B2 =(v1,...,vn) be two orthonormal bases. For any sequence of vectors, (w1,...,wn), in E, let detB1 (w1,...,wn) be the determinant of the matrix whose columns are the coordinates of the wj¡¯s over the basis B1 and similarly for detB2 (w1,...,wn). 
Prove that if B1 and B2 have the same orientation, then 
detB1 (w1,...,wn) = detB2 (w1,...,wn). 
Given any oriented vector space, E, for any sequence of vectors, (w1,...,wn), in E, the common value, detB(w1,...,wn), for all positive orthonormal bases, B, of E is denoted 
¦ËE(w1,...,wn) 
and called a volume form of (w1,...,wn). 
(3) Given any Euclidean oriented vector space, E, of dimension n for any n . 1 vectors, w1,...,wn.1, in E, check that the map 
x ¡ú ¦ËE(w1,...,wn.1,x) 
is a linear form. Then prove that there is a unique vector, denoted w1 ¡Á¡¤ ¡¤¡¤¡Á wn.1, such that 
¦ËE(w1,...,wn.1,x)=(w1 ¡Á¡¤ ¡¤¡¤¡Á wn.1) ¡¤ x, 
for all x ¡Ê E. The vector w1 ¡Á¡¤ ¡¤¡¤¡Á wn.1 is called the cross-product of (w1,...,wn.1). Itis a generalization of the cross-product in R3 (when n = 3). 
Problem 11.9. Given p vectors (u1,...,up) in a Euclidean space E of dimension n ¡Ý p, the Gram determinant (or Gramian) of the vectors (u1,...,up) is the determinant 
lu1l2 Ju1,u2¡¢ ... Ju1,up¡¢ Ju2,u1¡¢lu2l2 ... Ju2,up¡¢ 
Gram(u1,...,up)= .. . . 
.
. ...
.
.. . 
LLLLLLLLL 
LLLLLLLL2LJ¡¢J¡¢llu,uu,uu... 12ppp
(1) Prove that 
Gram(u1,...,un)= ¦ËE(u1,...,un)2 . 

Hint. If(e1,...,en) is an orthonormal basis and A is the matrix of the vectors (u1,...,un) over this basis, 
det(A)2 = det(ATA) = det(Ai ¡¤ Aj), 
where Ai denotes the ith column of the matrix A, and (Ai ¡¤ Aj) denotes the n ¡Á n matrix with entries Ai ¡¤ Aj. 
(2) Prove that 
lu1 ¡Á¡¤ ¡¤¡¤¡Á un.1l2 = Gram(u1,...,un.1). 

Hint. Letting w = u1 ¡Á¡¤ ¡¤¡¤¡Á un.1, observe that 
¦ËE(u1,...,un.1,w)= Jw, w¡¢ = lwl2 , 
and show that 
lwl4 = ¦ËE(u1,...,un.1,w)2 = Gram(u1,...,un.1,w) = Gram(u1,...,un.1)lwl2 . 
Problem 11.10. Let . : E ¡Á E ¡ú R be a bilinear form on a real vector space E of .nite dimension n. Given any basis (e1,...,en) of E, let A =(aij) be the matrix de.ned such that 
aij = .(ei,ej), 
1 ¡Ü i, j ¡Ü n. We call A the matrix of . w.r.t. the basis (e1,...,en). 
(1) For any two vectors x and y, if X and Y denote the column vectors of coordinates of x and y w.r.t. the basis (e1,...,en), prove that 

.(x, y)= XTAY. 

(2) Recall that A is a symmetric matrix if A = AT . Prove that . is symmetric if A is a symmetric matrix. 

(3) If (f1,...,fn) is another basis of E and P is the change of basis matrix from (e1,...,en) to (f1,...,fn), prove that the matrix of . w.r.t. the basis (f1,...,fn) is 


P TAP. 
The common rank of all matrices representing . is called the rank of .. 
11.11. PROBLEMS 
Problem 11.11. Let .: E ¡Á E ¡ú R be a symmetric bilinear form on a real vector space E of .nite dimension n. Two vectors x and y are said to be conjugate or orthogonal w.r.t. . if .(x, y) = 0. The main purpose of this problem is to prove that there is a basis of vectors that are pairwise conjugate w.r.t. .. 
(1) Prove that if .(x, x)=0 for all x ¡Ê E, then . is identically null on E. 
Otherwise, we can assume that there is some vector x ¡Ê E such that .(x, x) = 0. 
Use induction to prove that there is a basis of vectors (u1,...,un) that are pairwise 

conjugate w.r.t. .. Hint. For the induction step, proceed as follows. Let (u1,e2,...,en) be a basis of E, with .(u1,u1) = 0. Prove that there are scalars ¦Ë2,...,¦Ën such that each of the vectors 
vi = ei + ¦Ëiu1 
is conjugate to u1 w.r.t. ., where 2 ¡Ü i ¡Ü n, and that (u1,v2,...,vn) is a basis. 
(2) Let (e1,...,en) be a basis of vectors that are pairwise conjugate w.r.t. . and assume that they are ordered such that 
¦Èi =0 if1 ¡Ü i ¡Ü r,
.(ei,ei)= 
0 if r +1 ¡Ü i ¡Ü n, 
where r is the rank of .. Show that the matrix of . w.r.t. (e1,...,en) is a diagonal matrix, and that 
r
n 
.(x, y)= ¦Èixiyi, i=1 
aa 
where x = in =1 xiei and y = in =1 yiei. Prove that for every symmetric matrix A, there is an invertible matrix P such that 
P TAP = D, 
where D is a diagonal matrix. 
(3) Prove that there is an integer p,0 ¡Ü p ¡Ü r (where r is the rank of .), such that .(ui,ui) > 0 for exactly p vectors of every basis (u1,...,un) of vectors that are pairwise conjugate w.r.t. . (Sylvester¡¯s inertia theorem). 
Proceed as follows. Assume that in the basis (u1,...,un), for any x ¡Ê E, we have 
222 2
.(x, x)= ¦Á1x1 + ¡¤¡¤¡¤ + ¦Ápxp . ¦Áp+1xp+1 .¡¤ ¡¤¡¤. ¦Árxr, 
a 
where x = ni=1 xiui, and that in the basis (v1,...,vn), for any x ¡Ê E, we have 
222 2
.(x, x)= ¦Â1y1 + ¡¤¡¤¡¤ + ¦Âqyq . ¦Âq+1yq+1 .¡¤ ¡¤¡¤. ¦Âryr , 
a 
where x = in =1 yivi, with ¦Ái > 0, ¦Âi > 0, 1 ¡Ü i ¡Ü r. 
Assume that p>q and derive a contradiction. First consider x in the subspace F spanned by 
(u1,...,up,ur+1,...,un), 
and observe that .(x, x) ¡Ý 0 if x = 0. Next consider x in the subspace G spanned by 
(vq+1,...,vr), 
and observe that .(x, x) < 0 if x = 0. Prove that F ¡É G is nontrivial (i.e., contains some nonnull vector), and derive a contradiction. This implies that p ¡Ü q. Finish the proof. The pair (p, r . p) is called the signature of .. 
(4) A symmetric bilinear form . is de.nite if for every x ¡Ê E, if .(x, x) = 0, then x = 0. Prove that a symmetric bilinear form is de.nite i. its signature is either (n, 0) or (0,n). In 
other words, a symmetric de.nite bilinear form has rank n and is either positive or negative. Problem 11.12. Consider the n ¡Á n matrices Ri,j de.ned for all i, j with 1 ¡Ü i<j ¡Ü n and n ¡Ý 3, such that the only nonzero entries are 
Ri,j(i, j)= .1 Ri,j(i, i)=0 Ri,j(j, i)=1 Ri,j(j, j)=0 
Ri,j(k, k)=1, 1 ¡Ü k ¡Ü n, k = i, j. 
For example, 

.
. 
1 ... 
1 00 ¡¤¡¤¡¤ 0 .1 01 ¡¤¡¤¡¤ 00 
................... 

................... 

Ri,j 
= 

.. ..
.
.. ... 
.

.

.. 

.. 

0  0  ¡¤ ¡¤ ¡¤  1  0  
1  0  ¡¤ ¡¤ ¡¤  0  0  
1  
...  
1  

(1) 
Prove that the Ri,j are rotation matrices. Use the matrices Rij to form a basis of the n ¡Á n skew-symmetric matrices. 

(2) 
Consider the n ¡Á n symmetric matrices Si,j de.ned for all i, j with 1 ¡Ü i<j ¡Ü n and n ¡Ý 3, such that the only nonzero entries are Si,j(i, j)=1 Si,j(i, i)=0 Si,j(j, i)=1 


Si,j(j, j)=0 Si,j(k, k)=1, 1 ¡Ü k ¡Ü n, k = i, j, 
11.11. PROBLEMS 
and if i +2 ¡Ü j then Si,j(i +1,i +1) = .1, else if i> 1 and j = i + 1 then Si,j(1, 1) = .1, and if i = 1 and j = 2, then Si,j(3, 3) = .1. 
For example, 
.
. 
1 .. . 
1 00 ¡¤¡¤¡¤ 01 0 .1 ¡¤¡¤¡¤ 00 
Si,j 
= 

................... 

................... 

. 

.. ..
.
.. 

...

.

.. 

.. 

0  0  ¡¤ ¡¤ ¡¤  1  0  
1  0  ¡¤ ¡¤ ¡¤  0  0  
1  
...  
1  

Note that Si,j has a single diagonal entry equal to .1. Prove that the Si,j are rotations matrices. 
Use Problem 2.15 together with the Si,j to form a basis of the n ¡Á n symmetric matrices. 
(3) Prove that if n ¡Ý 3, the set of all linear combinations of matrices in SO(n) is the space Mn(R) of all n ¡Á n matrices. Prove that if n ¡Ý 3 and if a matrix A ¡Ê Mn(R) commutes with all rotations matrices, 
then A commutes with all matrices in Mn(R). What happens for n = 2? 
Problem 11.13. Let A be an n ¡Á n real invertible matrix. Prove that if A = Q1R1 and A = Q2R2 are two QR-decompositions of A where R1 and R2 are upper-triangular with positive diagonal entries, then Q1 = Q2 and R1 = R2. 
Problem 11.14. (1) Let H be the a.ne hyperplane in Rn given by the equation 
a1x1 + ¡¤¡¤¡¤ + anxn = c, 
with ai = 0 for some i, 1 ¡Ü i ¡Ü n. The linear hyperplane H0 parallel to H is given by the equation 
a1x1 + ¡¤¡¤¡¤ + anxn =0, 
and we say that a vector y ¡Ê Rn is orthogonal (or perpendicular) to H i. y is orthogonal to H0. Let h be the intersection of H with the line through the origin and perpendicular to H. Prove that the coordinates of h are given by 
c 
22 (a1,...,an). 
a1 + ¡¤¡¤¡¤ + an 
(2) For any point p ¡Ê H, prove that lhl¡Ülpl. Thus, it is natural to de.ne the distance d(O, H) from the origin O to the hyperplane H as d(O, H)= lhl. Prove that 
|c|
d(O, H)= . (a12 + ¡¤¡¤¡¤ + a2 n)21 
(3) Let S be a .nite set of n ¡Ý 3 points in the plane (R2). Prove that if for every pair of distinct points pi,pj ¡Ê S, there is a third point pk ¡Ê S (distinct from pi and pj) such that pi,pj,pk belong to the same (a.ne) line, then all points in S belong to a common (a.ne) line. Hint. Proceed by contradiction and use a minimality argument. This is either ¡Þ-hard or relatively easy, depending how you proceed! 
Problem 11.15. (The space of closed polygons in R2, after Hausmann and Knutson) 
An open polygon P in the plane is a sequence P =(v1,...,vn+1) of points vi ¡Ê R2 called vertices (with n ¡Ý 1). A closed polygon, for short a polygon, is an open polygon P =(v1,...,vn+1) such that vn+1 = v1. The sequence of edge vectors (e1,...,en) associated with the open (or closed) polygon P =(v1,...,vn+1) is de.ned by 
ei = vi+1 . vi,i =1, . . . , n. 
Thus, a closed or open polygon is also de.ned by a pair (v1, (e1,...,en)), with the vertices given by 
vi+1 = vi + ei,i =1, . . . , n. 
Observe that a polygon (v1, (e1,...,en)) is closed i. 
e1 + ¡¤¡¤¡¤ + en =0. 
Since every polygon (v1, (e1,...,en)) can be translated by .v1, so that v1 = (0, 0), we may assume that our polygons are speci.ed by a sequence of edge vectors. 
Recall that the plane R2 is isomorphic to C, via the isomorphism 
(x, y) ¡ú x + iy. 
We will represent each edge vector ek by the square of a complex number wk = ak +ibk. Thus, every sequence of complex numbers (w1,...,wn) de.nes a polygon (namely, (w12,...,wn2 )). This representation is many-to-one: the sequences (¡Àw1,..., ¡Àwn) describe the same poly-gon. To every sequence of complex numbers (w1,...,wn), we associate the pair of vectors (a, b), with a, b ¡Ê Rn, such that if wk = ak + ibk, then 
a =(a1,...,an),b =(b1,...,bn). 
The mapping (w1,...,wn) ¡ú (a, b) 
11.11. PROBLEMS 
is clearly a bijection, so we can also represent polygons by pairs of vectors (a, b) ¡Ê Rn ¡Á Rn . 
(1) 
Prove that a polygon P represented by a pair of vectors (a, b) ¡Ê Rn ¡Á Rn is closed i. a ¡¤ b = 0 and lal2 = lbl2. 

(2) 
Given a polygon P represented by a pair of vectors (a, b) ¡Ê Rn ¡Á Rn, the length l(P ) of the polygon P is de.ned by l(P )= |w1|2 + ¡¤¡¤¡¤ + |wn|2, with wk = ak + ibk. Prove that 


l(P )= lal22 + lbl22 . 
Deduce from (a) and (b) that every closed polygon of length 2 with n edges is represented by a n ¡Á 2 matrix A such that ATA = I. 
Remark: The space of all a n ¡Á 2 real matrices A such that ATA = I is a space known as the Stiefel manifold S(2,n). 
(3) Recall that in R2, the rotation of angle ¦È speci.ed by the matrix 
cos ¦È . sin ¦È 
R¦È = 
sin ¦È cos ¦È 
is expressed in terms of complex numbers by the map 
z ¡ú ze i¦È . 
Let P be a polygon represented by a pair of vectors (a, b) ¡Ê Rn ¡Á Rn . Prove that the polygon R¦È(P ) obtained by applying the rotation R¦È to every vertex wk 2 =(ak + ibk)2 of P is speci.ed by the pair of vectors 
(cos(¦È/2)a . sin(¦È/2)b, sin(¦È/2)a + cos(¦È/2)b)

.
. 

.... 

a1 b1 a2 b2 
.. 

.. 
.. 
an bn 
.... 

cos(¦È/2) sin(¦È/2) 

. sin(¦È/2) cos(¦È/2) 

.

= 

(4) The re.ection ¦Ñx about the x-axis corresponds to the map 
z ¡ú z, 

whose matrix is, 
10 
. 
0 .1 
Prove that the polygon ¦Ñx(P ) obtained by applying the re.ection ¦Ñx to every vertex wk 2 = (ak + ibk)2 of P is speci.ed by the pair of vectors 
.
. 

(a, .b)= 

.... 

a1 b1 a2 b2 
.. 

.. 
.. 
an bn 
.... 

10 

0 .1 

. 

(5) Let Q ¡Ê O(2) be any isometry such that det(Q)= .1 (a re.ection). Prove that there is a rotation R.¦È ¡Ê SO(2) such that 
Q = ¦Ñx . R.¦È. 
Prove that the isometry Q, which is given by the matrix 
cos ¦È sin ¦È 
Q = ,
sin ¦È . cos ¦È 
is the re.ection about the line corresponding to the angle ¦È/2 (the line of equation y = tan(¦È/2)x). 
Prove that the polygon Q(P ) obtained by applying the re.ection Q = ¦Ñx . R.¦È to every vertex wk 2 =(ak + ibk)2 of P , is speci.ed by the pair of vectors 
(cos(¦È/2)a + sin(¦È/2)b, sin(¦È/2)a . cos(¦È/2)b)

.
. 

.... 

a1 b1 a2 b2 
.. 

.. 
.. 
an bn 
.... 

cos(¦È/2) sin(¦È/2) 

sin(¦È/2) . cos(¦È/2) 

.

= 

(6) De.ne an equivalence relation ¡« on S(2,n) such that if A1,A2 ¡Ê S(2,n) are any n¡Á2 matrices such that AT = AT = I, then 
1 A12 A2 A1 ¡« A2 i. A2 = A1Q for some Q ¡Ê O(2). 
Prove that the quotient G(2,n)= S(2,n)/ ¡« is in bijection with the set of all 2-dimensional subspaces (the planes) of Rn . The space G(2,n) is called a Grassmannian manifold. 
Prove that up to translations and isometries in O(2) (rotations and re.ections), the n-sided closed polygons of length 2 are represented by planes in G(2,n). 
Problem 11.16. (1) Find two symmetric matrices, A and B, such that AB is not symmetric. 
(2) Find two matrices A and B such that 
AB A+B 
ee = e. 
Hint. Try 

.
.
.
. 
000 001 

.
0 01 0 .100 and use the Rodrigues formula. 
(3) Find some square matrices A, B such that AB = BA, yet 
AB A+B 
ee = e. 
.

.

.

A = ¦Ð 

00 .1 

and B = ¦Ð 

00 

, 

Hint. Look for 2 ¡Á 2 matrices with zero trace and use Problem 8.15. 

11.11. PROBLEMS 
Problem 11.17. Given a .eld K and any nonempty set I, let K(I) be the subset of the cartesian product KI consisting of all functions ¦Ë: I ¡ú K with .nite support, which means that ¦Ë(i) = 0 for all but .nitely many i ¡Ê I. We usually denote the function de.ned by ¦Ë as (¦Ëi)i¡ÊI , and call is a family indexed by I. We de.ne addition and multiplication by a scalar as follows: 
(¦Ëi)i¡ÊI +(¦Ìi)i¡ÊI =(¦Ëi + ¦Ìi)i¡ÊI , 
and ¦Á ¡¤ (¦Ìi)i¡ÊI =(¦Á¦Ìi)i¡ÊI . 
(1) 
Check that K(I) is a vector space. 

(2) 
If I is any nonempty subset, for any i ¡Ê I, we denote by ei the family (ej)j¡ÊI de.ned 


so that 1 if j = i 
ej = 
0 if j = i. 
Prove that the family (ei)i¡ÊI is linearly independent and spans K(I), so that it is a basis of K(I) called the canonical basis of K(I). When I is .nite, say of cardinality n, then prove that K(I) is isomorphic to Kn . 
(3) The function ¦É: I ¡ú K(I), such that ¦É(i)= ei for every i ¡Ê I, is clearly an injection. For any other vector space F , for any function f : I ¡ú F , prove that there is a unique 
linear map f : K(I) ¡ú F , such that f = f . ¦É, 
as in the following commutative diagram: 
¦É
I  . 
K(I) . 
        
f f
  j 
F 
We call the vector space K(I) the vector space freely generated by the set I. 
Problem 11.18. (Some pitfalls of in.nite dimension) Let E be the vector space freely generated by the set of natural numbers, N = {0, 1, 2,...}, and let (e0,e1,e2,...,en,...) be its canonical basis. We de.ne the function . such that 
. 
¦Äij if i, j ¡Ý 1,
.

.
1 if i = j =0,
.(ei,ej)= 
.
1/2j if i =0,j ¡Ý 1, . 1/2i if i ¡Ý 1,j =0, 
a 
and we extend . by bilinearity to a function .: E¡ÁE ¡ú K. This means that if u = i¡ÊN ¦Ëiei
a 
and v = j¡ÊN ¦Ìjej, then 
nn n 
.¦Ëiei,¦Ìjej = ¦Ëi¦Ìj.(ei,ej), i¡ÊN j¡ÊN i,j¡ÊN 
but remember that ¦Ëi = 0 and ¦Ìj =0 only for .nitely many indices i, j. 
(1) Prove that . is positive de.nite, so that it is an inner product on E. 
What would happen if we changed 1/2j to 1 (or any constant)? 

(2) Let H be the subspace of E spanned by the family (ei)i¡Ý1, a hyperplane in E. Find H¡Í and H¡Í¡Í, and prove that 
= H¡Í¡Í
H. 
(3) Let U be the subspace of E spanned by the family (e2i)i¡Ý1, and let V be the subspace of E spanned by the family (e2i.1)i¡Ý1. Prove that 
U¡Í 
= V V ¡Í 
= U 
U¡Í¡Í 
= U 
V ¡Í¡Í 
= V, 
yet 
(U ¡É V )¡Í = U¡Í + V ¡Í and 
(U + V )¡Í¡Í 
= U + V. If W is the subspace spanned by e0 and e1, prove that 
(W ¡É H)¡Í = W ¡Í + H¡Í . 
(4) Consider the dual space E. of E, and let (ei .)i¡ÊN be the family of dual forms of the basis (ei)i¡ÊN . Check that the family (ei .)i¡ÊN is linearly independent. 
(5) Let f ¡Ê E. be the linear form de.ned by 
f(ei)=1 forall i ¡Ê N. 
Prove that f is not in the subspace spanned by the ei . . If F is the subspace of E. spanned by the ei . and f, .nd F 0 and F 00, and prove that 
= F 00
F. 



