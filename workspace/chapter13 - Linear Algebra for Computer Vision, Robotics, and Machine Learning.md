Chapter 13 
Hermitian Spaces 
13.1 	Sesquilinear and Hermitian Forms, Pre-Hilbert Spaces and Hermitian Spaces 
In this chapter we generalize the basic results of Euclidean geometry presented in Chapter 11 to vector spaces over the complex numbers. Such a generalization is inevitable and not simply a luxury. For example, linear maps may not have real eigenvalues, but they always have complex eigenvalues. Furthermore, some very important classes of linear maps can be diagonalized if they are extended to the complexi.cation of a real vector space. This is the case for orthogonal matrices and, more generally, normal matrices. Also, complex vector spaces are often the natural framework in physics or engineering, and they are more convenient for dealing with Fourier series. However, some complications arise due to complex conjugation. 
Recall that for any complex number z ¡Ê C, if z = x + iy where x, y ¡Ê R, we let Rz = x, 
the real part of z, and Sz = y, the imaginary part of z. We also denote the conjugate of z = x + iy by z = x . iy, and the absolute value (or length, or modulus) of z by |z|. Recall that |z|2 = zz = x2 + y2 . 
There are many natural situations where a map .: E ¡Á E ¡ú C is linear in its .rst argument and only semilinear in its second argument, which means that .(u, ¦Ìv)= ¦Ì.(u, v), as opposed to .(u, ¦Ìv)= ¦Ì.(u, v). For example, the natural inner product to deal with functions f : R ¡ú C, especially Fourier series, is
 ¦Ð 
£¨f,gf(x)g(x)dx,
´¨ =
.¦Ð 
which is semilinear (but not linear) in g. Thus, when generalizing a result from the real case of a Euclidean space to the complex case, we always have to check very carefully that our proofs do not rely on linearity in the second argument. Otherwise, we need to revise our proofs, and sometimes the result is simply wrong! 
Before de.ning the natural generalization of an inner product, it is convenient to de.ne semilinear maps. 
449 

De.nition 13.1. Given two vector spaces E and F over the complex .eld C, a function 
f : E ¡ú F is semilinear if 
f(u + v)= f(u)+ f(v), f(¦Ëu)= ¦Ëf(u), 
for all u, v ¡Ê E and all ¦Ë ¡Ê C. 
Remark: Instead of de.ning semilinear maps, we could have de.ned the vector space E as the vector space with the same carrier set E whose addition is the same as that of E, but whose multiplication by a complex number is given by 
(¦Ë, u) ¡ú ¦Ëu. 
Then it is easy to check that a function f : E ¡ú C is semilinear i. f : E ¡ú C is linear. We can now de.ne sesquilinear forms and Hermitian forms. 
De.nition 13.2. Given a complex vector space E, a function .: E¡ÁE ¡ú Cis a sesquilinear form if it is linear in its .rst argument and semilinear in its second argument, which means that 
.(u1 + u2,v)= .(u1,v)+ .(u2,v), 
.(u, v1 + v2)= .(u, v1)+ .(u, v2), 
.(¦Ëu, v)= ¦Ë.(u, v), 
.(u, ¦Ìv)= ¦Ì.(u, v), 
for all u, v, u1,u2, v1,v2 ¡Ê E, and all ¦Ë, ¦Ì ¡Ê C. A function .: E ¡Á E ¡ú C is a Hermitian form if it is sesquilinear and if 
.(v, u)= .(u, v) 
for all all u, v ¡Ê E. 
Obviously, .(0,v)= .(u, 0) = 0. Also note that if .: E ¡Á E ¡ú C is sesquilinear, we have .(¦Ëu + ¦Ìv, ¦Ëu + ¦Ìv)= |¦Ë|2.(u, u)+ ¦Ë¦Ì.(u, v)+ ¦Ë¦Ì.(v, u)+ |¦Ì|2.(v, v), 
and if .: E ¡Á E ¡ú C is Hermitian, we have 
.(¦Ëu + ¦Ìv, ¦Ëu + ¦Ìv)= |¦Ë|2.(u, u)+2R(¦Ë¦Ì.(u, v)) + |¦Ì|2.(v, v). 
Note that restricted to real coe.cients, a sesquilinear form is bilinear (we sometimes say R-bilinear). 
De.nition 13.3. Given a sesquilinear form . : E ¡Á E ¡ú C, the function ¦µ: E ¡ú C de.ned such that ¦µ(u)= .(u, u) for all u ¡Ê E is called the quadratic form associated with .. 
13.1. HERMITIAN SPACES, PRE-HILBERT SPACES 
The standard example of a Hermitian form on Cn is the map . de.ned such that .((x1,...,xn), (y1,...,yn)) = x1y1 + x2y2 + ¡¤¡¤¡¤ + xnyn. 
This map is also positive de.nite, but before dealing with these issues, we show the following useful proposition. 
Proposition 13.1. Given a complex vector space E, the following properties hold: 
(1) 
A sesquilinear form .: E ¡Á E ¡ú C is a Hermitian form i. .(u, u) ¡Ê R for all u ¡Ê E. 

(2) 
If .: E ¡Á E ¡ú C is a sesquilinear form, then 4.(u, v)= .(u + v, u + v) . .(u . v, u . v) 


+ i.(u + iv, u + iv) . i.(u . iv, u . iv), 
and 
2.(u, v) = (1+ i)(.(u, u)+ .(v, v)) . .(u . v, u . v) . i.(u . iv, u . iv). 
These are called polarization identities. Proof. (1) If . is a Hermitian form, then .(v, u)= .(u, v) implies that .(u, u)= .(u, u), and thus .(u, u) ¡Ê R. If . is sesquilinear and .(u, u) ¡Ê R for all u ¡Ê E, then .(u + v, u + v)= .(u, u)+ .(u, v)+ .(v, u)+ .(v, v), which proves that .(u, v)+ .(v, u)= ¦Á, where ¦Á is real, and changing u to iu, we have i(.(u, v) . .(v, u)) = ¦Â, 
where ¦Â is real, and thus  
¦Á . i¦Â  ¦Á + i¦Â  
.(u, v) =  2  and  .(v, u) =  2  ,  
proving that . is Hermitian.  

(2) These identities are veri.ed by expanding the right-hand side, and we leave them as an exercise. 
Proposition 13.1 shows that a sesquilinear form is completely determined by the quadratic form ¦µ(u)= .(u, u), even if . is not Hermitian. This is false for a real bilinear form, unless it is symmetric. For example, the bilinear form .: R2 ¡Á R2 ¡ú R de.ned such that 
.((x1,y1), (x2,y2)) = x1y2 . x2y1 
is not identically zero, and yet it is null on the diagonal. However, a real symmetric bilinear form is indeed determined by its values on the diagonal, as we saw in Chapter 11. 
As in the Euclidean case, Hermitian forms for which .(u, u) ¡Ý 0 play an important role. 
De.nition 13.4. Given a complex vector space E, a Hermitian form .: E ¡Á E ¡ú C is positive if .(u, u) ¡Ý 0 for all u ¡Ê E, and positive de.nite if .(u, u) > 0 for all u 0.
=A pair£¨E, .´¨ where E is a complex vector space and . is a Hermitian form on E is called a pre-Hilbert space if . is positive, and a Hermitian (or unitary) space if . is positive de.nite. 
We warn our readers that some authors, such as Lang [42], de.ne a pre-Hilbert space as what we de.ne as a Hermitian space. We prefer following the terminology used in Schwartz 
[54] and Bourbaki [10]. The quantity .(u, v) is usually called the Hermitian product of u 
and v. We will occasionally call it the inner product of u and v. Given a pre-Hilbert space£¨E, .´¨, as in the case of a Euclidean space, we also denote 
.(u, v) by u ¡¤ v or£¨u, v´¨ or (u|v), 
 
and¦µ(u) by lul. 
Example 13.1. The complex vector space Cn under the Hermitian form 
.((x1,...,xn), (y1,...,yn)) = x1y1 + x2y2 + ¡¤¡¤¡¤ + xnyn 
is a Hermitian space. 
Example 13.2. LetÙ²2 denote the set of all countably in.nite sequences x =(xi)i¡ÊN of
o¡Þ o n
complex numbers such that |xi|2 is de.ned (i.e., the sequence |xi|2 converges as 
i=0 i=0 
n ¡ú¡Þ). It can be shown that the map .:Ù²2 ¡ÁÙ²2 ¡ú C de.ned such that 
¡Þ
C 
. ((xi)i¡ÊN, (yi)i¡ÊN)= xiyi 
i=0 
is well de.ned, andÙ²2 is a Hermitian space under .. Actually,Ù²2 is even a Hilbert space. 
Example 13.3. Let Cpiece[a, b] be the set of bounded piecewise continuous functions 
f :[a, b] ¡ú C under the Hermitian form
b 
£¨f, g´¨ = f(x)g(x)dx. 
a 
It is easy to check that this Hermitian form is positive, but it is not de.nite. Thus, under this Hermitian form, Cpiece[a, b] is only a pre-Hilbert space. 

13.1. HERMITIAN SPACES, PRE-HILBERT SPACES 
Example 13.4. Let C[a, b] be the set of complex-valued continuous functions f :[a, b] ¡ú C under the Hermitian form
b 
£¨f, g´¨ = f(x)g(x)dx. 
a 
It is easy to check that this Hermitian form is positive de.nite. Thus, C[a, b] is a Hermitian space. 
Example 13.5. Let E =Mn(C) be the vector space of complex n ¡Á n matrices. If we view a matrix A ¡Ê Mn(C) as a ¡°long¡± column vector obtained by concatenating together its columns, we can de.ne the Hermitian product of two matrices A, B ¡Ê Mn(C) as
n
C 
£¨A, B´¨ = aijbij, i,j=1 
which can be conveniently written as
£¨A, B´¨ = tr(ATB) = tr(B . A). 
Since this can be viewed as the standard Hermitian product on Cn2 , it is a Hermitian product on Mn(C). The corresponding norm 
lAlF = tr(A.A) 
is the Frobenius norm (see Section 8.2). 
If E is .nite-dimensional and if .: E ¡Á E ¡ú R is a sequilinear form on E, given any 
oo 
basis (e1,...,en) of E, we can write x = n xiei and y = n yjej, and we have 
i=1 j=1 
Cnnn
CC 
.(x, y)= .xiei,yjej = xiyj.(ei,ej). i=1 j=1 i,j=1 
If we let G =(gij) be the matrix given by gij = .(ej,ei), and if x and y are the column vectors associated with (x1,...,xn) and (y1,...,yn), then we can write 
TGT
.(x, y)= xy = y . Gx, 
where y corresponds to (y1,..., y). As in Section 11.1, we are committing the slight abuse of 
no 
notation of letting x denote both the vector x = ni=1 xiei and the column vector associated with (x1,...,xn) (and similarly for y). The ¡°correct¡± expression for .(x, y) is 
.(x, y)= y . Gx = xTGTy.
 Observe that in .(x, y)= y .Gx, the matrix involved is the transpose of the matrix (.(ei,ej)). The reason for this is that we want G to be positive de.nite when . is positive de.nite, not GT. 
Furthermore, observe that . is Hermitian i. G = G., and . is positive de.nite i. the matrix G is positive de.nite, that is, 
(Gx)T x = x . Gx > 0 for all x ¡Ê Cn,x =0. 
De.nition 13.5. The matrix G associated with a Hermitian product is called the Gram matrix of the Hermitian product with respect to the basis (e1,...,en). 
Conversely, if A is a Hermitian positive de.nite n ¡Á n matrix, it is easy to check that the Hermitian form£¨x, y´¨ = y . Ax 
is positive de.nite. If we make a change of basis from the basis (e1,...,en) to the basis (f1,...,fn), and if the change of basis matrix is P (where the jth column of P consists of 
'
the coordinates of fj over the basis (e1,...,en)), then with respect to coordinates x' and yover the basis (f1,...,fn), we have 
y . Gx =(y') . P . GP x', 
so the matrix of our inner product over the basis (f1,...,fn) is P .GP . We summarize these facts in the following proposition. 
Proposition 13.2. Let E be a .nite-dimensional vector space, and let (e1,...,en) be a basis of E. 
1. For any Hermitian inner product£¨.,.´¨ on E, if G =(gij) with gij =£¨ej,ei´¨ is the Gram matrix of the Hermitian product£¨.,.´¨ w.r.t. the basis (e1,...,en), then G is Hermitian positive de.nite. 

2. For any change of basis matrix P , the Gram matrix of£¨.,.´¨ with respect to the new basis is P .GP . 


3. If A is any n ¡Á n Hermitian positive de.nite matrix, then
£¨x, y´¨ = y . Ax 
is a Hermitian product on E. 
We will see later that a Hermitian matrix is positive de.nite i. its eigenvalues are all positive. 
The following result reminiscent of the .rst polarization identity of Proposition 13.1 can be used to prove that two linear maps are identical. 
Proposition 13.3. Given any Hermitian space E with Hermitian product£¨.,.´¨, for any linear map f : E ¡ú E, if£¨f(x),x´¨ =0 for all x ¡Ê E, then f =0. 

13.1. HERMITIAN SPACES, PRE-HILBERT SPACES 
Proof. Compute£¨f(x + y),x + y´¨ and£¨f(x . y),x . y´¨:£¨f(x + y),x + y´¨ =£¨f(x),x´¨ +£¨f(x),y´¨ +£¨f(y),x´¨ +£¨y, y´¨£¨f(x . y),x . y´¨ =£¨f(x),x´¨.£¨ f(x),y´¨.£¨ f(y),x´¨ +£¨y, y´¨; then subtract the second equation from the .rst to obtain£¨f(x + y),x + y´¨.£¨ f(x . y),x . y´¨ = 2(£¨f(x),y´¨ +£¨f(y),x´¨). If£¨f(u),u´¨ = 0 for all u ¡Ê E, we get£¨f(x),y´¨ +£¨f(y),x´¨ = 0 for all x, y ¡Ê E. Then the above equation also holds if we replace x by ix, and we obtain i£¨f(x),y´¨. i£¨f(y),x´¨ =0, for all x, y ¡Ê E, so we have
£¨f(x),y´¨ +£¨f(y),x´¨ =0
£¨f(x),y´¨.£¨ f(y),x´¨ =0, which implies that£¨f(x),y´¨ = 0 for all x, y ¡Ê E. Since£¨.,.´¨ is positive de.nite, we have 
f(x)=0 for all x ¡Ê E; that is, f = 0. 
One should be careful not to apply Proposition 13.3 to a linear map on a real Euclidean space because it is false! The reader should .nd a counterexample. 
The Cauchy¨CSchwarz inequality and the Minkowski inequalities extend to pre-Hilbert spaces and to Hermitian spaces. 
Proposition 13.4. Let£¨E, .´¨ be a pre-Hilbert space with associated quadratic form ¦µ. For all u, v ¡Ê E, we have the Cauchy¨CSchwarz inequality 
|.(u, v)|¡Ü ¦µ(u) ¦µ(v). 
Furthermore, if£¨E, .´¨ is a Hermitian space, the equality holds i. u and v are linearly de-pendent. 
We also have the Minkowski inequality 
¦µ(u + v) ¡Ü ¦µ(u)+ ¦µ(v). 
Furthermore, if£¨E, .´¨ is a Hermitian space, the equality holds i. u and v are linearly de-pendent, where in addition, if u =0 and v =0, then u = ¦Ëv for some real ¦Ë such that ¦Ë> 0. 
Proof. For all u, v ¡Ê E and all ¦Ì ¡Ê C, we have observed that 
.(u + ¦Ìv, u + ¦Ìv)= .(u, u)+2R(¦Ì.(u, v)) + |¦Ì|2.(v, v). 

Let .(u, v)= ¦Ñei¦È, where |.(u, v)| = ¦Ñ (¦Ñ ¡Ý 0). Let F : R ¡ú R be the function de.ned such that F (t) = ¦µ(u + tei¦È v), for all t ¡Ê R. The above shows that 
F (t)= .(u, u)+2t|.(u, v)| + t2.(v, v) = ¦µ(u)+2t|.(u, v)| + t2¦µ(v). Since . is assumed to be positive, we have F (t) ¡Ý 0 for all t ¡Ê R. If ¦µ(v) = 0, we must have .(u, v) = 0, since otherwise, F (t) could be made negative by choosing t negative and small enough. If ¦µ(v) > 0, in order for F (t) to be nonnegative, the equation 
¦µ(u)+2t|.(u, v)| + t2¦µ(v)=0 must not have distinct real roots, which is equivalent to |.(u, v)|2 ¡Ü ¦µ(u)¦µ(v). Taking the square root on both sides yields the Cauchy¨CSchwarz inequality. For the second part of the claim, if . is positive de.nite, we argue as follows. If u and v are linearly dependent, it is immediately veri.ed that we get an equality. Conversely, if |.(u, v)|2 = ¦µ(u)¦µ(v), then there are two cases. If ¦µ(v) = 0, since . is positive de.nite, we must have v = 0, so u and v are linearly dependent. Otherwise, the equation ¦µ(u)+2t|.(u, v)| + t2¦µ(v)=0 has a double root t0, and thus ¦µ(u + t0e i¦È v)=0. Since . is positive de.nite, we must have 
u + t0e i¦È v =0, 
which shows that u and v are linearly dependent. If we square the Minkowski inequality, we get ¦µ(u + v) ¡Ü ¦µ(u) + ¦µ(v)+2 ¦µ(u) ¦µ(v). However, we observed earlier that ¦µ(u + v) = ¦µ(u) + ¦µ(v)+2R(.(u, v)). 

13.1. HERMITIAN SPACES, PRE-HILBERT SPACES 
Thus, it is enough to prove that R(.(u, v)) ¡Ü ¦µ(u) ¦µ(v), but this follows from the Cauchy¨CSchwarz inequality |.(u, v)|¡Ü ¦µ(u) ¦µ(v) and the fact that Rz ¡Ü|z|. If . is positive de.nite and u and v are linearly dependent, it is immediately veri.ed that we get an equality. Conversely, if equality holds in the Minkowski inequality, we must have 
R(.(u, v)) =  ¦µ(u)  ¦µ(v),  
which implies that  
 
 
|.(u, v)| =  ¦µ(u)  ¦µ(v),  

since otherwise, by the Cauchy¨CSchwarz inequality, we would have R(.(u, v)) ¡Ü|.(u, v)| < ¦µ(u) ¦µ(v). Thus, equality holds in the Cauchy¨CSchwarz inequality, and R(.(u, v)) = |.(u, v)|. But then we proved in the Cauchy¨CSchwarz case that u and v are linearly dependent. Since we also just proved that .(u, v) is real and nonnegative, the coe.cient of proportionality between u and v is indeed nonnegative. As in the Euclidean case, if£¨E, .´¨ is a Hermitian space, the Minkowski inequality 
¦µ(u + v) ¡Ü ¦µ(u)+ ¦µ(v) 
shows that the map u ¡ú ¦µ(u) is a norm on E. The norm induced by . is called the Hermitian norm induced by .. We usually denote ¦µ(u) by lul, and the Cauchy¨CSchwarz inequality is written as 
|u ¡¤ v|¡Ülullvl. 
Since a Hermitian space is a normed vector space, it is a topological space under the topology induced by the norm (a basis for this topology is given by the open balls B0(u, ¦Ñ) of center u and radius ¦Ñ> 0, where 
B0(u, ¦Ñ)= {v ¡Ê E |lv . ul <¦Ñ}. 
If E has .nite dimension, every linear map is continuous; see Chapter 8 (or Lang [42, 43], Dixmier [18], or Schwartz [54, 55]). The Cauchy¨CSchwarz inequality 
|u ¡¤ v|¡Ülullvl 
shows that . : E ¡Á E ¡ú C is continuous, and thus, that ll is continuous. If£¨E, .´¨ is only pre-Hilbertian, lul is called a seminorm. In this case, the condition 
lul = 0 implies u =0 
is not necessarily true. However, the Cauchy¨CSchwarz inequality shows that if lul = 0, then u ¡¤ v = 0 for all v ¡Ê E. 
Remark: As in the case of real vector spaces, a norm on a complex vector space is induced by some positive de.nite Hermitian product£¨.,.´¨ i. it satis.es the parallelogram law: 
lu + vl2 + lu . vl2 = 2(lul2 + lvl2). 
This time the Hermitian product is recovered using the polarization identity from Proposition 
13.1: 4£¨u, v´¨ = lu + vl2 .lu . vl2 + i lu + ivl2 . i lu . ivl2 . 
It is easy to check that£¨u, u´¨ = lul2, and
£¨v, u´¨ =£¨u, v´¨
£¨iu, v´¨ = i£¨u, v´¨, 
so it is enough to check linearity in the variable u, and only for real scalars. This is easily done by applying the proof from Section 11.1 to the real and imaginary part of£¨u, v´¨; the details are left as an exercise. 
We will now basically mirror the presentation of Euclidean geometry given in Chapter 11 rather quickly, leaving out most proofs, except when they need to be seriously amended. 


13.2 Orthogonality, Duality, Adjoint of a Linear Map 
In this section we assume that we are dealing with Hermitian spaces. We denote the Her-mitian inner product by u ¡¤ v or£¨u, v´¨. The concepts of orthogonality, orthogonal family of vectors, orthonormal family of vectors, and orthogonal complement of a set of vectors are unchanged from the Euclidean case (De.nition 11.2). 
For example, the set C[.¦Ð, ¦Ð] of continuous functions f :[.¦Ð, ¦Ð] ¡ú C is a Hermitian space under the product
¦Ð 
£¨f, g´¨ = f(x)g(x)dx, 
.¦Ð 
and the family (eikx)k¡ÊZ is orthogonal. Propositions 11.4 and 11.5 hold without any changes. It is easy to show that
½Ð½Ð½Ð½Ð
½Ð½Ð½ÐC½Ð½Ð½Ð 
n
ui
C
C
2 
n
= 

luil2 +2R(ui ¡¤ uj). 
i=1 i=1 1¡Üi<j¡Ün 
13.2. ORTHOGONALITY, DUALITY, ADJOINT OF A LINEAR MAP 
Analogously to the case of Euclidean spaces of .nite dimension, the Hermitian product induces a canonical bijection (i.e., independent of the choice of bases) between the vector space E and the space E. . This is one of the places where conjugation shows up, but in this case, troubles are minor. 
Given a Hermitian space E, for any vector u ¡Ê E, let .lu : E ¡ú C be the map de.ned such that 
.l 
u(v)= u ¡¤ v, for all v ¡Ê E. Similarly, for any vector v ¡Ê E, let .rv : E ¡ú C be the map de.ned such that 
.rv(u)= u ¡¤ v, for all u ¡Ê E. 
Since the Hermitian product is linear in its .rst argument u, the map .rv is a linear form in E., and since it is semilinear in its second argument v, the map .lu is also a linear form in E. . Thus, we have two maps bl : E ¡ú E. and br : E ¡ú E., de.ned such that 
bl(u)= .lu, and br(v)= .rv. 
Proposition 13.5. The equations .lu = .ru and bl = br hold. Proof. Indeed, for all u, v ¡Ê E, we have 
bl(u)(v)= .lu(v) = u ¡¤ v = v ¡¤ u = .ru(v) = br(u)(v). 
Therefore, we use the notation .u for both .ul and .ur , and b for both bl and br . Theorem 13.6. Let E be a Hermitian space E. The map b: E ¡ú E. de.ned such that 
b(u)= .lu = .ur for all u ¡Ê E 
is semilinear and injective. When E is also of .nite dimension, the map b: E ¡ú E. is a canonical isomorphism. 
Proof. That b: E ¡ú E. is a semilinear map follows immediately from the fact that b = br , and that the Hermitian product is semilinear in its second argument. If .u = .v, then .u(w)= .v(w) for all w ¡Ê E, which by de.nition of .u and .v means that 
w ¡¤ u = w ¡¤ v for all w ¡Ê E, which by semilinearity on the right is equivalent to w ¡¤ (v . u)=0 forall w ¡Ê E, 
which implies that u = v, since the Hermitian product is positive de.nite. Thus, b: E ¡ú E. is injective. Finally, when E is of .nite dimension n, E. is also of dimension n, and then 
b: E ¡ú E. is bijective. Since b is semilinar, the map b: E ¡ú E. is an isomorphism. 
The inverse of the isomorphism b: E ¡ú E. is denoted by . : E. ¡ú E. 
As a corollary of the isomorphism b: E ¡ú E. we have the following result. 

Proposition 13.7. If E is a Hermitian space of .nite dimension, then every linear form f ¡Ê E. corresponds to a unique v ¡Ê E, such that 
f(u)= u ¡¤ v, for every u ¡Ê E. 
In particular, if f is not the zero form, the kernel of f, which is a hyperplane H, is precisely the set of vectors that are orthogonal to v. 
Remarks: 
1. The ¡°musical map¡± b: E ¡ú E. is not surjective when E has in.nite dimension. This result can be salvaged by restricting our attention to continuous linear maps and by assuming that the vector space E is a Hilbert space. 

2. 	Dirac¡¯s ¡°bra-ket¡± notation. Dirac invented a notation widely used in quantum me-chanics for denoting the linear form .u = b(u) associated to the vector u ¡Ê E via the duality induced by a Hermitian inner product. Dirac¡¯s proposal is to denote the vectors u in E by |u´¨, and call them kets; the notation |u´¨ is pronounced ¡°ket u.¡± Given two kets (vectors) |u´¨ and |v´¨, their inner product is denoted by


£¨u|v´¨ 
(instead of |u´¨¡¤|v´¨). The notation£¨u|v´¨ for the inner product of |u´¨ and |v´¨ anticipates duality. Indeed, we de.ne the dual (usually called adjoint) bra u of ket u, denoted by£¨u|, as the linear form whose value on any ket v is given by the inner product, so
£¨u|(|v´¨)=£¨u|v´¨. 
Thus, bra u =£¨u| is Dirac¡¯s notation for our b(u). Since the map b is semi-linear, we have£¨¦Ëu| = ¦Ë£¨u|. 
Using the bra-ket notation, given an orthonormal basis (|u1´¨,..., |un´¨), ket v (a vector) is written as 
n
C 
|v´¨ = £¨v|ui´¨|ui´¨, i=1
and the corresponding linear form bra v is written as
n	n
CC 
£¨v| = £¨v|ui´¨£¨ ui| = £¨ui|v´¨£¨ui|i=1i=1
over the dual basis (£¨u1|,...,£¨un|). As cute as it looks, we do not recommend using the Dirac notation. 

13.2. ORTHOGONALITY, DUALITY, ADJOINT OF A LINEAR MAP 
The existence of the isomorphism b: E ¡ú E. is crucial to the existence of adjoint maps. Indeed, Theorem 13.6 allows us to de.ne the adjoint of a linear map on a Hermitian space. Let E be a Hermitian space of .nite dimension n, and let f : E ¡ú E be a linear map. For every u ¡Ê E, the map 
v ¡ú u ¡¤ f(v) 
is clearly a linear form in E., and by Theorem 13.6, there is a unique vector in E denoted by f.(u), such that 
f.(u) ¡¤ v = u ¡¤ f(v), 
that is, f . (u) ¡¤ v = u ¡¤ f(v), for every v ¡Ê E. 
The following proposition shows that the map f. is linear. 
Proposition 13.8. Given a Hermitian space E of .nite dimension, for every linear map 
f : E ¡ú E there is a unique linear map f. : E ¡ú E such that 
f . (u) ¡¤ v = u ¡¤ f(v), for all u, v ¡Ê E. 
Proof. Careful inspection of the proof of Proposition 11.8 reveals that it applies unchanged. The only potential problem is in proving that f.(¦Ëu)= ¦Ëf.(u), but everything takes place in the .rst argument of the Hermitian product, and there, we have linearity. 
De.nition 13.6. Given a Hermitian space E of .nite dimension, for every linear map 
f : E ¡ú E, the unique linear map f. : E ¡ú E such that 
f . (u) ¡¤ v = u ¡¤ f(v), for all u, v ¡Ê E 
given by Proposition 13.8 is called the adjoint of f (w.r.t. to the Hermitian product). 
The fact that 
v ¡¤ u = u ¡¤ v 

implies that the adjoint f. of f is also characterized by 
f(u) ¡¤ v = u ¡¤ f . (v), 
for all u, v ¡Ê E. 
Given two Hermitian spaces E and F , where the Hermitian product on E is denoted by£¨.,.´¨1 and the Hermitian product on F is denoted by£¨.,.´¨2, given any linear map 
f : E ¡ú F , it is immediately veri.ed that the proof of Proposition 13.8 can be adapted to show that there is a unique linear map f. : F ¡ú E such that
£¨f(u),v´¨ =£¨u, f . (v)´¨
21 
for all u ¡Ê E and all v ¡Ê F . The linear map f. is also called the adjoint of f. 
As in the Euclidean case, the following properties immediately follow from the de.nition of the adjoint map. 
Proposition 13.9. (1) For any linear map f : E ¡ú F , we have 
f .. 
= f. 
(2) For any two linear maps f, g : E ¡ú F and any scalar ¦Ë ¡Ê R: 

(f + g) . = f . + g . (¦Ëf) . = ¦Ëf . . 

(3) If 	E, F, G are Hermitian spaces with respective inner products£¨.,.´¨1,£¨.,.´¨2, and£¨.,.´¨3, and if f : E ¡ú F and g : F ¡ú G are two linear maps, then 


(g . f) . = f . . g . . 
As in the Euclidean case, a linear map f : E ¡ú E (where E is a .nite-dimensional Hermitian space) is self-adjoint if f = f. . The map f is positive semide.nite i.
£¨f(x),x´¨¡Ý 0 all x ¡Ê E; 
positive de.nite i.£¨f(x),x´¨ > 0 all x ¡Ê E, x =0. 
An interesting corollary of Proposition 13.3 is that a positive semide.nite linear map must be self-adjoint. In fact, we can prove a slightly more general result. 
Proposition 13.10. Given any .nite-dimensional Hermitian space E with Hermitian prod-uct£¨.,.´¨, for any linear map f : E ¡ú E, if£¨f(x),x´¨¡Ê R for all x ¡Ê E, then f is self-adjoint. In particular, any positive semide.nite linear map f : E ¡ú E is self-adjoint. 
Proof. Since£¨f(x),x´¨¡Ê R for all x ¡Ê E, we have
£¨f(x),x´¨ =£¨f(x),x´¨ =£¨x, f(x)´¨ =£¨f . (x),x´¨, 
so we have£¨(f . f . )(x),x´¨ = 0 all x ¡Ê E, 
and Proposition 13.3 implies that f . f. = 0. 
Beware that Proposition 13.10 is false if E is a real Euclidean space. As in the Euclidean case, Theorem 13.6 can be used to show that any Hermitian space of .nite dimension has an orthonormal basis. The proof is unchanged. 
Proposition 13.11. Given any nontrivial Hermitian space E of .nite dimension n ¡Ý 1, there is an orthonormal basis (u1,...,un) for E. 

13.2. ORTHOGONALITY, DUALITY, ADJOINT OF A LINEAR MAP 
The Gram¨CSchmidt orthonormalization procedure also applies to Hermitian spaces of .nite dimension, without any changes from the Euclidean case! 
Proposition 13.12. Given a nontrivial Hermitian space E of .nite dimension n ¡Ý 1, from any basis (e1,...,en) for E we can construct an orthonormal basis (u1,...,un) for E with the property that for every k, 1 ¡Ü k ¡Ü n, the families (e1,...,ek) and (u1,...,uk) generate the same subspace. 
Remark: The remarks made after Proposition 11.10 also apply here, except that in the QR-decomposition, Q is a unitary matrix. 
As a consequence of Proposition 11.9 (or Proposition 13.12), given any Hermitian space of .nite dimension n, if (e1,...,en) is an orthonormal basis for E, then for any two vectors u = u1e1 + ¡¤¡¤¡¤ + unen and v = v1e1 + ¡¤¡¤¡¤ + vnen, the Hermitian product u ¡¤ v is expressed as 
n
C 
u ¡¤ v =(u1e1 + ¡¤¡¤¡¤ + unen) ¡¤ (v1e1 + ¡¤¡¤¡¤ + vnen)= uivi, 
i=1 
and the norm lul as 
n1/2
C 
lul = lu1e1 + ¡¤¡¤¡¤ + unenl = |ui|2 . 
i=1 
The fact that a Hermitian space always has an orthonormal basis implies that any Gram matrix G can be written as G = Q . Q, 
for some invertible matrix Q. Indeed, we know that in a change of basis matrix, a Gram matrix G becomes G ' = P .GP . If the basis corresponding to G ' is orthonormal, then G ' = I, so G =(P .1).P .1 . 
Proposition 11.11 also holds unchanged. 
Proposition 13.13. Given any nontrivial Hermitian space E of .nite dimension n ¡Ý 1, for any subspace F of dimension k, the orthogonal complement F ¡Í of F has dimension n . k, and E = F ¨’ F ¡Í . Furthermore, we have F ¡Í¡Í = F . 


13.3 	Linear Isometries (Also Called Unitary Transfor-mations) 
In this section we consider linear maps between Hermitian spaces that preserve the Hermitian norm. All de.nitions given for Euclidean spaces in Section 11.5 extend to Hermitian spaces, except that orthogonal transformations are called unitary transformation, but Proposition 
11.12 extends only with a modi.ed Condition (2). Indeed, the old proof that (2) implies 
(3) does not work, and the implication is in fact false! It can be repaired by strengthening Condition (2). For the sake of completeness, we state the Hermitian version of De.nition 
11.5. 
De.nition 13.7. Given any two nontrivial Hermitian spaces E and F of the same .nite dimension n, a function f : E ¡ú F is a unitary transformation, or a linear isometry, if it is linear and 
lf(u)l = lul, for all u ¡Ê E. 
Proposition 11.12 can be salvaged by strengthening Condition (2). 
Proposition 13.14. Given any two nontrivial Hermitian spaces E and F of the same .nite dimension n, for every function f : E ¡ú F , the following properties are equivalent: 
(1) 
f is a linear map and lf(u)l = lul, for all u ¡Ê E; 

(2) 
lf(v) . f(u)l = lv . ul and f(iu)= if(u), for all u, v ¡Ê E. 

(3) 
f(u) ¡¤ f(v)= u ¡¤ v, for all u, v ¡Ê E. 


Furthermore, such a map is bijective. 
Proof. The proof that (2) implies (3) given in Proposition 11.12 needs to be revised as follows. We use the polarization identity 
2.(u, v) = (1+ i)(lul2 + lvl2) .lu . vl2 . ilu . ivl2 . 
Since f(iv)= if(v), we get f(0) = 0 by setting v = 0, so the function f preserves distance and norm, and we get 
2.(f(u),f(v)) = (1+ i)(lf(u)l2 + lf(v)l2) .lf(u) . f(v)l2 . ilf(u) . if(v)l2 = (1+ i)(lf(u)l2 + lf(v)l2) .lf(u) . f(v)l2 
. ilf(u) . f(iv)l2 = (1+ i)(lul2 + lvl2) .lu . vl2 . ilu . ivl2 =2.(u, v), 
which shows that f preserves the Hermitian inner product as desired. The rest of the proof is unchanged. 
13.4. THE UNITARY GROUP, UNITARY MATRICES 
Remarks: 
(i) In the Euclidean case, we proved that the assumption 
lf(v) . f(u)l = lv . ul for all u, v ¡Ê E and f(0)=0 (2 ' ) 
implies (3). For this we used the polarization identity 
2u ¡¤ v = lul2 + lvl2 .lu . vl2 . 
In the Hermitian case the polarization identity involves the complex number i. In fact, the implication (2 ' ) implies (3) is false in the Hermitian case! Conjugation z ¡ú z satis.es (2 ' ) since 
|z2 . z1| = |z2 . z1| = |z2 . z1|, 
and yet, it is not linear! 
(ii) If we modify (2) by changing the second condition by now requiring that there be some 
¦Ó ¡Ê E such that 
f(¦Ó + iu)= f(¦Ó)+ i(f(¦Ó + u) . f(¦Ó)) 

for all u ¡Ê E, then the function g : E ¡ú E de.ned such that 
g(u)= f(¦Ó + u) . f(¦Ó) 
satis.es the old conditions of (2), and the implications (2) ¡ú (3) and (3) ¡ú (1) prove that g is linear, and thus that f is a.ne. In view of the .rst remark, some condition involving i is needed on f, in addition to the fact that f is distance-preserving. 


13.4 The Unitary Group, Unitary Matrices 
In this section, as a mirror image of our treatment of the isometries of a Euclidean space, we explore some of the fundamental properties of the unitary group and of unitary matrices. As an immediate corollary of the Gram¨CSchmidt orthonormalization procedure, we obtain the QR-decomposition for invertible matrices. In the Hermitian framework, the matrix of the adjoint of a linear map is not given by the transpose of the original matrix, but by its conjugate. 
De.nition 13.8. Given a complex m ¡Á n matrix A, the transpose AT of A is the n ¡Á m
3  
matrix AT = aT de.ned such that 
ij
T 
a = aji,
ij 
and the conjugate A of A is the m ¡Á n matrix A =(bij) de.ned such that 
bij = aij 
for all i, j,1 ¡Ü i ¡Ü m,1 ¡Ü j ¡Ü n. The adjoint A. of A is the matrix de.ned such that 
3 T 
A . =(AT)= A. 
Proposition 13.15. Let E be any Hermitian space of .nite dimension n, and let f : E ¡ú E be any linear map. The following properties hold: 
(1) 
The linear map f : E ¡ú E is an isometry i. 

f . f . = f . . f = id. 

(2) 
For every orthonormal basis (e1,...,en) of E, if the matrix of f is A, then the matrix of f. is the adjoint A. of A, and f is an isometry i. A satis.es the identities 


AA . = A . A = In, 
where In denotes the identity matrix of order n, i. the columns of A form an orthonor-mal basis of Cn, i. the rows of A form an orthonormal basis of Cn . 
Proof. (1) The proof is identical to that of Proposition 11.14 (1). 
(2) If (e1,...,en) is an orthonormal basis for E, let A =(aij) be the matrix of f, and let B =(bij) be the matrix of f. . Since f. is characterized by 
f . (u) ¡¤ v = u ¡¤ f(v) 
for all u, v ¡Ê E, using the fact that if w = w1e1 + ¡¤¡¤¡¤ + wnen, we have wk = w ¡¤ ek, for all k, 1 ¡Ü k ¡Ü n; letting u = ei and v = ej, we get 
bji = f . (ei) ¡¤ ej = ei ¡¤ f(ej)= f(ej) ¡¤ ei = aij, 
for all i, j,1 ¡Ü i, j ¡Ü n. Thus, B = A. . Now if X and Y are arbitrary matrices over the basis (e1,...,en), denoting as usual the jth column of X by Xj, and similarly for Y , a simple calculation shows that 
Y . X =(Xj ¡¤ Y i)1¡Üi,j¡Ün. 
Then it is immediately veri.ed that if X = Y = A, then A.A = AA. = In i. the column vectors (A1,...,An) form an orthonormal basis. Thus, from (1), we see that (2) is clear. 
Proposition 11.14 shows that the inverse of an isometry f is its adjoint f. . Proposition 
11.14 also motivates the following de.nition. 
De.nition 13.9. A complex n ¡Á n matrix is a unitary matrix if 
AA . = A . A = In. 
13.4. THE UNITARY GROUP, UNITARY MATRICES 
Remarks: 
(1) 
The conditions AA. = In, A.A = In, and A.1 = A. are equivalent. Given any two orthonormal bases (u1,...,un) and (v1,...,vn), if P is the change of basis matrix from (u1,...,un) to (v1,...,vn), it is easy to show that the matrix P is unitary. The proof of Proposition 13.14 (3) also shows that if f is an isometry, then the image of an orthonormal basis (u1,...,un) is an orthonormal basis. 

(2) 
Using the explicit formula for the determinant, we see immediately that 


det(A) = det(A). 
If f is a unitary transformation and A is its matrix with respect to any orthonormal basis, from AA. = I, we get 
det(AA . ) = det(A) det(A . ) = det(A)det(AT) = det(A)det(A)= | det(A)|2 , 
and so | det(A)| = 1. It is clear that the isometries of a Hermitian space of dimension n form a group, and that the isometries of determinant +1 form a subgroup. 
This leads to the following de.nition. 
De.nition 13.10. Given a Hermitian space E of dimension n, the set of isometries f : E ¡ú E forms a subgroup of GL(E, C) denoted by U(E), or U(n) when E = Cn , called the unitary group (of E). For every isometry f we have | det(f)| = 1, where det(f) denotes the determinant of f. The isometries such that det(f) = 1 are called rotations, or proper isometries, or proper unitary transformations, and they form a subgroup of the special linear group SL(E, C) (and of U(E)), denoted by SU(E), or SU(n) when E = Cn, called the special unitary group (of E). The isometries such that det(f) = 1 are called improper isometries, or improper unitary transformations, or .ip transformations. 
A very important example of unitary matrices is provided by Fourier matrices (up to a 
¡Ì 
factor of n), matrices that arise in the various versions of the discrete Fourier transform. For more on this topic, see the problems, and Strang [63, 66]. 
The group SU(2) turns out to be the group of unit quaternions, invented by Hamilton. This group plays an important role in the representation of rotations in SO(3) used in computer graphics and robotics; see Chapter 15. 
Now that we have the de.nition of a unitary matrix, we can explain how the Gram¨C Schmidt orthonormalization procedure immediately yields the QR-decomposition for matri-ces. 
De.nition 13.11. Given any complex n ¡Án matrix A,a QR-decomposition of A is any pair of n ¡Á n matrices (U, R), where U is a unitary matrix and R is an upper triangular matrix such that A = UR. 
Proposition 13.16. Given any n ¡Á n complex matrix A, if A is invertible, then there is a unitary matrix U and an upper triangular matrix R with positive diagonal entries such that A = UR. 
The proof is absolutely the same as in the real case! 
Remark: If A is invertible and if A = U1R1 = U2R2 are two QR-decompositions for A, then 
R1R2 .1 = U1 . U2. Then it is easy to show that there is a diagonal matrix D with diagonal entries such that |dii| = 1 for i =1,...,n, and U2 = U1D, R2 = D.R1. We have the following version of the Hadamard inequality for complex matrices. The proof is essentially the same as in the Euclidean case but it uses Proposition 13.16 instead of Proposition 11.16. 
Proposition 13.17. (Hadamard) For any complex n ¡Á n matrix A =(aij), we have 
nn1/2 nn1/2
nC nC 
| det(A)|¡Ü |aij|2 and | det(A)|¡Ü |aij|2 . 
i=1 j=1 j=1 i=1 
Moreover, equality holds i. either A has a zero row in the left inequality or a zero column in the right inequality, or A is unitary. 
We also have the following version of Proposition 11.18 for Hermitian matrices. The proof of Proposition 11.18 goes through because the Cholesky decomposition for a Hermitian positive de.nite A matrix holds in the form A = B.B, where B is upper triangular with positive diagonal entries. The details are left to the reader. 
Proposition 13.18. (Hadamard) For any complex n¡Án matrix A =(aij), if A is Hermitian positive semide.nite, then we have 
n
n 
det(A) ¡Ü aii. i=1 
Moreover, if A is positive de.nite, then equality holds i. A is a diagonal matrix. 
13.5 Hermitian Re.ections and QR-Decomposition 
If A is an n ¡Á n complex singular matrix, there is some (not necessarily unique) QR-decomposition A = QR with Q a unitary matrix which is a product of Householder re-.ections and R an upper triangular matrix, but the proof is more involved. One way to proceed is to generalize the notion of hyperplane re.ection. This is not really surprising since in the Hermitian case there are improper isometries whose determinant can be any unit complex number. Hyperplane re.ections are generalized as follows. 

13.5. HERMITIAN REFLECTIONS AND QR-DECOMPOSITION 
De.nition 13.12. Let E be a Hermitian space of .nite dimension. For any hyperplane H, for any nonnull vector w orthogonal to H, so that E = H ¨’ G, where G = Cw,a Hermitian re.ection about H of angle ¦È is a linear map of the form ¦ÑH, ¦È : E ¡ú E, de.ned such that 
¦ÑH, ¦È(u)= pH (u)+ e i¦È pG(u), 
for any unit complex number ei¦È = 1 (i.e. ¦È = k2¦Ð). For any nonzero vector w ¡Ê E, we denote by ¦Ñw,¦È the Hermitian re.ection given by ¦ÑH,¦È, where H is the hyperplane orthogonal to w. 
Since u = pH (u)+ pG(u), the Hermitian re.ection ¦Ñw, ¦È is also expressed as 
i¦È . 1)pG(u),¦Ñw, ¦È(u)= u +(e 
or as (u ¡¤ w)
¦Ñw, ¦È(u)= u +(e i¦È . 1) w. 
lwl2 
Note that the case of a standard hyperplane re.ection is obtained when ei¦È = .1, i.e., ¦È = ¦Ð. In this case, 
(u ¡¤ w)
¦Ñw, ¦Ð(u)= u . 2 w, 
lwl2 
and the matrix of such a re.ection is a Householder matrix, as in Section 12.1, except that w may be a complex vector. 
We leave as an easy exercise to check that ¦Ñw, ¦È is indeed an isometry, and that the inverse of ¦Ñw, ¦È is ¦Ñw, .¦È. If we pick an orthonormal basis (e1,...,en) such that (e1,...,en.1) is an orthonormal basis of H, the matrix of ¦Ñw, ¦È is 
In.1 0 
i¦È
0 e
We now come to the main surprise. Given any two distinct vectors u and v such that lul = lvl, there isn¡¯t always a hyperplane re.ection mapping u to v, but this can be done using two Hermitian re.ections! 
Proposition 13.19. Let E be any nontrivial Hermitian space. 
(1) 
For any two vectors u, v ¡Ê E such that u = v and lul = lvl, if u ¡¤ v = ei¦È|u ¡¤ v|, then the (usual) re.ection s about the hyperplane orthogonal to the vector v . e.i¦Èu is such that s(u)= ei¦Èv. 

(2) 
For any nonnull vector v ¡Ê E, for any unit complex number ei¦È =1, there is a Hermi-


tian re.ection ¦Ñv,¦È such that 
¦Ñv,¦È(v)= e i¦È v. 

As a consequence, for u and v as in (1), we have ¦Ñv,.¦È . s(u)= v. 
Proof. (1) Consider the (usual) re.ection about the hyperplane orthogonal to w = v .e.i¦Èu. We have 
(u ¡¤ (v . e.i¦Èu)) 
s(u)= u . 2(v . e .i¦È u). 
.i¦Èul2
lv . e
We need to compute 
.i¦È .i¦È .i¦È
.2u ¡¤ (v . eu) and (v . eu) ¡¤ (v . eu). Since u ¡¤ v = ei¦È|u ¡¤ v|, we have 
.i¦È i¦È 
eu ¡¤ v = |u ¡¤ v| and ev ¡¤ u = |u ¡¤ v|. Using the above and the fact that lul = lvl, we get .i¦È i¦È lul2
.2u ¡¤ (v . eu)=2e . 2u ¡¤ v, i¦È(lul2
=2e .|u ¡¤ v|), and 
.i¦È .i¦È .i¦È i¦È
(v . eu) ¡¤ (v . eu)= lvl2 + lul2 . eu ¡¤ v . ev ¡¤ u, = 2(lul2 .|u ¡¤ v|), and thus, (u ¡¤ (v . e.i¦Èu)) .i¦È i¦È(v . e .i¦È
.2(v . eu)= eu). 
.i¦È
l(v . eu)l2 But then, 
i¦È(v . e .i¦È i¦Èi¦È 
s(u)= u + eu)= u + ev . u = e v, and s(u)= ei¦Èv, as claimed. 
(2) This part is easier. Consider the Hermitian re.ection 
(u ¡¤ v)

¦Ñv,¦È(u)= u +(e i¦È . 1) v. 
lvl2 
We have (v ¡¤ v)
¦Ñv,¦È(v)= v +(e i¦È . 1) v, 
lvl2 
= v +(e i¦È . 1)v, i¦È 
= e v. Thus, ¦Ñv,¦È(v)= ei¦Èv. Since ¦Ñv,¦È is linear, changing the argument v to ei¦Èv, we get 
¦Ñv,.¦È(e i¦È v)= v, 
and thus, ¦Ñv,.¦È . s(u)= v. 

13.5. HERMITIAN REFLECTIONS AND QR-DECOMPOSITION 
Remarks: 
.i¦È.i¦È	i¦È
(1) If we use the vector v + eu instead of v . eu, we get s(u)= .ev. 

(2) Certain authors, 	such as Kincaid and Cheney [39] and Ciarlet [14], use the vector u + ei¦Èv instead of the vector v + e.i¦Èu. The e.ect of this choice is that they also get s(u)= .ei¦Èv. 

(3) If v = lul e1, where e1 is a basis vector, u ¡¤ e1 = a1, where a1 is just the coe.cient of u over the basis vector e1. Then, since u ¡¤ e1 = ei¦È|a1|, the choice of the plus sign in the vector lul e1 + e.i¦Èu has the e.ect that the coe.cient of this vector over e1 is lul + |a1|, and no cancellations takes place, which is preferable for numerical stability (we need to divide by the square norm of this vector). 


We now show that the QR-decomposition in terms of (complex) Householder matrices holds for complex matrices. We need the version of Proposition 13.19 and a trick at the end of the argument, but the proof is basically unchanged. 
Proposition 13.20. Let E be a nontrivial Hermitian space of dimension n. Given any orthonormal basis (e1,...,en), for any n-tuple of vectors (v1,...,vn), there is a sequence of n . 1 isometries h1,...,hn.1, such that hi is a (standard) hyperplane re.ection or the identity, and if (r1,...,rn) are the vectors given by 
½Ð½Ð 
rj = hn.1 .¡¤ ¡¤¡¤. h2 . h1(vj), 1 ¡Ü j ¡Ü n, 
½Ð½Ð 
then every rj is a linear combination of the vectors (e1,...,ej), (1 ¡Ü j ¡Ü n). Equivalently, the matrix R whose columns are the components of the rj over the basis (e1,...,en) is an upper triangular matrix. Furthermore, if we allow one more isometry hn of the form 
hn = ¦Ñen,.n .¡¤ ¡¤¡¤. ¦Ñe1,.1 
after h1,...,hn.1, we can ensure that the diagonal entries of R are nonnegative. 
Proof. The proof is very similar to the proof of Proposition 12.3, but it needs to be modi.ed a little bit since Proposition 13.19 is weaker than Proposition 12.2. We explain how to modify the induction step, leaving the base case and the rest of the proof as an exercise. 
As in the proof of Proposition 12.3, the vectors (e1,...,ek) form a basis for the subspace denoted as Uk' , the vectors (ek+1,...,en) form a basis for the subspace denoted as Uk '' , the 
' '' 	''' 
subspaces U and U are orthogonal, and E = U ¨’ Uk . Let
kk 	k 
uk+1 = hk .¡¤ ¡¤¡¤. h2 . h1(vk+1). 
We can write 
' '' 
uk+1 = uk+1 + uk+1, 
' ''''' 
where u ¡Ê U and u ¡Ê U Let
k+1 kk+1 k . 
'' 	i¦Èk+1 |u '' '' 
rk+1,k+1 =uk+1, and e k+1 ¡¤ ek+1| = uk+1 ¡¤ ek+1. 
If u '' k+1 = ei¦Èk+1 rk+1,k+1 ek+1, we let hk+1 = id. Otherwise, by Proposition 13.19(1) (with 
u = uk'' +1 and v = rk+1,k+1 ek+1), there is a unique hyperplane re.ection hk+1 such that 
'' i¦Èk+1
hk+1(uk+1)= erk+1,k+1 ek+1, 
where hk+1 is the re.ection about the hyperplane Hk+1 orthogonal to the vector 
.i¦Èk+1 '' 
wk+1 = rk+1,k+1 ek+1 . euk+1. 
At the end of the induction, we have a triangular matrix R, but the diagonal entries ei¦Èj rj, j of R may be complex. Letting 
hn = ¦Ñen, .¦Èn .¡¤ ¡¤¡¤. ¦Ñe1,.¦È1 , 
we observe that the diagonal entries of the matrix of vectors 
r ' = hn . hn.1 .¡¤ ¡¤¡¤. h2 . h1(vj)
j 
is triangular with nonnegative entries. 
.i¦Èk+1 u
Remark: For numerical stability, it is preferable to use wk+1 = rk+1,k+1 ek+1 + e'' k+1 
.i¦Èk+1 u '' 
instead of wk+1 = rk+1,k+1 ek+1 . ek+1. The e.ect of that choice is that the diagonal i(¦Èj +¦Ð)
entries in R will be of the form .ei¦Èj rj, j = erj, j. Of course, we can make these entries nonegative by applying 
hn = ¦Ñen,¦Ð.¦Èn .¡¤ ¡¤¡¤. ¦Ñe1,¦Ð.¦È1 
after hn.1. 
As in the Euclidean case, Proposition 13.20 immediately implies the QR-decomposition for arbitrary complex n ¡Á n-matrices, where Q is now unitary (see Kincaid and Cheney [39] and Ciarlet [14]). 
Proposition 13.21. For every complex n ¡Á n-matrix A, there is a sequence H1,...,Hn.1 of matrices, where each Hi is either a Householder matrix or the identity, and an upper triangular matrix R, such that 
R = Hn.1 ¡¤¡¤¡¤ H2H1A. 
As a corollary, there is a pair of matrices Q, R, where Q is unitary and R is upper triangular, such that A = QR (a QR-decomposition of A). Furthermore, R can be chosen so that its diagonal entries are nonnegative. This can be achieved by a diagonal matrix D with entries 
R
such that |dii| =1 for i =1,...,n, and we have A = QR 
with 
R¡¤¡¤ Hn.1D, R 
= D . R,
Q = H1 ¡¤ 
where R 
is upper triangular and has nonnegative diagonal entries. 
Proof. It is essentially identical to the proof of Proposition 12.4, and we leave the details as an exercise. For the last statement, observe that hn .¡¤ ¡¤¡¤. h1 is also an isometry. 

13.6. ORTHOGONAL PROJECTIONS AND INVOLUTIONS 


13.6 Orthogonal Projections and Involutions 
In this section we begin by assuming that the .eld K is not a .eld of characteristic 2. Recall that a linear map f : E ¡ú E is an involution i. f2 = id, and is idempotent i. f2 = f. We know from Proposition 5.7 that if f is idempotent, then 
E = Im(f) ¨’ Ker (f), 
and that the restriction of f to its image is the identity. For this reason, a linear involution is called a projection. The connection between involutions and projections is given by the following simple proposition. 
Proposition 13.22. For any linear map f : E ¡ú E, we have f2 = id i. 21 (id . f) is a projection i. 21 (id + f) is a projection; in this case, f is equal to the di.erence of the two projections 12 (id + f) and 12 (id . f). 
Proof. We have 
1(id . f) 2 = 1(id . 2f + f2)
24so 1(id . f) 2 = 1(id . f) i. f2 = id. 
22We also have 1(id + f) 2 = 1(id + 2f + f2),
24so 1(id + f) 2 = 1(id + f) i. f2 = id. 
22Obviously, f = 21 (id + f) . 12 (id . f). Proposition 13.23. For any linear map f : E ¡ú E, let U+ = Ker (12 (id . f)) and let U. = Im(12 (id . f)). If f2 = id, then 
11 
U+ = Ker (id . f) =Im (id+ f) ,
22and so, f(u)= u on U+ and f(u)= .u on U. . Proof. If f2 = id, then 
(id + f) . (id . f) = id . f2 = id . id = 0, which implies that 
11 
Im (id+ f) . Ker (id . f) . 
22
3 
Conversely, if u ¡Ê Ker 12 (id . f) , then f(u)= u, so 11 
(id + f)(u)= (u + u)= u,
22and thus 
11 
Ker (id . f) . Im (id+ f) . 
22Therefore, 
11 
U+ = Ker (id . f) =Im (id+ f) ,
22and so, f(u)= u on U+ and f(u)= .u on U. . 
We now assume that K = C. The involutions of E that are unitary transformations are characterized as follows. 
Proposition 13.24. Let f ¡Ê GL(E) be an involution. The following properties are equiva-lent: 
(a) The map f is unitary; that is, f ¡Ê U(E). 

(b) The subspaces U. = Im(12 (id . f)) and U+ = Im(12 (id + f)) are orthogonal. 
Furthermore, if E is .nite-dimensional, then (a) and (b) are equivalent to 



(c) The map is self-adjoint; that is, f = f. . Proof. If f is unitary, then from£¨f(u),f(v)´¨ =£¨u, v´¨ for all u, v ¡Ê E, we see that if u ¡Ê U+ 
and v ¡Ê U., we get£¨u, v´¨ =£¨f(u),f(v)´¨ =£¨u, .v´¨ =.£¨u, v´¨, 
so 2£¨u, v´¨ = 0, which implies£¨u, v´¨ = 0, that is, U+ and U. are orthogonal. Thus, (a) 
implies (b). 

Conversely, if (b) holds, since f(u)= u on U+ and f(u)= .u on U. , we see that£¨f(u),f(v)´¨ =£¨u, v´¨ if u, v ¡Ê U+ or if u, v ¡Ê U. . Since E = U+ ¨’ U. and since U+ and U. are orthogonal, we also have£¨f(u),f(v)´¨ =£¨u, v´¨ for all u, v ¡Ê E, and (b) implies (a). 
If E is .nite-dimensional, the adjoint f. of f exists, and we know that f.1 = f. . Since f is an involution, f2 = id, which implies that f. = f.1 = f. 
A unitary involution is the identity on U+ = Im(12 (id + f)), and f(v)= .v for all v ¡Ê U. = Im(21 (id . f)). Furthermore, E is an orthogonal direct sum E = U+ ¨’U. . We say that f is an orthogonal re.ection about U+ . In the special case where U+ is a hyperplane, we say that f is a hyperplane re.ection. We already studied hyperplane re.ections in the Euclidean case; see Chapter 12. 
If f : E ¡ú E is a projection (f2 = f), then 
(id . 2f)2 = id . 4f +4f2 = id . 4f +4f = id, 
so id . 2f is an involution. As a consequence, we get the following result. 
13.7. DUAL NORMS 
Proposition 13.25. If f : E ¡ú E is a projection (f2 = f), then Ker (f) and Im(f) are orthogonal i. f. = f. 
Proof. Apply Proposition 13.24 to g = id . 2f. Since id . g =2f we have 
1 
U+ = Ker (id . g) = Ker(f)
2
and U. = Im 1(id . g) = Im(f),
2which proves the proposition. 
A projection such that f = f. is called an orthogonal projection. 
If (a1 ...,ak) are k linearly independent vectors in Rn, let us determine the matrix P of the orthogonal projection onto the subspace of Rn spanned by (a1,...,ak). Let A be the n¡Ák matrix whose jth column consists of the coordinates of the vector aj over the canonical basis (e1,...,en). 
Any vector in the subspace (a1,...,ak) is a linear combination of the form Ax, for some x ¡Ê Rk . Given any y ¡Ê Rn, the orthogonal projection Py = Ax of y onto the subspace spanned by (a1,...,ak) is the vector Ax such that y . Ax is orthogonal to the subspace spanned by (a1,...,ak) (prove it). This means that y . Ax is orthogonal to every aj, which is expressed by 
AT(y . Ax) = 0; 
that is, ATAx = AT y. 
The matrix ATA is invertible because A has full rank k, thus we get 
x =(ATA).1AT y, 
and so Py = Ax = A(ATA).1AT y. 
Therefore, the matrix P of the projection onto the subspace spanned by (a1 ...,ak) is given by 
P = A(ATA).1AT . 
The reader should check that P 2 = P and P T = P . 


13.7 Dual Norms 
In the remark following the proof of Proposition 8.10, we explained that if (E, ll) and (F, ll) are two normed vector spaces and if we let L(E; F ) denote the set of all continuous (equivalently, bounded) linear maps from E to F , then, we can de.ne the operator norm (or subordinate norm) ll on L(E; F ) as follows: for every f ¡ÊL(E; F ), 
lf(x)l 
lfl = sup = sup lf(x)l . 
x¡ÊE lxl x¡ÊE 
x lxl=1
=0 
In particular, if F = C, then L(E; F )= E ' is the dual space of E, and we get the operator norm denoted by ll. given by 
lfl. = sup |f(x)|. 
x¡ÊE lxl=1 
The norm ll. is called the dual norm of ll on E ' . 
Let us now assume that E is a .nite-dimensional Hermitian space, in which case E ' = E. . Theorem 13.6 implies that for every linear form f ¡Ê E., there is a unique vector y ¡Ê E so that 
f(x)=£¨x, y´¨, 
for all x ¡Ê E, and so we can write 
lfl. = sup |£¨x, y´¨|. 
x¡ÊE lxl=1
The above suggests de.ning a norm llD on E. 
De.nition 13.13. If E is a .nite-dimensional Hermitian space and ll is any norm on E, for any y ¡Ê E we let 
lylD = sup |£¨x, y´¨|, 
x¡ÊE lxl=1
be the dual norm of ll (on E). If E is a real Euclidean space, then the dual norm is de.ned by 
lylD = sup £¨x, y´¨ 
x¡ÊE lxl=1
for all y ¡Ê E. 
Beware that ll is generally not the Hermitian norm associated with the Hermitian inner product. The dual norm shows up in convex programming; see Boyd and Vandenberghe [11], Chapters 2, 3, 6, 9. 
The fact that llD is a norm follows from the fact that ll. is a norm and can also be checked directly. It is worth noting that the triangle inequality for llD comes ¡°for free,¡± in the sense that it holds for any function p: E ¡ú R. 
Proposition 13.26. For any function p: E ¡ú R, if we de.ne pD by 
p D(x) = sup |£¨z, x´¨|, 
p(z)=1
13.7. DUAL NORMS 
then we have 
p D(x + y) ¡Ü p D(x)+ p D(y). Proof. We have 
p D(x + y) = sup |£¨z, x + y´¨|
p(z)=1
= sup (|£¨z, x´¨ +£¨z, y´¨|) 
p(z)=1 
¡Ü sup (|£¨z, x´¨| +|£¨z, y´¨|) 
p(z)=1 
¡Ü sup |£¨z, x´¨| + sup |£¨z, y´¨|
p(z)=1p(z)=1= p D(x)+ p D(y). 
De.nition 13.14. If p : E ¡ú R is a function such that 
(1) p(x) ¡Ý 0 for all x ¡Ê E, and p(x)=0 i. x = 0; 

(2) p(¦Ëx)= |¦Ë|p(x), for all x ¡Ê E and all ¦Ë ¡Ê C; 

(3) p is continuous, in the sense that for some basis (e1,...,en) of E, the function 


(x1,...,xn) ¡ú p(x1e1 + ¡¤¡¤¡¤ + xnen) 
from Cn to R is continuous, then we say that p is a pre-norm. Obviously, every norm is a pre-norm, but a pre-norm may not satisfy the triangle in-equality. 
Corollary 13.27. The dual norm of any pre-norm is actually a norm. Proposition 13.28. For all y ¡Ê E, we have 
lylD = sup |£¨x, y´¨| = sup R£¨x, y´¨. 
x¡ÊEx¡ÊE lxl=1lxl=1
Proof. Since E is .nite dimensional, the unit sphere Sn.1 = {x ¡Ê E |lxl =1} is compact, so there is some x0 ¡Ê Sn.1 such that 
lylD =|£¨x0,y´¨|. 
If£¨x0,y´¨ = ¦Ñei¦È, with ¦Ñ ¡Ý 0, then
|£¨e .i¦È x0,y´¨| = |e .i¦È£¨x0,y´¨| = |e .i¦È¦Ñei¦È| = ¦Ñ, 
so 
lylD = ¦Ñ =£¨e .i¦È x0,y´¨, (.) withe.i¦Èx0= lx0l = 1. On the other hand,
½Ð½Ð ½Ð½Ð 
R£¨x, y´¨ ¡Ü|£¨ x, y´¨|,  
so by (.) we get  
lylD  =  sup x¡ÊE  |£¨x, y´¨| =  sup x¡ÊE  R£¨x, y´¨,  
lxl=1 lxl=1 
as claimed.  

Proposition 13.29. For all x, y ¡Ê E, we have
|£¨x, y´¨| ¡Ü l xllylD
|£¨x, y´¨| ¡Ü l xlD lyl . 
Proof. If x = 0, then£¨x, y´¨ = 0 and these inequalities are trivial. If x = 0, since lx/ lxll = 1, 

by de.nition of lylD, we have|£¨x/ lxl ,y´¨| ¡Ü sup |£¨z, y´¨| = lylD , 
lzl=1
which yields
|£¨x, y´¨| ¡Ü l xllylD . The second inequality holds because|£¨x, y´¨| =|£¨y, x´¨|. It is not hard to show that for all y ¡Ê Cn , lylD = lyl¡Þ
1 
lylD = lyl
¡Þ 1 lylD 2 = lyl2 . 
Thus, the Euclidean norm is autodual. More generally, the following proposition holds. 
Proposition 13.30. If p, q ¡Ý 1 and 1/p +1/q =1, then for all y ¡Ê Cn, we have 
lylDp = lyl q . 
Proof. By H¡§older¡¯s inequality (Corollary 8.2), for all x, y ¡Ê Cn, we have|£¨x, y´¨| ¡Ü l xllyl ,
pq 
so lylD = sup |£¨x, y´¨| ¡Ü l yl . 
pq 
x¡ÊCn lxl =1
p 

13.7. DUAL NORMS 
For the converse, we consider the cases p =1, 1 <p< +¡Þ, and p =+¡Þ. First assume p = 1. The result is obvious for y = 0, so assume y = 0. Given y, if we pick xj =1 for some index j such that lyl¡Þ = max1¡Üi¡Ün |yi| = |yj|, and xk = 0 for k = j, then|£¨x, y´¨| = |yj| = lyl¡Þ, so lylD = lyl¡Þ.
1 
Now we turn to the case 1 <p< +¡Þ. Then we also have 1 <q< +¡Þ, and the equation 1/p +1/q = 1 is equivalent to pq = p + q, that is, p(q . 1) = q. Pick zj = yj|yj|q.2 for j =1,...,n, so that 
 

 1/p
 

 1/p
 

 1/pnnn
lzl ==|yj|(q.1)p=|yj|q. 
j=1 j=1 j=1 
Then if x = z/ lzl p, we have
nn 
zjyjyjyj|yj|q.2
j=1 j=1 
  
 C
C p||zjp 
o 
  
 C   
 
o 
  
 
|£¨x, y´¨| =
=

lzllzl 
pp
 
o
 1/q
C 
o 
n
 1/p
n 
|yj|qj=1 
j=1 
n 
|yj|q
j=1 
|yj|q= lyl
=

 

=

. 
q 
Thus lylDp = lyl q . 
Finally, if p = ¡Þ, then pick xj 

C
C 
= yj/|yj| if yj = 0, and xj =0 if yj = 0. Then
n
|£¨x, y´¨| =
yjyj/|yj|
    
  
= |yj| = lyl1 . 
      

yj =0 yj =0 
Thus lylD = lyl1.
¡Þ 
We can show that the dual of the spectral norm is the trace norm (or nuclear norm) also discussed in Section 20.5. Recall from Proposition 8.10 that the spectral norm lAl2 of a matrix A is the square root of the largest eigenvalue of A.A, that is, the largest singular value of A. 
Proposition 13.31. The dual of the spectral norm is given by 
lAlD ¡¤ 
2 = ¦Ò1 + ¡¤¡¤ + ¦Òr, 
where ¦Ò1 > ¡¤¡¤¡¤ >¦Òr > 0 are the singular values of A ¡Ê Mn(C) (which has rank r). 
Proof. In this case the inner product on Mn(C) is the Frobenius inner product£¨A, B´¨ = tr(B.A), and the dual norm of the spectral norm is given by 
lAlD = sup{|tr(A . B)||lBl=1}.
22 
If we factor A using an SVD as A = V ¦²U., where U and V are unitary and ¦² is a diagonal matrix whose r nonzero entries are the singular values ¦Ò1 > ¡¤¡¤¡¤ >¦Òr > 0, where r is the rank of A, then 
|tr(A . B)| = |tr(U¦²V . B)| = |tr(¦²V . BU)|, 
so if we pick B = VU., a unitary matrix such that lBl2 = 1, we get 
|tr(A . B)| = tr(¦²) = ¦Ò1 + ¡¤¡¤¡¤ + ¦Òr, 
and thus 
lAlD 2 ¡Ý ¦Ò1 + ¡¤¡¤¡¤ + ¦Òr. 
Since lBl2 = 1 and U and V are unitary, by Proposition 8.10 we have lV .BUl2 = lBl2 = 1. If Z = V .BU, by de.nition of the operator norm 
1= lZl2 = sup{lZxl2 |lxl2 =1}, 
so by picking x to be the canonical vector ej, we see that lZjl2 ¡Ü 1 where Zj is the jth column of Z, so |zjj|¡Ü 1, and since 
rrr
CC C 
|tr(¦²V . BU)| = |tr(¦²Z)| = ¦Òjzjj ¡Ü ¦Òj|zjj|¡Ü ¦Òj, 
j=1 j=1 j=1 
and we conclude that 
r
C 
|tr(¦²V . BU)|¡Ü ¦Òj. 
j=1 
The above implies that 
lAlD 2 ¡Ü ¦Ò1 + ¡¤¡¤¡¤ + ¦Òr, 
and since we also have lAlD ¡¤¡¤ + ¦Òr, we conclude that 
2 ¡Ý ¦Ò1 + ¡¤ 
lAlD 2 = ¦Ò1 + ¡¤¡¤¡¤ + ¦Òr, 
proving our proposition. 
De.nition 13.15. Given any complex matrix n ¡Á n matrix A of rank r, its nuclear norm (or trace norm) is given by 
lAlN = ¦Ò1 + ¡¤¡¤¡¤ + ¦Òr. 
The nuclear norm can be generalized to m ¡Á n matrices (see Section 20.5). The nuclear norm ¦Ò1 + ¡¤¡¤¡¤ + ¦Òr of an m ¡Á n matrix A (where r is the rank of A) is denoted by lAlN . The nuclear norm plays an important role in matrix completion. The problem is this. Given a matrix A0 with missing entries (missing data), one would like to .ll in the missing entries in A0 to obtain a matrix A of minimal rank. For example, consider the matrices 
12 1 . 12 
A0 = ,B0 = ,C0 = .
.. . 43 . 

13.7. DUAL NORMS 
All can be completed with rank 1. For A0, use any multiple of (1, 2) for the second row. For B0, use any numbers b and c such that bc = 4. For C0, the only possibility is d = 6. 
A famous example of this problem is the Net.ix competition. The ratings of m .lms by n viewers goes into A0. But the customers didn¡¯t see all the movies. Many ratings were missing. Those had to be predicted by a recommender system. The nuclear norm gave a good solution that needed to be adjusted for human psychology. 
Since the rank of a matrix is not a norm, in order to solve the matrix completion problem we can use the following ¡°convex relaxation.¡± Let A0 be an incomplete m ¡Á n matrix: 
Minimize lAlN subject to A = A0 in the known entries. 
The above problem has been extensively studied, in particular by Cand`es and Recht. Roughly, they showed that if A is an n ¡Á n matrix of rank r and K entries are known in A, then if K is large enough (K >Cn5/4r log n), with high probability, the recovery of A is perfect. See Strang [65] for details (Section III.5). 
We close this section by stating the following duality theorem. 
Theorem 13.32. If E is a .nite-dimensional Hermitian space, then for any norm ll on E, we have 
lylDD 
= lyl for all y ¡Ê E. Proof. By Proposition 13.29, we have
|£¨x, y´¨| ¡Ü l xlD lyl , 
so we get lylDD = sup |£¨x, y´¨|¡Ü l yl , for all y ¡Ê E. 
lxlD=1
It remains to prove that lyl¡ÜlylDD , for all y ¡Ê E. 
Proofs of this fact can be found in Horn and Johnson [36] (Section 5.5), and in Serre [57] (Chapter 7). The proof makes use of the fact that a nonempty, closed, convex set has a supporting hyperplane through each of its boundary points, a result known as Minkowski¡¯s lemma. For a geometric interpretation of supporting hyperplane see Figure 13.1. This result is a consequence of the Hahn¨CBanach theorem; see Gallier [25]. We give the proof in the case where E is a real Euclidean space. Some minor modi.cations have to be made when dealing with complex vector spaces and are left as an exercise. 
Since the unit ball B = {z ¡Ê E |lzl¡Ü 1} is closed and convex, the Minkowski lemma says for every x such that lxl = 1, there is an a.ne map g of the form 
g(z)=£¨z, w´¨.£¨ x, w´¨ 
with lwl = 1, such that g(x)=0 and g(z) ¡Ü 0 for all z such that lzl¡Ü 1. Then it is clear that 
sup £¨z, w´¨ =£¨x, w´¨, 
lzl=1

Figure 13.1: The orange tangent plane is a supporting hyperplane to the unit ball in R3 since this ball is entirely contained in ¡°one side¡± of the tangent plane. 
and so lwlD =£¨x, w´¨. 
It follows that 
lxlDD£¨x, w´¨ 
¡Ý£¨ w/ lwlD ,x´¨ =lwlD =1= lxl 
for all x such that lxl = 1. By homogeneity, this is true for all y ¡Ê E, which completes the proof in the real case. When E is a complex vector space, we have to view the unit ball B as a closed convex set in R2n and we use the fact that there is real a.ne map of the form 
g(z)=R£¨z, w´¨ .R£¨ x, w´¨ 
such that g(x)=0 and g(z) ¡Ü 0 for all z with lzl = 1, so that lwlD =R£¨x, w´¨. 
More details on dual norms and unitarily invariant norms can be found in Horn and Johnson [36] (Chapters 5 and 7). 


13.8 Summary 
The main concepts and results of this chapter are listed below: 
. Semilinear maps. 

. Sesquilinear forms; Hermitian forms. 

. Quadratic form associated with a sesquilinear form. 

. Polarization identities. 

13.8. SUMMARY 

. 	
Positive and positive de.nite Hermitian forms; pre-Hilbert spaces, Hermitian spaces. 

. 	
Gram matrix associated with a Hermitian product. 

. 	
The Cauchy¨CSchwarz inequality and the Minkowski inequality. 

. 	
Hermitian inner product, Hermitian norm. 

. 	
The parallelogram law. 

. 	
The musical isomorphisms b: E ¡ú E. and : E. ¡ú E; Theorem 13.6 (E is .nite-dimensional). 

. 	
The adjoint of a linear map (with respect to a Hermitian inner product). 

. 	
Existence of orthonormal bases in a Hermitian space (Proposition 13.11). 

. 	
Gram¨CSchmidt orthonormalization procedure. 

. 	
Linear isometries (unitary transformations). 

. 	
The unitary group, unitary matrices. 

. 	
The unitary group U(n). 

. 	
The special unitary group SU(n). 

. 	
QR-Decomposition for arbitrary complex matrices. 

. 	
The Hadamard inequality for complex matrices. 

. 	
The Hadamard inequality for Hermitian positive semide.nite matrices. 

. 	
Orthogonal projections and involutions; orthogonal re.ections. 

. 	
Dual norms. 

. 	
Nuclear norm (also called trace norm). 

. 	
Matrix completion. 


13.9 Problems 
Problem 13.1. Let (E,£¨.,.´¨) be a Hermitian space of .nite dimension. Prove that if 
f : E ¡ú E is a self-adjoint linear map (that is, f. = f), then£¨f(x),x´¨¡Ê R for all x ¡Ê E. 
Problem 13.2. Prove the polarization identities of Proposition 13.1. 
Problem 13.3. Let E be a real Euclidean space. Give an example of a nonzero linear map 
f : E ¡ú E such that£¨f(u),u´¨ = 0 for all u ¡Ê E. 
Problem 13.4. Prove Proposition 13.9. 
Problem 13.5. (1) Prove that every matrix in SU(2) is of the form 
a + ib c + id 
A = ,a 2 + b2 + c 2 + d2 =1, a,b,c,d ¡Ê R,
.c + id a . ib 
(2) Prove that the matrices 
10 i 0 010 i 
,, ,
01 0 .i .10 i 0 
all belong to SU(2) and are linearly independent over C. 
(3) Prove that the linear span of SU(2) over C is the complex vector space M2(C) of all complex 2 ¡Á 2 matrices. 
Problem 13.6. The purpose of this problem is to prove that the linear span of SU(n) over C is Mn(C) for all n ¡Ý 3. One way to prove this result is to adapt the method of Problem 11.12, so please review this problem. 
Every complex matrix A ¡Ê Mn(C) can be written as 
A + A. A . A. 
A =+ 
22 
where the .rst matrix is Hermitian and the second matrix is skew-Hermitian. Observe that if A =(zij) is a Hermitian matrix, that is A. = A, then zji = zij, so if zij = aij + ibij with aij,bij ¡Ê R, then aij = aji and bij = .bji. On the other hand, if A =(zij) is a skew-Hermitian matrix, that is A. = .A, then zji = .zij, so aij = .aji and bij = bji. 
The Hermitian and the skew-Hermitian matrices do not form complex vector spaces because they are not closed under multiplication by a complex number, but we can get around this problem by treating the real part and the complex part of these matrices separately and using multiplication by reals. 
13.9. PROBLEMS 
(1) Consider the matrices of the form 
.
. 
1 ... 
1 00 ¡¤¡¤¡¤ 0 i 01 ¡¤¡¤¡¤ 00 
Ri,j c 
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
i  0  ¡¤ ¡¤ ¡¤  0  0  
1  
...  

1 
).Ri,j
Prove that (Ri,j = I and det(Ri,j) = +1. Use the matrices Ri,j,Ri,j ¡Ê SU(n) and 
ccc c the matrices (Ri,j .(Ri,j).)/2 (from Problem 11.12) to form the real part of a skew-Hermitian . (Ri,j
matrix and the matrices (Rci,j c ).)/2 to form the imaginary part of a skew-Hermitian matrix. Deduce that the matrices in SU(n) span all skew-Hermitian matrices. 
(2) Consider matrices of the form 
Type 1 

.
. 
S1,2 = 
c 
........ 

0  .i  0  0  . . .  0  
i  0  0  0  . . .  0  
0  0  .1  0  . . .  0  

0  0  0  1  . . .  0  
. . .  . . .  . . .  . . .  ...  . . .  
0  0  0  0  . . .  1  

........ 

. 

Type 2 

.
. 
Si,i+1 c 
= 

.............. 

.1  
1  
.. .  
1  

0 .i 

i  0  
1  
...  
1  

.............. 

. 

Type 3 
1 
.
. 
.. . 1 
00 ¡¤¡¤¡¤ 0 .i 
0 .1 ¡¤¡¤¡¤ 00 
Si,j ..... 
= .. ... .
.
c .. .. 
................... 
00 ¡¤¡¤¡¤ 10 i 0 ¡¤¡¤¡¤ 00 1 .
.
. 1 
Prove that Si,j,Sci,j ¡Ê SU(n), and using diagonal matrices as in Problem 11.12, prove that the matrices Si,j can be used to form the real part of a Hermitian matrix and the matrices Sci,j can be used to form the imaginary part of a Hermitian matrix. 
(3) Use (1) and (2) to prove that the matrices in SU(n) span all Hermitian matrices. It follows that SU(n) spans Mn(C) for n ¡Ý 3. 
i 1 
A = . 
................... Problem13.7. Considerthecomplexmatrix .1 i 
Check that this matrix is symmetric but not Hermitian. Prove that 
det(¦ËI . A)= ¦Ë2 , 
and so the eigenvalues of A are 0, 0. 
Problem 13.8. Let (E,£¨.,.´¨) be a Hermitian space of .nite dimension and let f : E ¡ú E be a linear map. Prove that the following conditions are equivalent. 
(1) f . f. = f. . f (f is normal). 

(2)£¨
f(x),f(y)´¨ =£¨f.(x),f.(y)´¨ for all x, y ¡Ê E. 

(3) lf(x)l = lf.(x)l for all x ¡Ê E. 

(4) The map f can be diagonalized with respect to an orthonormal basis of eigenvectors. 

(5) There exist some linear maps g, h: E ¡ú E such that, g = g . ,£¨x, g(x)´¨¡Ý 0 for all x ¡Ê E, h.1 = h., and f = g . h = h . g. 

(6) There exist some linear map h: E ¡ú E such that h.1 = h. and f. = h . f. 

13.9. PROBLEMS 

(7) 
There is a polynomial P (with complex coe.cients) such that f. = P (f). 


Problem 13.9. Recall from Problem 12.7 that a complex n¡Án matrix H is upper Hessenberg if hjk =0 for all (j, k) such that j . k ¡Ý 0. Adapt the proof of Problem 12.7 to prove that given any complex n ¡Á n-matrix A, there are n . 2 ¡Ý 1 complex matrices H1,...,Hn.2, Householder matrices or the identity, such that 
B = Hn.2 ¡¤¡¤¡¤ H1AH1 ¡¤¡¤¡¤ Hn.2 
is upper Hessenberg. 
Problem 13.10. Prove that all y ¡Ê Cn , 
lylD 
1 = lyl¡Þ lylD = lyl
¡Þ 1 lylD 2 = lyl2 . 
Problem 13.11. The purpose of this problem is to complete each of the matrices A0,B0,C0 of Section 13.7 to a matrix A in such way that the nuclear norm lAlN is minimized. 
(1) Prove that the squares ¦Ò12 and ¦Ò22 of the singular values of 
12 
A = 
cd 
are the zeros of the equation 
¦Ë2 . (5 + c 2 + d2)¦Ë + (2c . d)2 =0. 
(2) Using the fact that 
 
lAlN = ¦Ò1 + ¦Ò2 =¦Ò12 + ¦Ò22 +2¦Ò1¦Ò2, 
prove that lAl2 =5+ c 2 + d2 +2|2c . d|.
N 
Consider the cases where 2c . d ¡Ý 0 and 2c . d ¡Ü 0, and show that in both cases we must have c = .2d, and that the minimum of f(c, d) = 5+ c2 + d2 +2|2c . d| is achieved by c = d = 0. Conclude that the matrix A completing A0 that minimizes lAlN is 
1  2  
A =  .  
0  0  

(3) Prove that the squares ¦Ò12 and ¦Ò22 of the singular values of 1 b 
A = 
c 4 

are the zeros of the equation 
¦Ë2 . (17 + b2 + c 2)¦Ë + (4 . bc)2 =0. 
(4) Prove that lAl2 = 17+ b2 + c 2 +2|4 . bc|.
N 
Consider the cases where 4 . bc ¡Ý 0 and 4 . bc ¡Ü 0, and show that in both cases we must have b2 = c2 . Then show that the minimum of f(c, d) = 17+ b2 + c2 +2|4 . bc| is achieved by b = c with .2 ¡Ü b ¡Ü 2. Conclude that the matrices A completing B0 that minimize lAl
N  
are given by  
A =  1 b  b 4  ,  .2 ¡Ü b ¡Ü 2.  

(5) Prove that the squares ¦Ò12 and ¦Ò22 of the singular values of 
12 
A = 
3 d 
are the zeros of the equation 
¦Ë2 . (14 + d2)¦Ë + (6 . d)2 =0 
(6) Prove that lAl2 = 14+ d2 +2|6 . d|.
N 
Consider the cases where 6 . d ¡Ý 0 and 6 . d ¡Ü 0, and show that the minimum of f(c, d)= 14 + d2 +2|6 . d| is achieved by d = 1. Conclude that the the matrix A completing C0 that minimizes lAlN is given by 
12 
A = . 
31 
Problem 13.12. Prove Theorem 13.32 when E is a .nite-dimensional Hermitian space. 


