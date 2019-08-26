Chapter 16 
Spectral Theorems in Euclidean and Hermitian Spaces 
16.1 	Introduction 
The goal of this chapter is to show that there are nice normal forms for symmetric matrices, skew-symmetric matrices, orthogonal matrices, and normal matrices. The spectral theorem for symmetric matrices states that symmetric matrices have real eigenvalues and that they can be diagonalized over an orthonormal basis. The spectral theorem for Hermitian matrices states that Hermitian matrices also have real eigenvalues and that they can be diagonalized over a complex orthonormal basis. Normal real matrices can be block diagonalized over an orthonormal basis with blocks having size at most two and there are re.nements of this normal form for skew-symmetric and orthogonal matrices. 
The spectral result for real symmetric matrices can be used to prove two characterizations of the eigenvalues of a symmetric matrix in terms of the Rayleigh ratio. The .rst charac-terization is the Rayleigh¨CRitz theorem and the second one is the Courant¨CFischer theorem. Both results are used in optimization theory and to obtain results about perturbing the eigenvalues of a symmetric matrix. 
In this chapter all vector spaces are .nite-dimensional real or complex vector spaces. 

16.2 	Normal Linear Maps: Eigenvalues and Eigenvec-tors 
We begin by studying normal maps, to understand the structure of their eigenvalues and eigenvectors. This section and the next three were inspired by Lang [41], Artin [3], Mac Lane and Birkho. [46], Berger [5], and Bertin [7]. 
De.nition 16.1. Given a Euclidean or Hermitian space E, a linear map f : E ¡ú E is normal if 
f . f . = f . . f. 
543 

A linear map f : E ¡ú E is self-adjoint if f = f. , skew-self-adjoint if f = .f., and orthogonal if f . f. = f. . f = id. 
Obviously, a self-adjoint, skew-self-adjoint, or orthogonal linear map is a normal linear map. Our .rst goal is to show that for every normal linear map f : E ¡ú E, there is an orthonormal basis (w.r.t. (.,.£©) such that the matrix of f over this basis has an especially nice form: it is a block diagonal matrix in which the blocks are either one-dimensional matrices (i.e., single entries) or two-dimensional matrices of the form 
¦Ë¦Ì ..¦Ì¦Ë 
This normal form can be further re.ned if f is self-adjoint, skew-self-adjoint, or orthog-onal. As a .rst step we show that f and f. have the same kernel when f is normal. 
Proposition 16.1. Given a Euclidean space E, if f : E ¡ú E is a normal linear map, then Ker f = Ker f. . 
Proof. First let us prove that 
(f(u),f(v)£© = (f . (u),f . (v)£© 
for all u, v ¡Ê E. Since f. is the adjoint of f and f . f. = f. . f, we have 
(f(u),f(u)£© = (u, (f . . f)(u)£©, = (u, (f . f . )(u)£©, = (f . (u),f . (u)£©. 
Since (.,.£© is positive de.nite, 
(f(u),f(u)£© = 0 i. f(u)=0, 
(f . (u),f . (u)£© = 0 i. f . (u)=0, 

and since (f(u),f(u)£© = (f . (u),f . (u)£©, 
we have f(u)=0 i. f . (u)=0. 
Consequently, Ker f = Ker f. . 
Assuming again that E is a Hermitian space, observe that Proposition 16.1 also holds. We deduce the following corollary. 
Proposition 16.2. Given a Hermitian space E, for any normal linear map f : E ¡ú E, we have Ker (f) ¡É Im(f) = (0). 
16.2. NORMAL LINEAR MAPS: EIGENVALUES AND EIGENVECTORS 
Proof. Assume v ¡Ê Ker (f) ¡É Im(f) = (0), which means that v = f(u) for some u ¡Ê E, and f(v) = 0. By Proposition 16.1, Ker (f) = Ker(f.), so f(v) = 0 implies that f.(v) = 0. Consequently, 
0= (f . (v),u£© = (v, f(u)£© = (v, v£©, 
and thus, v = 0. 
We also have the following crucial proposition relating the eigenvalues of f and f. . 
Proposition 16.3. Given a Hermitian space E, for any normal linear map f : E ¡ú E,a vector u is an eigenvector of f for the eigenvalue ¦Ë (in C) i. u is an eigenvector of f. for the eigenvalue ¦Ë. 
Proof. First it is immediately veri.ed that the adjoint of f . ¦Ë id is f. . ¦Ë id. Furthermore, f . ¦Ë id is normal. Indeed, 
(f . ¦Ë id) . (f . ¦Ë id) . =(f . ¦Ë id) . (f . . ¦Ë id), 
= f . f . . ¦Ëf . ¦Ëf . + ¦Ë¦Ë id, 
= f . . f . ¦Ëf . . ¦Ëf + ¦Ë¦Ë id, 
=(f . . ¦Ë id) . (f . ¦Ë id), 
=(f . ¦Ë id) . . (f . ¦Ë id). 
Applying Proposition 16.1 to f . ¦Ë id, for every nonnull vector u, we see that 
(f . ¦Ë id)(u)=0 i. (f . . ¦Ë id)(u)=0, 
which is exactly the statement of the proposition. 
The next proposition shows a very important property of normal linear maps: eigenvec-tors corresponding to distinct eigenvalues are orthogonal. 
Proposition 16.4. Given a Hermitian space E, for any normal linear map f : E ¡ú E, if u and v are eigenvectors of f associated with the eigenvalues ¦Ë and ¦Ì (in C) where ¦Ë 
= ¦Ì, then (u, v£© =0. 
Proof. Let us compute (f(u),v£© in two di.erent ways. Since v is an eigenvector of f for ¦Ì, by Proposition 16.3, v is also an eigenvector of f. for ¦Ì, and we have 
(f(u),v£© = (¦Ëu, v£© = ¦Ë(u, v£©, 
and 
(f(u),v£© = (u, f . (v)£© = (u, ¦Ìv£© = ¦Ì(u, v£©, 
where the last identity holds because of the semilinearity in the second argument. Thus ¦Ë(u, v£© = ¦Ì(u, v£©, that is, (¦Ë . ¦Ì)(u, v£© =0, which implies that (u, v£© = 0, since ¦Ë = ¦Ì. We can show easily that the eigenvalues of a self-adjoint linear map are real. 
Proposition 16.5. Given a Hermitian space E, all the eigenvalues of any self-adjoint linear 
map f : E ¡ú E are real. Proof. Let z (in C) be an eigenvalue of f and let u be an eigenvector for z. We compute (f(u),u£© in two di.erent ways. We have 
(f(u),u£© = (zu, u£© = z(u, u£©, and since f = f., we also have (f(u),u£© = (u, f . (u)£© = (u, f(u)£© = (u, zu£© = z(u, u£©. Thus, z(u, u£© = z(u, u£©, which implies that z = z, since u = 0, and z is indeed real. There is also a version of Proposition 16.5 for a (real) Euclidean space E and a self-adjoint map f : E ¡ú E since every real vector space E can be embedded into a complex vector space EC, and every linear map f : E ¡ú E can be extended to a linear map fC: EC ¡ú EC. De.nition 16.2. Given a real vector space E, let EC be the structure E ¡Á E under the addition operation (u1,u2)+(v1,v2)=(u1 + v1,u2 + v2), and let multiplication by a complex scalar z = x + iy be de.ned such that (x + iy) ¡¤ (u, v)=(xu . yv, yu + xv). The space EC is called the complexi.cation of E. It is easily shown that the structure EC is a complex vector space. It is also immediate that (0,v)= i(v, 0), and thus, identifying E with the subspace of EC consisting of all vectors of the form (u, 0), we can write 
(u, v)= u + iv. 


16.2. NORMAL LINEAR MAPS: EIGENVALUES AND EIGENVECTORS 
Observe that if (e1,...,en) is a basis of E (a real vector space), then (e1,...,en) is also a basis of EC (recall that ei is an abbreviation for (ei, 0)). 
A linear map f : E ¡ú E is extended to the linear map fC: EC ¡ú EC de.ned such that 
fC(u + iv)= f(u)+ if(v). 
For any basis (e1,...,en) of E, the matrix M(f) representing f over (e1,...,en) is iden-tical to the matrix M(fC) representing fC over (e1,...,en), where we view (e1,...,en) as a basis of EC. As a consequence, det(zI . M(f)) = det(zI . M(fC)), which means that f and fC have the same characteristic polynomial (which has real coe.cients). We know that every polynomial of degree n with real (or complex) coe.cients always has n complex roots (counted with their multiplicity), and the roots of det(zI . M(fC)) that are real (if any) are the eigenvalues of f. 
Next we need to extend the inner product on E to an inner product on EC. 
The inner product (.,.£© on a Euclidean space E is extended to the Hermitian positive de.nite form (.,.£©C on EC as follows: 
(u1 + iv1,u2 + iv2£©C = (u1,u2£© + (v1,v2£© + i((v1,u2£©.( u1,v2£©). 
It is easily veri.ed that (.,.£©C is indeed a Hermitian form that is positive de.nite, and it is clear that (.,.£©C agrees with (.,.£© on real vectors. Then given any linear map 
f : E ¡ú E, it is easily veri.ed that the map fC . de.ned such that 
fC. (u + iv)= f . (u)+ if . (v) 
for all u, v ¡Ê E is the adjoint of fC w.r.t. (.,.£©C. 
Proposition 16.6. Given a Euclidean space E, if f : E ¡ú E is any self-adjoint linear map, then every eigenvalue ¦Ë of fC is real and is actually an eigenvalue of f (which means that there is some real eigenvector u ¡Ê E such that f(u)= ¦Ëu). Therefore, all the eigenvalues of f are real. 
Proof. Let EC be the complexi.cation of E, (.,.£©C the complexi.cation of the inner product (.,.£© on E, and fC: EC ¡ú EC the complexi.cation of f : E ¡ú E. By de.nition of fC and (.,.£©C, if f is self-adjoint, we have 
(fC(u1 + iv1),u2 + iv2£©C = (f(u1)+ if(v1),u2 + iv2£©C = (f(u1),u2£© + (f(v1),v2£© 
+ i((u2,f(v1)£©.( f(u1),v2£©) = (u1,f(u2)£© + (v1,f(v2)£© 

+ i((f(u2),v1£©.( u1,f(v2)£©) = (u1 + iv1,f(u2)+ if(v2)£©C = (u1 + iv1,fC(u2 + iv2)£©C, 


which shows that fC is also self-adjoint with respect to (.,.£©C. 
As we pointed out earlier, f and fC have the same characteristic polynomial det(zI.fC)= det(zI . f), which is a polynomial with real coe.cients. Proposition 16.5 shows that the zeros of det(zI . fC) = det(zI . f) are all real, and for each real zero ¦Ë of det(zI . f), the linear map ¦Ëid . f is singular, which means that there is some nonzero u ¡Ê E such that f(u)= ¦Ëu. Therefore, all the eigenvalues of f are real. 
Proposition 16.7. Given a Hermitian space E, for any linear map f : E ¡ú E, if f is skew-self-adjoint, then f has eigenvalues that are pure imaginary or zero, and if f is unitary, then f has eigenvalues of absolute value 1. 
Proof. If f is skew-self-adjoint, f. = .f, and then by the de.nition of the adjoint map, for any eigenvalue ¦Ë and any eigenvector u associated with ¦Ë, we have 
¦Ë(u, u£© = (¦Ëu, u£© = (f(u),u£© = (u, f . (u)£© = (u, .f(u)£© = .(u, ¦Ëu£© = .¦Ë(u, u£©, 
and since u = 0 and (.,.£© is positive de.nite, (u, u£© = 0, so 
¦Ë = .¦Ë, 
which shows that ¦Ë = ir for some r ¡Ê R. 
If f is unitary, then f is an isometry, so for any eigenvalue ¦Ë and any eigenvector u associated with ¦Ë, we have 
|¦Ë|2(u, u£© = ¦Ë¦Ë(u, u£© = (¦Ëu, ¦Ëu£© = (f(u),f(u)£© = (u, u£©, 
and since u = 0, we obtain |¦Ë|2 = 1, which implies 
|¦Ë| =1. 


16.3 Spectral Theorem for Normal Linear Maps 
Given a Euclidean space E, our next step is to show that for every linear map f : E ¡ú E there is some subspace W of dimension 1 or 2 such that f(W ) . W . When dim(W ) = 1, the subspace W is actually an eigenspace for some real eigenvalue of f. Furthermore, when f is normal, there is a subspace W of dimension 1 or 2 such that f(W ) . W and f.(W ) . W . The di.culty is that the eigenvalues of f are not necessarily real. One way to get around this problem is to complexify both the vector space E and the inner product (.,.£© as we did in Section 16.2. 
Given any subspace W of a Euclidean space E, recall that the orthogonal complement W ¡Í of W is the subspace de.ned such that 
W ¡Í = {u ¡Ê E |(u, w£© =0, for all w ¡Ê W }. 
16.3. SPECTRAL THEOREM FOR NORMAL LINEAR MAPS 
Recall from Proposition 11.11 that E = W ¨’ W ¡Í (this can be easily shown, for example, by constructing an orthonormal basis of E using the Gram¨CSchmidt orthonormalization procedure). The same result also holds for Hermitian spaces; see Proposition 13.13. 
As a warm up for the proof of Theorem 16.12, let us prove that every self-adjoint map on a Euclidean space can be diagonalized with respect to an orthonormal basis of eigenvectors. 
Theorem 16.8. (Spectral theorem for self-adjoint linear maps on a Euclidean space) Given a Euclidean space E of dimension n, for every self-adjoint linear map f : E ¡ú E, there is an orthonormal basis (e1,...,en) of eigenvectors of f such that the matrix of f w.r.t. this 
basis is a diagonal matrix

.
. 
.... 

¦Ë1 ... ¦Ë2 ... 
,
.. .
.
.. ..
.
.. . 
....
... ¦Ën 
with ¦Ëi ¡Ê R. 
Proof. We proceed by induction on the dimension n of E as follows. If n = 1, the result is trivial. Assume now that n ¡Ý 2. From Proposition 16.6, all the eigenvalues of f are real, so pick some eigenvalue ¦Ë ¡Ê R, and let w be some eigenvector for ¦Ë. By dividing w by its norm, we may assume that w is a unit vector. Let W be the subspace of dimension 1 spanned by w. Clearly, f(W ) . W . We claim that f(W ¡Í) . W ¡Í, where W ¡Í is the orthogonal complement of W . 
Indeed, for any v ¡Ê W ¡Í, that is, if (v, w£© = 0, because f is self-adjoint and f(w)= ¦Ëw, we have 
(f(v),w£© = (v, f(w)£© 
= (v, ¦Ëw£© 
= ¦Ë(v, w£© =0 
since (v, w£© = 0. Therefore, f(W ¡Í) . W ¡Í . 
Clearly, the restriction of f to W ¡Í is self-adjoint, and we conclude by applying the induction hypothesis to W ¡Í (whose dimension is n . 1). 
We now come back to normal linear maps. One of the key points in the proof of Theorem 
16.8 is that we found a subspace W with the property that f(W ) . W implies that f(W ¡Í) . W ¡Í . In general, this does not happen, but normal maps satisfy a stronger property which ensures that such a subspace exists. 
The following proposition provides a condition that will allow us to show that a nor-mal linear map can be diagonalized. It actually holds for any linear map. We found the inspiration for this proposition in Berger [5]. 
Proposition 16.9. Given a Hermitian space E, for any linear map f : E ¡ú E and any 
JW 
subspace W of E, if f(W ) . W , then f. W ¡Í . W ¡Í . Consequently, if f(W ) . W and
JW JW 
f.(W ) . W , then fW ¡Í . W ¡Í and f. W ¡Í . W ¡Í . 
Proof. If u ¡Ê W ¡Í, then (w, u£© = 0 for all w ¡Ê W. 
However, (f(w),u£© = (w, f . (u)£©, 
and f(W ) . W implies that f(w) ¡Ê W . Since u ¡Ê W ¡Í, we get 
0= (f(w),u£© = (w, f . (u)£©, 
which shows that (w, f.(u)£© = 0 for all w ¡Ê W , that is, f.(u) ¡Ê W ¡Í . Therefore, we have f.(W ¡Í) . W ¡Í . 
JW 
We just proved that if f(W ) . W , then f. W ¡Í . W ¡Í . If we also have f.(W ) . W , then by applying the above fact to f. , we get f..(W ¡Í) . W ¡Í, and since f.. = f, this is just f(W ¡Í) . W ¡Í, which proves the second statement of the proposition. 
It is clear that the above proposition also holds for Euclidean spaces. 
Although we are ready to prove that for every normal linear map f (over a Hermitian space) there is an orthonormal basis of eigenvectors (see Theorem 16.13 below), we now return to real Euclidean spaces. 
Proposition 16.10. If f : E ¡ú E is a linear map and w = u + iv is an eigenvector of fC: EC ¡ú EC for the eigenvalue z = ¦Ë + i¦Ì, where u, v ¡Ê E and ¦Ë, ¦Ì ¡Ê R, then 
f(u)= ¦Ëu . ¦Ìv and f(v)= ¦Ìu + ¦Ëv. (.) 
As a consequence, 
fC(u . iv)= f(u) . if(v)=(¦Ë . i¦Ì)(u . iv), 
which shows that w = u . iv is an eigenvector of fC for z = ¦Ë . i¦Ì. 
Proof. Since fC(u + iv)= f(u)+ if(v) 
and fC(u + iv)=(¦Ë + i¦Ì)(u + iv)= ¦Ëu . ¦Ìv + i(¦Ìu + ¦Ëv), 
we have f(u)= ¦Ëu . ¦Ìv and f(v)= ¦Ìu + ¦Ëv. 
Using this fact, we can prove the following proposition. 


16.3. SPECTRAL THEOREM FOR NORMAL LINEAR MAPS 
Proposition 16.11. Given a Euclidean space E, for any normal linear map f : E ¡ú E, if w = u + iv is an eigenvector of fC associated with the eigenvalue z = ¦Ë + i¦Ì (where u, v ¡Ê E and ¦Ë, ¦Ì ¡Ê R), if ¦Ì =0 (i.e., z is not real) then (u, v£© =0 and (u, u£© = (v, v£©, which implies that u and v are linearly independent, and if W is the subspace spanned by u and v, then f(W )= W and f.(W )= W . Furthermore, with respect to the (orthogonal) basis (u, v), the restriction of f to W has the matrix 
¦Ë¦Ì 
.
.¦Ì¦Ë 
If ¦Ì =0, then ¦Ë is a real eigenvalue of f, and either u or v is an eigenvector of f for ¦Ë. If W is the subspace spanned by u if u =0, or spanned by v =0 if u =0, then f(W ) . W and f.(W ) . W . 
Proof. Since w = u + iv is an eigenvector of fC, by de.nition it is nonnull, and either u =0 or v = 0. Proposition 16.10 implies that u . iv is an eigenvector of fC for ¦Ë . i¦Ì. It is easy to check that fC is normal. However, if ¦Ì = 0, then ¦Ë + i¦Ì = ¦Ë . i¦Ì, and from Proposition 
16.4, the vectors u + iv and u . iv are orthogonal w.r.t. (.,.£©C, that is, 
(u + iv, u . iv£©C = (u, u£©.( v, v£© +2i(u, v£© =0. 
Thus we get (u, v£© = 0 and (u, u£© = (v, v£©, and since u = 0 or v = 0, u and v are linearly independent. Since 
f(u)= ¦Ëu . ¦Ìv and f(v)= ¦Ìu + ¦Ëv 
and since by Proposition 16.3 u + iv is an eigenvector of fC . for ¦Ë . i¦Ì, we have 
f . (u)= ¦Ëu + ¦Ìv and f . (v)= .¦Ìu + ¦Ëv, 
and thus f(W )= W and f.(W )= W , where W is the subspace spanned by u and v. When ¦Ì = 0, we have f(u)= ¦Ëu and f(v)= ¦Ëv, 
and since u = 0 or v = 0, either u or v is an eigenvector of f for ¦Ë. If W is the subspace spanned by u if u = 0, or spanned by v if u = 0, it is obvious that f(W ) . W and f.(W ) . W . Note that ¦Ë = 0 is possible, and this is why . cannot be replaced by =. 
The beginning of the proof of Proposition 16.11 actually shows that for every linear map 
f : E ¡ú E there is some subspace W such that f(W ) . W , where W has dimension 1 or 
2. In general, it doesn¡¯t seem possible to prove that W ¡Í is invariant under f. However, this 
happens when f is normal. We can .nally prove our .rst main theorem. 
Theorem 16.12. (Main spectral theorem) Given a Euclidean space E of dimension n, for every normal linear map f : E ¡ú E, there is an orthonormal basis (e1,...,en) such that the matrix of f w.r.t. this basis is a block diagonal matrix of the form 
.
. 
.... 

A1 ... A2 ... 
... 
... Ap 
... 
... 
... 
.... 

such that each block Aj is either a one-dimensional matrix (i.e., a real scalar) or a two-
dimensional matrix of the form  
Aj  =  ¦Ëj ¦Ìj  .¦Ìj ¦Ëj  ,  
where ¦Ëj, ¦Ìj ¡Ê R, with ¦Ìj > 0.  

Proof. We proceed by induction on the dimension n of E as follows. If n = 1, the result is trivial. Assume now that n ¡Ý 2. First, since C is algebraically closed (i.e., every polynomial has a root in C), the linear map fC: EC ¡ú EC has some eigenvalue z = ¦Ë + i¦Ì (where ¦Ë, ¦Ì ¡Ê R). Let w = u + iv be some eigenvector of fC for ¦Ë + i¦Ì (where u, v ¡Ê E). We can now apply Proposition 16.11. 
If ¦Ì = 0, then either u or v is an eigenvector of f for ¦Ë ¡Ê R. Let W be the subspace of dimension 1 spanned by e1 = u/½Ðu½Ð if u = 0, or by e1 = ½Ð otherwise. It is obvious 
v/½Ðvthat f(W ) . W and f.(W ) . WW. The orthogonal W ¡Í of W has dimension n . 1, and by Proposition 16.9, we have fW ¡Í . W ¡Í . But the restriction of f to W ¡Í is also normal, 
J 
and we conclude by applying the induction hypothesis to W ¡Í . 
If ¦Ì = 0, then (u, v£© = 0 and (u, u£© = (v, v£©, and if W is the subspace spanned by u/½Ðu½Ð and v/½Ðv½Ð, then f(W )= W and f.(W )= W . We also know that the restriction of f to W has the matrix 
¦Ë¦Ì .¦Ì¦Ë 
with respect to the basis (u/½Ðu½Ð, v/½Ðv½Ð). If ¦Ì< 0, we let ¦Ë1 = ¦Ë, ¦Ì1 = .¦Ì, e1 = u/½Ðu½Ð, and e2 = v/½Ðv½Ð. If ¦Ì> 0, we let ¦Ë1 = ¦Ë, ¦Ì1 = ¦Ì, e1 = v/½Ðv½Ð, and e2 = u/½Ðu½Ð. In all cases, it is easily veri.ed that the matrix of the restriction of f to W w.r.t. the orthonormal basis (e1,e2) is 
¦Ë1 .¦Ì1
A1 = , 
¦Ì1 ¦Ë1 
W
J
where ¦Ë1,¦Ì1 ¡Ê R, with ¦Ì1 > 0. However, W ¡Í has dimension n . 2, and by Proposition 16.9, fW ¡Í . W ¡Í . Since the restriction of f to W ¡Í is also normal, we conclude by applying the induction hypothesis to W ¡Í . 
After this relatively hard work, we can easily obtain some nice normal forms for the matrices of self-adjoint, skew-self-adjoint, and orthogonal linear maps. However, for the sake of completeness (and since we have all the tools to so do), we go back to the case of a Hermitian space and show that normal linear maps can be diagonalized with respect to an orthonormal basis. The proof is a slight generalization of the proof of Theorem 16.6. 

16.4. SELF-ADJOINT AND OTHER SPECIAL LINEAR MAPS 
Theorem 16.13. (Spectral theorem for normal linear maps on a Hermitian space) Given a Hermitian space E of dimension n, for every normal linear map f : E ¡ú E there is an orthonormal basis (e1,...,en) of eigenvectors of f such that the matrix of f w.r.t. this basis 
is a diagonal matrix

.
. 
.... 

¦Ë1 ... ¦Ë2 ... 
.. .
.
.. ..
.
.. . 
... ¦Ën 
....

, 

where ¦Ëj ¡Ê C. 
Proof. We proceed by induction on the dimension n of E as follows. If n = 1, the result is trivial. Assume now that n ¡Ý 2. Since C is algebraically closed (i.e., every polynomial has a root in C), the linear map f : E ¡ú E has some eigenvalue ¦Ë ¡Ê C, and let w be some unit eigenvector for ¦Ë. Let W be the subspace of dimension 1 spanned by w. Clearly, f(W ) . W . By Proposition 16.3, w is an eigenvector of f. for ¦Ë, and thus f.(W ) . W . By Proposition 16.9, we also have f(W ¡Í) . W ¡Í . The restriction of f to W ¡Í is still normal, and we conclude by applying the induction hypothesis to W ¡Í (whose dimension is n . 1). 
Theorem 16.13 implies that (complex) self-adjoint, skew-self-adjoint, and orthogonal lin-ear maps can be diagonalized with respect to an orthonormal basis of eigenvectors. In this latter case, though, an orthogonal map is called a unitary map. Proposition 16.5 also shows that the eigenvalues of a self-adjoint linear map are real, and Proposition 16.7 shows that the eigenvalues of a skew self-adjoint map are pure imaginary or zero, and that the eigenvalues of a unitary map have absolute value 1. 
Remark: There is a converse to Theorem 16.13, namely, if there is an orthonormal basis (e1,...,en) of eigenvectors of f, then f is normal. We leave the easy proof as an exercise. 
In the next section we specialize Theorem 16.12 to self-adjoint, skew-self-adjoint, and orthogonal linear maps. Due to the additional structure, we obtain more precise normal forms. 


16.4 	Self-Adjoint, Skew-Self-Adjoint, and Orthogonal Linear Maps 
We begin with self-adjoint maps. 
Theorem 16.14. Given a Euclidean space E of dimension n, for every self-adjoint linear map f : E ¡ú E, there is an orthonormal basis (e1,...,en) of eigenvectors of f such that the matrix of f w.r.t. this basis is a diagonal matrix 
.
. 
.... 

¦Ë1 ... ¦Ë2 ... 
,
..
. 
... ¦Ën 
....
... 
... 
... 

where ¦Ëi ¡Ê R. 
Proof. We already proved this; see Theorem 16.8. However, it is instructive to give a more direct method not involving the complexi.cation of (.,.£© and Proposition 16.5. Since C is algebraically closed, fC has some eigenvalue ¦Ë + i¦Ì, and let u + iv be some eigenvector of fC for ¦Ë+i¦Ì, where ¦Ë, ¦Ì ¡Ê R and u, v ¡Ê E. We saw in the proof of Proposition 
16.10 that f(u)= ¦Ëu . ¦Ìv and f(v)= ¦Ìu + ¦Ëv. 
Since f = f. , (f(u),v£© = (u, f(v)£© 
for all u, v ¡Ê E. Applying this to 
f(u)= ¦Ëu . ¦Ìv and f(v)= ¦Ìu + ¦Ëv, 
we get (f(u),v£© = (¦Ëu . ¦Ìv, v£© = ¦Ë(u, v£©. ¦Ì(v, v£© 
and (u, f(v)£© = (u, ¦Ìu + ¦Ëv£© = ¦Ì(u, u£© + ¦Ë(u, v£©, 
and thus we get ¦Ë(u, v£©. ¦Ì(v, v£© = ¦Ì(u, u£© + ¦Ë(u, v£©, 
that is, ¦Ì((u, u£© + (v, v£©)=0, 
which implies ¦Ì = 0, since either u =0 or v = 0. Therefore, ¦Ë is a real eigenvalue of f. Now going back to the proof of Theorem 16.12, only the case where ¦Ì = 0 applies, and the induction shows that all the blocks are one-dimensional. 
Theorem 16.14 implies that if ¦Ë1,...,¦Ëp are the distinct real eigenvalues of f, and Ei is the eigenspace associated with ¦Ëi, then 
E = E1 ¨’¡¤ ¡¤¡¤¨’ Ep, 
where Ei and Ej are orthogonal for all i = j. 
16.4. SELF-ADJOINT AND OTHER SPECIAL LINEAR MAPS 
Remark: Another way to prove that a self-adjoint map has a real eigenvalue is to use a little bit of calculus. We learned such a proof from Herman Gluck. The idea is to consider the real-valued function ¦µ: E ¡ú R de.ned such that 
¦µ(u)= (f(u),u£© 
for every u ¡Ê E. This function is C¡Þ , and if we represent f by a matrix A over some orthonormal basis, it is easy to compute the gradient vector
·É¦µ(X)= .¦µ(X),..., .¦µ(X)
.x1 .xn of ¦µ at X. Indeed, we .nd that·É¦µ(X)=(A + AT)X, where X is a column vector of size n. But since f is self-adjoint, A = AT, and thus·É¦µ(X)=2AX. The next step is to .nd the maximum of the function ¦µ on the sphere Sn.1 = {(x1,...,xn) ¡Ê Rn | x 21 + ¡¤¡¤¡¤ + x 2 =1}.
n 
Since Sn.1 is compact and ¦µ is continuous, and in fact C¡Þ, ¦µ takes a maximum at some X on Sn.1 . But then it is well known that at an extremum X of ¦µ we must have 
d¦µX (Y )=(·É¦µ(X),Y£© =0 for all tangent vectors Y to Sn.1 at X, and so·É¦µ(X) is orthogonal to the tangent plane at X, which means that
·É¦µ(X)= ¦ËX for some ¦Ë ¡Ê R. Since·É¦µ(X)=2AX, we get 2AX = ¦ËX, and thus ¦Ë/2 is a real eigenvalue of A (i.e., of f). Next we consider skew-self-adjoint maps. 
Theorem 16.15. Given a Euclidean space E of dimension n, for every skew-self-adjoint linear map f : E ¡ú E there is an orthonormal basis (e1,...,en) such that the matrix of f 
w.r.t. this basis is a block diagonal matrix of the form 
.
. 
.... 

A1 ... A2 ... 
.. .
.
....
.
.. . 
... Ap 
.... 

such that each block Aj is either 0 or a two-dimensional matrix of the form 0 .¦Ìj
Aj = , 
¦Ìj 0 
where ¦Ìj ¡Ê R, with ¦Ìj > 0. In particular, the eigenvalues of fC are pure imaginary of the form ¡Ài¦Ìj or 0. Proof. The case where n = 1 is trivial. As in the proof of Theorem 16.12, fC has some 
eigenvalue z = ¦Ë + i¦Ì, where ¦Ë, ¦Ì ¡Ê R. We claim that ¦Ë = 0. First we show that (f(w),w£© =0 for all w ¡Ê E. Indeed, since f = .f. , we get (f(w),w£© = (w, f . (w)£© = (w, .f(w)£© = .(w, f(w)£© = .(f(w),w£©, since (.,.£© is symmetric. This implies that (f(w),w£© =0. Applying this to u and v and using the fact that f(u)= ¦Ëu . ¦Ìv and f(v)= ¦Ìu + ¦Ëv, we get 0= (f(u),u£© = (¦Ëu . ¦Ìv, u£© = ¦Ë(u, u£©. ¦Ì(u, v£© and 0= (f(v),v£© = (¦Ìu + ¦Ëv, v£© = ¦Ì(u, v£© + ¦Ë(v, v£©, from which, by addition, we get ¦Ë((v, v£© + (v, v£©)=0. Since u =0 or v = 0, we have ¦Ë = 0. Then going back to the proof of Theorem 16.12, unless ¦Ì = 0, the case where u and v are orthogonal and span a subspace of dimension 2 applies, and the induction shows that all the blocks are two-dimensional or reduced to 0. 
Remark: One will note that if f is skew-self-adjoint, then ifC is self-adjoint w.r.t. (.,.£©C. By Proposition 16.5, the map ifC has real eigenvalues, which implies that the eigenvalues of fC are pure imaginary or 0. 
Finally we consider orthogonal linear maps. 

16.4. SELF-ADJOINT AND OTHER SPECIAL LINEAR MAPS 
Theorem 16.16. Given a Euclidean space E of dimension n, for every orthogonal linear map f : E ¡ú E there is an orthonormal basis (e1,...,en) such that the matrix of f w.r.t. 
such that each block Aj is either 1, .1, or a two-dimensional matrix of the form cos ¦Èj . sin ¦Èj
Aj = 
sin ¦Èj cos ¦Èj 
where 0 <¦Èj <¦Ð. In particular, the eigenvalues of fC are of the form cos ¦Èj ¡À i sin ¦Èj, 1, or .1. 
Proof. The case where n = 1 is trivial. It is immediately veri.ed that f . f. = f. . f = id implies that fC . fC . = fC . . fC = id, so the map fC is unitary. By Proposition 16.7, the eigenvalues of fC have absolute value 1. As a consequence, the eigenvalues of fC are of the form cos ¦È ¡À i sin ¦È, 1, or .1. The theorem then follows immediately from Theorem 16.12, where the condition ¦Ì> 0 implies that sin ¦Èj > 0, and thus, 0 <¦Èj <¦Ð. 
It is obvious that we can reorder the orthonormal basis of eigenvectors given by Theorem 
where each block Aj is a two-dimensional rotation matrix Aj = ¡ÀI2 of the form cos ¦Èj . sin ¦Èj
Aj = 
sin ¦Èj cos ¦Èj 
with 0 <¦Èj <¦Ð. 
The linear map f has an eigenspace E(1,f) = Ker(f . id) of dimension p for the eigen-value 1, and an eigenspace E(.1,f) = Ker(f + id) of dimension q for the eigenvalue .1. If det(f)=+1 (f is a rotation), the dimension q of E(.1,f) must be even, and the entries in .Iq can be paired to form two-dimensional blocks, if we wish. In this case, every rotation in SO(n) has a matrix of the form 
16.16, so that the matrix of f w.r.t. this basis is a block diagonal matrix of the form

. .  
A1  . . .  
. . . ......  . ..  . . .  . . . ......  
. . .  Ar  
.Iq  
. . .  Ip  

. 

.... 

this basis is a block diagonal matrix of the form
A1 ... A2 ... 
.. .
.
....
.
.. . 
... Ap 
. 
.... 

.
. 
.... 

A1 ... 
..
.
...
.
.. 
... Am ... In.2m 
.... 

where the .rst m blocks Aj are of the form 
cos ¦Èj . sin ¦Èj
Aj = 
sin ¦Èj cos ¦Èj 
with 0 <¦Èj ¡Ü ¦Ð. Theorem 16.16 can be used to prove a version of the Cartan¨CDieudonn¡äe theorem. 
Theorem 16.17. Let E be a Euclidean space of dimension n ¡Ý 2. For every isometry f ¡Ê O(E), if p = dim(E(1,f)) = dim(Ker (f . id)), then f is the composition of n . p re.ections, and n . p is minimal. 
Proof. From Theorem 16.16 there are r subspaces F1,...,Fr, each of dimension 2, such that 
E = E(1,f) ¨’ E(.1,f) ¨’ F1 ¨’¡¤ ¡¤¡¤¨’ Fr, 
and all the summands are pairwise orthogonal. Furthermore, the restriction ri of f to each 
H
Fi is a rotation ri = ¡Àid. Each 2D rotation ri can be written as the composition ri = si . si 
H
of two re.ections si and si about lines in Fi (forming an angle ¦Èi/2). We can extend si and 
H
si to hyperplane re.ections in E by making them the identity on Fi ¡Í . Then 
HH
s.¡¤ ¡¤¡¤. s
r . sr 1 . s1 
agrees with f on F1 ¨’ ¡¤ ¡¤¡¤ ¨’ Fr and is the identity on E(1,f) ¨’ E(.1,f). If E(.1,f) 
HH
has an orthonormal basis of eigenvectors (v1,...,vq), letting sbe the re.ection about the 
j 
hyperplane (vj)¡Í , it is clear that HH HH 
s .¡¤ ¡¤¡¤. s
q 1 
agrees with f on E(.1,f) and is the identity on E(1,f) ¨’ F1 ¨’¡¤ ¡¤¡¤¨’ Fr. But then 
HHHHH H
f = sq .¡¤ ¡¤¡¤. s1 . sr . sr .¡¤ ¡¤¡¤. s1 . s1, 
the composition of 2r + q = n . p re.ections. If f = st .¡¤ ¡¤¡¤. s1, 
for t re.ections si, it is clear that 
t
 
F = E(1,si) . E(1,f), 
i=1 
where E(1,si) is the hyperplane de.ning the re.ection si. By the Grassmann relation, if we intersect t ¡Ü n hyperplanes, the dimension of their intersection is at least n . t. Thus, n . t ¡Ü p, that is, t ¡Ý n . p, and n . p is the smallest number of re.ections composing f. 
As a corollary of Theorem 16.17, we obtain the following fact: If the dimension n of the Euclidean space E is odd, then every rotation f ¡Ê SO(E) admits 1 as an eigenvalue. 

16.5. NORMAL AND OTHER SPECIAL MATRICES 
Proof. The characteristic polynomial det(XI . f) of f has odd degree n and has real coef-.cients, so it must have some real root ¦Ë. Since f is an isometry, its n eigenvalues are of the form, +1, .1, and e¡Ài¦È, with 0 <¦È<¦Ð, so ¦Ë = ¡À1. Now the eigenvalues e¡Ài¦È appear in conjugate pairs, and since n is odd, the number of real eigenvalues of f is odd. This implies that +1 is an eigenvalue of f, since otherwise .1 would be the only real eigenvalue of f, and since its multiplicity is odd, we would have det(f)= .1, contradicting the fact that f is a rotation. 
When n = 3, we obtain the result due to Euler which says that every 3D rotation R has an invariant axis D, and that restricted to the plane orthogonal to D, it is a 2D rotation. Furthermore, if (a, b, c) is a unit vector de.ning the axis D of the rotation R and if the angle of the rotation is ¦È, if B is the skew-symmetric matrix 
.  .  
0  .c  b  
B = . c  0  .a. ,  
.b  a  0  

then the Rodigues formula (Proposition 11.15) states that 
R = I + sin ¦ÈB + (1 . cos ¦È)B2 . 
The theorems of this section and of the previous section can be immediately translated in terms of matrices. The matrix versions of these theorems is often used in applications so we brie.y present them in the section. 


16.5 Normal and Other Special Matrices 
First we consider real matrices. Recall the following de.nitions. 
De.nition 16.3. Given a real m ¡Á n matrix A, the transpose AT of A is the n ¡Á m matrix AT =(aT ) de.ned such that 
ij
T 
a
ij = aji for all i, j,1 ¡Ü i ¡Ü m,1 ¡Ü j ¡Ü n. A real n ¡Á n matrix A is 
. 	
normal if AAT = ATA, 

. 
symmetric if 


AT = A, 
. skew-symmetric if 
AT = .A, 
. orthogonal if 
AAT = ATA = In. 

Recall from Proposition 11.14 that when E is a Euclidean space and (e1,..., en) is an orthonormal basis for E, if A is the matrix of a linear map f : E ¡ú E w.r.t. the basis (e1,...,en), then AT is the matrix of the adjoint f. of f. Consequently, a normal linear map has a normal matrix, a self-adjoint linear map has a symmetric matrix, a skew-self-adjoint linear map has a skew-symmetric matrix, and an orthogonal linear map has an orthogonal matrix. 
Furthermore, if (u1,...,un) is another orthonormal basis for E and P is the change of basis matrix whose columns are the components of the ui w.r.t. the basis (e1,...,en), then P is orthogonal, and for any linear map f : E ¡ú E, if A is the matrix of f w.r.t (e1,...,en) and B is the matrix of f w.r.t. (u1,...,un), then 
B = P TAP. 
As a consequence, Theorems 16.12 and 16.14¨C16.16 can be restated as follows. 
Theorem 16.18. For every normal matrix A there is an orthogonal matrix P and a block diagonal matrix D such that A = PDP T, where D is of the form 
.
. 
D = 

.... 

D1 ... D2 ... 
.. .
.
....
.
.. . 
... Dp 
.... 

such that each block Dj is either a one-dimensional matrix (i.e., a real scalar) or a two-
dimensional matrix of the form  
Dj =  ¦Ëj ¦Ìj  .¦Ìj ¦Ëj  ,  
where ¦Ëj, ¦Ìj ¡Ê R, with ¦Ìj > 0.  

Theorem 16.19. For every symmetric matrix A there is an orthogonal matrix P and a diagonal matrix D such that A = PDP T, where D is of the form 
.
. 
D = 

.... 

¦Ë1 ... ¦Ë2 ... 
.. .
.
.. ..
.
.. . 
... ¦Ën 
....

, 

where ¦Ëi ¡Ê R. 
16.5. NORMAL AND OTHER SPECIAL MATRICES 
Theorem 16.20. For every skew-symmetric matrix A there is an orthogonal matrix P and a block diagonal matrix D such that A = PDP T, where D is of the form 
.
. 
D = 

.... 

D1 ... D2 ... 
.. .
.
....
.
.. . 
... Dp 
.... 

such that each block Dj is either 0 or a two-dimensional matrix of the form 0 .¦Ìj
Dj = , 
¦Ìj 0 
where ¦Ìj ¡Ê R, with ¦Ìj > 0. In particular, the eigenvalues of A are pure imaginary of the form ¡Ài¦Ìj, or 0. 
Theorem 16.21. For every orthogonal matrix A there is an orthogonal matrix P and a block diagonal matrix D such that A = PDP T, where D is of the form 
.
. 
D = 

.... 

D1 ... D2 ... 
.. .
.
....
.
.. . 
... Dp 
.... 

such that each block Dj is either 1, .1, or a two-dimensional matrix of the form 
cos ¦Èj . sin ¦Èj
Dj = 
sin ¦Èj cos ¦Èj 
where 0 <¦Èj <¦Ð. In particular, the eigenvalues of A are of the form cos ¦Èj ¡À i sin ¦Èj, 1, or .1. 
Theorem 16.21 can be used to show that the exponential map exp: so(n) ¡ú SO(n) is surjective; see Gallier [25]. 
We now consider complex matrices. 
WJ
De.nition 16.4. Given a complex m ¡Á n matrix A, the transpose AT matrix AT = aT de.ned such that 
ij 
T 
a
ij = aji for all i, j,1 ¡Ü i ¡Ü m,1 ¡Ü j ¡Ü n. The conjugate A of A is the m ¡Á n matrix A =(bij) de.ned such that bij = aij for all i, j,1 ¡Ü i ¡Ü m,1 ¡Ü j ¡Ü n. Given an m ¡Á n complex matrix A, the adjoint A. of A is the matrix de.ned such that A . =(AT)=(A)T . A complex n ¡Á n matrix A is 
of A is the n ¡Á m 

. 	
normal if 
AA . = A . A, 


. 
Hermitian if 


A . 
= A, 
. skew-Hermitian if 
A . 
= .A, 
. 	unitary if 
AA . = A . A = In. 

Recall from Proposition 13.15 that when E is a Hermitian space and (e1,..., en) is an orthonormal basis for E, if A is the matrix of a linear map f : E ¡ú E w.r.t. the basis (e1,...,en), then A. is the matrix of the adjoint f. of f. Consequently, a normal linear map has a normal matrix, a self-adjoint linear map has a Hermitian matrix, a skew-self-adjoint linear map has a skew-Hermitian matrix, and a unitary linear map has a unitary matrix. 
Furthermore, if (u1,...,un) is another orthonormal basis for E and P is the change of basis matrix whose columns are the components of the ui w.r.t. the basis (e1,...,en), then P is unitary, and for any linear map f : E ¡ú E, if A is the matrix of f w.r.t (e1,...,en) and B is the matrix of f w.r.t. (u1,...,un), then 
B = P . AP. 
Theorem 16.13 and Proposition 16.7 can be restated in terms of matrices as follows. 
Theorem 16.22. For every complex normal matrix A there is a unitary matrix U and a diagonal matrix D such that A = UDU. . Furthermore, if A is Hermitian, then D is a real matrix; if A is skew-Hermitian, then the entries in D are pure imaginary or zero; and if A is unitary, then the entries in D have absolute value 1. 


16.6 	Rayleigh¨CRitz Theorems and Eigenvalue Interlac-ing 
A fact that is used frequently in optimization problems is that the eigenvalues of a symmetric matrix are characterized in terms of what is known as the Rayleigh ratio, de.ned by 
xTAx 
R(A)(x)= T,x ¡Ê Rn ,x =0. 
xx 
16.6. RAYLEIGH¨CRITZ THEOREMS AND EIGENVALUE INTERLACING 
The following proposition is often used to prove the correctness of various optimization or approximation problems (for example PCA; see Section 21.4). It is also used to prove Proposition 16.25, which is used to justify the correctness of a method for graph-drawing (see Chapter 19). 
Proposition 16.23. (Rayleigh¨CRitz) If A is a symmetric n ¡Á n matrix with eigenvalues ¦Ë1 ¡Ü ¦Ë2 ¡Ü ¡¤ ¡¤¡¤ ¡Ü ¦Ën and if (u1,...,un) is any orthonormal basis of eigenvectors of A, where ui is a unit eigenvector associated with ¦Ëi, then 
xTAx 
max = ¦Ën 
xÊ± xT
=0 x 
(with the maximum attained for x = un), and 
xTAx 
xÊ±=0,x¡Ê{umax n.k+1,...,un}¡Í xTx = ¦Ën.k 
(with the maximum attained for x = un.k), where 1 ¡Ü k ¡Ü n . 1. Equivalently, if Vk is the subspace spanned by (u1,...,uk), then 
xTAx ¦Ëk = xÊ± max xT,k =1, . . . , n. 
=0,x¡ÊVk x 
Proof. First observe that 
xTAx 
max = max {x TAx | x T x =1}, 
xÊ± xTx
=0 x 
and similarly, 
xTAx Ê± max xT= max  x TAx | (x ¡Ê{un.k+1,...,un}¡Í) ¡Ä (x T x = 1) . 
x=0,x¡Ê{un.k+1,...,un}¡Í x x
Since A is a symmetric matrix, its eigenvalues are real and it can be diagonalized with respect to an orthonormal basis of eigenvectors, so let (u1,...,un) be such a basis. If we write 
n
 
x = xiui, i=1 
a simple computation shows that 
n
 
x TAx = ¦Ëixi 2 . i=1 
N 
If xTx = 1, then n x2 = 1, and since we assumed that ¦Ë1 ¡Ü ¦Ë2 ¡Ü ¡¤ ¡¤¡¤ ¡Ü ¦Ën, we get 
i=1 i 
nn 
x TAx = 2 ¡Ü ¦Ën x 2 = ¦Ën.
¦Ëixii i=1 i=1 
Thus, max x TAx | x T x =1 ¡Ü ¦Ën, 
x and since this maximum is achieved for en = (0, 0,..., 1), we conclude that 
max x TAx | x T x =1 = ¦Ën. 
x 
Next observe that x ¡Ê{un.k+1,...,un}¡Í and xTx = 1 i. xn.k+1 = ¡¤¡¤¡¤ = xn = 0 and Nn.k x2 = 1. Consequently, for such an x, we have 
i=1 i 
n.kn.k x TAx = ¦Ëix 2 i ¡Ü ¦Ën.k x 2 i = ¦Ën.k. i=1 i=1 
Thus, Ê± T¡ÍT|¡Ê{}¡Ä¡ÜAx ()(=1) ¦Ëmax xxu,...,uxx,..k+1knnnÊ± 
x and since this maximum is achieved for en.k = (0,..., 0, 1, 0,..., 0) with a 1 in position 
. k,weconcludethat n Ê± 
max x TAx | (x ¡Ê{un.k+1,...,un}¡Í) ¡Ä (x T x =1) = ¦Ën.k, 
x 
Ê± 
as claimed. 
For our purposes we need the version of Proposition 16.23 applying to min instead of max, whose proof is obtained by a trivial modi.cation of the proof of Proposition 16.23. 
Proposition 16.24. (Rayleigh¨CRitz) If A is a symmetric n ¡Á n matrix with eigenvalues ¦Ë1 ¡Ü ¦Ë2 ¡Ü ¡¤ ¡¤¡¤ ¡Ü ¦Ën and if (u1,...,un) is any orthonormal basis of eigenvectors of A, where ui is a unit eigenvector associated with ¦Ëi, then 
xTAx 
min = ¦Ë1 
x=0 xTx 
(with the minimum attained for x = u1), and 
xTAx 
min T= ¦Ëi 
x=0,x¡Ê{u1,...,ui.1}¡Í xx 
(with the minimum attained for x = ui), where 2 ¡Ü i ¡Ü n. Equivalently, if Wk = Vk¡Í.1 denotes the subspace spanned by (uk,...,un) (with V0 = (0)), then 
xTAx xTAx 
¦Ëk = min = min ,k =1, . . . , n. 
x=0,x¡ÊWk xTx x=0,x¡ÊV ¡Í xTx 
k.1 

16.6. RAYLEIGH¨CRITZ THEOREMS AND EIGENVALUE INTERLACING 
Propositions 16.23 and 16.24 together are known the Rayleigh¨CRitz theorem. 
As an application of Propositions 16.23 and 16.24, we prove a proposition which allows us to compare the eigenvalues of two symmetric matrices A and B = RTAR, where R is a rectangular matrix satisfying the equation RTR = I. 
De.nition 16.5. Given an n ¡Á n symmetric matrix A and an m ¡Á m symmetric B, with m ¡Ü n, if ¦Ë1 ¡Ü ¦Ë2 ¡Ü ¡¤¡¤¡¤ ¡Ü ¦Ën are the eigenvalues of A and ¦Ì1 ¡Ü ¦Ì2 ¡Ü ¡¤¡¤¡¤ ¡Ü ¦Ìm are the eigenvalues of B, then we say that the eigenvalues of B interlace the eigenvalues of A if 
¦Ëi ¡Ü ¦Ìi ¡Ü ¦Ën.m+i,i =1, . . . , m. 
For example, if n = 5 and m = 3, we have 
Firstweneedade.nition. ¡Ü¡Ü¦Ë¦Ë¦Ì113 (a)Theeigenvaluesof B interlacetheeigenvaluesof A. ¡Ü¡Ü¡Ü¡Ü¡Ü¡Ü(b)If ¦Ë¦Ë¦Ëaretheeigenvaluesof A andarethe ¡¤¡¤¡¤¡¤¡¤¡¤¦Ì¦Ì¦Ì1212nm eigenvaluesof B,andif ¦Ë,thenthereisaneigenvector of B witheigenvalue = ¦Ìviisuchthat Rv isaneigenvectorof A witheigenvalue ¦Ë¦Ì.iiProof. (a)Let()beanorthonormalbasisofeigenvectorsfor A,andlet()u,...,uv,...,v11nmbeanorthonormalbasisofeigenvectorsfor BLet Ubethesubspacespannedby()u,...,u. j 1jandlet Vbethesubspacespannedby().Forany i,thesubspace Vhasdimension v,...,vj 1ji T.i andthesubspace RUhasdimensionatmost i 1.Therefore,thereissomenonzero .i1 ¡ÍT¡Êwehave Rv (U)ByProposition16.24andusingthefactthat RRI,wehave =..i1Ê± 
¦Ë2 ¡Ü ¦Ì2 ¡Ü ¦Ë4 
¦Ë3 ¡Ü ¦Ì3 ¡Ü ¦Ë5. 
Proposition 16.25. Let A be an n ¡Á n symmetric matrix, R be an n ¡Á m matrix such that RTR = I (with m ¡Ü n), and let B = RTAR (an m ¡Á m matrix). The following properties hold: 
Ê± 
vector v ¡Ê Vi ¡É (RTUi.1)¡Í, and since 
v TRT uj =(Rv)T uj =0,j =1,...,i . 1, 
(Rv)TARv vTBv 
¦Ëi ¡Ü = . 
(Rv)TRv vTv 
On the other hand, by Proposition 16.23, 
xTBx xTBx 
¦Ìi == , 
x=0,x¡Ê{vi+1,max...,vn}¡Í xTx x=0,x¡Ê{v1,max...,vi} xTx 
so 
wTBw 
¡Ü ¦Ìi for all w ¡Ê Vi, 
wTw and since v ¡Ê Vi, we have vTBv 
¦Ëi ¡Ü T¡Ü ¦Ìi,i =1, . . . , m. 
vv We can apply the same argument to the symmetric matrices .A and .B, to conclude that 
.¦Ën.m+i ¡Ü.¦Ìi, 
that is, ¦Ìi ¡Ü ¦Ën.m+i,i =1, . . . , m. Therefore, ¦Ëi ¡Ü ¦Ìi ¡Ü ¦Ën.m+i,i =1, . . . , m, as desired. 
(b) If ¦Ëi = ¦Ìi, then 
(Rv)TARv vTBv 

¦Ëi = == ¦Ìi,
(Rv)TRv vTv so v must be an eigenvector for B and Rv must be an eigenvector for A, both for the eigenvalue ¦Ëi = ¦Ìi. 
Proposition 16.25 immediately implies the Poincar¡äe separation theorem. It can be used in situations, such as in quantum mechanics, where one has information about the inner 
T
products ui Auj. 
Proposition 16.26. (Poincar¡äe separation theorem) Let A be a n ¡Á n symmetric (or Her-mitian) matrix, let r be some integer with 1 ¡Ü r ¡Ü n, and let (u1,...,ur) be r orthonormal 
T
vectors. Let B =(ui Auj) (an r ¡Á r matrix), let ¦Ë1(A) ¡Ü ... ¡Ü ¦Ën(A) be the eigenvalues of A and ¦Ë1(B) ¡Ü ... ¡Ü ¦Ër(B) be the eigenvalues of B; then we have ¦Ëk(A) ¡Ü ¦Ëk(B) ¡Ü ¦Ëk+n.r(A),k =1, . . . , r. Observe that Proposition 16.25 implies that 
¦Ë1 + ¡¤¡¤¡¤ + ¦Ëm ¡Ü tr(RTAR) ¡Ü ¦Ën.m+1 + ¡¤¡¤¡¤ + ¦Ën. If P1 is the the n ¡Á (n . 1) matrix obtained from the identity matrix by dropping its last column, we have P1 TP1 = I, and the matrix B = P1 TAP1 is the matrix obtained from A by deleting its last row and its last column. In this case the interlacing result is 
¦Ë1 ¡Ü ¦Ì1 ¡Ü ¦Ë2 ¡Ü ¦Ì2 ¡Ü ¡¤ ¡¤¡¤ ¡Ü ¦Ìn.2 ¡Ü ¦Ën.1 ¡Ü ¦Ìn.1 ¡Ü ¦Ën, 
a genuine interlacing. We obtain similar results with the matrix Pn.r obtained by dropping the last n . r columns of the identity matrix and setting B = PnT.rAPn.r (B is the r ¡Á r matrix obtained from A by deleting its last n . r rows and columns). In this case we have the following interlacing inequalities known as Cauchy interlacing theorem: 
¦Ëk ¡Ü ¦Ìk ¡Ü ¦Ëk+n.r,k =1, . . . , r. (.) 


sults 
xTAx 
16.7.THECOURANT¨CFISCHERTHEOREM;PERTURBATIONRESULTS 16.7TheCourant¨CFischerTheorem;PerturbationRe-AnotherusefultooltoproveeigenvalueequalitiesistheCourant¨CFischercharacterizationof theeigenvaluesofasymmetricmatrix,alsoknownastheMin-max(andMax-min)theorem. Theorem16.27. (Courant¨CFischer) Let A beasymmetric matrixwitheigenvalues ¡Ánn Rn¡Ü¡Ü¡¤ ¡¤¡Ü V¦Ë¦Ë¦ËIfdenotesthesetofsubspacesof ofdimension k,then ¡¤ .12 kn¦Ëmin=maxk T¡ÊVW ¡ÊW,x=0 xxx.k+1nTAxx¦Ëmin=max .k byProposition16.23,wehave ¡Ý min¦Ë=maxmax .k TT¡ÊVÊ±¡ÊW ¡ÊÊ±=0VW,x=0xxxxx,xxkk Therefore,weneedtoprovethereverseinequality;thatis,wehavetoshowthat ¡Í¡Í¡ÊV¡É¡Ê¡ÉNowforany W ,ifwecanprovethat WV =(0),thenforanynonzero WVvk¡Ü¡Ü¦Ëmin=max .k Ê± ¡Í¡Í¡É¡Ý.Itremainstoprovethatdim(WV )1.However,dim(V)= k 1,sodim(V )= .k1..k1k1. k +1,andbyhypothesisdim(W )= kBytheGrassmannrelation, n . ¡Í¡Í¡Í¡Édim(W )+dim(V )=dim(WV )+dim(WV )+ , ¡Í Rn¡Üandsincedim(WV )dim()= +,weget n.k1eigenvaluesofasymmetricmatrixduetoHermannWeyl. 
T
W ¡ÊVk x¡ÊW,x=0 xx 
Ê± 
Ê± Proof. Letusconsiderthesecondequality,theproofofthe.rstequalitybeingsimilar.Let ()beanyorthonormalbasisofeigenvectorsof A,where isauniteigenvector u,...,uu1in
associated with ¦Ëi. Observe that the space Vk spanned by (u1,...,uk) has dimension k, and xTAx xTAx 
Ê± 
Ê± 

xTAx 
¦Ëk ¡Ü max , for all W ¡ÊVk.
T
x=0,x¡ÊW xx k.1 k.1, 
by Proposition 16.24 , we have 
xTAx vTAv xTAx x=0,x¡ÊV ¡Í xTxvTv x¡ÊW,x=0 xTx 
k.1 
k.1k.1k.1
k + n . k +1 ¡Ü dim(W ¡É Vk¡Í.1)+ n; that is, 1 ¡Ü dim(W ¡É Vk¡Í.1), as claimed. The Courant¨CFischer theorem yields the following useful result about perturbing the 
Proposition 16.28. Given two n ¡Án symmetric matrices A and B = A +¦¤A, if ¦Á1 ¡Ü ¦Á2 ¡Ü ¡¤¡¤¡¤ ¡Ü ¦Án are the eigenvalues of A and ¦Â1 ¡Ü ¦Â2 ¡Ü ¡¤ ¡¤¡¤ ¡Ü ¦Ân are the eigenvalues of B, then 
|¦Ák . ¦Âk|¡Ü ¦Ñ(¦¤A)¡Ü½Ð ¦¤A½Ð2 ,k =1, . . . , n. 
Proof. Let Vk be de.ned as in the Courant¨CFischer theorem and let Vk be the subspace spanned by the k eigenvectors associated with ¦Ë1,...,¦Ëk. By the Courant¨CFischer theorem applied to B, we have 
xTBx 
¦Âk = min max 
T
W ¡ÊVk x¡ÊW,x=0 xx xTBx 
¡Ü max 
T
x¡ÊVk xx xTAx xT¦¤ Ax 
= max + 
TT
x¡ÊVk xx xx xTAx xT¦¤Ax 
¡Ü max + max . 
Ê± 
x¡ÊVk xTx x¡ÊVk xTx 
By Proposition 16.23, we have 
xTAx 
¦Ák = max ,
T
x¡ÊVk xx 
so we obtain xTAx xT¦¤Ax 
¦Âk ¡Ü max + max 
TT
x¡ÊVk xx x¡ÊVk xx xT¦¤Ax 
= ¦Ák + max 
x¡ÊVk xTx xT¦¤Ax 
¡Ü ¦Ák + max .
T
x¡ÊRn xx 
Now by Proposition 16.23 and Proposition 8.9, we have 
max xT¦¤TAx = max ¦Ëi(¦¤A) ¡Ü ¦Ñ(¦¤A)¡Ü½Ð ¦¤A½Ð2 , 
x¡ÊRn xx i where ¦Ëi(¦¤A) denotes the ith eigenvalue of ¦¤A, which implies that ¦Âk ¡Ü ¦Ák + ¦Ñ(¦¤A) ¡Ü ¦Ák +½Ð¦¤A½Ð2 . By exchanging the roles of A and B, we also have ¦Ák ¡Ü ¦Âk + ¦Ñ(¦¤A) ¡Ü ¦Âk +½Ð¦¤A½Ð2 , and thus, |¦Ák . ¦Âk|¡Ü ¦Ñ(¦¤A)¡Ü½Ð ¦¤A½Ð2 ,k =1, . . . , n, as claimed. 
n 
(¦Ák . ¦Âk)2¡Ü½Ð ¦¤A½ÐF 2 , 
k=1 
proof; see Lax [44]. 
Ê± 
Ê± 
1. If i + j = k +1, then 
¦Ëi(A)+ ¦Ëj(B) ¡Ü ¦Ëk(A + B). 


2. If i + j = k + n, then 


16.7.THECOURANT¨CFISCHERTHEOREM;PERTURBATIONRESULTS Proposition16.28alsoholdsforHermitianmatrices. AprettyresultofWielandtandHo.manassertsthat ½Ð½Ð whereistheFrobeniusnorm.However,theproofissigni.cantlyharderthantheabove F TheCourant¨CFischertheoremcanalsobeusedtoprovesomefamousinequalitiesdueto HermannWeyl.Thesecanalsobeviewedasperturbationresults.Giventwosymmetric(or Hermitian)matrices A and B,let ¦Ë(A),¦Ë(B),and ¦Ë(AB)denotethe itheigenvalueof +iiiA,B,and AB,respectively,arrangedinnondecreasingorder. + Proposition16.29. (Weyl)Giventwosymmetric(orHermitian) matrices A and B¡Ánn , ¡Ü¡Üthefollowinginequalitieshold:Forall i,j,k with 1 i,j,k :n¡Ü¦Ë(AB) ¦Ë(A)+ ¦Ë(B)+ .kij¦Ë(AB)= min+ .k¡ÊH,x=0 xxxTTAxBx xx¦Ë(A)= ¦Ë(B)= maxmax,.ijÊ±¡Ê¡ÊF,x=0 G,x=0xxxxxx¡É¡É¡É.¡Édim(FGH)=dim(F )+dim(GH)dim(F +(GH)) ¡É.¡Ý.dim(GH)=dim(G)+dim(H)dim(GH)dim(G)+dim(H)+ n, 
Proof. Observe that the .rst set of inequalities is obtained form the second set by replacing A by .A and B by .B, so it is enough to prove the second set of inequalities. By the Courant¨CFischer theorem, there is a subspace H of dimension n . k + 1 such that 
xT(A + B)x T
Similarly, there exists a subspace F of dimension i and a subspace G of dimension j such that 
TT
We claim that F ¡É G ¡É H = (0). To prove this, we use the Grassmann relation twice. First, 
¡Ý dim(F ) + dim(G ¡É H) . n, 
and second, 
so 
dim(F ¡É G ¡É H) ¡Ý dim(F ) + dim(G) + dim(H) . 2n. 

However, dim(F ) + dim(G) + dim(H)= i + j + n . k +1 and i + j = k + n, so we have dim(F ¡É G ¡É H) ¡Ý i + j + n . k +1 . 2n = k + n + n . k +1 . 2n =1, which shows that F ¡É G ¡É H = (0). Then for any unit vector z ¡Ê F ¡É G ¡É H = (0), we have 
¦Ëk(A + B) ¡Ü z T(A + B)z, ¦Ëi(A) ¡Ý z TAz, ¦Ëj(B) ¡Ý z TBz, 
establishing the desired inequality ¦Ëk(A + B) ¡Ü ¦Ëi(A)+ ¦Ëj(B). 
In the special case i = j = k, we obtain 
¦Ë1(A)+ ¦Ë1(B) ¡Ü ¦Ë1(A + B),¦Ën(A + B) ¡Ü ¦Ën(A)+ ¦Ën(B). It follows that ¦Ë1 (as a function) is concave, while ¦Ën (as a function) is convex. If i = 1 and j = k, we obtain 
¦Ë1(A)+ ¦Ëk(B) ¡Ü ¦Ëk(A + B), and if i = k and j = n, we obtain 
¦Ëk(A + B) ¡Ü ¦Ëk(A)+ ¦Ën(B), 
and combining them, we get 
¦Ë1(A)+ ¦Ëk(B) ¡Ü ¦Ëk(A + B) ¡Ü ¦Ëk(A)+ ¦Ën(B). 
In particular, if B is positive semide.nite, since its eigenvalues are nonnegative, we obtain the following inequality known as the monotonicity theorem for symmetric (or Hermitian) matrices: if A and B are symmetric (or Hermitian) and B is positive semide.nite, then 
¦Ëk(A) ¡Ü ¦Ëk(A + B) k =1, . . . , n. 
The reader is referred to Horn and Johnson [36] (Chapters 4 and 7) for a very complete treatment of matrix inequalities and interlacing results, and also to Lax [44] and Serre [57]. 

16.8 Summary 
The main concepts and results of this chapter are listed below: 
. 	
Normal linear maps, self-adjoint linear maps, skew-self-adjoint linear maps, and or-thogonal linear maps. 

. 	
Properties of the eigenvalues and eigenvectors of a normal linear map. 

16.9. PROBLEMS 

. 	The complexi.cation of a real vector space, of a linear map, and of a Euclidean inner product. 

. 	The eigenvalues of a self-adjoint map in a Hermitian space are real. 

. 	The eigenvalues of a self-adjoint map in a Euclidean space are real. 

. 	Every self-adjoint linear map on a Euclidean space has an orthonormal basis of eigen-vectors. 

. 	Every normal linear map on a Euclidean space can be block diagonalized (blocks of size at most 2 ¡Á 2) with respect to an orthonormal basis of eigenvectors. 

. 	Every normal linear map on a Hermitian space can be diagonalized with respect to an orthonormal basis of eigenvectors. 

. 	The spectral theorems for self-adjoint, skew-self-adjoint, and orthogonal linear maps (on a Euclidean space). 

. 	The spectral theorems for normal, symmetric, skew-symmetric, and orthogonal (real) matrices. 

. 	The spectral theorems for normal, Hermitian, skew-Hermitian, and unitary (complex) matrices. 

. 	The Rayleigh ratio and the Rayleigh¨CRitz theorem. 

. 	Interlacing inequalities and the Cauchy interlacing theorem. 

. 	The Poincar¡äe separation theorem. 

. 	The Courant¨CFischer theorem. 

. 	Inequalities involving perturbations of the eigenvalues of a symmetric matrix. 

. 	The Weyl inequalities. 


16.9 Problems 
Problem 16.1. Prove that the structure EC introduced in De.nition 16.2 is indeed a com-plex vector space. 
Problem 16.2. Prove that the formula 
(u1 + iv1,u2 + iv2£©C = (u1,u2£© + (v1,v2£© + i((v1,u2£©.( u1,v2£©) 
de.nes a Hermitian form on EC that is positive de.nite and that (.,.£©C agrees with (.,.£©on real vectors. 
Problem 16.3. Given any linear map f : E ¡ú E, prove the map fC . de.ned such that 
fC. (u + iv)= f . (u)+ if . (v) 
for all u, v ¡Ê E is the adjoint of fC w.r.t. (.,.£©C. 
Problem 16.4. Let A be a real symmetric n ¡Á n matrix whose eigenvalues are nonnegative. Prove that for every p> 0, there is a real symmetric matrix S whose eigenvalues are nonnegative such that Sp = A. 
Problem 16.5. Let A be a real symmetric n ¡Á n matrix whose eigenvalues are positive. 
(1) Prove that there is a real symmetric matrix S such that A = eS . 
(2) Let S be a real symmetric n ¡Á n matrix. Prove that A = eS is a real symmetric n ¡Á n matrix whose eigenvalues are positive. 
Problem 16.6. Let A be a complex matrix. Prove that if A can be diagonalized with respect to an orthonormal basis, then A is normal. 
Problem 16.7. Let f : Cn ¡ú Cn be a linear map. 
(1) Prove that if f is diagonalizable and if ¦Ë1,...,¦Ën are the eigenvalues of f, then ¦Ë21,...,¦Ë2 n are the eigenvalues of f2, and if ¦Ë2 i = ¦Ë2 j implies that ¦Ëi = ¦Ëj, then f and f2 have the same eigenspaces. 

(2) Let f and g be two real self-adjoint linear maps f, g : Rn ¡ú Rn . Prove that if f and g have nonnegative eigenvalues (f and g are positve semide.nite) and if f2 = g2, then f = g. 


Problem 16.8. (1) Let so(3) be the space of 3 ¡Á 3 skew symmetric matrices 
.. 

.

..

    

a, b, c ¡Ê R 

.. 

. 

.

. 
0 .cb 
so(3) = 

0 .a

c 

.
.ba 0 
For any matrix 

.
. 
0  .c  b  
A = .c  0  .a. ¡Ê so(3),  
.b  a  0  
¡Ì  

if we let ¦È = a2 + b2 + c2, recall from Section 11.7 (the Rodrigues formula) that the expo-nential map exp: so(3) ¡ú SO(3) is given by 
A sin ¦È (1 . cos ¦È)
A2 
e = I3 + A + , if ¦È =0,
¦È¦È2 
with exp(03)= I3. 
(2) Prove that eA is an orthogonal matrix of determinant +1, i.e., a rotation matrix. 
(3) Prove that the exponential map exp: so(3) ¡ú SO(3) is surjective. For this proceed as follows: Pick any rotation matrix R ¡Ê SO(3); 
16.9. PROBLEMS 
(1) 
The case R = I is trivial. 

(2) 
If R = I and tr(R)= .1, then 


 	 
exp .1(R)=¦È (R . RT ) 1+2cos ¦È = tr(R). 
2 sin ¦È 
(Recall that tr(R)= r11 + r22 + r33, the trace of the matrix R). 
Show that there is a unique skew-symmetric B with corresponding ¦È satisfying 0 < 
B
¦È<¦Ð such that e= R. 
(3) If 	R = I and tr(R)= .1, then prove that the eigenvalues of R are 1, .1, .1, that R = RT, and that R2 = I. Prove that the matrix 
1 
S =(R . I)
2
is a symmetric matrix whose eigenvalues are .1, .1, 0. Thus S can be diagonalized with respect to an orthogonal matrix Q as 
.. 
.100 S = Q . 0 .10. QT . 0 00 
Prove that there exists a skew symmetric matrix 
.. 
0 .dc 
..
U = d 0 .b .cb 0 
so that 
U2 = S = 1(R . I). 

2
Observe that .. 
.(c2 + d2) bc bd 
U2 = . bc .(b2 + d2) cd ., 
bd cd .(b2 + c2) 

and use this to conclude that if U2 = S, then b2 + c2 + d2 = 1. Then show that 
.. .. 
. 0 .dc . exp .1(R)= (2k + 1)¦Ð . d 0 .b.,k ¡Ê Z ,
.	.
.cb 0 
where (b, c, d) is any unit vector such that for the corresponding skew symmetric matrix U, we have U2 = S. 
(4) To .nd a skew symmetric matrix U so that U2 = S = 12 (R . I) as in (3), we can solve the system 
.. 
b2 . 1 bc bd . bc c2 . 1 cd . = S. bd cdd2 . 1 
We immediately get b2,c2,d2, and then, since one of b, c, d is nonzero, say b, if we choose the positive square root of b2, we can determine c and d from bc and bd. 
Implement a computer program in Matlab to solve the above system. 
Problem 16.9. It was shown in Proposition 14.15 that the exponential map is a map exp: so(n) ¡ú SO(n), where so(n) is the vector space of real n ¡Á n skew-symmetric matrices. Use the spectral theorem to prove that the map exp: so(n) ¡ú SO(n) is surjective. 
Problem 16.10. Let u(n) be the space of (complex) n ¡Á n skew-Hermitian matrices (B. = .B) and let su(n) be its subspace consisting of skew-Hermitian matrice with zero trace (tr(B) = 0). 
(1) Prove that if B ¡Ê u(n), then eB ¡Ê U(n), and if if B ¡Ê su(n), then eB ¡Ê SU(n). Thus we have well-de.ned maps exp: u(n) ¡ú U(n) and exp: su(n) ¡ú SU(n). 
(2) 
Prove that the map exp: u(n) ¡ú U(n) is surjective. 

(3) 
Prove that the map exp: su(n) ¡ú SU(n) is surjective. 


Problem 16.11. Recall that a matrix B ¡Ê Mn(R) is skew-symmetric if BT = .B. Check that the set so(n) of skew-symmetric matrices is a vector space of dimension n(n . 1)/2, and thus is isomorphic to Rn(n.1)/2 . 
(1) Given a rotation matrix 
cos ¦È . sin ¦È 
R = ,
sin ¦È cos ¦È 
where 0 <¦È<¦Ð, prove that there is a skew symmetric matrix B such that 
R =(I . B)(I + B).1 . 
(2) Prove that the eigenvalues of a skew-symmetric matrix are either 0 or pure imaginary 
(that is, of the form i¦Ì for ¦Ì ¡Ê R.). Let C : so(n) ¡ú Mn(R) be the function (called the Cayley transform of B) given by 
C(B)=(I . B)(I + B).1 . 
Prove that if B is skew-symmetric, then I . B and I + B are invertible, and so C is well-de.ned. Prove that 
(I + B)(I . B)=(I . B)(I + B), 
and that (I + B)(I . B).1 =(I . B).1(I + B). 
16.9. PROBLEMS 
Prove that (C(B))TC(B)= I 
and that det C(B) = +1, 
so that C(B) is a rotation matrix. Furthermore, show that C(B) does not admit .1 as an eigenvalue. 
(3) Let SO(n) be the group of n ¡Á n rotation matrices. Prove that the map 
C : so(n) ¡ú SO(n) 
is bijective onto the subset of rotation matrices that do not admit .1 as an eigenvalue. Show that the inverse of this map is given by 
B =(I + R).1(I . R)=(I . R)(I + R).1 , 
where R ¡Ê SO(n) does not admit .1 as an eigenvalue. 
Problem 16.12. Please refer back to Problem 3.6. Let ¦Ë1,...,¦Ën be the eigenvalues of A (not necessarily distinct). Using Schur¡¯s theorem, A is similar to an upper triangular matrix B, that is, A = P BP .1 with B upper triangular, and we may assume that the diagonal entries of B in descending order are ¦Ë1,...,¦Ën. 
(1) If the Eij are listed according to total order given by 
i = h and j>k (i, j) < (h, k) i. 
or i<h. 
prove that RB is an upper triangular matrix whose diagonal entries are 
(¦Ën,...,¦Ë1,...,¦Ën,...,¦Ë1),
¦Ë¡­ ¡£ 
2
n
and that LB is an upper triangular matrix whose diagonal entries are 
(¦Ë1,...,¦Ë1 ...,¦Ën,...,¦Ën).
¦Ë¡­ ¡£ ¦Ë¡­ ¡£ 
nn 
Hint. Figure out what are RB(Eij)= EijB and LB(Eij)= BEij. 
(2) Use the fact that = R.1
LA = LP . LB . LP .1 ,RAP . RB . RP , 
to express adA = LA . RA in terms of LB . RB, and conclude that the eigenvalues of adA are ¦Ëi . ¦Ëj, for i =1,...,n, and for j = n, . . . , 1. 


