Chapter 16 
Spectral Theorems in Euclidean and Hermitian Spaces 
16.1 	Introduction 
The goal of this chapter is to show that there are nice normal forms for symmetric matrices, skew-symmetric matrices, orthogonal matrices, and normal matrices. The spectral theorem for symmetric matrices states that symmetric matrices have real eigenvalues and that they can be diagonalized over an orthonormal basis. The spectral theorem for Hermitian matrices states that Hermitian matrices also have real eigenvalues and that they can be diagonalized over a complex orthonormal basis. Normal real matrices can be block diagonalized over an orthonormal basis with blocks having size at most two and there are re.nements of this normal form for skew-symmetric and orthogonal matrices. 
The spectral result for real symmetric matrices can be used to prove two characterizations of the eigenvalues of a symmetric matrix in terms of the Rayleigh ratio. The .rst charac-terization is the Rayleigh�CRitz theorem and the second one is the Courant�CFischer theorem. Both results are used in optimization theory and to obtain results about perturbing the eigenvalues of a symmetric matrix. 
In this chapter all vector spaces are .nite-dimensional real or complex vector spaces. 

16.2 	Normal Linear Maps: Eigenvalues and Eigenvec-tors 
We begin by studying normal maps, to understand the structure of their eigenvalues and eigenvectors. This section and the next three were inspired by Lang [41], Artin [3], Mac Lane and Birkho. [46], Berger [5], and Bertin [7]. 
De.nition 16.1. Given a Euclidean or Hermitian space E, a linear map f : E �� E is normal if 
f . f . = f . . f. 
543 

A linear map f : E �� E is self-adjoint if f = f. , skew-self-adjoint if f = .f., and orthogonal if f . f. = f. . f = id. 
Obviously, a self-adjoint, skew-self-adjoint, or orthogonal linear map is a normal linear map. Our .rst goal is to show that for every normal linear map f : E �� E, there is an orthonormal basis (w.r.t. (.,.��) such that the matrix of f over this basis has an especially nice form: it is a block diagonal matrix in which the blocks are either one-dimensional matrices (i.e., single entries) or two-dimensional matrices of the form 
�˦� ..�̦� 
This normal form can be further re.ned if f is self-adjoint, skew-self-adjoint, or orthog-onal. As a .rst step we show that f and f. have the same kernel when f is normal. 
Proposition 16.1. Given a Euclidean space E, if f : E �� E is a normal linear map, then Ker f = Ker f. . 
Proof. First let us prove that 
(f(u),f(v)�� = (f . (u),f . (v)�� 
for all u, v �� E. Since f. is the adjoint of f and f . f. = f. . f, we have 
(f(u),f(u)�� = (u, (f . . f)(u)��, = (u, (f . f . )(u)��, = (f . (u),f . (u)��. 
Since (.,.�� is positive de.nite, 
(f(u),f(u)�� = 0 i. f(u)=0, 
(f . (u),f . (u)�� = 0 i. f . (u)=0, 

and since (f(u),f(u)�� = (f . (u),f . (u)��, 
we have f(u)=0 i. f . (u)=0. 
Consequently, Ker f = Ker f. . 
Assuming again that E is a Hermitian space, observe that Proposition 16.1 also holds. We deduce the following corollary. 
Proposition 16.2. Given a Hermitian space E, for any normal linear map f : E �� E, we have Ker (f) �� Im(f) = (0). 
16.2. NORMAL LINEAR MAPS: EIGENVALUES AND EIGENVECTORS 
Proof. Assume v �� Ker (f) �� Im(f) = (0), which means that v = f(u) for some u �� E, and f(v) = 0. By Proposition 16.1, Ker (f) = Ker(f.), so f(v) = 0 implies that f.(v) = 0. Consequently, 
0= (f . (v),u�� = (v, f(u)�� = (v, v��, 
and thus, v = 0. 
We also have the following crucial proposition relating the eigenvalues of f and f. . 
Proposition 16.3. Given a Hermitian space E, for any normal linear map f : E �� E,a vector u is an eigenvector of f for the eigenvalue �� (in C) i. u is an eigenvector of f. for the eigenvalue ��. 
Proof. First it is immediately veri.ed that the adjoint of f . �� id is f. . �� id. Furthermore, f . �� id is normal. Indeed, 
(f . �� id) . (f . �� id) . =(f . �� id) . (f . . �� id), 
= f . f . . ��f . ��f . + �˦� id, 
= f . . f . ��f . . ��f + �˦� id, 
=(f . . �� id) . (f . �� id), 
=(f . �� id) . . (f . �� id). 
Applying Proposition 16.1 to f . �� id, for every nonnull vector u, we see that 
(f . �� id)(u)=0 i. (f . . �� id)(u)=0, 
which is exactly the statement of the proposition. 
The next proposition shows a very important property of normal linear maps: eigenvec-tors corresponding to distinct eigenvalues are orthogonal. 
Proposition 16.4. Given a Hermitian space E, for any normal linear map f : E �� E, if u and v are eigenvectors of f associated with the eigenvalues �� and �� (in C) where �� 
= ��, then (u, v�� =0. 
Proof. Let us compute (f(u),v�� in two di.erent ways. Since v is an eigenvector of f for ��, by Proposition 16.3, v is also an eigenvector of f. for ��, and we have 
(f(u),v�� = (��u, v�� = ��(u, v��, 
and 
(f(u),v�� = (u, f . (v)�� = (u, ��v�� = ��(u, v��, 
where the last identity holds because of the semilinearity in the second argument. Thus ��(u, v�� = ��(u, v��, that is, (�� . ��)(u, v�� =0, which implies that (u, v�� = 0, since �� = ��. We can show easily that the eigenvalues of a self-adjoint linear map are real. 
Proposition 16.5. Given a Hermitian space E, all the eigenvalues of any self-adjoint linear 
map f : E �� E are real. Proof. Let z (in C) be an eigenvalue of f and let u be an eigenvector for z. We compute (f(u),u�� in two di.erent ways. We have 
(f(u),u�� = (zu, u�� = z(u, u��, and since f = f., we also have (f(u),u�� = (u, f . (u)�� = (u, f(u)�� = (u, zu�� = z(u, u��. Thus, z(u, u�� = z(u, u��, which implies that z = z, since u = 0, and z is indeed real. There is also a version of Proposition 16.5 for a (real) Euclidean space E and a self-adjoint map f : E �� E since every real vector space E can be embedded into a complex vector space EC, and every linear map f : E �� E can be extended to a linear map fC: EC �� EC. De.nition 16.2. Given a real vector space E, let EC be the structure E �� E under the addition operation (u1,u2)+(v1,v2)=(u1 + v1,u2 + v2), and let multiplication by a complex scalar z = x + iy be de.ned such that (x + iy) �� (u, v)=(xu . yv, yu + xv). The space EC is called the complexi.cation of E. It is easily shown that the structure EC is a complex vector space. It is also immediate that (0,v)= i(v, 0), and thus, identifying E with the subspace of EC consisting of all vectors of the form (u, 0), we can write 
(u, v)= u + iv. 


16.2. NORMAL LINEAR MAPS: EIGENVALUES AND EIGENVECTORS 
Observe that if (e1,...,en) is a basis of E (a real vector space), then (e1,...,en) is also a basis of EC (recall that ei is an abbreviation for (ei, 0)). 
A linear map f : E �� E is extended to the linear map fC: EC �� EC de.ned such that 
fC(u + iv)= f(u)+ if(v). 
For any basis (e1,...,en) of E, the matrix M(f) representing f over (e1,...,en) is iden-tical to the matrix M(fC) representing fC over (e1,...,en), where we view (e1,...,en) as a basis of EC. As a consequence, det(zI . M(f)) = det(zI . M(fC)), which means that f and fC have the same characteristic polynomial (which has real coe.cients). We know that every polynomial of degree n with real (or complex) coe.cients always has n complex roots (counted with their multiplicity), and the roots of det(zI . M(fC)) that are real (if any) are the eigenvalues of f. 
Next we need to extend the inner product on E to an inner product on EC. 
The inner product (.,.�� on a Euclidean space E is extended to the Hermitian positive de.nite form (.,.��C on EC as follows: 
(u1 + iv1,u2 + iv2��C = (u1,u2�� + (v1,v2�� + i((v1,u2��.( u1,v2��). 
It is easily veri.ed that (.,.��C is indeed a Hermitian form that is positive de.nite, and it is clear that (.,.��C agrees with (.,.�� on real vectors. Then given any linear map 
f : E �� E, it is easily veri.ed that the map fC . de.ned such that 
fC. (u + iv)= f . (u)+ if . (v) 
for all u, v �� E is the adjoint of fC w.r.t. (.,.��C. 
Proposition 16.6. Given a Euclidean space E, if f : E �� E is any self-adjoint linear map, then every eigenvalue �� of fC is real and is actually an eigenvalue of f (which means that there is some real eigenvector u �� E such that f(u)= ��u). Therefore, all the eigenvalues of f are real. 
Proof. Let EC be the complexi.cation of E, (.,.��C the complexi.cation of the inner product (.,.�� on E, and fC: EC �� EC the complexi.cation of f : E �� E. By de.nition of fC and (.,.��C, if f is self-adjoint, we have 
(fC(u1 + iv1),u2 + iv2��C = (f(u1)+ if(v1),u2 + iv2��C = (f(u1),u2�� + (f(v1),v2�� 
+ i((u2,f(v1)��.( f(u1),v2��) = (u1,f(u2)�� + (v1,f(v2)�� 

+ i((f(u2),v1��.( u1,f(v2)��) = (u1 + iv1,f(u2)+ if(v2)��C = (u1 + iv1,fC(u2 + iv2)��C, 


which shows that fC is also self-adjoint with respect to (.,.��C. 
As we pointed out earlier, f and fC have the same characteristic polynomial det(zI.fC)= det(zI . f), which is a polynomial with real coe.cients. Proposition 16.5 shows that the zeros of det(zI . fC) = det(zI . f) are all real, and for each real zero �� of det(zI . f), the linear map ��id . f is singular, which means that there is some nonzero u �� E such that f(u)= ��u. Therefore, all the eigenvalues of f are real. 
Proposition 16.7. Given a Hermitian space E, for any linear map f : E �� E, if f is skew-self-adjoint, then f has eigenvalues that are pure imaginary or zero, and if f is unitary, then f has eigenvalues of absolute value 1. 
Proof. If f is skew-self-adjoint, f. = .f, and then by the de.nition of the adjoint map, for any eigenvalue �� and any eigenvector u associated with ��, we have 
��(u, u�� = (��u, u�� = (f(u),u�� = (u, f . (u)�� = (u, .f(u)�� = .(u, ��u�� = .��(u, u��, 
and since u = 0 and (.,.�� is positive de.nite, (u, u�� = 0, so 
�� = .��, 
which shows that �� = ir for some r �� R. 
If f is unitary, then f is an isometry, so for any eigenvalue �� and any eigenvector u associated with ��, we have 
|��|2(u, u�� = �˦�(u, u�� = (��u, ��u�� = (f(u),f(u)�� = (u, u��, 
and since u = 0, we obtain |��|2 = 1, which implies 
|��| =1. 


16.3 Spectral Theorem for Normal Linear Maps 
Given a Euclidean space E, our next step is to show that for every linear map f : E �� E there is some subspace W of dimension 1 or 2 such that f(W ) . W . When dim(W ) = 1, the subspace W is actually an eigenspace for some real eigenvalue of f. Furthermore, when f is normal, there is a subspace W of dimension 1 or 2 such that f(W ) . W and f.(W ) . W . The di.culty is that the eigenvalues of f are not necessarily real. One way to get around this problem is to complexify both the vector space E and the inner product (.,.�� as we did in Section 16.2. 
Given any subspace W of a Euclidean space E, recall that the orthogonal complement W �� of W is the subspace de.ned such that 
W �� = {u �� E |(u, w�� =0, for all w �� W }. 
16.3. SPECTRAL THEOREM FOR NORMAL LINEAR MAPS 
Recall from Proposition 11.11 that E = W �� W �� (this can be easily shown, for example, by constructing an orthonormal basis of E using the Gram�CSchmidt orthonormalization procedure). The same result also holds for Hermitian spaces; see Proposition 13.13. 
As a warm up for the proof of Theorem 16.12, let us prove that every self-adjoint map on a Euclidean space can be diagonalized with respect to an orthonormal basis of eigenvectors. 
Theorem 16.8. (Spectral theorem for self-adjoint linear maps on a Euclidean space) Given a Euclidean space E of dimension n, for every self-adjoint linear map f : E �� E, there is an orthonormal basis (e1,...,en) of eigenvectors of f such that the matrix of f w.r.t. this 
basis is a diagonal matrix

.
. 
.... 

��1 ... ��2 ... 
,
.. .
.
.. ..
.
.. . 
....
... ��n 
with ��i �� R. 
Proof. We proceed by induction on the dimension n of E as follows. If n = 1, the result is trivial. Assume now that n �� 2. From Proposition 16.6, all the eigenvalues of f are real, so pick some eigenvalue �� �� R, and let w be some eigenvector for ��. By dividing w by its norm, we may assume that w is a unit vector. Let W be the subspace of dimension 1 spanned by w. Clearly, f(W ) . W . We claim that f(W ��) . W ��, where W �� is the orthogonal complement of W . 
Indeed, for any v �� W ��, that is, if (v, w�� = 0, because f is self-adjoint and f(w)= ��w, we have 
(f(v),w�� = (v, f(w)�� 
= (v, ��w�� 
= ��(v, w�� =0 
since (v, w�� = 0. Therefore, f(W ��) . W �� . 
Clearly, the restriction of f to W �� is self-adjoint, and we conclude by applying the induction hypothesis to W �� (whose dimension is n . 1). 
We now come back to normal linear maps. One of the key points in the proof of Theorem 
16.8 is that we found a subspace W with the property that f(W ) . W implies that f(W ��) . W �� . In general, this does not happen, but normal maps satisfy a stronger property which ensures that such a subspace exists. 
The following proposition provides a condition that will allow us to show that a nor-mal linear map can be diagonalized. It actually holds for any linear map. We found the inspiration for this proposition in Berger [5]. 
Proposition 16.9. Given a Hermitian space E, for any linear map f : E �� E and any 
JW 
subspace W of E, if f(W ) . W , then f. W �� . W �� . Consequently, if f(W ) . W and
JW JW 
f.(W ) . W , then fW �� . W �� and f. W �� . W �� . 
Proof. If u �� W ��, then (w, u�� = 0 for all w �� W. 
However, (f(w),u�� = (w, f . (u)��, 
and f(W ) . W implies that f(w) �� W . Since u �� W ��, we get 
0= (f(w),u�� = (w, f . (u)��, 
which shows that (w, f.(u)�� = 0 for all w �� W , that is, f.(u) �� W �� . Therefore, we have f.(W ��) . W �� . 
JW 
We just proved that if f(W ) . W , then f. W �� . W �� . If we also have f.(W ) . W , then by applying the above fact to f. , we get f..(W ��) . W ��, and since f.. = f, this is just f(W ��) . W ��, which proves the second statement of the proposition. 
It is clear that the above proposition also holds for Euclidean spaces. 
Although we are ready to prove that for every normal linear map f (over a Hermitian space) there is an orthonormal basis of eigenvectors (see Theorem 16.13 below), we now return to real Euclidean spaces. 
Proposition 16.10. If f : E �� E is a linear map and w = u + iv is an eigenvector of fC: EC �� EC for the eigenvalue z = �� + i��, where u, v �� E and ��, �� �� R, then 
f(u)= ��u . ��v and f(v)= ��u + ��v. (.) 
As a consequence, 
fC(u . iv)= f(u) . if(v)=(�� . i��)(u . iv), 
which shows that w = u . iv is an eigenvector of fC for z = �� . i��. 
Proof. Since fC(u + iv)= f(u)+ if(v) 
and fC(u + iv)=(�� + i��)(u + iv)= ��u . ��v + i(��u + ��v), 
we have f(u)= ��u . ��v and f(v)= ��u + ��v. 
Using this fact, we can prove the following proposition. 


16.3. SPECTRAL THEOREM FOR NORMAL LINEAR MAPS 
Proposition 16.11. Given a Euclidean space E, for any normal linear map f : E �� E, if w = u + iv is an eigenvector of fC associated with the eigenvalue z = �� + i�� (where u, v �� E and ��, �� �� R), if �� =0 (i.e., z is not real) then (u, v�� =0 and (u, u�� = (v, v��, which implies that u and v are linearly independent, and if W is the subspace spanned by u and v, then f(W )= W and f.(W )= W . Furthermore, with respect to the (orthogonal) basis (u, v), the restriction of f to W has the matrix 
�˦� 
.
.�̦� 
If �� =0, then �� is a real eigenvalue of f, and either u or v is an eigenvector of f for ��. If W is the subspace spanned by u if u =0, or spanned by v =0 if u =0, then f(W ) . W and f.(W ) . W . 
Proof. Since w = u + iv is an eigenvector of fC, by de.nition it is nonnull, and either u =0 or v = 0. Proposition 16.10 implies that u . iv is an eigenvector of fC for �� . i��. It is easy to check that fC is normal. However, if �� = 0, then �� + i�� = �� . i��, and from Proposition 
16.4, the vectors u + iv and u . iv are orthogonal w.r.t. (.,.��C, that is, 
(u + iv, u . iv��C = (u, u��.( v, v�� +2i(u, v�� =0. 
Thus we get (u, v�� = 0 and (u, u�� = (v, v��, and since u = 0 or v = 0, u and v are linearly independent. Since 
f(u)= ��u . ��v and f(v)= ��u + ��v 
and since by Proposition 16.3 u + iv is an eigenvector of fC . for �� . i��, we have 
f . (u)= ��u + ��v and f . (v)= .��u + ��v, 
and thus f(W )= W and f.(W )= W , where W is the subspace spanned by u and v. When �� = 0, we have f(u)= ��u and f(v)= ��v, 
and since u = 0 or v = 0, either u or v is an eigenvector of f for ��. If W is the subspace spanned by u if u = 0, or spanned by v if u = 0, it is obvious that f(W ) . W and f.(W ) . W . Note that �� = 0 is possible, and this is why . cannot be replaced by =. 
The beginning of the proof of Proposition 16.11 actually shows that for every linear map 
f : E �� E there is some subspace W such that f(W ) . W , where W has dimension 1 or 
2. In general, it doesn��t seem possible to prove that W �� is invariant under f. However, this 
happens when f is normal. We can .nally prove our .rst main theorem. 
Theorem 16.12. (Main spectral theorem) Given a Euclidean space E of dimension n, for every normal linear map f : E �� E, there is an orthonormal basis (e1,...,en) such that the matrix of f w.r.t. this basis is a block diagonal matrix of the form 
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
Aj  =  ��j ��j  .��j ��j  ,  
where ��j, ��j �� R, with ��j > 0.  

Proof. We proceed by induction on the dimension n of E as follows. If n = 1, the result is trivial. Assume now that n �� 2. First, since C is algebraically closed (i.e., every polynomial has a root in C), the linear map fC: EC �� EC has some eigenvalue z = �� + i�� (where ��, �� �� R). Let w = u + iv be some eigenvector of fC for �� + i�� (where u, v �� E). We can now apply Proposition 16.11. 
If �� = 0, then either u or v is an eigenvector of f for �� �� R. Let W be the subspace of dimension 1 spanned by e1 = u/��u�� if u = 0, or by e1 = �� otherwise. It is obvious 
v/��vthat f(W ) . W and f.(W ) . WW. The orthogonal W �� of W has dimension n . 1, and by Proposition 16.9, we have fW �� . W �� . But the restriction of f to W �� is also normal, 
J 
and we conclude by applying the induction hypothesis to W �� . 
If �� = 0, then (u, v�� = 0 and (u, u�� = (v, v��, and if W is the subspace spanned by u/��u�� and v/��v��, then f(W )= W and f.(W )= W . We also know that the restriction of f to W has the matrix 
�˦� .�̦� 
with respect to the basis (u/��u��, v/��v��). If ��< 0, we let ��1 = ��, ��1 = .��, e1 = u/��u��, and e2 = v/��v��. If ��> 0, we let ��1 = ��, ��1 = ��, e1 = v/��v��, and e2 = u/��u��. In all cases, it is easily veri.ed that the matrix of the restriction of f to W w.r.t. the orthonormal basis (e1,e2) is 
��1 .��1
A1 = , 
��1 ��1 
W
J
where ��1,��1 �� R, with ��1 > 0. However, W �� has dimension n . 2, and by Proposition 16.9, fW �� . W �� . Since the restriction of f to W �� is also normal, we conclude by applying the induction hypothesis to W �� . 
After this relatively hard work, we can easily obtain some nice normal forms for the matrices of self-adjoint, skew-self-adjoint, and orthogonal linear maps. However, for the sake of completeness (and since we have all the tools to so do), we go back to the case of a Hermitian space and show that normal linear maps can be diagonalized with respect to an orthonormal basis. The proof is a slight generalization of the proof of Theorem 16.6. 

16.4. SELF-ADJOINT AND OTHER SPECIAL LINEAR MAPS 
Theorem 16.13. (Spectral theorem for normal linear maps on a Hermitian space) Given a Hermitian space E of dimension n, for every normal linear map f : E �� E there is an orthonormal basis (e1,...,en) of eigenvectors of f such that the matrix of f w.r.t. this basis 
is a diagonal matrix

.
. 
.... 

��1 ... ��2 ... 
.. .
.
.. ..
.
.. . 
... ��n 
....

, 

where ��j �� C. 
Proof. We proceed by induction on the dimension n of E as follows. If n = 1, the result is trivial. Assume now that n �� 2. Since C is algebraically closed (i.e., every polynomial has a root in C), the linear map f : E �� E has some eigenvalue �� �� C, and let w be some unit eigenvector for ��. Let W be the subspace of dimension 1 spanned by w. Clearly, f(W ) . W . By Proposition 16.3, w is an eigenvector of f. for ��, and thus f.(W ) . W . By Proposition 16.9, we also have f(W ��) . W �� . The restriction of f to W �� is still normal, and we conclude by applying the induction hypothesis to W �� (whose dimension is n . 1). 
Theorem 16.13 implies that (complex) self-adjoint, skew-self-adjoint, and orthogonal lin-ear maps can be diagonalized with respect to an orthonormal basis of eigenvectors. In this latter case, though, an orthogonal map is called a unitary map. Proposition 16.5 also shows that the eigenvalues of a self-adjoint linear map are real, and Proposition 16.7 shows that the eigenvalues of a skew self-adjoint map are pure imaginary or zero, and that the eigenvalues of a unitary map have absolute value 1. 
Remark: There is a converse to Theorem 16.13, namely, if there is an orthonormal basis (e1,...,en) of eigenvectors of f, then f is normal. We leave the easy proof as an exercise. 
In the next section we specialize Theorem 16.12 to self-adjoint, skew-self-adjoint, and orthogonal linear maps. Due to the additional structure, we obtain more precise normal forms. 


16.4 	Self-Adjoint, Skew-Self-Adjoint, and Orthogonal Linear Maps 
We begin with self-adjoint maps. 
Theorem 16.14. Given a Euclidean space E of dimension n, for every self-adjoint linear map f : E �� E, there is an orthonormal basis (e1,...,en) of eigenvectors of f such that the matrix of f w.r.t. this basis is a diagonal matrix 
.
. 
.... 

��1 ... ��2 ... 
,
..
. 
... ��n 
....
... 
... 
... 

where ��i �� R. 
Proof. We already proved this; see Theorem 16.8. However, it is instructive to give a more direct method not involving the complexi.cation of (.,.�� and Proposition 16.5. Since C is algebraically closed, fC has some eigenvalue �� + i��, and let u + iv be some eigenvector of fC for ��+i��, where ��, �� �� R and u, v �� E. We saw in the proof of Proposition 
16.10 that f(u)= ��u . ��v and f(v)= ��u + ��v. 
Since f = f. , (f(u),v�� = (u, f(v)�� 
for all u, v �� E. Applying this to 
f(u)= ��u . ��v and f(v)= ��u + ��v, 
we get (f(u),v�� = (��u . ��v, v�� = ��(u, v��. ��(v, v�� 
and (u, f(v)�� = (u, ��u + ��v�� = ��(u, u�� + ��(u, v��, 
and thus we get ��(u, v��. ��(v, v�� = ��(u, u�� + ��(u, v��, 
that is, ��((u, u�� + (v, v��)=0, 
which implies �� = 0, since either u =0 or v = 0. Therefore, �� is a real eigenvalue of f. Now going back to the proof of Theorem 16.12, only the case where �� = 0 applies, and the induction shows that all the blocks are one-dimensional. 
Theorem 16.14 implies that if ��1,...,��p are the distinct real eigenvalues of f, and Ei is the eigenspace associated with ��i, then 
E = E1 ���� ������ Ep, 
where Ei and Ej are orthogonal for all i = j. 
16.4. SELF-ADJOINT AND OTHER SPECIAL LINEAR MAPS 
Remark: Another way to prove that a self-adjoint map has a real eigenvalue is to use a little bit of calculus. We learned such a proof from Herman Gluck. The idea is to consider the real-valued function ��: E �� R de.ned such that 
��(u)= (f(u),u�� 
for every u �� E. This function is C�� , and if we represent f by a matrix A over some orthonormal basis, it is easy to compute the gradient vector
�ɦ�(X)= .��(X),..., .��(X)
.x1 .xn of �� at X. Indeed, we .nd that�ɦ�(X)=(A + AT)X, where X is a column vector of size n. But since f is self-adjoint, A = AT, and thus�ɦ�(X)=2AX. The next step is to .nd the maximum of the function �� on the sphere Sn.1 = {(x1,...,xn) �� Rn | x 21 + ������ + x 2 =1}.
n 
Since Sn.1 is compact and �� is continuous, and in fact C��, �� takes a maximum at some X on Sn.1 . But then it is well known that at an extremum X of �� we must have 
d��X (Y )=(�ɦ�(X),Y�� =0 for all tangent vectors Y to Sn.1 at X, and so�ɦ�(X) is orthogonal to the tangent plane at X, which means that
�ɦ�(X)= ��X for some �� �� R. Since�ɦ�(X)=2AX, we get 2AX = ��X, and thus ��/2 is a real eigenvalue of A (i.e., of f). Next we consider skew-self-adjoint maps. 
Theorem 16.15. Given a Euclidean space E of dimension n, for every skew-self-adjoint linear map f : E �� E there is an orthonormal basis (e1,...,en) such that the matrix of f 
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

such that each block Aj is either 0 or a two-dimensional matrix of the form 0 .��j
Aj = , 
��j 0 
where ��j �� R, with ��j > 0. In particular, the eigenvalues of fC are pure imaginary of the form ��i��j or 0. Proof. The case where n = 1 is trivial. As in the proof of Theorem 16.12, fC has some 
eigenvalue z = �� + i��, where ��, �� �� R. We claim that �� = 0. First we show that (f(w),w�� =0 for all w �� E. Indeed, since f = .f. , we get (f(w),w�� = (w, f . (w)�� = (w, .f(w)�� = .(w, f(w)�� = .(f(w),w��, since (.,.�� is symmetric. This implies that (f(w),w�� =0. Applying this to u and v and using the fact that f(u)= ��u . ��v and f(v)= ��u + ��v, we get 0= (f(u),u�� = (��u . ��v, u�� = ��(u, u��. ��(u, v�� and 0= (f(v),v�� = (��u + ��v, v�� = ��(u, v�� + ��(v, v��, from which, by addition, we get ��((v, v�� + (v, v��)=0. Since u =0 or v = 0, we have �� = 0. Then going back to the proof of Theorem 16.12, unless �� = 0, the case where u and v are orthogonal and span a subspace of dimension 2 applies, and the induction shows that all the blocks are two-dimensional or reduced to 0. 
Remark: One will note that if f is skew-self-adjoint, then ifC is self-adjoint w.r.t. (.,.��C. By Proposition 16.5, the map ifC has real eigenvalues, which implies that the eigenvalues of fC are pure imaginary or 0. 
Finally we consider orthogonal linear maps. 

16.4. SELF-ADJOINT AND OTHER SPECIAL LINEAR MAPS 
Theorem 16.16. Given a Euclidean space E of dimension n, for every orthogonal linear map f : E �� E there is an orthonormal basis (e1,...,en) such that the matrix of f w.r.t. 
such that each block Aj is either 1, .1, or a two-dimensional matrix of the form cos ��j . sin ��j
Aj = 
sin ��j cos ��j 
where 0 <��j <��. In particular, the eigenvalues of fC are of the form cos ��j �� i sin ��j, 1, or .1. 
Proof. The case where n = 1 is trivial. It is immediately veri.ed that f . f. = f. . f = id implies that fC . fC . = fC . . fC = id, so the map fC is unitary. By Proposition 16.7, the eigenvalues of fC have absolute value 1. As a consequence, the eigenvalues of fC are of the form cos �� �� i sin ��, 1, or .1. The theorem then follows immediately from Theorem 16.12, where the condition ��> 0 implies that sin ��j > 0, and thus, 0 <��j <��. 
It is obvious that we can reorder the orthonormal basis of eigenvectors given by Theorem 
where each block Aj is a two-dimensional rotation matrix Aj = ��I2 of the form cos ��j . sin ��j
Aj = 
sin ��j cos ��j 
with 0 <��j <��. 
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
cos ��j . sin ��j
Aj = 
sin ��j cos ��j 
with 0 <��j �� ��. Theorem 16.16 can be used to prove a version of the Cartan�CDieudonn��e theorem. 
Theorem 16.17. Let E be a Euclidean space of dimension n �� 2. For every isometry f �� O(E), if p = dim(E(1,f)) = dim(Ker (f . id)), then f is the composition of n . p re.ections, and n . p is minimal. 
Proof. From Theorem 16.16 there are r subspaces F1,...,Fr, each of dimension 2, such that 
E = E(1,f) �� E(.1,f) �� F1 ���� ������ Fr, 
and all the summands are pairwise orthogonal. Furthermore, the restriction ri of f to each 
H
Fi is a rotation ri = ��id. Each 2D rotation ri can be written as the composition ri = si . si 
H
of two re.ections si and si about lines in Fi (forming an angle ��i/2). We can extend si and 
H
si to hyperplane re.ections in E by making them the identity on Fi �� . Then 
HH
s.�� ����. s
r . sr 1 . s1 
agrees with f on F1 �� �� ���� �� Fr and is the identity on E(1,f) �� E(.1,f). If E(.1,f) 
HH
has an orthonormal basis of eigenvectors (v1,...,vq), letting sbe the re.ection about the 
j 
hyperplane (vj)�� , it is clear that HH HH 
s .�� ����. s
q 1 
agrees with f on E(.1,f) and is the identity on E(1,f) �� F1 ���� ������ Fr. But then 
HHHHH H
f = sq .�� ����. s1 . sr . sr .�� ����. s1 . s1, 
the composition of 2r + q = n . p re.ections. If f = st .�� ����. s1, 
for t re.ections si, it is clear that 
t
 
F = E(1,si) . E(1,f), 
i=1 
where E(1,si) is the hyperplane de.ning the re.ection si. By the Grassmann relation, if we intersect t �� n hyperplanes, the dimension of their intersection is at least n . t. Thus, n . t �� p, that is, t �� n . p, and n . p is the smallest number of re.ections composing f. 
As a corollary of Theorem 16.17, we obtain the following fact: If the dimension n of the Euclidean space E is odd, then every rotation f �� SO(E) admits 1 as an eigenvalue. 

16.5. NORMAL AND OTHER SPECIAL MATRICES 
Proof. The characteristic polynomial det(XI . f) of f has odd degree n and has real coef-.cients, so it must have some real root ��. Since f is an isometry, its n eigenvalues are of the form, +1, .1, and e��i��, with 0 <��<��, so �� = ��1. Now the eigenvalues e��i�� appear in conjugate pairs, and since n is odd, the number of real eigenvalues of f is odd. This implies that +1 is an eigenvalue of f, since otherwise .1 would be the only real eigenvalue of f, and since its multiplicity is odd, we would have det(f)= .1, contradicting the fact that f is a rotation. 
When n = 3, we obtain the result due to Euler which says that every 3D rotation R has an invariant axis D, and that restricted to the plane orthogonal to D, it is a 2D rotation. Furthermore, if (a, b, c) is a unit vector de.ning the axis D of the rotation R and if the angle of the rotation is ��, if B is the skew-symmetric matrix 
.  .  
0  .c  b  
B = . c  0  .a. ,  
.b  a  0  

then the Rodigues formula (Proposition 11.15) states that 
R = I + sin ��B + (1 . cos ��)B2 . 
The theorems of this section and of the previous section can be immediately translated in terms of matrices. The matrix versions of these theorems is often used in applications so we brie.y present them in the section. 


16.5 Normal and Other Special Matrices 
First we consider real matrices. Recall the following de.nitions. 
De.nition 16.3. Given a real m �� n matrix A, the transpose AT of A is the n �� m matrix AT =(aT ) de.ned such that 
ij
T 
a
ij = aji for all i, j,1 �� i �� m,1 �� j �� n. A real n �� n matrix A is 
. 	
normal if AAT = ATA, 

. 
symmetric if 


AT = A, 
. skew-symmetric if 
AT = .A, 
. orthogonal if 
AAT = ATA = In. 

Recall from Proposition 11.14 that when E is a Euclidean space and (e1,..., en) is an orthonormal basis for E, if A is the matrix of a linear map f : E �� E w.r.t. the basis (e1,...,en), then AT is the matrix of the adjoint f. of f. Consequently, a normal linear map has a normal matrix, a self-adjoint linear map has a symmetric matrix, a skew-self-adjoint linear map has a skew-symmetric matrix, and an orthogonal linear map has an orthogonal matrix. 
Furthermore, if (u1,...,un) is another orthonormal basis for E and P is the change of basis matrix whose columns are the components of the ui w.r.t. the basis (e1,...,en), then P is orthogonal, and for any linear map f : E �� E, if A is the matrix of f w.r.t (e1,...,en) and B is the matrix of f w.r.t. (u1,...,un), then 
B = P TAP. 
As a consequence, Theorems 16.12 and 16.14�C16.16 can be restated as follows. 
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
Dj =  ��j ��j  .��j ��j  ,  
where ��j, ��j �� R, with ��j > 0.  

Theorem 16.19. For every symmetric matrix A there is an orthogonal matrix P and a diagonal matrix D such that A = PDP T, where D is of the form 
.
. 
D = 

.... 

��1 ... ��2 ... 
.. .
.
.. ..
.
.. . 
... ��n 
....

, 

where ��i �� R. 
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

such that each block Dj is either 0 or a two-dimensional matrix of the form 0 .��j
Dj = , 
��j 0 
where ��j �� R, with ��j > 0. In particular, the eigenvalues of A are pure imaginary of the form ��i��j, or 0. 
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
cos ��j . sin ��j
Dj = 
sin ��j cos ��j 
where 0 <��j <��. In particular, the eigenvalues of A are of the form cos ��j �� i sin ��j, 1, or .1. 
Theorem 16.21 can be used to show that the exponential map exp: so(n) �� SO(n) is surjective; see Gallier [25]. 
We now consider complex matrices. 
WJ
De.nition 16.4. Given a complex m �� n matrix A, the transpose AT matrix AT = aT de.ned such that 
ij 
T 
a
ij = aji for all i, j,1 �� i �� m,1 �� j �� n. The conjugate A of A is the m �� n matrix A =(bij) de.ned such that bij = aij for all i, j,1 �� i �� m,1 �� j �� n. Given an m �� n complex matrix A, the adjoint A. of A is the matrix de.ned such that A . =(AT)=(A)T . A complex n �� n matrix A is 
of A is the n �� m 

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

Recall from Proposition 13.15 that when E is a Hermitian space and (e1,..., en) is an orthonormal basis for E, if A is the matrix of a linear map f : E �� E w.r.t. the basis (e1,...,en), then A. is the matrix of the adjoint f. of f. Consequently, a normal linear map has a normal matrix, a self-adjoint linear map has a Hermitian matrix, a skew-self-adjoint linear map has a skew-Hermitian matrix, and a unitary linear map has a unitary matrix. 
Furthermore, if (u1,...,un) is another orthonormal basis for E and P is the change of basis matrix whose columns are the components of the ui w.r.t. the basis (e1,...,en), then P is unitary, and for any linear map f : E �� E, if A is the matrix of f w.r.t (e1,...,en) and B is the matrix of f w.r.t. (u1,...,un), then 
B = P . AP. 
Theorem 16.13 and Proposition 16.7 can be restated in terms of matrices as follows. 
Theorem 16.22. For every complex normal matrix A there is a unitary matrix U and a diagonal matrix D such that A = UDU. . Furthermore, if A is Hermitian, then D is a real matrix; if A is skew-Hermitian, then the entries in D are pure imaginary or zero; and if A is unitary, then the entries in D have absolute value 1. 


16.6 	Rayleigh�CRitz Theorems and Eigenvalue Interlac-ing 
A fact that is used frequently in optimization problems is that the eigenvalues of a symmetric matrix are characterized in terms of what is known as the Rayleigh ratio, de.ned by 
xTAx 
R(A)(x)= T,x �� Rn ,x =0. 
xx 
16.6. RAYLEIGH�CRITZ THEOREMS AND EIGENVALUE INTERLACING 
The following proposition is often used to prove the correctness of various optimization or approximation problems (for example PCA; see Section 21.4). It is also used to prove Proposition 16.25, which is used to justify the correctness of a method for graph-drawing (see Chapter 19). 
Proposition 16.23. (Rayleigh�CRitz) If A is a symmetric n �� n matrix with eigenvalues ��1 �� ��2 �� �� ���� �� ��n and if (u1,...,un) is any orthonormal basis of eigenvectors of A, where ui is a unit eigenvector associated with ��i, then 
xTAx 
max = ��n 
xʱ xT
=0 x 
(with the maximum attained for x = un), and 
xTAx 
xʱ=0,x��{umax n.k+1,...,un}�� xTx = ��n.k 
(with the maximum attained for x = un.k), where 1 �� k �� n . 1. Equivalently, if Vk is the subspace spanned by (u1,...,uk), then 
xTAx ��k = xʱ max xT,k =1, . . . , n. 
=0,x��Vk x 
Proof. First observe that 
xTAx 
max = max {x TAx | x T x =1}, 
xʱ xTx
=0 x 
and similarly, 
xTAx ʱ max xT= max  x TAx | (x ��{un.k+1,...,un}��) �� (x T x = 1) . 
x=0,x��{un.k+1,...,un}�� x x
Since A is a symmetric matrix, its eigenvalues are real and it can be diagonalized with respect to an orthonormal basis of eigenvectors, so let (u1,...,un) be such a basis. If we write 
n
 
x = xiui, i=1 
a simple computation shows that 
n
 
x TAx = ��ixi 2 . i=1 
N 
If xTx = 1, then n x2 = 1, and since we assumed that ��1 �� ��2 �� �� ���� �� ��n, we get 
i=1 i 
nn 
x TAx = 2 �� ��n x 2 = ��n.
��ixii i=1 i=1 
Thus, max x TAx | x T x =1 �� ��n, 
x and since this maximum is achieved for en = (0, 0,..., 1), we conclude that 
max x TAx | x T x =1 = ��n. 
x 
Next observe that x ��{un.k+1,...,un}�� and xTx = 1 i. xn.k+1 = ������ = xn = 0 and Nn.k x2 = 1. Consequently, for such an x, we have 
i=1 i 
n.kn.k x TAx = ��ix 2 i �� ��n.k x 2 i = ��n.k. i=1 i=1 
Thus, ʱ T��T|��{}�ġ�Ax ()(=1) ��max xxu,...,uxx,..k+1knnnʱ 
x and since this maximum is achieved for en.k = (0,..., 0, 1, 0,..., 0) with a 1 in position 
. k,weconcludethat n ʱ 
max x TAx | (x ��{un.k+1,...,un}��) �� (x T x =1) = ��n.k, 
x 
ʱ 
as claimed. 
For our purposes we need the version of Proposition 16.23 applying to min instead of max, whose proof is obtained by a trivial modi.cation of the proof of Proposition 16.23. 
Proposition 16.24. (Rayleigh�CRitz) If A is a symmetric n �� n matrix with eigenvalues ��1 �� ��2 �� �� ���� �� ��n and if (u1,...,un) is any orthonormal basis of eigenvectors of A, where ui is a unit eigenvector associated with ��i, then 
xTAx 
min = ��1 
x=0 xTx 
(with the minimum attained for x = u1), and 
xTAx 
min T= ��i 
x=0,x��{u1,...,ui.1}�� xx 
(with the minimum attained for x = ui), where 2 �� i �� n. Equivalently, if Wk = Vk��.1 denotes the subspace spanned by (uk,...,un) (with V0 = (0)), then 
xTAx xTAx 
��k = min = min ,k =1, . . . , n. 
x=0,x��Wk xTx x=0,x��V �� xTx 
k.1 

16.6. RAYLEIGH�CRITZ THEOREMS AND EIGENVALUE INTERLACING 
Propositions 16.23 and 16.24 together are known the Rayleigh�CRitz theorem. 
As an application of Propositions 16.23 and 16.24, we prove a proposition which allows us to compare the eigenvalues of two symmetric matrices A and B = RTAR, where R is a rectangular matrix satisfying the equation RTR = I. 
De.nition 16.5. Given an n �� n symmetric matrix A and an m �� m symmetric B, with m �� n, if ��1 �� ��2 �� ������ �� ��n are the eigenvalues of A and ��1 �� ��2 �� ������ �� ��m are the eigenvalues of B, then we say that the eigenvalues of B interlace the eigenvalues of A if 
��i �� ��i �� ��n.m+i,i =1, . . . , m. 
For example, if n = 5 and m = 3, we have 
Firstweneedade.nition. �ܡܦ˦˦�113 (a)Theeigenvaluesof B interlacetheeigenvaluesof A. �ܡܡܡܡܡ�(b)If �˦˦�aretheeigenvaluesof A andarethe �������������̦̦�1212nm eigenvaluesof B,andif ��,thenthereisaneigenvector of B witheigenvalue = ��viisuchthat Rv isaneigenvectorof A witheigenvalue �˦�.iiProof. (a)Let()beanorthonormalbasisofeigenvectorsfor A,andlet()u,...,uv,...,v11nmbeanorthonormalbasisofeigenvectorsfor BLet Ubethesubspacespannedby()u,...,u. j 1jandlet Vbethesubspacespannedby().Forany i,thesubspace Vhasdimension v,...,vj 1ji T.i andthesubspace RUhasdimensionatmost i 1.Therefore,thereissomenonzero .i1 ��T��wehave Rv (U)ByProposition16.24andusingthefactthat RRI,wehave =..i1ʱ 
��2 �� ��2 �� ��4 
��3 �� ��3 �� ��5. 
Proposition 16.25. Let A be an n �� n symmetric matrix, R be an n �� m matrix such that RTR = I (with m �� n), and let B = RTAR (an m �� m matrix). The following properties hold: 
ʱ 
vector v �� Vi �� (RTUi.1)��, and since 
v TRT uj =(Rv)T uj =0,j =1,...,i . 1, 
(Rv)TARv vTBv 
��i �� = . 
(Rv)TRv vTv 
On the other hand, by Proposition 16.23, 
xTBx xTBx 
��i == , 
x=0,x��{vi+1,max...,vn}�� xTx x=0,x��{v1,max...,vi} xTx 
so 
wTBw 
�� ��i for all w �� Vi, 
wTw and since v �� Vi, we have vTBv 
��i �� T�� ��i,i =1, . . . , m. 
vv We can apply the same argument to the symmetric matrices .A and .B, to conclude that 
.��n.m+i ��.��i, 
that is, ��i �� ��n.m+i,i =1, . . . , m. Therefore, ��i �� ��i �� ��n.m+i,i =1, . . . , m, as desired. 
(b) If ��i = ��i, then 
(Rv)TARv vTBv 

��i = == ��i,
(Rv)TRv vTv so v must be an eigenvector for B and Rv must be an eigenvector for A, both for the eigenvalue ��i = ��i. 
Proposition 16.25 immediately implies the Poincar��e separation theorem. It can be used in situations, such as in quantum mechanics, where one has information about the inner 
T
products ui Auj. 
Proposition 16.26. (Poincar��e separation theorem) Let A be a n �� n symmetric (or Her-mitian) matrix, let r be some integer with 1 �� r �� n, and let (u1,...,ur) be r orthonormal 
T
vectors. Let B =(ui Auj) (an r �� r matrix), let ��1(A) �� ... �� ��n(A) be the eigenvalues of A and ��1(B) �� ... �� ��r(B) be the eigenvalues of B; then we have ��k(A) �� ��k(B) �� ��k+n.r(A),k =1, . . . , r. Observe that Proposition 16.25 implies that 
��1 + ������ + ��m �� tr(RTAR) �� ��n.m+1 + ������ + ��n. If P1 is the the n �� (n . 1) matrix obtained from the identity matrix by dropping its last column, we have P1 TP1 = I, and the matrix B = P1 TAP1 is the matrix obtained from A by deleting its last row and its last column. In this case the interlacing result is 
��1 �� ��1 �� ��2 �� ��2 �� �� ���� �� ��n.2 �� ��n.1 �� ��n.1 �� ��n, 
a genuine interlacing. We obtain similar results with the matrix Pn.r obtained by dropping the last n . r columns of the identity matrix and setting B = PnT.rAPn.r (B is the r �� r matrix obtained from A by deleting its last n . r rows and columns). In this case we have the following interlacing inequalities known as Cauchy interlacing theorem: 
��k �� ��k �� ��k+n.r,k =1, . . . , r. (.) 


sults 
xTAx 
16.7.THECOURANT�CFISCHERTHEOREM;PERTURBATIONRESULTS 16.7TheCourant�CFischerTheorem;PerturbationRe-AnotherusefultooltoproveeigenvalueequalitiesistheCourant�CFischercharacterizationof theeigenvaluesofasymmetricmatrix,alsoknownastheMin-max(andMax-min)theorem. Theorem16.27. (Courant�CFischer) Let A beasymmetric matrixwitheigenvalues ��nn Rn�ܡܡ� ���� V�˦˦�Ifdenotesthesetofsubspacesof ofdimension k,then �� .12 kn��min=maxk T��VW ��W,x=0 xxx.k+1nTAxx��min=max .k byProposition16.23,wehave �� min��=maxmax .k TT��Vʱ��W ��ʱ=0VW,x=0xxxxx,xxkk Therefore,weneedtoprovethereverseinequality;thatis,wehavetoshowthat �͡͡�V�ɡʡ�Nowforany W ,ifwecanprovethat WV =(0),thenforanynonzero WVvk�ܡܦ�min=max .k ʱ �͡͡ɡ�.Itremainstoprovethatdim(WV )1.However,dim(V)= k 1,sodim(V )= .k1..k1k1. k +1,andbyhypothesisdim(W )= kBytheGrassmannrelation, n . �͡͡͡�dim(W )+dim(V )=dim(WV )+dim(WV )+ , �� Rn��andsincedim(WV )dim()= +,weget n.k1eigenvaluesofasymmetricmatrixduetoHermannWeyl. 
T
W ��Vk x��W,x=0 xx 
ʱ 
ʱ Proof. Letusconsiderthesecondequality,theproofofthe.rstequalitybeingsimilar.Let ()beanyorthonormalbasisofeigenvectorsof A,where isauniteigenvector u,...,uu1in
associated with ��i. Observe that the space Vk spanned by (u1,...,uk) has dimension k, and xTAx xTAx 
ʱ 
ʱ 

xTAx 
��k �� max , for all W ��Vk.
T
x=0,x��W xx k.1 k.1, 
by Proposition 16.24 , we have 
xTAx vTAv xTAx x=0,x��V �� xTxvTv x��W,x=0 xTx 
k.1 
k.1k.1k.1
k + n . k +1 �� dim(W �� Vk��.1)+ n; that is, 1 �� dim(W �� Vk��.1), as claimed. The Courant�CFischer theorem yields the following useful result about perturbing the 
Proposition 16.28. Given two n ��n symmetric matrices A and B = A +��A, if ��1 �� ��2 �� ������ �� ��n are the eigenvalues of A and ��1 �� ��2 �� �� ���� �� ��n are the eigenvalues of B, then 
|��k . ��k|�� ��(��A)�ܽ� ��A��2 ,k =1, . . . , n. 
Proof. Let Vk be de.ned as in the Courant�CFischer theorem and let Vk be the subspace spanned by the k eigenvectors associated with ��1,...,��k. By the Courant�CFischer theorem applied to B, we have 
xTBx 
��k = min max 
T
W ��Vk x��W,x=0 xx xTBx 
�� max 
T
x��Vk xx xTAx xT�� Ax 
= max + 
TT
x��Vk xx xx xTAx xT��Ax 
�� max + max . 
ʱ 
x��Vk xTx x��Vk xTx 
By Proposition 16.23, we have 
xTAx 
��k = max ,
T
x��Vk xx 
so we obtain xTAx xT��Ax 
��k �� max + max 
TT
x��Vk xx x��Vk xx xT��Ax 
= ��k + max 
x��Vk xTx xT��Ax 
�� ��k + max .
T
x��Rn xx 
Now by Proposition 16.23 and Proposition 8.9, we have 
max xT��TAx = max ��i(��A) �� ��(��A)�ܽ� ��A��2 , 
x��Rn xx i where ��i(��A) denotes the ith eigenvalue of ��A, which implies that ��k �� ��k + ��(��A) �� ��k +�Ц�A��2 . By exchanging the roles of A and B, we also have ��k �� ��k + ��(��A) �� ��k +�Ц�A��2 , and thus, |��k . ��k|�� ��(��A)�ܽ� ��A��2 ,k =1, . . . , n, as claimed. 
n 
(��k . ��k)2�ܽ� ��A��F 2 , 
k=1 
proof; see Lax [44]. 
ʱ 
ʱ 
1. If i + j = k +1, then 
��i(A)+ ��j(B) �� ��k(A + B). 


2. If i + j = k + n, then 


16.7.THECOURANT�CFISCHERTHEOREM;PERTURBATIONRESULTS Proposition16.28alsoholdsforHermitianmatrices. AprettyresultofWielandtandHo.manassertsthat �н� whereistheFrobeniusnorm.However,theproofissigni.cantlyharderthantheabove F TheCourant�CFischertheoremcanalsobeusedtoprovesomefamousinequalitiesdueto HermannWeyl.Thesecanalsobeviewedasperturbationresults.Giventwosymmetric(or Hermitian)matrices A and B,let ��(A),��(B),and ��(AB)denotethe itheigenvalueof +iiiA,B,and AB,respectively,arrangedinnondecreasingorder. + Proposition16.29. (Weyl)Giventwosymmetric(orHermitian) matrices A and B��nn , �ܡ�thefollowinginequalitieshold:Forall i,j,k with 1 i,j,k :n�ܦ�(AB) ��(A)+ ��(B)+ .kij��(AB)= min+ .k��H,x=0 xxxTTAxBx xx��(A)= ��(B)= maxmax,.ijʱ�ʡ�F,x=0 G,x=0xxxxxx�ɡɡ�.��dim(FGH)=dim(F )+dim(GH)dim(F +(GH)) ��.��.dim(GH)=dim(G)+dim(H)dim(GH)dim(G)+dim(H)+ n, 
Proof. Observe that the .rst set of inequalities is obtained form the second set by replacing A by .A and B by .B, so it is enough to prove the second set of inequalities. By the Courant�CFischer theorem, there is a subspace H of dimension n . k + 1 such that 
xT(A + B)x T
Similarly, there exists a subspace F of dimension i and a subspace G of dimension j such that 
TT
We claim that F �� G �� H = (0). To prove this, we use the Grassmann relation twice. First, 
�� dim(F ) + dim(G �� H) . n, 
and second, 
so 
dim(F �� G �� H) �� dim(F ) + dim(G) + dim(H) . 2n. 

However, dim(F ) + dim(G) + dim(H)= i + j + n . k +1 and i + j = k + n, so we have dim(F �� G �� H) �� i + j + n . k +1 . 2n = k + n + n . k +1 . 2n =1, which shows that F �� G �� H = (0). Then for any unit vector z �� F �� G �� H = (0), we have 
��k(A + B) �� z T(A + B)z, ��i(A) �� z TAz, ��j(B) �� z TBz, 
establishing the desired inequality ��k(A + B) �� ��i(A)+ ��j(B). 
In the special case i = j = k, we obtain 
��1(A)+ ��1(B) �� ��1(A + B),��n(A + B) �� ��n(A)+ ��n(B). It follows that ��1 (as a function) is concave, while ��n (as a function) is convex. If i = 1 and j = k, we obtain 
��1(A)+ ��k(B) �� ��k(A + B), and if i = k and j = n, we obtain 
��k(A + B) �� ��k(A)+ ��n(B), 
and combining them, we get 
��1(A)+ ��k(B) �� ��k(A + B) �� ��k(A)+ ��n(B). 
In particular, if B is positive semide.nite, since its eigenvalues are nonnegative, we obtain the following inequality known as the monotonicity theorem for symmetric (or Hermitian) matrices: if A and B are symmetric (or Hermitian) and B is positive semide.nite, then 
��k(A) �� ��k(A + B) k =1, . . . , n. 
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

. 	Every normal linear map on a Euclidean space can be block diagonalized (blocks of size at most 2 �� 2) with respect to an orthonormal basis of eigenvectors. 

. 	Every normal linear map on a Hermitian space can be diagonalized with respect to an orthonormal basis of eigenvectors. 

. 	The spectral theorems for self-adjoint, skew-self-adjoint, and orthogonal linear maps (on a Euclidean space). 

. 	The spectral theorems for normal, symmetric, skew-symmetric, and orthogonal (real) matrices. 

. 	The spectral theorems for normal, Hermitian, skew-Hermitian, and unitary (complex) matrices. 

. 	The Rayleigh ratio and the Rayleigh�CRitz theorem. 

. 	Interlacing inequalities and the Cauchy interlacing theorem. 

. 	The Poincar��e separation theorem. 

. 	The Courant�CFischer theorem. 

. 	Inequalities involving perturbations of the eigenvalues of a symmetric matrix. 

. 	The Weyl inequalities. 


16.9 Problems 
Problem 16.1. Prove that the structure EC introduced in De.nition 16.2 is indeed a com-plex vector space. 
Problem 16.2. Prove that the formula 
(u1 + iv1,u2 + iv2��C = (u1,u2�� + (v1,v2�� + i((v1,u2��.( u1,v2��) 
de.nes a Hermitian form on EC that is positive de.nite and that (.,.��C agrees with (.,.��on real vectors. 
Problem 16.3. Given any linear map f : E �� E, prove the map fC . de.ned such that 
fC. (u + iv)= f . (u)+ if . (v) 
for all u, v �� E is the adjoint of fC w.r.t. (.,.��C. 
Problem 16.4. Let A be a real symmetric n �� n matrix whose eigenvalues are nonnegative. Prove that for every p> 0, there is a real symmetric matrix S whose eigenvalues are nonnegative such that Sp = A. 
Problem 16.5. Let A be a real symmetric n �� n matrix whose eigenvalues are positive. 
(1) Prove that there is a real symmetric matrix S such that A = eS . 
(2) Let S be a real symmetric n �� n matrix. Prove that A = eS is a real symmetric n �� n matrix whose eigenvalues are positive. 
Problem 16.6. Let A be a complex matrix. Prove that if A can be diagonalized with respect to an orthonormal basis, then A is normal. 
Problem 16.7. Let f : Cn �� Cn be a linear map. 
(1) Prove that if f is diagonalizable and if ��1,...,��n are the eigenvalues of f, then ��21,...,��2 n are the eigenvalues of f2, and if ��2 i = ��2 j implies that ��i = ��j, then f and f2 have the same eigenspaces. 

(2) Let f and g be two real self-adjoint linear maps f, g : Rn �� Rn . Prove that if f and g have nonnegative eigenvalues (f and g are positve semide.nite) and if f2 = g2, then f = g. 


Problem 16.8. (1) Let so(3) be the space of 3 �� 3 skew symmetric matrices 
.. 

.

..

    

a, b, c �� R 

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
A = .c  0  .a. �� so(3),  
.b  a  0  
��  

if we let �� = a2 + b2 + c2, recall from Section 11.7 (the Rodrigues formula) that the expo-nential map exp: so(3) �� SO(3) is given by 
A sin �� (1 . cos ��)
A2 
e = I3 + A + , if �� =0,
�Ȧ�2 
with exp(03)= I3. 
(2) Prove that eA is an orthogonal matrix of determinant +1, i.e., a rotation matrix. 
(3) Prove that the exponential map exp: so(3) �� SO(3) is surjective. For this proceed as follows: Pick any rotation matrix R �� SO(3); 
16.9. PROBLEMS 
(1) 
The case R = I is trivial. 

(2) 
If R = I and tr(R)= .1, then 


 	 
exp .1(R)=�� (R . RT ) 1+2cos �� = tr(R). 
2 sin �� 
(Recall that tr(R)= r11 + r22 + r33, the trace of the matrix R). 
Show that there is a unique skew-symmetric B with corresponding �� satisfying 0 < 
B
��<�� such that e= R. 
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
. 0 .dc . exp .1(R)= (2k + 1)�� . d 0 .b.,k �� Z ,
.	.
.cb 0 
where (b, c, d) is any unit vector such that for the corresponding skew symmetric matrix U, we have U2 = S. 
(4) To .nd a skew symmetric matrix U so that U2 = S = 12 (R . I) as in (3), we can solve the system 
.. 
b2 . 1 bc bd . bc c2 . 1 cd . = S. bd cdd2 . 1 
We immediately get b2,c2,d2, and then, since one of b, c, d is nonzero, say b, if we choose the positive square root of b2, we can determine c and d from bc and bd. 
Implement a computer program in Matlab to solve the above system. 
Problem 16.9. It was shown in Proposition 14.15 that the exponential map is a map exp: so(n) �� SO(n), where so(n) is the vector space of real n �� n skew-symmetric matrices. Use the spectral theorem to prove that the map exp: so(n) �� SO(n) is surjective. 
Problem 16.10. Let u(n) be the space of (complex) n �� n skew-Hermitian matrices (B. = .B) and let su(n) be its subspace consisting of skew-Hermitian matrice with zero trace (tr(B) = 0). 
(1) Prove that if B �� u(n), then eB �� U(n), and if if B �� su(n), then eB �� SU(n). Thus we have well-de.ned maps exp: u(n) �� U(n) and exp: su(n) �� SU(n). 
(2) 
Prove that the map exp: u(n) �� U(n) is surjective. 

(3) 
Prove that the map exp: su(n) �� SU(n) is surjective. 


Problem 16.11. Recall that a matrix B �� Mn(R) is skew-symmetric if BT = .B. Check that the set so(n) of skew-symmetric matrices is a vector space of dimension n(n . 1)/2, and thus is isomorphic to Rn(n.1)/2 . 
(1) Given a rotation matrix 
cos �� . sin �� 
R = ,
sin �� cos �� 
where 0 <��<��, prove that there is a skew symmetric matrix B such that 
R =(I . B)(I + B).1 . 
(2) Prove that the eigenvalues of a skew-symmetric matrix are either 0 or pure imaginary 
(that is, of the form i�� for �� �� R.). Let C : so(n) �� Mn(R) be the function (called the Cayley transform of B) given by 
C(B)=(I . B)(I + B).1 . 
Prove that if B is skew-symmetric, then I . B and I + B are invertible, and so C is well-de.ned. Prove that 
(I + B)(I . B)=(I . B)(I + B), 
and that (I + B)(I . B).1 =(I . B).1(I + B). 
16.9. PROBLEMS 
Prove that (C(B))TC(B)= I 
and that det C(B) = +1, 
so that C(B) is a rotation matrix. Furthermore, show that C(B) does not admit .1 as an eigenvalue. 
(3) Let SO(n) be the group of n �� n rotation matrices. Prove that the map 
C : so(n) �� SO(n) 
is bijective onto the subset of rotation matrices that do not admit .1 as an eigenvalue. Show that the inverse of this map is given by 
B =(I + R).1(I . R)=(I . R)(I + R).1 , 
where R �� SO(n) does not admit .1 as an eigenvalue. 
Problem 16.12. Please refer back to Problem 3.6. Let ��1,...,��n be the eigenvalues of A (not necessarily distinct). Using Schur��s theorem, A is similar to an upper triangular matrix B, that is, A = P BP .1 with B upper triangular, and we may assume that the diagonal entries of B in descending order are ��1,...,��n. 
(1) If the Eij are listed according to total order given by 
i = h and j>k (i, j) < (h, k) i. 
or i<h. 
prove that RB is an upper triangular matrix whose diagonal entries are 
(��n,...,��1,...,��n,...,��1),
�ˡ� �� 
2
n
and that LB is an upper triangular matrix whose diagonal entries are 
(��1,...,��1 ...,��n,...,��n).
�ˡ� �� �ˡ� �� 
nn 
Hint. Figure out what are RB(Eij)= EijB and LB(Eij)= BEij. 
(2) Use the fact that = R.1
LA = LP . LB . LP .1 ,RAP . RB . RP , 
to express adA = LA . RA in terms of LB . RB, and conclude that the eigenvalues of adA are ��i . ��j, for i =1,...,n, and for j = n, . . . , 1. 


