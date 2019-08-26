Chapter 5 
Direct Sums, Rank-Nullity Theorem, A.ne Maps 
In this chapter all vector spaces are de.ned over an arbitrary .eld K. For the sake of concreteness, the reader may safely assume that K = R. 
5.1 Direct Products 
There are some useful ways of forming new vector spaces from older ones. De.nition 5.1. Given p ◎ 2 vector spaces E1,...,Ep, the product F = E1 ≠· ··≠ Ep can be made into a vector space by de.ning addition and scalar multiplication as follows: (u1,...,up)+(v1,...,vp)=(u1 + v1,...,up + vp) ┡(u1,...,up)=(┡u1, . . . , ┡up), for all ui,vi Å Ei and all ┡ Å R. The zero vector of E1 ≠· ··≠ Ep is the p-tuple (0,..., 0),
   
p 
where the ith zero is the zero vector of Ei. With the above addition and multiplication, the vector space F = E1 ≠···≠ Ep is called the direct product of the vector spaces E1,...,Ep. 
As a special case, when E1 = ··· = Ep = R, we .nd again the vector space F = Rp. The projection maps pri : E1 ≠· ··≠ Ep ∪ Ei given by 
pri(u1,...,up)= ui are clearly linear. Similarly, the maps ini : Ei ∪ E1 ≠· ··≠ Ep given by ini(ui) = (0,..., 0,ui, 0,..., 0) 
125 

are injective and linear. If dim(Ei)= ni and if (ei 1,...,eini ) is a basis of Ei for i =1,...,p, then it is easy to see that the n1 + ··· + np vectors 
(e1 1, 0, . . . , 0),  . . . ,  (e1 n1 , 0, . . . , 0),  
.  .  .  
.  .  .  
.  .  .  
(0, . . . , 0, ei 1, 0, . . . , 0),  . . . ,  (0, . . . , 0, ei ni , 0, . . . , 0),  
.  .  .  
.  .  .  
.  .  .  
(0, . . . , 0, ep 1),  . . . ,  (0, . . . , 0, ep np )  

form a basis of E1 ≠· ··≠ Ep, and so 
dim(E1 ≠· ··≠ Ep) = dim(E1)+ ··· + dim(Ep). 

5.2 Sums and Direct Sums 
Let us now consider a vector space E and p subspaces U1,...,Up of E. We have a map 
a: U1 ≠· ··≠ Up ∪ E 
given by a(u1,...,up)= u1 + ··· + up, 
with ui Å Ui for i =1,...,p. It is clear that this map is linear, and so its image is a subspace of E denoted by 
U1 + ··· + Up 
and called the sum of the subspaces U1,...,Up. By de.nition, 
U1 + ··· + Up = {u1 + ··· + up | ui Å Ui, 1 ● i ● p}, 
and it is immediately veri.ed that U1 + ··· + Up is the smallest subspace of E containing U1,...,Up. This also implies that U1 + ··· + Up does not depend on the order of the factors Ui; in particular, 
U1 + U2 = U2 + U1. 
De.nition 5.2. For any vector space E and any p ◎ 2 subspaces U1,...,Up of E, if the map a: U1 ≠· ··≠ Up ∪ E de.ned above is injective, then the sum U1 + ··· + Up is called a direct sum and it is denoted by 
U1 쮵· ··쮵 Up. 
The space E is the direct sum of the subspaces Ui if 
E = U1 쮵· ··쮵 Up. 
5.2. SUMS AND DIRECT SUMS 
If the map a is injective, then by Proposition 2.16 we have Ker a = {(0,..., 0)} where 
p each 0 is the zero vector of E, which means that if ui Å Ui for i =1,...,p and if 
u1 + ··· + up =0, 
then (u1,...,up) = (0,..., 0), that is, u1 =0,...,up = 0. 
Proposition 5.1. If the map a: U1 ≠···≠ Up ∪ E is injective, then every u Å U1 + ··· + Up has a unique expression as a sum 
u = u1 + ··· + up, 
with ui Å Ui, for i =1,...,p. 
Proof. If u = v1 + ··· + vp = w1 + ··· + wp, 
with vi,wi Å Ui, for i =1,...,p, then we have 
w1 . v1 + ··· + wp . vp =0, 
and since vi,wi Å Ui and each Ui is a subspace, wi .vi Å Ui. The injectivity of a implies that wi .vi = 0, that is, wi = vi for i =1,...,p, which shows the uniqueness of the decomposition of u. 
Proposition 5.2. If the map a: U1 ≠···≠ Up ∪ E is injective, then any p nonzero vectors u1,...,up with ui Å Ui are linearly independent. 
Proof. To see this, assume that 
┡1u1 + ··· + ┡pup =0 
for some ┡i Å R. Since ui Å Ui and Ui is a subspace, ┡iui Å Ui, and the injectivity of a implies that ┡iui = 0, for i =1,...,p. Since ui = 0, we must have ┡i = 0 for i =1,...,p; that is, u1,...,up with ui Å Ui and ui= 0 are linearly independent.  
Observe that if a is injective, then we must have Ui ℃ Uj = (0) whenever i However, 
= j. this condition is generally not su.cient if p ◎ 3. For example, if E = R2 and U1 the line spanned by e1 = (1, 0), U2 is the line spanned by d = (1, 1), and U3 is the line spanned by e2 = (0, 1), then U1℃U2 = U1℃U3 = U2℃U3 = {(0, 0)}, but U1+U2 = U1+U3 = U2+U3 = R2 , so U1 + U2 + U3 is not a direct sum. For example, d is expressed in two di.erent ways as 
d = (1, 1) = (1, 0) + (0, 1) = e1 + e2. 
See Figure 5.1. 


As in the case of a sum, U1 쮵 U2 = U2 쮵 U1. Observe that when the map a is injective, then it is a linear isomorphism between U1 ≠· ··≠ Up and U1 쮵· ··쮵 Up. The di.erence is that U1 ≠· ··≠ Up is de.ned even if the spaces Ui are not assumed to be subspaces of some common space. 
If E is a direct sum E = U1 쮵···쮵Up, since any p nonzero vectors u1,...,up with ui Å Ui are linearly independent, if we pick a basis (uk)kÅIj in Uj for j =1,...,p, then (ui)iÅI with I = I1 ″· ··″ Ip is a basis of E. Intuitively, E is split into p independent subspaces. 
Conversely, given a basis (ui)iÅI of E, if we partition the index set I as I = I1 ″···″ Ip, then each subfamily (uk)kÅIj spans some subspace Uj of E, and it is immediately veri.ed that we have a direct sum 
E = U1 쮵· ··쮵 Up. 
De.nition 5.3. Let f : E ∪ E be a linear map. For any subspace U of E, if f(U) . U we say that U is invariant under f. 
Assume that E is .nite-dimensional, a direct sum E = U1 쮵· ··쮵 Up, and that each Uj is invariant under f. If we pick a basis (ui)iÅI as above with I = I1 ″ · ·· ″ Ip and with each (uk)kÅIj a basis of Uj, since each Uj is invariant under f, the image f(uk) of every basis vector uk with k Å Ij belongs to Uj, so the matrix A representing f over the basis (ui)iÅI is a block diagonal matrix of the form 
.
. 
A = 

.... 

A1 
A2 
....

,

..
. 
Ap 

5.2. SUMS AND DIRECT SUMS 
with each block Aj a dj ≠ dj-matrix with dj = dim(Uj) and all other entries equal to 0. If 
dj = 1 for j =1,...,p, the matrix A is a diagonal matrix. 
There are natural injections from each Ui to E denoted by ini : Ui ∪ E. 
Now, if p = 2, it is easy to determine the kernel of the map a: U1 ≠ U2 ∪ E. We have 
a(u1,u2)= u1 + u2 = 0 i. u1 = .u2,u1 Å U1,u2 Å U2, 
which implies that Ker a = {(u, .u) | u Å U1 ℃ U2}. 
Now, U1 ℃ U2 is a subspace of E and the linear map u ∪ (u, .u) is clearly an isomorphism between U1 ℃ U2 and Ker a, so Ker a is isomorphic to U1 ℃ U2. As a consequence, we get the following result: 
Proposition 5.3. Given any vector space E and any two subspaces U1 and U2, the sum U1 + U2 is a direct sum i. U1 ℃ U2 = (0). 
An interesting illustration of the notion of direct sum is the decomposition of a square matrix into its symmetric part and its skew-symmetric part. Recall that an n ≠ n matrix A Å Mn is symmetric if AT = A, skew -symmetric if AT = .A. It is clear that s 
S(n)= {A Å Mn | AT = A} and Skew(n)= {A Å Mn | AT = .A} 
are subspaces of Mn, and that S(n) ℃ Skew(n) = (0). Observe that for any matrix A Å Mn, the matrix H(A)=(A + AT)/2 is symmetric and the matrix S(A)=(A . AT)/2 is skew-symmetric. Since 
A + AT A . AT 
A = H(A)+ S(A)= + ,
22 we see that Mn = S(n)+ Skew(n), and since S(n) ℃ Skew(n) = (0), we have the direct sum 
Mn = S(n) 쮵 Skew(n). 
Remark: The vector space Skew(n) of skew-symmetric matrices is also denoted by so(n). It is the Lie algebra of the group SO(n). 
Proposition 5.3 can be generalized to any p ◎ 2 subspaces at the expense of notation. The proof of the following proposition is left as an exercise. 
Proposition 5.4. Given any vector space E and any p ◎ 2 subspaces U1,...,Up, the fol-lowing properties are equivalent: 
(1) 
The sum U1 + ··· + Up is a direct sum. 

(2) 
We have 


 p 
령 
Ui ℃Uj= (0),i =1, . . . , p. j=1,j 
=i 
(3) We have  Ui ℃  i.1령 j=1 Uj  = (0),  i = 2, . . . , p.  
Because of the isomorphism  

U1 ≠· ··≠ Up ≒ U1 쮵· ··쮵 Up, 
we have 
Proposition 5.5. If E is any vector space, for any (.nite-dimensional) subspaces U1,..., Up of E, we have dim(U1 쮵· ··쮵 Up) = dim(U1)+ ··· + dim(Up). If E is a direct sum E = U1 쮵· ··쮵 Up, since every u Å E can be written in a unique way as u = u1 + ··· + up with ui Å Ui for i =1 ...,p, we can de.ne the maps ┪i : E ∪ Ui, called projections, by ┪i(u)= ┪i(u1 + ··· + up)= ui. It is easy to check that these maps are linear and satisfy the following properties: 
 
┪i if i = j
┪j . ┪i =
0 if i = j, ┪1 + ··· + ┪p = idE. For example, in the case of the direct sum Mn = S(n) 쮵 Skew(n), the projection onto S(n) is given by A + AT 
┪1(A)= H(A)= ,
2 and the projection onto Skew(n) is given by A . AT 
┪2(A)= S(A)= . 
2 Clearly, H(A)+S(A)= A, H(H(A)) = H(A), S(S(A)) = S(A), and H(S(A)) = S(H(A)) = 
0. 
A function f such that f . f = f is said to be idempotent. Thus, the projections ┪i are idempotent. Conversely, the following proposition can be shown: 

5.3. THE RANK-NULLITY THEOREM; GRASSMANN’S RELATION 
Proposition 5.6. Let E be a vector space. For any p ◎ 2 linear maps fi : E ∪ E, if 
fi if i = jfj . fi = 
0 if i = j, f1 + ··· + fp = idE, then if we let Ui = fi(E), we have a direct sum E = U1 쮵· ··쮵 Up. 
We also have the following proposition characterizing idempotent linear maps whose proof is also left as an exercise. 
Proposition 5.7. For every vector space E, if f : E ∪ E is an idempotent linear map, i.e., f . f = f, then we have a direct sum E = Ker f 쮵 Im f, so that f is the projection onto its image Im f. We are now ready to prove a very crucial result relating the rank and the dimension of the kernel of a linear map. 


5.3 The Rank-Nullity Theorem; Grassmann’s Relation 
We begin with the following theorem which shows that given a linear map f : E ∪ F , its domain E is the direct sum of its kernel Ker f with some isomorphic copy of its image Im f. 
Theorem 5.8. (Rank-nullity theorem) Let f : E ∪ F be a linear map with .nite image. For any choice of a basis (f1,...,fr) of Im f, let (u1,...,ur) be any vectors in E such that fi = f(ui), for i =1,...,r. If s: Im f ∪ E is the unique linear map de.ned by s(fi)= ui, for i =1,...,r, then s is injective, f . s = id, and we have a direct sum 
E = Ker f 쮵 Im s 
as illustrated by the following diagram: 
f 
 Im f . F. 
Ker f ∪ E = Ker f 쮵 Im s 묏∪ 
s
See Figure 5.2. As a consequence, if E is .nite-dimensional, then dim(E) = dim(Ker f) + dim(Im f) = dim(Ker f) + rk(f). 

Proof. The vectors u1,...,ur must be linearly independent since otherwise we would have a nontrivial linear dependence 
┡1u1 + ··· + ┡rur =0, 
and by applying f, we would get the nontrivial linear dependence 
0= ┡1f(u1)+ ··· + ┡rf(ur)= ┡1f1 + ··· + ┡rfr, 
contradicting the fact that (f1,...,fr) is a basis. Therefore, the unique linear map s given by s(fi)= ui, for i =1,...,r, is a linear isomorphism between Im f and its image, the subspace spanned by (u1,...,ur). It is also clear by de.nition that f . s = id. For any u Å E, let 
h = u . (s . f)(u). 
Since f . s = id, we have 
f(h)= f(u . (s . f)(u)) = f(u) . (f . s . f)(u) = f(u) . (id . f)(u)= f(u) . f(u)=0, 
which shows that h Å Ker f. Since h = u . (s . f)(u), it follows that 
u = h + s(f(u)), 
with h Å Ker f and s(f(u)) Å Im s, which proves that 
E = Ker f + Im s. 

5.3. THE RANK-NULLITY THEOREM; GRASSMANN’S RELATION 
Now if u Å Ker f ℃ Im s, then u = s(v) for some v Å F and f(u) = 0 since u Å Ker f. Since u = s(v) and f . s = id, we get 
0= f(u)= f(s(v)) = v, 
and so u = s(v)= s(0) = 0. Thus, Ker f ℃ Im s = (0), which proves that we have a direct sum 
E = Ker f 쮵 Im s. 
The equation 
dim(E) = dim(Ker f) + dim(Im f) = dim(Ker f) + rk(f) 
is an immediate consequence of the fact that the dimension is an additive property for direct sums, that by de.nition the rank of f is the dimension of the image of f, and that dim(Im s) = dim(Im f), because s is an isomorphism between Im f and Im s. 
Remark: The statement E = Ker f 쮵 Im s holds if E has in.nite dimension. It still holds if Im(f) also has in.nite dimension. 
De.nition 5.4. The dimension dim(Ker f) of the kernel of a linear map f is called the nullity of f. 
We now derive some important results using Theorem 5.8. 
Proposition 5.9. Given a vector space E, if U and V are any two .nite-dimensional sub-spaces of E, then 
dim(U) + dim(V ) = dim(U + V ) + dim(U ℃ V ), 
an equation known as Grassmann’s relation. 
Proof. Recall that U + V is the image of the linear map 
a: U ≠ V ∪ E 
given by a(u, v)= u + v, 
and that we proved earlier that the kernel Ker a of a is isomorphic to U ℃ V . By Theorem 5.8, 
dim(U ≠ V ) = dim(Ker a) + dim(Im a), 
but dim(U ≠ V ) = dim(U) + dim(V ), dim(Ker a) = dim(U ℃ V ), and Im a = U + V , so the Grassmann relation holds. 
The Grassmann relation can be very useful to .gure out whether two subspace have a nontrivial intersection in spaces of dimension > 3. For example, it is easy to see that in R5 , there are subspaces U and V with dim(U) = 3 and dim(V ) = 2 such that U ℃ V = (0); for example, let U be generated by the vectors (1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), and V be generated by the vectors (0, 0, 0, 1, 0) and (0, 0, 0, 0, 1). However, we claim that if dim(U)=3 and dim(V ) = 3, then dim(U ℃ V ) ◎ 1. Indeed, by the Grassmann relation, we have 
dim(U) + dim(V ) = dim(U + V ) + dim(U ℃ V ), 
namely 3+3 = 6 = dim(U + V ) + dim(U ℃ V ), 
and since U + V is a subspace of R5, dim(U + V ) ● 5, which implies 
6 ● 5 + dim(U ℃ V ), 
that is 1 ● dim(U ℃ V ). 
As another consequence of Proposition 5.9, if U and V are two hyperplanes in a vector space of dimension n, so that dim(U)= n . 1 and dim(V )= n . 1, the reader should show that 
dim(U ℃ V ) ◎ n . 2, 
and so, if U = V , then dim(U ℃ V )= n . 2. 
Here is a characterization of direct sums that follows directly from Theorem 5.8. 
Proposition 5.10. If U1,...,Up are any subspaces of a .nite dimensional vector space E, then 
dim(U1 + ··· + Up) ● dim(U1)+ ··· + dim(Up), 
and 
dim(U1 + ··· + Up) = dim(U1)+ ··· + dim(Up) i. the Uis form a direct sum U1 쮵· ··쮵 Up. 
Proof. If we apply Theorem 5.8 to the linear map 
a: U1 ≠· ··≠ Up ∪ U1 + ··· + Up 
given by a(u1,...,up)= u1 + ··· + up, we get 
dim(U1 + ··· + Up) = dim(U1 ≠· ··≠ Up) . dim(Ker a) = dim(U1)+ ··· + dim(Up) . dim(Ker a), 
so the inequality follows. Since a is injective i. Ker a = (0), the Uis form a direct sum i. the second equation holds. 

5.3. THE RANK-NULLITY THEOREM; GRASSMANN’S RELATION 
Another important corollary of Theorem 5.8 is the following result: 
Proposition 5.11. Let E and F be two vector spaces with the same .nite dimension dim(E) = dim(F )= n. For every linear map f : E ∪ F , the following properties are equivalent: 
(a) 
f is bijective. 

(b) 
f is surjective. 

(c) 
f is injective. 

(d) 
Ker f = (0). 


Proof. Obviously, (a) implies (b). If f is surjective, then Im f = F , and so dim(Im f)= n. By Theorem 5.8, 
dim(E) = dim(Ker f) + dim(Im f), 
and since dim(E)= n and dim(Im f)= n, we get dim(Ker f) = 0, which means that 
Ker f = (0), and so f is injective (see Proposition 2.16). This proves that (b) implies (c). 
If f is injective, then by Proposition 2.16, Ker f = (0), so (c) implies (d). 
Finally, assume that Ker f = (0), so that dim(Ker f)=0 and f is injective (by Proposi-tion 2.16). By Theorem 5.8, 
dim(E) = dim(Ker f) + dim(Im f), 
and since dim(Ker f)=0, we get 
dim(Im f) = dim(E) = dim(F ), 
which proves that f is also surjective, and thus bijective. This proves that (d) implies (a) and concludes the proof. 
One should be warned that Proposition 5.11 fails in in.nite dimension. A linear map may be injective without being surjective and vice versa. 
Here are a few applications of Proposition 5.11. Let A be an n ≠ n matrix and assume that A some right inverse B, which means that B is an n ≠ n matrix such that 
AB = I. 
The linear map associated with A is surjective, since for every u Å Rn, we have A(Bu)= u. By Proposition 5.11, this map is bijective so B is actually the inverse of A; in particular BA = I. 
Similarly, assume that A has a left inverse B, so that 
BA = I. 

This time the linear map associated with A is injective, because if Au = 0, then BAu = B0 = 0, and since BA = I we get u = 0. Again, by Proposition 5.11, this map is bijective so B is actually the inverse of A; in particular AB = I. 
Now assume that the linear system Ax = b has some solution for every b. Then the linear map associated with A is surjective and by Proposition 5.11, A is invertible. 
Finally assume that the linear system Ax = b has at most one solution for every b. Then the linear map associated with A is injective and by Proposition 5.11, A is invertible. 
We also have the following basic proposition about injective or surjective linear maps. 
Proposition 5.12. Let E and F be vector spaces, and let f : E ∪ F be a linear map. If 
f : E ∪ F is injective, then there is a surjective linear map r : F ∪ E called a retraction, such that r . f = idE. See Figure 5.3. If f : E ∪ F is surjective, then there is an injective linear map s: F ∪ E called a section, such that f . s = idF . See Figure 5.2. 

Proof. Let (ui)iÅI be a basis of E. Since f : E ∪ F is an injective linear map, by Proposition 2.17, (f(ui))iÅI is linearly independent in F . By Theorem 2.9, there is a basis (vj)jÅJ of F , where I . J, and where vi = f(ui), for all i Å I. By Proposition 2.17, a linear map r : F ∪ E can be de.ned such that r(vi)= ui, for all i Å I, and r(vj)= w for all j Å (J . I), where w is any given vector in E, say w = 0. Since r(f(ui)) = ui for all i Å I, by Proposition 2.17, we have r . f = idE. 
Now assume that f : E ∪ F is surjective. Let (vj)jÅJ be a basis of F . Since f : E ∪ F is surjective, for every vj Å F , there is some uj Å E such that f(uj)= vj. Since (vj)jÅJ is a basis of F , by Proposition 2.17, there is a unique linear map s: F ∪ E such that s(vj)= uj. Also since f(s(vj)) = vj, by Proposition 2.17 (again), we must have f . s = idF . 

5.4. AFFINE MAPS 
Remark: Proposition 5.12 also holds if E or F has in.nite dimension. 
The converse of Proposition 5.12 is obvious. 
The notion of rank of a linear map or of a matrix important, both theoretically and practically, since it is the key to the solvability of linear equations. We have the following simple proposition. 
Proposition 5.13. Given a linear map f : E ∪ F , the following properties hold: 
(i) 
rk(f) + dim(Ker f) = dim(E). 

(ii) 
rk(f) ● min(dim(E), dim(F )). 


Proof. Property (i) follows from Proposition 5.8. As for (ii), since Im f is a subspace of F , we have rk(f) ● dim(F ), and since rk(f) +dim(Ker f) = dim(E), we have rk(f) ● dim(E). 
The rank of a matrix is de.ned as follows. 
De.nition 5.5. Given a m ≠ n-matrix A =(aij), the rank rk(A) of the matrix A is the maximum number of linearly independent columns of A (viewed as vectors in Rm). 
In view of Proposition 2.10, the rank of a matrix A is the dimension of the subspace of Rm generated by the columns of A. Let E and F be two vector spaces, and let (u1,...,un) be a basis of E, and (v1,...,vm) a basis of F . Let f : E ∪ F be a linear map, and let M(f) be its matrix w.r.t. the bases (u1,...,un) and (v1,...,vm). Since the rank rk(f) of f is the dimension of Im f, which is generated by (f(u1),...,f(un)), the rank of f is the maximum number of linearly independent vectors in (f(u1),...,f(un)), which is equal to the number of linearly independent columns of M(f), since F and Rm are isomorphic. Thus, we have rk(f) = rk(M(f)), for every matrix representing f. 
We will see later, using duality, that the rank of a matrix A is also equal to the maximal number of linearly independent rows of A. 
5.4 A.ne Maps 
We showed in Section 2.7 that every linear map f must send the zero vector to the zero vector; that is, 
f(0) = 0. 
Yet for any .xed nonzero vector u Å E (where E is any vector space), the function tu given by 
tu(x)= x + u, for all x Å E 
shows up in practice (for example, in robotics). Functions of this type are called translations. They are not linear for u = 0, since tu(0) = 0+ u = u. 
More generally, functions combining linear maps and translations occur naturally in many applications (robotics, computer vision, etc.), so it is necessary to understand some basic properties of these functions. For this, the notion of a.ne combination turns out to play a key role. 
Recall from Section 2.7 that for any vector space E, given any family (ui)iÅI of vectors ui Å E, an a.ne combination of the family (ui)iÅI is an expression of the form
령령 
┡iui with┡i =1, iÅIiÅI 
where (┡i)iÅI is a family of scalars. 
A linear combination places no restriction on the scalars involved, but an a.ne com-bination is a linear combination with the restriction that the scalars ┡i must add up to 1. Nevertheless, a linear combination can always be viewed as an a.ne combination using the following trick involving 0. For any family (ui)iÅI of vectors in E and for any family of 
w 
scalars (┡i)iÅI , we can write the linear combination iÅI ┡iui as an a.ne combination as follows:
령령 령 
┡iui =┡iui +1 .┡i 0. 
iÅIiÅIiÅI 
A.ne combinations are also called barycentric combinations. 
Although this is not obvious at .rst glance, the condition that the scalars ┡i add up to 1 ensures that a.ne combinations are preserved under translations. To make this precise, consider functions f : E ∪ F , where E and F are two vector spaces, such that there is some linear map h: E ∪ F and some .xed vector b Å F (a translation vector), such that 
f(x)= h(x)+ b, for all x Å E. 
The map f given by x1 8/5 .6/5 x1 1 
∪ + 
x23/10 2/5 x2 1 
is an example of the composition of a linear map with a translation. We claim that functions of this type preserve a.ne combinations. 
Proposition 5.14. For any two vector spaces E and F , given any function f : E ∪ F de.ned such that 
f(x)= h(x)+ b, for all x Å E, 
where h: E ∪ F is a linear map and b is some .xed vector in F , for every a.ne combination 
ww 
(with ┡i =1), we have 
iÅI ┡iuiiÅI 
령령 
f┡iui =┡if(ui). 
iÅIiÅI 
In other words, f preserves a.ne combinations. 


5.4. AFFINE MAPS 
w 
Proof. By de.nition of f, using the fact that h is linear and the fact that iÅI ┡i = 1, we 
have  
령  령  
f ┡iui  = h ┡iui  + b  
iÅ  iÅI 

령 
=┡ih(ui)+1b iÅI 
령령 
=┡ih(ui)+┡i b iÅIiÅI
령 
=┡i(h(ui)+ b) iÅI
령 
=┡if(ui), iÅI 
as claimed. 
w 
Observe how the fact that iÅI ┡i = 1 was used in a crucial way in Line 3. Surprisingly, the converse of Proposition 5.14 also holds. 
Proposition 5.15. For any two vector spaces E and F , let f : E ∪ F be any function that 
ww 
preserves a.ne combinations, i.e., for every a.ne combination (with = 
iÅI ┡iuiiÅI ┡i 
1), we have 
령령 
f┡iui =┡if(ui). iÅIiÅI 
Then for any a Å E, the function h: E ∪ F given by h(x)= f(a + x) . f(a) is a linear map independent of a, and f(a + x)= h(x)+ f(a), for all x Å E. In particular, for a =0, if we let c = f(0), then f(x)= h(x)+ c, for all x Å E. Proof. First, let us check that h is linear. Since f preserves a.ne combinations and since a + u + v =(a + u)+(a + v) . a is an a.ne combination (1 + 1 . 1 = 1), we have h(u + v)= f(a + u + v) . f(a) = f((a + u)+(a + v) . a) . f(a) = f(a + u)+ f(a + v) . f(a) . f(a) = f(a + u) . f(a)+ f(a + v) . f(a) = h(u)+ h(v). 
This proves that 
h(u + v)= h(u)+ h(v), u,v Å E. 
Observe that a + ┡u = ┡(a + u) + (1 . ┡)a is also an a.ne combination (┡ +1 . ┡ = 1), so we have 
h(┡u)= f(a + ┡u) . f(a) = f(┡(a + u) + (1 . ┡)a) . f(a) = ┡f(a + u) + (1 . ┡)f(a) . f(a) = ┡(f(a + u) . f(a)) = ┡h(u). 
This proves that h(┡u)= ┡h(u),u Å E, ┡ Å R. 
Therefore, h is indeed linear. For any b Å E, since b + u =(a + u) . a + b is an a.ne combination (1 . 1+1=1), we have 
f(b + u) . f(b)= f((a + u) . a + b) . f(b) = f(a + u) . f(a)+ f(b) . f(b) = f(a + u) . f(a), 
which proves that for all a, b Å E, 
f(b + u) . f(b)= f(a + u) . f(a),u Å E. 
Therefore h(x)= f(a + u) . f(a) does not depend on a, and it is obvious by the de.nition of h that f(a + x)= h(x)+ f(a), for all x Å E. 
For a = 0, we obtain the last part of our proposition. 
We should think of a as a chosen origin in E. The function f maps the origin a in E to the origin f(a) in F . Proposition 5.15 shows that the de.nition of h does not depend on the origin chosen in E. Also, since 
f(x)= h(x)+ c, for all x Å E 
for some .xed vector c Å F , we see that f is the composition of the linear map h with the translation tc (in F ). The unique linear map h as above is called the linear map associated with f, and it is 
.∪ 
sometimes denoted by f . In view of Propositions 5.14 and 5.15, it is natural to make the following de.nition. 

5.4. AFFINE MAPS 
De.nition 5.6. For any two vector spaces E and F , a function f : E ∪ F is an a.ne 
w 
map if f preserves a.ne combinations, i.e., for every a.ne combination iÅI ┡iui (with
w 
iÅI ┡i = 1), we have 
령령 
f┡iui =┡if(ui). 
iÅIiÅI 
Equivalently, a function f : E ∪ F is an a.ne map if there is some linear map h: E ∪ F 
.∪ 
(also denoted by f ) and some .xed vector c Å F such that 
f(x)= h(x)+ c, for all x Å E. 
Note that a linear map always maps the standard origin 0 in E to the standard origin 0 in F . However an a.ne map usually maps 0 to a nonzero vector c = f(0). This is the “translation component” of the a.ne map. 
When we deal with a.ne maps, it is often fruitful to think of the elements of E and F not only as vectors but also as points. In this point of view, points can only be combined using a.ne combinations, but vectors can be combined in an unrestricted fashion using linear combinations. We can also think of u + v as the result of translating the point u by the translation tv. These ideas lead to the de.nition of a.ne spaces. 
The idea is that instead of a single space E, an a.ne space consists of two sets E and 
.∪ .∪ 
E , where E is just an unstructured set of points, and E is a vector space. Furthermore, the 
.∪ .∪ 
vector space E acts on E. We can think of E as a set of translations speci.ed by vectors, 
.∪ 
and given any point a Å E and any vector (translation) u Å E , the result of translating a by u is the point (not vector) a + u. Formally, we have the following de.nition. 
De.nition 5.7. An a.ne space is either the degenerate space reduced to the empty set, or 
 .∪  .∪ 
a tripleE, E, +consisting of a nonempty set E (of points), a vector space E (of trans-
.∪ 
lations, or free vectors), and an action +: E ≠ E ∪ E, satisfying the following conditions. 
(A1) a +0 = a, for every a Å E. 
.∪ 
(A2) (a + u)+ v = a +(u + v), for every a Å E, and every u, v Å E . 
.∪ 
(A3) For any two points a, b Å E, there is a unique u Å E such that a + u = b. 
.
.∪ ∪ 
The unique vector u Å E such that a + u = b is denoted by ab, or sometimes by ab, or even by b . a. Thus, we also write 
.
∪ 
b = a + ab 
(or b = a + ab, or even b = a +(b . a)). 

It is important to note that adding or rescaling points does not make sense! However, 

.∪ 
using the fact that E acts on E is a special way (this action is transitive and faithful), it is possible to de.ne rigorously the notion of a.ne combinations of points and to de.ne a.ne spaces, a.ne maps, etc. However, this would lead us to far a.eld, and for our purposes it is enough to stick to vector spaces and we will not distinguish between vector addition + and translation of a point by a vector +. Still, one should be aware that a.ne combinations really apply to points, and that points are not vectors! 
If E and F are .nite dimensional vector spaces with dim(E)= n and dim(F )= m, then it is useful to represent an a.ne map with respect to bases in E in F . However, the translation part c of the a.ne map must be somehow incorporated. There is a standard trick to do this which amounts to viewing an a.ne map as a linear map between spaces of dimension n + 1 and m + 1. We also have the extra .exibility of choosing origins a Å E and b Å F . 
Let (u1,...,un) be a basis of E,(v1,...,vm) be a basis of F , and let a Å E and b Å F be any two .xed vectors viewed as origins. Our a.ne map f has the property that if v = f(u), then 
v . b = f(a + u . a) . b = f(a) . b + h(u . a), 
where the last equality made use of the fact that h(x)= f(a + x) . f(a). Ifwelet y = v . b, x = u . a, and d = f(a) . b, then 
y = h(x)+ d, x Å E. 
Over the basis U =(u1,...,un), we write 
x = x1u1 + ··· + xnun, 
and over the basis V =(v1,...,vm), we write 
y = y1v1 + ··· + ymvm, d = d1v1 + ··· + dmvm. 
Then since y = h(x)+ d, 
if we let A be the m ≠ n matrix representing the linear map h, that is, the jth column of A consists of the coordinates of h(uj) over the basis (v1,...,vm), then we can write 
yV = AxU + dV . 
where xU =(x1,...,xn)T , yV =(y1,...,ym)T, and dV =(d1,...,dm)T . The above is the ma-trix representation of our a.ne map f with respect to (a, (u1,...,un)) and (b, (v1,...,vm)). 
The reason for using the origins a and b is that it gives us more .exibility. In particular, we can choose b = f(a), and then f behaves like a linear map with respect to the origins a and b = f(a). 

5.4. AFFINE MAPS 
When E = F , if there is some a Å E such that f(a)= a (a is a .xed point of f), then we can pick b = a. Then because f(a)= a, we get 
v = f(u)= f(a + u . a)= f(a)+ h(u . a)= a + h(u . a), 
that is v . a = h(u . a). 
With respect to the new origin a, if we de.ne x and y by 
x = u . a y = v . a, 
then we get y = h(x). 
Therefore, f really behaves like a linear map, but with respect to the new origin a (not the standard origin 0). This is the case of a rotation around an axis that does not pass through the origin. 
Remark: A pair (a, (u1,...,un)) where (u1,...,un) is a basis of E and a is an origin chosen in E is called an a.ne frame. 
We now describe the trick which allows us to incorporate the translation part d into the matrix A. We de.ne the (m + 1) ≠ (n + 1) matrix A' obtained by .rst adding d as the (n + 1)th column and then (0,..., 0, 1) as the (m + 1)th row: 
n 
Ad 
A' 
= . 
0n 1 
It is clear that 
y Adx 
= 
10n 11 
i. y = Ax + d. 
This amounts to considering a point x Å Rn as a point (x, 1) in the (a.ne) hyperplane Hn+1 in Rn+1 of equation xn+1 = 1. Then an a.ne map is the restriction to the hyperplane Hn+1 of the linear map fffrom Rn+1 to Rm+1 corresponding to the matrix A' which maps Hn+1 into Hm+1 (ff(Hn+1) . Hm+1). Figure 5.4 illustrates this process for n = 2. 
For example, the map 
x1 11 x1 3 

∪ + 
x213 x2 0 

de.nes an a.ne map f which is represented in R3 by 
... ... 
x1 113 x1 
.x2. ∪ .130..x2.. 
1 001 1 

It is easy to check that the point a = (6, .3) is .xed by f, which means that f(a)= a, so by translating the coordinate frame to the origin a, the a.ne map behaves like a linear map. The idea of considering Rn as an hyperplane in Rn+1 can be used to de.ne projective maps. 


5.5 Summary 
The main concepts and results of this chapter are listed below: 
. 
Direct products, sums, direct sums. 

. 
Projections. 

. 
The fundamental equation 

dim(E) = dim(Ker f) + dim(Im f) = dim(Ker f) + rk(f) 
(The rank-nullity theorem; Theorem 5.8). 


. 
Grassmann’s relation 


dim(U) + dim(V ) = dim(U + V ) + dim(U ℃ V ). 

5.6. PROBLEMS 
. Characterizations of a bijective linear map f : E ∪ F . 

. Rank of a matrix. 

. A.ne Maps. 




5.6 Problems 
Problem 5.1. Let V and W be two subspaces of a vector space E. Prove that if V ″ W is 
a subspace of E, then either V . W or W . V . Problem 5.2. Prove that for every vector space E, if f : E ∪ E is an idempotent linear map, i.e., f . f = f, then we have a direct sum 
E = Ker f 쮵 Im f, so that f is the projection onto its image Im f. Problem 5.3. Let U1,...,Up be any p ◎ 2 subspaces of some vector space E and recall 
that the linear map 
a: U1 ≠· ··≠ Up ∪ E 
is given by a(u1,...,up)= u1 + ··· + up, with ui Å Ui for i =1,...,p. 
(1) If we let Zi . U1 ≠· ··≠ Up be given by 
In general, for any given i, the condition Ui ℃ p Uj = (0) does not necessarily 
j=1,j=i 
imply that Zi = (0). Thus, let 
     

령
p
Zi = u1,...,ui.1, . uj,ui+1,...,up 
j=1,j=i 
p령 p령  
uj Å Ui ℃  Uj  ,  
j=1,j=i  j=1,j=i  
for i = 1, . . . , p, then prove that  
Ker a = Z1 = · · · = Zp.  

w 

     

Z = u1,...,ui.1,ui,ui+1,...,up 
p령 p령  
ui = .  uj, ui Å Ui ℃  Uj  , 1 ● i ● p  .  
j=1,j=i  j=1,j=i  

Since Ker a = Z1 = ··· = Zp, we have Z = Ker a. Prove that if 
p
령 
Ui ℃ Uj = (0) 1 ● i ● p, 
j=1,j=i 
then Z = Ker a = (0). 
(2) Prove that U1 + ··· + Up is a direct sum i. 
p
령 
Ui ℃ Uj = (0) 1 ● i ● p. 
j=1,j=i 
Problem 5.4. Assume that E is .nite-dimensional, and let fi : E ∪ E be any p ◎ 2 linear maps such that 
f1 + ··· + fp = idE. 
Prove that the following properties are equivalent: 
(1) fi 2 = fi,1 ● i ● p. 
(2) fj . fi = 0, for all i = j,1 ● i, j ● p. Hint. Use Problem 5.2. 
Let U1,...,Up be any p ◎ 2 subspaces of some vector space E. Prove that U1 + ··· + Up is a direct sum i. 
i.1
령 
Ui ℃ Uj = (0),i =2, . . . , p. j=1 
Problem 5.5. Given any vector space E, a linear map f : E ∪ E is an involution if f . f = id. 
(1) Prove that an involution f is invertible. What is its inverse? 

(2) Let E1 and E.1 be the subspaces of E de.ned as follows: 
E1 = {u Å E | f(u)= u}



E.1 = {u Å E | f(u)= .u}. Prove that we have a direct sum 
E = E1 쮵 E.1. Hint. For every u Å E, write 
u + f(u) u . f(u) 
u =+ . 
22 
5.6. PROBLEMS 
(3) If E is .nite-dimensional and f is an involution, prove that there is some basis of E with respect to which the matrix of f is of the form 
Ik 0 Ik,n.k = ,
0 .In.k 
where Ik is the k ≠ k identity matrix (similarly for In.k) and k = dim(E1). Can you give a geometric interpretation of the action of f (especially when k = n . 1)? 
Problem 5.6. An n ≠ n matrix H is upper Hessenberg if hjk = 0 for all (j, k) such that j . k ◎ 0. An upper Hessenberg matrix is unreduced if hi+1i = 0 for i =1,...,n . 1. Prove that if H is a singular unreduced upper Hessenberg matrix, then dim(Ker (H)) = 1. 
Problem 5.7. Let A be any n ≠ k matrix. 
(1) 
Prove that the k ≠ k matrix ATA and the matrix A have the same nullspace. Use this to prove that rank(ATA) = rank(A). Similarly, prove that the n ≠ n matrix AAT and the matrix AT have the same nullspace, and conclude that rank(AAT) = rank(AT). 

We will prove later that rank(AT) = rank(A). 

(2) 
Let a1,...,ak be k linearly independent vectors in Rn (1 ● k ● n), and let A be the n ≠ k matrix whose ith column is ai. Prove that ATA has rank k, and that it is invertible. Let P = A(ATA).1AT (an n ≠ n matrix). Prove that 


P 2 
= P P T 
= P. 
What is the matrix P when k = 1? 
(3) Prove that the image of P is the subspace V spanned by a1,...,ak, or equivalently the set of all vectors in Rn of the form Ax, with x Å Rk . Prove that the nullspace U of P is the set of vectors u Å Rn such that ATu = 0. Can you give a geometric interpretation of U? 
Conclude that P is a projection of Rn onto the subspace V spanned by a1,...,ak, and that 
Rn 
= U 쮵 V. 
Problem 5.8. A rotation R┍ in the plane R2 is given by the matrix 
cos ┍ . sin ┍ 
R┍ = . 
sin ┍ cos ┍ 
(1) Use Matlab to show the action of a rotation R┍ on a simple .gure such as a triangle or a rectangle, for various values of ┍, including ┍ = ┪/6, ┪/4, ┪/3, ┪/2. 
(2) 
Prove that R┍ is invertible and that its inverse is R.┍. 

(3) 
For any two rotations R┒ and R┑, prove that 


R┑ . R┒ = R┒ . R┑ = R┒+┑. 
Use (2)-(3) to prove that the rotations in the plane form a commutative group denoted SO(2). 
Problem 5.9. Consider the a.ne map R┍,(a1,a2) in R2 given by 
y1 cos ┍ . sin ┍x1 a1 
=+ . 
y2 sin ┍ cos ┍x2 a2 
(1) Prove that if ┍ = k2┪, with k Å Z, then R┍,(a1,a2) has a unique .xed point (c1,c2), that is, there is a unique point (c1,c2) such that 
c1 c1 = R┍,(a1,a2) , 
c2 c2 
and this .xed point is given by 
c1 1 cos(┪/2 . ┍/2) . sin(┪/2 . ┍/2) a1 
= . 
c2 2 sin(┍/2) sin(┪/2 . ┍/2) cos(┪/2 . ┍/2) a2 
(2) In this question we still assume that ┍ = k2┪, with k Å Z. By translating the coordinate system with origin (0, 0) to the new coordinate system with origin (c1,c2), which means that if (x1,x2) are the coordinates with respect to the standard origin (0, 0) and if (x1' ,x 2' ) are the coordinates with respect to the new origin (c1,c2), we have 
x1 = x1 ' + c1 
x2 = x2 ' + c2 
and similarly for (y1,y2) and (y1' ,y 2' ), then show that 
y1 x1 = R┍,(a1,a2)
y2 x2 
becomes 
y ' x ' 
1 = R┍ 1 . 
y ' x ' 
22 
Conclude that with respect to the new origin (c1,c2), the a.ne map R┍,(a1,a2) becomes the rotation R┍. We say that R┍,(a1,a2) is a rotation of center (c1,c2). 
(3) 
Use Matlab to show the action of the a.ne map R┍,(a1,a2) on a simple .gure such as a triangle or a rectangle, for ┍ = ┪/3 and various values of (a1,a2). Display the center (c1,c2) of the rotation. 

What kind of transformations correspond to ┍ = k2┪, with k Å Z? 

(4) 
Prove that the inverse of R┍,(a1,a2) is of the form R.┍,(b1,b2), and .nd (b1,b2) in terms of ┍ and (a1,a2). 


(5) Given two a.ne maps R┒,(a1,a2) and R┑,(b1,b2), prove that 
R┑,(b1,b2) . R┒,(a1,a2) = R┒+┑,(t1,t2) 
for some (t1,t2), and .nd (t1,t2) in terms of ┑,(a1,a2) and (b1,b2). 

5.6. PROBLEMS 
Even in the case where (a1,a2) = (0, 0), prove that in general 
R┑,(b1,b2) . R┒ = R┒ . R┑,(b1,b2). 
Use (4)-(5) to show that the a.ne maps of the plane de.ned in this problem form a nonabelian group denoted SE(2). 
Prove that R┑,(b1,b2) . R┒,(a1,a2) is not a translation (possibly the identity) i. ┒ + ┑ = k2┪, for all k Å Z. Find its center of rotation when (a1,a2) = (0, 0). 
If ┒ + ┑ = k2┪, then R┑,(b1,b2) . R┒,(a1,a2) is a pure translation. Find the translation vector of R┑,(b1,b2) . R┒,(a1,a2). 
Problem 5.10. (A.ne subspaces) A subset A of Rn is called an a.ne subspace if either A = ., or there is some vector a Å Rn and some subspace U of Rn such that 
A = a + U = {a + u | u Å U}. 
We de.ne the dimension dim(A) of A as the dimension dim(U) of U. 
(1) If A = a + U, why is a ÅA? 
What are a.ne subspaces of dimension 0? What are a.ne subspaces of dimension 1 

(begin with R2)? What are a.ne subspaces of dimension 2 (begin with R3)? Prove that any nonempty a.ne subspace is closed under a.ne combinations. 
(2) Prove that if A = a + U is any nonempty a.ne subspace, then A = b + U for any b ÅA. 
(3) Let A be any nonempty subset of Rn closed under a.ne combinations. For any 
a ÅA, prove that Ua = {x . a Å Rn | x Å A} 
is a (linear) subspace of Rn such that 
A = a + Ua. 
Prove that Ua does not depend on the choice of a ÅA; that is, Ua = Ub for all a, b ÅA. In fact, prove that 
Ua = U = {y . x Å Rn | x, y Å A}, for all a ÅA, 
and so A = a + U, for any a ÅA. 
Remark: The subspace U is called the direction of A. 
(4) Two nonempty a.ne subspaces A and B are said to be parallel i. they have the same direction. Prove that that if A = B and A and B are parallel, then A℃B = .. 
Remark: The above shows that a.ne subspaces behave quite di.erently from linear sub-spaces. 
Problem 5.11. (A.ne frames and a.ne maps) For any vector v =(v1,...,vn) Å Rn, let v Å Rn+1
fbe the vector vf=(v1,...,vn, 1). Equivalently, vf=(fv1,..., vfn+1) Å Rn+1 is the vector de.ned by 
vi if 1 ● i ● n,
f=
vi 
1 if i = n +1. 
(1) 
For any m + 1 vectors (u0,u1,...,um) with ui Å Rn and m ● n, prove that if the m vectors (u1 . u0,...,um . u0) are linearly independent, then the m + 1 vectors (uf0,..., ufm) are linearly independent. 

(2) 
Prove that if the m + 1 vectors (uf0,..., ufm) are linearly independent, then for any choice of i, with 0 ● i ● m, the m vectors uj . ui for j Å{0,...,m} with j . i = 0 are linearly independent. 


Any m + 1 vectors (u0,u1,...,um) such that the m + 1 vectors (uf0,..., ufm) are linearly independent are said to be a.nely independent. 
From (1) and (2), the vector (u0,u1,...,um) are a.nely independent i. for any any choice of i, with 0 ● i ● m, the m vectors uj . ui for j Å{0,...,m} with j . i = 0 are linearly independent. If m = n, we say that n + 1 a.nely independent vectors (u0,u1,...,un) form an a.ne frame of Rn . 
(3) if (u0,u1,...,un) is an a.ne frame of Rn, then prove that for every vector v Å Rn , there is a unique (n + 1)-tuple (┡0,┡1,...,┡n) Å Rn+1, with ┡0 + ┡1 + ··· + ┡n = 1, such that 
v = ┡0u0 + ┡1u1 + ··· + ┡nun. 
The scalars (┡0,┡1,...,┡n) are called the barycentric (or a.ne) coordinates of v w.r.t. the a.ne frame (u0,u1,...,un). 
If we write ei = ui . u0, for i =1,...,n, then prove that we have 
v = u0 + ┡1e1 + ··· + ┡nen, 
and since (e1,...,en) is a basis of Rn (by (1) & (2)), the n-tuple (┡1,...,┡n) consists of the standard coordinates of v . u0 over the basis (e1,...,en). 
Conversely, for any vector u0 Å Rn and for any basis (e1,...,en) of Rn, let ui = u0 + ei for i =1,...,n. Prove that (u0,u1,...,un) is an a.ne frame of Rn, and for any v Å Rn, if 
v = u0 + x1e1 + ··· + xnen, 
with (x1,...,xn) Å Rn (unique), then 
v = (1 . (x1 + ··· + xx))u0 + x1u1 + ··· + xnun, 
so that (1 . (x1 + ··· + xx)),x1, ··· ,xn), are the barycentric coordinates of v w.r.t. the a.ne frame (u0,u1,...,un). 
The above shows that there is a one-to-one correspondence between a.ne frames (u0,..., un) and pairs (u0, (e1,...,en)), with (e1,...,en) a basis. Given an a.ne frame (u0,...,un), 

5.6. PROBLEMS 
we obtain the basis (e1,...,en) with ei = ui . u0, for i =1,...,n; given the pair (u0, (e1,..., en)) where (e1,...,en) is a basis, we obtain the a.ne frame (u0,...,un), with ui = u0 + ei, for i =1,...,n. There is also a one-to-one correspondence between barycentric coordinates 
w.r.t. the a.ne frame (u0,...,un) and standard coordinates w.r.t. the basis (e1,...,en). The barycentric cordinates (┡0,┡1,...,┡n) of v (with ┡0 + ┡1 + ··· + ┡n = 1) yield the standard coordinates (┡1,...,┡n) of v . u0; the standard coordinates (x1,...,xn) of v . u0 yield the barycentric coordinates (1 . (x1 + ··· + xn),x1,...,xn) of v. 
(4) 
Let (u0,...,un) be any a.ne frame in Rn and let (v0,...,vn) be any vectors in Rm . Prove that there is a unique a.ne map f : Rn ∪ Rm such that 

f(ui)= vi,i =0, . . . , n. 

(5) 
Let (a0,...,an) be any a.ne frame in Rn and let (b0,...,bn) be any n + 1 points in Rn . Prove that there is a unique (n + 1) ≠ (n + 1) matrix 


B  w  
A =  
0  1  

corresponding to the unique a.ne map f such that 
f
f(ai)= bi,i =0, . . . , n, 
in the sense that 
bi,Afai 
i =0, . . . , n, 

=

and that A is given by 

 

 .1
f
ff
f
ff
b0 b1 bna0a1 an
Make sure to prove that the bottom row of A is (0,..., 0, 1). 
In the special case where (a0,...,an) is the canonical a.ne frame with ai = ei+1 for i =0,...,n . 1 and an = (0,..., 0) (where ei is the ith canonical basis vector), show that
A =

···

···

. 

.
. 
10 ··· 00 01 ··· 00 
.. 
.
...... 

...... 

 
fffff
 
a0a1 an
a1 an
···

= 

.. 

.
. 00
.. 

00 ··· 10 
11 ··· 11 

and

.
. 
10 ··· 00 
01 ··· 00

...... 

...... 

f 
a0
 .1 
.. 
.
···

= 

.. 

.
. 00 
.

.. 

00 ··· 10 
.1 .1 ··· .11 

.  . .  .  .  .  
x1  x2  x3  1  0  0  x1 . x3  x2 . x3  x3  
A = .y1  y2  y3. . 0  1  0. = .y1 . y3  y2 . y3  y3 ..  
1  1  1  .1  .1  1  0  0  1  

(6) Recall that a nonempty a.ne subspace A of Rn is any nonempty subset of Rn closed under a.ne combinations. For any a.ne map f : Rn ∪ Rm, for any a.ne subspace A of Rn, and any a.ne subspace B of Rm, prove that f(A) is an a.ne subspace of Rm, and that f.1(B) is an a.ne subspace of Rn . 



