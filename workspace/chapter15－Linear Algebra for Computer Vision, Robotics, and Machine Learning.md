Chapter 15 

Unit Quaternions and Rotations in SO(3) 
This chapter is devoted to the representation of rotations in SO(3) in terms of unit quater-nions. Since we already de.ned the unitary groups SU(n), the quickest way to introduce the unit quaternions is to de.ne them as the elements of the group SU(2). 
The skew .eld H of quaternions and the group SU(2) of unit quaternions are discussed in Section 15.1. In Section 15.2, we de.ne a homomorphism r : SU(2) ¡ú SO(3) and prove that its kernel is {.I,I}. We compute the rotation matrix Rq associated with the rotation rq induced by a unit quaternion q in Section 15.3. In Section 15.4, we prove that the homomor-phism r : SU(2) ¡ú SO(3) is surjective by providing an algorithm to construct a quaternion from a rotation matrix. In Section 15.5 we de.ne the exponential map exp: su(2) ¡ú SU(2) where su(2) is the real vector space of skew-Hermitian 2 ¡Á 2 matrices with zero trace. We prove that exponential map exp: su(2) ¡ú SU(2) is surjective and give an algorithm for .nding a logarithm. We discuss quaternion interpolation and prove the famous slerp inter-polation formula due to Ken Shoemake in Section 15.6. This formula is used in robotics and computer graphics to deal with interpolation problems. In Section 15.7, we prove that there is no ¡°nice¡± section s: SO(3) ¡ú SU(2) of the homomorphism r : SU(2) ¡ú SO(3), in the sense that any section of r is neither a homomorphism nor continuous. 
15.1 	The Group SU(2) of Unit Quaternions and the Skew Field H of Quaternions 
De.nition 15.1. The unit quaternions are the elements of the group SU(2), namely the group of 2 ¡Á 2 complex matrices of the form 
.  .  
¦Á .¦Â  ¦Â ¦Á  ¦Á, ¦Â ¡Ê C, ¦Á¦Á + ¦Â¦Â = 1.  

The quaternions are the elements of the real vector space H = RSU(2). 

521 

Let 1, i, j, k be the matrices 
10 i 0 010 i 
1 = , i = , j = , k = ,
01 0 .i .10 i 0 
then H is the set of all matrices of the form 
X = a1 + bi + cj + dk, a,b,c,d ¡Ê R. 
Indeed, every matrix in H is of the form 
a + ib c + id 
X = , a,b,c,d ¡Ê R.
.(c . id) a . ib 
It is easy (but a bit tedious) to verify that the quaternions 1, i, j, k satisfy the famous identities discovered by Hamilton: 
i2 = j2 = k2 = ijk = .1, 
ij = .ji = k, 
jk = .kj = i, 
ki = .ik = j. 
Thus, the quaternions are a generalization of the complex numbers, but there are three square roots of .1 and multiplication is not commutative. 
Given any two quaternions X = a1+ bi + cj + dk and Y = a'1+ b'i+ c'j + d'k, Hamilton¡¯s famous formula 
XY =(aa' . bb' . cc' . dd')1 +(ab' + ba' + cd' . dc')i ''
+(ac+ ca+ db' . bd')j +(ad' + da' + bc' . cb')k 
looks mysterious, but it is simply the result of multiplying the two matrices 
a + ib c + id a' + ib' c' + id' 
X = and Y = .
.(c . id) a . ib .(c' . id') a' . ib' 
It is worth noting that this formula was discovered independently by Olinde Rodrigues in 1840, a few years before Hamilton (Veblen and Young [72]). However, Rodrigues was working with a di.erent formalism, homogeneous transformations, and he did not discover the quaternions. 
If 
a + ib c + id 
X = , a,b,c,d ¡Ê R,
.(c . id) a . ib it is immediately veri.ed that 
XX . = X . X =(a 2 + b2 + c 2 + d2)1. 
15.2. REPRESENTATION OF ROTATION IN SO(3) BY QUATERNIONS IN SU(2)523 

Also observe that 

a . ib .(c + id)
X . =	= a1 . bi . cj . dk. 
c . id a + ib 
This implies that if X 

= 0, then X 	is invertible and its inverse is given by 2 + b22 + d2).1X . 
X.1 =(a + c. 
As a consequence, it can be veri.ed that H is a skew .eld (a noncommutative .eld). It is also a real vector space of dimension 4 with basis (1, i, j, k); thus as a vector space, H is isomorphic to R4 . 
De.nition 15.2. A concise notation for the quaternion X de.ned by ¦Á = a + ib and ¦Â = c + id is 
X =[a, (b, c, d)]. We call a the scalar part of X and (b, c, d) the vector part of X. With this notation, X. =[a, .(b, c, d)], which is often denoted by X. The quaternion X is called the conjugate of q. If q is a unit quaternion, then q is the multiplicative inverse of q. 


15.2 	Representation of Rotations in SO(3) by Quater-nions in SU(2) 
The key to representation of rotations in SO(3) by unit quaternions is a certain group homomorphism called the adjoint representation of SU(2). To de.ne this mapping, .rst we de.ne the real vector space su(2) of skew Hermitian matrices. 
De.nition 15.3. The (real) vector space su(2) of 2 ¡Á 2 skew Hermitian matrices with zero trace is given by 
 

su(2) =

ix y + iz 

.y + iz .ix

    

(x, y, z) ¡Ê R3
 

. 

Observe that for every matrix A ¡Ê su(2), we have A. = .A, that is, A is skew Hermitian, and that tr(A) = 0. 
De.nition 15.4. The adjoint representation of the group SU(2) is the group homomorphism Ad: SU(2) ¡ú GL(su(2)) de.ned such that for every q ¡Ê SU(2), with 
q =  ¦Á .¦Â  ¦Â ¦Á  ¡Ê SU(2),  
we have  

Adq(A)= qAq . ,A ¡Ê su(2), where q . is the inverse of q (since SU(2) is a unitary group) and is given by ¦Á .¦Â 
q . = . 
¦Â¦Á 
One needs to verify that the map Adq is an invertible linear map from su(2) to itself, and that Ad is a group homomorphism, which is easy to do. 
In order to associate a rotation ¦Ñq (in SO(3)) to q, we need to embed R3 into H as the pure quaternions, by 
ix y + iz 
¦×(x, y, z)= , (x, y, z) ¡Ê R3 .
.y + iz .ix 
Then q de.nes the map ¦Ñq (on R3) given by 
¦Ñq(x, y, z)= ¦×.1(q¦×(x, y, z)q . ). 
Therefore, modulo the isomorphism ¦×, the linear map ¦Ñq is the linear isomorphism Adq. In fact, it turns out that ¦Ñq is a rotation (and so is Adq), which we will prove shortly. So, the representation of rotations in SO(3) by unit quaternions is just the adjoint representation of SU(2); its image is a subgroup of GL(su(2)) isomorphic to SO(3). 
Technically, it is a bit simpler to embed R3 in the (real) vector spaces of Hermitian matrices with zero trace, 
xz . iy x, y, z ¡Ê R . 
z + iy .x 
Since the matrix ¦×(x, y, z) is skew-Hermitian, the matrix .i¦×(x, y, z) is Hermitian, and we have 
xz . iy
.i¦×(x, y, z)= = x¦Ò3 + y¦Ò2 + z¦Ò1, 
z + iy .x where ¦Ò1,¦Ò2,¦Ò3 are the Pauli spin matrices 
01 0 .i 10 
¦Ò1 = ,¦Ò2 = ,¦Ò3 = . 
10 i 00 .1 
Matrices of the form x¦Ò3 + y¦Ò2 + z¦Ò1 are Hermitian matrices with zero trace. 
It is easy to see that every 2 ¡Á 2 Hermitian matrix with zero trace must be of this form. (observe that (i¦Ò1, i¦Ò2, i¦Ò3) forms a basis of su(2). Also, i = i¦Ò3, j = i¦Ò2, k = i¦Ò1.) 
Now, if A = x¦Ò3 + y¦Ò2 + z¦Ò1 is a Hermitian 2 ¡Á 2 matrix with zero trace, we have 
(qAq . ) . = qA . q . = qAq . , 
so qAq. is also Hermitian, and 
tr(qAq . ) = tr(Aq . q) = tr(A), 
and qAq. also has zero trace. Therefore, the map A ¡ú qAq. preserves the Hermitian matrices with zero trace. We also have 
xz . iy 22
det(x¦Ò3 + y¦Ò2 + z¦Ò1)=det = .(x + y + z 2), 
z + iy .x 
15.2. REPRESENTATION OF ROTATION IN SO(3) BY QUATERNIONS IN SU(2)525 

and det(qAq . ) = det(q) det(A) det(q . ) = det(A)= .(x 2 + y 2 + z 2). We can embed R3 into the space of Hermitian matrices with zero trace by .(x, y, z)= x¦Ò3 + y¦Ò2 + z¦Ò1. Note that 
..1
. = .i¦× and = i¦×.1 . De.nition 15.5. The unit quaternion q ¡Ê SU(2) induces a map rq on R3 by rq(x, y, z)= ..1(q.(x, y, z)q . )= ..1(q(x¦Ò3 + y¦Ò2 + z¦Ò1)q . ). The map rq is clearly linear since . is linear. 
Proposition 15.1. For every unit quaternion q ¡Ê SU(2), the linear map rq is orthogonal, 
that is, rq ¡Ê O(3). 
Proof. Since

.£ü (x, y, z)£ü2 = .(x 2 + y 2 + z 2) = det(x¦Ò3 + y¦Ò2 + z¦Ò1) = det(.(x, y, z)), we have
.£ü rq
(x, y, z)£ü2 
= det(.(rq
. )
  
  
(x, y, z))) = det(q(x¦Ò3 + y¦Ò2 + z¦Ò1)q (x, y, z)2
,

= det(x¦Ò3 + y¦Ò2 + z¦Ò1)= .
and we deduce that rq is an isometry. Thus, rq ¡Ê O(3). 
In fact, rq is a rotation, and we can show this by .nding the .xed points of rq. Let q be a unit quaternion of the form 
¦Á¦Â 
q = 
.¦Â¦Á 
with ¦Á = a + ib, ¦Â = c + id, and a2 + b2 + c2 + d2 =1(a, b, c, d ¡Ê R). If b = c = d = 0, then q = I and rq is the identity so we may assume that (b, c, d)= (0, 0, 0). 
Proposition 15.2. If (b, c, d) = (0, 0, 0), then the .xed points of rq are solutions (x, y, z) of the linear system 
.dy + cz =0 cx . by =0 dx . bz =0. 
This linear system has the nontrivial solution (b, c, d) and has rank 2. Therefore, rq has the eigenvalue 1 with multiplicity 1, and rq is a rotation whose axis is determined by (b, c, d). 
Proof. We have rq(x, y, z)=(x, y, z) i. ..1(q(x¦Ò3 + y¦Ò2 + z¦Ò1)q . )=(x, y, z) i. q(x¦Ò3 + y¦Ò2 + z¦Ò1)q . = .(x, y, z), and since .(x, y, z)= x¦Ò3 + y¦Ò2 + z¦Ò1 = A with 
xz . iy
A = , 
z + iy .x we see that rq(x, y, z)=(x, y, z) i. qAq . = A i. qA = Aq. We have ¦Á¦Â xz . iy ¦Áx + ¦Âz + i¦Ây ¦Áz . i¦Áy . ¦Âx 
qA == 
.¦Â¦Á z + iy .x .¦Âx + ¦Áz + i¦Áy .¦Âz + i¦Ây . ¦Áx and 
xz . iy¦Á¦Â ¦Áx . ¦Âz + i¦Ây ¦Âx + ¦Áz . i¦Áy
Aq == . 
z + iy .x .¦Â¦Á ¦Áz + i¦Áy + ¦Âx ¦Âz + i¦Ây . ¦Áx By equating qA and Aq, we get 
i(¦Â . ¦Â)y +(¦Â + ¦Â)z =0 
2¦Âx + i(¦Á . ¦Á)y +(¦Á . ¦Á)z =0 
2¦Âx + i(¦Á . ¦Á)y +(¦Á . ¦Á)z =0 
i(¦Â . ¦Â)y +(¦Â + ¦Â)z =0. 
The .rst and the fourth equation are identical and the third equation is obtained by conju-gating the second, so the above system reduces to 
i(¦Â . ¦Â)y +(¦Â + ¦Â)z =0 2¦Âx + i(¦Á . ¦Á)y +(¦Á . ¦Á)z =0. 
Replacing ¦Á by a + ib and ¦Â by c + id, we get 
.dy + cz =0 cx . by + i(dx . bz)=0, 

15.2. REPRESENTATION OF ROTATION IN SO(3) BY QUATERNIONS IN SU(2)527 

which yields the equations 
.dy + cz =0 cx . by =0 dx . bz =0. 
This linear system has the nontrivial solution (b, c, d) and the matrix of this system is 
.  .  
0  .d  c  
.c  .b  0 . .  
d  0  .b  

Since (b, c, d) =(0, 0, 0), this matrix always has a 2 ¡Á 2 submatrix which is nonsingular, so it has rank 2, and consequently its kernel is the one-dimensional space spanned by (b, c, d). Therefore, rq has the eigenvalue 1 with multiplicity 1. If we had det(rq)= .1, then the eigenvalues of rq would be either (.1, 1, 1) or (.1,ei¦È,e.i¦È) with ¦È = k2¦Ð (with k ¡Ê Z), contradicting the fact that 1 is an eigenvalue with multiplicity 1. Therefore, rq is a rotation; in fact, its axis is determined by (b, c, d). 
In summary, q ¡ú rq is a map r from SU(2) to SO(3). 
Theorem 15.3. The map r : SU(2) ¡ú SO(3) is homomorphism whose kernel is {I, .I}. 
Proof. This map is a homomorphism, because if q1,q2 ¡Ê SU(2), then 
rq2 (rq1 (x, y, z)) = ..1(q2.(rq1 (x, y, z))q2. ) 
= ..1(q2.(..1(q1.(x, y, z)q1. ))q2. ) 
= ..1((q2q1).(x, y, z)(q2q1) . ) 
= rq2q1 (x, y, z). 
The computation that showed that if (b, c, d) = (0, 0, 0), then rq has the eigenvalue 1 with multiplicity 1 implies the following: if rq = I3, namely rq has the eigenvalue 1 with multi-plicity 3, then (b, c, d) = (0, 0, 0). But then a = ¡À1, and so q = ¡ÀI2. Therefore, the kernel of the homomorphism r : SU(2) ¡ú SO(3) is {I, .I}. 
Remark: Perhaps the quickest way to show that r maps SU(2) into SO(3) is to observe that the map r is continuous. Then, since it is known that SU(2) is connected, its image by r lies in the connected component of I, namely SO(3). 
The map r is surjective, but this is not obvious. We will return to this point after .nding the matrix representing rq explicitly. 


15.3 Matrix Representation of the Rotation rq Given a unit quaternion q of the form 
¦Á¦Â 
q = 
.¦Â¦Á 
with ¦Á = a + ib, ¦Â = c + id, and a2 + b2 + c2 + d2 =1(a, b, c, d ¡Ê R), to .nd the matrix representing the rotation rq we need to compute 
¦Á¦Â xz . iy ¦Á .¦Â 
q(x¦Ò3 + y¦Ò2 + z¦Ò1)q . = . 
.¦Â¦Á z + iy .x ¦Â¦Á 
First we have 
xz . iy ¦Á .¦Â x¦Á + z¦Â . iy¦Â .x¦Â + z¦Á . iy¦Á 
= . 
z + iy .x ¦Â¦Á z¦Á + iy¦Á . x¦Â .z¦Â . iy¦Â . x¦Á 
Next, we have 
¦Á¦Â x¦Á + z¦Â . iy¦Â .x¦Â + z¦Á . iy¦Á A1 A2 .¦Â¦Á z¦Á + iy¦Á . x¦Â .z¦Â . iy¦Â . x¦Á = A3 A4 , 
with 
A1 =(¦Á¦Á . ¦Â¦Â)x + i(¦Á¦Â . ¦Á¦Â)y +(¦Á¦Â + ¦Á¦Â)z A2 = .2¦Á¦Âx . i(¦Á2 + ¦Â2)y +(¦Á2 . ¦Â2)z 
A3 = .2¦Á¦Âx + i(¦Á2 + ¦Â 2 )y +(¦Á2 . ¦Â 2 )z 
A4 = .(¦Á¦Á . ¦Â¦Â)x . i(¦Á¦Â . ¦Á¦Â)y . (¦Á¦Â + ¦Á¦Â)z. 

Since ¦Á = a + ib and ¦Â = c + id, with a, b, c, d ¡Ê R, we have 
2 . d2
¦Á¦Á . ¦Â¦Â = a 2 + b2 . c i(¦Á¦Â . ¦Á¦Â) = 2(bc . ad) ¦Á¦Â + ¦Á¦Â = 2(ac + bd) .¦Á¦Â = .ac + bd . i(ad + bc) .i(¦Á2 + ¦Â2) = 2(ab + cd) . i(a 2 . b2 + c 2 . d2) ¦Á2 . ¦Â22 . b2 . c 2 + d2 
= a + i2(ab . cd). 
Using the above, we get 
(¦Á¦Á . ¦Â¦Â)x + i(¦Á¦Â . ¦Á¦Â)y +(¦Á¦Â + ¦Á¦Â)z =(a 2 + b2 . c 2 . d2)x + 2(bc . ad)y + 2(ac + bd)z, 
15.3. MATRIX REPRESENTATION OF THE ROTATION RQ 
and 
. 2¦Á¦Âx . i(¦Á2 + ¦Â2)y +(¦Á2 . ¦Â2)z 
= 2(.ac + bd)x + 2(ab + cd)y +(a 2 . b2 . c 2 + d2)z 
. i[2(ad + bc)x +(a 2 . b2 + c 2 . d2)y + 2(.ab + cd)z]. 
If we write 
x ' z ' . iy ' 
q(x¦Ò3 + y¦Ò2 + z¦Ò1)q . = , 
z ' + iy ' .x ' 
we obtain 
x ' =(a 2 + b2 . c 2 . d2)x + 2(bc . ad)y + 2(ac + bd)z 
y ' = 2(ad + bc)x +(a 2 . b2 + c 2 . d2)y + 2(.ab + cd)z 
z ' = 2(.ac + bd)x + 2(ab + cd)y +(a 2 . b2 . c 2 + d2)z. 
In summary, we proved the following result. 
Proposition 15.4. The matrix representing rq is 
.. 
a2 + b2 . c2 . d2 2bc . 2ad 2ac +2bd Rq = . 2bc +2ad a2 . b2 + c2 . d2 .2ab +2cd . . .2ac +2bd 2ab +2cd a2 . b2 . c2 + d2 
Since a2 + b2 + c2 + d2 =1, this matrix can also be written as 
.. 
2a2 +2b2 . 12bc . 2ad 2ac +2bd Rq = . 2bc +2ad 2a2 +2c2 . 1 .2ab +2cd .. .2ac +2bd 2ab +2cd 2a2 +2d2 . 1 
The above is the rotation matrix in Euler form induced by the quaternion q, which is the matrix corresponding to ¦Ñq. This is because ..1
. = .i¦×, = i¦×.1 , 
so 
rq(x, y, z)= ..1(q.(x, y, z)q . )= i¦×.1(q(.i¦×(x, y, z))q . ) 
= ¦×.1(q¦×(x, y, z)q . )= ¦Ñq(x, y, z), 
and so rq = ¦Ñq. 
We showed that every unit quaternion q ¡Ê SU(2) induces a rotation rq ¡Ê SO(3), but it is not obvious that every rotation can be represented by a quaternion. This can shown in various ways. 
One way to is use the fact that every rotation in SO(3) is the composition of two re.ec-tions, and that every re.ection ¦Ò of R3 can be represented by a quaternion q, in the sense that 
¦Ò(x, y, z)= ...1(q.(x, y, z)q . ). 
Note the presence of the negative sign. This is the method used in Gallier [25] (Chapter 9). 



15.4 	An Algorithm to Find a Quaternion Representing a Rotation 
Theorem 15.5. The homomorphim r : SU(2) ¡ú SO(3) is surjective. 
Here is an algorithmic method to .nd a unit quaternion q representing a rotation matrix R, which provides a proof of Theorem 15.5. Let 
a + ib 	c + id 2 + b22 
q = 	,a + c + d2 =1, a,b,c,d ¡Ê R. 
.(c . id) a . ib First observe that the trace of Rq is given by 
tr(Rq)=3a 2 . b2 . c 2 . d2 , 
but since a2 + b2 + c2 + d2 = 1, we get tr(Rq)=4a2 . 1, so 
2 tr(Rq)+1 
a = . 
4 
If R ¡Ê SO(3) is any rotation matrix and if we write 
.. 
r11 r12 r13 
..
R = 	r21 r22 r23 r31 r32 r33, 
we are looking for a unit quaternion q ¡Ê SU(2) such that Rq = R. Therefore, we must have 
2 tr(R)+1 
a = . 
4 We also know that tr(R) = 1+2cos ¦È, where ¦È ¡Ê [0,¦Ð] is the angle of the rotation R, so we get 
cos ¦È +1 ¦È
2	2 
a = = cos ,
22 
which implies that 
¦È 
|a| = cos (0 ¡Ü ¦È ¡Ü ¦Ð). 
2 
Note that we may assume that ¦È ¡Ê [0,¦Ð], because if ¦Ð ¡Ü ¦È ¡Ü 2¦Ð, then ¦È . 2¦Ð ¡Ê [.¦Ð, 0], and then the rotation of angle ¦È . 2¦Ð and axis determined by the vector (b, c, d) is the same as the rotation of angle 2¦Ð . ¦È ¡Ê [0,¦Ð] and axis determined by the vector .(b, c, d). There are two cases. 
15.4. AN ALGORITHM TO FIND A QUATERNION REPRESENTING A ROTATION531 

Case 1 . tr(R)= .1, or equivalently ¦È = ¦Ð. In this case a = 0. Pick 
 
tr(R)+1 
a =. 
2 
Then by equating R . RT and Rq . RqT, we get 
4ab = r32 . r23 4ac = r13 . r31 4ad = r21 . r12, 
which yields 
r32 . r23 r13 . r31 r21 . r12
b = ,c = ,d = . 
4a 4a 4a 
Case 2 . tr(R)= .1, or equivalently ¦È = ¦Ð. In this case a = 0. By equating R + RT and Rq + RqT, we get 
4bc = r21 + r12 4bd = r13 + r31 4cd = r32 + r23. 
By equating the diagonal terms of R and Rq, we also get 
1+ r11
b2 
= 
2 2 1+ r22 
c = 
2 1+ r33
d2 
= . 
2 
Since q = 0 and a = 0, at least one of b, c, d is nonzero. If b = 0, let 
¡Ì 
1+ r11
b =  ¡Ì  ,  
2  
and determine c, d using  

4bc = r21 + r12 4bd = r13 + r31. 
If c = 0, let 
¡Ì 
1+ r22 
c = ¡Ì , 
2 
and determine b, d using 
4bc = r21 + r12 
4cd = r32 + r23. 
If d = 0, let 
¡Ì 
1+ r33
d = ¡Ì , 
2 and determine b, c using 
4bd = r13 + r31 
4cd = r32 + r23. 
It is easy to check that whenever we computed a square root, if we had chosen a negative sign instead of a positive sign, we would obtain the quaternion .q. However, both q and .q determine the same rotation rq. 
The above discussion involving the cases tr(R)= .1 and tr(R)= .1 is reminiscent of the procedure for .nding a logarithm of a rotation matrix using the Rodrigues formula (see Section 11.7). This is not surprising, because if 
.. 
0 .u3 u2 
..
B = u3 0 .u1 .u2 u1 0 
222
and if we write ¦È = u1 + u2 + u3 (with 0 ¡Ü ¦È ¡Ü ¦Ð), then the Rodrigues formula says that 
sin ¦È (1 . cos ¦È) 
e B = I + B + B2 ,¦È =0,
¦È¦È2 
with e0 = I. It is easy to check that tr(eB) = 1+2cos ¦È. Then it is an easy exercise to check that the quaternion q corresponding to the rotation R = eB (with B = 0) is given by 
¦È ¦Èu1 u2 u3 
q = cos , sin ,, . 
22 ¦È¦È¦È 
So the method for .nding the logarithm of a rotation R is essentially the same as the method for .nding a quaternion de.ning R. 
Remark: Geometrically, the group SU(2) is homeomorphic to the 3-sphere S3 in R4 , 
S3 222 
= {(x, y, z, t) ¡Ê R4 | x + y + z + t2 =1}. 
However, since the kernel of the surjective homomorphism r : SU(2) ¡ú SO(3) is {I, .I}, as a topological space, SO(3) is homeomorphic to the quotient of S3 obtained by identifying antipodal points (x, y, z, t) and .(x, y, z, t). This quotient space is the (real) projective space RP3, and it is more complicated than S3 . The space S3 is simply-connected, but RP3 is not. 

15.5. THE EXPONENTIAL MAP EXP: SU(2) ¡ú SU(2) 


15.5 The Exponential Map exp: su(2) ¡ú SU(2) 
Given any matrix A ¡Ê su(2), with 
iu1 u2 + iu3
A = ,
.u2 + iu3 .iu1 
it is easy to check that A2 = .¦È2 10 ,
01 
222
with ¦È = u1 + u2 + u3. Then we have the following formula whose proof is very similar to the proof of the formula given in Proposition 8.22. 
Proposition 15.6. For every matrix A ¡Ê su(2), with 
iu1 u2 + iu3
A = ,
.u2 + iu3 .iu1 
222
if we write ¦È = u1 + u2 + u3, then 
sin ¦È 
e A = cos ¦ÈI + A, ¦È =0,
¦È 
0
and e= I. 
Therefore, by the discussion at the end of the previous section, eA is a unit quaternion representing the rotation of angle 2¦È and axis (u1,u2,u3) (or I when ¦È = k¦Ð, k ¡Ê Z). The above formula shows that we may assume that 0 ¡Ü ¦È ¡Ü ¦Ð. Proposition 15.6 shows that the exponential yields a map exp: su(2) ¡ú SU(2). It is an analog of the exponential map exp: so(3) ¡ú SO(3). 
Remark: Because so(3) and su(2) are real vector spaces of dimension 3, they are isomorphic, and it is easy to construct an isomorphism. In fact, so(3) and su(2) are isomorphic as Lie algebras, which means that there is a linear isomorphism preserving the the Lie bracket [A, B]= AB . BA. However, as observed earlier, the groups SU(2) and SO(3) are not isomorphic. 
An equivalent, but often more convenient, formula is obtained by assuming that u = (u1,u2,u3) is a unit vector, equivalently det(A) = 1, in which case A2 = .I, so we have 
¦ÈA 
e = cos ¦ÈI + sin ¦ÈA. 
Using the quaternion notation, this is read as 
e ¦ÈA = [cos ¦È, sin ¦Èu]. 
Proposition 15.7. The exponential map exp: su(2) ¡ú SU(2) is surjective 

Proof. We give an algorithm to .nd the logarithm A ¡Ê su(2) of a unit quaternion 
¦Á¦Â 
q = 
.¦Â¦Á with ¦Á = a + bi and ¦Â = c + id. If q = I (i.e. a = 1), then A = 0. If q = .I (i.e. a = .1), then i 0 
A = ¡À¦Ð. 
0 .i Otherwise, a = ¡À1 and (b, c, d) = (0, 0, 0), and we are seeking some A = ¦ÈB ¡Ê su(2) with det(B)=1 and 0 <¦È<¦Ð, such that, by Proposition 15.6, q = e ¦ÈB = cos ¦ÈI + sin ¦ÈB. Let 
iu1 u2 + iu3
B = ,
.u2 + iu3 .iu1 with u =(u1,u2,u3) a unit vector. We must have ¦ÈB . (e ¦ÈB) .. 
a = cos ¦È, e = q . q. Since 0 <¦È<¦Ð, we have sin ¦È = 0, and iu1 u2 + iu3 ¦Á . ¦Á 2¦Â 
2 sin ¦È = . 
.u2 + iu3 .iu1 .2¦Â¦Á . ¦Á 
Thus, we get 
11 
u1 = b, u2 + iu3 =(c + id);
sin ¦È sin ¦È that is, 
cos ¦È = a (0 <¦È<¦Ð) 1 
(u1,u2,u3)= (b, c, d). 
sin ¦È 
Since a2 +b2 +c2 +d2 = 1 and a = cos ¦È, the vector (b, c, d)/ sin ¦È is a unit vector. Furthermore if the quaternion q is of the form q = [cos ¦È, sin ¦Èu] where u =(u1,u2,u3) is a unit vector (with 0 <¦È<¦Ð), then 
iu1 u2 + iu3
A = ¦È (. log)
.u2 + iu3 .iu1 is a logarithm of q. 
15.6. QUATERNION INTERPOLATION . 
Observe that not only is the exponential map exp: su(2) ¡ú SU(2) surjective, but the above proof shows that it is injective on the open ball 
{¦ÈB ¡Ê su(2) | det(B)=1, 0 ¡Ü ¦È<¦Ð}. 
Also, unlike the situation where in computing the logarithm of a rotation matrix R ¡Ê SO(3) we needed to treat the case where tr(R)= .1 (the angle of the rotation is ¦Ð) in a special way, computing the logarithm of a quaternion (other than ¡ÀI) does not require any case analysis; no special case is needed when the angle of rotation is ¦Ð. 
15.6 Quaternion Interpolation . 
We are now going to derive a formula for interpolating between two quaternions. This formula is due to Ken Shoemake, once a Penn student and my TA! Since rotations in SO(3) can be de.ned by quaternions, this has applications to computer graphics, robotics, and computer vision. 
First we observe that multiplication of quaternions can be expressed in terms of the inner product and the cross-product in R3 . Indeed, if q1 =[a, u1] and q2 =[a2,u2], itcanbe veri.ed that 
q1q2 =[a1,u1][a2,u2]=[a1a2 . u1 ¡¤ u2,a1u2 + a2u1 + u1 ¡Á u2]. (. mult) 
We will also need the identity 
u ¡Á (u ¡Á v)=(u ¡¤ v)u . (u ¡¤ u)v. 
Given a quaternion q expressed as q = [cos ¦È, sin ¦Èu], where u is a unit vector, we can interpolate between I and q by .nding the logs of I and q, interpolating in su(2), and then exponentiating. We have 
00 iu1 u2 + iu3
A = log(I)= ,B = log(q)= ¦È,
00 .u2 + iu3 .iu1 
and so q = eB . Since SU(2) is a compact Lie group and since the inner product on su(2) given by 
(X, Y£© = tr(XTY ) 
is Ad(SU(2))-invariant, it induces a biinvariant Riemannian metric on SU(2), and the curve 
¦Ë ¡ú e ¦ËB ,¦Ë ¡Ê [0, 1] 
is a geodesic from I to q in SU(2). We write q¦Ë = e¦ËB . Given two quaternions q1 and q2, because the metric is left invariant, the curve 
¦Ë ¡ú Z(¦Ë)= q1(q1 .1 q2)¦Ë ,¦Ë ¡Ê [0, 1] 
is a geodesic from q1 to q2. Remarkably, there is a closed-form formula for the interpolant Z(¦Ë). 
Say q1 = [cos ¦È, sin ¦Èu] and q2 = [cos ., sin .v], and assume that q1 = q2 and q1 = .q2. First, we compute q.1q2. Since q.1 = [cos ¦È, . sin ¦Èu], we have 
q .1 q2 = [cos ¦È cos . + sin ¦È sin .(u ¡¤ v), . sin ¦È cos .u + cos ¦È sin .v . sin ¦È sin .(u ¡Á v)]. De.ne ¦¸ by cos¦¸ = cos ¦È cos . + sin ¦È sin .(u ¡¤ v). (. ¦¸) Since q1 = q2 and q1 = .q2, we have 0 < ¦¸ <¦Ð, so we get .1 (. sin ¦È cos .u + cos ¦È sin .v . sin ¦È sin .(u ¡Á v) 
q1 q2 = cos ¦¸, sin ¦¸ ,
sin ¦¸ where the term multiplying sin ¦¸ is a unit vector because q1 and q2 are unit quaternions, so q1 .1 q2 is also a unit quaternion. By (. log), we have .1 q2)¦Ë
(q1 
(. sin ¦È cos .u + cos ¦È sin .v . sin ¦È sin .(u ¡Á v) 

= cos ¦Ë¦¸, sin ¦Ë¦¸ . 
sin ¦¸ Next we need to compute q1(q1 .1 q2)¦Ë . The scalar part of this product is sin ¦Ë¦¸ sin ¦Ë¦¸ 
s = cos ¦È cos ¦Ë¦¸ + sin2 ¦È cos .(u ¡¤ u) . sin ¦È sin . cos ¦È(u ¡¤ v)
sin¦¸ sin¦¸ sin ¦Ë¦¸ 
+ sin2 ¦È sin .(u ¡¤ (u ¡Á v)). 
sin ¦¸ Since u ¡¤ (u ¡Á v) = 0, the last term is zero, and since u ¡¤ u = 1 and sin ¦È sin .(u ¡¤ v) = cos¦¸ . cos ¦È cos ., we get sin ¦Ë¦¸ sin ¦Ë¦¸ 
s = cos ¦È cos ¦Ë¦¸ + sin2 ¦È cos . . cos ¦È(cos ¦¸ . cos ¦È cos .)
sin¦¸ sin¦¸ sin ¦Ë¦¸ sin ¦Ë¦¸ 
= cos ¦È cos ¦Ë¦¸ + (sin2 ¦È + cos2 ¦È) cos . . cos ¦È cos ¦¸ 
sin¦¸ sin¦¸ (cos ¦Ë¦¸ sin ¦¸ . sin ¦Ë¦¸ cos ¦¸) cos ¦È sin ¦Ë¦¸ 
= + cos . 
sin¦¸ sin¦¸ 
sin(1 . ¦Ë)¦¸ sin ¦Ë¦¸ 

= cos ¦È + cos .. 
sin¦¸ sin¦¸ 
15.6. QUATERNION INTERPOLATION ¢à 
The vector part of the product q1(q1 .1 q2)¦Ë is given by 
sin ¦Ë¦¸ sin ¦Ë¦¸ 
¦Í = . cos ¦È sin ¦È cos .u + cos 2 ¦È sin .v 
sin¦¸ sin¦¸ 
sin ¦Ë¦¸ 

. cos ¦È sin ¦È sin .(u ¡Á v) + cos ¦Ë¦¸ sin ¦Èu 
sin ¦¸ sin ¦Ë¦¸ sin ¦Ë¦¸ 
. sin2 ¦È cos .(u ¡Á u) + cos ¦È sin ¦È sin .(u ¡Á v)
sin¦¸ sin¦¸ sin ¦Ë¦¸ . sin2 ¦È sin .(u ¡Á (u ¡Á v)). 
sin ¦¸ 
We have u ¡Á u = 0, the two terms involving u ¡Á v cancel out, 
u ¡Á (u ¡Á v)=(u ¡¤ v)u . (u ¡¤ u)v, 
and u ¡¤ u = 1, so we get 
sin ¦Ë¦¸ sin ¦Ë¦¸ 
¦Í = . cos ¦È sin ¦È cos .u + cos ¦Ë¦¸ sin ¦Èu + cos 2 ¦È sin .v 
sin¦¸ sin¦¸ 
sin ¦Ë¦¸ sin ¦Ë¦¸ 

+ sin2 ¦È sin .v . sin2 ¦È sin .(u ¡¤ v)u. 
sin¦¸ sin¦¸ 
Using sin ¦È sin .(u ¡¤ v) = cos¦¸ . cos ¦È cos ., 
we get 
sin ¦Ë¦¸ sin ¦Ë¦¸ 
¦Í = . cos ¦È sin ¦È cos .u + cos ¦Ë¦¸ sin ¦Èu + sin .v 
sin¦¸ sin¦¸ sin ¦Ë¦¸ 
. sin ¦È(cos ¦¸ . cos ¦È cos .)u 
sin ¦¸ sin ¦Ë¦¸ sin ¦Ë¦¸ 
= cos ¦Ë¦¸ sin ¦Èu + sin .v . sin ¦È cos ¦¸ u 
sin¦¸ sin¦¸ (cos ¦Ë¦¸ sin ¦¸ . sin ¦Ë¦¸ cos ¦¸) sin ¦Ë¦¸ 
= sin ¦Èu + sin .v 
sin¦¸ sin¦¸ sin(1 . ¦Ë)¦¸ sin ¦Ë¦¸ 
= sin ¦Èu + sin . v. 
sin¦¸ sin¦¸ 
Putting the scalar part and the vector part together, we obtain 
.1 sin(1 . ¦Ë)¦¸ sin ¦Ë¦¸ 
q1(q1 q2)¦Ë = cos ¦È + cos .,
sin¦¸ sin¦¸ sin(1 . ¦Ë)¦¸ sin ¦Ë¦¸ 
sin ¦Èu + sin .v ,
sin¦¸ sin¦¸ sin(1 . ¦Ë)¦¸ sin ¦Ë¦¸ 
= [cos ¦È, sin ¦Èu] + [cos ., sin .v]. 
sin¦¸ sin¦¸ 
This yields the celebrated slerp interpolation formula 

.1 sin(1 . ¦Ë)¦¸ sin ¦Ë¦¸ 
Z(¦Ë)= q1(q1 q2)¦Ë = q1 + q2,
sin¦¸ sin¦¸ 
with cos¦¸ = cos ¦È cos . + sin ¦È sin .(u ¡¤ v). 


15.7 	Nonexistence of a ¡°Nice¡± Section from SO(3) to SU(2) 
We conclude by discussing the problem of a consistent choice of sign for the quaternion q representing a rotation R = ¦Ñq ¡Ê SO(3). We are looking for a ¡°nice¡± section s: SO(3) ¡ú SU(2), that is, a function s satisfying the condition 
¦Ñ . s = id, 
where ¦Ñ is the surjective homomorphism ¦Ñ: SU(2) ¡ú SO(3). 
Proposition 15.8. Any section s: SO(3) ¡ú SU(2) of ¦Ñ is neither a homomorphism nor continuous. 
Intuitively, this means that there is no ¡°nice and simple ¡± way to pick the sign of the quaternion representing a rotation. The following proof is due to Marcel Berger. 
Proof. Let ¦£ be the subgroup of SU(2) consisting of all quaternions of the form q = [a, (b, 0, 0)]. Then, using the formula for the rotation matrix Rq corresponding to q (and the fact that a2 + b2 = 1), we get 
.	. 
10 0 Rq = .02a2 . 1 .2ab . . 02ab 2a2 . 1 
Since a2 + b2 = 1, we may write a = cos ¦È, b = sin ¦È, and we see that 
.	. 
10 0 
.	.
Rq = 	0 cos2¦È . sin 2¦È, 0 sin 2¦È cos 2¦È 
a rotation of angle 2¦È around the x-axis. Thus, both ¦£ and its image are isomorphic to SO(2), which is also isomorphic to U(1) = {w ¡Ê C ||w| =1}. By identifying i and i, and identifying ¦£ and its image to U(1), if we write w = cos ¦È + i sin ¦È ¡Ê ¦£, the restriction of the map ¦Ñ to ¦£ is given by ¦Ñ(w)= w2 . 
15.7. NONEXISTENCE OF A ¡°NICE¡± SECTION FROM SO(3) TO SU(2) 
We claim that any section s of ¦Ñ is not a homomorphism. Consider the restriction of s to U(1). Then since ¦Ñ . s = id and ¦Ñ(w)= w2, for .1 ¡Ê ¦Ñ(¦£) ¡Ö U(1), we have 
.1= ¦Ñ(s(.1)) = (s(.1))2 . 
On the other hand, if s is a homomorphism, then 
(s(.1))2 = s((.1)2)= s(1) = 1, 
contradicting (s(.1))2 = .1. 
We also claim that s is not continuous. Assume that s(1) = 1, the case where s(1) = .1 being analogous. Then s is a bijection inverting ¦Ñ on ¦£ whose restriction to U(1) must be given by 
s(cos ¦È + i sin ¦È) = cos(¦È/2) + i sin(¦È/2), .¦Ð ¡Ü ¦È < ¦Ð. 
If ¦È tends to ¦Ð, that is z = cos ¦È + i sin ¦È tends to .1 in the upper-half plane, then s(z) tends to i, but if ¦È tends to .¦Ð, that is z tends to .1 in the lower-half plane, then s(z) tends to .i, which shows that s is not continuous. 
Another way (due to Jean Dieudonn¡äe) to prove that a section s of ¦Ñ is not a homomor-phism is to prove that any unit quaternion is the product of two unit pure quaternions. Indeed, if q =[a, u] is a unit quaternion, if we let q1 = [0,u1], where u1 is any unit vector orthogonal to u, then 
q1q =[.u1 ¡¤ u, au1 + u1 ¡Á u] = [0, au1 + u1 ¡Á u]= q2 
is a nonzero unit pure quaternion. This is because if a = 0 then au1 +u1 ¡Áu = 0 (since u1 ¡Áu is orthogonal to au1 = 0), and if a = 0 then u = 0, so u1 ¡Á u = 0 (since u1 is orthogonal to u). But then, q1 .1 = [0, .u1] is a unit pure quaternion and we have 
.1 
q = qq2,
1 
a product of two pure unit quaternions. 
We also observe that for any two pure quaternions q1,q2, there is some unit quaternion q such that 
.1 
q2 = qq1q. This is just a restatement of the fact that the group SO(3) is transitive. Since the kernel of ¦Ñ: SU(2) ¡ú SO(3) is {I, .I}, the subgroup s(SO(3)) would be a normal subgroup of index 2 in SU(2). Then we would have a surjective homomorphism ¦Ç from SU(2) onto the quotient group SU(2)/s(SO(3)), which is isomorphic to {1, .1}. Now, since any two pure quaternions are conjugate of each other, ¦Ç would have a constant value on the unit pure quaternions. Since k = ij, we would have 
¦Ç(k)= ¦Ç(ij)=(¦Ç(i))2 =1. 
Consequently, ¦Ç would map all pure unit quaternions to 1. But since every unit quaternion is the product of two pure quaternions, ¦Ç would map every unit quaternion to 1, contradicting the fact that it is surjective onto {.1, 1}. 


15.8 Summary 
The main concepts and results of this chapter are listed below: 
. 
The group SU(2) of unit quaternions. 

. 
The skew .eld H of quaternions. 

. 
Hamilton¡¯s identities. 

. 
The (real) vector space su(2) of 2 ¡Á 2 skew Hermitian matrices with zero trace. 

. 
The adjoint representation of SU(2). 

. 
The (real) vector space su(2) of 2 ¡Á 2 Hermitian matrices with zero trace. 

. 
The group homomorphism r : SU(2) ¡ú SO(3); Ker (r)= {+I, .I}. 

. 
The matrix representation Rq of the rotation rq induced by a unit quaternion q. 

. 
Surjectivity of the homomorphism r : SU(2) ¡ú SO(3). 

. 
The exponential map exp: su(2) ¡ú SU(2). 

. 
Surjectivity of the exponential map exp: su(2) ¡ú SU(2). 

. 
Finding a logarithm of a quaternion. 

. 
Quaternion interpolation. 

. 
Shoemake¡¯s slerp interpolation formula. 

. 
Sections s: SO(3) ¡ú SU(2) of r : SU(2) ¡ú SO(3). 



15.9 Problems 
Problem 15.1. Verify the quaternion identities i2 = j2 = k2 = ijk = .1, ij = .ji = k, jk = .kj = i, ki = .ik = j. Problem 15.2. Check that for every quaternion X = a1 + bi + cj + dk, we have 2 + b22 + d2)1.
XX . = X . X =(a + c Conclude that if X = 0, then X is invertible and its inverse is given by 22 + d2).1X . 
X.1 =(a + b2 + c. 
15.9. PROBLEMS 
Problem 15.3. Given any two quaternions X = a1+bi+cj+dk and Y = a ' 1+b ' i+c ' j+d ' k, prove that 
XY =(aa ' . bb ' . cc ' . dd ' )1 +(ab ' + ba ' + cd ' . dc ' )i 
+(ac ' + ca ' + db ' . bd ' )j +(ad ' + da ' + bc ' . cb ' )k. 
Also prove that if X =[a, U] and Y =[a ' ,U ' ], the quaternion product XY can be expressed as XY =[aa ' . U ¡¤ U ' , aU ' + a ' U + U ¡Á U ' ]. 
Problem 15.4. Let Ad: SU(2) ¡ú GL(su(2)) be the map de.ned such that for every q ¡Ê SU(2), 
Adq(A)= qAq . ,A ¡Ê su(2), where q . is the inverse of q (since SU(2) is a unitary group) Prove that the map Adq is an invertible linear map from su(2) to itself and that Ad is a group homomorphism. 
Problem 15.5. Prove that every Hermitian matrix with zero trace is of the form x¦Ò3 + 
y¦Ò2 + z¦Ò1, with  
¦Ò1 =  0 1  1 0  ,  ¦Ò2 =  0 i  .i 0  ,  ¦Ò3 =  1 0  0 .1  .  
Check that i = i¦Ò3, j = i¦Ò2, and that k = i¦Ò1.  
Problem 15.6. If  .  .  
0  .u3  u2  
B = . u3  0  .u1. ,  
.u2  u1  0  

222
and if we write ¦È = u1 + u2 + u3 (with 0 ¡Ü ¦È ¡Ü ¦Ð), then the Rodrigues formula says that B sin ¦È (1 . cos ¦È)
B2 
e = I + B + ,¦È =0,
¦È¦È2 with e0 = I. Check that tr(eB) = 1+2cos ¦È. Prove that the quaternion q corresponding to the rotation R = eB (with B = 0) is given by 
¦È ¦Èu1 u2 u3 
q = cos , sin ,, . 
22 ¦È¦È¦È Problem 15.7. For every matrix A ¡Ê su(2), with 
iu1 u2 + iu3
A = ,
.u2 + iu3 .iu1 
222
prove that if we write ¦È = u, then 
1 + u2 + u3sin ¦È 
e A = cos ¦ÈI + A, ¦È =0,
¦È and e0 = I. Conclude that eA is a unit quaternion representing the rotation of angle 2¦È and axis (u1,u2,u3) (or I when ¦È = k¦Ð, k ¡Ê Z). 
Problem 15.8. Write a Matlab program implementing the method of Section 15.4 for .nding a unit quaternion corresponding to a rotation matrix. 
Problem 15.9. Show that there is a very simple method for producing an orthonormal frame in R4 whose .rst vector is any given nonnull vector (a, b, c, d). 
Problem 15.10. Let i, j, and k, be the unit vectors of coordinates (1, 0, 0), (0, 1, 0), and (0, 0, 1) in R3 . 
(1) Describe geometrically the rotations de.ned by the following quaternions: p = (0,i),q = (0,j). Prove that the interpolant Z(¦Ë)= p(p.1q)¦Ë is given by 
Z(¦Ë) = (0, cos(¦Ë¦Ð/2)i + sin(¦Ë¦Ð/2)j) . Describe geometrically what this rotation is. 
(2) Repeat Question (1) with the rotations de.ned by the quaternions 
  
¡Ì 
13 
p =, i,q = (0,j). 
22 
Prove that the interpolant Z(¦Ë) is given by 
  
¡Ì 
13 
Z(¦Ë) =cos(¦Ë¦Ð/2), cos(¦Ë¦Ð/2)i + sin(¦Ë¦Ð/2)j. 
22 Describe geometrically what this rotation is. 
(3) Repeat Question (1) with the rotations de.ned by the quaternions 
11 1 

p = ¡Ì , ¡Ì i, q =0, ¡Ì (i + j) . 
22 2Prove that the interpolant Z(¦Ë) is given by 
11 
Z(¦Ë)= ¡Ì cos(¦Ë¦Ð/3) .¡Ì sin(¦Ë¦Ð/3), 
26 
¡Ì¡Ì 
2 
(1/ 2 cos(¦Ë¦Ð/3) + 1/ 6 sin(¦Ë¦Ð/3))i + ¡Ì sin(¦Ë¦Ð/3)j. 
6 Problem 15.11. Prove that w ¡Á (u ¡Á v)=(w ¡¤ v)u . (u ¡¤ w)v. Conclude that u ¡Á (u ¡Á v)=(u ¡¤ v)u . (u ¡¤ u)v. 



