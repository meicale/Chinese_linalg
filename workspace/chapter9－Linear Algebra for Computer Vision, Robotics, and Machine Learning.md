Chapter 9 
Iterative Methods for Solving Linear Systems 
9.1 Convergence of Sequences of Vectors and Matrices 
In Chapter 7 we discussed some of the main methods for solving systems of linear equations. These methods are direct methods, in the sense that they yield exact solutions (assuming in.nite precision!). 
Another class of methods for solving linear systems consists in approximating solutions using iterative methods. The basic idea is this: Given a linear system Ax = b (with A a square invertible matrix in Mn(C)), .nd another matrix B �� Mn(C) and a vector c �� Cn , such that 
1. 
The matrix I . B is invertible 

2. 
The unique solution xuof the system Ax = b is identical to the unique solution uuof the 


system 
u = Bu + c, 

and then starting from any vector u0, compute the sequence (uk) given by 
uk+1 = Buk + c, k �� N. 
Under certain conditions (to be clari.ed soon), the sequence (uk) converges to a limit uuwhich is the unique solution of u = Bu + c, and thus of Ax = b. 
Consequently, it is important to .nd conditions that ensure the convergence of the above sequences and to have tools to compare the ��rate�� of convergence of these sequences. Thus, we begin with some general results about the convergence of sequences of vectors and ma-trices. 
Let (E, 11) be a normed vector space. Recall from Section 8.7 that a sequence (uk) of vectors uk �� E converges to a limit u �� E, if for every��> 0, there some natural number N such that 
1uk . u1�ܷ�, for all k �� N. 
319 

We write u = lim uk. 
kh���� 
If E is a .nite-dimensional vector space and dim(E)= n, we know from Theorem 8.5 that any two norms are equivalent, and if we choose the norm 11��, we see that the convergence of the sequence of vectors uk is equivalent to the convergence of the n sequences of scalars formed by the components of these vectors (over any basis). The same property applies to the .nite-dimensional vector space Mm,n(K) of m �� n matrices (with K = R or K = C), 
(k)
which means that the convergence of a sequence of matrices Ak =(a) is equivalent to the 
ij 
convergence of the m �� n sequences of scalars (a(k)), with i, j .xed (1 �� i �� m,1 �� j �� n).
ij The .rst theorem below gives a necessary and su.cient condition for the sequence (Bk) of powers of a matrix B to converge to the zero matrix. Recall that the spectral radius ��(B) of a matrix B is the maximum of the moduli |��i| of the eigenvalues of B. 
Theorem 9.1. For any square matrix B, the following conditions are equivalent: 
(1) limkh���� Bk =0, 

(2) limkh���� Bkv =0, for all vectors v, 

(3) ��(B) < 1, 

(4) 1B1 < 1, for some subordinate matrix norm 11. 


Proof. Assume (1) and let 11 be a vector norm on E and 11 be the corresponding matrix norm. For every vector v �� E, because 11 is a matrix norm, we have 
1Bk v1��1Bk11v1, 
and since limkh���� Bk = 0 means that limkh���� 1Bk1 = 0, we conclude that limkh���� 1Bkv1 = 0, that is, limkh���� Bkv = 0. This proves that (1) implies (2). 
Assume (2). If we had ��(B) �� 1, then there would be some eigenvector u ( 
= 0) and some eigenvalue �� such that Bu = ��u, |��| = ��(B) �� 1, 
but then the sequence (Bku) would not converge to 0, because Bku = ��ku and |��k| = |��|k �� 
1. It follows that (2) implies (3). 
Assume that (3) holds, that is, ��(B) < 1. By Proposition 8.12, we can .nd��> 0 small enough that ��(B)+��< 1, and a subordinate matrix norm 11 such that 
1B1�� ��(B)+��, 
which is (4). Finally, assume (4). Because 11 is a matrix norm, 
1Bk1��1B1k , 
and since 1B1 < 1, we deduce that (1) holds. 

9.2. CONVERGENCE OF ITERATIVE METHODS 
The following proposition is needed to study the rate of convergence of iterative methods. 
Proposition 9.2. For every square matrix B �� Mn(C) and every matrix norm 11, we have 
1Bk11/k
lim = ��(B). 
kh���� 
Proof. We know from Proposition 8.6 that ��(B) ��1B1, and since ��(B)=(��(Bk))1/k, we deduce that 
��(B) ��1Bk11/k for all k �� 1, 
and so 
1Bk11/k
��(B) �� lim . 
kh���� 
Now let us prove that for every��> 0, there is some integer N(��) such that 
1Bk11/k �� ��(B)+�� for all k �� N(��), 
which proves that 
1Bk11/k �� ��(B),
lim 
kh���� 
and our proposition. For any given��> 0, let B�� be the matrix 
B B�� = ��(B)+��. 
Since 1B��1 < 1, Theorem 9.1 implies that limkh���� B��k= 0. Consequently, there is some integer N(��) such that for all k �� N(��), we have 
1Bk1 
1Bk1 =(��(B)+��)k �� 1, 
which implies that 1Bk11/k �� ��(B)+��, 
as claimed. 
We now apply the above results to the convergence of iterative methods. 
9.2 Convergence of Iterative Methods 
Recall that iterative methods for solving a linear system Ax = b (with A �� Mn(C) invertible) consists in .nding some matrix B and some vector c, such that I . B is invertible, and the unique solution xuof Ax = b is equal to the unique solution u 
of u = Bu + c. Then starting from any vector u0, compute the sequence (uk) given by 
uk+1 = Buk + c, k �� N, 
and say that the iterative method is convergent i. 
lim uk u,
= u
kh���� 
for every initial vector u0. 
Here is a fundamental criterion for the convergence of any iterative methods based on a matrix B, called the matrix of the iterative method. 
Theorem 9.3. Given a system u = Bu + c as above, where I . B is invertible, the following statements are equivalent: 
(1) 
The iterative method is convergent. 

(2) 
��(B) < 1. 


(3) 1B1 < 1, for some subordinate matrix norm 11. Proof. De.ne the vector ek (error vector) by 
u,
ek = uk . u
where u 
is the unique solution of the system u = Bu + c. Clearly, the iterative method is convergent i. 
lim ek =0. 
kh���� 
We claim that ek = Bk e0,k �� 0, 
where e0 = u0 . u
. 
This is proven by induction on k. The base case k = 0 is trivial. By the induction hypothesis, ek = Bke0, and since uk+1 = Buk + c, we get 
u,
uk+1 . u 
= Buk + c . u
and because u 
= Bu 
+ c and ek = Bke0 (by the induction hypothesis), we obtain = Bk+1 
uk+1 . u 
= Buk . Bu 
= B(uk . u
)= Bek = BBk e0 e0, 
proving the induction step. Thus, the iterative method converges i. 
Bk
lim e0 =0. 
kh���� 
Consequently, our theorem follows by Theorem 9.1. 
The next proposition is needed to compare the rate of convergence of iterative methods. It shows that asymptotically, the error vector ek = Bke0 behaves at worst like (��(B))k . 
9.2. CONVERGENCE OF ITERATIVE METHODS 
Proposition 9.4. Let 11 be any vector norm, let B �� Mn(C) be a matrix such that I . B is invertible, and let u 
be the unique solution of u = Bu + c. 
(1) If (uk) is any sequence de.ned iteratively by 
uk+1 = Buk + c, k �� N, 
then   
1uk . u
11/k
lim sup= ��(B). 
kh����ߵu0.һuߵ=1 
(2) Let B1 and B2 be two matrices such that I . B1 and I . B2 are invertibe, assume that both u = B1u + c1 and u = B2u + c2 have the same unique solution u
, and consider any two sequences (uk) and (vk) de.ned inductively by 
uk+1 = B1uk + c1 
vk+1 = B2vk + c2, 
with u0 = v0. If ��(B1) <��(B2), then for any��> 0, there is some integer N(��), such that 

for all k �� N(��), we have  
sup  1vk . u
1  1/k ��  ��(B2)  . 

һߵ=11uk . u
1��(B1)+��
ߵu0.u
Proof. Let 11 be the subordinate matrix norm. Recall that 
uk . u 
= Bk e0, with e0 = u0 . u
. For every k �� N, we have (��(B1))k = ��(B1 k) ��1B1 k1 = sup1B1 k e01, 
ߵe0ߵ=1 
which implies 
1 e011/k 11/k
��(B1) = sup1Bk = 1B1 k , 
ߵe0ߵ=1 
and Statement (1) follows from Proposition 9.2. Because u0 = v0, we have uk . u 
= B1 k e0 
vk . u 
= B2 k e0, with e0 = u0 . u 
= v0 . u
. Again, by Proposition 9.2, for every��> 0, there is some natural number N(��) such that if k �� N(��), then 
sup1B1 k e011/k �� ��(B1)+��. 
ߵe0ߵ=1 Furthermore, for all k �� N(��), there exists a vector e0 = e0(k) such that 2 e011/k 11/k �� ��(B2),
1e01 = 1 and 1Bk = 1B2 k 
which implies Statement (2). 

In light of the above, we see that when we investigate new iterative methods, we have to deal with the following two problems: 
1. 
Given an iterative method with matrix B, determine whether the method is conver-gent. This involves determining whether ��(B) < 1, or equivalently whether there is a subordinate matrix norm such that 1B1 < 1. By Proposition 8.11, this implies that I . B is invertible (since 1. B1 = 1B1, Proposition 8.11 applies). 

2. 
Given two convergent iterative methods, compare them. The iterative method which is faster is that whose matrix has the smaller spectral radius. 


We now discuss three iterative methods for solving linear systems: 
1. 
Jacobi��s method 

2. 
Gauss�CSeidel��s method 

3. 
The relaxation method. 


9.3 	Description of the Methods of Jacobi, Gauss�CSeidel, and Relaxation 
The methods described in this section are instances of the following scheme: Given a linear system Ax = b, with A invertible, suppose we can write A in the form 
A = M . N, 
with M invertible, and ��easy to invert,�� which means that M is close to being a diagonal or a triangular matrix (perhaps by blocks). Then Au = b is equivalent to 
Mu = Nu + b, 
that is, u = M.1Nu + M.1b. 
Therefore, we are in the situation described in the previous sections with B = M.1N and c = M.1b. In fact, since A = M . N, we have 
B = M.1N = M.1(M . A)= I . M.1A, 	(.) 
which shows that I . B = M.1A is invertible. The iterative method associated with the matrix B = M.1N is given by 
uk+1 = M.1Nuk + M.1b, k �� 0, 	(.) 
9.3. METHODS OF JACOBI, GAUSS�CSEIDEL, AND RELAXATION 
starting from any arbitrary vector u0. From a practical point of view, we do not invert M, and instead we solve iteratively the systems 
Muk+1 = Nuk + b, k �� 0. 
Various methods correspond to various ways of choosing M and N from A. The .rst two methods choose M and N as disjoint submatrices of A, but the relaxation method allows some overlapping of M and N. 
To describe the various choices of M and N, it is convenient to write A in terms of three submatrices D, E, F , as 
A = D . E . F, 
where the only nonzero entries in D are the diagonal entries in A, the only nonzero entries in E are entries in A below the the diagonal, and the only nonzero entries in F are entries in A above the diagonal. More explicitly, if 
.
. 
a11 a12 a13 ������ a1n.1 a1n a21 a22 a23 
������ a2n.1 a2n a31 a32 a33 ������ a3n.1 a3n 
. . .  . . .  . . .  .. .  . . .  . . .  
an.1 1  an.1 2  an.1 3  �� �� ��  an.1 n.1  an.1 n  
an 1  an 2  an 3  �� �� ��  an n.1  an n  

............ 

............ 

,

A = 

then 

.
. 
............ 

00 ������ 00 00 ������ 00 
a11 a22 
00 ������ 00
a33 
... ..
.
..... .
.
... .. 
000 ������ an.1 n.1 0 000 ������ 0 
ann 
............ 

D = 

, 

.
. 
000 ������ 00 0
............ 
00 ������ 0

a21 
............ 

0 .. ..
..
. . ....
..
.. .. 
0 ������ 0

a31 a32 
.E = 

, 

..
.0 0 
������ 0
an 1 an 2 an 3 ann.1 
an.11 an.12 an.13 
.
. 
............ 

0  a12  a13  �� �� ��  a1n.1  a1n  
0  0  a23  �� �� ��  a2n.1  a2n  
0  0  0  ...  a3n.1  a3n  

. . .  . . .  . . .  ...  ...  . . .  
0  0  0  �� �� ��  0  an.1 n  
0  0  0  �� �� ��  0  0  

............ 

.F = 

. 

In Jacobi��s method, we assume that all diagonal entries in A are nonzero, and we pick M = D N = E + F, so that by (.), B = M.1N = D.1(E + F )= I . D.1A. As a matter of notation, we let J = I . D.1A = D.1(E + F ), which is called Jacobi��s matrix. The corresponding method, Jacobi��s iterative method, com-putes the sequence (uk) using the recurrence 
uk+1 = D.1(E + F )uk + D.1b, k �� 0. 
In practice, we iteratively solve the systems Duk+1 =(E + F )uk + b, k �� 0. 
If we write uk =(u1k,...,ukn), we solve iteratively the following system: k+1 kk
a11u = .a12u������ .a1nu+ b1
12 n k+1 
a22u2 = .a21u1 k ������ .a2nunk + b2 
. .. 
. .. 
.
. .. 
k+1 kk
an.1 n.1un.1 = .an.11u1 ������ .an.1 nun + bn.1 k+1 kk k
annun = .an 1u1 .an 2u2 .ann.1un.1 + bn 
In Matlab one step of Jacobi iteration is achieved by the following function: 

9.3. METHODS OF JACOBI, GAUSS�CSEIDEL, AND RELAXATION 
function v = Jacobi2(A,b,u) n = size(A,1); v = zeros(n,1); 
for i = 1:n v(i,1) = u(i,1) + (-A(i,:)*u + b(i))/A(i,i); end end 
In order to run m iteration steps, run the following function: 
function u = jacobi(A,b,u0,m) u = u0; for j= 1:m 
u = Jacobi2(A,b,u); end end 
Example 9.1. Consider the linear system 
.
....
. 
2 .10 0 x1 25 
... 

.12 .1 

0 

0 .12 .1 

... 
... 

x2 
x3 
... 

= 

... 

.24 

21 

...

. 

00 .12 x4 .15 We check immediately that the solution is 
x1 = 11,x2 = .3,x3 =7,x4 = .4. It is easy to see that the Jacobi matrix is 
.
. 
0100 

J = 

1 

2 

... 

1010 

0101 

...

. 

0010 After 10 Jacobi iterations, we .nd the approximate solution x1 = 10.2588,x2 = .2.5244,x3 =5.8008,x4 = .3.7061. After 20 iterations, we .nd the approximate solution x1 = 10.9110,x2 = .2.9429,x3 =6.8560,x4 = .3.9647. After 50 iterations, we .nd the approximate solution 
x1 = 10.9998,x2 = .2.9999,x3 =6.9998,x4 = .3.9999, 
and After 60 iterations, we .nd the approximate solution x1 = 11.0000,x2 = .3.0000,x3 =7.0000,x4 = .4.0000, correct up to at least four decimals. It can be shown (see Problem 9.6) that the eigenvalues of J are 
  
�� 2�� 3�� 4�� 
cos, cos , cos , cos ,
55 5 5 so the spectral radius of J = B is 
  
�� 
��(J)=cos=0.8090 < 1. 
5By Theorem 9.3, Jacobi��s method converges for the matrix of this example. Observe that we can try to ��speed up�� the method by using the new value u1 k+1 instead kk+2 k+1 k+1
of uin solving for u using the second equations, and more generally, use u ,...,u 
12 1 i.1 kk k+1
instead of u1,...,uin solving for u in the ith equation. This observation leads to the 
i.1 i 
system 
k+1 kk
= ������ .a1nu+ b1
a11u1 .a12u2 n 
k+1 k+1 k
= ������ .a2nu+ b2
a22u2 .a21u1 n 
... 
... ,
... 
k+1 k+1 �� k
an.1 n.1un.1 = .an.11u1 ���� .an.1 nun + bn.1 k+1 k+1 k+1 k+1 
annun = .an 1u1 .an 2u2 .ann.1un.1 + bn 
which, in matrix form, is written Duk+1 = Euk+1 + Fuk + b. Because D is invertible and E is lower triangular, the matrix D . E is invertible, so the above equation is equivalent to uk+1 =(D . E).1Fuk +(D . E).1b, k �� 0. The above corresponds to choosing M and N to be M = D . E N = F, and the matrix B is given by B = M.1N =(D . E).1F. Since M = D . E is invertible, we know that I . B = M.1A is also invertible. 
9.3. METHODS OF JACOBI, GAUSS�CSEIDEL, AND RELAXATION 
The method that we just described is the iterative method of Gauss�CSeidel, and the matrix B is called the matrix of Gauss�CSeidel and denoted by L1, with 
L1 =(D . E).1F. 
One of the advantages of the method of Gauss�CSeidel is that is requires only half of the memory used by Jacobi��s method, since we only need 
k+1 k+1 kk 
u ,...,u i.1 ,u i+1,...,u 
1 n 
to compute uik+1 . We also show that in certain important cases (for example, if A is a tridiagonal matrix), the method of Gauss�CSeidel converges faster than Jacobi��s method (in this case, they both converge or diverge simultaneously). 
In Matlab one step of Gauss�CSeidel iteration is achieved by the following function: 
function u = GaussSeidel3(A,b,u) 
n = size(A,1); 
fori =1:n 

u(i,1) = u(i,1) + (-A(i,:)*u + b(i))/A(i,i); end end 
It is remarkable that the only di.erence with Jacobi2 is that the same variable u is used on both sides of the assignment. In order to run m iteration steps, run the following function: 
function u = GaussSeidel1(A,b,u0,m) u = u0; for j= 1:m 
u = GaussSeidel3(A,b,u); end end 
Example 9.2. Consider the same linear system 
.
....
. 
2 .10 0 x1 25 
.12 .1 

0 

0 .12 .1 

... 
... 

x2 
x3 
... 

= 

... 

.24 

21 

00 .12 x4 .15 as in Example 9.1, whose solution is x1 = 11,x2 = .3,x3 =7,x4 = .4. After 10 Gauss�CSeidel iterations, we .nd the approximate solution 
x1 = 10.9966,x2 = .3.0044,x3 =6.9964,x4 = .4.0018. 
After 20 iterations, we .nd the approximate solution 
x1 = 11.0000,x2 = .3.0001,x3 =6.9999,x4 = .4.0000. 
After 25 iterations, we .nd the approximate solution 
x1 = 11.0000,x2 = .3.0000,x3 =7.0000,x4 = .4.0000, 
correct up to at least four decimals. We observe that for this example, Gauss�CSeidel��s method converges about twice as fast as Jacobi��s method. It will be shown in Proposition 9.8 that for a tridiagonal matrix, the spectral radius of the Gauss�CSeidel matrix L1 is given by 
��(L1)=(��(J))2 , 
so our observation is consistent with the theory. 
The new ingredient in the relaxation method is to incorporate part of the matrix D into 
N: we de.ne M and N by D 
M = . E 
�� 1 . �� 
N = D + F, 
�� 
where �� = 0 is a real parameter to be suitably chosen. Actually, we show in Section 9.4 that for the relaxation method to converge, we must have �� �� (0, 2). Note that the case �� =1 corresponds to the method of Gauss�CSeidel. 
If we assume that all diagonal entries of D are nonzero, the matrix M is invertible. The matrix B is denoted by L�� and called the matrix of relaxation, with 
D .1 1 . �� 
L�� = . ED + F =(D . ��E).1((1 . ��)D + ��F ). 
�ئ� 
The number �� is called the parameter of relaxation. 
When ��> 1, the relaxation method is known as successive overrelaxation, abbreviated as SOR. 
At .rst glance the relaxation matrix L�� seems at lot more complicated than the Gauss�C Seidel matrix L1, but the iterative system associated with the relaxation method is very similar to the method of Gauss�CSeidel, and is quite simple. Indeed, the system associated with the relaxation method is given by 
D 1 . �� 
. Euk+1 = D + Fuk + b,
�ئ� 
which is equivalent to 
(D . ��E)uk+1 = ((1 . ��)D + ��F )uk + ��b, 
9.3. METHODS OF JACOBI, GAUSS�CSEIDEL, AND RELAXATION 
and can be written 
Duk+1 = Duk . ��(Duk . Euk+1 . Fuk . b). 
Explicitly, this is the system 
k+1 kk 	kk 
��	. b1)
a11u	1 = a11u1 . ��(a11u1 + ���� + a1n.1un.1 + a1nun k+1 kk+1 kk 
a22u2 = a22u2 . ��(a21u1 + ������ + a2n.1un.1 + a2nun . b2) . 
. 
. 
k+1 kk+1 	k+1 k 
annun = annun . ��(an 1u1 ++ ������ + ann.1un.1 + annun . bn). 
In Matlab one step of relaxation iteration is achieved by the following function: 
function u = relax3(A,b,u,omega) 
n = size(A,1); 
fori =1:n 

u(i,1) = u(i,1) + omega*(-A(i,:)*u + b(i))/A(i,i); end end 
Observe that function relax3 is obtained from the function GaussSeidel3 by simply insert-ing �� in front of the expression (.A(i, :) . u + b(i))/A(i, i). In order to run m iteration steps, run the following function: 
function u = relax(A,b,u0,omega,m) u = u0; for j= 1:m 
u = relax3(A,b,u,omega); end end 
Example 9.3. Consider the same linear system as in Examples 9.1 and 9.2, whose solution is x1 = 11,x2 = .3,x3 =7,x4 = .4. 
After 10 relaxation iterations with �� =1.1, we .nd the approximate solution 
x1 = 11.0026,x2 = .2.9968,x3 =7.0024,x4 = .3.9989. 
After 10 iterations with �� =1.2, we .nd the approximate solution 
x1 = 11.0014,x2 = .2.9985,x3 =7.0010,x4 = .3.9996. 
After 10 iterations with �� =1.3, we .nd the approximate solution 
x1 = 10.9996,x2 = .3.0001,x3 =6.9999,x4 = .4.0000. 
After 10 iterations with �� =1.27, we .nd the approximate solution x1 = 11.0000,x2 = .3.0000,x3 =7.0000,x4 = .4.0000, 
correct up to at least four decimals. We observe that for this example the method of relax-ation with �� =1.27 converges faster than the method of Gauss�CSeidel. This observation will be con.rmed by Proposition 9.10. 
What remains to be done is to .nd conditions that ensure the convergence of the relax-ation method (and the Gauss�CSeidel method), that is: 
1. 
Find conditions on ��, namely some interval I . R so that �� �� I implies ��(L��) < 1; we will prove that �� �� (0, 2) is a necessary condition. 

2. 
Find if there exist some optimal value ��0 of �� �� I, so that 


��(L��0 ) = inf ��(L��). 
�ء�I 
We will give partial answers to the above questions in the next section. 
It is also possible to extend the methods of this section by using block decompositions of the form A = D . E . F , where D, E, and F consist of blocks, and D is an invertible block-diagonal matrix. See Figure 9.1. 

Figure 9.1: A schematic representation of a block decomposition A = D . E . F , where D = ��4 = ��3 i=1Ei, and F = ��3 
i=1Di, E i=1Fi. 
9.4. CONVERGENCE OF THE METHODS 
9.4 	Convergence of the Methods of Gauss�CSeidel and Relaxation 
We begin with a general criterion for the convergence of an iterative method associated with a (complex) Hermitian positive de.nite matrix, A = M . N. Next we apply this result to the relaxation method. 
Proposition 9.5. Let A be any Hermitian positive de.nite matrix, written as A = M . N, with M invertible. Then M. + N is Hermitian, and if it is positive de.nite, then ��(M.1N) < 1, 
so that the iterative method converges. Proof. Since M = A + N and A is Hermitian, A. = A, so we get M . + N = A . + N . + N = A + N + N . = M + N . =(M . + N) . , which shows that M. + N is indeed Hermitian. Because A is Hermitian positive de.nite, the function 
v �� (v . Av)1/2 
from Cn to R is a vector norm 11, and let 11 also denote its subordinate matrix norm. We prove that 
1M.1N1 < 1, which by Theorem 9.1 proves that ��(M.1N) < 1. By de.nition 1M.1N1 = 1I . M.1A1 = sup1v . M.1Av1, 
ߵvߵ=1 
which leads us to evaluate 1v . M.1Av1 when 1v1 = 1. If we write w = M.1Av, using the facts that 1v1 = 1, v = A.1Mw, A. = A, and A = M . N, we have 
1v . w12 =(v . w) . A(v . w) = 1v12 . v . Aw . w . Av + w . Aw . M . 
=1 . ww . w . Mw + w . Aw =1 . w . (M . + N)w. 
Now since we assumed that M. + N is positive de.nite, if w = 0, then w .(M. + N)w> 0, and we conclude that 
if 1v1 =1, then 1v . M.1Av1 < 1. 
Finally, the function 
v ��1v . M.1Av1 
is continuous as a composition of continuous functions, therefore it achieves its maximum on the compact subset {v �� Cn |1v1 =1}, which proves that 
sup1v . M.1Av1 < 1, 
ߵvߵ=1 
and completes the proof. 
Now as in the previous sections, we assume that A is written as A = D . E . F , with D invertible, possibly in block form. The next theorem provides a su.cient condition (which turns out to be also necessary) for the relaxation method to converge (and thus, for the method of Gauss�CSeidel to converge). This theorem is known as the Ostrowski-Reich theorem. 
Theorem 9.6. If A = D . E . F is Hermitian positive de.nite, and if 0 <��< 2, then the relaxation method converges. This also holds for a block decomposition of A. 
Proof. Recall that for the relaxation method, A = M . N with D 
M = . E 
�� 1 . �� 
N = D + F, 
�� and because D. = D, E. = F (since A is Hermitian) and �� = 0 is real, we have D. 1 . �� 2 . �� 
M . . E . 
+ N =+ D + F = D. 
�ئ� �� If D consists of the diagonal entries of A, then we know from Section 7.8 that these entries are all positive, and since �� �� (0, 2), we see that the matrix ((2 .��)/��)D is positive de.nite. If D consists of diagonal blocks of A, because A is positive, de.nite, by choosing vectors z obtained by picking a nonzero vector for each block of D and padding with zeros, we see that each block of D is positive de.nite, and thus D itself is positive de.nite. Therefore, in all cases, M. + N is positive de.nite, and we conclude by using Proposition 9.5. 
Remark: What if we allow the parameter �� to be a nonzero complex number �� �� C? In this case, we get 
D. 1 . �� 11 
M . + N = . E . + D + F =+ . 1 D. 
�ئ� �ئ� 
But, 
11 �� + �� . �ئ� 1 . (�� . 1)(�� . 1) 1 .|�� . 1|2 
+ . 1= = = ,
�ئ� �ئ� |��|2 |��|2 
9.5. CONVERGENCE METHODS FOR TRIDIAGONAL MATRICES 
so the relaxation method also converges for �� �� C, provided that 
|�� . 1| < 1. 
This condition reduces to 0 <��< 2 if �� is real. 
Unfortunately, Theorem 9.6 does not apply to Jacobi��s method, but in special cases, Proposition 9.5 can be used to prove its convergence. On the positive side, if a matrix is strictly column (or row) diagonally dominant, then it can be shown that the method of Jacobi and the method of Gauss�CSeidel both converge. The relaxation method also converges if �� �� (0, 1], but this is not a very useful result because the speed-up of convergence usually occurs for ��> 1. 
We now prove that, without any assumption on A = D . E . F , other than the fact that A and D are invertible, in order for the relaxation method to converge, we must have �� �� (0, 2). 
Proposition 9.7. Given any matrix A = D . E . F , with A and D invertible, for any �� =0, we have 
��(L��) ��|�� . 1|, 
.1 D 1.��
where L�� = �� . E �� D + F . Therefore, the relaxation method (possibly by blocks) does not converge unless �� �� (0, 2). If we allow �� to be complex, then we must have 
|�� . 1| < 1 
for the relaxation method to converge. 
Proof. Observe that the product ��1 ������ ��n of the eigenvalues of L��, which is equal to det(L��), is given by 
1 . �� 
det  D + F  
��1 �� �� �� ��n = det(L��) =  �� D  = (1 . ��)n .  
det  . E  
��  
It follows that  

��(L��) ��|��1 ������ ��n|1/n = |�� . 1|. The proof is the same if �� �� C. 
9.5 	Convergence of the Methods of Jacobi, Gauss�CSeidel, and Relaxation for Tridiagonal Matrices 
We now consider the case where A is a tridiagonal matrix, possibly by blocks. In this case, we obtain precise results about the spectral radius of J and L��, and as a consequence, about the convergence of these methods. We also obtain some information about the rate of convergence of these methods. We begin with the case �� = 1, which is technically easier to deal with. The following proposition gives us the precise relationship between the spectral radii ��(J) and ��(L1) of the Jacobi matrix and the Gauss�CSeidel matrix. 
Proposition 9.8. Let A be a tridiagonal matrix (possibly by blocks). If ��(J) is the spectral radius of the Jacobi matrix and ��(L1) is the spectral radius of the Gauss�CSeidel matrix, then we have 
��(L1)=(��(J))2 . 
Consequently, the method of Jacobi and the method of Gauss�CSeidel both converge or both diverge simultaneously (even when A is tridiagonal by blocks); when they converge, the method of Gauss�CSeidel converges faster than Jacobi��s method. 
Proof. We begin with a preliminary result. Let A(��) with a tridiagonal matrix by block of 
the form 

. 

A1 ��.1C1 00 ������ 0 ��B1 A2 ��.1C2 0 ������ 0 
....
....
0. . . ������ . 
. 
...
. ...
. ������ ...0 0 ������ 0 ��Bp.2 Ap.1 ��.1Cp.1 0 ������ ������ 0 ��Bp.1 Ap 
. 

........ 

........ 

A(��)= 

, 

then det(A(��)) = det(A(1)),�� =0. 
To prove this fact, form the block diagonal matrix 
P (��) = diag(��I1,�� 2I2,...,��pIp), 
where Ij is the identity matrix of the same dimension as the block Aj. Then it is easy to see that 
A(��)= P (��)A(1)P (��).1 , 
and thus, det(A(��)) = det(P (��)A(1)P (��).1) = det(A(1)). 
Since the Jacobi matrix is J = D.1(E + F ), the eigenvalues of J are the zeros of the characteristic polynomial pJ (��) = det(��I . D.1(E + F )), 
and thus, they are also the zeros of the polynomial 
qJ (��) = det(��D . E . F ) = det(D)pJ (��). 
Similarly, since the Gauss�CSeidel matrix is L1 =(D . E).1F , the zeros of the characteristic polynomial 
pL1 (��) = det(��I . (D . E).1F ) 
9.5. CONVERGENCE METHODS FOR TRIDIAGONAL MATRICES 
are also the zeros of the polynomial 
qL1 (��) = det(��D . ��E . F ) = det(D . E)pL1 (��). 
Since A = D . E . F is tridiagonal (or tridiagonal by blocks), ��2D . ��2E . F is also tridiagonal (or tridiagonal by blocks), and by using our preliminary result with �� = �� = 0, we get 
qL1 (��2) = det(��2D . ��2E . F ) = det(��2D . ��E . ��F )= ��n qJ (��). 
By continuity, the above equation also holds for �� = 0. But then we deduce that: 
1. 
For any �� = 0, if �� is an eigenvalue of L1, then ��1/2 and .��1/2 are both eigenvalues of J, where ��1/2 is one of the complex square roots of ��. 

2. 
For any �� = 0, if �� and .�� are both eigenvalues of J, then ��2 is an eigenvalue of L1. 


The above immediately implies that ��(L1)=(��(J))2 . 
We now consider the more general situation where �� is any real in (0, 2). 
Proposition 9.9. Let A be a tridiagonal matrix (possibly by blocks), and assume that the eigenvalues of the Jacobi matrix are all real. If �� �� (0, 2), then the method of Jacobi and the method of relaxation both converge or both diverge simultaneously (even when A is tridiagonal by blocks). When they converge, the function �� �� ��(L��) (for �� �� (0, 2)) has a unique minimum equal to ��0 . 1 for 
��0 = 1+һ 1 . 2(��(J))2 , 
where 1 <��0 < 2 if ��(J) > 0. We also have ��(L1)=(��(J))2, as before. 
Proof. The proof is very technical and can be found in Serre [57] and Ciarlet [14]. As in the proof of the previous proposition, we begin by showing that the eigenvalues of the matrix L�� are the zeros of the polynomial 
�� + �� . 1 D 
qL�� (��) = det D . ��E . F = det . EpL�� (��),
�ئ� 
where pL�� (��) is the characteristic polynomial of L��. Then using the preliminary fact from Proposition 9.8, it is easy to show that 
��2 + �� . 1 
qL�� (��2)= ��n qJ ,
�˦� 
for all �� �� C, with �� = 0. This time we cannot extend the above equation to �� = 0. This leads us to consider the equation 
��2 + �� . 1 
= ��, 
�˦� 
which is equivalent to 
��2 . ���ئ� + �� . 1=0, 
for all �� = 0. Since �� = 0, the above equivalence does not hold for �� = 1, but this is not a problem since the case �� = 1 has already been considered in the previous proposition. Then we can show the following: 
1. For any �� = 0, if �� is an eigenvalue of L��, then 
�� + �� . 1 �� + �� . 1 
, . 
��1/2�ئ�1/2�� 
are eigenvalues of J. 
2. For every 	�� = 0, if �� and .�� are eigenvalues of J, then ��+(��, ��) and ��.(��, ��) are eigenvalues of L��, where ��+(��, ��) and ��.(��, ��) are the squares of the roots of the equation 
��2 . ���ئ� + �� . 1=0. 
It follows that ��(L��) = max {max(|��+(��, ��)|, |��.(��, ��)|)}, 
�� | pJ (��)=0
and since we are assuming that J has real roots, we are led to study the function 
M(��, ��) = max{|��+(��, ��)|, |��.(��, ��)|}, 
where �� �� R and �� �� (0, 2). Actually, because M(.��, ��)= M(��, ��), it is only necessary to consider the case where �� �� 0. Note that for �� = 0, the roots of the equation 
��2 . ���ئ� + �� . 1=0. 
are �� ���� �� ��2��2 . 4�� +4 
. 
2 In turn, this leads to consider the roots of the equation 
��2��2 . 4�� +4=0, 
which are 
�� 
2(1 �� 1 . ��2) 
��2 , 
for �� = 0. Since we have 
�� 	�̡� 
2(1+ 1 . ��2) 2(1+ 1 . ��2)(1 . 1 . ��2)2 
= 	�� = �� 
��2 
��2(1 . 1 . ��2)1 . 1 . ��2 
9.5. CONVERGENCE METHODS FOR TRIDIAGONAL MATRICES 
and 
�� �̡� 
2(1 . 1 . ��2) 2(1+ 1 . ��2)(1 . 1 . ��2)2 
= �� = �� ,
��2 
��2(1+ 1 . ��2) 1+ 1 . ��2 
these roots are  
2  2  
��0(��) = 1 + �� 1 . ��2 ,  ��1(��) = 1 . �� 1 . ��2 .  

Observe that the expression for ��0(��) is exactly the expression in the statement of our proposition! The rest of the proof consists in analyzing the variations of the function M(��, ��) by considering various cases for ��. In the end, we .nd that the minimum of ��(L��) is obtained for ��0(��(J)). The details are tedious and we omit them. The reader will .nd complete proofs in Serre [57] and Ciarlet [14]. 
Combining the results of Theorem 9.6 and Proposition 9.9, we obtain the following result which gives precise information about the spectral radii of the matrices J, L1, and L��. 
Proposition 9.10. Let A be a tridiagonal matrix (possibly by blocks) which is Hermitian positive de.nite. Then the methods of Jacobi, Gauss�CSeidel, and relaxation, all converge for �� �� (0, 2). There is a unique optimal relaxation parameter 
��0 = 1+һ 1 . 2(��(J))2 , 
such that 
��(L��0 ) = inf ��(L��)= ��0 . 1. 
0<��<2 
Furthermore, if ��(J) > 0, then 
��(L��0 ) <��(L1)=(��(J))2 <��(J), and if ��(J)=0, then ��0 =1 and ��(L1)= ��(J)=0. Proof. In order to apply Proposition 9.9, we have to check that J = D.1(E + F ) has real 
eigenvalues. However, if �� is any eigenvalue of J and if u is any corresponding eigenvector, then D.1(E + F )u = ��u implies that (E + F )u = ��Du, and since A = D . E . F , the above shows that (D . A)u = ��Du, that is, Au = (1 . ��)Du. Consequently, 
u . Au = (1 . ��)u . Du, and since A and D are Hermitian positive de.nite, we have u .Au > 0 and u .Du > 0 if u = 0, which proves that �� �� R. The rest follows from Theorem 9.6 and Proposition 9.9. 
Remark: It is preferable to overestimate rather than underestimate the relaxation param-eter when the optimum relaxation parameter is not known exactly. 
9.6 Summary 
The main concepts and results of this chapter are listed below: 
. 	
Iterative methods. Splitting A as A = M . N. 

. 	
Convergence of a sequence of vectors or matrices. 

. 	
A criterion for the convergence of the sequence (Bk) of powers of a matrix B to zero in terms of the spectral radius ��(B). 

. 	
A characterization of the spectral radius ��(B) as the limit of the sequence (1Bk11/k). 

. 	
A criterion of the convergence of iterative methods. 

. 	
Asymptotic behavior of iterative methods. 

. 	
Splitting A as A = D.E.F , and the methods of Jacobi, Gauss�CSeidel, and relaxation (and SOR). 

. 	
The Jacobi matrix, J = D.1(E + F ). 

. 	
The Gauss�CSeidel matrix, L1 =(D . E).1F . 

. 	
The matrix of relaxation, L�� =(D . ��E).1((1 . ��)D + ��F ). 

. 	
Convergence of iterative methods: a general result when A = M . N is Hermitian positive de.nite. 

. 	
A su.cient condition for the convergence of the methods of Jacobi, Gauss�CSeidel, and relaxation. The Ostrowski-Reich theorem: A is Hermitian positive de.nite and �� �� (0, 2). 

. 	
A necessary condition for the convergence of the methods of Jacobi , Gauss�CSeidel, and relaxation: �� �� (0, 2). 

. 	
The case of tridiagonal matrices (possibly by blocks). Simultaneous convergence or divergence of Jacobi��s method and Gauss�CSeidel��s method, and comparison of the spectral radii of ��(J) and ��(L1): ��(L1)=(��(J))2 . 

. 	
The case of tridiagonal Hermitian positive de.nite matrices (possibly by blocks). The methods of Jacobi, Gauss�CSeidel, and relaxation, all converge. 

. 	
In the above case, there is a unique optimal relaxation parameter for which ��(L��0 ) < ��(L1)=(��(J))2 <��(J) (if ��(J) = 0). 


9.7. PROBLEMS 341 
9.7 Problems 
Problem 9.1. Consider the matrix 
.
. 
. 

12 .2 
11 1

.

A = 

. 

22 1 Prove that ��(J)=0 and ��(L1)=2, so ��(J) < 1 <��(L1), where J is Jacobi��s matrix and L1 is the matrix of Gauss�CSeidel. Problem 9.2. Consider the matrix 
.
. 
A =

. 

2 .11 
2 22

.

. 

.1 .12 
�� 
Prove that ��(J)= 5/2 and ��(L1)=1/2, so ��(L1) <��(J), where where J is Jacobi��s matrix and L1 is the matrix of Gauss�CSeidel. Problem 9.3. Consider the following linear system: 
.
....
. 
2 .10 0 x1 19 
... 

.12 .10 

0 

.12 .1 

... 
... 

x2 
x3 
... 

= 

... 

19 

.3 

...

. 

00 .12 x4 .12 
(1) Solve the above system by Gaussian elimination. 
kkkk
(2) Compute the sequences of vectors uk =(u1,u2,u3,u4) for k =1,..., 10, using the methods of Jacobi, Gauss�CSeidel, and relaxation for the following values of ��: �� = 1.1, 1.2,..., 1.9. In all cases, the initial vector is u0 = (0, 0, 0, 0). 
c
Problem 9.4. Recall that a complex or real n �� n matrix A is strictly row diagonally 
n
dominant if |aii| > 
j=1,j��=i |aij| for i =1,...,n. 
(1) Prove that if A is strictly row diagonally dominant, then Jacobi��s method converges. 
(2) Prove that if A is strictly row diagonally dominant, then Gauss�CSeidel��s method converges. 
Problem 9.5. Prove that the converse of Proposition 9.5 holds. That is, if A is a Hermitian positive de.nite matrix writen as A = M . N with M invertible, if the Hermitan matrix M. + N is positive de.nite, and if ��(M.1N) < 1, then A is positive de.nite. 
Problem 9.6. Consider the following tridiagonal n �� n matrix: 

. .  
2  .1  0  
.1  2  .1 

...... 

...... 

A = 

1 

(n + 1)2 
...

.

.

.

.

. 

. 

. 

.12 .1 
0 .12 

(1) Prove that the eigenvalues of the Jacobi matrix J are given by k�� 
��k = cos ,k =1, . . . , n. 
n +1 
Hint. 

First show that the Jacobi matrix is 01 0 10 1 
.
. 

...... 

...... 

J = 

1 

2 

...

.

.

.

.

. 

. 

. 

1 01 
0 10 

Then the eigenvalues and the eigenvectors of J are solutions of the system of equations y0 =0 yk+1 + yk.1 =2��yk,k =1,...,n yn+1 =0. It is well known that the general solution to the above recurrence is given by yk = ��z1 k + ��z2 k ,k =0,...,n +1, (with ��, �� = 0) where z1 and z2 are the zeros of the equation z 2 . 2��z +1=0. It follows that z2 = z1 .1 and z1 + z2 =2��. The boundary condition y0 = 0 yields �� + �� = 0, so yk = ��(z1 k . z1 .k), and the boundary condition yn+1 = 0 yields 2(n+1)
z=1.
1 Deduce that we may assume that the n possible values (z1)k for z1 are given by 
k��i 
(z1)k = en+1 ,k =1, . . . , n, 
and .nd 2��k =(z1)k +(z1).k 1 . 
(k)(k)
Show that an eigenvector (y1 ,...,yn ) associated wih the eigenvalue ��k is given by 
(k) kj�� 
y= sin ,j =1, . . . , n. 
j 
n +1 
(2) Find the spectral radius ��(J), ��(L1), and ��(L��0 ), as functions of h =1/(n + 1). 


