Chapter 18 
Graphs and Graph Laplacians; Basic Facts 
In this chapter and the next we present some applications of linear algebra to graph theory. Graphs (undirected and directed) can be de.ned in terms of various matrices (incidence and adjacency matrices), and various connectivity properties of graphs are captured by properties of these matrices. Another very important matrix is associated with a (undirected) graph: the graph Laplacian. The graph Laplacian is symmetric positive de.nite, and its eigenvalues capture some of the properties of the underlying graph. This is a key fact that is exploited in graph clustering methods, the most powerful being the method of normalized cuts due to Shi and Malik [58]. This chapter and the next constitute an introduction to algebraic and spectral graph theory. We do not discuss normalized cuts, but we discuss graph drawings. Thorough presentations of algebraic graph theory can be found in Godsil and Royle [28] and Chung [13]. 
We begin with a review of basic notions of graph theory. Even though the graph Laplacian is fundamentally associated with an undirected graph, we review the de.nition of both directed and undirected graphs. For both directed and undirected graphs, we de.ne the degree matrix D, the incidence matrix B, and the adjacency matrix A. Then we de.ne a weighted graph. This is a pair (V, W ), where V is a .nite set of nodes and W is a m ℅ m symmetric matrix with nonnegative entries and zero diagonal entries (where m = |V |). 
For every node vi ﹋ V , the degree d(vi) (or di) of vi is the sum of the weights of the edges adjacent to vi: m
m 
di = d(vi)= wij. j=1 
The degree matrix is the diagonal matrix 
D = diag(d1,...,dm). 
The notion of degree is illustrated in Figure 18.1. Then we introduce the (unnormalized) 
graph Laplacian L of a directed graph G in an ※old-fashion§ way, by showing that for any 

607 


orientation of a graph G, BBT = D . A = L 
is an invariant. We also de.ne the (unnormalized) graph Laplacian L of a weighted graph G =(V, W ) as L = D . W . We show that the notion of incidence matrix can be generalized to weighted graphs in a simple way. For any graph G考 obtained by orienting the underlying graph of a weighted graph G =(V, W ), there is an incidence matrix B考 such that 
B考(B考)T = D . W = L. 
We also prove that 
m
1 m
xTLx = wij(xi . xj)2 for all x ﹋ Rm . 
2 
i,j=1 
Consequently, xTLx does not depend on the diagonal entries in W , and if wij ≡ 0 for all i, j ﹋{1,...,m}, then L is positive semide.nite. Then if W consists of nonnegative entries, the eigenvalues 0 = 竹1 ≒ 竹2 ≒ ... ≒ 竹m of L are real and nonnegative, and there is an orthonormal basis of eigenvectors of L. We show that the number of connected components of the graph G =(V, W ) is equal to the dimension of the kernel of L, which is also equal to the dimension of the kernel of the transpose (B考)T of any incidence matrix B考 obtained by orienting the underlying graph of G. 
We also de.ne the normalized graph Laplacians Lsym and Lrw, given by 
Lsym = D.1/2LD.1/2 = I . D.1/2WD.1/2 Lrw = D.1L = I . D.1W, 
and prove some simple properties relating the eigenvalues and the eigenvectors of L, Lsym and Lrw. These normalized graph Laplacians show up when dealing with normalized cuts. 
Next, we turn to graph drawings (Chapter 19). Graph drawing is a very attractive appli-cation of so-called spectral techniques, which is a fancy way of saying that that eigenvalues and eigenvectors of the graph Laplacian are used. Furthermore, it turns out that graph clustering using normalized cuts can be cast as a certain type of graph drawing. 
Given an undirected graph G =(V, E), with |V | = m, we would like to draw G in Rn for n (much) smaller than m. The idea is to assign a point 老(vi) in Rn to the vertex vi ﹋ V , for every vi ﹋ V , and to draw a line segment between the points 老(vi) and 老(vj). Thus, a graph drawing is a function 老: V ↙ Rn . 
We de.ne the matrix of a graph drawing 老 (in Rn) as a m ℅ n matrix R whose ith row consists of the row vector 老(vi) corresponding to the point representing vi in Rn . Typically, we want n<m; in fact n should be much smaller than m. 
Since there are in.nitely many graph drawings, it is desirable to have some criterion to decide which graph is better than another. Inspired by a physical model in which the edges are springs, it is natural to consider a representation to be better if it requires the springs to be less extended. We can formalize this by de.ning the energy of a drawing R by 
m 
E(R)= ��老(vi) . 老(vj)��2 , 
{vi,vj }﹋E
where 老(vi) is the ith row of R and��老(vi) . 老(vj)��2 is the square of the Euclidean length of the line segment joining 老(vi) and 老(vj). 
Then ※good drawings§ are drawings that minimize the energy function E. Of course, the trivial representation corresponding to the zero matrix is optimum, so we need to impose extra constraints to rule out the trivial solution. 
We can consider the more general situation where the springs are not necessarily identical. This can be modeled by a symmetric weight (or sti.ness) matrix W =(wij), with wij ≡ 0. In this case, our energy function becomes 
m 
E(R)= wij��老(vi) . 老(vj)��2 . 
{vi,vj }﹋E 
Following Godsil and Royle [28], we prove that 
E(R) = tr(RTLR), 
where L = D . W, 
is the familiar unnormalized Laplacian matrix associated with W , and where D is the degree matrix associated with W . 
It can be shown that there is no loss in generality in assuming that the columns of R are pairwise orthogonal and that they have unit length. Such a matrix satis.es the equation RTR = I and the corresponding drawing is called an orthogonal drawing. This condition also rules out trivial drawings. 
Then we prove the main theorem about graph drawings (Theorem 19.2), which essentially says that the matrix R of the desired graph drawing is constituted by the n eigenvectors of L associated with the smallest nonzero n eigenvalues of L. We give a number examples of graph drawings, many of which are borrowed or adapted from Spielman [60]. 
18.1 	Directed Graphs, Undirected Graphs, Incidence Matrices, Adjacency Matrices, Weighted Graphs 
De.nition 18.1. A directed graph is a pair G =(V, E), where V = {v1,...,vm} is a set of nodes or vertices, and E . V ℅ V is a set of ordered pairs of distinct nodes (that is, pairs (u, v) ﹋ V ℅ V with u Given any edge e =(u, v), we let s(e)= u be the 
= v), called edges. source of e and t(e)= v be the target of e. 
Remark: 	Since an edge is a pair (u, v) with u v, self-loops are not allowed. Also, there 
= is at most one edge from a node u to a node v. Such graphs are sometimes called simple graphs. 
An example of a directed graph is shown in Figure 18.2. 

De.nition 18.2. For every node v ﹋ V , the degree d(v) of v is the number of edges leaving or entering v: 
d(v)= |{u ﹋ V | (v, u) ﹋ E or (u, v) ﹋ E}|. 
We abbreviate d(vi) as di. The degree matrix, D(G), is the diagonal matrix 
D(G) = diag(d1,...,dm). 
18.1. DIRECTED GRAPHS, UNDIRECTED GRAPHS, WEIGHTED GRAPHS 611 For example, for graph G1, we have 
.
. 
D(G1)= 
..... 

20000 
04000 

00300 

00030 
00002 

..... 

. 

Unless confusion arises, we write D instead of D(G). 
De.nition 18.3. Given a directed graph G =(V, E), for any two nodes u, v ﹋ V ,a path from u to v is a sequence of nodes (v0,v1,...,vk) such that v0 = u, vk = v, and (vi,vi+1) is an edge in E for all i with 0 ≒ i ≒ k . 1. The integer k is the length of the path. A path is closed if u = v. The graph G is strongly connected if for any two distinct nodes u, v ﹋ V , there is a path from u to v and there is a path from v to u. 
Remark: The terminology walk is often used instead of path, the word path being reserved to the case where the nodes vi are all distinct, except that v0 = vk when the path is closed. 
The binary relation on V ℅ V de.ned so that u and v are related i. there is a path from u to v and there is a path from v to u is an equivalence relation whose equivalence classes are called the strongly connected components of G. 
De.nition 18.4. Given a directed graph G =(V, E), with V = {v1,...,vm}, if E = {e1,...,en}, then the incidence matrix B(G) of G is the m ℅ n matrix whose entries bij are given by 
. .. .. 

+1 if s(ej)= vi 
.1 if t(ej)= vi
bij 
= 

0 otherwise. 

Here is the incidence matrix of the graph G1: 
.
. 
B = 

..... 

1100000 
.10 .1 .11 0 0 
0 .11 0 0 0 1 

00010 .1 .1 
0000 .11 0 

..... 

. 

Observe that every column of an incidence matrix contains exactly two nonzero entries, +1 and .1. Again, unless confusion arises, we write B instead of B(G). 
When a directed graph has m nodes v1,...,vm and n edges e1,...,en, a vector x ﹋ Rm can be viewed as a function x: V ↙ R assigning the value xi to the node vi. Under this interpretation, Rm is viewed as RV . Similarly, a vector y ﹋ Rn can be viewed as a function in RE . This point of view is often useful. For example, the incidence matrix B can be interpreted as a linear map from RE to RV , the boundary map, and BT can be interpreted as a linear map from RV to RE, the coboundary map. 

Remark: Some authors adopt the opposite convention of sign in de.ning the incidence matrix, which means that their incidence matrix is .B. 
Undirected graphs are obtained from directed graphs by forgetting the orientation of the edges. 
De.nition 18.5. A graph (or undirected graph) is a pair G =(V, E), where V = {v1,...,vm}is a set of nodes or vertices, and E is a set of two-element subsets of V (that is, subsets {u, v}, with u, v ﹋ V and u = v), called edges. 
Remark: Since an edge is a set {u, v}, we have u = v, so self-loops are not allowed. Also, for every set of nodes {u, v}, there is at most one edge between u and v. As in the case of directed graphs, such graphs are sometimes called simple graphs. 
An example of a graph is shown in Figure 18.3. 
De.nition 18.6. For every node v ﹋ V , the degree d(v) of v is the number of edges incident to v: 
d(v)= |{u ﹋ V |{u, v}﹋ E}|. 
The degree matrix D(G) (or simply, D) is de.ned as in De.nition 18.2. 
De.nition 18.7. Given a (undirected) graph G =(V, E), for any two nodes u, v ﹋ V ,a path from u to v is a sequence of nodes (v0,v1,...,vk) such that v0 = u, vk = v, and {vi,vi+1} is an edge in E for all i with 0 ≒ i ≒ k . 1. The integer k is the length of the path. A path is closed if u = v. The graph G is connected if for any two distinct nodes u, v ﹋ V , there is a path from u to v. 
Remark: The terminology walk or chain is often used instead of path, the word path being reserved to the case where the nodes vi are all distinct, except that v0 = vk when the path is closed. 

18.1. DIRECTED GRAPHS, UNDIRECTED GRAPHS, WEIGHTED GRAPHS 613 
The binary relation on V ℅V de.ned so that u and v are related i. there is a path from u to v is an equivalence relation whose equivalence classes are called the connected components of G. 
The notion of incidence matrix for an undirected graph is not as useful as in the case of directed graphs 
De.nition 18.8. Given a graph G =(V, E), with V = {v1,...,vm}, if E = {e1,...,en}, then the incidence matrix B(G) of G is the m ℅ n matrix whose entries bij are given by 
 

+1 if ej = {vi,vk} for some k 
bij =
0 otherwise. 
Unlike the case of directed graphs, the entries in the incidence matrix of a graph (undi-rected) are nonnegative. We usually write B instead of B(G). 
De.nition 18.9. If G =(V, E) is a directed or an undirected graph, given a node u ﹋ V , any node v ﹋ V such that there is an edge (u, v) in the directed case or {u, v} in the undirected case is called adjacent to u, and we often use the notation 
u ‵ v. 
Observe that the binary relation ‵ is symmetric when G is an undirected graph, but in general it is not symmetric when G is a directed graph. 
The notion of adjacency matrix is basically the same for directed or undirected graphs. 
De.nition 18.10. Given a directed or undirected graph G =(V, E), with V = {v1,...,vm}, the adjacency matrix A(G) of G is the symmetric m ℅ m matrix (aij) such that 
(1) If G is directed, then 
 

1 if there is some edge (vi,vj) ﹋ E or some edge (vj,vi) ﹋ E 
aij =
0 otherwise. 
(2) Else if G is undirected, then 
 

1 if there is some edge {vi,vj}﹋ E 
aij =
0 otherwise. 
As usual, unless confusion arises, we write A instead of A(G). Here is the adjacency matrix of both graphs G1 and G2: 
.
. 
A = 

..... 

0  1  1  0  0  
1  0  1  1  1  
1  1  0  1  0  

01101 
01010 

..... 

. 

If G =(V, E) is an undirected graph, the adjacency matrix A of G can be viewed as a linear map from RV to RV , such that for all x ﹋ Rm, we have 
m 
(Ax)i = xj; j‵i 
that is, the value of Ax at vi is the sum of the values of x at the nodes vj adjacent to vi. The adjacency matrix can be viewed as a di.usion operator. This observation yields a geometric interpretation of what it means for a vector x ﹋ Rm to be an eigenvector of A associated with some eigenvalue 竹; we must have 
m 
竹xi = xj,i =1, . . . , m, j‵i 
which means that the the sum of the values of x assigned to the nodes vj adjacent to vi is equal to 竹 times the value of x at vi. 
De.nition 18.11. Given any undirected graph G =(V, E), an orientation of G is a function 考 : E ↙ V ℅ V assigning a source and a target to every edge in E, which means that for every edge {u, v}﹋ E, either 考({u, v})=(u, v) or 考({u, v})=(v, u). The oriented graph G考 obtained from G by applying the orientation 考 is the directed graph G考 =(V, E考), with E考 = 考(E). 
The following result shows how the number of connected components of an undirected graph is related to the rank of the incidence matrix of any oriented graph obtained from G. 
Proposition 18.1. Let G =(V, E) be any undirected graph with m vertices, n edges, and c connected components. For any orientation 考 of G, if B is the incidence matrix of the oriented graph G考 , then c = dim(Ker (BT)), and B has rank m . c. Furthermore, the nullspace of BT has a basis consisting of indicator vectors of the connected components of G; that is, vectors (z1,...,zm) such that zj =1 i. vj is in the ith component Ki of G, and zj =0 otherwise. 
Proof. (After Godsil and Royle [28], Section 8.3). The fact that rank(B)= m . c will be proved last. 
Let us prove that the kernel of BT has dimension c. A vector z ﹋ Rm belongs to the kernel of BT i. BTz = 0 i. zTB = 0. In view of the de.nition of B, for every edge {vi,vj}of G, the column of B corresponding to the oriented edge 考({vi,vj}) has zero entries except for a +1 and a .1 in position i and position j or vice-versa, so we have 
zi = zj. 
An easy induction on the length of the path shows that if there is a path from vi to vj in G (unoriented), then zi = zj. Therefore, z has a constant value on any connected component of 
G. It follows that every vector z ﹋ Ker (BT) can be written uniquely as a linear combination 
z = 竹1z 1 + ﹞﹞﹞ + 竹cz c , 

18.1. DIRECTED GRAPHS, UNDIRECTED GRAPHS, WEIGHTED GRAPHS 615 

where the vector zi corresponds to the ith connected component Ki of G and is de.ned such that 
i 1 i. vj ﹋ Ki 
z = 
j 
0 otherwise. 
This shows that dim(Ker (BT)) = c, and that Ker (BT) has a basis consisting of indicator vectors. 
Since BT is a n ℅ m matrix, we have 
m = dim(Ker (BT)) + rank(BT), 
and since we just proved that dim(Ker (BT)) = c, we obtain rank(BT)= m . c. Since B and BT have the same rank, rank(B)= m . c, as claimed. 
De.nition 18.12. Following common practice, we denote by 1 the (column) vector (of dimension m) whose components are all equal to 1. 
Since every column of B contains a single +1 and a single .1, the rows of BT sum to zero, which can be expressed as BT1 =0. 
According to Proposition 18.1, the graph G is connected i. B has rank m.1 i. the nullspace of BT is the one-dimensional space spanned by 1. 
In many applications, the notion of graph needs to be generalized to capture the intuitive idea that two nodes u and v are linked with a degree of certainty (or strength). Thus, we assign a nonnegative weight wij to an edge {vi,vj}; the smaller wij is, the weaker is the link (or similarity) between vi and vj, and the greater wij is, the stronger is the link (or similarity) between vi and vj. 
De.nition 18.13. A weighted graph is a pair G =(V, W ), where V = {v1,...,vm} is a set of nodes or vertices, and W is a symmetric matrix called the weight matrix, such that wij ≡ 0 for all i, j ﹋{1,...,m}, and wii = 0 for i =1,...,m. We say that a set {vi,vj} is an edge i. wij > 0. The corresponding (undirected) graph (V, E) with E = {{vi,vj}| wij > 0}, is called the underlying graph of G. 
Remark: Since wii = 0, these graphs have no self-loops. We can think of the matrix W as a generalized adjacency matrix. The case where wij ﹋{0, 1} is equivalent to the notion of a graph as in De.nition 18.5. 
We can think of the weight wij of an edge {vi,vj} as a degree of similarity (or a.nity) in an image, or a cost in a network. An example of a weighted graph is shown in Figure 18.4. The thickness of an edge corresponds to the magnitude of its weight. 

De.nition 18.14. Given a weighted graph G =(V, W ), for every node vi ﹋ V , the degree d(vi) of vi is the sum of the weights of the edges adjacent to vi: 
m
m 
d(vi)= wij. j=1 
Note that in the above sum, only nodes vj such that there is an edge {vi,vj} have a nonzero contribution. Such nodes are said to be adjacent to vi, and we write vi ‵ vj. The degree matrix D(G) (or simply, D) is de.ned as before, namely by D(G) = diag(d(v1),...,d(vm)). 
The weight matrix W can be viewed as a linear map from RV to itself. For all x ﹋ Rm , we have 
m 
(Wx)i = wijxj; j‵i 
that is, the value of Wx at vi is the weighted sum of the values of x at the nodes vj adjacent to vi. 
Observe that W 1 is the (column) vector (d(v1),...,d(vm)) consisting of the degrees of the nodes of the graph. 
We now de.ne the most important concept of this chapter: the Laplacian matrix of a graph. Actually, as we will see, it comes in several .avors. 


18.2 Laplacian Matrices of Graphs 
Let us begin with directed graphs, although as we will see, graph Laplacians are funda-mentally associated with undirected graph. The key proposition below shows how given an 
18.2. 	LAPLACIAN MATRICES OF GRAPHS 
undirected graph G, for any orientation 考 of G, B考(B考)T relates to the adjacency matrix A (where B考 is the incidence matrix of the directed graph G考). We reproduce the proof in Gallier [24] (see also Godsil and Royle [28]). 
Proposition 18.2. Given any undirected graph G, for any orientation 考 of G, if B考is the incidence matrix of the directed graph G考 , A is the adjacency matrix of G考, and D is the degree matrix such that Dii = d(vi), then 
  
B考(B考)T = D . A. 
Consequently, L = B考(B考)T is independent of the orientation 考 of G, and D.A is symmetric 
  
and positive semide.nite; that is, the eigenvalues of D . A are real and nonnegative. 
Proof. 	The entry B考(B考)T i , and the jth row b考j .
ij is the inner product of the ith row b考 of B考 If i = j, then as +1 if s(ek)= vi b考 
= if t(ek)= vi
ik 	.1 0 otherwise 
. .. 	.. 考考考考weseethat bbd().If ij,then bb=0i.thereissomeedge with()= ﹞﹞= v= esevikki
ij andt(ek)=vjiorivice-versa (which are mutually exclusive cases, since G考 arises by orienting an undirected graph), in which case, b考i ﹞ bj考 = .1. Therefore, 
B考(B考)T = D . A, 
as claimed. For every x ﹋ Rm, we have 
x TLx = x TB考(B考)T x = ((B考)T x)T(B考)T x =(B考)T x22 ≡ 0, 
since the Euclidean norm���� is positive (de.nite). Therefore, L = B考(B考)T is positive 
2 
semide.nite. It is well-known that a real symmetric matrix is positive semide.nite i. its eigenvalues are nonnegative. 
De.nition 18.15. The matrix L = B考(B考)T = D . A is called the (unnormalized) graph Laplacian of the graph G考 . The (unnormalized) graph Laplacian of an undirected graph G =(V, E) is de.ned by 
L = D . A. 
For example, the graph Laplacian of graph G1 is 
2  .1  .1  0  0  
L = ..... .1 .1 0  4 .1 .1  .1 3 .1  .1 .1 3  .1 0 .1 ..... .  
0  .1  0  .1  2  

.
. 
Observe that each row of L sums to zero (because (B考)T1 = 0). Consequently, the vector 1 is in the nullspace of L. 
Remarks: 
1. 
With the unoriented version of the incidence matrix (see De.nition 18.8), it can be 

shown that BBT = D + A. 

2. 
As pointed out by Evangelos Chatzipantazis, Proposition 18.2 in which the incidence matrix B考 is replaced by the incidence matrix B of any arbitrary directed graph G does not hold. The problem is that such graphs may have both edges (vi,vj) and (vj,vi) between two distinct nodes vi and vj, and as a consequence, the inner product bi ﹞ bj = .2 instead of .1. A simple counterexample is given by the directed graph with three vertices and four edges whose incidence matrix is given by 


.. 
1 .10 .1 
..
B = .11 .10 . 0011 
We have  
.  .  .  .  .  .  
3  .2  .1  3  0  0  0  1  1  
BBT = ..2  3  .1. = .0  3  0. . .1  0  1. = D . A.  
.1  .1  2  0  0  2  1  1  0  

The natural generalization of the notion of graph Laplacian to weighted graphs is this: De.nition 18.16. Given any weighted graph G =(V, W ) with V = {v1,...,vm}, the (unnormalized) graph Laplacian L(G) of G is de.ned by L(G)= D(G) . W, where D(G) = diag(d1,...,dm) is the degree matrix of G (a diagonal matrix), with 
m
m 
di = wij. j=1 
As usual, unless confusion arises, we write D instead of D(G) and L instead of L(G). 
The graph Laplacian can be interpreted as a linear map from RV to itself. For all x ﹋ RV , we have 
m 
(Lx)i = wij(xi . xj). j‵i 
It is clear from the equation L = D . W that each row of L sums to 0, so the vector 1 is the nullspace of L, but it is less obvious that L is positive semide.nite. One way to prove it is to generalize slightly the notion of incidence matrix. 

18.2. LAPLACIAN MATRICES OF GRAPHS 
De.nition 18.17. Given a weighted graph G =(V, W ), with V = {v1,...,vm}, if {e1,..., en} are the edges of the underlying graph of G (recall that {vi,vj} is an edge of this graph i. wij > 0), for any oriented graph G考 obtained by giving an orientation to the underlying graph of G, the incidence matrix B考 of G考 is the m ℅ n matrix whose entries bij are given 
by 

bij 
= 

. .. .. 

﹟ 
+ wij if s(ej)= vi 
﹟
.

wij 
if t(ej)= vi 
0 otherwise. 

For example, given the weight matrix 

.
. 
W = 

... 

0363 
3003 

6003 
3330 

...

, 

the incidence matrix B corresponding to the orientation of the underlying graph of W where an edge (i, j) is oriented positively i. i<j is 
.
. 
B = 

... 

1.7321 2.4495 1.73210 0 
.1.73210 0 1.7321 0 

0 .2.44950 0 1.7321 
00 .1.7321 .1.7321 .1.7321 

...

. 

The reader should verify that BBT = D . W . This is true in general, see Proposition 18.3. 
It is easy to see that Proposition 18.1 applies to the underlying graph of G. For any oriented graph G考 obtained from the underlying graph of G, the rank of the incidence matrix B考 is equal to m.c, where c is the number of connected components of the underlying graph of G, and we have (B考)T1 = 0. We also have the following version of Proposition 18.2 whose proof is immediately adapted. 
Proposition 18.3. Given any weighted graph G =(V, W ) with V = {v1,...,vm}, if B考 is the incidence matrix of any oriented graph G考 obtained from the underlying graph of G and D is the degree matrix of G, then 
B考(B考)T = D . W = L. 
Consequently, B考(B考)T is independent of the orientation of the underlying graph of G and L = D . W is symmetric and positive semide.nite; that is, the eigenvalues of L = D . W are real and nonnegative. 
Another way to prove that L is positive semide.nite is to evaluate the quadratic form xTLx. 
= 
j=1 
(
Proposition 18.4. For any m ℅ m symmetric matrix W 

where D is the degree matrix associated with W (that is, di 
m
m 
D . W

= 

(wij), if we let L = 
m 
wij), then we have 
1
TLx = 
wij(xi . xj)2 for all x ﹋ Rm 
.

x 

2 
i,j=1 
Consequently, xTLx does not depend on the diagonal entries in W , and if wij ≡ 0 for all i, j ﹋{1,...,m}, then L is positive semide.nite. 
Proof. We have 
x TLx = x TDx . x TWx 
mm= dixi . wijxixj i=1 i,j=1 m
m 
2 m 
m 

 

 

m
m
22
dixi . 2 wijxixj + dix 
mm
1 

= 

i
2
i=1 m
m 
i,j=1 i=1 
1 

wij(xi . xj)2 
= 

. 

2 
i,j=1 
Obviously, the quantity on the right-hand side does not depend on the diagonal entries in W , and if wij ≡ 0 for all i, j, then this quantity is nonnegative. 
Proposition 18.4 immediately implies the following facts: For any weighted graph G = (V, W ), 
1. 
The eigenvalues 0 = 竹1 ≒ 竹2 ≒ ... ≒ 竹m of L are real and nonnegative, and there is an orthonormal basis of eigenvectors of L. 

2. 
The smallest eigenvalue 竹1 of L is equal to 0, and 1 is a corresponding eigenvector. 


It turns out that the dimension of the nullspace of L (the eigenspace of 0) is equal to the number of connected components of the underlying graph of G. 
Proposition 18.5. Let G =(V, W ) be a weighted graph. The number c of connected com-ponents K1,...,Kc of the underlying graph of G is equal to the dimension of the nullspace of L, which is equal to the multiplicity of the eigenvalue 0. Furthermore, the nullspace of L has a basis consisting of indicator vectors of the connected components of G, that is, vectors (f1,...,fm) such that fj =1 i. vj ﹋ Ki and fj =0 otherwise. 
Proof. Since L = BBT for the incidence matrix B associated with any oriented graph obtained from G, and since L and BT have the same nullspace, by Proposition 18.1, the dimension of the nullspace of L is equal to the number c of connected components of G and the indicator vectors of the connected components of G form a basis of Ker (L). 

18.3. NORMALIZED LAPLACIAN MATRICES OF GRAPHS 
Proposition 18.5 implies that if the underlying graph of G is connected, then the second eigenvalue 竹2 of L is strictly positive. 
Remarkably, the eigenvalue 竹2 contains a lot of information about the graph G (assuming that G =(V, E) is an undirected graph). This was .rst discovered by Fiedler in 1973, and for this reason, 竹2 is often referred to as the Fiedler number. For more on the properties of the Fiedler number, see Godsil and Royle [28] (Chapter 13) and Chung [13]. More generally, the spectrum (0,竹2,...,竹m) of L contains a lot of information about the combinatorial structure of the graph G. Leverage of this information is the object of spectral graph theory. 


18.3 Normalized Laplacian Matrices of Graphs 
It turns out that normalized variants of the graph Laplacian are needed, especially in appli-cations to graph clustering. These variants make sense only if G has no isolated vertices. 
De.nition 18.18. Given a weighted graph G =(V, W ), a vertex u ﹋ V is isolated if it is not incident to any other vertex. This means that every row of W contains some strictly positive entry. 
If G has no isolated vertices, then the degree matrix D contains positive entries, so it is invertible and D.1/2 makes sense; namely 
D.1/2 = diag(d.1/2,...,d.1/2),
1 m 
and similarly for any real exponent 汐. 
De.nition 18.19. Given any weighted directed graph G =(V, W ) with no isolated vertex and with V = {v1,...,vm}, the (normalized) graph Laplacians Lsym and Lrw of G are de.ned by 
Lsym = D.1/2LD.1/2 = I . D.1/2WD.1/2 Lrw = D.1L = I . D.1W. 
Observe that the Laplacian Lsym = D.1/2LD.1/2 is a symmetric matrix (because L and D.1/2 are symmetric) and that D1/2
Lrw = D.1/2Lsym. 
The reason for the notation Lrw is that this matrix is closely related to a random walk on the graph G. 
Example 18.1. As an example, the matrices Lsym and Lrw associated with the graph G1 
are

.
. 
..... 

1.0000 .0.3536 .0.40820 0 .0.3536 1.0000 .0.2887 .0.2887 .0.3536 .0.4082 .0.2887 1.0000 .0.3333 0 
0 .0.2887 .0.3333 1.0000 .0.4082 
0 .0.3536 0 .0.4082 1.0000 

..... 

Lsym = 
and

.
. 
Lrw = 
..... 

1.0000 .0.5000 .0.50000 0 .0.2500 1.0000 .0.2500 .0.2500 .0.2500 .0.3333 .0.3333 1.0000 .0.3333 0 
0 .0.3333 .0.3333 1.0000 .0.3333 
0 .0.5000 0 .0.5000 1.0000 

..... 

. 

Since the unnormalized Laplacian L can be written as L = BBT, where B is the incidence matrix of any oriented graph obtained from the underlying graph of G =(V, W ), if we let 
Bsym = D.1/2B, 
we get 
BT
Lsym = Bsymsym. In particular, for any singular decomposition Bsym = U曳V T of Bsym (with U an m ℅ m orthogonal matrix, 曳 a ※diagonal§ m℅n matrix of singular values, and V an n℅n orthogonal matrix), the eigenvalues of Lsym are the squares of the top m singular values of Bsym, and the vectors in U are orthonormal eigenvectors of Lsym with respect to these eigenvalues (the squares of the top m diagonal entries of 曳). Computing the SVD of Bsym generally yields more accurate results than diagonalizing Lsym, especially when Lsym has eigenvalues with 
m 
high multiplicity. 
There are simple relationships between the eigenvalues and the eigenvectors of Lsym, and Lrw. There is also a simple relationship with the generalized eigenvalue problem Lx = 竹Dx. 
Proposition 18.6. Let G =(V, W ) be a weighted graph without isolated vertices. The graph Laplacians, L, Lsym, and Lrw satisfy the following properties: 
(1) The matrix Lsym is symmetric and positive semide.nite. In fact, 
2
m
1 

xi 
-xj
﹟. 
di dj 
x 

TLsymx = 
for all x ﹋ Rm 
.

wij 
2 

i,j=1 
(2) 
The normalized graph Laplacians Lsym and Lrw have the same spectrum (0 = 糸1 ≒ 糸2 ≒ ... ≒ 糸m), and a vector u =0 is an eigenvector of Lrw for 竹 i. D1/2u is an eigenvector of Lsym for 竹. 

(3) 
The graph Laplacians L and Lsym are symmetric and positive semide.nite. 

(4) 
A vector u =0 is a solution of the generalized eigenvalue problem Lu = 竹Du i. D1/2u is an eigenvector of Lsym for the eigenvalue 竹 i. u is an eigenvector of Lrw for the eigenvalue 竹. 

(5) 
The graph Laplacians, L and Lrw have the same nullspace. For any vector u, we have u ﹋ Ker (L) i. D1/2u ﹋ Ker (Lsym). 

18.3. NORMALIZED LAPLACIAN MATRICES OF GRAPHS 

(6) 
The vector 1 is in the nullspace of Lrw, and D1/21 is in the nullspace of Lsym. 

(7) 
For every eigenvalue 糸i of the normalized graph Laplacian Lsym, we have 0 ≒ 糸i ≒ 2. Furthermore, 糸m =2 i. the underlying graph of G contains a nontrivial connected bipartite component. 

(8) 
If 	m ≡ 2 and if the underlying graph of G is not a complete graph,1 then 糸2 ≒ 1. Furthermore the underlying graph of G is a complete graph i. 糸2 = mm .1 . 

(9) 
If m ≡ 2 and if the underlying graph of G is connected, then 糸2 > 0. 

(10) 
If m ≡ 2 and if the underlying graph of G has no isolated vertices, then 糸m ≡ mm .1 . 


Proof. (1) We have Lsym = D.1/2LD.1/2, and D.1/2 is a symmetric invertible matrix (since it is an invertible diagonal matrix). It is a well-known fact of linear algebra that if B is an invertible matrix, then a matrix S is symmetric, positive semide.nite i. BSBT is symmetric, positive semide.nite. Since L is symmetric, positive semide.nite, so is Lsym = D.1/2LD.1/2 . The formula 
m
1 mxi xj 2 
x TLsymx = wij ﹟. -for all x ﹋ Rm 
2 di dj
i,j=1 
follows immediately from Proposition 18.4 by replacing x by D.1/2x, and also shows that Lsym is positive semide.nite. 
(2) Since 
Lrw = D.1/2LsymD1/2 , 
the matrices Lsym and Lrw are similar, which implies that they have the same spectrum. In fact, since D1/2 is invertible, 
Lrwu = D.1Lu = 竹u 
i. 
D.1/2Lu = 竹D1/2 u 
i. 
D1/2
D.1/2LD.1/2D1/2 u = Lsymu = 竹D1/2 u, 
which shows that a vector u = 0 is an eigenvector of Lrw for 竹 i. D1/2u is an eigenvector of Lsym for 竹. 
(3) 
We already know that L and Lsym are positive semide.nite. 

(4) 
Since D.1/2 is invertible, we have 


Lu = 竹Du 
i. 
D.1/2Lu = 竹D1/2 u 
1Recall that an undirected graph is complete if for any two distinct nodes u, v, there is an edge {u, v}. 
i. 
D1/2
D.1/2LD.1/2D1/2 u = Lsymu = 竹D1/2 u, 
which shows that a vector u = 0 is a solution of the generalized eigenvalue problem Lu = 竹Du i. D1/2u is an eigenvector of Lsym for the eigenvalue 竹. The second part of the statement follows from (2). 
(5) 
Since D.1 is invertible, we have Lu = 0 i. D.1Lu = Lrwu = 0. Similarly, since D.1/2 is invertible, we have Lu = 0 i. D.1/2LD.1/2D1/2u = 0 i. D1/2u ﹋ Ker (Lsym). 

(6) 
Since L1 = 0, we get Lrw1 = D.1L1 = 0. That D1/21 is in the nullspace of Lsym follows from (2). Properties (7)每(10) are proven in Chung [13] (Chapter 1). 


The eigenvalues the matrices Lsym and Lrw from Example 18.1 are 
0, 7257, 1.1667, 1.5, 1.6076. 
On the other hand, the eigenvalues of the unormalized Laplacian for G1 are 
0, 1.5858, 3, 4.4142, 5. 
Remark: Observe that although the matrices Lsym and Lrw have the same spectrum, the matrix Lrw is generally not symmetric, whereas Lsym is symmetric. 
A version of Proposition 18.5 also holds for the graph Laplacians Lsym and Lrw. This fol-lows easily from the fact that Proposition 18.1 applies to the underlying graph of a weighted graph. The proof is left as an exercise. 
Proposition 18.7. Let G =(V, W ) be a weighted graph. The number c of connected com-ponents K1,...,Kc of the underlying graph of G is equal to the dimension of the nullspace of both Lsym and Lrw, which is equal to the multiplicity of the eigenvalue 0. Furthermore, the nullspace of Lrw has a basis consisting of indicator vectors of the connected components of G, that is, vectors (f1,...,fm) such that fj =1 i. vj ﹋ Ki and fj =0 otherwise. For Lsym, a basis of the nullpace is obtained by multiplying the above basis of the nullspace of Lrw by D1/2 . 
A particularly interesting application of graph Laplacians is graph clustering. 
18.4 Graph Clustering Using Normalized Cuts 
In order to explain this problem we need some de.nitions. 
De.nition 18.20. Given any subset of nodes A . V , we de.ne the volume vol(A) of A as the sum of the weights of all edges adjacent to nodes in A: 
m m
m 
vol(A)= wij. vi﹋Aj=1 
18.4. GRAPH CLUSTERING USING NORMALIZED CUTS 
Given any two subsets A, B . V (not necessarily distinct), we de.ne links(A, B) by 

m 
links(A, B)= wij. vi﹋A,vj ﹋B 
The quantity links(A, A) = links(A, A) (where A = V . A denotes the complement of A in V ) measures how many links escape from A (and A). We de.ne the cut of A as 
cut(A) = links(A, A). 
The notion of volume is illustrated in Figure 18.5 and the notions of cut is illustrated in Figure 18.6. 


The above concepts play a crucial role in the theory of normalized cuts. This beautiful and deeply original method .rst published in Shi and Malik [58], has now come to be a ※textbook chapter§ of computer vision and machine learning. It was invented by Jianbo Shi and Jitendra Malik and was the main topic of Shi＊s dissertation. This method was extended to K ≡ 3 clusters by Stella Yu in her dissertation [75] and is also the subject of Yu and Shi [76]. 
Given a set of data, the goal of clustering is to partition the data into di.erent groups according to their similarities. When the data is given in terms of a similarity graph G, where the weight wij between two nodes vi and vj is a measure of similarity of vi and vj, the problem can be stated as follows: Find a partition (A1,...,AK ) of the set of nodes V into di.erent groups such that the edges between di.erent groups have very low weight (which indicates that the points in di.erent clusters are dissimilar), and the edges within a group have high weight (which indicates that points within the same cluster are similar). 
The above graph clustering problem can be formalized as an optimization problem, using the notion of cut mentioned earlier. If we want to partition V into K clusters, we can do so by .nding a partition (A1,...,AK ) that minimizes the quantity 
KK
mm
11 
cut(A1,...,AK ) = cut(Ai) = links(Ai,Ai). 
22 
i=1 i=1 
For K = 2, the mincut problem is a classical problem that can be solved e.ciently, but in practice, it does not yield satisfactory partitions. Indeed, in many cases, the mincut solution separates one vertex from the rest of the graph. What we need is to design our cost function in such a way that it keeps the subsets Ai ※reasonably large§ (reasonably balanced). 
An example of a weighted graph and a partition of its nodes into two clusters is shown in Figure 18.7. 

A way to get around this problem is to normalize the cuts by dividing by some measure of each subset Ai. A solution using the volume vol(Ai) of Ai (for K = 2) was proposed and 
18.5. SUMMARY 
investigated in a seminal paper of Shi and Malik [58]. Subsequently, Yu (in her dissertation [75]) and Yu and Shi [76] extended the method to K> 2 clusters. The idea is to minimize the cost function 
KK
m
m 
Ncut(A1,...,AK )= links(Ai,Ai) 
= 

cut(Ai,Ai) 
vol(Ai) vol(Ai)
i=1 i=1 
. 

The next step is to express our optimization problem in matrix form, and this can be done in terms of Rayleigh ratios involving the graph Laplacian in the numerators. This theory is very beautiful, but we do not have the space to present it here. The interested reader is referred to Gallier [23]. 
18.5 Summary 
The main concepts and results of this chapter are listed below: 
. 
Directed graphs, undirected graphs. 

. 
Incidence matrices, adjacency matrices. 

. 
Weighted graphs. 

. 
Degree matrix. 

. 
Graph Laplacian (unnormalized). 

. 
Normalized graph Laplacian. 

. 
Spectral graph theory. 

. 
Graph clustering using normalized cuts. 


18.6 Problems 
Problem 18.1. Find the unnormalized Laplacian of the graph representing a triangle and of the graph representing a square. 
Problem 18.2. Consider the complete graph Km on m ≡ 2 nodes. 
(1) Prove that the normalized Laplacian Lsym of K is 
.
. 
1 .1/(m . 1) ... .1/(m . 1) .1/(m . 1) .1/(m . 1) 1 ... .1/(m . 1) .1/(m . 1) ..
...
Lsym = 
...... 

...... 

.

. 

.

.

..

. 

. 

.

. 

. 

.1/(m . 1) .1/(m . 1) ... 1 .1/(m . 1) .1/(m . 1) .1/(m . 1) ... .1/(m . 1) 1 
(2) Prove that the characteristic polynomial of Lsym is
竹 . 11/(m . 1) ... 1/(m . 1) 1/(m . 1) 

           

           

1/(m . 1) 竹 . 1 ... 1/(m . 1) 1/(m . 1)

 
 
m.1 
= 竹竹 . 

m 

m . 1

..
...
. 

. 

.

.. 

.

. 

. 

.

. 

. 

1/(m . 1) 1/(m . 1) ... 竹 . 11/(m . 1) 1/(m . 1) 1/(m . 1) ... 1/(m . 1) 竹 . 1
Hint. First subtract the second column from the .rst, factor 竹 . m/(m . 1), and then add the .rst row to the second. Repeat this process. You will end up with the determinant
    

竹 . 1/(m . 1) 1 

1/(m . 1) 竹 . 1

    

. 

Problem 18.3. Consider the complete bipartite graph Km,n on m + n ≡ 3 nodes, with edges between each of the .rst m ≡ 1 nodes to each of the last n ≡ 1 nodes. Prove that the eigenvalues of the normalized Laplacian Lsym of Km,n are 0 with multiplicity m + n . 2 and 1 with multiplicity 2. 
Problem 18.4. Let G be a graph with a set of nodes V with m ≡ 2 elements, without isolated nodes, and let Lsym = D.1/2LD.1/2 be its normalized Laplacian (with L its unnor-malized Laplacian). 
(1) For any y ﹋ RV , consider the Rayleigh ratio 
yTLsym y
R = . 
yTy 
Prove that if x = D.1/2y, then 
m 

(x(u) . x(v))2 
xTLx u‵v
R ==

m

. 

(D1/2x)TD1/2xdvx(v)2 
v 
(2) Prove that the second eigenvalue 糸2 of Lsym is given by 
m 

(x(u) . x(v))2 
u‵v
糸2 = min
m

. 

1TDx=0,x珧 dvx(v)2
=0 
v 
(3) Prove that the largest eigenvalue 糸m of Lsym is given by 
m 

(x(u) . x(v))2 
u‵v
糸m = max 
x珧
=0 
m

. 

dvx(v)2 
v 
18.6. PROBLEMS 
Problem 18.5. Let G be a graph with a set of nodes V with m ≡ 2 elements, without isolated nodes. If 0 = 糸1 ≒ 糸1 ≒ ... ≒ 糸m are the eigenvalues of Lsym, prove the following properties: 
(1) We have 糸1 + 糸2 + ﹞﹞﹞ + 糸m = m. 

(2) We have 糸2 ≒ m/(m . 1), with equality holding i. G = Km, the complete graph on m nodes. 

(3) We have 糸m ≡ m/(m . 1). 

(4) If G is not a complete graph, then 糸2 ≒ 1 


Hint. If a and b are nonadjacent nodes, consider the function x given by 
db if v = a 
x(v)= 	.da if v = b 0 if v = a, b, 
..... 	(mm 
anduseProblem18.4(2). 	珧 
(5) Prove 	that 糸m ≒ 2. Prove that 糸m = 2 i. the underlying graph of G contains a nontrivial connected bipartite component. 
(
Hint. Use Problem 18.4(3). 
(6) Prove that if G is connected, then 糸2 > 0. 
Problem 18.6. Let G be a graph with a set of nodes V with m ≡ 2 elements, without isolated nodes. Let vol(G)=v﹋V dv and let 
v dvx(v) 
x =. 
vol(G) 
Prove that (x(u) . x(v))2 
糸2 = min u‵v 	. 
x=0 dv(x(v) . x)2 
v 
Problem 18.7. Let G be a connected bipartite graph. Prove that if 糸 is an eigenvalue of Lsym, then 2 . 糸 is also an eigenvalue of Lsym. 
Problem 18.8. Prove Proposition 18.7. 


