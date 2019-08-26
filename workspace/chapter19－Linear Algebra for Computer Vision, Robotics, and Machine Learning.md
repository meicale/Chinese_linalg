Chapter 19 
Spectral Graph Drawing 
19.1 Graph Drawing and Energy Minimization 
Let G =(V, E) be some undirected graph. It is often desirable to draw a graph, usually in the plane but possibly in 3D, and it turns out that the graph Laplacian can be used to design surprisingly good methods. Say |V | = m. The idea is to assign a point ¦Ñ(vi) in Rn to the vertex vi ¡Ê V , for every vi ¡Ê V , and to draw a line segment between the points ¦Ñ(vi) and ¦Ñ(vj) i. there is an edge {vi,vj}. 
De.nition 19.1. Let G =(V, E) be some undirected graph with m vertices. A graph drawing is a function ¦Ñ: V ¡ú Rn, for some n ¡Ý 1. The matrix of a graph drawing ¦Ñ (in Rn) is a m ¡Á n matrix R whose ith row consists of the row vector ¦Ñ(vi) corresponding to the point representing vi in Rn . 
For a graph drawing to be useful we want n ¡Ü m; in fact n should be much smaller than m, typically n =2 or n = 3. 
De.nition 19.2. A graph drawing is balanced i. the sum of the entries of every column of the matrix of the graph drawing is zero, that is, 
1TR =0. 
If a graph drawing is not balanced, it can be made balanced by a suitable translation. We may also assume that the columns of R are linearly independent, since any basis of the column space also determines the drawing. Thus, from now on, we may assume that n ¡Ü m. 
Remark: A graph drawing ¦Ñ: V ¡ú Rn is not required to be injective, which may result in degenerate drawings where distinct vertices are drawn as the same point. For this reason, we prefer not to use the terminology graph embedding, which is often used in the literature. This is because in di.erential geometry, an embedding always refers to an injective map. The term graph immersion would be more appropriate. 
631 

As explained in Godsil and Royle [28], we can imagine building a physical model of G by connecting adjacent vertices (in Rn) by identical springs. Then it is natural to consider a representation to be better if it requires the springs to be less extended. We can formalize this by de.ning the energy of a drawing R by 
E(R)= £ü¦Ñ(vi) . ¦Ñ(vj)£ü2 , 
{vi,vj }¡ÊE
where ¦Ñ(vi) is the ith row of R and£ü¦Ñ(vi) . ¦Ñ(vj)£ü2 is the square of the Euclidean length of the line segment joining ¦Ñ(vi) and ¦Ñ(vj). 
Then, ¡°good drawings¡± are drawings that minimize the energy function E. Of course, the trivial representation corresponding to the zero matrix is optimum, so we need to impose extra constraints to rule out the trivial solution. 
We can consider the more general situation where the springs are not necessarily identical. This can be modeled by a symmetric weight (or sti.ness) matrix W =(wij), with wij ¡Ý 0. Then our energy function becomes 
 
E(R)=wij£ü¦Ñ(vi) . ¦Ñ(vj)£ü2 . 
{vi,vj }¡ÊE 
It turns out that this function can be expressed in terms of the Laplacian L = D . W . The following proposition is shown in Godsil and Royle [28]. We give a slightly more direct proof. 
Proposition 19.1. Let G =(V, W ) be a weighted graph, with |V | = m and W an m ¡Á m symmetric matrix, and let R be the matrix of a graph drawing ¦Ñ of G in Rn (a m¡Án matrix). If L = D . W is the unnormalized Laplacian matrix associated with W , then 
E(R) = tr(RTLR). 
Proof. Since ¦Ñ(vi) is the ith row of R (and ¦Ñ(vj) is the jth row of R), if we denote the kth column of R by Rk, using Proposition 18.4, we have 
 
E(R)=wij£ü¦Ñ(vi) . ¦Ñ(vj)£ü2 
{vi,vj }¡ÊE 
n
  
= wij(Rik . Rjk)2 
k=1{vi,vj }¡ÊE 
n1 m
  
= wij(Rik . Rjk)2 
2 
k=1 i,j=1 n
 
=(Rk)TLRk = tr(RTLR), k=1 
as claimed. 

19.1. GRAPH DRAWING AND ENERGY MINIMIZATION 
Since the matrix RTLR is symmetric, it has real eigenvalues. Actually, since L is positive semide.nite, so is RTLR. Then the trace of RTLR is equal to the sum of its positive eigenvalues, and this is the energy E(R) of the graph drawing. 
If R is the matrix of a graph drawing in Rn, then for any n ¡Á n invertible matrix M, the map that assigns ¦Ñ(vi)M to vi is another graph drawing of G, and these two drawings convey the same amount of information. From this point of view, a graph drawing is determined by the column space of R. Therefore, it is reasonable to assume that the columns of R are pairwise orthogonal and that they have unit length. Such a matrix satis.es the equation RTR = I. 
De.nition 19.3. If the matrix R of a graph drawing satis.es the equation RTR = I, then the corresponding drawing is called an orthogonal graph drawing. 
This above condition also rules out trivial drawings. The following result tells us how to .nd minimum energy orthogonal balanced graph drawings, provided the graph is connected. Recall that 
L1 =0, 
as we already observed. 
Theorem 19.2. Let G =(V, W ) be a weighted graph with |V | = m. If L = D . W is the (unnormalized) Laplacian of G, and if the eigenvalues of L are 0= ¦Ë1 <¦Ë2 ¡Ü ¦Ë3 ¡Ü ... ¡Ü ¦Ëm, then the minimal energy of any balanced orthogonal graph drawing of G in Rn is equal to ¦Ë2 +¡¤¡¤¡¤+¦Ën+1 (in particular, this implies that n<m). The m¡Án matrix R consisting of any unit eigenvectors u2,...,un+1 associated with ¦Ë2 ¡Ü ... ¡Ü ¦Ën+1 yields a balanced orthogonal graph drawing of minimal energy; it satis.es the condition RTR = I. 
Proof. We present the proof given in Godsil and Royle [28] (Section 13.4, Theorem 13.4.1). The key point is that the sum of the n smallest eigenvalues of L is a lower bound for tr(RTLR). This can be shown using a Rayleigh ratio argument; see Proposition 16.25 (the Poincar¡äe separation theorem). Then any n eigenvectors (u1,...,un) associated with ¦Ë1,...,¦Ën achieve this bound. Because the .rst eigenvalue of L is ¦Ë1 = 0 and because 
¡Ì 
we are assuming that ¦Ë2 > 0, we have u1 = 1/m. Since the uj are pairwise orthogonal 
¡Ì 
for i =2,...,n and since ui is orthogonal to u1 = 1/m, the entries in ui add up to 0. Consequently, for any f with 2 ¡Ü f ¡Ü n, by deleting u1 and using (u2,...,uR), we obtain a balanced orthogonal graph drawing in RR.1 with the same energy as the orthogonal graph drawing in RR using (u1,u2,...,uR). Conversely, from any balanced orthogonal drawing in 
RR.1 
using (u2,...,uR), we obtain an orthogonal graph drawing in RR using (u1,u2,...,uR) with the same energy. Therefore, the minimum energy of a balanced orthogonal graph drawing in Rn is equal to the minimum energy of an orthogonal graph drawing in Rn+1, and this minimum is ¦Ë2 + ¡¤¡¤¡¤ + ¦Ën+1. 
Since 1 spans the nullspace of L, using u1 (which belongs to Ker L) as one of the vectors in R would have the e.ect that all points representing vertices of G would have the same .rst coordinate. This would mean that the drawing lives in a hyperplane in Rn, which is undesirable, especially when n = 2, where all vertices would be collinear. This is why we omit the .rst eigenvector u1. 
Observe that for any orthogonal n ¡Á n matrix Q, since 
tr(RTLR) = tr(QTRTLRQ), 
the matrix RQ also yields a minimum orthogonal graph drawing. This amounts to applying the rigid motion QT to the rows of R. 
In summary, if ¦Ë2 > 0, an automatic method for drawing a graph in R2 is this: 
1. 
Compute the two smallest nonzero eigenvalues ¦Ë2 ¡Ü ¦Ë3 of the graph Laplacian L (it is possible that ¦Ë3 = ¦Ë2 if ¦Ë2 is a multiple eigenvalue); 

2. 
Compute two unit eigenvectors u2,u3 associated with ¦Ë2 and ¦Ë3, and let R =[u2 u3] be the m ¡Á 2 matrix having u2 and u3 as columns. 

3. 
Place vertex vi at the point whose coordinates is the ith row of R, that is, (Ri1,Ri2). 


This method generally gives pleasing results, but beware that there is no guarantee that distinct nodes are assigned distinct images since R can have identical rows. This does not seem to happen often in practice. 


19.2 Examples of Graph Drawings 
We now give a number of examples using Matlab. Some of these are borrowed or adapted from Spielman [60]. 
Example 1. Consider the graph with four nodes whose adjacency matrix is 
.
. 
...
1 0110 
We use the following program to compute u2 and u3: 
A =[0110;100 1;10 01;011 0]; 
D = diag(sum(A)); 
L =D -A; 
[v, e] = eigs(L); 
gplot(A, v(:,[3 2])) 
hold on; 
gplot(A, v(:,[3 2]),¡¯o¡¯) 

0110 1
100 

A = 

. 

100 

19.2. EXAMPLES OF GRAPH DRAWINGS 

The graph of Example 1 is shown in Figure 19.1. The function eigs(L) computes the six largest eigenvalues of L in decreasing order, and corresponding eigenvectors. It turns out that ¦Ë2 = ¦Ë3 = 2 is a double eigenvalue. 
Example 2. Consider the graph G2 shown in Figure 18.3 given by the adjacency matrix 
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

We use the following program to compute u2 and u3: 
A =[0110 0;10 11 1; 11 01 0; 01 101; 01010]; 
D = diag(sum(A)); 
L =D -A; 
[v, e] = eig(L); 
gplot(A, v(:, [2 3])) 
hold on 
gplot(A, v(:, [2 3]),¡¯o¡¯) 

The function eig(L) (with no s at the end) computes the eigenvalues of L in increasing order. The result of drawing the graph is shown in Figure 19.2. Note that node v2 is assigned to the point (0, 0), so the di.erence between this drawing and the drawing in Figure 18.3 is that the drawing of Figure 19.2 is not convex. 
Example 3. Consider the ring graph de.ned by the adjacency matrix A given in the Matlab program shown below: 

Figure 19.2: Drawing of the graph from Example 2. 
A = diag(ones(1, 11),1); 
A =A +A¡¯; 
A(1, 12) = 1; A(12, 1) = 1; 
D = diag(sum(A)); 
L =D -A; 
[v, e] = eig(L); 
gplot(A, v(:, [2 3])) 
hold on 
gplot(A, v(:, [2 3]),¡¯o¡¯) 


Observe that we get a very nice ring; see Figure 19.3. Again ¦Ë2 =0.2679 is a double eigenvalue (and so are the next pairs of eigenvalues, except the last, ¦Ë12 = 4). 
19.2. EXAMPLES OF GRAPH DRAWINGS 
Example 4. In this example adapted from Spielman, we generate 20 randomly chosen points in the unit square, compute their Delaunay triangulation, then the adjacency matrix of the corresponding graph, and .nally draw the graph using the second and third eigenvalues of the Laplacian. 
A = zeros(20,20); 
xy = rand(20, 2); 
trigs = delaunay(xy(:,1), xy(:,2)); 
elemtrig = ones(3) -eye(3); 
for i = 1:length(trigs), 

A(trigs(i,:),trigs(i,:)) = elemtrig; end A = double(A >0); gplot(A,xy) D = diag(sum(A)); L =D -A; [v, e] = eigs(L, 3, ¡¯sm¡¯); figure(2) gplot(A, v(:, [2 1])) hold on gplot(A, v(:, [2 1]),¡¯o¡¯) 
The Delaunay triangulation of the set of 20 points and the drawing of the corresponding graph are shown in Figure 19.4. The graph drawing on the right looks nicer than the graph on the left but is is no longer planar. 

Example 5. Our last example, also borrowed from Spielman [60], corresponds to the skeleton of the ¡°Buckyball,¡± a geodesic dome invented by the architect Richard Buckminster Fuller (1895¨C1983). The Montr¡äeal Biosph`ere is an example of a geodesic dome designed by Buckminster Fuller. 
A = full(bucky); 
D = diag(sum(A)); 
L =D -A; 
[v, e] = eig(L); 
gplot(A, v(:, [2 3])) 
hold on; 
gplot(A,v(:, [2 3]), ¡¯o¡¯) 

Figure 19.5 shows a graph drawing of the Buckyball. This picture seems a bit squashed for two reasons. First, it is really a 3-dimensional graph; second, ¦Ë2 =0.2434 is a triple eigenvalue. (Actually, the Laplacian of L has many multiple eigenvalues.) What we should really do is to plot this graph in R3 using three orthonormal eigenvectors associated with ¦Ë2. 

A 3D picture of the graph of the Buckyball is produced by the following Matlab program, and its image is shown in Figure 19.6. It looks better! 
[x, y] = gplot(A, v(:, [2 3])); [x, z] = gplot(A, v(:, [2 4])); plot3(x,y,z) 


19.3 Summary 
The main concepts and results of this chapter are listed below: 
. 
Graph drawing. 

19.3. SUMMARY 

. 
Matrix of a graph drawing. 

. 
Balanced graph drawing. 

. 
Energy E(R) of a graph drawing. 

. 
Orthogonal graph drawing. 

. 
Delaunay triangulation. 

. 
Buckyball. 





