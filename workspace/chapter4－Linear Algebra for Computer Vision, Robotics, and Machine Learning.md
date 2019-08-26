Chapter 4 
Haar Bases, Haar Wavelets, Hadamard Matrices 
In this chapter, we discuss two types of matrices that have applications in computer science and engineering: 
(1) Haar matrices and the corresponding Haar wavelets, a fundamental tool in signal pro-cessing and computer graphics. 
2) Hadamard matrices which have applications in error correcting codes, signal processing, and low rank approximation. 
4.1 	Introduction to Signal Compression Using Haar Wavelets 
We begin by considering Haar wavelets in R4 . Wavelets play an important role in audio and video signal processing, especially for compressing long signals into much smaller ones that still retain enough information so that when they are played, we can¡¯t see or hear any di.erence. 
Consider the four vectors w1,w2,w3,w4 given by 
..
..
..
.. 
11 	1 0 

w1 = 
... 

1 

1 

...

w2 = 
... 

1 

.1 

...

w3 = 
... 

.1 

0 

...

w4 = 
... 

0 

1 

...

. 

1 .1	0 .1 
Note that these vectors are pairwise orthogonal, so they are indeed linearly independent (we will see this in a later chapter). Let W = {w1,w2,w3,w4} be the Haar basis, and let U = {e1,e2,e3,e4} be the canonical basis of R4 . The change of basis matrix W = PW,U from 
101 

U to W is given by 

.
. 
11 1 0 
11 .10 

W = 

... 

...

,

1 .10 1 
1 .10 .1 

and we easily .nd that the inverse of W is given by 

.
..
. 
1/4000 1111 

W .1 
= 
... 

01/40 0 

0 01/20 

... 
... 

...

11 

.1 .1 

. 

1 

.1 

00 

0 0 01/2 001 .1 
So the vector v = (6, 4, 5, 1) over the basis U becomes c =(c1,c2,c3,c4) over the Haar basis W, with 
..
...
..
... 
c1 1/4000 1111 6 4 
... 

c2 
c3 
... 

= 

... 

01/40 0 

0 01/20 

... 
... 

... 
... 

4 

5 

... 

= 

... 

1 

1 

...

. 

11 

.1 .1 

1 

.1 

00 

c4 0 0 01/2 001 .11 2 
Given a signal v =(v1,v2,v3,v4), we .rst transform v into its coe.cients c =(c1,c2,c3,c4) over the Haar basis by computing c = W .1v. Observe that 
v1 + v2 + v3 + v4 
c1 = 
4 
is the overall average value of the signal v. The coe.cient c1 corresponds to the background of the image (or of the sound). Then, c2 gives the coarse details of v, whereas, c3 gives the details in the .rst part of v, and c4 gives the details in the second half of v. 
Reconstruction of the signal consists in computing v = Wc. The trick for good compres-sion is to throw away some of the coe.cients of c (set them to zero), obtaining a compressed signal ¦Åc, and still retain enough crucial information so that the reconstructed signal ¦Åv = W¦Åc looks almost as good as the original signal v. Thus, the steps are: 
input v .¡ú coe.cients c = W .1 v .¡ú compressed ¦Åc .¡ú compressed ¦Åv = W¦Åc. 
This kind of compression scheme makes modern video conferencing possible. 
It turns out that there is a faster way to .nd c = W .1v, without actually using W .1 . This has to do with the multiscale nature of Haar wavelets. 
Given the original signal v = (6, 4, 5, 1) shown in Figure 4.1, we compute averages and half di.erences obtaining Figure 4.2. We get the coe.cients c3 = 1 and c4 = 2. Then again we compute averages and half di.erences obtaining Figure 4.3. We get the coe.cients c1 = 4 and c2 = 1. Note that the original signal v can be reconstructed from the two signals 
4.2. HAAR MATRICES, SCALING PROPERTIES OF HAAR WAVELETS 


in Figure 4.2, and the signal on the left of Figure 4.2 can be reconstructed from the two signals in Figure 4.3. In particular, the data from Figure 4.2 gives us 
v1 + v2 v1 . v2
5+1= + = v1
22 
v1 + v2 v1 . v2
5 . 1= . = v2
22 
v3 + v4 v3 . v4
3+2= + = v3
22 
v3 + v4 v3 . v4
3 . 2= . = v4. 
22 

4.2 	Haar Bases and Haar Matrices, Scaling Properties of Haar Wavelets 
The method discussed in Section 4.2 can be generalized to signals of any length 2n . The previous case corresponds to n = 2. Let us consider the case n = 3. The Haar basis 

.
. 
W = 

........... 

11 1 0 1 0 0 0 11 1 0 .10 0 0 11 .10 0 1 0 0 11 .10 0 .10 0 
1  .1  0  1  0  0  1  0  
1  .1  0  1  0  0  .1  0  
1  .1  0  .1  0  0  0  1  

........... 

. 

1 .10 .10 0 0 .1 The columns of this matrix are orthogonal, and it is easy to see that W .1 = diag(1/8, 1/8, 1/4, 1/4, 1/2, 1/2, 1/2, 1/2)W T. A pattern is beginning to emerge. It looks like the second Haar basis vector w2 is the ¡°mother¡± of all the other basis vectors, except the .rst, whose purpose is to perform aver-aging. Indeed, in general, given 
,

_
_
w2 = (1,..., 1, .1,..., .1) 
2n 
the other Haar basis vectors are obtained by a ¡°scaling and shifting process.¡± Starting from w2, the scaling process generates the vectors 
w3,w5,w9,...,w2j +1,...,w2n.1+1, 
such that w2j+1+1 is obtained from w2j +1 by forming two consecutive blocks of 1 and .1 of half the size of the blocks in w2j +1, and setting all other entries to zero. Observe that w2j +1 has 2j blocks of 2n.j elements. The shifting process consists in shifting the blocks of 1 and .1 in w2j +1 to the right by inserting a block of (k . 1)2n.j zeros from the left, with 0 ¡Ü j ¡Ü n . 1 and 1 ¡Ü k ¡Ü 2j. Note that our convention is to use j as the scaling index and k as the shifting index. Thus, we obtain the following formula for w2j +k: 
w2j +k(i)= 
. .
. .
. 

0  1 ¡Ü i ¡Ü (k . 1)2n.j  
1  (k . 1)2n.j + 1 ¡Ü i ¡Ü (k . 1)2n.j + 2n.j.1  
.1  (k . 1)2n.j + 2n.j.1 + 1 ¡Ü i ¡Ü k2n.j  
0  k2n.j + 1 ¡Ü i ¡Ü 2n ,  

4.2. HAAR MATRICES, SCALING PROPERTIES OF HAAR WAVELETS 
with 0 ¡Ü j ¡Ü n . 1 and 1 ¡Ü k ¡Ü 2j. Of course 
w1 = (1,..., 1) . 2n 
The above formulae look a little better if we change our indexing slightly by letting k vary from 0 to 2j . 1, and using the index j instead of 2j. 
De.nition 4.1. The vectors of the Haar basis of dimension 2n are denoted by 
w1,h00,h10,h11 ,h20,h12 ,h22,h32 ,...,hkj ,...,h2nn..1 1.1, 
where 
. 
01 ¡Ü i ¡Ü k2n.j
.

.
1 k2n.j +1 ¡Ü i ¡Ü k2n.j +2n.j.1 
hj 
(i)= 
k.
.1 k2n.j +2n.j.1 +1 ¡Ü i ¡Ü (k + 1)2n.j . 0(k + 1)2n.j +1 ¡Ü i ¡Ü 2n , with 0 ¡Ü j ¡Ü n . 1 and 0 ¡Ü k ¡Ü 2j . 1. The2n ¡Á 2n matrix whose columns are the vectors w1,h00,h10,h11,h20,h21,h22,h23,...,hjk,...,hn2n..1 1.1, (in that order), is called the Haar matrix of dimension 2n, and is denoted by Wn. It turns out that there is a way to understand these formulae better if we interpret a vector u =(u1,...,um) as a piecewise linear function over the interval [0, 1). De.nition 4.2. Given a vector u =(u1,...,um), the piecewise linear function plf(u) is de.ned such that 
i . 1 i 
plf(u)(x)= ui, ¡Ü x< , 1 ¡Ü i ¡Ü m. 
mm In words, the function plf(u) has the value u1 on the interval [0, 1/m), the value u2 on [1/m, 2/m), etc., and the value um on the interval [(m . 1)/m, 1). For example, the piecewise linear function associated with the vector u = (2.4, 2.2, 2.15, 2.05, 6.8, 2.8, .1.1, .1.3) is shown in Figure 4.4. Then each basis vector hjk corresponds to the function ¦×kj = plf(hjk). In particular, for all n, the Haar basis vectors h00 = w2 = (1,..., 1, .1,..., .1) 2n 

yield the same piecewise linear function ¦× given by 
¦×(x)= 

. .. .. 

1  if  0 ¡Ü x < 1/2  
.1  if  1/2 ¡Ü x < 1  
0  otherwise,  

whose graph is shown in Figure 4.5. It is easy to see that ¦×

jk
is given by the simple expression 


¦×

jk
(x)= ¦×(2jx . k), 0 ¡Ü j ¡Ü n . 1, 0 ¡Ü k ¡Ü 2j . 1. 
The above formula makes it clear that ¦×

jk
is obtained from ¦× by scaling and shifting. 

De.nition 4.3. The function ¦Õ00 = plf(w1) is the piecewise linear function with the constant 
value 1 on [0, 1), and the functions ¦×

jk 
= 

plf(h

jk
) together with ¦Õ00 
are known as the Haar 

wavelets. 
Rather than using W .1 to convert a vector u to a vector c of coe.cients over the Haar basis, and the matrix W to reconstruct the vector u from its Haar coe.cients c, we can use faster algorithms that use averaging and di.erencing. 

4.2. HAAR MATRICES, SCALING PROPERTIES OF HAAR WAVELETS 
If c is a vector of Haar coe.cients of dimension 2n, we compute the sequence of vectors u0,u1 ,..., un as follows: 
0 
u = c 
j+1 j
u= uuj+1(2i . 1) = uj(i)+ uj(2j + i) uj+1(2i)= uj(i) . uj(2j + i), 
n
for j =0,...,n . 1 and i =1,..., 2j. The reconstructed vector (signal) is u = u. 
nn.10
If u is a vector of dimension 2n, we compute the sequence of vectors c,c,...,cas follows: 
n 
c = u 
jj+1 
c= ccj(i)=(cj+1(2i . 1) + cj+1(2i))/2 cj(2j + i)=(cj+1(2i . 1) . cj+1(2i))/2, 
for j = n . 1,..., 0 and i =1,..., 2j. The vector over the Haar basis is c = c0 . 
We leave it as an exercise to implement the above programs in Matlab using two variables u and c, and by building iteratively 2j. Here is an example of the conversion of a vector to its Haar coe.cients for n = 3. 
Given the sequence u = (31, 29, 23, 17, .6, .8, .2, .4), we get the sequence 
c 3 = (31, 29, 23, 17, .6, .8, 2, .4)
 
2 31+29 23+17 .6 . 8 .2 . 4 31 . 29 23 . 17 .6 . (.8) 
c =,,,,,, ,
222222 2 
 
.2 . (.4) 
2
= (30, 20, .7, .3, 1, 3, 1, 1)

  
1 30 + 20 .7 . 3 30 . 20 .7 . (.3) 
c =,,, , 1, 3, 1, 1
222 2 
= (25, .5, 5, .2, 1, 3, 1, 1)
  
25 . 5 25 . (.5) 
c 0 =,, 5, .2, 1, 3, 1, 1
22 
= (10, 15, 5, .2, 1, 3, 1, 1) 
so c = (10, 15, 5, .2, 1, 3, 1, 1). Conversely, given c = (10, 15, 5, .2, 1, 3, 1, 1), we get the 

sequence 
u 0 = (10, 15, 5, .2, 1, 3, 1, 1) 
u 1 = (10 + 15, 10 . 15, 5, .2, 1, 3, 1, 1) = (25, .5, 5, .2, 1, 3, 1, 1) 
u 2 = (25+5, 25 . 5, .5+(.2), .5 . (.2), 1, 3, 1, 1) 

= (30, 20, .7, .3, 1, 3, 1, 1) u 3 = (30+1, 30 . 1, 20 + 3, 20 . 3, .7+1, .7 . 1, .3+1, .3 . 1) = (31, 29, 23, 17, .6, .8, .2, .4), 
which gives back u = (31, 29, 23, 17, .6, .8, .2, .4). 


4.3 Kronecker Product Construction of Haar Matrices 
There is another recursive method for constructing the Haar matrix Wn of dimension 2n that makes it clearer why the columns of Wn are pairwise orthogonal, and why the above algorithms are indeed correct (which nobody seems to prove!). If we split Wn into two 2n ¡Á 2n.1 matrices, then the second matrix containing the last 2n.1 columns of Wn has a very simple structure: it consists of the vector 
(1, .1, 0,..., 0) 
_

2n
_ 

and 2n.1 . 1 shifted copies of it, as illustrated below for n = 3: 
.
. 
........... 

1000 .10 0 0 0100 0 .10 0 
0  0  1  0  
0  0  .1  0  
0  0  0  1  
0  0  0  .1  

........... 

. 

Observe that this matrix can be obtained from the identity matrix I2n.1 , in our example 
.
. 
I4 = 
... 

1000 
0100 

0010 
0001 

...

, 

by forming the 2n ¡Á 2n.1 matrix obtained by replacing each 1 by the column vector 
1 .1 
4.3. KRONECKER PRODUCT CONSTRUCTION OF HAAR MATRICES 
and each zero by the column vector 
0 
. 
0 
Now the .rst half of Wn, that is the matrix consisting of the .rst 2n.1 columns of Wn, can be obtained from Wn.1 by forming the 2n ¡Á 2n.1 matrix obtained by replacing each 1 by the column vector 
1 
,
1 each .1 by the column vector 
.1 
,
.1 and each zero by the column vector 
0 
. 
0 
For n = 3, the .rst half of W3 is the matrix 
.
. 
........... 

1  1  1  0  
1  1  1  0  
1  1  .1  0  
1  1  .1  0  
1  .1  0  1  
1  .1  0  1  
1  .1  0  .1  
1  .1  0  .1  

........... 

which is indeed obtained from 

.
. 
... 

11 1 0 
11 .10 

1 .10 1 
1 .10 .1 

...

W2 = 
using the process that we just described. 
These matrix manipulations can be described conveniently using a product operation on matrices known as the Kronecker product. 
De.nition 4.4. Given a m¡Án matrix A =(aij) and a p¡Áq matrix B =(bij), the Kronecker product (or tensor product) A . B of A and B is the mp ¡Á nq matrix 
.
. 
A . B = 

.... 

a11Ba12B ¡¤¡¤¡¤ a1nB a21Ba22B ¡¤¡¤¡¤ a2nB 
.. .
.
. ...
.
.. . 
.... 

. 

am1Bam2B ¡¤¡¤¡¤ amnB 
It can be shown that . is associative and that 
(A . B)(C . D)= AC . BD (A . B)T = AT . BT , 
whenever AC and BD are well de.ned. Then it is immediately veri.ed that Wn is given by the following neat recursive equations: 
11 
Wn = Wn.1 . I2n.1 . ,
1 .1 
with W0 = (1). If we let 10 20 
B1 =2 = 
01 02 
and for n ¡Ý 1, Bn 0 
Bn+1 =2 ,
0 I2n 
then it is not hard to use the Kronecker product formulation of Wn to obtain a rigorous proof of the equation 
Wn TWn = Bn, for all n ¡Ý 1. 
The above equation o.ers a clean justi.cation of the fact that the columns of Wn are pairwise orthogonal. 
Observe that the right block (of size 2n ¡Á 2n.1) shows clearly how the detail coe.cients in the second half of the vector c are added and subtracted to the entries in the .rst half of the partially reconstructed vector after n . 1 steps. 


4.4 Multiresolution Signal Analysis with Haar Bases 
An important and attractive feature of the Haar basis is that it provides a multiresolution analysis of a signal. Indeed, given a signal u, if c =(c1,...,c2n ) is the vector of its Haar coef-.cients, the coe.cients with low index give coarse information about u, and the coe.cients with high index represent .ne information. For example, if u is an audio signal corresponding to a Mozart concerto played by an orchestra, c1 corresponds to the ¡°background noise,¡± c2 to the bass, c3 to the .rst cello, c4 to the second cello, c5,c6,c7,c7 to the violas, then the violins, etc. This multiresolution feature of wavelets can be exploited to compress a signal, that is, to use fewer coe.cients to represent it. Here is an example. 
Consider the signal 
u = (2.4, 2.2, 2.15, 2.05, 6.8, 2.8, .1.1, .1.3), 
whose Haar transform is c = (2, 0.2, 0.1, 3, 0.1, 0.05, 2, 0.1). 
4.4. MULTIRESOLUTION SIGNAL ANALYSIS WITH HAAR BASES 
The piecewise-linear curves corresponding to u and c are shown in Figure 4.6. Since some of the coe.cients in c are small (smaller than or equal to 0.2) we can compress c by replacing them by 0. We get 
c2 = (2, 0, 0, 3, 0, 0, 2, 0), 
and the reconstructed signal is 
u2 = (2, 2, 2, 2, 7, 3, .1, .1). 
The piecewise-linear curves corresponding to u2 and c2 are shown in Figure 4.7. 


An interesting (and amusing) application of the Haar wavelets is to the compression of audio signals. It turns out that if your type load handel in Matlab an audio .le will be loaded in a vector denoted by y, and if you type sound(y), the computer will play this piece of music. You can convert y to its vector of Haar coe.cients c. The length of y is 73113, 

216
so .rst tuncate the tail of y to get a vector of length 65536 = . A plot of the signals corresponding to y and c is shown in Figure 4.8. Then run a program that sets all coe.cients of c whose absolute value is less that 0.05 to zero. This sets 37272 coe.cients to 0. The resulting vector c2 is converted to a signal y2. A plot of the signals corresponding to y2 and c2 is shown in Figure 4.9. When you type sound(y2), you .nd that the music doesn¡¯t di.er 

much from the original, although it sounds less crisp. You should play with other numbers greater than or less than 0.05. You should hear what happens when you type sound(c). It plays the music corresponding to the Haar transform c of y, and it is quite funny. 

4.5. HAAR TRANSFORM FOR DIGITAL IMAGES 


4.5 Haar Transform for Digital Images 
Another neat property of the Haar transform is that it can be instantly generalized to matrices (even rectangular) without any extra e.ort! This allows for the compression of digital images. But .rst we address the issue of normalization of the Haar coe.cients. As we observed earlier, the 2n ¡Á 2n matrix Wn of Haar basis vectors has orthogonal columns, but its columns do not have unit length. As a consequence, Wn T is not the inverse of Wn, but rather the matrix 
W .1 W T 
n = Dnn 
with 
 
2.n , 2.n , 2.(n.1), 2.(n.1), 2.(n.2)
Dn = diag,..., 2.(n.2),..., 
20 21 22 
 
2.1 
,..., 2.1 . 
2n.1
De.nition 4.5. The orthogonal matrix 
1 
Hn = WnDn 2 
whose columns are the normalized Haar basis vectors, with 
  
1 2. n 2 , 2. n , 2. n.1 , 2. n.1 , 2. n.2 
2 

22222 2
Dn = diag,..., 2. n.2 ,..., 2. 21 ,..., 2. 1 
20 21 22 2n.1
is called the normalized Haar transform matrix. Given a vector (signal) u, we call c = Hn Tu the normalized Haar coe.cients of u. 
Because Hn is orthogonal, Hn .1 = Hn T . 
Then a moment of re.ection shows that we have to slightly modify the algorithms to compute Hn Tu and Hnc as follows: When computing the sequence of ujs, use 
¡Ì 
uj+1(2i . 1) = (uj(i)+ uj(2j + i))/ 2 
¡Ì 
uj+1(2i)=(uj(i) . uj(2j + i))/ 2, 
j
and when computing the sequence of cs, use 
¡Ì 
cj(i)=(cj+1(2i . 1) + cj+1(2i))/ 2 
¡Ì 
cj(2j + i)=(cj+1(2i . 1) . cj+1(2i))/ 2. 
¡Ì 
Note that things are now more symmetric, at the expense of a division by 2. However, for long vectors, it turns out that these algorithms are numerically more stable. 
¡Ì 
Remark: Some authors (for example, Stollnitz, Derose and Salesin [62]) rescale c by 1/ 2n 
¡Ì 
and u by 2n . This is because the norm of the basis functions ¦×j is not equal to 1 (under 
 1 k 
the inner product (f, g£© =0 f(t)g(t)dt). The normalized basis functions are the functions 
¡Ì 
2j¦×kj . Let us now explain the 2D version of the Haar transform. We describe the version using the matrix Wn, the method using Hn being identical (except that Hn .1 = Hn T, but this does not hold for Wn .1). Givena2m ¡Á 2n matrix A, we can .rst convert the rows of A to their Haar coe.cients using the Haar transform Wn .1, obtaining a matrix B, and then convert the columns of B to their Haar coe.cients, using the matrix W .1 . Because columns and rows 
m 
are exchanged in the .rst step, B = A(Wn .1)T , 
= W .1
and in the second step C m B, thus, we have 
C = W .1A(W .1)T = DmW TAWn Dn.
mn m 
In the other direction, given a 2m ¡Á 2n matrix C of Haar coe.cients, we reconstruct the matrix A (the image) by .rst applying Wm to the columns of C, obtaining B, and then Wn T to the rows of B. Therefore 
A = WmCW n T . 
Of course, we don¡¯t actually have to invert Wm and Wn and perform matrix multiplications. We just have to use our algorithms using averaging and di.erencing. Here is an example. 
If the data matrix (the image) is the 8 ¡Á 8 matrix 
.
. 
A = 

........... 

64  2  3  61  60  6  7  57  
9  55  54  12  13  51  50  16  
17  47  46  20  21  43  42  24  
40  26  27  37  36  30  31  33  

32  34  35  29  28  38  39  25  
41  23  22  44  45  19  18  48  
49  15  14  52  53  11  10  56  
8  58  59  5  4  62  63  1  

........... 

, 

then applying our algorithms, we .nd that 

.
. 
C = 

........... 

32.500 0 0000 000 0 0000 000 04 .44 .4 000 04 .44 .4 
0  0  0.5  0.5  27  .25  23  .21  
0  0  .0.5  .0.5  .11  9  .7  5  
0  0  0.5  0.5  .5  7  .9  11  
0  0  .0.5  .0.5  21  .23  25  .27  

........... 

. 

4.5. HAAR TRANSFORM FOR DIGITAL IMAGES 
As we can see, C has more zero entries than A; it is a compressed version of A. We can further compress C by setting to 0 all entries of absolute value at most 0.5. Then we get 
.
. 
C2 = 
........... 

32.50000 000 00000 000 0 000 4 .44 .4 0 000 4 .44 .4 
0  0  0  0  27  .25  23  .21  
0  0  0  0  .11  9  .7  5  
0  0  0  0  .5  7  .9  11  
0  0  0  0  21  .23  25  .27  

........... 

. 

We .nd that the reconstructed image is 

.
. 
A2 = 
........... 

63.5  1.5  3.5  61.5  59.5  5.5  7.5  57.5  
9.5  55.5  53.5  11.5  13.5  51.5  49.5  15.5  
17.5  47.5  45.5  19.5  21.5  43.5  41.5  23.5  
39.5  25.5  27.5  37.5  35.5  29.5  31.5  33.5  

31.5  33.5  35.5  29.5  27.5  37.5  39.5  25.5  
41.5  23.5  21.5  43.5  45.5  19.5  17.5  47.5  
49.5  15.5  13.5  51.5  53.5  11.5  9.5  55.5  
7.5  57.5  59.5  5.5  3.5  61.5  63.5  1.5  

........... 

, 

which is pretty close to the original image matrix A. 
It turns out that Matlab has a wonderful command, image(X) (also imagesc(X), which often does a better job), which displays the matrix X has an image in which each entry is shown as a little square whose gray level is proportional to the numerical value of that entry (lighter if the value is higher, darker if the value is closer to zero; negative values are treated as zero). The images corresponding to A and C are shown in Figure 4.10. The compressed images corresponding to A2 and C2 are shown in Figure 4.11. The compressed versions appear to be indistinguishable from the originals! 


If we use the normalized matrices Hm and Hn, then the equations relating the image matrix A and its normalized Haar transform C are 
C = HTAHn
m
CHT .
A = Hmn 
The Haar transform can also be used to send large images progressively over the internet. Indeed, we can start sending the Haar coe.cients of the matrix C starting from the coarsest coe.cients (the .rst column from top down, then the second column, etc.), and at the receiving end we can start reconstructing the image as soon as we have received enough data. 
Observe that instead of performing all rounds of averaging and di.erencing on each row and each column, we can perform partial encoding (and decoding). For example, we can perform a single round of averaging and di.erencing for each row and each column. The result is an image consisting of four subimages, where the top left quarter is a coarser version of the original, and the rest (consisting of three pieces) contain the .nest detail coe.cients. We can also perform two rounds of averaging and di.erencing, or three rounds, etc. The second round of averaging and di.erencing is applied to the top left quarter of the image. Generally, the kth round is applied to the 2m+1.k ¡Á 2n+1.k submatrix consisting of the .rst 
2m+1.k 
rows and the .rst 2n+1.k columns (1 ¡Ü k ¡Ü n) of the matrix obtained at the end of the previous round. This process is illustrated on the image shown in Figure 4.5. The result of performing one round, two rounds, three rounds, and nine rounds of averaging is shown in Figure 4.13. Since our images have size 512 ¡Á 512, nine rounds of averaging yields the Haar transform, displayed as the image on the bottom right. The original image has completely disappeared! We leave it as a fun exercise to modify the algorithms involving averaging and di.erencing to perform k rounds of averaging/di.erencing. The reconstruction algorithm is 
4.5. HAAR TRANSFORM FOR DIGITAL IMAGES 

a little tricky. 
A nice and easily accessible account of wavelets and their uses in image processing and computer graphics can be found in Stollnitz, Derose and Salesin [62]. A very detailed account is given in Strang and and Nguyen [66], but this book assumes a fair amount of background in signal processing. 
22n
We can .nd easily a basis of 2n ¡Á 2n = vectors wij (2n ¡Á 2n matrices) for the linear map that reconstructs an image from its Haar coe.cients, in the sense that for any 2n ¡Á 2n matrix C of Haar coe.cients, the image matrix A is given by 
2n2n
22 
A = cijwij. i=1 j=1 
Indeed, the matrix wij is given by the so-called outer product 
wij = wi(wj)T . 
22n
Similarly, there is a basis of 2n ¡Á 2n = vectors hij (2n ¡Á 2n matrices) for the 2D Haar transform, in the sense that for any 2n ¡Á 2n matrix A, its matrix C of Haar coe.cients is given by 
2n2n
22 
C = aijhij. i=1 j=1 
  
If the columns of W .1 are w1,...,w2n , then   )T
hij = w(w.
ij

We leave it as exercise to compute the bases (wij) and (hij) for n = 2, and to display the corresponding images using the command imagesc. 


4.6 Hadamard Matrices 
There is another famous family of matrices somewhat similar to Haar matrices, but these matrices have entries +1 and .1 (no zero entries). 
De.nition 4.6. A real n ¡Á n matrix H is a Hadamard matrix if hij = ¡À1 for all i, j such that 1 ¡Ü i, j ¡Ü n and if 
HTH = nIn. 
Thus the columns of a Hadamard matrix are pairwise orthogonal. Because H is a square matrix, the equation HTH = nIn shows that H is invertible, so we also have HHT = nIn. 
4.6. HADAMARD MATRICES 
The following matrices are example of Hadamard matrices: 
.
. 
11 1 1 

... 

...

11 

1 

.11 .1 

11 

.1 .1 

H2 = ,H4 = 
,

1 

.1 

1 .1 .11 

and

.
. 
........... 

11 1 1 1 1 1 1 1 .11 .11 .11 .1 11 .1 .11 1 .1 .1 1 .1 .11 1 .1 .11 
1  1  1  1  .1  .1  .1  .1  
1  .1  1  .1  .1  1  .1  1  
1  1  .1  .1  .1  .1  1  1  
1  .1  .1  1  .1  1  1  .1  

........... 

H8 = 
. 

A natural question is to determine the positive integers n for which a Hadamard matrix of dimension n exists, but surprisingly this is an open problem. The Hadamard conjecture is that for every positive integer of the form n =4k, there is a Hadamard matrix of dimension 
n. What is known is a necessary condition and various su.cient conditions. 
Theorem 4.1. If H is an n ¡Á n Hadamard matrix, then either n =1, 2, or n =4k for some positive integer k. 
Sylvester introduced a family of Hadamard matrices and proved that there are Hadamard matrices of dimension n =2m for all m ¡Ý 1 using the following construction. 
Proposition 4.2. (Sylvester, 1867) If H is a Hadamard matrix of dimension n, then the block matrix of dimension 2n, 
HH 
,
H .H 
is a Hadamard matrix. 
If we start with 
11 
H2 = ,
1 .1 
we obtain an in.nite family of symmetric Hadamard matrices usually called Sylvester¨C Hadamard matrices and denoted by H2m . The Sylvester¨CHadamard matrices H2,H4 and H8 are shown on the previous page. 
In 1893, Hadamard gave examples of Hadamard matrices for n = 12 and n = 20. At the present, Hadamard matrices are known for all n =4k ¡Ü 1000, except for n = 668, 716, and 
892. 
Hadamard matrices have various applications to error correcting codes, signal processing, and numerical linear algebra; see Seberry, Wysocki and Wysocki [56] and Tropp [69]. For example, there is a code based on H32 that can correct 7 errors in any 32-bit encoded block, and can detect an eighth. This code was used on a Mariner spacecraft in 1969 to transmit pictures back to the earth. 
For every m ¡Ý 0, the piecewise a.ne functions plf((H2m )i) associated with the 2m rows of the Sylvester¨CHadamard matrix H2m are functions on [0, 1] known as the Walsh functions. It is customary to index these 2m functions by the integers 0, 1,..., 2m .1 in such a way that the Walsh function Wal(k, t) is equal to the function plf((H2m )i) associated with the Row i of H2m that contains k changes of signs between consecutive groups of +1 and consecutive groups of .1. For example, the .fth row of H8, namely
 	 
1 .1 .111 .1 .11, 
has .ve consecutive blocks of +1s and .1s, four sign changes between these blocks, and thus is associated with Wal(4,t). In particular, Walsh functions corresponding to the rows of H8 (from top down) are: 
Wal(0,t), Wal(7,t), Wal(3,t), Wal(4,t), 
Wal(1,t), Wal(6,t), Wal(2,t), Wal(5,t). 
Because of the connection between Sylvester¨CHadamard matrices and Walsh functions, Sylvester¨CHadamard matrices are called Walsh¨CHadamard matrices by some authors. For every m, the 2m Walsh functions are pairwise orthogonal. The countable set of Walsh functions Wal(k, t) for all m ¡Ý 0 and all k such that 0 ¡Ü k ¡Ü 2m . 1 can be ordered in such a way that it is an orthogonal Hilbert basis of the Hilbert space L2([0, 1)]; see Seberry, Wysocki and Wysocki [56]. 
The Sylvester¨CHadamard matrix H2m plays a role in various algorithms for dimension reduction and low-rank matrix approximation. There is a type of structured dimension-reduction map known as the subsampled randomized Hadamard transform, for short SRHT; see Tropp [69] and Halko, Martinsson and Tropp [33]. For f¡¶ n =2m, an SRHT matrix is an f ¡Á n matrix of the form 
 
n
¦µ=fRHD, 
where 
1. 	
D is a random n ¡Á n diagonal matrix whose entries are independent random signs. 

2. 	
H = n.1/2Hn, a normalized Sylvester¨CHadamard matrix of dimension n. 

3. 	
R is a random f ¡Á n matrix that restricts an n-dimensional vector to f coordinates, chosen uniformly at random. 


It is explained in Tropp [69] that for any input x such that lxl2 = 1, the probability that 
i 
|(HDx)i|¡Ý n.1 log(n) for any i is quite small. Thus HD has the e.ect of ¡°.attening¡± the input x. The main result about the SRHT is that it preserves the geometry of an entire subspace of vectors; see Tropp [69] (Theorem 1.3). 

4.7. SUMMARY 


4.7 Summary 
The main concepts and results of this chapter are listed below: 
. 
Haar basis vectors and a glimpse at Haar wavelets. 

. 
Kronecker product (or tensor product) of matrices. 

. 
Hadamard and Sylvester¨CHadamard matrices. 

. 
Walsh functions. 



4.8 Problems 
Problem 4.1. (Haar extravaganza) Consider the matrix 
.
. 
W3,3 = 
........... 

10001 0 0 0 1000 .10 0 0 01000 1 0 0 0100 0 .10 0 
0  0  1  0  0  0  1  0  
0  0  1  0  0  0  .1  0  
0  0  0  1  0  0  0  1  
0  0  0  1  0  0  0  .1  

........... 

. 

(1) Show that given any vector c =(c1,c2,c3,c4,c5,c6,c7,c8), the result W3,3c of applying W3,3 to c is 
W3,3c =(c1 + c5,c1 . c5,c2 + c6,c2 . c6,c3 + c7,c3 . c7,c4 + c8,c4 . c8), 
the last step in reconstructing a vector from its Haar coe.cients. 
(2) Prove that the inverse of W3,3 is (1/2)W3T ,3. Prove that the columns and the rows of W3,3 are orthogonal. 
(3) Let W3,2 and W3,1 be the following matrices: 
.
.
.
. 
........... 
101 00000 11000000 0 0 0 0 0 
W3,2 
= 

........... 

........... 

........... 

10 

.1 

0 0000 

1 

.1 

00000 

010 10000 

0 0 10000 

01 0 .10000 

0 0 01000 

W3,1 
=

, 

. 

000 01000 

0 0 00100 

000 00100 

0 0 00010 

000 00010 

0 0 00001 

0 000 00001 00000001 Show that given any vector c =(c1,c2,c3,c4,c5,c6,c7,c8), the result W3,2c of applying W3,2 to c is 
W3,2c =(c1 + c3,c1 . c3,c2 + c4,c2 . c4,c5,c6,c7,c8), 
the second step in reconstructing a vector from its Haar coe.cients, and the result W3,1c of applying W3,1 to c is 
W3,1c =(c1 + c2,c1 . c2,c3,c4,c5,c6,c7,c8), 
the .rst step in reconstructing a vector from its Haar coe.cients. Conclude that 
W3,3W3,2W3,1 = W3, 
the Haar matrix 

.
. 
........... 

11 1 0 1 0 0 0 11 1 0 .10 0 0 11 .10 0 1 0 0 11 .10 0 .10 0 
1  .1  0  1  0  0  1  0  
1  .1  0  1  0  0  .1  0  
1  .1  0  .1  0  0  0  1  
1  .1  0  .1  0  0  0  .1  

........... 

W3 = 
. 

Hint. First check that 
W2 04,4
W3,2W3,1 = ,
04,4 I4 
where

.
. 
... 

11 1 0 
11 .10 

1 .10 1 
1 .10 .1 

...

W2 = 
. 

(4) Prove that the columns and the rows of W3,2 and W3,1 are orthogonal. Deduce from this that the columns of W3 are orthogonal, and the rows of W3 .1 are orthogonal. Are the rows of W3 orthogonal? Are the columns of W3 .1 orthogonal? Find the inverse of W3,2 and the inverse of W3,1. 
Problem 4.2. This is a continuation of Problem 4.1. 
(1) For any n ¡Ý 2, the 2n ¡Á 2n matrix Wn,n is obtained form the two rows 
1, 0,..., 0, 1, 0,..., 0 
_

2n.1
_
_

2n.1
_ 

1, 0,..., 0, .1, 0,..., 0 

_

2n.1
_
_

2n.1
_ 

by shifting them 2n.1 . 1 times over to the right by inserting a zero on the left each time. 
4.8. PROBLEMS 
Given any vector c =(c1,c2,...,c2n ), show that Wn,nc is the result of the last step in the process of reconstructing a vector from its Haar coe.cients c. Prove that W .1 = (1/2)W T ,
n,n n,nand that the columns and the rows of Wn,n are orthogonal. 
(2) Given a m ¡Á n matrix A =(aij) and a p ¡Á q matrix B =(bij), the Kronecker product (or tensor product) A . B of A and B is the mp ¡Á nq matrix 
.
. 
A . B = 

.... 

a11Ba12B ¡¤¡¤¡¤ a1nB a21Ba22B ¡¤¡¤¡¤ a2nB 
.. .
.
. ...
.
.. . 
am1Bam2B ¡¤¡¤¡¤ amnB
.... 

. 

It can be shown (and you may use these facts without proof) that . is associative and that (A . B)(C . D)= AC . BD (A . B)T = AT . BT , whenever AC and BD are well de.ned. 
Check that  
Wn,n =  I2n.1 .  1 1  I2n.1 .  1 .1  ,  
and that  
Wn =  Wn.1 .  1 1  I2n.1 .  1 .1  .  
Use the above to reprove that  

W T 
=2I2n .
Wn,nn,n 
Let 
10 20 
B1 =2 = 
01 02 
and for n ¡Ý 1, 
Bn 0 
Bn+1 =2 . 
0 I2n Prove that 
Wn TWn = Bn, for all n ¡Ý 1. 
(3) The matrix Wn,i is obtained from the matrix Wi,i (1 ¡Ü i ¡Ü n . 1) as follows: 
Wi,i 02i ,2n.2i 
Wn,i = . 
02n.2i ,2i I2n.2i 
It consists of four blocks, where 02i ,2n.2i and 02n.2i ,2i are matrices of zeros and I2n.2i is the identity matrix of dimension 2n . 2i . 
Explain what Wn,i does to c and prove that 
Wn,nWn,n.1 ¡¤¡¤¡¤ Wn,1 = Wn, 
where Wn is the Haar matrix of dimension 2n . 
Hint. Use induction on k, with the induction hypothesis 

Wk 02k ,2n.2k Wn,kWn,k.1 ¡¤¡¤¡¤ Wn,1 = . 
02n.2k ,2k I2n.2k 
Prove that the columns and rows of Wn,k are orthogonal, and use this to prove that the columns of Wn and the rows of Wn .1 are orthogonal. Are the rows of Wn orthogonal? Are the columns of Wn .1 orthogonal? Prove that 
1 W T 
W .12 k,k 02k ,2n.2k 
= .
n,k 02n.2k ,2k I2n.2k 
Problem 4.3. Prove that if H is a Hadamard matrix of dimension n, then the block matrix of dimension 2n, HH 
,
H .H 
is a Hadamard matrix. 
Problem 4.4. Plot the graphs of the eight Walsh functions Wal(k, t) for k =0, 1,..., 7. 
Problem 4.5. Describe a recursive algorithm to compute the product H2m x of the Sylvester¨C Hadamard matrix H2m by a vector x ¡Ê R2m that uses m recursive calls. 



