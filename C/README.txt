   trimmer.c 
   This program implements the O(n^2) dynamic programming algorithm of Tsourakakis et al. (ACM JEA '11)
   and some demos. 
   author: Charalampos E. Tsourakakis  

   Paper: Approximation Algorithms for Speeding up Dynamic Programming and Denoising aCGH data
   C. Tsourakakis, R. Peng, M. Tsiarli, G. Miller, R. Schwartz
   Howto compile: gcc -o trimmer trimmer.c and if it doesn't work (because of the math.h library) 
                  compile using gcc -lm -Wall  -o trimmer trimmer.c   
	Tested in Ubuntu  10.x and Windows 7
