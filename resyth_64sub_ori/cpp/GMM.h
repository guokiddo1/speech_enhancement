
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
using std::vector;

// This structure is for one mixture
struct Mix 
{
   float *Mean;  // mean vector
   float *variance; // variance vector
   float weight;        // weight
   vector<int> Sample_Index; // This vector stores the indices of data points which belongs to this mixture  
   int Total_Sample;   // This counts the total number of samples of this mixture
};

class GMM   // This class creates a GMM
{
 public:
   GMM(int,int,int);   // Constructor
   void Read_Feature(double*);  // Read in features from a file
   void Random_Init();                 // Randomly assign initial centers
   void Iterate();                     // Do k-mean iterations
   void EM();                          // Do EM
   float* GetMean(int i){return Mixtures[i].Mean;}
   float* GetVar(int i){return Mixtures[i].variance;}
   float GetWeight(int i){return Mixtures[i].weight;}
   ~GMM();

 private:
	int DIM;        // This is data dimensionality
	int Num_Mix;     // Number of GMM mixtures
	int NumData;   // This is the number of data to be used
	int Min_samples;  // This is the minimum number of samples in each cluster
	float Min_weight; // This is the weight floor for k-mean
	int Max_Iter; // This is the maximum number of iterations for k-mean
	float threshold; // This is the stop criteria for k-mean
	float Infinity ;
	float Perturb; // When a bad cell appears, replace its center by the center
                          // of the most popular cell plus this perturbation
	float threshold_EM; // stop criteria for EM
   	float **Feature;  // Each row is a feature
	Mix *Mixtures;        // Each element is a mixture
	float **E_Prob; // This matrix stores Pr(class i|data j) in E step
	void Recluster();             // Recluster samples in k-mean
	int Check_Badcell();        // Check if there is a mixture with less than Min_samples in k-mean 
   void Adjust_Badmean(int);   // If there is a bad cell, add a new center in the most popular cluster
   float  Update_Mean();           // Update centers after reclustering
   void Update_variance();         // update variances after k-mean
   void Update_weight();           // update weights after k-mean
   float Distance(int, int);     // Compute L2-norm
   int find_closest_mixture(int);  // find closest mixture to a pattern in k-mean
   float Compute_Gaussian(int,float *); // Compute gaussian for a mixture
   void E_step();                       // This performs E-step in EM
   float M_step();                       // This performs M-step in EM
};
// ---------------------------------------------------End  Class Definition ---------------------------------------------------- 

