#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "momet.h"
#include "hv.h"

#define INF 1e40
#define EPS 1e-6

Momet::Momet() {
	
}

Momet::~Momet() {
	
}

double Momet::hypervolume(vector<vector<double> > PFKnown, vector<double> reference) {
	PFKnown = removeDominated(PFKnown);
	
	int n = PFKnown.size();
	int d = PFKnown[0].size();
	
	double *ref = new double[d];
	double *data = new double[n * d];
	for (int i = 0; i < d; i++)
		ref[i] = reference[i];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < d; j++)
			data[i * d + j] = PFKnown[i][j];
	
	double hv = fpli_hv(data, d, n, ref);

	delete [] data;
	delete [] ref;
	
	return hv;
}

double Momet::bestAchievement(vector<vector<double> >& PFknown, int &pos, vector<double> zref, vector<double> weights) {
	pos = 1;
	double bestAch = chebyAchievement(PFknown.at(0), zref, weights);
	for (unsigned int i = 1; i < PFknown.size(); ++i) {
		double currentAch = chebyAchievement(PFknown.at(i), zref, weights);
		if (currentAch < bestAch) {
			bestAch = currentAch;
			pos = i + 1;
		}
	}

	return bestAch;
}

double Momet::chebyAchievement(vector<double> x, vector<double> zref, vector<double> weights) {
	double maxDeviation = weights.at(0) * (x.at(0) - zref.at(0));
	for (unsigned int i = 1; i < zref.size(); ++i) {
		double deviation = weights.at(i) * (x.at(i) - zref.at(i));
		if (deviation > maxDeviation)
			maxDeviation = deviation;
	}

	return maxDeviation;
}

vector<double> Momet::nadirPoint(vector<vector<double> >& PFknown) {
	vector<double> nadirEstimate(PFknown.front());
	int PFKnownSize = PFknown.size();
	int nObjs = nadirEstimate.size();
	for (int i = 1; i < PFKnownSize; ++i)
		for (int k = 0; k < nObjs; ++k)
			if (nadirEstimate.at(k) < PFknown.at(i).at(k))
				nadirEstimate.at(k) = PFknown.at(i).at(k);

	return nadirEstimate;
}

double Momet::errorRatio(vector<vector<double> > PFknown, vector<vector<double> > PFtrue) {
	PFknown = removeDominated(PFknown);
	
	double errors = 0.0;
	bool dominated = false;
	int PFKnownSize = PFknown.size();
	int PFTrueSize = PFtrue.size();

	for (int i = 0; i < PFKnownSize; ++i) {
		dominated = false;
		for (int j = 0; j < PFTrueSize && !dominated; ++j)
			if (dominates(PFtrue.at(j), PFknown.at(i)))
				dominated = true;

		if (dominated)
			errors++;
	}

	return errors / PFKnownSize;
}

double Momet::gDistanceSphere(vector<vector<double> >& PFknown) {
	int PFSize = PFknown.size();
	int numObjs = PFknown[0].size();
	double sumTotal = 0.0;
	for (unsigned i = 0; i < PFknown.size(); ++i) {
		double sumVector = 0.0;
		for (int k = 0; k < numObjs; ++k)
			sumVector += PFknown[i][k] * PFknown[i][k];

		sumTotal += sumVector;
	}

	return (sumTotal / ((double) PFSize)) - 1.0;
}

double Momet::gDistancePlane(vector<vector<double> >& PFknown) {
	int PFSize = PFknown.size();
	int numObjs = PFknown[0].size();
	double sumTotal = 0.0;
	for (unsigned i = 0; i < PFknown.size(); ++i) {
		double sumVector = 0.0;
		for (int k = 0; k < numObjs; ++k)
			sumVector += PFknown[i][k];

		sumTotal += sumVector;
	}

	return (sumTotal / ((double) PFSize)) - 1.0;
}

double Momet::gValueDTLZ(vector<vector<double> >& ParetoSetApprox, int K_distVars) {
	int nVars = ParetoSetApprox[0].size();
	int PSSize = ParetoSetApprox.size();
	double gValue, gValueTotal;

	// K_distVars: num. of distance-related variables (Huband's notation).
	// L_posVars : num. of position-related variables (Huband's notation).
	// nVars = K_distVars + L_posVars : total num. of variables.
	int L_posVars = nVars - K_distVars;

	gValueTotal = 0.0;
	for (int i = 0; i < PSSize; ++i) {
		/* Compute g */
		gValue = 0.0;
		for (int j = L_posVars; j < nVars; ++j)
			gValue += ParetoSetApprox[i][j];

		gValue = 1.0 + (9.0 / K_distVars) * gValue;

		gValueTotal += gValue;
	}

	return gValueTotal / ((double) PSSize);
}

double Momet::spacing(vector<vector<double> > PFKnown) {
	PFKnown = removeDominated(PFKnown);
	
	double mean = 0.0;
	double sum = 0.0;
	int PFKnownSize = PFKnown.size();

	// Compute the mean of the distance between neighboring vectors.
	for (int i = 0; i < PFKnownSize; ++i)
		mean += nearestNeighborDistance(PFKnown, i);

	mean /= PFKnownSize;

	// Compute the standard deviation of interdistances.
	for (int i = 0; i < PFKnownSize; ++i)
		sum += pow(mean - nearestNeighborDistance(PFKnown, i), 2);

	return sqrt(sum / (PFKnownSize - 1));
}

double Momet::coverage(vector<vector<double> > a, vector<vector<double> > b) {
	a = removeDominated(a);
	b = removeDominated(b);
	
	double ndominatedOfB = 0;
	for (unsigned int i = 0; i < b.size(); i++) {
		bool isDominated = false;
		for (unsigned int j = 0; j < a.size() && !isDominated; j++) {
			if (weaklyDominates(a[j], b[i]))
				isDominated = true;
		}

		if (isDominated)
			ndominatedOfB++;
	}

	return ndominatedOfB / (double) b.size();
}

double Momet::addEpsilonIndicator(vector<vector<double> > a, vector<vector<double> > b) {
	a = removeDominated(a);
	b = removeDominated(b);
	
	int nObjs = a[0].size();
	double distance, maxDist, minDist;

	// Create an |a|x|b| matrix with the distances from points in A (rows) to points in B (cols).
	vector < vector<double> > distMatrix(a.size());
	for (unsigned int i = 0; i < a.size(); ++i)
		distMatrix[i].resize(b.size());

	// Step 1: Compute the distance matrix epsilon_{z1,z2}.
	for (unsigned int i = 0; i < a.size(); ++i) {
		for (unsigned int j = 0; j < b.size(); ++j) {
			// Find the maximum difference among objectives for all points i in A & j in B.
			maxDist = a[i][0] - b[j][0];
			for (int k = 1; k < nObjs; ++k) {
				distance = a[i][k] - b[j][k];
				if (maxDist < distance)
					maxDist = distance;
			}

			distMatrix[i][j] = maxDist;
		}
	}

	//Step 2: For each point z2 in B compute its distance to the nearest point z1 in A.
	vector<double> epsilonZ2(b.size());
	for (unsigned int i = 0; i < b.size(); ++i) {

		// Find the distance from point i in B to the nearest point z1 in A.
		minDist = distMatrix[0][i];
		for (unsigned int j = 1; j < a.size(); ++j)
			if (minDist > distMatrix[j][i])
				minDist = distMatrix[j][i];

		epsilonZ2[i] = minDist;
	}

	//Step 3: Find the point z2 in B which is farthest to its nearest point z1 in A.
	maxDist = epsilonZ2[0];
	for (unsigned int i = 1; i < b.size(); ++i)
		if (maxDist < epsilonZ2[i])
			maxDist = epsilonZ2[i];

	return maxDist;
}

double Momet::multEpsilonIndicator(vector<vector<double> > a, vector<vector<double> > b) {
	a = removeDominated(a);
	b = removeDominated(b);
	
	int nObjs = a[0].size();
	double distance, maxDist, minDist;
	double delta = 0;
	
	for (unsigned int i = 0; i < a.size(); i++) // To ensure we have positive values
		for (int j = 0; j < nObjs; j++)
			if (a[i][j] < EPS)
				delta = max(delta, fabs(a[i][j]) + EPS);
	for (unsigned int i = 0; i < b.size(); i++) // To ensure we have positive values
		for (int j = 0; j < nObjs; j++)
			if (b[i][j] < EPS)
				delta = max(delta, fabs(b[i][j]) + EPS);

	// Create an |a|x|b| matrix with the distances from points in A (rows) to points in B (cols).
	vector < vector<double> > distMatrix(a.size());
	for (unsigned int i = 0; i < a.size(); ++i)
		distMatrix[i].resize(b.size());

	// Step 1: Compute the distance matrix epsilon_{z1,z2}.
	for (unsigned int i = 0; i < a.size(); ++i) {
		for (unsigned int j = 0; j < b.size(); ++j) {
			// Find the maximum difference among objectives for all points i in A & j in B.
			maxDist = (a[i][0] + delta) / (b[j][0] + delta);
			for (int k = 1; k < nObjs; ++k) {
				distance = (a[i][k] + delta) / (b[j][k] + delta);
				if (maxDist < distance)
					maxDist = distance;
			}

			distMatrix[i][j] = maxDist;
		}
	}

	// Step 2: For each point z2 in B compute its distance to the nearest point z1 in A.
	vector<double> epsilonZ2(b.size());
	for (unsigned int i = 0; i < b.size(); ++i) {
		// Find the distance from point i in B to the nearest point z1 in A.
		minDist = distMatrix[0][i];
		for (unsigned int j = 1; j < a.size(); ++j)
			if (minDist > distMatrix[j][i])
				minDist = distMatrix[j][i];

		epsilonZ2[i] = minDist;
	}

	// Step 3: Find the point z2 in B which is farthest to its nearest point z1 in A.
	maxDist = epsilonZ2[0];
	for (unsigned int i = 1; i < b.size(); ++i)
		if (maxDist < epsilonZ2[i])
			maxDist = epsilonZ2[i];

	return maxDist;
}

double Momet::nearestNeighborDistance(vector<vector<double> >& PFKnown, int i) {
	double minDistance = INF;
	double sum = 0.0;
	int PFKnownSize = PFKnown.size();
	int nObjs = PFKnown.at(i).size();

	for (int j = 0; j < PFKnownSize; ++j) {
		if (j != i) {
			sum = 0.0;
			for (int k = 0; k < nObjs; ++k)
				sum += fabs(PFKnown[i][k] - PFKnown[j][k]);

			if (minDistance > sum)
				minDistance = sum;
		}
	}

	return minDistance;
}

double Momet::euclideanDistance(vector<double>& x, vector<double>& y) {
	double sum = 0.0;
	int nObjs = x.size();

	for (int i = 0; i < nObjs; ++i)
		sum += pow(x[i] - y[i], 2);

	return sqrt(sum);
}

bool Momet::dominates(vector<double>& x, vector<double>& y) {
	unsigned counter = 0;
	unsigned dimension = x.size();

	for (unsigned i = 0; i < dimension; ++i) {
		if (x[i] > y[i])
			return false;
		else if (x[i] == y[i])
			++counter;
	}

	// If the vectors are equal
	if (counter == dimension)
		return false;
	else
		return true;
}

bool Momet::weaklyDominates(vector<double>& x, vector<double>& y) {
	unsigned int dimension = x.size();

	for (unsigned int i = 0; i < dimension; ++i)
		if (x[i] > y[i])
			return false;

	return true;
}

vector<vector<double> > Momet::removeDominated(vector<vector<double> > x) {
	vector<vector<double> > nd;
	for (unsigned int i = 0; i < x.size(); i++) {
		bool isNd = true;
		for (unsigned int j = i+1; j < x.size(); j++)
			if (weaklyDominates(x[j], x[i])) {
				isNd = false;
				break;
			}
		
		if (isNd)
			nd.push_back(x[i]);
	}
		
	return nd;
}

double Momet::nearestDistance(vector<double>& ind, vector<vector<double> >& PFKnown) {
	double minDistance = euclideanDistance(ind, PFKnown[0]);
	for (unsigned int j = 1; j < PFKnown.size(); ++j) {
		double distance = euclideanDistance(ind, PFKnown[j]);
		if (distance < minDistance)
			minDistance = distance;
	}

	return minDistance;
}

double Momet::genDistance(vector<vector<double> > PFknown, vector<vector<double> > PFtrue) {
	PFknown = removeDominated(PFknown);
	
	double sum = 0.0, p = 2;
	int PFKnownSize = PFknown.size();

	for (int i = 0; i < PFKnownSize; ++i) {
		double minDistance = nearestDistance(PFknown[i], PFtrue);
		sum += pow(minDistance, p);
	}

	return pow(sum, 1.0 / p) / PFKnownSize;
}

double Momet::genDistanceP(vector<vector<double> > PFknown, vector<vector<double> > PFtrue) {
	PFknown = removeDominated(PFknown);
	
	double sum = 0.0, p = 2;
	int PFKnownSize = PFknown.size();

	for (int i = 0; i < PFKnownSize; ++i) {
		double minDistance = nearestDistance(PFknown[i], PFtrue);
		sum += pow(minDistance, p);
	}
	
	return pow(sum / PFKnownSize, 1.0 / p);
}

double Momet::invertedGenDistance(vector<vector<double> > PFknown, vector<vector<double> > PFtrue) {
	PFknown = removeDominated(PFknown);
	
	double sum = 0.0, p = 2;
	int PFTrueSize = PFtrue.size();

	for (int i = 0; i < PFTrueSize; ++i) {
		double minDistance = nearestDistance(PFtrue[i], PFknown);
		sum += pow(minDistance, p);
	}

	return pow(sum, 1.0 / p) / PFTrueSize;
}

double Momet::invertedGenDistanceP(vector<vector<double> > PFknown, vector<vector<double> > PFtrue) {
	PFknown = removeDominated(PFknown);
	
	double sum = 0.0, p = 2;
	int PFTrueSize = PFtrue.size();

	for (int i = 0; i < PFTrueSize; ++i) {
		double minDistance = nearestDistance(PFtrue[i], PFknown);
		sum += pow(minDistance, p);
	}

	return pow(sum / PFTrueSize, 1.0 / p);
}

double Momet::deltaP(vector<vector<double> > PFknown, vector<vector<double> > PFtrue) {
	return max(genDistanceP(PFknown, PFtrue), invertedGenDistanceP(PFknown, PFtrue));
}
