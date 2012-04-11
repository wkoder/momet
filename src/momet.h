/*
 * momet.h
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#ifndef MOMET_H_
#define MOMET_H_

#include <stdlib.h>
#include <vector>

using namespace std;

class Momet {
public:
	Momet();
	virtual ~Momet();


	double errorRatio(vector<vector<double> > PFknown, vector<vector<double> > PFtrue);
	double spacing(vector<vector<double> > PFKnown);
	double coverage(vector<vector<double> > a, vector<vector<double> > b);
	double addEpsilonIndicator(vector<vector<double> > a, vector<vector<double> > b);
	double multEpsilonIndicator(vector<vector<double> > a, vector<vector<double> > b);
	double hypervolume(vector<vector<double> > PFKnown, vector<double> reference);

	double genDistance(vector<vector<double> > PFknown, vector<vector<double> > PFtrue);
	double invertedGenDistance(vector<vector<double> > PFknown, vector<vector<double> > PFtrue);
	double genDistanceP(vector<vector<double> > PFknown, vector<vector<double> > PFtrue);
	double invertedGenDistanceP(vector<vector<double> > PFknown, vector<vector<double> > PFtrue);
	double deltaP(vector<vector<double> > PFknown, vector<vector<double> > PFtrue);

	vector<double> nadirPoint(vector<vector<double> >& PFknown);
	double bestAchievement(vector<vector<double> >& PFknown, int &pos, vector<double> zref, vector<double> weights);
	double gDistanceSphere(vector<vector<double> >& PFknown);
	double gDistancePlane(vector<vector<double> >& PFknown);
	double gValueDTLZ(vector<vector<double> >& ParetoSetApprox, int K_distVars);

private:
	bool dominates(vector<double>& x, vector<double>& y);
	bool weaklyDominates(vector<double>& x, vector<double>& y);
	double nearestNeighborDistance(vector<vector<double> >& PFKnown, int i);
	double nearestDistance(vector<double>& ind, vector<vector<double> >& PFKnown);
	double euclideanDistance(vector<double>& x, vector<double>& y);
	double chebyAchievement(vector<double> x, vector<double> zref, vector<double> weights);
	vector<vector<double> > removeDominated(vector<vector<double> > x);
};

#endif /* MOMET_H_ */
