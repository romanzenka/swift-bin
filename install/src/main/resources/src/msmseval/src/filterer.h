/***************************************************************************************************
File:		filterer.cpp
Created:	2006/14/05
Author:		Jason W H Wong <jason.wong@chem.ox.ac.uk>
	
Purpose:	Implements all the methods to evaluate tandem mass spectra

License:	The stand-alone msmsEval application and its source code is free to academic, 
		government and non-profit users for non-commercial use. Use of the msmsEval
		application and source code in part or in whole for commercial use by non-profit 
		and for-profit entities requires a license agreement. Parties interested in 
		commercial should contact jason.wong@chem.ox.ac.uk for permission to use
		the msmsEval application or its source code.

		The use of msmsEval in any research must cite as appropriate:
		Wong, J.W.H., Sullivan, M.J., Cartwright, H.M. and Cagney, G. (2007) Empirical modelling
		of tandem mass spectra quality for high-throughput proteomics. BMC Bioinformatics. 8:51.

Disclaimer:	msmsEval and all associated documents ARE PROVIDED "AS IS" WITHOUT ANY WARRANTY 
		OF ANY KIND, EITHER EXPRESS, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, 
		ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
		FREEDOM FROM INFRINGEMENT, OR THAT msmsEval and all associated documents WILL 
		BE ERROR FREE. The authur MAKES NO REPRESENTATIONS THAT THE USE OF msmsEval or 
		any associated documents WILL NOT INFRINGE ANY PATENT OR PROPRIETARY RIGHTS OF 
		THIRD PARTIES. In no event will the authors be liable for any damages, including 
		but not limited to direct, indirect, special or consequential damages, arising out 
		of, resulting from, or in any way connected with the use of msmsEval or any associated 
		documents.
***************************************************************************************************/
#ifndef FILTERER_H
#define FILTERER_H

#include <vector>
#include <string>
#include "fileHandle.h"
#include <fstream>

using namespace std;


class GreaterThanIntn{
public:
	bool operator()(const peakData &a, const peakData &b)
	{return a.intn> b.intn;}
};

class GreaterThanIntnS{
public:
	bool operator()(const peakData *a, const peakData *b)
	{return a->intn> b->intn;}
};

class LessThanMZ{
public:
	bool operator()(const peakData &a, const peakData &b)
	{return a.exact_mz< b.exact_mz;}
};


class GreaterThanInt{
public:
	bool operator()(const int &a, const int &b)
	{return a> b;}
};

class filterer
{

public:

	filterer();
	virtual ~filterer();

	void analyze(vector<scan*> &sec, parameters &params);
	void run_EM(vector <scan*> &sec, const char* path, parameters &params, bool quiet);
	void sortGB(vector <scan*> &sec);
	void filterData(vector <scan*> &sec, double cutoff, bool pF);

	double setCutoff(double p);

private:

	int goodSegments(scan *s1);
	void neutralPeaks(scan *s1);
	void evalCharge(scan *s1);
	void calculateParents(scan *s1);
	bool AAdiff(double mz1, double mz2);
	int commonDiff(scan *s1);
    void complement(scan *s1);
	void calNormTIC(vector<scan*> &sec);
	void calIntnRatio(scan *s1);

	inline int round_int(double x);

	double discrim_m(scan &s, parameters &p);
	
	inline double cal_pF(double Fp,double Fm,double p, double m);
	inline double cal_mF(double Fp,double Fm,double p, double m);
	inline double cal_Fp(double f, double avg, double stdev);
	inline double cal_Fm(double m, double avg, double stdev);
	double cal_p(vector <scan*> &sec);
	double cal_pavg(double p, vector <scan*> &sec);
	double cal_pstdev(double p, double avg, vector <scan*> &sec);
	double cal_mavg(double m, vector <scan*> &sec);
	double cal_mstdev(double m, double avg, vector <scan*> &sec);
	double c_norm(double x);

	double cutoff;
	
	parameters em_params;
	double final_p_avg;
	double final_p_std;
};
#endif
