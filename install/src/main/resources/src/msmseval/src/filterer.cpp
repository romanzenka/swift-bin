/***************************************************************************************************
File:		filterer.cpp
Created:	2006/19/12
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

Disclaimer:	msmsEval and all associated documents ARE PROVIDED "AS IS" WITH ANY WARRANTY 
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

#include "filterer.h"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>

#include <errno.h>

#define PI	 3.1415926535897932
#define A 	 71.03711
#define R 	 156.10111
#define N 	 114.04293
#define D 	 115.02694
#define C 	 160.00919
#define E 	 129.04259
#define Q 	 128.05858
#define G 	 57.02146
#define H 	 137.05891
#define IL 	 113.08406
#define K 	 128.09496
#define M 	 131.04049
#define F 	 147.06841
#define P 	 97.05276
#define S 	 87.03203
#define T 	 101.04768
#define W 	 186.07931
#define Y 	 163.06333
#define V 	 99.06841


filterer::filterer(){
	//tout.open("dump",ios::in|ios::trunc);
}

filterer::~filterer(){
	//tout.close();
}

double filterer::setCutoff(double p)
{
	cutoff = p;
	return cutoff;
}

void filterer::analyze(vector<scan*> &sec, parameters &params)
{
	for (unsigned int i = 0; i < sec.size();++i){
		sec.at(i)->TICFore = 0;
		sec.at(i)->TICBack = 0;
		sec.at(i)->sumPeaksIntn = 0;
		if (i % 500 == 0){
			cout <<".";
			cout.flush();
		}
		if (sec.at(i)->mz.size() != 0){		
			peakData p;
			for (unsigned int j = 0; j < sec.at(i)->mz.size(); ++j){
				p.has_a = false;
				p.has_comp = false;
				p.has_h2o_nh3 = false;
				p.has_iso = false;
				p.posi = j;
				p.intn = sec.at(i)->intn.at(j);
				p.exact_mz = sec.at(i)->mz.at(j);
				p.mz = round_int(sec.at(i)->mz.at(j));
				sec.at(i)->peaks.push_back(p);
			}
			sec.at(i)->num_peaks = sec.at(i)->peaks.size();
			sec.at(i)->charge[0] = 0;
			sec.at(i)->charge[1] = 0;
			sec.at(i)->TIC = 0;
			float sumP = 0;
			double maxFore = 0;
			double maxBack = 0;
			double maxPeak = 0;
			double secondPeak = 0;
			for (unsigned int k = 0; k < sec.at(i)->peaks.size(); ++k){
				sec.at(i)->TIC += sec.at(i)->peaks.at(k).intn;
				sumP += sec.at(i)->peaks.at(k).intn;
				if (sec.at(i)->peaks.at(k).mz <= sec.at(i)->preMZ){
					sec.at(i)->TICFore+=sec.at(i)->peaks.at(k).intn;
					if (sec.at(i)->peaks.at(k).intn > maxFore)
						maxFore = sec.at(i)->peaks.at(k).intn;
				}
				else{
					sec.at(i)->TICBack+=sec.at(i)->peaks.at(k).intn;
					if (sec.at(i)->peaks.at(k).intn > maxBack)
						maxBack = sec.at(i)->peaks.at(k).intn;
				}
				if (sec.at(i)->peaks.at(k).intn>maxPeak){
					secondPeak = maxPeak;
					maxPeak=sec.at(i)->peaks.at(k).intn;
				}
				else if (sec.at(i)->peaks.at(k).intn>secondPeak && sec.at(i)->peaks.at(k).intn<maxPeak){
					secondPeak = sec.at(i)->peaks.at(k).intn;
				}
			}
			sec.at(i)->TICFore -= maxFore;
			sec.at(i)->TICBack -= maxBack;
			sec.at(i)->sumPeaksIntn=sumP;
			
			if (sec.at(i)->TICFore*0.125 > sec.at(i)->TICBack){
				sec.at(i)->charge[0]=1;
				sec.at(i)->charge[1]=0;
				for (unsigned int j = 0; j < sec.at(i)->peaks.size(); ++j){
					if (sec.at(i)->peaks.at(j).mz > sec.at(i)->preMZ){
						sec.at(i)->peaks.erase(sec.at(i)->peaks.begin()+j);
						j--;
					}
				}
			}
			sec.at(i)->num_peaks = sec.at(i)->peaks.size();
			sec.at(i)->basePeakIntensity= maxPeak;
			sec.at(i)->secBasePeakIntn = secondPeak;
			
			calIntnRatio(sec.at(i));
			goodSegments(sec.at(i));
			neutralPeaks(sec.at(i));
			calculateParents(sec.at(i));
			complement(sec.at(i));
			int num_r_peaks = 0;
			for (unsigned int j = 0; j < sec.at(i)->peaks.size(); ++j){
				if (sec.at(i)->peaks.at(j).rank<4)
					num_r_peaks++;
			}
			sec.at(i)->num_r_peaks = num_r_peaks;
			sec.at(i)->iso_ratio = 0;
			if (num_r_peaks != 0)
				sec.at(i)->iso_ratio = (double)sec.at(i)->num_iso/(double)sec.at(i)->num_r_peaks;
			if (num_r_peaks != 0)
				sec.at(i)->h2o_ratio = (double)sec.at(i)->num_h2o/(double)sec.at(i)->num_r_peaks;
		}
//tout << endl;
	}
	calNormTIC(sec);	
	for (unsigned int i = 0; i < sec.size();++i){		
		sec.at(i)->discrim=discrim_m(*sec.at(i),params);
	
//		if (sec.at(i)->discrim>cutoff)
//			sec.at(i)->good = true;
//		else
//			sec.at(i)->good = false;
//		evalCharge(sec.at(i));
	}
	cout <<"Done"<<endl;
}

bool filterer::AAdiff(double mz1, double mz2)
{
	//double diff = fabs(mz1-mz2);
	bool matched = false;
	if ((fabs(mz1-mz2) >= A-0.3 && fabs(mz1-mz2) <= A+0.3) ||
		(fabs(mz1-mz2) >= R-0.3 && fabs(mz1-mz2) <= R+0.3) ||
		(fabs(mz1-mz2) >= N-0.3 && fabs(mz1-mz2) <= N+0.3) ||
		(fabs(mz1-mz2) >= D-0.3 && fabs(mz1-mz2) <= D+0.3) ||
		(fabs(mz1-mz2) >= C-0.3 && fabs(mz1-mz2) <= C+0.3) ||
		(fabs(mz1-mz2) >= E-0.3 && fabs(mz1-mz2) <= E+0.3) ||
		(fabs(mz1-mz2) >= Q-0.3 && fabs(mz1-mz2) <= Q+0.3) ||
		(fabs(mz1-mz2) >= G-0.3 && fabs(mz1-mz2) <= G+0.3) ||
		(fabs(mz1-mz2) >= H-0.3 && fabs(mz1-mz2) <= H+0.3) ||
		(fabs(mz1-mz2) >= IL-0.3 && fabs(mz1-mz2) <= IL+0.3) ||
		(fabs(mz1-mz2) >= K-0.3 && fabs(mz1-mz2) <= K+0.3) ||
		(fabs(mz1-mz2) >= M-0.3 && fabs(mz1-mz2) <= M+0.3) ||
		(fabs(mz1-mz2) >= F-0.3 && fabs(mz1-mz2) <= F+0.3) ||
		(fabs(mz1-mz2) >= P-0.3 && fabs(mz1-mz2) <= P+0.3) ||
		(fabs(mz1-mz2) >= S-0.3 && fabs(mz1-mz2) <= S+0.3) ||
		(fabs(mz1-mz2) >= T-0.3 && fabs(mz1-mz2) <= T+0.3) ||
		(fabs(mz1-mz2) >= W-0.3 && fabs(mz1-mz2) <= W+0.3) ||
		(fabs(mz1-mz2) >= Y-0.3 && fabs(mz1-mz2) <= Y+0.3)){
			matched = true;
	}
	return matched;
}

int filterer::goodSegments(scan *s1)
{
  int goodSegs = 0;
  s1->SNR = 0;
  int segs = 0;
  double endmass = s1->mz.back();
  if (s1->charge[0] ==1)
	endmass = s1->preMZ+57;
  if (!s1->mz.empty()){
    unsigned int index = 0;
    float spos = 0;
    while (index < s1->mz.size() && s1->mz.at(index) < 200){
      spos = s1->mz.at(index);
      index++;
    }
    vector <peakData *> tmpData;

    	while (index < s1->peaks.size() && spos < endmass){
		while (index < s1->peaks.size() && s1->peaks.at(index).mz < spos +57){
			tmpData.push_back(&s1->peaks.at(index));
			index++;
		}
		spos+=57;
		segs++;
		if (!tmpData.empty()){
			sort(tmpData.begin(),tmpData.end(),GreaterThanIntnS());
			int rank = 1;
			tmpData.at(0)->rank =rank;
			for (unsigned int j = 1; j < tmpData.size(); ++j){
				bool is_iso = false;
				for (unsigned int k = 0; k < j; ++k){
					if (tmpData.at(j)->exact_mz < tmpData.at(k)->exact_mz+2.4 && tmpData.at(j)->exact_mz > tmpData.at(k)->exact_mz)
						is_iso = true;
				}
				if (is_iso){
					tmpData.at(j)->rank =999;
					tmpData.erase(tmpData.begin()+j);
					j--;
				}
				else{
					tmpData.at(j)->rank =rank;
					rank++;
				}
			}
			if (tmpData.size() > 5){
				tmpData.resize(5);
			}
			if (tmpData.size() >= 3  && tmpData.front()->intn > tmpData.back()->intn*3){
				goodSegs++;
			}
			else if ((tmpData.size() ==2|| tmpData.size() ==1)&& tmpData.front()->intn > s1->basePeakIntensity*0.4)
				goodSegs++;
      		}
      		tmpData.clear();
    	}
  }

  if (segs <= 5 || endmass < 500)
	s1->goodSegs = 0;
  else if (segs <= 10)
	s1->goodSegs = ((double)goodSegs/((double)segs*2));
  else
  	s1->goodSegs = ((double)goodSegs/(double)segs);
  return goodSegs;
}

void filterer::neutralPeaks(scan *s1)
{
  s1->num_a = 0;
  s1->num_h2o = 0;
  s1->num_iso = 0;
  s1->AA_diff_ratio = 0;
  sort(s1->peaks.begin(),s1->peaks.end(),GreaterThanIntn());
  for (unsigned int i = 0; i < s1->peaks.size(); ++i){
    if (s1->peaks.at(i).rank < 4){
      for (unsigned int j = 0; j < s1->peaks.size() && s1->peaks.at(i).intn != 0; ++j){
	if (s1->peaks.at(j).exact_mz > s1->peaks.at(i).exact_mz && s1->peaks.at(j).exact_mz <= s1->peaks.at(i).exact_mz+1.4){
	  s1->peaks.at(i).has_iso = true;
	  s1->num_iso++;
	}
	if (s1->peaks.at(j).exact_mz > s1->peaks.at(i).exact_mz-18.5 && s1->peaks.at(j).exact_mz <= s1->peaks.at(i).exact_mz-17.5){
	  s1->peaks.at(i).has_h2o_nh3 = true;
	  s1->num_h2o++;
	}
      }
    }
  }

  bool found = false;
  double  sumRank = 0;
  int hits = 0;
  int sumR1 = 0;
  // Calculating differences
  for (unsigned int k = 0; k < s1->peaks.size(); ++k){
      if (s1->peaks.at(k).rank <=1){
	sumR1++;
      for (unsigned int j = k+1; j < s1->peaks.size()&& !found; ++j){
	if (s1->peaks.at(j).rank <=2){
	  bool diffAA = AAdiff(s1->peaks.at(j).exact_mz,s1->peaks.at(k).exact_mz);
	  if (diffAA){
		hits++;
		if (hits < 3)
			sumRank++;
		else
			sumRank--;
	  }
	}
	s1->AA_diff_ratio += sumRank;
	hits = 0;
	sumRank = 0;
      }
     }
      found = false;
  }
	if (sumR1 <= 2)
	    s1->AA_diff_ratio = 0;
	else	
	    s1->AA_diff_ratio/=(double)sumR1;
	if (s1->AA_diff_ratio < 0 )
	    s1->AA_diff_ratio = 0;
		
  sort(s1->peaks.begin(),s1->peaks.end(),LessThanMZ());
}

void filterer::complement(scan *s1)
{
  s1->match2 = 0;
  s1->match3 = 0;
	s1->normMatch1 =0;
	s1->normMatch2 =0;
	s1->normMatch3 =0;
	for (unsigned int j = 0; j < s1->peaks.size(); ++j){
		if (s1->peaks.at(j).rank < 10){
			for (unsigned int k = 0; k < s1->peaks.size(); ++k){
				if (s1->parent1-s1->peaks.at(j).exact_mz > s1->peaks.at(k).exact_mz-1.4 && s1->parent1-s1->peaks.at(j).exact_mz < s1->peaks.at(k).exact_mz+1.4){
					s1->normMatch1 += 1.0/max(s1->peaks.at(k).rank,s1->peaks.at(j).rank);
				}
				if (s1->parent2-s1->peaks.at(j).exact_mz > s1->peaks.at(k).exact_mz-1.4 && s1->parent2-s1->peaks.at(j).exact_mz < s1->peaks.at(k).exact_mz+1.4){
					s1->normMatch2 += 1.0/max(s1->peaks.at(k).rank,s1->peaks.at(j).rank);
					s1->match2++;
				}
				if (s1->parent3-s1->peaks.at(j).exact_mz*2 > s1->peaks.at(k).exact_mz-1.4 && s1->parent3-s1->peaks.at(j).exact_mz*2 < s1->peaks.at(k).exact_mz+1.4){
					s1->normMatch3 += 1.0/max(s1->peaks.at(k).rank,s1->peaks.at(j).rank);
					s1->match3++;
				}
			}
		}
	}
    if (s1->TICFore*0.125 > s1->TICBack){
	s1->comple = s1->normMatch1;
    }
    else
	s1->comple = max(s1->normMatch2,s1->normMatch3);
}

void filterer::calNormTIC(vector<scan*> &sec)
{
	double meanTIC = 0;
	for (unsigned int i = 0;i < sec.size(); ++i){
		meanTIC+=sec.at(i)->TIC;
	}
	meanTIC/=(double)sec.size();
	for (unsigned int i = 0;i < sec.size(); ++i){
		sec.at(i)->normTIC = sec.at(i)->TIC/meanTIC;
	}
}

void filterer::calIntnRatio(scan *s1)
{
	if (s1->intn.size() < 5){
		s1->intnRatio1 = 1;
		s1->intnRatio20 = 1;
		return;
	}
	double peaks_over_1 = 0;
	double peaks_over_20 = 0;
	for (unsigned int i = 0; i < s1->intn.size(); ++i){
		double cur_intn = s1->intn.at(i)/s1->basePeakIntensity;
		if (cur_intn > 0.01)
			peaks_over_1++;
		if (cur_intn > 0.2)
			peaks_over_20++;
	}
	s1->intnRatio1 = peaks_over_1/(double)s1->intn.size();
	s1->intnRatio20 =  peaks_over_20/(double)s1->intn.size();
}

void filterer::evalCharge(scan *s1)
{
//  if (!s1->good){
//    s1->charge[0] = 0;
//    s1->charge[1] = -1;
//  }
//  else{
    if (s1->TICFore*0.125 > s1->TICBack){
      s1->charge[0]=1;
      s1->charge[1]=0;
    }
    else{
    	if (s1->match2 > 2 || s1->match3 > 2){
    		if (s1->preMZ < 700 && s1->normMatch3 >= s1->normMatch2*1.5){
	  			s1->charge[0]=3;
	  			s1->charge[1]=(int)s1->match3;
			}
			if (s1->preMZ < 700 && s1->normMatch2 < s1->normMatch3*4){
	  			s1->charge[0]=4;
	  			s1->charge[1]=(int)s1->match3;
			}
			else if (s1->normMatch2 >= s1->normMatch3*2){
	  			s1->charge[0]=2;
	  			s1->charge[1]=(int)s1->match2;
			}
			else if (s1->normMatch3 >= s1->normMatch2*2){
	  			s1->charge[0]=3;
	  			s1->charge[1]=(int)s1->match3;
			}
			else {
	  			s1->charge[0]=4;
	  			s1->charge[1]=(int)max(s1->match2,s1->match3);
			}
		}
		else{
			s1->charge[0]=4;
			s1->charge[1]=(int)max(s1->match2,s1->match3);
      	}
    }
  //}
}

void filterer::calculateParents(scan *s1)
{
	int pCounts[5000];
	double pSum[5000];
	for (unsigned int k = 0; k < 5000; ++k){
		pCounts[k] = 0;
		pSum[k] = 0;
	}

	for (unsigned int k = 0; k <s1->peaks.size(); ++k){
		if (s1->peaks.at(k).rank < 10){
			for (unsigned int l = 0; l <s1->peaks.size(); ++l){
				double sumval = s1->peaks.at(k).exact_mz+s1->peaks.at(l).exact_mz;
				double sum2val = s1->peaks.at(k).exact_mz*2+s1->peaks.at(l).exact_mz;
				if ((sumval > s1->preMZ-3 &&sumval < s1->preMZ+3) || (sumval > s1->preMZ*2-1-3 && sumval < s1->preMZ*2-1+3)){
					int roundVal = round_int(sumval);
					if (roundVal < 5000){
					pCounts[roundVal]++;
					pSum[roundVal]+=sumval;
					}
				}
				if ((sumval > s1->preMZ*3-1-3 &&sumval < s1->preMZ*3-1+3)){
					int roundVal = round_int(sum2val);
					if (roundVal < 5000){
					pCounts[roundVal]++;
					pSum[roundVal]+=sum2val;
					}
				}
			}
		}
	}
  int maxP1 = 0;
  int maxP2 = 0;
  int maxP3 = 0;
  int pos1 = 0;
  int pos2 = 0;
  int pos3 = 0;
  unsigned int pre = round_int(s1->preMZ);
  for (unsigned int k =pre-3; k <= pre+3; ++k){
    if (k > 0 && k < 5000 && pCounts[k] > maxP1){
      maxP1 = pCounts[k];
      pos1 = k;
    }
  }
  pre = round_int(s1->preMZ*2-1);
  for (unsigned int k =pre-3; k <= pre+3; ++k){
    if (k > 0 && k < 5000 && pCounts[k] > maxP2){
      maxP2 = pCounts[k];
      pos2 = k;
    }
  }
  pre = round_int(s1->preMZ*3-1);
  for (unsigned int k =pre-3; k <= pre+3; ++k){
    if (k > 0 && k < 5000 && pCounts[k] > maxP3){
      maxP3 = pCounts[k];
      pos3 = k;
    }
  }
  if (maxP1 == 0)
    s1->parent1 = s1->preMZ;					// stored as parent+H
  else
    s1->parent1 = pSum[pos1]/pCounts[pos1];
  if (maxP2 == 0)
    s1->parent2 = s1->preMZ*2-1;					// stored as parent+H
  else
    s1->parent2 = pSum[pos2]/pCounts[pos2];
  if (maxP3 == 0)
    s1->parent3 = s1->preMZ*3-1;					// stored as parent+H
  else
    s1->parent3 = pSum[pos3]/pCounts[pos3];
}

int filterer::commonDiff(scan *s1)
{
	int common = 0;
	if (s1->charge[0] == 1){
		vector <int> dists;
		dists.assign(300,0);
		vector <peakData *> r;
		for (unsigned int i = 0; i < s1->peaks.size(); ++i){
			if (s1->peaks.at(i).rank == 1 || s1->peaks.at(i).rank == 2){
				r.push_back(&s1->peaks.at(i));
			}
		}
		for (unsigned int i = 0 ; i < r.size(); ++i){
			for (unsigned int j = i; j < r.size() && j < i+5; ++j){
				double diff = fabs(r.at(i)->exact_mz-r.at(j)->exact_mz);
				if (diff > 20 && diff < 300){
					dists[round_int(diff)]++;
				}
			}
		}
		for (unsigned int i = 0; i < dists.size(); ++i){
			if (dists.at(i) >= 2){
				common += dists.at(i)-1;
				cout << dists.at(i) << " ";
			}
		}
		cout << endl;
	}
	s1->commonDiff = common;
	return common;
}

inline int filterer::round_int (double x)
{
	int i = (int)round(x);
  return (i);
}

double filterer::discrim_m(scan &s, parameters &p)
{
	return p.p_num_peaks*s.num_peaks+p.p_normTIC*s.normTIC+p.p_goodSegs*s.goodSegs+
	       p.p_intnRatio1*s.intnRatio1+p.p_intnRatio20*s.intnRatio20+
	       p.p_comple*s.comple+p.p_iso_ratio*s.iso_ratio+p.p_h2o_ratio*s.h2o_ratio+
	       p.p_AA_diff_ratio*s.AA_diff_ratio+p.p_const;
}

void filterer::run_EM(vector <scan*> &sec,const char* path, parameters & params, bool quiet)
{
	em_params = params;
	cout << "Performing EM optimization";
	double pavg = em_params.em_p_avg;
	double pstdev = em_params.em_p_std;
	double p = em_params.em_p_prior;

	double mavg = em_params.em_m_avg;
	double mstdev = em_params.em_m_std;
	double m = 1-p;

	vector <double> pavg_list;
	vector <double> pstdev_list;
	vector <double> mavg_list;
	vector <double> mstdev_list;
	vector <double> p_list;
	vector <double> lsqr_list;

	unsigned int iters = 0;
	bool negged = false;
	bool stab_neg = false;
	bool stab_pos = false;

	double minma[2];
	double maxma[2];
	while (iters < 500 && (!stab_neg &&!stab_pos)){
		minma[0]=10;
		minma[1]=10;
		maxma[0]=0;
		maxma[1]=0;
		iters++;
		//lsqr_list.push_back(least_sqr(p,pavg,pstdev,mavg,mstdev));
		pavg_list.push_back(pavg);
		pstdev_list.push_back(pstdev);
		mavg_list.push_back(mavg);
		mstdev_list.push_back(mstdev);
		p_list.push_back(p);
		for (unsigned int i = 0; i < sec.size(); ++i){
			sec.at(i)->Fp = cal_Fp(sec.at(i)->discrim,pavg,pstdev);
			sec.at(i)->Fm = cal_Fp(sec.at(i)->discrim,mavg,mstdev);
			sec.at(i)->pF = cal_pF(sec.at(i)->Fp,sec.at(i)->Fm,p,m);
			sec.at(i)->mF = 1-sec.at(i)->pF;
 	 //		if (sec.at(i)->discrim > 5){
  	//				sec.at(i)->pF = 1;
  	//				sec.at(i)->mF = 0;
 	//		}
			if (sec.at(i)->pF < minma[0]){
				minma[0]=sec.at(i)->pF;
				minma[1]=sec.at(i)->discrim;
			}
			if (sec.at(i)->pF > maxma[0]){
				maxma[0]=sec.at(i)->pF;
				maxma[1]=sec.at(i)->discrim;
			}
		}
		for (unsigned int i = 0; i < sec.size(); ++i){
			if (sec.at(i)->discrim < minma[1] && sec.at(i)->pF > minma[0]){
				sec.at(i)->pF = minma[0];
				sec.at(i)->mF = 1-minma[0];
			}
			if (sec.at(i)->discrim > maxma[1] && sec.at(i)->pF < maxma[0]){
				sec.at(i)->pF = maxma[0];
				sec.at(i)->mF = 1-maxma[0];
			}
			if (sec.at(i)->discrim < -12){
				sec.at(i)->pF = 0.0001;
				sec.at(i)->mF = 0.9999;
			}
		}
		if (!negged){
			p = cal_p(sec);
			pavg = cal_pavg(p,sec);
			pstdev = cal_pstdev(p,pavg,sec);
		}
		if (!stab_neg){
			m = 1-p;
			mavg = cal_mavg(m,sec);
			mstdev = cal_mstdev(m,mavg,sec);
		}
		if (pavg < em_params.em_p_avg_llimit && iters > 1){
			negged = true;
			double maxp = 0;
			int mpos = 0;
			for (unsigned int i = 1; i < pavg_list.size();++i){
				if (pavg_list.at(i) > maxp){

					maxp = pavg_list.at(i);
					mpos = i;
				}
			}
			p = p_list.back();
			pavg = pavg_list.back();
			pstdev = pstdev_list.back();
			m = 1-p;
			mavg = mavg_list.at(mpos);
			mstdev = mstdev_list.at(mpos);
		}
		if (mavg > em_params.em_m_avg_ulimit && iters > 4){
			m = 1-p;
			mavg = mavg_list.back();
			mstdev = mstdev_list.back();
		}		
		if (iters > 2 && fabs(mavg-mavg_list.back())< 0.0001)
			stab_neg = true;
		if (iters > 2 && (!negged||stab_neg) && fabs(pavg-pavg_list.back())< 0.0001)
			stab_pos = true;
	}
	//lsqr_list.push_back(least_sqr(p,pavg,pstdev,mavg,mstdev));
	pavg_list.push_back(pavg);
	pstdev_list.push_back(pstdev);
	mavg_list.push_back(mavg);
	mstdev_list.push_back(mstdev);
	p_list.push_back(p);
	cout << "\nConverged in "<<iters << " iterations" << endl;
	if (!quiet){
		ofstream out;
		out.open(path,ios::in|ios::trunc);
		out << "p(+),Average+,Stdev+,p(-),Average-,Stdev-"<<endl;
		for (unsigned int i = 0; i < iters+1; ++i){
			out << p_list.at(i)<<","<<pavg_list.at(i)<<","<<pstdev_list.at(i)<<","<<1-p_list.at(i)<<","<<mavg_list.at(i)<<","<<mstdev_list.at(i)<<endl;
		}
	}
//	for (unsigned int i = 0; i < sec.size(); ++i){
//		if (sec.at(i)->pF > cutoff)
//			sec.at(i)->good = true;
//		else
//			sec.at(i)->good = false;
//		evalCharge(sec.at(i));
//	}
	final_p_avg = pavg;
	final_p_std = pstdev;
}

double filterer::cal_pF(double Fp,double Fm,double p, double m)
{
	double pF = (Fp*p)/((Fp*p)+(Fm*m));
	return pF;
}

double filterer::cal_mF(double Fp,double Fm,double p, double m)
{
	double mF = (Fm*m)/((Fm*m)+(Fp*p));
	return mF;
}

double filterer::cal_Fp(double f, double avg, double stdev)
{
	double Fp = (1/(stdev*sqrt(2*PI)))*exp(-(pow(f-avg,2)/(2*pow(stdev,2))));;
	return Fp;
}

double filterer::cal_Fm(double m, double avg, double stdev)
{
	double beta = stdev*0.77969680123368;
	double mu = avg+(0.577215664901532*beta);	
	double exp_t = exp((m-mu)/beta);
	double Fm = (exp_t*exp(-exp_t))/beta;
	return Fm;
}


double filterer::cal_p(vector <scan*> &sec)
{
	double sumpFF = 0;
	double num_gscan = 0;
	for (unsigned int i = 0; i < sec.size(); ++i){
		if (sec.at(i)->mz.size() >= 5){
			num_gscan++;
			sumpFF+=sec.at(i)->pF;
		}
	}
	if (num_gscan == 0)
		return 0;
	double prior = sumpFF/num_gscan;
	if (prior < em_params.em_prior_llimit)
		prior = em_params.em_prior_llimit;
	if (prior > em_params.em_prior_ulimit)
		prior = em_params.em_prior_ulimit;
	return prior;
}

double filterer::cal_pavg(double p, vector <scan*> &sec)
{
	double sumpFF = 0;
	double num_gscan = 0;
	for (unsigned int i = 0; i < sec.size(); ++i){
		if (sec.at(i)->mz.size() >= 5){
			num_gscan++;
			sumpFF+=(sec.at(i)->discrim*sec.at(i)->pF);
		}
	}
	if (num_gscan == 0 || p == 0)
		return 0;
	return (1/(num_gscan*p))*sumpFF;
}	
double filterer::cal_pstdev(double p, double avg, vector <scan*> &sec)
{
	double sumpFF = 0;
	double num_gscan = 0;
	for (unsigned int i = 0; i < sec.size(); ++i){
		if (sec.at(i)->mz.size() >= 5){
			num_gscan++;
			sumpFF+=(pow(sec.at(i)->discrim-avg,2)*sec.at(i)->pF);
		}
	}
	if (num_gscan == 0 || p == 0)
		return 0;
	double dev = sqrt((1/(num_gscan*p))*sumpFF);
	if (dev < em_params.em_p_std_llimit)
		dev = em_params.em_p_std_llimit;
	if (dev > em_params.em_p_std_ulimit)
		dev = em_params.em_p_std_ulimit;
	return dev;
}

double filterer::cal_mavg(double m, vector <scan*> &sec)
{
	double sumpFF = 0;
	double num_gscan = 0;
	for (unsigned int i = 0; i < sec.size(); ++i){
		if (sec.at(i)->mz.size() >= 5){
			num_gscan++;
			sumpFF+=(sec.at(i)->discrim*(1-sec.at(i)->pF));
		}
	}
	if (num_gscan == 0 || m == 0)
		return 0;
	return (1/(num_gscan*m))*sumpFF;
}

double filterer::cal_mstdev(double m, double avg, vector <scan*> &sec)
{
	double sumpFF = 0;
	double num_gscan = 0;
	for (unsigned int i = 0; i < sec.size(); ++i){
		if (sec.at(i)->mz.size() >= 5){
			num_gscan++;
			sumpFF+=(pow(sec.at(i)->discrim-avg,2)*(1-sec.at(i)->pF));
		}
	}
	if (num_gscan == 0 || m == 0)
		return 0;
	double dev = sqrt((1/(num_gscan*m))*sumpFF);
	if (dev < em_params.em_m_std_llimit)
		dev = em_params.em_m_std_llimit;
	if (dev > em_params.em_m_std_ulimit)
		dev = em_params.em_m_std_ulimit;
	return dev;
}

double filterer::c_norm(double x)
{
	int neg = (x<0);
	if(neg) x *= -1;
	double k(1/(1+0.2316419*x));
	double y=((((1.330274429*k-1.821255978)*k+1.781477937)*k-0.356563782)*k+0.319381530)*k;
	y = 1.0 - 0.398942280401*exp(-0.5*x*x)*y;
	return (1-neg)*y + neg*(1-y);
}

void filterer::sortGB(vector <scan*> &sec)
{
	for (unsigned int i = 0; i < sec.size(); ++i){
		double z_score = (sec.at(i)->discrim-final_p_avg)/final_p_std;
		sec.at(i)->PP = c_norm(z_score);
	//	if (sec.at(i)->PP > cutoff)
	//		sec.at(i)->good = true;
	//	else
	//		sec.at(i)->good = false;
		evalCharge(sec.at(i));
	}
}

void filterer::filterData(vector <scan*> &sec, double cutoff, bool pF)
{
	for (unsigned int i = 0; i < sec.size(); ++i){
		if (pF){
			if (sec.at(i)->pF > cutoff)
				sec.at(i)->good = true;
			else
				sec.at(i)->good = false;
		}
		else{
			if (sec.at(i)->PP > cutoff)
				sec.at(i)->good = true;
			else
				sec.at(i)->good = false;
		}
	}
}
			
