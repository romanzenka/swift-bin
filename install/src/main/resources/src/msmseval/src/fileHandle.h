/***************************************************************************************************
File:		fileHandle.h
Created:	2006/14/05
Author:		Jason W H Wong <jason.wong@chem.ox.ac.uk>
	
Purpose:	Implements reading of mzXML files and writing of DTA files.

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

#ifndef FILEHANDLE_H
#define FILEHANDLE_H

#include <vector>
#include <string>

using namespace std;

struct peakData{
  double mz;
  double exact_mz;
  double intn;
  bool has_comp;
  bool has_iso;
  bool has_h2o_nh3;
  bool has_a;
  peakData *comp;
  int rank;
  double score;
  double prob;
  int posi;
};

struct AA{
  string name;
  float mass;
  peakData *peak;
  bool hasPeak;
};

struct scan
{
  long scanNum;
  int level;
 // char*	scanType;
  long	numPakets;
 // double	RT;
  double	lowMass;
  double	highMass ;
  double	TIC;
  double	basePeakMass;
  double	basePeakIntensity;
//  long	channel;
//  long	uniformTime;
//  double	frequency;
//  double	collision;
  double	preIntn;
  double	preMZ;
  vector <float> mz;
  vector <float> intn;
  vector <peakData> peaks;
  double sumPeaksIntn;
  int charge[2];
  double TICFore;
  double TICBack;
  bool good;
  double parentMass;
  double parent1;
  double parent2;
  double parent3;
  int match2;
  int match3;
  int num_a;
  int num_h2o;
  int num_iso;
  double SNR;
  int num_r_peaks;

  int commonDiff;
  double secBasePeakIntn;
  double normMatch1;
  double normMatch2;
  double normMatch3;

  int num_peaks;
  double normTIC;
  double goodSegs;
  double intnRatio1;
  double intnRatio20;
  double comple;
  double iso_ratio;
  double h2o_ratio;
  double AA_diff_ratio;

  double discrim;
  double Fp;
  double pF;
  double Fm;
  double mF;
  
  double PP;
};

struct parameters
{
  double p_num_peaks;
  double p_normTIC;
  double p_goodSegs;
  double p_intnRatio1;
  double p_intnRatio20;
  double p_comple;
  double p_iso_ratio;
  double p_h2o_ratio;
  double p_AA_diff_ratio;
  double p_const;
  
  double em_p_avg;
  double em_p_std;
  double em_p_prior;
  double em_m_avg;
  double em_m_std;
  
  double em_prior_ulimit;
  double em_prior_llimit;
  double em_p_avg_ulimit;
  double em_p_avg_llimit;
  double em_p_std_ulimit;
  double em_p_std_llimit;
  double em_m_avg_ulimit;
  double em_m_avg_llimit;
  double em_m_std_ulimit;
  double em_m_std_llimit;
};

class fileHandle
{

 public:

  fileHandle();
  virtual ~fileHandle();

  int importFile(char *_filename, int isMZXML);
  int importRegParmas(char *_filename);
  int importFound(char *_filename);
  int exportSpecEval(const char* path);
  int exportDTA(char *foldername,unsigned  int begin,unsigned  int end, bool guessmulti, bool extracharge);
  int exportBadCSV(char *_filename, unsigned int start,unsigned  int end);

  string filename;
  vector<scan> all;
  
  int m_tSpectraTotal;

  parameters param;

  vector <scan*> pri;   // MS data
  vector <scan*> sec;	// MSMS data
  vector <scan*> ter;	// MS3 data
  
  vector <int> found_scans;
};
#endif
