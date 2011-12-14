/***************************************************************************************************
File:		fileHandle.cpp
Created:	2006/19/12
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/timeb.h>
#include <ctime>
#include <algorithm>
#include <set>
#include "loadmspectrum.h"


#include "fileHandle.h"

fileHandle::fileHandle(){
	param.p_num_peaks = 0;
	param.p_normTIC = 0;
	param.p_goodSegs = 0;
	param.p_intnRatio1 = 0;
	param.p_intnRatio20 = 0;
	param.p_comple = 0;
	param.p_iso_ratio = 0;
	param.p_h2o_ratio = 0;
	param.p_AA_diff_ratio = 0;
	param.p_const = 0;
	param.em_p_avg = 5.0;
	param.em_p_std = 3.26216318;
	param.em_p_prior = 0.2;
	param.em_m_avg = -5.0;
	param.em_m_std = 2.2;
	param.em_prior_ulimit = 1;
	param.em_prior_llimit = 0;
	param.em_p_avg_ulimit = 10000;
	param.em_p_avg_llimit = -10000;
	param.em_p_std_ulimit = 10000;
	param.em_p_std_llimit = -10000;
	param.em_m_avg_ulimit = 10000;
	param.em_m_avg_llimit = -10000;
	param.em_m_std_ulimit = 10000;
	param.em_m_std_llimit = -10000;
}


fileHandle::~fileHandle(){

}

int fileHandle::importFile(char *_filename, int isMZXML)
{

  string strValue = _filename;
  filename = _filename;
  string::size_type pos = filename.find_last_of("/");
  if (pos != string::npos){
     filename.erase(0,pos+1);
  }
	FILE *pStream;
	pStream = fopen(_filename,"r");
	char *pValue = new char[1024];
	size_t tRead = 256;
	if(pStream)	{
		tRead = fread(pValue,1,256,pStream);
	}
	if(pStream)	{
		if(pValue[0] == 1 && pValue[1] == -95 || (pValue[3] == 'F' && pValue[5] == 'i' && pValue[7] == 'n') )	{
			cout << "\nFailed to read spectrum file: " << strValue.c_str() << "\n";
			cout << "Most likely cause: using a Finnigan raw spectrum.\nUse ASCII formatted dta, pkl or mgf files ONLY!\n\n";
			cout.flush();
			delete pValue;
			return false;
		}
		else	{
			size_t d = 0;
			while(d < tRead)	{
				if(pValue[d] == '\0')	{
					cout << "\nFailed to read spectrum file: " << strValue.c_str() << "\n";
					cout << "Most likely cause: using a binary spectrum file.\nUse ASCII formatted dta, pkl or mgf files ONLY!\n\n";
					cout.flush();
					delete pValue;
					return false;
				}
				d++;
			}
		}
		fclose(pStream);
	}
	delete pValue;


// Start section of code from XTandem
  if(isMZXML == 1 || (strstr(_filename,"mzXML") || strstr(_filename,"MZXML") || strstr(_filename,"XML") || strstr(_filename,"xml") && strstr(_filename,"mzdata.xml")==0))	{
    loadmzxml ldMzxml(all);
    if(ldMzxml.open(strValue))	{
      ldMzxml.get();
	  m_tSpectraTotal = all.size();
    }

  }

  else	{
   loadmzdata ldMzdata(all);
    if(ldMzdata.open(strValue))  {
      ldMzdata.get();
	  m_tSpectraTotal = all.size();
    }
  }

// end section of code from XTandem.
  for (unsigned int i = 0; i < all.size();++i)
  {
	all.at(i).num_peaks = 0;
	all.at(i).normTIC=0;
	all.at(i).goodSegs=0;
	all.at(i).intnRatio1=0;
	all.at(i).intnRatio20=0;
	all.at(i).comple=0;
	all.at(i).iso_ratio=0;
	all.at(i).h2o_ratio=0;
	all.at(i).AA_diff_ratio=0;
	all.at(i).TIC = 0;
	all.at(i).charge[0] = 1;
	all.at(i).charge[1] = 0;
	all.at(i).discrim=0;
	all.at(i).Fp=0;
	all.at(i).pF=0;
	all.at(i).Fm=0;
	all.at(i).mF=0;
	all.at(i).TICFore = 0;
	all.at(i).TICBack = 0;
	all.at(i).sumPeaksIntn = 0;
//  double	preIntn;
//  double	preMZ;
  	all.at(i).good=false;
//  double parentMass;
//  double parent1;
//  double parent2;
  	//double parent3;
  	all.at(i).match2=0;
 	all.at(i).match3=0;
 	all.at(i).num_a=0;
  	all.at(i).num_h2o=0;
  	all.at(i).num_iso=0;
  	all.at(i).SNR=0;
  	all.at(i).num_r_peaks=0;
  	all.at(i).commonDiff=0;
  	all.at(i).secBasePeakIntn=0;
//  double normMatch1;
//  double normMatch2;
//  double normMatch3;
  	
    if (all.at(i).level == 2)
      sec.push_back(&all.at(i));
  }

  return 1;
}

int fileHandle::importRegParmas(char *_filename)
{
	cout << "Reading discriminant parameters..."<<endl<<flush;
	string strValue = _filename;
	ifstream in((const char*)_filename);
	char pValue[256]="";
	if(in.is_open()){
		in.getline(pValue,256);
		if(strstr(pValue,"msmsEval")==NULL)	{
			cout << "\nFailed to read parameter file: " << strValue.c_str() << "\n";
			cout << "Please use the example file that came with the source as a template!\n\n";
			cout.flush();
			return false;
		}
		else	{
			in.getline(pValue,256);
			while(!in.eof())	{
				char tmp[256];
				strcpy(tmp,pValue);
				char * pch = strtok(tmp," \t");
				if (pch != NULL && strstr(pch,"#") == NULL){
					strcpy(pValue,pch);
					pch = strtok(NULL," \t");
					if (pch == NULL){
						cout << "Error in parameter file at: "<< pValue<<endl;
						return 0;
					}
					double val = atof(pch);
					if(strcmp(pValue,"P_NUM_PEAKS")==0)
						param.p_num_peaks = val;
					else if (strcmp(pValue,"P_NORM_TIC")==0)
						param.p_normTIC = val;
					else if (strcmp(pValue,"P_GOOD_SEGS")==0)
						param.p_goodSegs = val;
					else if (strcmp(pValue,"P_INTN_RATIO_1")==0)
						param.p_intnRatio1 = val;
					else if (strcmp(pValue,"P_INTN_RATIO_20")==0)
						param.p_intnRatio20 = val;
					else if (strcmp(pValue,"P_COMPLEMENT")==0)
						param.p_comple = val;
					else if (strcmp(pValue,"P_ISO_RATIO")==0)
						param.p_iso_ratio = val;
					else if (strcmp(pValue,"P_H2O_RATIO")==0)
						param.p_h2o_ratio = val;
					else if (strcmp(pValue,"P_AA_DIFF_RATIO")==0)
						param.p_AA_diff_ratio = val;
					else if (strcmp(pValue,"P_CONST")==0)
						param.p_const = val;
					else if (strcmp(pValue,"EM_P_AVG")==0)
						param.em_p_avg = val;
					else if (strcmp(pValue,"EM_P_STD")==0)
						param.em_p_std = val;
					else if (strcmp(pValue,"EM_P_PRIOR")==0)
						param.em_p_prior = val;
					else if (strcmp(pValue,"EM_M_AVG")==0)
						param.em_m_avg = val;
					else if (strcmp(pValue,"EM_M_STD")==0)
						param.em_m_std = val;
					else if (strcmp(pValue,"EM_PRIOR_ULIMIT")==0)
						param.em_prior_ulimit = val;
					else if (strcmp(pValue,"EM_PRIOR_LLIMIT")==0)
						param.em_prior_llimit = val;
					else if (strcmp(pValue,"EM_P_AVG_ULLIMIT")==0)
						param.em_p_avg_ulimit = val;
					else if (strcmp(pValue,"EM_P_AVG_LLIMIT")==0)
						param.em_p_avg_llimit = val;
					else if (strcmp(pValue,"EM_P_STD_ULIMIT")==0)
						param.em_p_std_ulimit = val;
					else if (strcmp(pValue,"EM_P_STD_LLIMIT")==0)
						param.em_p_std_llimit = val;
					else if (strcmp(pValue,"EM_M_AVG_ULIMIT")==0)
						param.em_m_avg_ulimit = val;
					else if (strcmp(pValue,"EM_M_AVG_LLIMIT")==0)
						param.em_m_avg_llimit = val;
					else if (strcmp(pValue,"EM_M_STD_ULIMIT")==0)
						param.em_m_std_ulimit = val;
					else if (strcmp(pValue,"EM_M_STD_LLIMIT")==0)
						param.em_m_std_llimit = val;
					else
						cout << "Unrecognized variable with value: "<<val<<endl<<flush;
				}
				in.getline(pValue,256);
			}
		}
		if (param.p_const == 0)
			cout << "Warning: P_CONST = 0, possibly due to absence of new line at end of parameter file"<<endl<<flush;
	}
	else
		cout << "Cannot open parameter file."<<endl;
	in.close();
	return 1;
}

int fileHandle::importFound(char* _filename)
{
	ifstream in((const char*)_filename);
	char buffer[65536]="";
	if (in.is_open()){
		if (! in.eof()){
			in.getline (buffer,65536);
			if (!iswalnum(buffer[0]) && !iswpunct(buffer[0]) && !iswspace(buffer[0]) && !iswdigit(buffer[0])){
				cout << "BAD FILE"<<endl;
				return -1;
			}
		}
		in.getline (buffer,65536);
		while (!in.eof()){
			if (iswdigit(buffer[0])){
				char * pch;
				pch = strtok(buffer,",");
				while (pch != NULL){
					int scanNum;
					scanNum=atoi(pch);
					found_scans.push_back(scanNum);
					pch = strtok(NULL,",");
				}
			}
			in.getline (buffer,65536);
		}
	}
	else
		{ cout << "Error opening file"; return 0; }
  return 1;
}

int fileHandle::exportSpecEval(const char* path)
{
	ofstream out;
	out.open(path,ios::in|ios::trunc);

	out << "Scan #,Parent, NPeaks,NormTIC,GoodSegs,IntnRatio1%,IntnRatio20%,complements,IsoRatio,H2ORatio,AADiffRatio,discriminant,P(+|D),Z_prob"<<endl;
	for (unsigned int i = 0; i < sec.size(); ++i){
		out<< sec.at(i)->scanNum <<','<<sec.at(i)->preMZ<<','<<sec.at(i)->mz.size()<<','<<sec.at(i)->normTIC<<','<<sec.at(i)->goodSegs<<','<<sec.at(i)->intnRatio1<<','<<sec.at(i)->intnRatio20<<','<<sec.at(i)->comple<<','<<sec.at(i)->iso_ratio<<','<<sec.at(i)->h2o_ratio<<','<<sec.at(i)->AA_diff_ratio<<','<<sec.at(i)->discrim<<','<<sec.at(i)->pF<<','<<sec.at(i)->PP<<endl;
	}
	out.close();
	return 1;
}

int fileHandle::exportDTA(char *foldername, unsigned int start,unsigned  int end, bool guessmulti, bool extracharge)
{
	  if (sec.empty())
	    return 0;
	  ofstream out;
	  char t[256]="";
	  string flname;
	  string fname;
	  strcat(foldername,"/");
	  cout << "Exporting DTA to: "<<foldername<< " From: "<< start << " To: " <<end<<endl;
	  int n1Files = 0;
	  int n2Files = 0;
	  int n3Files = 0;
	  int n4Files = 0;
	  int n5Files = 0;
	  int rej = 0;
	  int found_file = 0;
	  if (!guessmulti){
	  	for (unsigned int i = start; i < end; ++i){
			if (sec.at(i)->charge[0] != 1){
	  			sec.at(i)->charge[0] = 4;
			}
	  	}
	  }
	  for (unsigned int i = start; i < end; ++i){
	  	bool found = false;
	  	for (unsigned int j = 0; j < found_scans.size() && !found; ++j){
	  		if (found_scans.at(j) == sec.at(i)->scanNum)
	  			found = true;
	  	}
	  	if (found){
	  		found_file++;
	  		rej++;
	  	}
	  	else{
		  	if (i %250 == 0){
		  		cout << ".";
				cout.flush();
			}
		  	if (sec.at(i)->good){
		    		sprintf(t, "%d", sec.at(i)->scanNum);
		    		if (sec.at(i)->charge[0] == 4){
		        		fname = filename+'.'+t+'.'+t+'.'+"2.dta";
					flname = foldername+fname;
					out.open(flname.c_str(),ios::in|ios::trunc);
					n2Files++;
					float parent = (sec.at(i)->preMZ*2)-1;
					out<<parent<<" 2"<<endl;
					for (unsigned int j =0; j < sec.at(i)->peaks.size(); ++j){
		 	  			out<<sec.at(i)->peaks.at(j).exact_mz<<" "<<sec.at(i)->peaks.at(j).intn<<endl;
			  		}
					out.close();
					fname = filename+'.'+t+'.'+t+'.'+"3.dta";
					flname = foldername+fname;
					out.open(flname.c_str(),ios::in|ios::trunc);
					n3Files++;
					parent = (sec.at(i)->preMZ*3)-2;
					out<<parent<<" 3"<<endl;
					for (unsigned int j =0; j < sec.at(i)->peaks.size(); ++j){
			  			out<<sec.at(i)->peaks.at(j).exact_mz<<" "<<sec.at(i)->peaks.at(j).intn<<endl;
					}
					out.close();
					if (extracharge){
			        		fname = filename+'.'+t+'.'+t+'.'+"4.dta";
						flname = foldername+fname;
						out.open(flname.c_str(),ios::in|ios::trunc);
						n4Files++;
						float parent = (sec.at(i)->preMZ*4)-3;
						out<<parent<<" 4"<<endl;
						for (unsigned int j =0; j < sec.at(i)->peaks.size(); ++j){
			 	  			out<<sec.at(i)->peaks.at(j).exact_mz<<" "<<sec.at(i)->peaks.at(j).intn<<endl;
				  		}
						out.close();
						fname = filename+'.'+t+'.'+t+'.'+"5.dta";
						flname = foldername+fname;
						out.open(flname.c_str(),ios::in|ios::trunc);
						n5Files++;
						parent = (sec.at(i)->preMZ*5)-4;
						out<<parent<<" 5"<<endl;
						for (unsigned int j =0; j < sec.at(i)->peaks.size(); ++j){
				  			out<<sec.at(i)->peaks.at(j).exact_mz<<" "<<sec.at(i)->peaks.at(j).intn<<endl;
						}
						out.close();
					}					
		    	}
		   		else{
					fname = filename+'.'+t+'.'+t+'.';
					sprintf(t, "%d", sec.at(i)->charge[0]);
					fname = fname+t+".dta";
					flname = foldername+fname;
					out.open(flname.c_str(),ios::in|ios::trunc);
					if (sec.at(i)->charge[0] == 1)
				   		n1Files++;
					else if (sec.at(i)->charge[0] == 2)
				  		n2Files++;
					else if (sec.at(i)->charge[0]  ==3)
				   		n3Files++;
					float parent = (sec.at(i)->preMZ*sec.at(i)->charge[0])-(sec.at(i)->charge[0]-1);
					out<<parent<<" "<<sec.at(i)->charge[0]<<endl;
					for (unsigned int j =0; j < sec.at(i)->peaks.size(); ++j){
				  		out<<sec.at(i)->peaks.at(j).exact_mz<<" "<<sec.at(i)->peaks.at(j).intn<<endl;
					}
					out.close();
				}
		  	}
	  		else
				rej++;
		  }
	  }
	  cout <<endl;
	  cout<<"# 1+ files generated: "<<n1Files<<endl;
	  cout<<"# 2+ files generated: "<<n2Files<<endl;
	  cout<<"# 3+ files generated: "<<n3Files<<endl;
	  cout<<"# 4+ files generated: "<<n4Files<<endl;
	  cout<<"# 5+ files generated: "<<n5Files<<endl;
	  cout<<"Found scans filtered: "<<found_file<<endl;
	  cout<<"Total scans filtered: "<<rej<<endl;
	  return 1;
}

int fileHandle::exportBadCSV(char *_filename, unsigned int start,unsigned  int end)
{
	  if (sec.empty())
	    return 0;
	  char t[256]="";
	  cout << "Exporting \"Bad\" CSV file to: "<<_filename <<endl;
	  ofstream out;
	  out.open(_filename,ios::in|ios::trunc);
	  int bad_scans = 0;
	  for (unsigned int i = start; i < end; ++i){
		if (!sec.at(i)->good){
			bad_scans++;
			out<<sec.at(i)->scanNum;
			bool found = false;
			for (unsigned int j = i+1; j < end && !found; ++j){
				if (!sec.at(j)->good){
					found = true;
				}
			}
			if (found)
				out<< ", ";
			else
				out<<endl;
		}
	  }
	  cout << "Total bad scans: "<<bad_scans<<endl<<endl;
}
