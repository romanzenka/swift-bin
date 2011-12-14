/***************************************************************************************************
File:		main.cpp
Created:	2007/9/8
Author:		Jason W H Wong <jason.wong@chem.ox.ac.uk>
	
Purpose:	The interface for msmsEval.

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

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <fstream>
#include "filterer.h"
#include "fileHandle.h"

using namespace std;

int usage();

int main(int argc, char *argv[])
{
  filterer fil;
  int c;

  char *in_filename = NULL;
  char *out_folder = ".";
  char *params_filename = NULL;
  char *eval_filename = NULL;
  char *em_filename = NULL;
  char *found_filename = NULL;
  char *bad_filename = NULL;
  int sflag = 0;
  int tflag = 0;
  int rflag = 0;
  int bflag = 0;
  int eflag = 0;
  int dflag = 0;
  int mflag = 0;
  int oflag = 0;
  int fflag = 0;
  int pflag = 0;
  int lflag = 0;
  bool xflag = false;
  bool cflag = false;
  float pVal = 0.01;
  float pF = 0.9;
  unsigned int first_scan = 0;
  unsigned int last_scan = 0;
  bool qflag = false;
  int isMZXML = 0;
  int zflag = 0;
  
  if (argc == 1) usage();

  while ((c = getopt(argc,argv,"hqcxr:p:s:t:b:e:d:o:m:f:z:l:")) != -1) {
    switch(c){
    case 's': //source file
      sflag = 1;
      in_filename = optarg;
      break;
    case 't': // destination dir
      tflag = 1;
      out_folder = optarg;
      break;
    case 'l': // bad file
      lflag = 1;
      bad_filename = optarg;
      break;
    case 'd': // params file
      dflag = 1;
      params_filename = optarg;
      break;
    case 'f': // found file
      fflag = 1;
      found_filename = optarg;
      break;
    case 'o': // eval file
      oflag = 1;
      eval_filename = optarg;
      break;
    case 'm': // em file
      mflag = 1;
      em_filename = optarg;
      break;
    case 'r': // removed probability
      rflag = 1;
      pVal = atof(optarg);
      break;
    case 'p': // pF probability
      pflag = 1;
      pF = atof(optarg);
      break;
    case 'b': // start scan
      bflag = 1;
      first_scan = (unsigned int)atoi(optarg);
      break;
    case 'e': // last scan
      eflag = 1;
      last_scan = (unsigned int)atoi(optarg);
      break;
    case 'z': // file format
      zflag = 1;
      if (strstr(optarg,"mzXML") || strstr(optarg,"mzxml") || strstr(optarg,"MZXML") )
         isMZXML = 1;
      else if ( strstr(optarg,"mzdata") || strstr(optarg,"mzData") || strstr(optarg,"MZDATA") )
         isMZXML = 2;
      break;
    case 'q':
      qflag = true;
      break;
    case 'c':
      cflag = true;
      break;
    case 'x':
      xflag = true;
      break;
    case 'h':  //help
      usage();
    default: 
      usage();
    }
  }
  
  if (sflag == 0){
  	cout << "Please enter input file name (-s)."<<endl;
  	exit(0);
  }
  if (dflag == 0){
  	cout << "Please enter discriminant parameters file name (-d)."<<endl;
  	exit(0);
  }
  if (pflag == 1 && rflag == 1){
  	cout << "Specify only either -p (filter by p(+|D)) or -r (filter by fraction of spectra removed)."<<endl;
  	exit(0);
  }
 
  if (pVal > 1)
	pVal = 1;
  else if (pVal <= 0)
	pVal = 0;

  cout << "\nmsmsEval 1.3 (03-Aug-07) \nby Jason W. H. Wong, Cartwright Lab, PTCL, University of Oxford, 2006-7"<<endl<<endl;

 	fileHandle fHandle;
 	filterer filter;
 	double cut = filter.setCutoff(pVal);
 	cout<<"source: "<<in_filename<<endl;
 	cout<<"params: "<<params_filename<<endl;
	if (tflag != 0){
		cout <<"dest: "<<out_folder;
	}
	if (pflag != 0){
		cout <<" (pF cutoff="<<pF<<")"; 
	} else {
		cout <<" (p cutoff="<<cut<<")"; 
	}
 	cout << endl;
    int s = fHandle.importFile(in_filename,isMZXML);
    if (s == 0)
		return -1;
    s = fHandle.importRegParmas(params_filename);
    if (s == 0)
		return -1;
	if (fflag){
		s = fHandle.importFound(found_filename);
   		if (s == 0)
			return -1;
	}
    string temp_name;
    string tmpin = in_filename;

    temp_name+=tmpin;
    temp_name+="_eval.csv";
    if (oflag)
    	temp_name = eval_filename;

    string em;
    em+=tmpin;
    em+="_em.csv";
    if (mflag)
    	em = em_filename;

 if (s!= 0 && fHandle.sec.size() != 0){
    cout << "Spectra to be analyzed: "<< fHandle.sec.size()<<" ";
    cout.flush();
    filter.analyze(fHandle.sec,fHandle.param);
    filter.run_EM(fHandle.sec,em.c_str(),fHandle.param,qflag);
    filter.sortGB(fHandle.sec);
    if (pflag && !rflag)
    	filter.filterData(fHandle.sec,pF,true);
    else
    	filter.filterData(fHandle.sec,pVal,false);
    if (tflag != 0){
    	if (first_scan < last_scan || first_scan < 0)
       		first_scan = 0;
    	if (last_scan > fHandle.sec.size()|| last_scan == 0)
       		last_scan = fHandle.sec.size();
    	fHandle.exportDTA(out_folder,first_scan,last_scan,cflag,xflag);
    }
    if (lflag != 0){
    	if (first_scan < last_scan || first_scan < 0)
       		first_scan = 0;
    	if (last_scan > fHandle.sec.size()|| last_scan == 0)
       		last_scan = fHandle.sec.size();
    	fHandle.exportBadCSV(bad_filename,first_scan,last_scan);
    }
    if (!qflag)
    	fHandle.exportSpecEval(temp_name.c_str());
    cout << "Analysis results output to: "<< temp_name.c_str()<<endl;
 }
 else
   cout << "\nFile is empty or does not exist."<<endl<<endl;
}

int usage() {
	printf("\nmsmsEval 1.3 (03-Aug-07) \nby Jason W. H. Wong, Cartwright Lab, PTCL, Univeristy of Oxford, 2006-7\n");
	printf("\nUsage:msmsEval <-s source_file> <-d param_file>[-t target_dir][-l bad_filename]\n");
	printf("       	          [-f found_scans][-o eval_output_file][-m em_output_file]\n");
	printf("       	          [-p pF probability][-r remove probability]\n");
	printf("                  [-b first_scan][-e last_scan][-z file format][-c][-q][-x] \n");
	printf("	----------------------------------------------------------------\n");
	printf("       	(source_file)   path of source file (mzXML only)\n");
	printf("       	(param_file)   path of discriminant/em parameters file\n");
	printf("       	-t   where \"target_dir\" is path of target directory for DTA export\n");
        printf("       	-l   outputs CSV file contain all \"Bad\" scans\n");
	printf("       	-f   where \"found_scans\" is a file with scans that are already annotated\n");
	printf("       	-o   where \"eval_output_file\" is path of target file name for the\n"); 
	printf("             evaulation summary\n");
	printf("       	-m   where \"em_output_file\" is path of target file name for EM\n"); 
	printf("             algorithm summary\n");
	printf("        -p   where \"pF probability\" is a float specifying the target p(+|D)\n"); 
	printf("             cutoff for DTA files generation (default: 0.9)\n");
	printf("        -r   where \"remove fraction\" is a float specifying the  \n"); 
	printf("             estimated fraction of identifiable spectra a user is willing to \n"); 
	printf("             potentially 'scrifice' when generating DTA files (default: 0.01)\n");	
    printf("       	-b   specify \"first_scan\" for exporting a limited range of scans\n");
    printf("       	-e   specify \"last_scan\" for exporting a limited range of scans\n");
    printf("       	-z   force \"file_format\" regardless of file extension.\n");
    printf("                 either: \"mzZML\" or \"mzData\"\n");
    printf("       	-c   if specified, msmsEval will try to guess the charge of a spectrum\n"); 
    printf("             for DTA output\n");
    printf("       	-x   generate extra charge states (+4/+5) when making dta files\n");
    printf("       	-q   suppress summary outputs\n");
	printf("\nUsage examples: \n");
	printf("msmsEval -s /home/guest/raw_data/example.mzXML -d ./msmsEval.params \n");
	printf("msmsEval -s /home/guest/raw_data/example.mzXML -t ./home/guest/filtered_dtas/ \n");
	printf("         -d ./msmsEval.params -p 0.9 -c -z mzXML\n\n");
	printf("BY default, the evaluation file and the results of the em algorithm will be \n");
	printf("outputed in the same directory as the source file named '<source file>_eval.csv'\n");
	printf("and '<source_file>_em.csv' respectively\n\n");	

	exit(1);
}
