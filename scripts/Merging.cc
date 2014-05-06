#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <map>

#include <TFile.h>

using namespace std;

int main() {


  system("source GetList.sh");


  ifstream dsFile("log_dataset", ios::in );

  if(dsFile) {
    string line;
    while( std::getline(dsFile, line) ) {

      ifstream ifile( ("log_"+line).c_str(), ios::in );
      
      if( ifile ) {
	
	map<int, std::pair<int,string> > files;
	map<int, std::pair<int,string> >::iterator iter;

	string f;
	
	while( std::getline(ifile, f) ) {

	  int posB = f.find("Ntuple");
	  int posE = f.find("_",posB+6);

	  int pos1 = f.find("_",posE+1);
	  int pos2 = f.find("_",pos1+1);

	  string nJob = f.substr(posE+1, pos1-posE-1);
	  string jIt = f.substr(posE+1, pos1-posE-1);
	  int NJob = atoi( nJob.c_str() );
	  int jobIt = atoi( jIt.c_str() );

	  cout<<f<<"   "<<NJob<<endl;

	  string fname = "/tmp/mmarionn/"+line+"/"+f;
	  TFile* tmpF= TFile::Open( (fname).c_str() );
	  if(tmpF==NULL) {delete tmpF; continue;}
	  delete tmpF;
	  

	  std::pair<int,string> k(jobIt,f);
	  
	  iter = files.find( NJob );

	  if( iter != files.end() ) {
	    cout<<" file in double! "<<(*iter).first<<"   "
		<<(*iter).second.second<<"     "<<(*iter).second.first
		<<"   "<<jobIt<<endl;
	    
	    if( (*iter).second.first < jobIt )
	      {
	    files.erase( NJob );
	    files[ NJob ] = k;
	      }
	  }
	  else
	    files[ NJob ] = k;
      
	}
	
	string hadd="hadd /tmp/mmarionn/"+line+".root ";
	
	for(iter = files.begin();iter != files.end(); iter++) {
	  string nf = " /tmp/mmarionn/"+line+"/"+(*iter).second.second;
	  
	  hadd += nf;
	  
	}
	cout<<"\n\n =============> Mergin "<<line<<" <===================  "<<endl;
	system( hadd.c_str() );
      }
      else { cout<<" Raté "<<line<<endl; }
      
    }
  }
  

}
