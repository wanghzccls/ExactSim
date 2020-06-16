#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "SimStruct.h"
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>


int mkpath(string s, mode_t mode=0755){
	size_t pre=0, pos;
	string dir;
	int mdret;
	if(s[s.size()-1]!='/'){
		s+='/';
	}
	while((pos=s.find_first_of('/',pre))!=string::npos){
		dir=s.substr(0,pos++);
		pre=pos;
		if(dir.size()==0) continue;
		if((mdret=::mkdir(dir.c_str(),mode)) && errno!=EEXIST){
			return mdret;
		}
	}
	return mdret;
}



void usage() {
    cout << "./EXSim -d <dataset> -f <filelabel> -algo <Algorithm> [-e epsilon (default 0.001)] [-qn querynum (default 50)]" << endl;
}

//查看参数是否合规
int check_inc(int i, int max) {
    if (i == max) {
        usage();
        exit(1);
    }
    return i + 1;
}


int main(int argc, char **argv){
    int i = 1;
    char *endptr;
    string filename;
    string outputname;
    string filelabel;
    string queryname;
    int querynum = -1;
    double eps = 0.001;
    double cvalue = 0.6;
    double level = 1000;
    string algo = "ExactSim";
    if(argc < 1){
        usage();
        exit(1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-d")) {
            i = check_inc(i, argc);
            filename = argv[i];
        } else if (!strcmp(argv[i], "-f")) {
            i = check_inc(i, argc);
            filelabel = argv[i];
        } else if (!strcmp(argv[i], "-algo")) {
            i = check_inc(i, argc);
            algo = argv[i];
        } else if (!strcmp(argv[i], "-e")) {
            i = check_inc(i, argc);
            eps = strtod(argv[i], &endptr);
            if ((eps == 0 || eps > 1) && endptr) {
                cerr << "Invalid eps argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-c")) {
            i = check_inc(i, argc);
            cvalue = strtod(argv[i], &endptr);
            if ((cvalue == 0 || cvalue > 1) && endptr) {
                cerr << "Invalid c argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-qn")) {
            i = check_inc(i, argc);
            querynum = strtod(argv[i], &endptr);
            if ((querynum < 0) && endptr) {
                cerr << "Invalid walknum argument" << endl;
                exit(1);
            }
        } else if (!strcmp(argv[i], "-o")) {
            i = check_inc(i, argc);
            outputname = argv[i];
        } else {
            usage();
            exit(1);
        }
        i++;
    }
    
   
    
    SimStruct sim = SimStruct(filename, filelabel, eps, cvalue);
    if(querynum == -1){
        querynum = max((int)sim.vert,50); 
    }

    
    //Generate query nodes
    //
    if(algo == "GEN_QUERY"){
        ofstream data_idx("query/" + filelabel + ".query");
        for(int i = 0; i < querynum; i++){
		data_idx << sim.R.generateRandom() % sim.vert << "\n";
        }
        data_idx.close();
    }

  
    //ExactSim Algorithm
    //
    if(algo == "ExactSim"){
	queryname = "query/" + filelabel + ".query";    
        ifstream query(queryname);
        cout<<"querynum="<<querynum<<endl;
	for(int i=0;i<querynum;i++){
	    int nodeId;
            query >> nodeId;
            cout<<i<<": "<<nodeId<<endl;
	    stringstream ss,ss_dir;
	    ss_dir<<"results/"<<filelabel<<"/"<<eps<<"/";
	    mkpath(ss_dir.str());
	    ss<<ss_dir.str()<<nodeId<<".txt";
    	    outputname = ss.str();
	    sim.ExactSim(nodeId, outputname, level);
	}
	query.close();

	cout<< "\n====ExactSim done!====\n"<<endl;
	double avg_querytime=sim.avg_time/(double)querynum;
	cout<<"Average Query Time: "<<avg_querytime<<" s\n"<<endl;
    }
    
    return 0;
}
