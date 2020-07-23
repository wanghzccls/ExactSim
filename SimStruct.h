#ifndef SIMSTRUCT_H
#define SIMSTRUCT_H

#define INT_MAX 32767

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include <chrono>
#include <cstdlib>
#include <errno.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <unistd.h>


void process_mem_usage(double& vm_usage, double& resident_set){
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}



class SimStruct{
   public:
  	Graph g; //Graph class
  	Random R;    
 	int vert; //number of nodes
 	string filelabel;
  	double sqrtC;                                                      
  	double C_value;
  	double eps; 
  	double delta;
  	double avg_time;
  	double avg_peakmem; 
  	int *vecpos;
  	double *reserve;
  	double *residue;
  	double *newReserve;
  	double *newResidue;
  	bool* isInArray;
  	long power_walk_num;
  	int power_ite_num;
  	unordered_map<int, unordered_map<int,double> > P_lw;//transition probablity used in function cal_Zw()
  	unordered_map<int, unordered_map<int,unordered_map<int,double> > > P_ly;//<final_node<level,<levelnode,pr>>>
  	unordered_map<int, unordered_map<int, vector<pair<int,double> > > > P_ally;//<level <node <final_node, pr> > >
  	unordered_map<int, unordered_map<int, double> > Fl;//used in function cal_Zw()
  	unordered_map<int, unordered_map<int, double> > answer;
  	vector<vector<int> > forwardNode;
  	vector<int> pathnum_set;
  	vector<int> new_pathnum_set;
  	vector<int> candidate_set_ford;
  	vector<int> new_candidate_set_ford;   
	

  SimStruct(string name, string file_label, double epsilon, double cvalue) {
    	R = Random();
    	C_value = cvalue;
    	sqrtC = sqrt(C_value);
    	eps=epsilon;
    	delta=0.01;

	filelabel = file_label;
    	g.inputGraph(name);
    	vert = g.n;
    	power_walk_num=(long)((C_value*0.00001)*(log(vert/delta)/log(2))/(eps*eps));
    	power_ite_num=(int)((log(eps*(1-C_value))/log(C_value)-1)+1);

    	cout<<"eps: "<<eps<<endl;
    	cout<<"decay c: "<<C_value<<endl;
    	cout<<"Sampling times: "<<power_walk_num<<endl;
    	cout<<"Iteration times: "<<power_ite_num<<endl;

    	avg_time = 0;
    	avg_peakmem=0.0;

    	reserve = new double[vert];
    	residue = new double[vert];
    	newReserve = new double[vert];
    	newResidue = new double[vert];
    	isInArray = new bool[vert];
    	vecpos=new int[vert];
    	for(int i = 0; i < vert; i++){
      		isInArray[i] = false;
      		residue[i] = 0;
      		newResidue[i] = 0;
      		reserve[i] = 0;
      		newReserve[i] = 0;
      		vecpos[i]=0;
    	}
    	srand(unsigned(time(0)));
    	cout << "====init done!====\n" << endl;
  }

  ~SimStruct() {
    	delete[] reserve;
    	delete[] residue;
    	delete[] newReserve;
    	delete[] newResidue;
    	delete[] isInArray;
    	delete[] vecpos;
    	vector<vector<int> >().swap(forwardNode);
    	vector<int>().swap(pathnum_set);
    	vector<int>().swap(new_pathnum_set);
    	vector<int>().swap(candidate_set_ford);
    	vector<int>().swap(new_candidate_set_ford);
  }


  //calculate transition probablity
  unordered_map<int, unordered_map<int, double> > forwardPushByLevel(int u, int maxLevel, int &deterlevel, long walk_num, vector<vector<int> > &forwardNode){
	answer.clear();
	forwardNode.clear();
        pathnum_set.clear();
  	new_pathnum_set.clear();
    	
	residue[u] = 1;
    	candidate_set_ford.push_back(u);
	pathnum_set.push_back(1);
	int tempLevel = 0;
	bool deter_flag=false;
	
	while(1){
      		long candidateSetSize=candidate_set_ford.size();
		if(candidateSetSize<=0){
			deterlevel=tempLevel-1;
			break;
		}
		long current_pathnum=0;
		long residue_count=0;
		if(tempLevel>1){
			for(int j=0;j<candidateSetSize;j++){
				int tempNode=candidate_set_ford[j];
				int inSize=g.getInSize(tempNode);
				current_pathnum+=pathnum_set[j]*inSize;
				if(current_pathnum*current_pathnum>=walk_num){
					deterlevel=tempLevel;
					deter_flag=true;
					break;;
				}
			}
		}
      		for(int j = 0; j < candidateSetSize; j++){
			int tempNode = candidate_set_ford[j];
			answer[tempLevel][tempNode]=residue[tempNode];
			if(tempLevel==maxLevel||deter_flag==true){
				residue[tempNode] = 0;  
				continue;
			}			
			int inSize = g.getInSize(tempNode);
			for(int k = 0; k < inSize; k++){
	  			int newNode = g.getInVert(tempNode, k);
	  			newResidue[newNode] += residue[tempNode]/(double)inSize;
	  			if(!isInArray[newNode]){
	    				isInArray[newNode] = true;
	    				new_candidate_set_ford.push_back(newNode);
					new_pathnum_set.push_back(pathnum_set[j]);	
					vecpos[newNode] = residue_count;
					residue_count+=1;
				}
				else{
					new_pathnum_set[vecpos[newNode]]+=pathnum_set[j];
				}
				
			}
			residue[tempNode] = 0;  
		}
		
		long new_candidateSetSize=new_candidate_set_ford.size();
		for(int j=0;j<new_candidateSetSize;j++){
			int tempNode=new_candidate_set_ford[j];
			residue[tempNode] = newResidue[tempNode];
			newResidue[tempNode] = 0;
			isInArray[tempNode]=false;
			vecpos[tempNode]=0;
		}
		
		forwardNode.push_back(candidate_set_ford);
		
      		candidate_set_ford.swap(new_candidate_set_ford);
		new_candidate_set_ford.clear();
		pathnum_set.swap(new_pathnum_set);
		new_pathnum_set.clear();

		tempLevel++;
		if(tempLevel > maxLevel){
			deterlevel=maxLevel;
			break;
      		}
		if(deter_flag==true){
			break;
		}
    	}		
	
	candidate_set_ford.clear();    	
	pathnum_set.clear();

	return answer;
  }


  //calculate Z(w,q)
  void cal_Zw(int w,int &deterlevel, int maxlevel, long total_walknum, double &pl_value){
   	pl_value=0.0;
	P_lw.clear();//unordered_map <int,unordered_map<int,double> >
	P_ly.clear();//unordered_map <int,unordered_map<int,double> >
  	P_ally.clear();//unordered_map <int,unordered_map<int, vector<pair<int,double> > > >
	Fl.clear();//unordered_map <int,unordered_map<int,double> >
	forwardNode.clear();

  	P_lw=forwardPushByLevel(w, maxlevel, deterlevel, total_walknum, forwardNode);

	if(deterlevel<=1){
		deterlevel=1;
		pl_value=C_value/g.getInSize(w);
		return;
	}
	for(int level=1;level<=deterlevel;level++){
		for(int ilevel=level;ilevel>=1;ilevel--){
			for(int itey=0;itey<forwardNode[ilevel].size();itey++){
				int y=forwardNode[ilevel][itey];
				if(ilevel==level){
					P_ally[level-ilevel][y].push_back(make_pair(y,1));
					P_ly[y][level-ilevel][y]=1;
				}
				else{
					unordered_map<int,int> isinvec;
					int cnt=0;
					for(int k=0;k<g.getInSize(y);k++){
						int nextNode=g.getInVert(y,k);
						for(int veci=0;veci<P_ally[level-ilevel-1][nextNode].size();veci++){
							double prob=P_ally[level-ilevel-1][nextNode][veci].second/g.getInSize(y);
							int finalNode=P_ally[level-ilevel-1][nextNode][veci].first;
							if(isinvec.find(finalNode)==isinvec.end()){
								P_ally[level-ilevel][y].push_back(make_pair(finalNode,prob));
								isinvec[finalNode]=P_ally[level-ilevel][y].size()-1;
								cnt+=1;
								P_ly[finalNode][level-ilevel][y]=prob;
							}
							else{
								int pos=isinvec[finalNode];
								P_ally[level-ilevel][y][pos].second+=prob;
								P_ly[finalNode][level-ilevel][y]+=prob;
							}
						}
					}
					isinvec.clear();
				}
			}
		}
		for(int itex=0;itex<forwardNode[level].size();itex++){
			int x=forwardNode[level][itex];
			double Pr_repeatmeet=0.0;
			for(int jlevel=1;jlevel<=level-1;jlevel++){
				for(int itey=0;itey<forwardNode[jlevel].size();itey++){
					int y=forwardNode[jlevel][itey];
					if(P_ly[x][level-jlevel].find(y)!=P_ly[x][level-jlevel].end()){
						double tmp_ply=P_ly[x][level-jlevel][y];
						Pr_repeatmeet+=pow(C_value,level-jlevel)*tmp_ply*tmp_ply*Fl[jlevel][y];
					}
				}
			}
			Fl[level][x]=pow(C_value,level)*P_lw[level][x]*P_lw[level][x]-Pr_repeatmeet;
			if(Fl[level][x]<0.01*eps){
				Fl[level][x]=0.0;
			}
			pl_value+=Fl[level][x];
		}
		P_ly.clear();
		P_ally.clear();
	}

  	return;
  }
    
  
 
  double sampleD(int nodeId,double ppr_nodeId,long walk_num, int maxlevel, int &deterlevel, double para_r){ 
	double meet = 0;
        double pl=0.0;
	deterlevel=1;    	

	if(walk_num==0.0){
		deterlevel=1;
		return 1-C_value/g.getInSize(nodeId); 
	}
    	else{
		bool first_lflag=false; 
		long firstlevelnum=0;
		int winSize=g.getInSize(nodeId);
		if(winSize*winSize>=para_r*walk_num){
			deterlevel=1;
			first_lflag=true;
		}
		else{
			for(int li=0;li<winSize;li++){
				int tmpN=g.getInVert(nodeId,li);
				firstlevelnum+=g.getInSize(tmpN);
				if(firstlevelnum*firstlevelnum>=para_r*walk_num){
					deterlevel=1;
					first_lflag=true;
					break;
				}
			}
			if(first_lflag==false){
				pl=0.0;
				cal_Zw(nodeId,deterlevel,maxlevel,(long)(para_r*walk_num),pl);
			}
		}
	        for(long i = 0; i < walk_num; i++){
			bool flag=true;	  
			int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
			int count =0;
			if(g.getInSize(nodeId)<2){
				break;
			}
			u_newNode = g.getInVert(nodeId, R.generateRandom() % g.getInSize(nodeId));
			v_newNode = g.getInVert(nodeId, R.generateRandom() % g.getInSize(nodeId));
			for(int l=1;l<deterlevel;l++){
				if (u_newNode ==  v_newNode){
					flag=false;
					break;
				}
				if((g.getInSize(u_newNode)==0)||(g.getInSize(v_newNode)==0)){
					flag=false;
					break;
				}
				int u_tmpNode,v_tmpNode;
				u_tmpNode = g.getInVert(u_newNode, (int)(R.generateRandom() % g.getInSize(u_newNode)));
				v_tmpNode = g.getInVert(v_newNode, (int)(R.generateRandom() % g.getInSize(v_newNode)));
				u_newNode=u_tmpNode;
				v_newNode=v_tmpNode;		
			}	
			if ((u_newNode != v_newNode)&&(flag==true)){
				while(R.drand() < C_value){
	    				int length = g.getInSize(u_newNode);
	    				if(length == 0)
	      					break;
	    				int r = R.generateRandom() % length;
	    				u_nextNode = g.getInVert(u_newNode, r);
	    				length = g.getInSize(v_newNode);
	    				if(length == 0)
	      					break;
	    				r = R.generateRandom() % length;
	    				v_nextNode = g.getInVert(v_newNode, r);
	    				if(u_nextNode == v_nextNode){
	      					meet+=1.0;
						break;
	    				}
	    				u_newNode = u_nextNode;
	    				v_newNode = v_nextNode;
	  			}   
			}
         	}

		meet/=walk_num;

		if(deterlevel<=1){
			return (1-C_value/g.getInSize(nodeId) - C_value * meet);
		}
		else{
			return (1-pl-(pow(C_value,deterlevel))*meet);
		}
	}
  }


  
  void ExactSim(int u, string outputFile, double para_r){
    if(g.getInSize(u)==0){
	ofstream fout;
    	fout.open(outputFile);
    	fout.setf(ios::fixed,ios::floatfield);
    	fout.precision(15);
    	if(!fout) cout<<"Fail to open the writed file"<<endl;
    	fout<<u<<" "<<1.0<<endl;
	fout.close();
    	return;
    }    

    clock_t t_start = clock();
    
    int targetlevel=power_ite_num;//iteration nums
    double* scores = new double[vert];
    double* nextScores = new double[vert];
    double* d = new double[vert];//matrix D
    double** ulevelscores =new double*[targetlevel+1];
    double* simrank = new double[vert];
    double alpha = 1 - sqrtC;//teleport probability in PPR computation
    long nr_pre = power_walk_num;//thesampling number in total
    int deterlevel;//l(k)
    int max_deterlevel=5;//upper bound of l(k)  
    double avg_ppr=0.0;//Average value of each node's PPR
    int ppr_nnz=0;//nonzero of PPR 
    for(int i = 0; i < vert; i++){
      	simrank[i]=0;
      	scores[i] = 0;
      	nextScores[i] = 0;
    }
    

    //PPR computation for source node u
    //
    scores[u]=1.0;
    for(int i=0;i<targetlevel;i++){
        for(int j=0;j<vert;j++){
	    for(int k=0;k<g.getInSize(j);k++){
		int tempNode = g.getInVert(j,k);
		nextScores[tempNode]+=(1-alpha)*scores[j]/g.getInSize(j);	
	    }
	}
	nextScores[u]+=alpha;
	for(int j=0;j<vert;j++){
		scores[j]=nextScores[j];
		nextScores[j]=0.0;
		if((i==targetlevel-1)&&(scores[j]>0.0)&&(g.getInSize(j)>1)){
			avg_ppr+=scores[j];	
			ppr_nnz+=1;
		}
	}
    }
    if(ppr_nnz>0){
    	avg_ppr/=ppr_nnz;       
    } 




    //Matrix D's estimation
    //
    for(int j = 0; j < vert; j++){
      	if(g.getInSize(j) == 0){
		d[j] = 1;
		continue;
      	}
      	else if((scores[j]==0)||(g.getInSize(j)==1)){
		d[j] = 1 - C_value/(double) g.getInSize(j);
      		continue;
      	}
      	else{
      		long nrj_pre = (long) nr_pre*scores[j]/(avg_ppr*ppr_nnz);//the sample number from node vj
		d[j]=sampleD(j, scores[j], nrj_pre, max_deterlevel, deterlevel, para_r); 	
      	}
    }




    //SimRank vectors' linearization
    //    
    for(int i=0;i<vert;i++){
	nextScores[i]=0.0;
    }
    nextScores[u]=1.0;
    for(int level=0;level<=targetlevel;level++){
	ulevelscores[level] = new double[vert];
	for(int i=0;i<vert;i++){
		ulevelscores[level][i]=d[i]*nextScores[i];
		scores[i]=0.0;
	}
	for(int j=0;j<vert;j++){
		for(int k=0;k<g.getOutSize(j);k++){
			int tempNode=g.getOutVert(j,k);
			if(g.getInSize(tempNode)>0){
				scores[j]+=nextScores[tempNode]/g.getInSize(tempNode);
			}
		}
	}
	for(int j=0;j<vert;j++){
		nextScores[j]=scores[j];
	}
    }    

    for(int j=0;j<vert;j++){
	simrank[j]=ulevelscores[targetlevel][j];
    }
    
    for(int level=1;level<=targetlevel;level++){
	for(int i = 0; i < vert; i++){
		nextScores[i] = 0;
      		scores[i]=0;
	}
	for(int j=0;j<vert;j++){
		for(int k=0;k<g.getInSize(j);k++){
			int tempNode=g.getInVert(j,k);
			if(g.getInSize(j)>0){
				nextScores[j]+=simrank[tempNode]/g.getInSize(j);
			}
		}
		int tmplevel=targetlevel-level;
		nextScores[j]=nextScores[j]*C_value;
		scores[j]=nextScores[j]+ulevelscores[tmplevel][j];
	}
	for(int i=0;i<vert;i++){
		simrank[i]=scores[i];
	}
    }
    
    clock_t t_end =clock();
    double querytime=(t_end-t_start)/(double) CLOCKS_PER_SEC;
    cout<<"querytime = "<<querytime<<endl;
    avg_time+=querytime;




    //Output results
    //
    cout<<"Outfile address:"<<outputFile<<endl;
    ofstream fout;
    fout.open(outputFile);
    fout.setf(ios::fixed,ios::floatfield);
    fout.precision(15);
    if(!fout) cout<<"Fail to open the writed file"<<endl;
    for(int j = 0; j < vert; j++){
        if(j == u){
	    simrank[j]=1.0;
	    fout<<j<<" "<<simrank[j]<<endl;
	}
	else if(simrank[j]>0.0){
		fout<<j<<" "<<simrank[j]<<endl;
	}
    }
    fout.close();
    

    delete[] scores;
    delete[] nextScores;   
    delete[] d;
    delete[] simrank;
    for(int i=0;i<=targetlevel;i++){
	delete[] ulevelscores[i]; 
    }
  }
  


};




#endif
