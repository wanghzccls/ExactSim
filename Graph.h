#ifndef GRAPH_H
#define GRAPH_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;

class Graph{
public:
	int n;	//number of nodes
	int m;	//number of edges
	int** inAdjList;
	int** outAdjList;
	int* indegree;
	int* outdegree;

	Graph(){
		
	}

	~Graph(){

	}

	void inputGraph(string filename){
		m=0;
		ifstream infile(filename.c_str());
		infile >> n;
		cout<<"Nodes n="<<n<<endl;

		indegree=new int[n];
		outdegree=new int[n];
		for(int i=0;i<n;i++){
			indegree[i]=0;
			outdegree[i]=0;
		}
		//read graph and get degree info
		int from;
		int to;
		//while(infile.good())
		while(infile>>from>>to){
			//infile >> from >> to;
			//cout << "from:" << from << " to: " << to << endl;
			outdegree[from]++;
			indegree[to]++;
		}

		inAdjList=new int*[n];
		outAdjList=new int*[n];
		for(int i=0;i<n;i++){
			inAdjList[i]=new int[indegree[i]];
			outAdjList[i]=new int[outdegree[i]];
			//cout << i << "indegree: " << indegree[i] << endl;
		}
		int* pointer_in=new int[n];
		int* pointer_out=new int[n];
		for(int i=0;i<n;i++){
			pointer_in[i]=0;
			pointer_out[i]=0;
		}
		infile.clear();
		infile.seekg(0);

		clock_t t1=clock();
		infile >> n;
		while(infile>>from>>to){
			outAdjList[from][pointer_out[from]]=to;
			pointer_out[from]++;
			inAdjList[to][pointer_in[to]]=from;
			pointer_in[to]++;

			m++;
		}
		infile.close();
		clock_t t2=clock();
		cout<<"Edges m="<<m<<endl;
		cout<<"reading in graph takes "<<(t2-t1)/(1.0*CLOCKS_PER_SEC)<<" s."<<endl;
 		cout<<"====Graphs reading done!====\n"<<endl;

		delete[] pointer_in;
		delete[] pointer_out;
	}
	int getInSize(int vert){
		return indegree[vert];
	}
	int getInVert(int vert, int pos){
		return inAdjList[vert][pos];
	}
	int getOutSize(int vert){
		return outdegree[vert];
	}
	int getOutVert(int vert, int pos){
		return outAdjList[vert][pos];
	}

};


#endif
