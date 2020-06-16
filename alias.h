#ifndef ALIAS_H
#define ALIAS_H

#include <random>
#include <algorithm>
#include <stack>

using namespace std;

class Alias{
public:
	double* p;
	int* h;
	pair<int, int>* map1;
	int n;
	Alias(vector<pair<pair<int, int>, double> > pi){
		double sum = 0;
		n = pi.size();
		stack<int> small;
		stack<int> big;
		p = new double[n];
		h = new int[n];
		map1 = new pair<int, int>[n];
		for(int i = 0; i < n; i++){
			sum += pi[i].second;
			map1[i] = pi[i].first;
		}
		for(int i = 0; i < n; i++){
			p[i] = pi[i].second * n / sum;
			if(p[i] > 1)
				big.push(i);
			else
				small.push(i);
		}
		while(!small.empty() && !big.empty()){
			int smallId = small.top();
			small.pop();
			int bigId = big.top();
			h[smallId] = bigId;
			p[bigId] -= (1-p[smallId]);
			if(p[bigId] < 1){
				small.push(bigId);
				big.pop();
			}
		}
	}

	~Alias(){
		delete[] p;
		delete[] h;
		delete[] map1;
	}
	pair<int, int> generateRandom(Random& R){
		int firstId = R.drand() * n;
		pair<int, int> answer = R.drand() < p[firstId] ? map1[firstId] : map1[h[firstId]];
		return answer;
	}
};

#endif