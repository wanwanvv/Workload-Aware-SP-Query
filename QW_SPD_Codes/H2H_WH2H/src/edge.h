#pragma once
using namespace std;
class edge
{
public:
	edge() {};
	edge(const int& u, const int& v, const int& weight) {
		this->u = u;
		this->v = v;
		this->weight = weight;
	}
	int u;
	int v;
	int weight;
};

