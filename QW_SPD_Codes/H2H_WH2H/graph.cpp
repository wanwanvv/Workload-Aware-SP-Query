#include "src/graph.h"

graph::graph(const char* graph_file, int maxListSize, float gamma, int blocksize, int maxDegree, bool directed) {
	ifstream ifs(graph_file);
	ifs >> numOfVertices >> numOfEdges;//read number of vertices and edges from graph file

	edge temp;
	int u{ 0 }, v{ 0 };
	int w;
	
	if (directed) {
		es_directed.resize(numOfVertices);
		//es_dijkstra.resize(numOfVertices);
		degree.assign(numOfVertices, make_pair(0, 0));
		//neighbors.resize(numOfVertices);
		//read graph data
		for (; ifs >> u >> v >> w;) {
			temp.u = u;
			temp.v = v;
			temp.weight = w;
			degree[u].second++;
			degree[v].first++;
			//neighbors[u].push_back(v);
			//neighbors[v].push_back(u);
			es_directed[u].first.push_back(temp);
			es_directed[v].second.push_back(edge(u, v, w));
			//es_dijkstra[u].first.push_back(temp);
			//es_dijkstra[v].second.push_back(edge(u, v, w));
		}
		for (int i = 0; i < numOfVertices; i++) {
			//neighbors[i].sort();
			//neighbors[i].unique();
			es_directed[i].first.sort([](edge e1, edge e2) {return e1.v < e2.v; });
			es_directed[i].second.sort([](edge e1, edge e2) {return e1.u < e2.u; });
			//es_dijkstra[i].first.sort([](edge e1, edge e2) {return e1.v < e2.v; });
			//es_dijkstra[i].second.sort([](edge e1, edge e2) {return e1.u < e2.u; });
		}
		/*dis_dijkstra.resize(numOfVertices);
		for (int i = 0; i < numOfVertices; i++) {
			dis_dijkstra[i].resize(numOfVertices);
			dijkstra(i, 0);
		}*/
		/*for (int i = 0; i < numOfVertices; i++) {
			for (int j = 0; j < numOfVertices; j++) {
				cout << dis_dijkstra[i][j] << " ";
			}
			cout << endl;
		}*/
	}
	else {
		es.resize(numOfVertices);
		//read graph data
		for (; ifs >> u >> v >> w;) {
			temp.u = u;
			temp.v = v;
			temp.weight = w;
			es[u].push_back(temp);
		}
		for (int i = 0; i < numOfVertices; i++) {
			es[i].sort([](edge e1, edge e2) {return e1.v < e2.v; });
		}
	}

	ifs.close();

	//initialize variables
	MAXLISTSIZE = maxListSize;
	GAMMA = gamma; //Make the weight of degree as 1.
	BLOCKSIZE = blocksize;
	MAXDEGREE = maxDegree;
	isRemoved.assign(numOfVertices, 0);
	level.assign(numOfVertices, -1);
	nodes.resize(numOfVertices);
	vertexListIterator.assign(MAXLISTSIZE, NULL);
	vertexIterator.assign(numOfVertices, NULL);
}

void graph::readHFVertexData(const char* HFVertex_file) {
	ifstream ifsHF(HFVertex_file);
	assert(ifsHF.is_open());

	HFVertex.assign(numOfVertices, make_pair(0, 0));
	allHFVertex.reserve(numOfVertices);

	//Read data of high frequency vertices.
	int frq;
	int u;
	for (int i = 0, j = 0; ifsHF >> u >> frq; ) {
		HFVertex[u] = make_pair(u, frq);
		allHFVertex.push_back(u);
	}
	ifsHF.close();

	//Set value of other vertices as 0.
	for (int i = 0; i < numOfVertices; i++) {
		if (HFVertex[i].second == 0) {
			HFVertex[i].first = i;
		}
	}

	//block division
	int maxAvgValue = 0;
	int curValue = 0;
	for (int i = 0; i < allHFVertex.size() && i < BLOCKSIZE; i++) {
		maxAvgValue += HFVertex[allHFVertex[i]].second;
	}
	maxAvgValue /= BLOCKSIZE;
	for (int i = 0; i < allHFVertex.size(); i++)
	{
		if (i % BLOCKSIZE == 0) {
			curValue = 0;
			for (int j = i; j < allHFVertex.size() && j < i + BLOCKSIZE; j++) {
				curValue += HFVertex[allHFVertex[j]].second;
			}
			curValue = curValue / (float)BLOCKSIZE / (float)maxAvgValue * MAXDEGREE;
		}
		HFVertex[allHFVertex[i]].second = curValue;
	}
}

void graph::addVertexToVertexList_readData_HFVertices(const int& nodeNum, const int& value) {
	//If the iterator of current vertex is not NULL, delete current vertex.
	if (vertexIterator[nodeNum] && vertexListIterator[value]) {
		(**vertexListIterator[value])->erase(*vertexIterator[nodeNum]);//Erase current vertex from the list with ID 'value'.
		if ((**vertexListIterator[value])->size() == 0) {//If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[value]);
			delete vertexListIterator[value];
			vertexListIterator[value] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curValue = calValue(nodeNum);
	if (vertexListIterator[curValue]) {//If the list with ID 'value' is not empty, push the current vertex into it.
		(**vertexListIterator[curValue])->push_back(make_pair(nodeNum, curValue));
		list<pair<int, int>>::iterator it = (**vertexListIterator[curValue])->end();
		it--;
		*vertexIterator[nodeNum] = it;
	}
	else {//Create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curValue] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, curValue));
		int curMaxValue = curValue - 1;	//curMaxValue is the maximum value of current inverted List.
		while (curMaxValue >= 0 && !vertexListIterator[curMaxValue]) {//Get the maximum value.
			curMaxValue--;
		}
		if (curMaxValue > 0) {
			list<list<pair<int, int>>*>::iterator it = vertexList.begin();
			for (; it != vertexList.end() && (*it)->front().first != (**vertexListIterator[curMaxValue])->begin()->first; it++);
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curValue] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*(vertexListIterator[curValue]) = vertexList.begin();
		}
		auto it = tempList->end();
		*(vertexIterator[nodeNum]) = --it;
	}
}

void graph::addVertexToVertexList_readData_HFVertices_directed(const int& nodeNum, const int& value) {
	//If the iterator of current vertex is not NULL, delete current vertex.
	if (vertexIterator[nodeNum] && vertexListIterator[value]) {
		(**vertexListIterator[value])->erase(*vertexIterator[nodeNum]);//Erase current vertex from the list with ID 'value'.
		if ((**vertexListIterator[value])->size() == 0) {//If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[value]);
			delete vertexListIterator[value];
			vertexListIterator[value] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curValue = calValue_directed(nodeNum);
	if (vertexListIterator[curValue]) {//If the list with ID 'value' is not empty, push the current vertex into it.
		(**vertexListIterator[curValue])->push_back(make_pair(nodeNum, curValue));
		list<pair<int, int>>::iterator it = (**vertexListIterator[curValue])->end();
		it--;
		*vertexIterator[nodeNum] = it;
	}
	else {//Create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curValue] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, curValue));
		int curMaxValue = curValue - 1;	//curMaxValue is the maximum value of current inverted List.
		while (curMaxValue >= 0 && !vertexListIterator[curMaxValue]) {//Get the maximum value.
			curMaxValue--;
		}
		if (curMaxValue > 0) {
			list<list<pair<int, int>>*>::iterator it = vertexList.begin();
			for (; it != vertexList.end() && (*it)->front().first != (**vertexListIterator[curMaxValue])->begin()->first; it++);
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curValue] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*(vertexListIterator[curValue]) = vertexList.begin();
		}
		auto it = tempList->end();
		*(vertexIterator[nodeNum]) = --it;
	}
}

void graph::addVertexToVertexList_readData(const int& nodeNum, const int& degree) {
	if (vertexIterator[nodeNum]) {//If the iterator of current vertex is not NULL, delete current vertex.
		(**vertexListIterator[degree])->erase(*vertexIterator[nodeNum]);//Erase current vertex from the list with ID 'degree'.
		if ((**vertexListIterator[degree])->size() == 0) { //If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[degree]);
			delete vertexListIterator[degree];
			vertexListIterator[degree] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curDegree = es[nodeNum].size();
	if (!vertexListIterator[curDegree]) {//If the list with ID 'degree' is empty, create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curDegree] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, 0));
		int curMaxDegree = curDegree - 1;//curMaxDegree is the maximum degree of current inverted List.
		while (curMaxDegree >= 0 && !vertexListIterator[curMaxDegree]) {//Get the maximum degree.
			curMaxDegree--;
		}
		if (curMaxDegree > 0) {
			list<list<pair<int, int>>*>::iterator it = vertexList.begin();
			for (; (*it)->front().first != (**vertexListIterator[curMaxDegree])->begin()->first && it != vertexList.end(); it++);
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curDegree] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*vertexListIterator[curDegree] = vertexList.begin();
		}
		*(vertexIterator[nodeNum]) = --tempList->end();
		(**vertexListIterator[curDegree]) = tempList;//Save the iterator of the new list.
	}
	else {//If the list with ID 'degree' is not empty, push the current vertex into it.
		if (sequence == 0) {//If the process is to read data, push vertex into inverted list by ID directly.
			(**vertexListIterator[curDegree])->push_back(make_pair(nodeNum, 0));
			*vertexIterator[nodeNum] = --(**vertexListIterator[curDegree])->end();
		}
		else {
			(**vertexListIterator[curDegree])->push_back(make_pair(nodeNum, 0));
			list<pair<int, int>>::iterator it = (**vertexListIterator[curDegree])->end();
			it--;
			*vertexIterator[nodeNum] = it;
		}
	}
}

void graph::addVertexToVertexList_readData_directed(const int& nodeNum, const int& degree) {
	if (vertexIterator[nodeNum]) {//If the iterator of current vertex is not NULL, delete current vertex.
		(**vertexListIterator[degree])->erase(*vertexIterator[nodeNum]);//Erase current vertex from the list with ID 'degree'.
		if ((**vertexListIterator[degree])->size() == 0) { //If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[degree]);
			delete vertexListIterator[degree];
			vertexListIterator[degree] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curDegree = calDegree(nodeNum);
	if (!vertexListIterator[curDegree]) {//If the list with ID 'degree' is empty, create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curDegree] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, 0));
		int curMaxDegree = curDegree - 1;//curMaxDegree is the maximum degree of current inverted List.
		while (curMaxDegree >= 0 && !vertexListIterator[curMaxDegree]) {//Get the maximum degree.
			curMaxDegree--;
		}
		if (curMaxDegree > 0) {
			list<list<pair<int, int>>*>::iterator it = vertexList.begin();
			for (; (*it)->front().first != (**vertexListIterator[curMaxDegree])->begin()->first && it != vertexList.end(); it++);
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curDegree] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*vertexListIterator[curDegree] = vertexList.begin();
		}
		*(vertexIterator[nodeNum]) = --tempList->end();
		(**vertexListIterator[curDegree]) = tempList;//Save the iterator of the new list.
	}
	else {//If the list with ID 'degree' is not empty, push the current vertex into it.
		if (sequence == 0) {//If the process is to read data, push vertex into inverted list by ID directly.
			(**vertexListIterator[curDegree])->push_back(make_pair(nodeNum, 0));
			*vertexIterator[nodeNum] = --(**vertexListIterator[curDegree])->end();
		}
		else {
			(**vertexListIterator[curDegree])->push_back(make_pair(nodeNum, 0));
			list<pair<int, int>>::iterator it = (**vertexListIterator[curDegree])->end();
			it--;
			*vertexIterator[nodeNum] = it;
		}
	}
}

void graph::addVertexToVertexList(const int& nodeNum, const int degree) {
	if (vertexIterator[nodeNum] && vertexListIterator[degree]) {//If the iterator of current vertex is not NULL, delete current vertex.
		(**vertexListIterator[degree])->erase(*vertexIterator[nodeNum]); //Erase current vertex from the list with ID 'degree'.
		if ((**vertexListIterator[degree])->size() == 0) {//If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[degree]);
			delete vertexListIterator[degree];
			vertexListIterator[degree] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curDegree = es[nodeNum].size();
	if (vertexListIterator[curDegree]) {//If the list with ID 'degree' is not empty, push the current vertex into it.
		(**vertexListIterator[curDegree])->push_back(make_pair(nodeNum, 0));
		list<pair<int, int>>::iterator it = (**vertexListIterator[curDegree])->end();
		it--;
		*vertexIterator[nodeNum] = it;
	}
	else {//If the list with ID 'degree' is empty, create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curDegree] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, 0));
		int curMaxDegree = curDegree - 1;//curMaxDegree is the maximum degree of current inverted List.
		while (curMaxDegree >= 0 && !vertexListIterator[curMaxDegree]) {//Get the maximum degree.
			curMaxDegree--;
		}
		if (curMaxDegree > 0) {
			list<list<pair<int, int>>*>::iterator it = *vertexListIterator[curMaxDegree];
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curDegree] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*vertexListIterator[curDegree] = vertexList.begin();
		}
		*(vertexIterator[nodeNum]) = --tempList->end();
		(**vertexListIterator[curDegree]) = tempList;//Save the iterator of the new list.
	}
}

void graph::addVertexToVertexList_directed(const int& nodeNum, const int degree) {
	if (vertexIterator[nodeNum] && vertexListIterator[degree]) {//If the iterator of current vertex is not NULL, delete current vertex.
		(**vertexListIterator[degree])->erase(*vertexIterator[nodeNum]); //Erase current vertex from the list with ID 'degree'.
		if ((**vertexListIterator[degree])->size() == 0) {//If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[degree]);
			delete vertexListIterator[degree];
			vertexListIterator[degree] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curDegree = calDegree(nodeNum);
	if (vertexListIterator[curDegree]) {//If the list with ID 'degree' is not empty, push the current vertex into it.
		(**vertexListIterator[curDegree])->push_back(make_pair(nodeNum, 0));
		list<pair<int, int>>::iterator it = (**vertexListIterator[curDegree])->end();
		it--;
		*vertexIterator[nodeNum] = it;
	}
	else {//If the list with ID 'degree' is empty, create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curDegree] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, 0));
		int curMaxDegree = curDegree - 1;//curMaxDegree is the maximum degree of current inverted List.
		while (curMaxDegree >= 0 && !vertexListIterator[curMaxDegree]) {//Get the maximum degree.
			curMaxDegree--;
		}
		if (curMaxDegree > 0) {
			list<list<pair<int, int>>*>::iterator it = *vertexListIterator[curMaxDegree];
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curDegree] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*vertexListIterator[curDegree] = vertexList.begin();
		}
		*(vertexIterator[nodeNum]) = --tempList->end();
		(**vertexListIterator[curDegree]) = tempList;//Save the iterator of the new list.
	}
}

void graph::addVertexToVertexList_HFVertices(const int& nodeNum, const int value) {
	if (vertexIterator[nodeNum] && vertexListIterator[value]) {//If the iterator of current vertex is not NULL, delete current vertex.
		(**vertexListIterator[value])->erase(*vertexIterator[nodeNum]);//Erase current vertex from the list with ID 'value'.
		if ((**vertexListIterator[value])->size() == 0) {//If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[value]);
			delete vertexListIterator[value];
			vertexListIterator[value] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curValue = calValue(nodeNum);
	if (vertexListIterator[curValue]) {//If the list with ID 'value' is not empty, push the current vertex into it.
		(**vertexListIterator[curValue])->push_back(make_pair(nodeNum, curValue));
		list<pair<int, int>>::iterator it = (**vertexListIterator[curValue])->end();
		it--;
		*vertexIterator[nodeNum] = it;
	}
	else {//If the list with ID 'value' is empty, create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curValue] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, curValue));
		int curMaxValue = curValue - 1;//curMaxValue is the maximum value of current inverted List.
		while (curMaxValue >= 0 && !vertexListIterator[curMaxValue]) {//Get the maximum value.
			curMaxValue--;
		}
		if (curMaxValue > 0) {
			list<list<pair<int, int>>*>::iterator it = vertexList.begin();
			for (; it != vertexList.end() && (*it)->front().first != (**vertexListIterator[curMaxValue])->begin()->first; it++);
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curValue] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*(vertexListIterator[curValue]) = vertexList.begin();
		}
		auto it = tempList->end();
		*(vertexIterator[nodeNum]) = --it;
	}
}

void graph::addVertexToVertexList_HFVertices_directed(const int& nodeNum, const int value) {
	if (vertexIterator[nodeNum] && vertexListIterator[value]) {//If the iterator of current vertex is not NULL, delete current vertex.
		(**vertexListIterator[value])->erase(*vertexIterator[nodeNum]);//Erase current vertex from the list with ID 'value'.
		if ((**vertexListIterator[value])->size() == 0) {//If the list becomes empty, erase the list from the first dimension.
			vertexList.erase(*vertexListIterator[value]);
			delete vertexListIterator[value];
			vertexListIterator[value] = NULL;
		}
		delete vertexIterator[nodeNum];
		vertexIterator[nodeNum] = NULL;
	}
	vertexIterator[nodeNum] = new list<pair<int, int>>::iterator;
	int curValue = calValue_directed(nodeNum);
	if (vertexListIterator[curValue]) {//If the list with ID 'value' is not empty, push the current vertex into it.
		(**vertexListIterator[curValue])->push_back(make_pair(nodeNum, curValue));
		list<pair<int, int>>::iterator it = (**vertexListIterator[curValue])->end();
		it--;
		*vertexIterator[nodeNum] = it;
	}
	else {//If the list with ID 'value' is empty, create a new list, and then push the current vertex into it.
		list<pair<int, int>>* tempList = new list<pair<int, int>>;
		vertexListIterator[curValue] = new list<list<pair<int, int>>*>::iterator;
		tempList->push_back(make_pair(nodeNum, curValue));
		int curMaxValue = curValue - 1;//curMaxValue is the maximum value of current inverted List.
		while (curMaxValue >= 0 && !vertexListIterator[curMaxValue]) {//Get the maximum value.
			curMaxValue--;
		}
		if (curMaxValue > 0) {
			list<list<pair<int, int>>*>::iterator it = vertexList.begin();
			for (; it != vertexList.end() && (*it)->front().first != (**vertexListIterator[curMaxValue])->begin()->first; it++);
			vertexList.insert(++it, tempList);//Insert the new list into the correct position to make sure that the inverted list is in order.
			*vertexListIterator[curValue] = --it;//Save the iterator of the new list.
		}
		else {//Push the new list at the beginning of the inverted list.
			vertexList.push_front(tempList);
			*(vertexListIterator[curValue]) = vertexList.begin();
		}
		auto it = tempList->end();
		*(vertexIterator[nodeNum]) = --it;
	}
}

void graph::removeDegreeOne_HFVertices(const int& number) {
	int curNum = number;
	edge temp;
	temp.u = es[curNum].front().u;
	temp.v = es[curNum].front().v;
	temp.weight = es[curNum].front().weight;
	int v_originalDegree = es[temp.v].size();
	isRemoved[curNum] = 1;
	nodes[curNum].curNum = curNum;//Record the ID of current removing vertex.
	nodes[curNum].neighbor.push_back(temp.v);
	nodes[curNum].neighborDis.push_back(temp.weight);//Save distance to the neighbor.
	vertexList.front()->erase(*vertexIterator[curNum]);
	nodes[curNum].level = sequence;//Set sequence.
	level[curNum] = sequence++;
	curNum = temp.v;
	//Process edge (u,v).
	for (auto it = es[curNum].begin(); it != es[curNum].end(); it++) {
		if ((*it).v == temp.u) {
			es[curNum].erase(it);
			break;
		}
	}

	//update vertex v.
	int curMinValue = 0;
	while (!vertexListIterator[curMinValue] || (**(vertexListIterator[curMinValue]))->size() == 0) {//Get the minimum value from inverted list.
		curMinValue++;
	}

	if (es[curNum].size() == 1) {//If the degree of vertex is one ,insert into list directly.
		isAltered[curNum] = 0;
		addVertexToVertexList_HFVertices(curNum, (*vertexIterator[curNum])->second);
	}
	else {
		if (!isAltered[curNum]) {
			if (calValue(curNum) >= curMinValue) {
				isAltered[curNum] = 1;
				changedNeighbors.push(make_pair(curNum, v_originalDegree));
			}
			else {//If the degree of the vertex is smaller than the minimum value, update the list.
				if (vertexList.front()->size() == 0) {
					vertexList.pop_front();//Delete the list with value 1.
					delete vertexListIterator[1];
					vertexListIterator[1] = NULL;
				}
				addVertexToVertexList_HFVertices(curNum, (*vertexIterator[curNum])->second);
			}
		}
	}
}

void graph::removeDegreeOne_HFVertices_directed(const int& number) {

	int curNum = number;
	nodes[curNum].neighbor_in.reserve(1);
	nodes[curNum].neighbor_out.reserve(1);
	nodes[curNum].neighborDis_in.reserve(1);
	nodes[curNum].neighborDis_out.reserve(1);

	int v = es_directed[curNum].first.front().v;
	//Record the original number of neighbors of the neighbor vertex.
	int originalDegree = calDegree(v);

	//Remove edge (v, curNum) in es_directed[v] and record the neighbors in nodes..
	for (list<edge>::iterator it = es_directed[v].first.begin(); it != es_directed[v].first.end(); ) { //find the edge (v,curNum) and (curNum,v)
		if (it->v == curNum) {
			nodes[curNum].neighborDis_in.push_back(it->weight);//Save distance from the neighbor.
			nodes[curNum].neighbor_in.push_back(it->u);
			es_directed[v].first.erase(it++);
			degree[v].second--;//update the degree
			break;
		}
		else if (it->v < curNum) {
			it++;
		}
		else {
			break;//There is no edge (v, curNum) in es_directed[v].first.
		}
	}

	//Remove edge (curNum, v) in es_directed[v].second and record the neighbors in nodes..
	for (list<edge>::iterator it = es_directed[v].second.begin(); it != es_directed[v].second.end(); ) { //find the edge (v,curNum) and (curNum,v)
		if (it->u == curNum) {
			nodes[curNum].neighborDis_out.push_back(it->weight);//Save distance from the neighbor.
			nodes[curNum].neighbor_out.push_back(it->v);
			es_directed[v].second.erase(it++);
			degree[v].first--;//update the degree
			break;
		}
		else if (it->u < curNum) {
			it++;
		}
		else {
			break;//There is no edge (curNum, v) in es_directed[v].second.
		}
	}

	isRemoved[curNum] = 1;
	nodes[curNum].curNum = curNum;//Record the ID of current removing vertex.
	vertexList.front()->erase(*vertexIterator[curNum]);
	nodes[curNum].level = sequence;//Set sequence.
	level[curNum] = sequence++;

	//Update neighbors.
	curNum = es_directed[curNum].first.front().v;

	int curValue = calValue_directed(curNum);
	//update the neighbor vertex.
	int curMinValue = 0;
	while (!vertexListIterator[curMinValue] || (**(vertexListIterator[curMinValue]))->size() == 0) {//Get the minimum value from inverted list.
		curMinValue++;
	}
	if (!isAltered[curNum]) {
		if (curValue != 1 && curValue >= curMinValue) {
			isAltered[curNum] = 1;
			changedNeighbors.push(make_pair(curNum, originalDegree));
		}
		else {//If the value of the vertex is smaller than the minimum value, update the list.
			addVertexToVertexList_HFVertices_directed(curNum, (*vertexIterator[curNum])->second);
			if (curValue == 1) {//If value of the neighbor becomes one, process the neighbor directly.
				removeDegreeOne_HFVertices_directed(curNum);
			}
		}
	}
}

void graph::removeDegreeOne(const int& number) {
	int curNum = number;
	edge temp;
	temp.u = es[curNum].front().u;
	temp.v = es[curNum].front().v;
	temp.weight = es[curNum].front().weight;
	int v_originalDegree = es[temp.v].size();
	isRemoved[curNum] = 1;
	nodes[curNum].curNum = curNum;//Record the ID of current removing vertex.
	nodes[curNum].neighbor.push_back(temp.v);
	nodes[curNum].neighborDis.push_back(temp.weight);//Save distances to all the neighbor.
	vertexList.front()->erase(*vertexIterator[curNum]);
	nodes[curNum].level = sequence;//Set sequence.
	level[curNum] = sequence++;
	curNum = temp.v;
	//Process edge (u,v).
	for (list<edge>::iterator it = es[curNum].begin(); it != es[curNum].end(); it++) {
		if (it->v == temp.u) {
			es[curNum].erase(it);
			break;
		}
	}

	//update vertex v.
	int curMinDegree = 1;
	for (auto it = vertexList.begin(); it != vertexList.end(); it++) {
		if ((*it)->size() != 0) {
			curMinDegree = es[(*it)->front().first].size();
			break;
		}
	}
	if (!isAltered[curNum]) {
		if (es[curNum].size() != 1 && es[curNum].size() >= curMinDegree) {
			isAltered[curNum] = 1;
			changedNeighbors.push(make_pair(curNum, v_originalDegree));
		}
		else {//If the degree of the vertex is smaller than the minimum degree, update the list.
			addVertexToVertexList(curNum, v_originalDegree);
			if (es[curNum].size() == 1) {//If degree of v becomes one, process v directly.
				removeDegreeOne(curNum);
			}
		}
	}
}

void graph::removeDegreeOne_directed(const int& number) {
	
	int curNum = number;
	nodes[curNum].neighbor_in.reserve(1);
	nodes[curNum].neighbor_out.reserve(1);
	nodes[curNum].neighborDis_in.reserve(1);
	nodes[curNum].neighborDis_out.reserve(1);

	int v = es_directed[curNum].first.front().v;
	//Record the original number of neighbors of the neighbor vertex.
	int originalDegree = calDegree(v);

	//Remove edge (v, curNum) in es_directed[v] and record the neighbors in nodes..
	for (list<edge>::iterator it = es_directed[v].first.begin(); it != es_directed[v].first.end(); ) { //find the edge (v,curNum) and (curNum,v)
		if (it->v == curNum) {
			nodes[curNum].neighborDis_in.push_back(it->weight);//Save distance from the neighbor.
			nodes[curNum].neighbor_in.push_back(it->u);
			es_directed[v].first.erase(it++);
			degree[v].second--;//update the degree
			break;
		}
		else if (it->v < curNum) {
			it++;
		}
		else {
			break;//There is no edge (v, curNum) in es_directed[v].first.
		}
	}

	//Remove edge (curNum, v) in es_directed[v].second and record the neighbors in nodes..
	for (list<edge>::iterator it = es_directed[v].second.begin(); it != es_directed[v].second.end(); ) { //find the edge (v,curNum) and (curNum,v)
		if (it->u == curNum) {
			nodes[curNum].neighborDis_out.push_back(it->weight);//Save distance from the neighbor.
			nodes[curNum].neighbor_out.push_back(it->v);
			es_directed[v].second.erase(it++);
			degree[v].first--;//update the degree
			break;
		}
		else if (it->u < curNum) {
			it++;
		}
		else {
			break;//There is no edge (curNum, v) in es_directed[v].second.
		}
	}

	isRemoved[curNum] = 1;
	nodes[curNum].curNum = curNum;//Record the ID of current removing vertex.
	vertexList.front()->erase(*vertexIterator[curNum]);
	nodes[curNum].level = sequence;//Set sequence.
	level[curNum] = sequence++;
	
	//Update neighbors.
	curNum = es_directed[curNum].first.front().v;

	int curDegree = calDegree(curNum);
	//update the neighbor vertex.
	int curMinDegree = 0;
	while (!vertexListIterator[curMinDegree] || (**(vertexListIterator[curMinDegree]))->size() == 0) {//Get the minimum value from inverted list.
		curMinDegree++;
	}
	
	if (!isAltered[curNum]) {
		if (curDegree != 1 && curDegree >= curMinDegree) {
			isAltered[curNum] = 1;
			changedNeighbors.push(make_pair(curNum, originalDegree));
		}
		else {//If the degree of the vertex is smaller than the minimum degree, update the list.
			addVertexToVertexList_directed(curNum, originalDegree);
			if (curDegree == 1) {//If degree of the neighbor becomes one, process the neighbor directly.
				removeDegreeOne_directed(curNum);
			}
		}
	}
}

void graph::removeByDegree(const int& number, const int& degree) {
	int edgeNum = degree * (degree - 1) / 2;//Calculate the number of edges to be updated.
	vector<edge> edgeProcess;//Record edges after removing the vertex.
	edgeProcess.reserve(degree * (degree - 1));
	vector<pair<int, int>> originalEdge(degree); //Record edges before removing the vertex.
	vector<int> originalDegree(degree);//Record the degree of neighbors before removing the vertex.
	edge temp;

	//Process edges to be added into the graph after removing the vertex.
	nodes[number].curNum = number;//Record the ID of current removed vertex.
	nodes[number].neighbor.reserve(degree);
	nodes[number].neighborDis.reserve(degree);
	int j = 0;
	for (list<edge>::iterator it = es[number].begin(); it != es[number].end(); it++, j++) {
		nodes[number].neighbor.push_back(it->v);//Record neighbors in nodes.
		nodes[number].neighborDis.push_back(it->weight);//Record distance to neighbors in nodes.
		originalEdge[j].first = it->v;
		originalEdge[j].second = it->weight;
		originalDegree[j] = es[it->v].size();
	}

	//Calculate the edge weight of added edges.
	for (int i = 0; i < degree; i++) {
		for (int j = 0; j < degree; j++) {
			if (j != i) {
				temp.u = originalEdge[i].first;
				temp.v = originalEdge[j].first;
				temp.weight = originalEdge[i].second + originalEdge[j].second;
				edgeProcess.push_back(temp);
			}
		}
	}

	//Remove current vertex and corresponding edges.
	isRemoved[number] = 1;
	nodes[number].level = sequence;//Set sequence.
	level[number] = sequence++;//Set sequence++.
	int calEdgeNum = 0;
	//Add edges into es[u] which is similar with merge-sort.
	for (int i = 0; i < degree; i++) {
		int u = originalEdge[i].first;
		flag = 0;
		list<edge>::iterator it = es[u].begin();
		for (; it != es[u].end() && calEdgeNum < (i + 1) * (degree - 1);) { //Find edge (u,v).
			if (it->v == number) {//Delete edge £¨u,curNum£©.
				es[u].erase(it++);
				flag = 1;
				continue;
			}
			if (it->v == edgeProcess[calEdgeNum].v) {//If the edge is the same as an edge of vertex u, compare the edge weight.
				if (it->weight > edgeProcess[calEdgeNum].weight) {
					it->weight = edgeProcess[calEdgeNum].weight;
				}
				it++;
				calEdgeNum++;
			}
			else if (it->v < edgeProcess[calEdgeNum].v) {
				it++;
			}
			else {
				es[u].insert(it, edgeProcess[calEdgeNum++]);
			}
		}
		while (flag == 0) {//guarantee all the edges have been deleted.
			if (it->v == number) {//Delete edge £¨u,curNum£©.
				es[u].erase(it);
				break;
			}
			it++;
		}
		while (calEdgeNum < (i + 1) * (degree - 1)) {//Add other edges into es[u].
			es[u].push_back(edgeProcess[calEdgeNum++]);
		}
	}

	int curMaxDegree = es[vertexList.front()->front().first].size();//Get the minimum degree from inverted list.
	//Update vertices and add them into inverted list.
	for (int i = 0; i < degree; i++) {
		int processNum = originalEdge[i].first;
		if (es[processNum].size() != originalDegree[i] && !isAltered[processNum]) {//If the degree of a neighbor has been changed, update the neighbor.
			if (es[processNum].size() != 1 && es[processNum].size() >= curMaxDegree) {
				isAltered[processNum] = 1;
				changedNeighbors.push(make_pair(processNum, originalDegree[i]));
			}
			else {//If the degree of neighbor is 1 or smaller than the minimum degree in inverted list, update the neighbor.
				addVertexToVertexList(processNum, originalDegree[i]);
			}
		}
	}
}

void graph::removeByDegree_directed(const int& number, const int& degree) {
	int edgeNum = degree * (degree - 1) / 2;//Calculate the number of edges to be updated.
	vector<edge> outEdgeProcess;//Record out-edges after removing the vertex.
	vector<edge> inEdgeProcess;//Record in-edges after removing the vertex.
	vector<edge> originalEdge_in;//Record in-edges before removing the vertex.
	originalEdge_in.reserve(degree);
	vector<edge> originalEdge_out;//Record in-edges before removing the vertex.
	originalEdge_out.reserve(degree);
	list<pair<int, int>> originalDegree;//Record the degree of neighbors before removing the vertex. First is the ID of neighbors. Second is the degree of neighbors.
	//originalDegree.reserve(2 * degree);
	vector<int> numberOfOutEdges;//Record the number of out-edges when removing the vertex.
	numberOfOutEdges.reserve(degree + 1);
	vector<int> numberOfInEdges;//Record the number of in-edges when removing the vertex.
	numberOfInEdges.reserve(degree + 1);
	edge temp;

	//Process edges to be added into the graph after removing the vertex.
	nodes[number].curNum = number;//Record the ID of current removed vertex.
	nodes[number].neighbor_in.reserve(degree);
	nodes[number].neighbor_out.reserve(degree);
	nodes[number].neighborDis_in.reserve(degree);
	nodes[number].neighborDis_out.reserve(degree);

	//Record distances, neighbors, original degree of neighbors, and original edges from es_directed[number].
	for (list<edge>::iterator it = es_directed[number].first.begin(); it != es_directed[number].first.end(); it++) {
		nodes[number].neighbor_out.push_back(it->v);//Record neighbors in nodes.
		nodes[number].neighborDis_out.push_back(it->weight);//Record distance of the out-edge.
		int degree = calDegree(it->v);
		originalDegree.push_back(make_pair(it->v, degree));
		originalEdge_out.push_back(edge(it->u, it->v, it->weight));
	}
	for (list<edge>::iterator it = es_directed[number].second.begin(); it != es_directed[number].second.end(); it++) {
		nodes[number].neighbor_in.push_back(it->u);//Record neighbors in nodes.
		nodes[number].neighborDis_in.push_back(it->weight);//Record distance of the in-edge.
		int degree = calDegree(it->u);
		originalDegree.push_back(make_pair(it->u, degree));
		originalEdge_in.push_back(edge(it->u, it->v, it->weight));
	}

	outEdgeProcess.reserve(originalEdge_in.size() * originalEdge_out.size());
	int numOfCurEdges = 0;
	//Calculate the out-edge weight of added edges.
	numberOfOutEdges.push_back(0);
	for (int i = 0; i < originalEdge_in.size(); i++) {
		for (int j = 0; j < originalEdge_out.size(); j++) {
			temp.u = originalEdge_in[i].u;
			temp.v = originalEdge_out[j].v;
			if (temp.u != temp.v) {
				temp.weight = originalEdge_in[i].weight + originalEdge_out[j].weight;	
			}
			else {//Exclude that the in-edge and out-edge are the same.
				temp.weight = INT_MAX;
			}
			outEdgeProcess.push_back(temp);
			numOfCurEdges++;
		}
		if (numOfCurEdges > 0 && numOfCurEdges != numberOfOutEdges.back())
			numberOfOutEdges.push_back(numOfCurEdges);
	}

	inEdgeProcess.reserve(originalEdge_in.size() * originalEdge_out.size());
	numOfCurEdges = 0;
	//Calculate the in-edge weight of added edges.
	numberOfInEdges.push_back(0);
	for (int i = 0; i < originalEdge_out.size(); i++) {
		for (int j = 0; j < originalEdge_in.size(); j++) {
			temp.u = originalEdge_in[j].u;
			temp.v = originalEdge_out[i].v;
			if (temp.u != temp.v) {
				temp.weight = originalEdge_in[j].weight + originalEdge_out[i].weight;
			}
			else {//Exclude that the in-edge and out-edge are the same.
				temp.weight = INT_MAX;
			}
			inEdgeProcess.push_back(temp);
			numOfCurEdges++;
		}
		if (numOfCurEdges > 0 && numOfCurEdges != numberOfInEdges.back())
			numberOfInEdges.push_back(numOfCurEdges);
	}

	//Remove current vertex and corresponding edges.
	isRemoved[number] = 1;
	nodes[number].level = sequence;//Set sequence.
	level[number] = sequence++;//Set sequence++.

	//Add out-edges into es_directed[u].first which is similar with merge-sort.
	for (int i = 0; i < numberOfOutEdges.size() - 1; i++) {
		int u = outEdgeProcess[numberOfOutEdges[i]].u;
		updateEdges_out(u, numberOfOutEdges[i + 1], number, numberOfOutEdges[i], outEdgeProcess);
	}

	//Add in-edges into es_directed[v].second which is similar with merge-sort.
	for (int i = 0; i < numberOfInEdges.size() - 1; i++) {
		int v = inEdgeProcess[numberOfInEdges[i]].v;
		updateEdges_in(v, numberOfInEdges[i + 1], number, numberOfInEdges[i], inEdgeProcess);
	}

	int curMinDegree = degree;
	if (vertexList.front()->size() != 0) {
		int minDegreeVertex = vertexList.front()->front().first;
		curMinDegree = calDegree(minDegreeVertex);
	}

	//Update vertices and add them into inverted list.
	originalDegree.sort([](pair<int, int> a, pair<int, int> b) {return a.first < b.first; });
	originalDegree.unique();//Every vertex should be processed once in directed graph.
	for (auto it = originalDegree.begin(); it != originalDegree.end(); it++) {
		int curDegree = calDegree(it->first); 
		if (curDegree != it->second && !isAltered[it->first]) {//If the degree of a neighbor has been changed, update the neighbor.
				// changedNeighbors.push(*it);
				// isAltered[it->first] = 1;
			if (curDegree < curMinDegree) {//The current processed vertex with smallest degree should be removed first in the next loop.
				// flag = 3;
				addVertexToVertexList_directed(it->first, it->second);
				// if (curDegree == 1) {//If degree of v becomes one, process v directly.
				// 	removeDegreeOne_directed(it->first);
				// }
			}
			else{
				changedNeighbors.push(*it);
				isAltered[it->first] = 1;
			}
		}
	}
}

void graph::updateEdges_out(const int& u, const int& outNum, const int& number, const int& calNum, const vector<edge>& edgeProcess) {
	int calEdgeNum = calNum;
	list<edge>::iterator it = es_directed[u].first.begin();
	flag = 0;
	while (it != es_directed[u].first.end() && it->u == u && calEdgeNum < outNum) { //Find edge (u,v).
		if (it->v == number) {//Delete edge £¨u,curNum£©.
			es_directed[u].first.erase(it++);
			flag = 1;
			continue;
		}
		if (it->v == edgeProcess[calEdgeNum].v) {//If the edge is the same as an edge of vertex u, compare the edge weight.
			if (it->weight > edgeProcess[calEdgeNum].weight) {
				it->weight = edgeProcess[calEdgeNum].weight;
			}
			it++;
			calEdgeNum++;
		}
		else if ((it->u == u && it->v < edgeProcess[calEdgeNum].v)) {
			it++;
		}
		else {
			if (edgeProcess[calEdgeNum].u != edgeProcess[calEdgeNum].v)
				es_directed[u].first.insert(it, edgeProcess[calEdgeNum++]);
			else
				calEdgeNum++;
		}
	}

	while (it != es_directed[u].first.end() && flag == 0 && it->u == u) {//Guarantee all the edges have been deleted.
		if (it->v == number) {//Delete edge £¨u, curNum£©.
			es_directed[u].first.erase(it++);
			break;
		}
		it++;
	}
	while (calEdgeNum < outNum) {//Add other edges into es_directed[u].
		if (edgeProcess[calEdgeNum].weight == INT_MAX) {
			calEdgeNum++;
			continue;
		}
		if (it != es_directed[u].first.end()) {
			es_directed[u].first.insert(it++, edgeProcess[calEdgeNum++]);
		}
		else {
			es_directed[u].first.push_back(edgeProcess[calEdgeNum++]);
		}
	}
}

void graph::updateEdges_in(const int& v, const int& inNum, const int& number, const int& calNum, const vector<edge>& edgeProcess) {
	int calEdgeNum = calNum;
	list<edge>::iterator it = es_directed[v].second.begin();
	flag = 0;
	while (it != es_directed[v].second.end() && it->v == v &&  calEdgeNum < inNum) { //Find edge (u,v).
		if (it->u == number) {//Delete edge £¨u,curNum£©.
			es_directed[v].second.erase(it++);
			flag = 1;
			continue;
		}
		if (it->u == edgeProcess[calEdgeNum].u) {//If the edge is the same as an edge of vertex u, compare the edge weight.
			if (it->weight > edgeProcess[calEdgeNum].weight) {
				it->weight = edgeProcess[calEdgeNum].weight;
			}
			it++;
			calEdgeNum++;
		}
		else if ((it->v == v && it->u < edgeProcess[calEdgeNum].u)) {
			it++;
		}
		else {
			if (edgeProcess[calEdgeNum].u != edgeProcess[calEdgeNum].v)
				es_directed[v].second.insert(it, edgeProcess[calEdgeNum++]);
			else
				calEdgeNum++;
		}
	}
	while (it != es_directed[v].second.end() && flag == 0 && it->v == v) {//Guarantee the edge have been deleted.
		if (it->u == number) {//Delete edge £¨u, curNum£©.
			es_directed[v].second.erase(it++);
			break;
		}
		it++;
	}
	while (calEdgeNum < inNum) {//Add other edges into es_directed[u].
		if (edgeProcess[calEdgeNum].weight == INT_MAX) {
			calEdgeNum++;
			continue;
		}
		if (it != es_directed[v].second.end()) {
			es_directed[v].second.insert(it++, edgeProcess[calEdgeNum++]);
		}
		else {
			es_directed[v].second.push_back(edgeProcess[calEdgeNum++]);
		}
	}
}

void graph::removeByValue(const int& number, const int& degree) {
	int edgeNum = degree * (degree - 1) / 2;//Calculate the number of edges to be updated.
	vector<edge> edgeProcess;//Record edges after removing the vertex.
	edgeProcess.reserve((degree - 1) * (degree - 1));
	vector<pair<int, int>> originalEdge(degree); //Record edges before removing the vertex.
	vector<int> originalDegree(degree);//Record the degree of neighbors before removing the vertex.
	edge temp;

	//Process edges to be added into the graph after removing the vertex.
	nodes[number].curNum = number;//Record the ID of current removed vertex.
	nodes[number].neighbor.reserve(degree);
	nodes[number].neighborDis.reserve(degree);
	int j = 0;
	for (auto it = es[number].begin(); it != es[number].end(); it++, j++) {
		nodes[number].neighbor.push_back((*it).v);//Record neighbors in nodes.
		nodes[number].neighborDis.push_back((*it).weight);//Record distance to neighbors in nodes.
		originalEdge[j].first = (*it).v;
		originalEdge[j].second = (*it).weight;
		originalDegree[j] = es[(*it).v].size();
	}

	//Calculate the edge weight of added edges.
	for (int i = 0; i < degree; i++) {
		for (int j = 0; j < degree; j++) {
			if (j != i) {
				temp.u = originalEdge[i].first;
				temp.v = originalEdge[j].first;
				temp.weight = originalEdge[i].second + originalEdge[j].second;
				edgeProcess.push_back(temp);
			}
		}
	}

	//Remove current vertex and corresponding edges.
	isRemoved[number] = 1;
	nodes[number].level = sequence;//Set sequence.
	level[number] = sequence++;//Set sequence++.
	int calEdgeNum = 0;

	//Add edges into es[u] which is similar with merge-sort.
	for (int i = 0; i < degree; i++) {
		int u = originalEdge[i].first;
		flag = 0;
		auto it = es[u].begin();
		for (; it != es[u].end() && calEdgeNum < (i + 1) * (degree - 1);) { //Find edge (u,v).
			if ((*it).v == number) {//Delete edge £¨u,curNum£©.
				es[u].erase(it++);
				flag = 1;
				continue;
			}
			if ((*it).v == edgeProcess[calEdgeNum].v) {//If the edge is the same as an edge of vertex u, compare the edge weight.
				if ((*it).weight > edgeProcess[calEdgeNum].weight) {
					(*it).weight = edgeProcess[calEdgeNum].weight;
				}
				it++;
				calEdgeNum++;
			}
			else if ((*it).v < edgeProcess[calEdgeNum].v) {
				it++;
			}
			else {
				edge* temp = new edge(edgeProcess[calEdgeNum].u, edgeProcess[calEdgeNum].v, edgeProcess[calEdgeNum].weight);
				calEdgeNum++;
				es[u].insert(it, *temp);
			}
		}
		while (flag == 0) {//guarantee all the edges have been deleted.
			if ((*it).v == number) {//Delete edge £¨u,curNum£©.
				es[u].erase(it);
				flag = 1;
				break;
			}
			it++;
		}
		while (calEdgeNum < (i + 1) * (degree - 1)) {//Add other edges into es[u].
			edge* temp = new edge(edgeProcess[calEdgeNum].u, edgeProcess[calEdgeNum].v, edgeProcess[calEdgeNum].weight);
			calEdgeNum++;
			es[u].push_back(*temp);
		}
	}

	int curMaxValue = 0;
	//Get the minimum value from inverted list.
	while (!vertexListIterator[curMaxValue] || (**(vertexListIterator[curMaxValue]))->size() == 0) {
		curMaxValue++;
	}

	//Update vertices and add them into inverted list.
	for (int i = 0; i < degree; i++) {
		int processNum = originalEdge[i].first;
		if (es[processNum].size() == 1) {//If the degree of neighbor is 1, insert the neighbor into the inverted list.
			isAltered[processNum] = 0;
			addVertexToVertexList_HFVertices(processNum, (*vertexIterator[processNum])->second);
		}
		else {
			if (!isAltered[processNum]) {
				if (es[processNum].size() != originalDegree[i]) {//If the degree of a neighbor has been changed, update the neighbor.
					if (calValue(processNum) >= curMaxValue) {
						isAltered[processNum] = 1;
						changedNeighbors.push(make_pair(processNum, (*vertexIterator[processNum])->second));
					}
					else {//If the degree of neighbor is smaller than the minimum degree in inverted list, update the neighbor.
						addVertexToVertexList_HFVertices(processNum, (*vertexIterator[processNum])->second);
					}
				}
			}
		}
	}
}

void graph::removeByValue_directed(const int& number, const int& degree) {
	int edgeNum = degree * (degree - 1) / 2;//Calculate the number of edges to be updated.
	vector<edge> outEdgeProcess;//Record out-edges after removing the vertex.
	vector<edge> inEdgeProcess;//Record in-edges after removing the vertex.
	vector<edge> originalEdge_in;//Record in-edges before removing the vertex.
	originalEdge_in.reserve(degree);
	vector<edge> originalEdge_out;//Record in-edges before removing the vertex.
	originalEdge_out.reserve(degree);
	list<pair<int, int>> originalDegree;//Record the degree of neighbors before removing the vertex. First is the ID of neighbors. Second is the degree of neighbors.
	//originalDegree.reserve(2 * degree);
	vector<int> numberOfOutEdges;//Record the number of out-edges when removing the vertex.
	numberOfOutEdges.reserve(degree + 1);
	vector<int> numberOfInEdges;//Record the number of in-edges when removing the vertex.
	numberOfInEdges.reserve(degree + 1);
	edge temp;

	//Process edges to be added into the graph after removing the vertex.
	nodes[number].curNum = number;//Record the ID of current removed vertex.
	nodes[number].neighbor_in.reserve(degree);
	nodes[number].neighbor_out.reserve(degree);
	nodes[number].neighborDis_in.reserve(degree);
	nodes[number].neighborDis_out.reserve(degree);

	//Record distances, neighbors, original degree of neighbors, and original edges from es_directed[number].
	for (list<edge>::iterator it = es_directed[number].first.begin(); it != es_directed[number].first.end(); it++) {
		nodes[number].neighbor_out.push_back(it->v);//Record neighbors in nodes.
		nodes[number].neighborDis_out.push_back(it->weight);//Record distance of the out-edge.
		int degree = calDegree(it->v);
		originalDegree.push_back(make_pair(it->v, degree));
		originalEdge_out.push_back(edge(it->u, it->v, it->weight));
	}
	for (list<edge>::iterator it = es_directed[number].second.begin(); it != es_directed[number].second.end(); it++) {
		nodes[number].neighbor_in.push_back(it->u);//Record neighbors in nodes.
		nodes[number].neighborDis_in.push_back(it->weight);//Record distance of the in-edge.
		int degree = calDegree(it->u);
		originalDegree.push_back(make_pair(it->u, degree));
		originalEdge_in.push_back(edge(it->u, it->v, it->weight));
	}

	outEdgeProcess.reserve(originalEdge_in.size() * originalEdge_out.size());
	int numOfCurEdges = 0;
	//Calculate the out-edge weight of added edges.
	numberOfOutEdges.push_back(0);
	for (int i = 0; i < originalEdge_in.size(); i++) {
		for (int j = 0; j < originalEdge_out.size(); j++) {
			temp.u = originalEdge_in[i].u;
			temp.v = originalEdge_out[j].v;
			if (temp.u != temp.v) {
				temp.weight = originalEdge_in[i].weight + originalEdge_out[j].weight;
			}
			else {//Exclude that the in-edge and out-edge are the same.
				temp.weight = INT_MAX;
			}
			outEdgeProcess.push_back(temp);
			numOfCurEdges++;
		}
		if (numOfCurEdges > 0 && numOfCurEdges != numberOfOutEdges.back())
			numberOfOutEdges.push_back(numOfCurEdges);
	}

	inEdgeProcess.reserve(originalEdge_in.size() * originalEdge_out.size());
	numOfCurEdges = 0;
	//Calculate the in-edge weight of added edges.
	numberOfInEdges.push_back(0);
	for (int i = 0; i < originalEdge_out.size(); i++) {
		for (int j = 0; j < originalEdge_in.size(); j++) {
			temp.u = originalEdge_in[j].u;
			temp.v = originalEdge_out[i].v;
			if (temp.u != temp.v) {
				temp.weight = originalEdge_in[j].weight + originalEdge_out[i].weight;
			}
			else {//Exclude that the in-edge and out-edge are the same.
				temp.weight = INT_MAX;
			}
			inEdgeProcess.push_back(temp);
			numOfCurEdges++;
		}
		if (numOfCurEdges > 0 && numOfCurEdges != numberOfInEdges.back())
			numberOfInEdges.push_back(numOfCurEdges);
	}

	//Remove current vertex and corresponding edges.
	isRemoved[number] = 1;
	nodes[number].level = sequence;//Set sequence.
	level[number] = sequence++;//Set sequence++.

	//Add out-edges into es_directed[u].first which is similar with merge-sort.
	for (int i = 0; i < numberOfOutEdges.size() - 1; i++) {
		int u = outEdgeProcess[numberOfOutEdges[i]].u;
		updateEdges_out(u, numberOfOutEdges[i + 1], number, numberOfOutEdges[i], outEdgeProcess);
	}

	//Add in-edges into es_directed[v].second which is similar with merge-sort.
	for (int i = 0; i < numberOfInEdges.size() - 1; i++) {
		int v = inEdgeProcess[numberOfInEdges[i]].v;
		updateEdges_in(v, numberOfInEdges[i + 1], number, numberOfInEdges[i], inEdgeProcess);
	}

	int curMinValue = calValue_directed(number);
	//Get the minimum value from inverted list.
	if (vertexList.front()->size() != 0) {
		int minValueVertex = vertexList.front()->front().first;
		curMinValue = calValue_directed(minValueVertex);
	}

	/*while (!vertexListIterator[curMaxValue] || (**(vertexListIterator[curMaxValue]))->size() == 0) {
		curMaxValue++;
	}*/

	//Update vertices and add them into inverted list.
	originalDegree.sort([](pair<int, int> a, pair<int, int> b) {return a.first < b.first; });
	originalDegree.unique();//Every vertex should be processed once in directed graph.
	for (auto it = originalDegree.begin(); it != originalDegree.end(); it++) {
		int curDegree = calDegree(it->first);
		if (curDegree != it->second && !isAltered[it->first]) {//If the degree of a neighbor has been changed, update the neighbor.
			// changedNeighbors.push(*it);
			// isAltered[it->first] = 1;
			if (calValue_directed(it->first) < curMinValue) {//The current processed vertex with smallest value should be removed first in the next loop.
				//flag = 3;
				addVertexToVertexList_HFVertices_directed(it->first, (*vertexIterator[it->first])->second);
			}
			else{
				changedNeighbors.push(*it);
				isAltered[it->first] = 1;
			}
		}
	}
}

inline int graph::calValue(const int& num) {
	if (es[num].size() == 1) {
		return 1;
	}
	int value = round(GAMMA * HFVertex[num].second + (1 - GAMMA) * es[num].size());
	return 2 > value ? 2 : value;//Make sure that vertices with degree one should be processed first.
}

inline int graph::calValue_directed(const int& num) {
	int degree = calDegree(num);
	if (degree == 1) {
		return 1;
	}
	int value = round(GAMMA * HFVertex[num].second + (1 - GAMMA) * degree);
	return 2 > value ? 2 : value;//Make sure that vertices with degree one should be processed first.
}

inline int graph::calDegree(const int& num) {
	int degree = es_directed[num].first.size() + es_directed[num].second.size();//Here minDegree is converted by in-degree plus out-degree.
	if (degree == 2 && es_directed[num].first.front().v == es_directed[num].second.front().u) {//Vertices with just one neighbor should be processed at first.
		degree = 1;
	}
	return degree;
}

void graph::removeVertices_H2H() {
	//Push vertices to the inverted list.
	for (int i = 0; i < numOfVertices; i++) {
		addVertexToVertexList_readData(i, es[i].size());
	}
	isAltered.assign(numOfVertices, 0);

	//Vertices with degree one must be processed at first.
	while (vertexListIterator[1]) {
		list<pair<int, int>>::iterator it = vertexList.front()->begin();
		removeDegreeOne(it->first);
		if (vertexList.front()->size() == 0) {
			vertexList.pop_front();//Erase the list with degree one.
			delete vertexListIterator[1];
			vertexListIterator[1] = NULL;
		}
	}

	int lastDeleteDegree{ -1 };
	int curFlag = 0;
	//Remove vertices by minimumDegree, the last vertex compose the root of tree decomposition.
	//In each while loop, just one vertex will be removed.
	while (sequence < numOfVertices - 1) {
		//Vertices with degree one must be processed at first.
		while (vertexListIterator[1] && sequence < numOfVertices - 1) {
			list<pair<int, int>>::iterator it = vertexList.front()->begin();
			removeDegreeOne(it->first);
			if (vertexList.front()->size() == 0) {
				vertexList.pop_front();//Erase the list with degree one.
				delete vertexListIterator[1];
				vertexListIterator[1] = NULL;
			}
		}
		if (sequence == numOfVertices - 1) {//The last vertex compose the root of tree decomposition, quit the procedure.
			break;
		}
		curFlag = 0;
		for (list<list<pair<int, int>>*>::iterator listIt = vertexList.begin(); listIt != vertexList.end() && curFlag == 0; listIt++) {
			//The minimum degree has been changed.
				//We update neighbors in batch instead of one-by-one.
				//That speeds up the update procedure and does not affect the overall time complexity.
			if (isAltered[(*listIt)->begin()->first]) {
				while (!changedNeighbors.empty()) {
					int curNum = changedNeighbors.front().first;
					addVertexToVertexList(curNum, changedNeighbors.front().second);
					isAltered[curNum] = 0;//Exclude vertices in changedNeighbors with the same ID.
					changedNeighbors.pop();
				}
				curFlag = 2;
				break;
			}
			for (list<pair<int, int>>::iterator curIt = (*listIt)->begin(); curIt != (*listIt)->end(); curIt++) {
				if (sequence < numOfVertices - 1) {
					lastDeleteDegree = es[curIt->first].size();
					removeByDegree(curIt->first, lastDeleteDegree);
					(*listIt)->erase(curIt);
					curFlag = 1;
					break;
				}
			}
			//If the current processed list becomes empty, clear the list.
			if (curFlag == 1 && vertexListIterator[lastDeleteDegree] && (*listIt)->size() == 0) {
				vertexList.erase(listIt);//
				delete vertexListIterator[lastDeleteDegree];
				vertexListIterator[lastDeleteDegree] = NULL;
				break;
			}
		}
	}
	//Set the last vertex as the root node.
	level[vertexList.front()->front().first] = numOfVertices - 1;
	//Reorder the neighbors and distances in nodes.
	for (int i = 0; i < numOfVertices; i++) {
		nodes[i].reOrder(level);
	}
}

void graph::removeVertices_H2H_directed() {
	//Push vertices to the inverted list.
	for (int i = 0; i < numOfVertices; i++) {
		int curDegree = calDegree(i);
		addVertexToVertexList_readData_directed(i, curDegree);
	}
	isAltered.assign(numOfVertices, 0);

	//Vertices with degree one must be processed at first.
	while (vertexListIterator[1]) {
		list<pair<int, int>>::iterator it = vertexList.front()->begin();
		removeDegreeOne_directed(it->first);
		if (vertexList.front()->size() == 0) {
			vertexList.pop_front();//Erase the list with degree one.
			delete vertexListIterator[1];
			vertexListIterator[1] = NULL;
		}
	}

	int lastDeleteDegree{ -1 };
	int curFlag = 0;
	//Remove vertices by minimumDegree, the last vertex compose the root of tree decomposition.
	//In each while loop, just one vertex will be removed.
	while (sequence < numOfVertices - 1) {
		//Vertices with degree one must be processed at first.
		while (vertexListIterator[1] && sequence < numOfVertices - 1) {
			list<pair<int, int>>::iterator it = vertexList.front()->begin();
			removeDegreeOne_directed(it->first);
			if (vertexList.front()->size() == 0) {
				vertexList.pop_front();//Erase the list with degree one.
				delete vertexListIterator[1];
				vertexListIterator[1] = NULL;
			}
		}
		if (sequence == numOfVertices - 1) {//The last vertex compose the root of tree decomposition, quit the procedure.
			break;
		}
		curFlag = 0;
		//if (sequence == 14) {
		//	sequence = sequence;
		//}
		for (list<list<pair<int, int>>*>::iterator listIt = vertexList.begin(); listIt != vertexList.end() && curFlag == 0; listIt++) {
			//The minimum degree has been changed.
			//We update neighbors in batch instead of one-by-one.
			//That speeds up the update procedure and does not affect the overall time complexity.
			//The current vertex with smallest degree should be removed first in the next loop.
			if (/*flag == 3 ||*/ isAltered[(*listIt)->begin()->first]) {
				while (!changedNeighbors.empty()) {
					int curNum = changedNeighbors.front().first;
					addVertexToVertexList_directed(curNum, changedNeighbors.front().second);
					isAltered[curNum] = 0;//Exclude vertices in changedNeighbors with the same ID.
					changedNeighbors.pop();
				}
				//flag = 0;
				curFlag = 2;
				break;
			}
			for (list<pair<int, int>>::iterator curIt = (*listIt)->begin(); curIt != (*listIt)->end(); curIt++) {
				if (sequence < numOfVertices - 1) {
					lastDeleteDegree = calDegree(curIt->first); 
					removeByDegree_directed(curIt->first, lastDeleteDegree);
					(*listIt)->erase(curIt);
					curFlag = 1;
					break;
				}
			}
			//If the current processed list becomes empty, clear the list.
			if (curFlag == 1 && vertexListIterator[lastDeleteDegree] && (*listIt)->size() == 0) {
				vertexList.erase(listIt);//
				delete vertexListIterator[lastDeleteDegree];
				vertexListIterator[lastDeleteDegree] = NULL;
				break;
			}
		}
	}
	//Set the last vertex as the root node.
	level[vertexList.front()->front().first] = numOfVertices - 1;
	//Reorder the neighbors and distances in nodes.
	for (int i = 0; i < numOfVertices; i++) {
		nodes[i].reOrder_directed(level);
	}
}

void graph::removeVertices_WH2H() {
	//Push vertices to the inverted list.
	int curValue = 0;
	for (int i = 0; i < numOfVertices; i++) {
		curValue = calValue(i);
		addVertexToVertexList_readData_HFVertices(i, curValue);
	}
	isAltered.assign(numOfVertices, 0);

	//Vertices with degree one must be processed at first.
	while (vertexListIterator[1]) {
		list<pair<int, int>>::iterator it = vertexList.front()->begin();
		removeDegreeOne_HFVertices(it->first);
		if (vertexList.front()->size() == 0) {
			vertexList.pop_front();//Erase the list with degree one.
			delete vertexListIterator[1];
			vertexListIterator[1] = NULL;
		}
	}
	int lastDeleteValue{ -1 };
	int lastDeleteDegree{ -1 };
	int curFlag = 0;
	//Remove vertices by minimumDegree, the last vertex compose the root of tree decomposition.
	//In each while loop, just one vertex will be removed.
	while (sequence < numOfVertices - 1) {
		//Vertices with degree one must be processed at first.
		while (vertexListIterator[1] && sequence < numOfVertices - 1) {
			list<pair<int, int>>::iterator it = vertexList.front()->begin();
			removeDegreeOne_HFVertices(it->first);
			if (vertexList.front()->size() == 0) {
				vertexList.pop_front();//Erase the list with degree one.
				delete vertexListIterator[1];
				vertexListIterator[1] = NULL;
			}
		}
		if (sequence == numOfVertices - 1) {//The last vertex compose the root of tree decomposition, quit the procedure.
			break;
		}

		curFlag = 0;
		int nodeNum = 0;
		for (list<list<pair<int, int>>*>::iterator listIt = vertexList.begin(); curFlag == 0 && listIt != vertexList.end(); listIt++) {
			if (isAltered[(*listIt)->begin()->first]) {//The minimum value has been changed.
				while (!changedNeighbors.empty()) {
					int curNum = changedNeighbors.front().first;
					if (level[curNum] != -1) {
						isAltered[curNum] = 0;
						changedNeighbors.pop();
						continue;
					}
					addVertexToVertexList_HFVertices(curNum, (*vertexIterator[curNum])->second);
					isAltered[curNum] = 0;//Exclude vertices in changedNeighbors with the same ID.
					changedNeighbors.pop();
				}
				curFlag = 2;
				break;
			}
			for (list<pair<int, int>>::iterator curIt = (*listIt)->begin(); curIt != (*listIt)->end(); curIt++) {
				lastDeleteDegree = es[curIt->first].size();
				lastDeleteValue = curIt->second;
				nodeNum = curIt->first;
				removeByValue(curIt->first, lastDeleteDegree);
				if (*listIt)
					(*listIt)->erase(curIt);
				curFlag = 1;
				break;
			}
			//If the current processed list becomes empty, clear the list.
			if (curFlag == 1 && (*listIt)->size() == 0 && vertexListIterator[lastDeleteValue] != NULL) {
				vertexList.erase(listIt);//
				delete vertexListIterator[lastDeleteValue];
				vertexListIterator[lastDeleteValue] = NULL;
				break;
			}
		}
	}
	//Set the last vertex as the root node.
	level[vertexList.front()->front().first] = numOfVertices - 1;
	//Reorder the neighbors and distances in nodes.
	for (int i = 0; i < numOfVertices; i++) {
		nodes[i].reOrder(level);
	}
}

void graph::removeVertices_WH2H_directed() {
	//Push vertices to the inverted list.
	int curValue = 0;
	for (int i = 0; i < numOfVertices; i++) {
		curValue = calValue_directed(i);
		addVertexToVertexList_readData_HFVertices_directed(i, curValue);
	}
	isAltered.assign(numOfVertices, 0);

	//Vertices with degree one must be processed at first.
	while (vertexListIterator[1]) {
		list<pair<int, int>>::iterator it = vertexList.front()->begin();
		removeDegreeOne_HFVertices_directed(it->first);
		if (vertexList.front()->size() == 0) {
			vertexList.pop_front();//Erase the list with degree one.
			delete vertexListIterator[1];
			vertexListIterator[1] = NULL;
		}
	}
	int lastDeleteValue{ -1 };
	int lastDeleteDegree{ -1 };
	int curFlag = 0;
	//Remove vertices by minimumDegree, the last vertex compose the root of tree decomposition.
	//In each while loop, just one vertex will be removed.
	while (sequence < numOfVertices - 1) {
		//Vertices with degree one must be processed at first.
		while (vertexListIterator[1] && sequence < numOfVertices - 1) {
			list<pair<int, int>>::iterator it = vertexList.front()->begin();
			removeDegreeOne_HFVertices_directed(it->first);
			if (vertexList.front()->size() == 0) {
				vertexList.pop_front();//Erase the list with degree one.
				delete vertexListIterator[1];
				vertexListIterator[1] = NULL;
			}
		}
		if (sequence == numOfVertices - 1) {//The last vertex compose the root of tree decomposition, quit the procedure.
			break;
		}

		curFlag = 0;
		int nodeNum = 0;
		for (list<list<pair<int, int>>*>::iterator listIt = vertexList.begin(); curFlag == 0 && listIt != vertexList.end(); listIt++) {
			//The minimum degree has been changed.
			//We update neighbors in batch instead of one-by-one.
			//That speeds up the update procedure and does not affect the overall time complexity.
			//The current vertex with smallest degree should be removed first in the next loop.
			if (/*3 == flag ||*/ isAltered[(*listIt)->begin()->first]) {//The minimum value has been changed.
				while (!changedNeighbors.empty()) {
					int curNum = changedNeighbors.front().first;
					//if (level[curNum] != -1) {
					//	isAltered[curNum] = 0;
					//	changedNeighbors.pop();
					//	continue;
					//}
					addVertexToVertexList_HFVertices_directed(curNum, (*vertexIterator[curNum])->second);
					isAltered[curNum] = 0;//Exclude vertices in changedNeighbors with the same ID.
					changedNeighbors.pop();
				}
				curFlag = 2;
				//flag = 0;
				break;
			}
			for (list<pair<int, int>>::iterator curIt = (*listIt)->begin(); curIt != (*listIt)->end(); curIt++) {
				if (sequence < numOfVertices - 1) {
					lastDeleteDegree = calDegree(curIt->first); //es[curIt->first].size();
					lastDeleteValue = curIt->second;
					nodeNum = curIt->first;
					removeByValue_directed(curIt->first, lastDeleteDegree);
					//if (*listIt)
					(*listIt)->erase(curIt);
					curFlag = 1;
					break;
				}
			}
			//If the current processed list becomes empty, clear the list.
			if (curFlag == 1 && (*listIt)->size() == 0 && vertexListIterator[lastDeleteValue] != NULL) {
				vertexList.erase(listIt);//
				delete vertexListIterator[lastDeleteValue];
				vertexListIterator[lastDeleteValue] = NULL;
				break;
			}
		}
	}
	//Set the last vertex as the root node.
	level[vertexList.front()->front().first] = numOfVertices - 1;
	//Reorder the neighbors and distances in nodes.
	for (int i = 0; i < numOfVertices; i++) {
		nodes[i].reOrder_directed(level);
	}
}

int graph::dijkstra(int x, int y) {

	const int INF = 0x3f3f3f3f;
	const int N = numOfVertices;

	vector<vector<int>> graph(N);
	for (int i = 0; i < numOfVertices; i++) {
		graph[i].resize(numOfVertices);
	}
	vector<int> check;
	check.assign(N, 0);
	vector<int> dis;
	dis.assign(N, 0);
	int n;


	int m, a, b, c;
	int i, p, j;
	n = numOfVertices;
	m = numOfEdges;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			if (i == j)
				graph[i][j] = 0;
			else
				graph[i][j] = INF;
		}

	for (int i = 0; i < numOfVertices; i++) {
		for (list<edge>::iterator it = es_dijkstra[i].first.begin(); it != es_dijkstra[i].first.end(); it++) {
			graph[it->u][it->v] = it->weight;
		}
		for (list<edge>::iterator it = es_dijkstra[i].second.begin(); it != es_dijkstra[i].second.end(); it++) {
			graph[it->u][it->v] = it->weight;
		}
	}

	int u0 = x;
	int inf, mid;
	for (i = 0; i < n; i++)
		dis[i] = graph[u0][i];
	check[u0] = -1;
	for (i = 0; i < n; i++)
	{
		inf = INF;
		mid = u0;
		for (j = 0; j < n; j++)
		{
			if (check[j] != -1 && dis[j] < inf)
			{
				inf = dis[j];
				mid = j;
			}
		}
		if (inf != INF)
		{
			check[mid] = -1;
			for (j = 0; j < n; j++)
			{
				if (check[j] != -1 && dis[j] > graph[mid][j] + dis[mid])
					dis[j] = graph[mid][j] + dis[mid];
			}
		}
	}
	for (int i = 0; i < n; i++) {
		dis_dijkstra[x][i] = dis[i];
	}
	return dis[y];
}
