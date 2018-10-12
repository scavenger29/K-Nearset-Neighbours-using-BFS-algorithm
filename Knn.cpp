#include <bits/stdc++.h>
#include <time.h>
#include <chrono>
using namespace std;

int comDim; // For using the comparator

class point{
public:
	int index;
	vector<double> dimVals; // This vector will store coordinates of the data
	point(int id, vector<double> &vals){
		index  = id;
		dimVals = vals;
	}
};

vector<point> dataPoints; // making this array global

bool comparator(int a, int b){
	if(dataPoints[a].dimVals[comDim] < dataPoints[b].dimVals[comDim]){
		return true;
	}
	else if(dataPoints[a].dimVals[comDim] > dataPoints[b].dimVals[comDim]){
		return false;
	}
	else{
		return (dataPoints[a].index < dataPoints[b].index);
	}
}

bool comparator1(pair<double, int> a, pair<double, int> b){
	if(a.first < b.first) return true;
	else if(a.first>b.first) return false;
	else{
		int dimension = dataPoints[a.second].dimVals.size();
		for(int i=0;i<dimension;i++){
			if(dataPoints[a.second].dimVals[i]<dataPoints[b.second].dimVals[i]) return true;
			else if(dataPoints[a.second].dimVals[i]>dataPoints[b.second].dimVals[i]) return false;
		}
		return false;
	}
}

bool comparator2(point a, point b){
	int dimension = a.dimVals.size();
	for(int i=0;i<dimension;i++){
		if(a.dimVals[i]>b.dimVals[i]) return false;
		else if(a.dimVals[i] < b.dimVals[i]) return true;
	}
	return true;
}

double findMedian(vector<point> &v, int dim){
	int k = v.size();
	if(k%2==0){
		return (v[k/2].dimVals[dim] + v[(k-1)/2].dimVals[dim])/2.0;
	}
	else{
		return v[k/2].dimVals[dim];
	}
}

class node{
public:
	bool isLeaf;
	int splitDim;
	double splitVal;
	int pointIndex; // Not null if the leaf node is the leafNode
	node *left, *right;
	int numPoints;
	vector< pair<double, double> > MBR; // will store min and max of each dimension. // Will correspond to the Minimum bounding region
	
	node(bool b, int a, double d, int p, int c, vector< pair<double, double> > &v){
		isLeaf = b;
		splitDim = a;
		splitVal = d;
		pointIndex = p;
		numPoints = c;
		MBR = v;
		left = NULL;
		right = NULL;
	}

	void dfs(){
		cout << "************" << endl;
		cout << isLeaf << " " << splitDim << " " << splitVal << " " << numPoints << endl;
		if(isLeaf){
			cout << this->pointIndex << endl;
			cout << dataPoints[this->pointIndex].dimVals[0] << " " << dataPoints[this->pointIndex].dimVals[1] << endl;
		}
		else{
			if((this->left) != NULL)	(this->left)->dfs();
			if((this->right) != NULL)	(this->right)->dfs();
		}
	}

};

node *constructNode(vector< vector<int> > &allDataPointsSorted, int currentLevel, int d){
	if(allDataPointsSorted.size()==0) return NULL;

	// current spliting dimension is currentLevel%d
	// allDataPointsSorted stores the points in this particular grid in sorted manner as per all possible dimensions - d
	int splitDim = currentLevel%d;
	//double median = findMedian(allDataPointsSorted[splitDim], splitDim);

	int numPoints = allDataPointsSorted[0].size();
	int medianIndex;
	if(numPoints%2) medianIndex = numPoints/2;
	else medianIndex = numPoints/2 - 1;
	double median = dataPoints[allDataPointsSorted[splitDim][medianIndex]].dimVals[splitDim];
	
	set<int> indices;
	for(int i=0;i<=medianIndex;i++){
		indices.insert(allDataPointsSorted[splitDim][i]);
	}

	vector< pair<double, double> > v;
	for(int i=0;i<d;i++){
		v.push_back({dataPoints[allDataPointsSorted[i][0]].dimVals[i], dataPoints[allDataPointsSorted[i][numPoints-1]].dimVals[i]});
	}
	bool isLeaf = false;
	if(allDataPointsSorted[0].size()==1){
		isLeaf = true;
	}

	if(isLeaf){
		
		node *nd = new node(isLeaf, splitDim, median, allDataPointsSorted[0][0], numPoints, v);
		 
		allDataPointsSorted.clear();
		return nd;
	}
	else{
		vector< vector<int> > allDataPointsSortedLeft, allDataPointsSortedRight;
		// less than equal to goes to left
		// more than goes to right
		for(int i=0;i<d;i++){
			vector<int> temp1, temp2;
			for(int j=0;j<allDataPointsSorted[i].size();j++){
				if(indices.count(allDataPointsSorted[i][j])==1){
					temp1.push_back(allDataPointsSorted[i][j]);
				}
				else temp2.push_back(allDataPointsSorted[i][j]);
			}
			if(temp1.size()!=0) allDataPointsSortedLeft.push_back(temp1);
			if(temp2.size()!=0) allDataPointsSortedRight.push_back(temp2);
		}
		allDataPointsSorted.clear();
		node *nd = new node(isLeaf, splitDim, median, -1, numPoints, v);
		
		nd->left = constructNode(allDataPointsSortedLeft, currentLevel+1, d);
		nd->right = constructNode(allDataPointsSortedRight, currentLevel+1, d);
		return nd;
	}
}

double L2distance(point &a, point &b, int d){
	double sum = 0;
	for(int i=0;i<d;i++){
		sum += (a.dimVals[i]-b.dimVals[i])*((a.dimVals[i]-b.dimVals[i]));
	}
	return sqrt(sum);
}

double MBR_dist(point &queryPoint, vector< pair<double, double> > &MBR, int d)
{
	double eucld_dist=0;
	double delt_i = 0;
	for(int i=0;i<d;i++){
		if(queryPoint.dimVals[i] < MBR[i].first)
			delt_i = MBR[i].first - queryPoint.dimVals[i];
		else if (queryPoint.dimVals[i] > MBR[i].second)
			delt_i = queryPoint.dimVals[i] - MBR[i].second;
		else
			delt_i = 0;
		eucld_dist += delt_i * delt_i ;
	}
	eucld_dist = sqrt(eucld_dist);
	return eucld_dist;
}

class kdTree{
public:
	node *root;

	kdTree(){
		root = NULL;
	}

	void constructKDtree(vector< vector<int> > &allDataPointsSorted, int currentLevel, int d){
		root = constructNode(allDataPointsSorted, currentLevel, d);
	}

	void KDtreeDfs(){
		root->dfs();
	}

	vector< pair<double, int> > kNN(int k, point queryPoint, int d){
		
		priority_queue< pair<double, int> > maxheap;
		priority_queue <  pair<double, node*>, vector<  pair<double, node*> >, greater<  pair<double, node*> > > minheap;

		double eucld_dist = MBR_dist(queryPoint, root->MBR, d);

		node *alt_root = root;
		minheap.push({eucld_dist, alt_root});

		
		while(minheap.size()!=0)
		{
			
			double dist = minheap.top().first;
			node *nd = minheap.top().second;
			minheap.pop();
			if(maxheap.size()==k && (dist> maxheap.top().first)){
				break;
			}
			if(maxheap.size()<k){
				if (nd->pointIndex != -1){
					maxheap.push({dist,nd->pointIndex});
				}
				else{
					node *left_nd = nd->left;
					node *right_nd = nd->right;
					double left_dist = MBR_dist(queryPoint, left_nd->MBR, d);
					double right_dist = MBR_dist(queryPoint, right_nd->MBR, d);
					minheap.push({left_dist, left_nd});
					minheap.push({right_dist, right_nd});	
				}
			}
			else{
				double top_dist = maxheap.top().first;
				double top_index = maxheap.top().second;
				if(dist > top_dist || (dist==top_dist && top_index< (nd->pointIndex))){
					continue;
				}
				else{
					if (nd->pointIndex != -1){
							maxheap.pop();
							maxheap.push({dist,nd->pointIndex});
					}
					else{
						node *left_nd = nd->left;
						node *right_nd = nd->right;
						double left_dist = MBR_dist(queryPoint, left_nd->MBR, d);
						double right_dist = MBR_dist(queryPoint, right_nd->MBR, d);
						if(left_dist <= top_dist)
							minheap.push({left_dist, left_nd});
						if(right_dist <= top_dist)
							minheap.push({right_dist, right_nd});
					}
				}	
			}
		}
		vector< pair<double, int> > ans;
		while(!maxheap.empty()){
			ans.push_back(maxheap.top());
			maxheap.pop();
		}
		return ans;
	}

	vector< pair<double, int> > sequentialScan(int k, point queryPoint, int d){
		priority_queue< pair<double, int> > maxHeap;
		for(int i=0;i<dataPoints.size();i++){
			double dist = L2distance(queryPoint, dataPoints[i], d);
			if(maxHeap.size()<k){
				maxHeap.push({dist, i});
			}
			else{
				if(maxHeap.top().first > dist){
					maxHeap.pop();
					maxHeap.push({dist, i});
				}
				else if(maxHeap.top().first == dist && maxHeap.top().second > i){
					maxHeap.pop();
					maxHeap.push({dist, i});
				}
			}
		}
		vector< pair<double, int> > ans;
		while(!maxHeap.empty()){
			ans.push_back(maxHeap.top());
			maxHeap.pop();
		}
		return ans; // It has points in decreasing order of distances.
	}
};


int main(int argc, char *argv[]){
	using namespace std::chrono;

	int d; 
	int numDataPoints;
	// d is the dimension of the kd tree
	vector< vector<int> > allDataPointsSorted;
	ifstream inFile;
	
	inFile.open(argv[1]);
	inFile >> d;
	inFile >> numDataPoints;

	for(int i=0;i<numDataPoints;i++){
		vector<double> temp;
		
		for(int j=0;j<d;j++){
			double dou;
			inFile >> dou;
			temp.push_back(dou);
		}
		dataPoints.push_back(point(i, temp));
	}
	inFile.close();

	sort(dataPoints.begin(), dataPoints.end(), comparator2);

	for(int i=0;i<numDataPoints;i++){
		dataPoints[i].index = i;
	}
	
	for(int i=0;i<d;i++){
		vector<int> temp;
		for(int i=0;i<numDataPoints;i++){
			temp.push_back(i);
		}
		comDim = i;
		sort(temp.begin(), temp.end(), comparator);
		
		allDataPointsSorted.push_back(temp);
	}


    // cout << "Total running time of the code is " << difftime(tfinish1,tinit) << " seconds" << endl;

	// Now i have d vectors, each sorted along the one direction.
	kdTree myKDtree = kdTree();
	myKDtree.constructKDtree(allDataPointsSorted, 0, d);


    // cout << "Total running time of the code is " << difftime(tfinish,tinit) << " seconds" << endl;

	cout << 0 << endl;

	// myKDtree.KDtreeDfs();
	
	// allDataPointsSorted will get cleared, so if we need it anywhere in future, let me know, I may need to recode somethings
	string queryFile;
	int k;
	cin >> queryFile;
	cin >> k; // for the knn query
	int number=0; //number of querypoints
	inFile.open(queryFile);
	
	inFile >> d;
	inFile >> number;
	// now take input for the query point
	ofstream ansFile;
	ansFile.open("results.txt");
	double time_knn = 0;
	double time_seq = 0;
	while(number--){
		vector<double> qDimVals;
		for(int i=0;i<d;i++){
			double a;
			inFile >> a;
			qDimVals.push_back(a);
		}
		
		point queryPoint(-1, qDimVals);
		high_resolution_clock::time_point t0 = high_resolution_clock::now();
		 vector< pair<double, int> > seqScanAns = myKDtree.sequentialScan(k, queryPoint, d);
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		vector< pair<double, int> > bestFirstAns = myKDtree.kNN(k, queryPoint, d);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double,std::milli> time_span = t2 - t1;
		duration<double,std::milli> time_span_sequential = t1 - t0;
		//cout<<time_span.count()<<" "<<number<<endl;
		time_knn+=time_span.count();
		time_seq +=  time_span_sequential.count();
		// sort(seqScanAns.begin(), seqScanAns.end(), comparator1);
		// sort(bestFirstAns.begin(), bestFirstAns.end(), comparator1);

		
		
		// ansFile.open("results1.txt");

		// for(int i=seqScanAns.size()-1;i>=0;i--){
		// 	for(int j=0;j<d;j++){
		// 		ansFile << dataPoints[seqScanAns[i].second].dimVals[j];
		// 		if(j!=d-1) ansFile<< " ";
		// 	}
		// 	if(i!=0) ansFile << endl;
		// }
		// ansFile.close();
		

		
		for(int i=bestFirstAns.size()-1;i>=0;i--){
			for(int j=0;j<d;j++){
				ansFile << dataPoints[bestFirstAns[i].second].dimVals[j];
				if(j!=d-1) ansFile<< " ";
			}
			if(i!=0) ansFile << endl;
		}
		ansFile << endl;
	}
	inFile.close();
	ansFile.close();
     cout << "Total running time of the knn code is " << time_knn<< " milliseconds" << endl;
     cout << "Total running time of the seq code is " << time_seq<< " milliseconds" << endl;
	cout << 1 << endl;
}