#include <bits/stdc++.h>
#include "omp.h"
#include <sys/time.h>
using namespace std;
int DIM;
float epsilon;
int total_data_points;
vector<vector<float> > input_set; //contains the input data points
vector<int> vec_0_N;
int N,M,K;
float mincost;
vector<vector<float> > opti_centers;
vector<vector<int> > all_subset_points;
vector<vector<int> > k_indices_of_subset_points;

void print( vector<int> l){
    for(vector<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}

void print( vector<float> l){
    for(vector<float>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}

float calc_distance(vector<float> &v1,vector<float> &v2)
{
	int sz = v1.size();
	float distsq = 0;
	for(int i=0;i<sz;i++)
	{
		distsq += ( (v1[i]-v2[i]) * (v1[i]-v2[i]) );
	}
	//cout << "distsq " <<  distsq <<  endl;
	return distsq;
}

vector<int> d2sampling(vector<vector<float> > &centers)
{
	map<pair<float,int>, int> mymap;
	vector<int> ans;
	for(int i=0;i<total_data_points;i++)
	{
		float mindist = INT_MAX;
		for(int j=0;j<centers.size();j++)
		{
			mindist = min(calc_distance(centers[j],input_set[i]),mindist );
		}
		mymap.insert(make_pair(make_pair(mindist,i),0 ));
	}
	std::map<pair<float,int>, int>::reverse_iterator it;
	int count = 0;
	for (it=mymap.rbegin(); it!=mymap.rend(); ++it)
	{
		if(count >= N)
			break;
		ans.push_back(it->first.second);
		//cout << "dist " << it->first.first << endl;
		count++;
	}
	//cout << "finished" << endl;
	//print(ans);
	return ans;
}

float cost(vector<vector<float> > &centers)
{
	float cst = 0;
	for(int i=0;i<total_data_points;i++)
	{
		float mindist = INT_MAX;
		for(int j=0;j<centers.size();j++)
		{
			mindist = min(calc_distance(centers[j],input_set[i]),mindist );
		}
		cst += mindist;
	}
	//cout << cst << " cost " << endl;
	return cst;
}


void generate_all_subsets(vector<int> data,int size, int left, int index)
{
	if(left==0){
			all_subset_points.push_back(data);
	}
    for(int i=index;i<size;i++){
    	data.push_back(vec_0_N[i]);
    	generate_all_subsets(data,size,left-1,i+1);
    	data.pop_back();
    }
}
///////////////////////////////////////////////
void subsets(vector<int> data,vector<vector<float> > &centers,vector<int> index_of_n_sampled_points,int size, int left, int index)
{
	if(left==0){
			vector<float> mean;
			for(int i=0;i<DIM;i++)
				mean.push_back(0.0);
			int data_size = data.size();
			for(int i = 0 ; i < data_size;i++)
			{
				for(int j=0;j<DIM;j++)
				{
					mean[j]+= ((input_set[index_of_n_sampled_points[data[i]]][j] *1.0)/(data_size));
				}
			}
			centers.push_back(mean);
			if(centers.size() == K)
			{
				float new_cost = cost(centers);
				//cout << "new cost " << centers.size() << endl;
				if(new_cost < mincost)
				{
					mincost = new_cost;
					for(int i=0;i<opti_centers.size();i++)
					{
						opti_centers[i].erase(opti_centers[i].begin(),opti_centers[i].end());
					}
					opti_centers.erase(opti_centers.begin(),opti_centers.end());
					for(int i=0;i<centers.size();i++)
					{
						opti_centers.push_back(centers[i]);
					}
					cout << mincost <<  " " << new_cost  << endl;
				}

			}
			else
			{
				data.erase(data.begin(),data.end());
				subsets(data,centers,d2sampling(centers),N,M,0);

			}
			centers.pop_back();

		return;
	}
    for(int i=index;i<size;i++){
    	data.push_back(vec_0_N[i]);
    	subsets(data,centers,index_of_n_sampled_points,size,left-1,i+1);
    	data.pop_back();
    }
}

///////////////////////////////////////////////


void iterative_subset_open_mp(int i,int num_points)
{
	vector<int> tmp1;
	tmp1.push_back(i);
	//cout << i << " thread no. " << omp_get_thread_num() << endl;
	queue<vector<int> > k_indices_of_subset_points_t;
	k_indices_of_subset_points_t.push(tmp1);
	while((k_indices_of_subset_points_t.front()).size() < K)
	{
			vector<int> tmp = k_indices_of_subset_points_t.front();
			k_indices_of_subset_points_t.pop();
			//k_indices_of_subset_points.erase(k_indices_of_subset_points.begin());
			for(int j=0;j<num_points;j++)
			{
				tmp.push_back(j);
				k_indices_of_subset_points_t.push(tmp);
				tmp.pop_back();
			}
	}
	#pragma omp critical
	{
		while(!k_indices_of_subset_points_t.empty())
		{
			k_indices_of_subset_points.push_back(k_indices_of_subset_points_t.front());
			k_indices_of_subset_points_t.pop();
		}
	}
}



vector<int> generate_k_random(int n, int k)
{
	vector<int> v,v1;
	epsilon = 1/20;
	for(int i=0;i<n;i++)
		v1.push_back(i);
	for(int i=0;i<k;i++)
	{
		int randint = rand()%(n-i)+i;
		v.push_back(v1[randint]);
		iter_swap(v1.begin()+i,v1.begin()+randint);
	}
	return v;
}
int countn;
map<vector<vector<float> >,vector<int> > dp;
void main_openmp(vector<int> indices,vector<int> index_of_n_sampled_points)
{
	vector<vector<float> > centers;
	vector<float> mean;
	vector<int> data = all_subset_points[indices[0]];
	for(int i=0;i<DIM;i++)
		mean.push_back(0.0);
	int data_size = data.size();
	for(int i = 0 ; i < data_size;i++)
	{
		for(int j=0;j<DIM;j++)
		{
			mean[j]+= ((input_set[index_of_n_sampled_points[data[i]]][j] *1.0)/(data_size));
		}
	}
	centers.push_back(mean);
	for(int i=1;i<K;i++)
	{
		vector<float> mean1;
		vector<int> data1 = all_subset_points[indices[i]];
		for(int i=0;i<DIM;i++)
			mean1.push_back(0.0);
		int data_size1 = data1.size();

		map<vector<vector<float> >,vector<int> >::iterator it = dp.find(centers);
		if(it != dp.end())
		{
			index_of_n_sampled_points = it->second;
		}
		else
		{
			index_of_n_sampled_points = d2sampling(centers);
			#pragma omp critical
			{
				dp[centers] = index_of_n_sampled_points;
				//countn++;
				//cout << countn << endl;
			}
		}
		for(int i = 0 ; i < data_size;i++)
		{
			for(int j=0;j<DIM;j++)
			{
				mean1[j]+= ((input_set[index_of_n_sampled_points[data1[i]]][j] *1.0)/(data_size1));
			}
		}
		centers.push_back(mean1);
	}
	float new_cost = cost(centers);
	#pragma omp critical
	{			//cout << "new cost " << centers.size() << endl;
		if(new_cost < mincost)
		{
			mincost = new_cost;
			for(int i=0;i<opti_centers.size();i++)
			{
				opti_centers[i].erase(opti_centers[i].begin(),opti_centers[i].end());
			}
			opti_centers.erase(opti_centers.begin(),opti_centers.end());
			for(int i=0;i<centers.size();i++)
			{
				opti_centers.push_back(centers[i]);
			}
			cout << mincost <<  " " << new_cost  << endl;
		}
	}
}

int main()
{
	freopen("test1.txt","r", stdin);
	//ifstream myfile ("digitdata.txt");
	//int row = 1001,col = 157;
	int row = 12,col = 2;
	countn = 0;
	mincost = INT_MAX;
	K=3;
	for(int i=0;i<row;i++)
	{
		vector<float> v1;
		for(int j=0;j<col;j++)
		{
			int a1;
			cin >> a1;
			//cout <<  i << " " << j << " "  << a1  << endl;
			v1.push_back((float)a1);
		}
		input_set.push_back(v1);
	}
	struct timeval start, end;
	gettimeofday(&start, NULL);
	DIM = input_set[0].size();
	N= 10;
	M= 3;
	//N = 400,M = 50;
	total_data_points = input_set.size();
	vector<int> index_of_n_sampled_points =  generate_k_random(total_data_points,N);
	//print(index_of_n_sampled_points);
	/*cout << v1.size() << endl;// " " << v1[0].size() << endl;
	for(int i =0;i<v1.size();i++)
		cout << v1[i] << " ";
	cout << endl;*/
	/*for(int i=0;i<row-1;i++)
	{
		for(int j=0;j<col;j++)
		{
			cout << i << " " << j << " " << input_set[i][j] << "\n " ;
		}
		cout << endl;
	}*/
	//return 0;

    for(int i=0;i<N;i++)
        vec_0_N.push_back(i);
    //generate_all_subsets(data_subset,vec_0_N.size(),M,0);

	//generate_subset_indices(data_subset,vec_total_subset.size(),K,0);
	cout << "\n" << all_subset_points.size() << "\n";
	//iterative_subset();
	//cout << "\n" << k_indices_of_subset_points.size() << "\n";
	//int *a = new int[K];
	//for(int i=0;i<K;i++)
	//	a[i] = 0;
	//iterative_subset_1(a,all_subset_points.size(),0);
	int i1;
	//#pragma omp parallel for num_threads(24)
	//for(i1=0;i1<all_subset_points.size();i1++)
	//{
		//iterative_subset_open_mp(i1,all_subset_points.size());
	//}

	cout << "Done" << endl;
    int size = vec_0_N.size();
    int i;
    #pragma omp parallel for num_threads(8)
    for(i=0;i<size;i++){
        vector<int> data;
        vector<vector<float> > centers;
        data.push_back(vec_0_N[i]);
        cout << "loop counter = "<< i << " thread id = "<< omp_get_thread_num()<< endl;
        subsets(data,centers,index_of_n_sampled_points,vec_0_N.size(),M-1,i+1);
    	//subsets(data,centers,index_of_n_sampled_points,size,left-1,i+1);
    	data.pop_back();
    }

    //for(i1=0;i1<k_indices_of_subset_points.size();i1++)
	//{
		//main_openmp(k_indices_of_subset_points[i1],index_of_n_sampled_points);
	//}
	/*for(int i=0;i<k_indices_of_subset_points.size();i++)
	{
		for(int j=0;j<k_indices_of_subset_points[i].size();j++)
		{
			cout << k_indices_of_subset_points[i][j] << " " ;
		}
		cout << endl;
	}
	for(int i=0;i<K;i++)
	{
		cout << a[i] << " ";
	}
	cout << endl;*/
	//subsets(data_subset,centers,index_of_n_sampled_points,vec_0_N.size(),M,0);
	gettimeofday(&end, NULL);

	float delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
         		end.tv_usec - start.tv_usec) / 1.e6;
    cout<<delta<<endl;
	cout << "final centers are" << endl;
	for(int i=0;i<opti_centers.size();i++)
	{
		print(opti_centers[i]);
	}

	return 0;
}
