#include <bits/stdc++.h>
using namespace std;
int DIM;
float epsilon;
int total_data_points;
vector<vector<float> > input_set; //contains the input data points
vector<int> vec_0_N;
int N,M,K;
float mincost;
vector<vector<float> > opti_centers;
/*vector<int> helper(vector<vector<int> > list_points,vector<vector<int> > centers);
{
	int num_centers = centers.size();
}*/
/*void print( list<int> l){
    for(list<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}*/
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
int main()
{
	freopen("test.txt","r", stdin);
	//ifstream myfile ("digitdata.txt");
	//int row = 1001,col = 157;
	int row = 15,col = 2;
	
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
	/*for(int i=0;i<row;i++)
	{
		vector<float> v1;
		for(int j=0;j<=col;j++)
			{
			if(i == 0)
			{
				string a;
				cin >> a;
			}
			else
			{
				if(j==0)
				{
					string a;
					cin >> a;
				}
				else
				{
					int a1;
					cin >> a1;
					cout <<  i << " " << j << " "  << a1  << endl;
					v1.push_back((float)a1);
				}
			}
		}
		if(i!=0)
			input_set.push_back(v1);
	}*/
	//return 0;
	DIM = input_set[0].size();
	N= 12;
	M= 4;
	//N = 400,M = 50;
	total_data_points = input_set.size();
	vector<int> index_of_n_sampled_points =  generate_k_random(total_data_points,N);
	print(index_of_n_sampled_points);
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
	vector<int> data_subset;
	vector<vector<float> > centers;
	subsets(data_subset,centers,index_of_n_sampled_points,vec_0_N.size(),M,0);
	cout << "final centers are" << endl;
	for(int i=0;i<opti_centers.size();i++)
	{
		print(opti_centers[i]);
	}

	return 0;
}