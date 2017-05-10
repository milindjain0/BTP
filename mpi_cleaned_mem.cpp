#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/time.h>
#include "omp.h"

#define MASTER 0    /* taskid of first task */

using namespace std;

int DIM;    //dimension of the input dataset
int N,      //number of points sampled
    M,      //size of subset to be taken from N sampled points
    K,      //Number of cluseters
    total_data_points;      //total number of data points
float epsilon;
vector<vector<float> > input_set;   //contains the input data points

vector< vector<int> > all_subset_points;
vector<int> all_subset_points1;

vector<vector<int> > task_indices;
vector<int> task_indices1;

vector<int> vec_0_N;
float mincost;  //Cost of current optimal center
vector<vector<float> > opti_centers;    //Best centers found

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

//Calculate distance between two points.
float calc_distance(vector<float> &v1,vector<float> &v2)
{
	int sz = v1.size();
	float distsq = 0;
	for(int i=0;i<sz;i++)
	{
		distsq += ( (v1[i]-v2[i]) * (v1[i]-v2[i]) );
	}
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
		count++;
	}
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
	return cst;
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



void p_subset_master(int i)
{
	vector<int>* tmp1=new vector<int>();
    tmp1->push_back(i);
    int num_points = N;
	queue< vector<int>* > k_indices_of_subset_points_t;
	k_indices_of_subset_points_t.push(tmp1);
	while( !k_indices_of_subset_points_t.empty() && (k_indices_of_subset_points_t.front())->size() < M )
	{
			vector<int>* tmp=new vector<int>();
            vector<int> * tt = k_indices_of_subset_points_t.front();
            tmp->insert(tmp->end(),tt->begin(),tt->end());
			k_indices_of_subset_points_t.pop();
			for(int j=tmp->back()+1;j<num_points;j++)
			{

                int flag = 0;
                for(int k = 0;k <tmp->size(); k++)
                {
                    if(tmp->at(k) == j)
                    {
                        flag = 1;
                        break;
                    }
                }
                if(flag == 1)
                    continue;
                vector<int>* tmp2=new vector<int>();
				tmp->push_back(j);
                tmp2->insert(tmp2->end(),tmp->begin(),tmp->end());
				k_indices_of_subset_points_t.push(tmp2);
				tmp->pop_back();
			}
	}
	#pragma omp critical
	{
		while(!k_indices_of_subset_points_t.empty())
		{
            vector<int>* v = k_indices_of_subset_points_t.front();
            if(v->size() == M)
            {
                all_subset_points1.insert(all_subset_points1.end(),
                v->begin(),
                v->end());
            }
                //all_subset_points.push_back(*v);
			k_indices_of_subset_points_t.pop();
		}
	}
}

void generate_all_subsets(vector<int> data,int size, int left, int index)
{
	if(left==0){
			//all_subset_points1.push_back(data);
            all_subset_points1.insert(all_subset_points1.end(),data.begin(),data.end());

	}
    for(int i=index;i<size;i++){
    	data.push_back(vec_0_N[i]);
    	generate_all_subsets(data,size,left-1,i+1);
    	data.pop_back();
    }
}


map<vector<vector<float> >,vector<int> > dp;
void main_openmp(vector<int> indices,vector<int> index_of_n_sampled_points)
{


	vector<vector<float> > centers;
	vector<float> mean;
	vector<int> data = all_subset_points[indices.at(0)];

    // if(data[0] == 85 && data[1] == 10 && indices[2] == 0 )
    // {
    //      cout << "data is \n";
    //      print(data);
    // }
    for(int i=0;i<DIM;i++)
		mean.push_back(0.0);

	int data_size = data.size();
	for(int k = 0 ; k < data_size;k++)
	{
		for(int j=0;j<DIM;j++)
		{
			mean[j]+= ((input_set[index_of_n_sampled_points[data.at(k)]][j] *1.0)/(data_size));
		}
	}
	centers.push_back(mean);

	for(int i=1;i<K;i++)
	{
		vector<float> mean1;
		vector<int> data1 = all_subset_points[indices.at(i)];
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
			//#pragma omp critical
			//{
				dp[centers] = index_of_n_sampled_points;
			//}
		}
		for(int k = 0 ; k < data_size1;k++)
		{
			for(int j=0;j<DIM;j++)
			{
				mean1[j]+= ((input_set[index_of_n_sampled_points[data1.at(k)]][j] *1.0)/(data_size1));
			}
		}
		centers.push_back(mean1);
	}
	float new_cost = cost(centers);
    // if(data[0] == 1 && data[1] == 2 && data[2] == 3 && data[3] == 5 )
    // {
    //     std::cout << "cost at optimal is " << new_cost << '\n';
    //     std::cout << "centers are" << '\n';
    //     for(int i=0;i<centers.size();i++)
    //         print(centers[i]);
    //     std::cout << "centers are" << '\n';
    //
    // }
	//#pragma omp critical
	//{			//cout << "new cost " << centers.size() << endl;
    // if(new_cost < 12)
    //     cout << "solution found" << endl;
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
	//}
}


void subset_indices_master(int i,int num_points)
{
    //std::cout << "num_points= " << num_points << " i = " <<  i  << '\n';
	vector<int>* tmp1=new vector<int>();
    tmp1->push_back(i);
	queue< vector<int>* > k_indices_of_subset_points_t;
	k_indices_of_subset_points_t.push(tmp1);
	while( !k_indices_of_subset_points_t.empty() && (k_indices_of_subset_points_t.front())->size() < K )
	{
			vector<int>* tmp=new vector<int>();
            vector<int> * tt = k_indices_of_subset_points_t.front();
            tmp->insert(tmp->end(),tt->begin(),tt->end());
            //print(tt);
            //print(tmp);
			k_indices_of_subset_points_t.pop();
			//k_indices_of_subset_points.erase(k_indices_of_subset_points.begin());
			for(int j=0;j<num_points;j++)
			{
                vector<int>* tmp2=new vector<int>();
				tmp->push_back(j);
                tmp2->insert(tmp2->end(),tmp->begin(),tmp->end());
				k_indices_of_subset_points_t.push(tmp2);
				tmp->pop_back();
			}
	}
	#pragma omp critical
	{
		while(!k_indices_of_subset_points_t.empty())
		{
            vector<int>* v = k_indices_of_subset_points_t.front();
            if(v->size() == K)
                task_indices1.insert(task_indices1.end(),
                                        v->begin(),
                                        v->end());
			k_indices_of_subset_points_t.pop();
		}
	}
}


int main(int argc, char *argv[])
{
    int	   numtasks,       /* number of tasks in partition */
	       taskid;         /* a task identifier */
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    freopen("test1.txt","r", stdin);
	//ifstream myfile ("digitdata.txt");
	//int row = 1001,col = 157;
	int row = 12,col = 2;

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
    vector<int> index_of_n_sampled_points;
    int sz;
    int sz_all_subset_point;
    int size_per_process;
    int size_task_indices;
    for(int i=0;i<N;i++)
		vec_0_N.push_back(i);
    if(taskid == MASTER)
    {
        index_of_n_sampled_points =  generate_k_random(total_data_points,N);
        std::cout << "sampled points are" << '\n';
        print(index_of_n_sampled_points);
        sz = index_of_n_sampled_points.size();
    }
    else
    {
        sz = N;
        index_of_n_sampled_points.resize(N);
    }
    MPI_Bcast(
        &index_of_n_sampled_points[0],
        N,
        MPI_INT,
        MASTER,
        MPI_COMM_WORLD);

    // for(int i= 0; i < numtasks; i++)
    // {
    //     if(taskid == i)
    //     {
    //         cout << "my id = " << taskid << endl;
    //         print(index_of_n_sampled_points);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    if(taskid == MASTER)
    {
        int i1;
        #pragma omp parallel for num_threads(24)
    	for(i1=0;i1<N;i1++)
    	{
    		p_subset_master(i1);
        }
        sz_all_subset_point = all_subset_points1.size()/M;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(
        &sz_all_subset_point,
        1,
        MPI_INT,
        MASTER,
        MPI_COMM_WORLD);

    if(taskid != MASTER)
            all_subset_points1.resize(sz_all_subset_point*M);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(
        &all_subset_points1[0],
        sz_all_subset_point*M,
        MPI_INT,
        MASTER,
        MPI_COMM_WORLD);

    if (taskid == MASTER)
    {

        //Task Division
        cout << "size all_subset_points =  " << sz_all_subset_point << endl;
        int i2;
        #pragma omp parallel for num_threads(24)
    	for(i2=0;i2<sz_all_subset_point;i2++)
    	{
    		subset_indices_master(i2,sz_all_subset_point);
    	}

        size_task_indices = task_indices1.size()/K;


        size_per_process = (int)ceil(size_task_indices*1.0/numtasks);
        cout << "size_per_process =  " << size_per_process << endl;
        cout << "num task =  " << numtasks << endl;
        cout << "size task indices =  " << task_indices1.size()/K << endl;
        for(int i=1;i<numtasks;i++)
        {
            int send_size = min(size_per_process*K,size_task_indices*K - i*K*size_per_process);
            MPI_Send(
                &send_size,
                1,
                MPI_INT,
                i,
                MASTER,
                MPI_COMM_WORLD);
            MPI_Send(
                &task_indices1[i*K*size_per_process],
                send_size,
                MPI_INT,
                i,
                MASTER,
                MPI_COMM_WORLD);
        }

    }
    else
    {
        MPI_Recv(
            &size_task_indices,
            1,
            MPI_INT,
            MASTER,
            MASTER,
            MPI_COMM_WORLD,
            &status);
        cout << "size_task_indices slaves " << size_task_indices << endl;
        task_indices1.resize(size_task_indices);
        size_task_indices = size_task_indices/K;
        //printf("rank %d receiving received %d\n", taskid, sz_all_subset_point);
        MPI_Recv(
            &task_indices1[0],
            size_task_indices*K,
            MPI_INT,
            MASTER,
            MASTER,
            MPI_COMM_WORLD,
            &status);
        //MPI_Barrier(MPI_COMM_WORLD);
    }

    for(int i=0;i<sz_all_subset_point;i++)
    {
        vector<int> tmp;
        for(int j=0;j<M;j++)
        {
            tmp.push_back(all_subset_points1[i*M+j]);
        }
        if(tmp[0] == 1 && tmp[1] == 2 && tmp[2] == 3 && tmp[3] == 5 )
            cout << "found in " << taskid << endl;
        all_subset_points.push_back(tmp);
    }



    if(taskid == MASTER)
    {
        for(int i=0;i<size_per_process;i++)
        {
            vector<int> tmp;
            for(int j=0;j<K;j++)
            {
                tmp.push_back(task_indices1[i*K+j]);
            }
            task_indices.push_back(tmp);
        }
        //  for(int i=0;i<size_per_process;i++)
        //  {
        //      for(int j=0;j<K;j++)
        //          cout << task_indices[i][j] << " ";
        //      cout << endl;
        //  }
        int i1;
        for(i1=0;i1<size_per_process;i1++)
        {
            main_openmp(task_indices[i1],index_of_n_sampled_points);
        }
    }
    else
    {
        for(int i=0;i<size_task_indices;i++)
        {
            vector<int> tmp;

            for(int j=0;j<K;j++)
            {
                tmp.push_back(task_indices1[i+j]);
            }
            if(tmp[0] == 85 && tmp[1] == 18 && tmp[2] == 0)
                cout << "found in slave" << endl;
            task_indices.push_back(tmp);
        }
        int i1;
        for(i1=0;i1<size_task_indices;i1++)
        {
            main_openmp(task_indices[i1],index_of_n_sampled_points);
        }

    }

    float localres[2];
    float globalres[2];
    localres[0] = mincost;


    localres[1] = (float)taskid;

    MPI_Allreduce(localres, globalres, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);



    if(taskid == (int)globalres[1])
    {
        cout << globalres[0] << " ,  " << globalres[1] << endl;
        for(int i=0;i<opti_centers.size();i++)
        {
            for(int j=0;j<opti_centers[i].size();j++)
            {
                cout << opti_centers[i][j] << " ";
            }
            cout << endl;
        }
    }
    MPI_Finalize();
    return 0;
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

	vector<int> data_subset;
	vector<vector<float> > centers;
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
