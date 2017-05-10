#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/time.h>
#include "omp.h"

#define MASTER 0    /* taskid of first task */

using namespace std;
float* heuristic1;
int DIM;    //dimension of the input dataset
int N,      //number of points sampled
    M,      //size of subset to be taken from N sampled points
    K,      //Number of cluseters
    total_data_points;      //total number of data points
float epsilon;
vector<vector<float> > input_set;   //contains the input data points

vector< vector<int> > all_subset_points;
vector<int> all_subset_points1;
vector<int> all_subset_points2;

vector<vector<int> > task_indices;
vector<int> task_indices1;

vector<int> vec_0_N;
float mincost;  //Cost of current optimal center
float maxrad;
vector<vector<float> > opti_centers;    //Best centers found

float **cost_map;
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

float maxdistfromset(vector<int>& index_of_n_sampled_points,vector<int>& m){
    float thread_maxrad[4]={0};
    int i1;
    #pragma omp parallel for num_threads(4)

    for(i1=0;i1<m.size();i1++)
    {
        for(int j1=i1+1;j1<m.size();j1++)
        {
            std::vector<float> v1 = input_set[index_of_n_sampled_points[m[i1]]];
            std::vector<float> v2 = input_set[index_of_n_sampled_points[m[j1]]];
            thread_maxrad[omp_get_thread_num()] = max(calc_distance(v1,v2),thread_maxrad[omp_get_thread_num()]);
        }
    }
    float maxrad1 = 0;
    for(int i=0;i<4;i++)
    {
        maxrad1= max(maxrad1,thread_maxrad[i]);
    }
    return maxrad1;
}

vector<int> d2sampling(vector<vector<float> > &centers)
{
	map<pair<float,int>, int> mymap;
	vector<int> ans;
    int i;
    #pragma omp parallel for num_threads(4)
	for(i=0;i<total_data_points;i++)
	{
		float mindist = FLT_MAX;
		for(int j=0;j<centers.size();j++)
		{
			mindist = min(calc_distance(centers[j],input_set[i]),mindist );
		}
        #pragma omp critical
        {
            mymap.insert(make_pair(make_pair(mindist,i),0 ));
        }
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
    int i;
    //cout << "in cost \n";
    float cost_per_thread[4] = {};
    //cout << "size of centers " << centers.size() << endl;
    #pragma omp parallel for num_threads(4)
	for(i=0;i<total_data_points;i++)
	{
		float mindist = FLT_MAX;
		for(int j=0;j<centers.size();j++)
		{
            float r;
            r =calc_distance(centers[j],input_set[i]);
			mindist = min(r,mindist );
            //cout << "md= " << mindist << endl;
		}
        cost_per_thread[omp_get_thread_num()]+=mindist;
		//cst += mindist;
        //cout << "cst = " << cst << endl;
	}
    //cout << "77t888 cost function " << omp_get_num_threads() << endl;
    for(int i=0;i<4;i++)
        cst+= cost_per_thread[i];
    //std::cout << " cost = " << cst << '\n';
	return cst;
}

vector<int> generate_k_random(int n, int k)
{
	vector<int> v,v1;
	epsilon = 1/20;
	for(int i=0;i<n;i++)
		v1.push_back(i);
    //std::random_shuffle(v1.begin(),v1.end()) ;
	for(int i=0;i<k;i++)
	{
		int randint = rand()%(n-i)+i;
		v.push_back(v1[randint]);
		iter_swap(v1.begin()+i,v1.begin()+randint);
	}
	return v;
}



void p_subset_master(int i,std::vector<int>& index_of_n_sampled_points)
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
            tmp->erase(tmp->begin(),tmp->end());
            tt->erase(tt->begin(),tt->end());
	}
	#pragma omp critical
	{
		while(!k_indices_of_subset_points_t.empty())
		{
            vector<int>* v = k_indices_of_subset_points_t.front();
            if(v->size() == M)
            {
                #ifndef random1
                int randint = rand()%100;
                if(randint == 1){
                    float md=maxdistfromset(index_of_n_sampled_points,*v);
                    if(md<maxrad/K)
                    all_subset_points2.insert(all_subset_points2.end(),
                    v->begin(),
                    v->end());
                }
                #else
                float md=maxdistfromset(index_of_n_sampled_points,*v);
                if(md<maxrad/K)
                all_subset_points2.insert(all_subset_points2.end(),
                v->begin(),
                v->end());
                #endif
            }
                //all_subset_points.push_back(*v);
            v->erase(v->begin(),v->end());
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

    // for(int i=0;i<total_data_points;i++)
    // {
    //     for(int j=0;j<K;j++)
    //         cost_map[i][j] = -1;
    // }
    //MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    //cout << " taskid= " << taskid  << endl;
	vector<vector<float> > centers;
	vector<float> mean;
	vector<int>* data = &all_subset_points[indices.at(0)];
    for(int i=0;i<DIM;i++)
		mean.push_back(0.0);

	int data_size = data->size();
	for(int k = 0 ; k < data_size;k++)
	{
		for(int j=0;j<DIM;j++)
		{
			mean[j]+= ((input_set[index_of_n_sampled_points[data->at(k)]][j] *1.0)/(data_size));
		}
	}
	centers.push_back(mean);
    //h1[0]=cost(centers);
//    if(h1[0]*0.5>heuristic1[0]){
        //cout  << "here0" <<endl;
//        return;
//    }

    // if(cost(centers)*0.5 > mincost)
    // {
    //     cout << "here" << endl;
    //     return;
    // }
	for(int i=1;i<K;i++)
	{
		vector<float> mean1;
		vector<int>* data1 = &all_subset_points[indices.at(i)];
   		for(int j=0;j<DIM;j++)
			mean1.push_back(0.0);
		int data_size1 = data1->size();

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
				mean1[j]+= ((input_set[index_of_n_sampled_points[data1->at(k)]][j] *1.0)/(data_size1));
			}
		}
		centers.push_back(mean1);
        //h1[i]=cost(centers);
        //if(h1[i]*0.5>heuristic1[i]){
            //cout << "here " << i << endl;
            //return;
        //}

	}
	float new_cost = cost(centers);
	if(new_cost < mincost)
	{
        cout << "mincost= " << mincost <<  " new_cost= " <<  new_cost  << endl;
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
        //for(int i=0;i<K;i++)
            //heuristic1[i] = h1[i];
	}
	//}
}


void subset_indices_master(int i,int num_points,vector<int> index_of_n_sampled_points)
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
            tmp->erase(tmp->begin(),tmp->end());
            tt->erase(tt->begin(),tt->end());

	}
//	#pragma omp critical
//	{
		while(!k_indices_of_subset_points_t.empty())

		{
            vector<int> v = *(k_indices_of_subset_points_t.front());
            if(v.size() == K)
            {
                //cout << "here" << endl;
                #ifndef random2
                int randint = rand()%100;
                if(randint == 1)
                    main_openmp(v,index_of_n_sampled_points);
                #else
                    main_openmp(v,index_of_n_sampled_points);
                #endif
            }
                //task_indices1.insert(task_indices1.end(),
                    //                    v->begin(),
                    //                    v->end());
            (k_indices_of_subset_points_t.front())->erase((k_indices_of_subset_points_t.front())->begin(),(k_indices_of_subset_points_t.front())->end());
			k_indices_of_subset_points_t.pop();
		}
	//}
}


int main(int argc, char *argv[])
{
    srand (time(NULL));
    int	   numtasks,       /* number of tasks in partition */
	       taskid;         /* a task identifier */
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    freopen(argv[1],"r", stdin);
	int row = atoi(argv[2]),col =atoi(argv[3]);

	mincost = FLT_MAX;
	K=atoi(argv[4]);
    heuristic1=new float[K];
    for(int i=0;i<K;i++)
        heuristic1[i] = FLT_MAX;
    for(int i=0;i<row;i++){
		vector<float> v1;
		for(int j=0;j<col;j++)
		{
			float a1;
			cin >> a1;
            //cout <<  i << " " << j << " "  << a1  << endl;
			v1.push_back((float)a1);
		}
        //print(v1);
		input_set.push_back(v1);
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	DIM = input_set[0].size();
	N=  atoi(argv[5]);
	M=  atoi(argv[6]);
	//N = 400,M = 50;
	total_data_points = input_set.size();
    //cost_map = new float*[total_data_points];
    maxrad = 0;
    int i3;
    float thread_maxrad[4]={};
    #pragma omp parallel for num_threads(4)
    for(i3=0;i3<total_data_points;i3++)
    {
        for(int j1=i3+1;j1<total_data_points;j1++)
        {
            thread_maxrad[omp_get_thread_num()] = max(calc_distance(input_set[i3],input_set[j1]),thread_maxrad[omp_get_thread_num()]);
        }
    }
    for(int i=0;i<4;i++)
    {
        maxrad = max(maxrad,thread_maxrad[i]);
    }
    // for(int i=0;i<total_data_points;i++)
    // {
    //     cost_map[i] = new float[K];
    // }
    vector<int> index_of_n_sampled_points;
    int sz;
    int sz_all_subset_point;
    int size_per_process;
    int size_task_indices;
    for(int i=0;i<N;i++)
		vec_0_N.push_back(i);
    if(taskid == MASTER)
    {
        index_of_n_sampled_points =  generate_k_random(total_data_points,N);;
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

    // for(int i= 0; i < numtasks; i{}++)
    // {
    //     if(taskid == i)
    //     {
    //         cout << "my id = " << taskid << endl;
    //         print(index_of_n_sampled_points);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

        int subset_per_process = N/numtasks;
        if(subset_per_process==0)
            subset_per_process = 1;
        int i1;
        int loop_end;
        if(taskid == numtasks-1)
            loop_end = N;
        else
            loop_end = (taskid+1)*subset_per_process;
        #pragma omp parallel for num_threads(4)
    	for(i1=taskid*subset_per_process;i1<loop_end;i1++)
    	{
    		p_subset_master(i1,index_of_n_sampled_points);
        }
        int sz_all_subset_point2 = all_subset_points2.size()/M;
        int slaves_send_size = sz_all_subset_point2*M;
        cout << "taskid=  " <<  taskid  << "slaves_send_size "<< slaves_send_size << endl;
        int *recvcounts = new int[numtasks];
        MPI_Gather(
        &slaves_send_size,
        1,
        MPI_INT,
        &recvcounts[0],
        1,
        MPI_INT,
        MASTER,
        MPI_COMM_WORLD);

    int *displs = new int[numtasks];
    int cnt_all_subset_point=0;
    if(taskid ==MASTER)
    {
        int cnter=0;
        for(int i=0;i<numtasks;i++)
        {
            displs[i]=cnter;
            cnter+=(recvcounts[i]);
            cnt_all_subset_point+=recvcounts[i];
        }
    }
    // if(taskid == MASTER)
    // {
    //     cout << "displs" << endl;
    //     for(int i=0;i<numtasks;i++)
    //     cout<< displs[i] << " ";
    //     cout << endl;
    // }
    cout << "cnt_all_subset_point= " << cnt_all_subset_point << endl;
    if(taskid == MASTER)
            all_subset_points1.resize(cnt_all_subset_point);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gatherv(
    &all_subset_points2[0],
    slaves_send_size,
    MPI_INT,
    &all_subset_points1[0],
    &recvcounts[0],
    &displs[0],
    MPI_INT,
    MASTER,
    MPI_COMM_WORLD);
    sz_all_subset_point = all_subset_points1.size()/M;
    MPI_Barrier(MPI_COMM_WORLD);
    print(all_subset_points1);
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

        //Task Division
        //cout << "size all_subset_points =  " << sz_all_subset_point << endl;

        int task_per_process = sz_all_subset_point/numtasks;
        if(task_per_process ==0)
            task_per_process=1;
        int i2;

        int loop_end1;
        if(taskid == numtasks-1)
            loop_end1 = sz_all_subset_point;
        else
            loop_end1 = (taskid+1)*task_per_process;
        for(int i=0;i<numtasks;i++)
        {
            if(taskid == i)
            {
                cout << "taskid  " << taskid<< endl;
                cout << "loop_start " << taskid*task_per_process<< endl;
                cout << "loop_end " << loop_end1<< endl;

                cout << endl;
                cout << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);


        }

        //size_task_indices = task_indices1.size()/K;

    // for(int i=0;i<numtasks;i++)
    // {
    //     if(taskid == i)
    //     {
    //         print(task_indices1);
    //         cout << endl;
    //         cout << endl;
    //         cout << endl;
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    //
    //
    // }

    all_subset_points2.erase(all_subset_points2.begin(),all_subset_points2.end());

    for(int i=0;i<sz_all_subset_point;i++)
    {
        vector<int> tmp;
        for(int j=0;j<M;j++)
        {
            tmp.push_back(all_subset_points1[i*M+j]);
        }
        all_subset_points.push_back(tmp);

    }
    all_subset_points1.erase(all_subset_points1.begin(),all_subset_points1.end());
    cout << "before subset_indices_master" << endl;
    //#pragma omp parallel for num_threads(4)
    for(i2=taskid*task_per_process;i2<loop_end1;i2++)
    {
        subset_indices_master(i2,sz_all_subset_point,index_of_n_sampled_points);
    }
    cout << "after subset_indices_master" << endl;

    //if(taskid == MASTER)
    //{
        // for(int i=0;i<size_task_indices;i++)
        // {
        //     vector<int> tmp;
        //     for(int j=0;j<K;j++)
        //     {
        //         tmp.push_back(task_indices1[i*K+j]);
        //     }
        //     task_indices.push_back(tmp);
        // }

        //  for(int i=0;i<size_per_process;i++)
        //  {
        //      for(int j=0;j<K;j++)
        //          cout << task_indices[i][j] << " ";
        //      cout << endl;
        //  }
        // int i1;
        // for(i1=0;i1<size_task_indices;i1++)
        // {
        //     main_openmp(task_indices[i1],index_of_n_sampled_points);
        // }
    //}
    // else
    // {
    //     for(int i=0;i<size_task_indices;i++)
    //     {
    //         vector<int> tmp;
    //
    //         for(int j=0;j<K;j++)
    //         {
    //             tmp.push_back(task_indices1[i+j]);
    //         }
    //         task_indices.push_back(tmp);
    //     }
    //     int i1;
    //     for(i1=0;i1<size_task_indices;i1++)
    //     {
    //         main_openmp(task_indices[i1],index_of_n_sampled_points);
    //     }
    //
    // }
    task_indices1.erase(task_indices1.begin(),task_indices1.end());

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
    MPI_Barrier(MPI_COMM_WORLD);
    if(taskid == MASTER)
    {
        gettimeofday(&end, NULL);
        float delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
        end.tv_usec - start.tv_usec) / 1.e6;
        cout << endl << "TIME TAKEN = " << delta << " sec" << endl;
        //cout<<delta<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(taskid==MASTER)
    {
        cout << "row = " << row << endl;
        cout << "col = " << col << endl;
        cout << "N = " << N << endl;
        cout << "M = " << M << endl;
        cout << "K = " << K << endl;
    }
    MPI_Finalize();
    return 0;
}
