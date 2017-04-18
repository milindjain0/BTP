 #include <bits/stdc++.h>
#include "omp.h"
#include <sys/time.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
using namespace std;
int DIM;
float epsilon;
int total_data_points;
float **input_set; //contains the input data points
float *input_set_gpu;
vector<int> vec_0_N;
int N,M,K;
float mincost;
vector<vector<float> > opti_centers;
float** final_centers;
vector<vector<int> > all_subset_points;
vector<vector<int> > k_indices_of_subset_points;


void print( vector<int> l){
    for(vector<int>::iterator it=l.begin(); it!=l.end() ; ++it)
            cout << " " << *it;
    cout<<endl;
}

void print( float *l){
    for(int i=0;i<DIM;i++)
        cout << l[i] << " ";
    cout << endl;
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

float calc_distance(vector<float> &v1,float* v2)
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

__global__ void kernel(float *data, int row, int col,int dim) {

  printf("Element (%d, %d) = %d\n", row, col, data[(row*dim)+col]);

}
__global__ void gpu_fun(float *v1,float* v2,float* tmp_mindist,int row,int dim)
{
        int i =threadIdx.x;
        int j = row*dim +i;
        float val =  v1[i] - v2[j];
        //printf("inside kernel %d\n",v2[j]);
        //printf(" %d %d\n", i,val*val);// << i << endl;
		tmp_mindist[i] =  val*val;
}


vector<int> d2sampling(vector<vector<float> > &centers)
{
	map<pair<float,int>, int> mymap;
	vector<int> ans;
	for(int i=0;i<total_data_points;i++)
	{
		float mindist = INT_MAX;
        //  print(input_set[i]);
		for(int j=0;j<centers.size();j++)
		{
            float *centers1 = new float[DIM];
            for(int k=0;k<DIM;k++)
                centers1[k] = centers[j][k];
            float* centers_gpu;
            float* mindist_gpu;
            float * mindist_tmp= new float[DIM];
            *mindist_tmp = 0;
            cudaMalloc(&mindist_gpu, sizeof(float)*DIM);
            cudaMemcpy(mindist_gpu, mindist_tmp, sizeof(float)*DIM, cudaMemcpyHostToDevice);

            const size_t a_size = sizeof(float) * DIM;
            cudaMalloc(&centers_gpu, a_size);
            cudaMemcpy(centers_gpu, centers1, a_size, cudaMemcpyHostToDevice);

            //mindist = min(calc_distance(centers[j],input_set[i]),mindist );
            gpu_fun<<<1,DIM>>>(centers_gpu,input_set_gpu,mindist_gpu,i,DIM);
            cudaMemcpy(mindist_tmp, mindist_gpu, sizeof(float)*DIM, cudaMemcpyDeviceToHost);
            float res=0;
            for(int k=0;k<DIM;k++)
                res+=mindist_tmp[k];
            //cout << "res " <<  res << " mindist " << calc_distance(centers[j],input_set[i])<<endl;
            mindist = min(res,mindist);
            cudaFree(centers_gpu);
            cudaFree(mindist_gpu);

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

__device__ static float atomicMin(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fminf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}


__global__ void cost(float* data,float *centers,float *ans,int dim,int num_centers) {

    //printf("%s\n", "here");
    int id = threadIdx.x;
    int blkid = blockIdx.x;
    float val=0;
    for(int i=0;i<dim;i++)
    {
         float tmp = data[id*dim+i]-centers[blkid*dim+i];
         val = val + tmp*tmp ;
    }
    atomicMin(&ans[id], val);
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
    float *gpu_ans;
    cudaMalloc((void**)&gpu_ans,sizeof(float)*total_data_points);
    float* dev_centers;
    cudaMalloc((void **)&dev_centers, K*DIM* sizeof(float));
	vector<vector<float> > centers;
    float **center_new  = new float*[K];
    for(int i = 0;i < K; i++)
        center_new[i] = new float[DIM];
	vector<float> mean;
	vector<int> data = all_subset_points[indices[0]];
    float* init_ans = new float[total_data_points];
    for(int i = 0 ; i < total_data_points ; i++)
        init_ans[i] = INT_MAX;
    float* ans = new float[total_data_points];
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
    for(int i=0;i<DIM;i++)
    {
        center_new[0][i] = mean[i];
    }
	centers.push_back(mean);
	for(int i=1;i<K;i++)
	{
		vector<float> mean1;
		vector<int> data1 = all_subset_points[indices[i]];
		for(int j=0;j<DIM;j++)
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
		for(int j = 0 ; j < data_size1;j++)
		{
			for(int k=0;k<DIM;k++)
			{
				mean1[k]+= ((input_set[index_of_n_sampled_points[data1[j]]][k] *1.0)/(data_size1));
			}
		}
		centers.push_back(mean1);
        for(int x=0;x<DIM;x++)
        {
            center_new[i][x] = mean1[x];
        }
	}
    float new_cost = 0;

    for(int i = 0 ; i < K ; i++)
        cudaMemcpy(dev_centers + i*DIM, center_new[i], DIM*sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_ans,init_ans,sizeof(float)*total_data_points,cudaMemcpyHostToDevice);
    cost<<<K,total_data_points>>>(input_set_gpu,dev_centers,gpu_ans,DIM,K);
    cudaMemcpy( ans, gpu_ans, total_data_points* sizeof(float), cudaMemcpyDeviceToHost );
    for(int i=0;i<total_data_points;i++)
    {
        new_cost += ans[i];
    }
    //cout << new_cost << "  " << cost(centers) << endl;
	#pragma omp critical
	{			//cout << "new cost " << centers.size() << endl;
		if(new_cost < mincost)
		{
			mincost = new_cost;
            for(int i=0;i<K;i++)
                for(int j=0;j<DIM;j++)
                    final_centers[i][j] = center_new[i][j];
			cout << mincost <<  " " << new_cost  << endl;
		}
	}
    cudaFree(dev_centers);
    cudaFree(gpu_ans);

}

int main()
{
	freopen("test1.txt","r", stdin);
	//ifstream myfile ("digitdata.txt");
	//int row = 1001,col = 157;
	int row = 12,col = 2;
    input_set = new float*[row];
    for(int i=0;i<row;i++)
        input_set[i] = new float[col];
	countn = 0;
	mincost = INT_MAX;
	K=3;
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			int a1;
			cin >> a1;
            input_set[i][j] = a1;
			//v1.push_back((float)a1);
		}
		//input_set.push_back(v1);
	}
    cudaMalloc((void **)&input_set_gpu, row*col* sizeof(float));
    for(int i = 0 ; i < row ; i++)    {
        // cudaMalloc((void **)&hd_array[i], length[i] * sizeof(int));
        cudaMemcpy(input_set_gpu + i*col, input_set[i], col*sizeof(float), cudaMemcpyHostToDevice);
    }



	struct timeval start, end;
	gettimeofday(&start, NULL);
	//DIM = input_set[0].size();
    DIM = col;
	N= 10;
	M= 3;
    final_centers = new float*[K];
    for(int i = 0 ; i < K ; i++)
    {
        final_centers[i] = new float[DIM];
    }

    total_data_points = row;
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
	vector<int> data_subset;

	generate_all_subsets(data_subset,vec_0_N.size(),M,0);

	//generate_subset_indices(data_subset,vec_total_subset.size(),K,0);
	cout << "\n" << all_subset_points.size() << "\n";
	//iterative_subset();
	//cout << "\n" << k_indices_of_subset_points.size() << "\n";
	//int *a = new int[K];
	//for(int i=0;i<K;i++)
	//	a[i] = 0;
	//iterative_subset_1(a,all_subset_points.size(),0);
	int i1;
	#pragma omp parallel for num_threads(8)
	for(i1=0;i1<all_subset_points.size();i1++)
	{
		iterative_subset_open_mp(i1,all_subset_points.size());
	}

	cout << "Done" << endl;
	#pragma omp parallel for num_threads(8)
	for(i1=0;i1<k_indices_of_subset_points.size();i1++)
	{
		main_openmp(k_indices_of_subset_points[i1],index_of_n_sampled_points);
	}
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
	for(int i=0;i<K;i++)
	{
		print(final_centers[i]);
	}
    cudaFree(input_set_gpu);

	return 0;
}
