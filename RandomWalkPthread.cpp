/*
 *@description - This file implements the random walk algorithm using
 	 	 	 	 pthread to approximate the PageRank algorithm.
 	 	 	 	 It assumes the unweighted web graph file is in coordinate
                 format, i.e., each line of the input file should be in the
                 format of
                 <j> <i>
                 where j is the destination node, i is the source node.
                 The node numbers are continuous intergers starting from 0.
 *@author: Yao Zhu (yzhucs@gmail.com).
 */

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <set>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

using namespace std;

// deallocate a pointer array.
template<typename T>
void deallocate_pointer_array(T **ptr, int len) {
	if (ptr == NULL) {
		return;
	}
	for (int i = 0; i < len; i++) {
		if (ptr[i] != NULL) {
			delete [] ptr[i];
			ptr[i] = NULL;
		}
	}
}

/*
 *@description - read in the webgraph data. (src[m], dest[m]) is a directed edge
 	 	 	 	 in the webgraph. This function also returns the maximum node
 	 	 	 	 number, which is number of nodes - 1 (we assume the node
 	 	 	 	 number starts from 0).
 *@param - webgraph_file, name of the webgraph file.
 *@param[out] - src, the array of source nodes.
 *@param[out] - dest, the array of the destination nodes.
 *@param[out] - nnz, number of edges.
 *@return - the maximum node number.
 */
int read_webgraph(const char *webgraph_file, int **src, int **dest, int *nnz) {
	std::ifstream wbfile(webgraph_file);
	string line;
	*nnz = 0;
	// count the number of edges.
	while (std::getline(wbfile, line)) {
		// check the 1st character of this line.
		if (line[0] >= '0' && line[0] <= '9') {
			(*nnz)++;
		}
	}
	// allocate storage in coordinate format.
	*src = new int[*nnz];
	*dest = new int[*nnz];
	if (*src == NULL || *dest == NULL) {
		cerr << "cannot allocate src or dest in read_webgraph()." << endl;
		exit(-1);
	}
	// return to the beginning of the input file.
	wbfile.clear();
	wbfile.seekg(0, ios::beg);
	int i = 0;
	int max_nn = -1;	// maximum node number.
	while (std::getline(wbfile, line)) {
		// check the 1st character of this line.
		if (line[0] >= '0' && line[0] <= '9') {
			std::istringstream iss(line);
			iss >> (*dest)[i] >> (*src)[i];
			if ((*src)[i] > max_nn) {
				max_nn = (*src)[i];
			}
			if ((*dest)[i] > max_nn) {
				max_nn = (*dest)[i];
			}
			i++;
		}
	}
	wbfile.close();
	return max_nn;
}

/*
 *@description - construct the CSR format of the sparse PageRank matrix
 	 	 	 	 (src, dest) from the coordinate format given by (src, dest).
 *@param - (src, dest) defines the web graph edges.
 *@param - nnz, number of nonzero elements. it's the length of
 	 	   src, dest, and col_ind.
 *@param - N, number of nodes, it's the length of row_ptr and out_degree.
 *@param - (col_ind, row_ptr) represents the CSR format of the
 	 	   matrix (src, dest).
 */
void coord2csr(const int* const src, const int* const dest,
			   const int nnz, const int N, int *col_ind, int *row_ptr) {
	if (src == NULL || dest == NULL) {
		cerr << "none of src and dest can be NULL in coord2csr()." << endl;
		exit(-1);
	}
	if (col_ind == NULL || row_ptr == NULL) {
		cerr << "none of col_ind, and row_ptr can be NULL in "
				"coord2csr()." << endl;
		exit(-1);
	}
	// we need the out degree to construct the CSR format of (src, dest) matrix.
	int *degree = new int[N];
	std::fill(degree, degree + N, 0);
	for (int l = 0; l < nnz; l++) {
		int i = src[l];
		degree[i]++;
	}
	// compute row_ptr as the cumsum of degree.
	// note row_ptr[N] = nnz. the node numbers are in [0..N-1].
	row_ptr[0] = 0;
	for (int i = 1; i < N+1; i++) {
		row_ptr[i] = row_ptr[i-1] + degree[i-1];
	}
	// construct col_ind according to row_ptr.
	for (int l = 0; l < nnz; l++) {
		int i = src[l];
		int j = dest[l];
		col_ind[row_ptr[i]] = j;
		row_ptr[i]++;
	}
	// recompute row_ptr as the cumsum of degree.
	row_ptr[0] = 0;
	for (int i = 1; i < N+1; i++) {
		row_ptr[i] = row_ptr[i-1] + degree[i-1];
	}
	// deallocate degree.
	if (degree != NULL) {
		delete [] degree;
	}
}

/*
 *@description - get the id of the thread that has the node with the given
 	 	 	 	 node number. This partition assumes all the remnant elements
 	 	 	 	 are given to the last thread.
 *@param - nn, node number.
 *@param - N, total number of nodes.
 *@param - nthread, number of threads.
 */
int nn2tid(const int nn, const int N, const int nthread) {
	int quota = int(N / nthread);
	if (nn >= quota * nthread) {
		return (nthread - 1);
	} else {
		return int(nn / quota);
	}
}

/*
 *@description - get the number of rows a thread t has from the row partition
 				 scheme.
 *@param - N, total number of rows (one row correspond to one node).
 *@param - nthread, total number of processes.
 *@param - tid,  id of the thread.
 */
int get_Nlocal(const int N, const int nthread, const int tid) {
	int quota = int(N / nthread);
	if (tid == nthread - 1) {	// the last thread.
		return (N - (nthread - 1) * quota);
	} else {
		return quota;
	}
}

/*
 *@description - get the first node number that belongs to a thread t.
 */
int get_start_nn(const int N, const int nthread, const int tid) {
	int quota = int(N / nthread);
	return (quota * tid);
}

// print out_walk for test. we use a lock to synchronize the output from
// different threads.
pthread_mutex_t print_lock;		// lock for printing purpose.
void print_out_walk(const int* const out_walk, const int start_nn,
					const int Nlocal, int tid) {
	pthread_mutex_lock(&print_lock);
	cerr << "thread: " << tid << "------" << endl;
	for (int i = start_nn; i < start_nn + Nlocal; i++) {
		cerr << "(" << i << ", " << out_walk[i] << ") " << endl;
	}
	pthread_mutex_unlock(&print_lock);
}

//---------global variables shared by all threads.
#define TOP_NUM		100
int N;								// total number of nodes.
int num_threads;					// total number of threads.
int num_iter;						// number of desired steps.

// define a structure for data sharing among threads.
class thread_share_t {
public:
	int tid_;
	int *col_ind_;
	int *row_ptr_;
	int *out_walk_;
	int **in_walk_;
	pthread_mutex_t *in_walk_lock_;
public:
	thread_share_t() {
		tid_ = -1;
		col_ind_ = NULL;
		row_ptr_ = NULL;
		out_walk_ = NULL;
		in_walk_ = NULL;
		in_walk_lock_ = NULL;
	};
};

pthread_barrier_t in_walk_barr;
pthread_barrier_t out_walk_barr;
pthread_barrier_t timing_barr;
double tElapsed = 0.0;

// space used for merging local top ranked nodes into a globally top
// ranked nodes.
int *local_top_list_num_marks = NULL;	// of size num_threads*TOP_NUM. thread
										// t stores its results at
										// [TOP_NUM*t..TOP_NUM*(t+1)-1]
int *local_top_list_nn = NULL;
int *local_top_list_pointer = NULL;		// of size num_threads,
										// local_top_list_pointer[t] should be
										// in the range of [TOP_NUM*t..TOP_NUM*(t+1)-1]
int global_top_list_num_marks[TOP_NUM];
int global_top_list_nn[TOP_NUM];

// the class definition for using partial_sort().
class node_mark {
public:
    int num_marks_;
    int nn_;
};
bool comp_node_mark(const node_mark& nm_i, const node_mark& nm_j) {
	return (nm_i.num_marks_ > nm_j.num_marks_);
}

// this function returns the winner thread id according to the score array.
// the score of thread t is at location score[thread_pointer[t]].
// the winner achieves the maximal score.
int get_winner_tid(const int* const score, const int* const thread_pointer,
				   int num_threads) {
	if (score == NULL || thread_pointer == NULL) {
		cerr << "none of score and thread_pointer can be NULL in "
				"get_winner_tid()." << endl;
		exit(-1);
	}
	int max_score = INT_MIN;
	int winner_tid = -1;
	for (int tid = 0; tid < num_threads; tid++) {
		if (score[thread_pointer[tid]] > max_score) {
			max_score = score[thread_pointer[tid]];
			winner_tid = tid;
		}
	}
	if (winner_tid == -1) {
		cerr << "error in get_winner_tid()." << endl;
		exit(-1);
	} else {
		return winner_tid;
	}
}

/*
 *@description - this function randomly selects a neighbor for node i.
 *@param - i, the current node number.
 *@param - (col_ind, row_ptr) represents the CSR format of the graph. the
 	 	   neighbors of i are {col_ind[row_ptr[i]]..col_ind[row_ptr[i+1]-1]}
 *@param[out] - seed, the seed for the random number generator rand_r.
 	 	 	 	seed will be updated through each call to rand_r().
 */
int get_random_neighbor(const int i, const int* const col_ind,
						const int* const row_ptr, unsigned int *seed) {
	if (col_ind == NULL || row_ptr == NULL) {
		cerr << "none of col_ind and row_ptr can be NULL in "
				"get_random_neighbor()." << endl;
		exit(-1);
	}
	if (seed == NULL) {
		cerr << "seed cannot be NULL in get_random_neighbor()." << endl;
		exit(-1);
	}
	if (row_ptr[i+1] - row_ptr[i] == 0) {
		// a dangling node.
		return i;
	} else {
		int displs = floor(float(rand_r(seed)) / RAND_MAX *
						  (row_ptr[i+1]+1-row_ptr[i]));
		displs += row_ptr[i];
		if (displs == row_ptr[i+1] || displs == row_ptr[i+1] + 1) {
			return i;
		} else {
			return col_ind[displs];
		}
	}
}

/*
 *@description - the function implements the thread for advancing the random
 				 walks originating from the nodes a thread has.
 *@param - threadid, id of the thread.
 */
void* random_walk_thread(void *p_data) {
	thread_share_t *p_share_data = static_cast<thread_share_t*>(p_data);
	// unpack the shared data.
	int tid = p_share_data->tid_;				// record its thread id.
	int *col_ind = p_share_data->col_ind_;
	int *row_ptr = p_share_data->row_ptr_;
	int *out_walk = p_share_data->out_walk_;
	int **in_walk = p_share_data->in_walk_;
	pthread_mutex_t *in_walk_lock = p_share_data->in_walk_lock_;

	if (col_ind == NULL || row_ptr == NULL || out_walk == NULL ||
		in_walk == NULL || in_walk_lock == NULL) {
		cerr << "none of col_ind, row_ptr, out_walk, in_walk, in_walk_lock can "
				"be NULL in random_walk_thread() for thread " << tid << endl;
		exit(-1);
	}

	// check the availability of global variables.
	if (local_top_list_num_marks == NULL || local_top_list_nn == NULL ||
		local_top_list_pointer == NULL) {
		cerr << "none of local_top_list_num_marks, local_top_list_nn, "
				"local_top_list_pointer can be NULL in random_walk_thread() "
				"for thread " << tid << endl;
		exit(-1);
	}

	unsigned int seed = (unsigned int)tid;	// local seed for rand_r().

	// get the starting node belong to this thread.
	int start_nn = get_start_nn(N, num_threads, tid);
	// get the number of nodes belong to this thread.
	int Nlocal = get_Nlocal(N, num_threads, tid);

	struct timeval start_tv, end_tv;
	// barrier for time profiling.
	pthread_barrier_wait(&timing_barr);
	if (tid == 0) {
		gettimeofday(&start_tv, NULL);
	}

	// advance the random walks for num_iter steps.
	for (int iter = 0; iter < num_iter; iter++) {
		// for each node this thread possesses.
		for (int i = start_nn; i < start_nn + Nlocal; i++) {
			for (int cnt = 0; cnt < out_walk[i]; cnt++) {
				// get a random neighbor as the destination for one walk.
				int j = get_random_neighbor(i, col_ind, row_ptr, &seed);
				// try to update in_walk[j].
				pthread_mutex_lock(&in_walk_lock[j]);
				in_walk[j][0]++;
				pthread_mutex_unlock(&in_walk_lock[j]);
			}
		}
		// barrier to wait for all work updating in_walk complete.
	    int rc = pthread_barrier_wait(&in_walk_barr);
	    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	    {
	        cerr << "could not wait on in_walk_barrier in "
	        		"random_walk_thread()." << endl;
	        exit(-1);
	    }

		// copy in_walk to out_walk and zero out in_walk.
	    for (int i = start_nn; i < start_nn + Nlocal; i++) {
	    	out_walk[i] = in_walk[i][0];
	    	in_walk[i][0] = 0;
	    }

		// barrier to wait for all work copying in_walk complete.
	    rc = pthread_barrier_wait(&out_walk_barr);
	    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	    {
	        cerr << "could not wait on out_walk_barrier in "
	        		"random_walk_thread()." << endl;
	        exit(-1);
	    }
	}

	// use partial sort to find the local top TOP_NUM nodes.
	vector<node_mark> nm_list;
	nm_list.resize(Nlocal);
	node_mark curr_nm;
	for (int i = start_nn; i < start_nn + Nlocal; i++) {
		curr_nm.num_marks_ = out_walk[i];
		curr_nm.nn_ = i;
		nm_list.push_back(curr_nm);		// push_back is copying.
	}
	// find the top min(TOP_NUM, Nlocal) nodes.
	partial_sort(nm_list.begin(), nm_list.begin() + min(TOP_NUM, Nlocal),
				 nm_list.end(), comp_node_mark);

	// set up the data for global merge.
	for (int i = 0; i < TOP_NUM; i++) {
		local_top_list_num_marks[TOP_NUM * tid + i] = nm_list[i].num_marks_;
		local_top_list_nn[TOP_NUM * tid + i] = nm_list[i].nn_;
	}
	local_top_list_pointer[tid] = TOP_NUM * tid;

	// barrier for local partial_sort to complete and for time profiling.
	pthread_barrier_wait(&timing_barr);

	if (tid == 0) {
		// merge the local top ranked nodes to a global list.
		for (int i = 0; i < TOP_NUM; i++) {
			int winner_tid = get_winner_tid(local_top_list_num_marks,
											local_top_list_pointer,
											num_threads);
			global_top_list_num_marks[i] =
				local_top_list_num_marks[local_top_list_pointer[winner_tid]];
			global_top_list_nn[i] =
				local_top_list_nn[local_top_list_pointer[winner_tid]];
			// advance the pointer of winner_tid.
			local_top_list_pointer[winner_tid]++;
		}
		gettimeofday(&end_tv, NULL);
		tElapsed = (end_tv.tv_sec + end_tv.tv_usec/1000000.0) -
				   (start_tv.tv_sec + start_tv.tv_usec/1000000.0);
	}

	pthread_exit(0);
}

/*
 *@description - it assumes the node number starts from 0.
 *@param - argv[1] - name of the file storing web graph.
 *@param - argv[2] - name of output file.
 *@param - argv[3] - num_iter, the number of desired steps for RW.
 *@param - argv[4] - num_threads, the number of threads.
 */
int parallel_random_walk(int argc, char *argv[]) {
	if (argc != 5) {
		cerr << "to run this program must supply the following command "
				"line arguments (in order)" << endl;
		cerr << "argv[1]---web graph file." << endl;
		cerr << "argv[2]---output file." << endl;
		cerr << "argv[3]---number of desired steps for random walks." << endl;
		cerr << "argv[4]---number of threads." << endl;
		exit(-1);
	}

	// CSR representation of the (src, dest) matrix.
	int *col_ind = NULL;
	int *row_ptr = NULL;
	// variables for storing number of marks for each node. we use 2D array to
	// store in_walk to avoid false sharing.
	int *out_walk = NULL;		// number of marks for current snapshot.
	const int stride = 32;		// tune this parameter to avoid false sharing.
	int **in_walk = NULL;		// number of marks for next snapshot.
	// locks for in_walk.
	pthread_mutex_t *in_walk_lock = NULL;

	//--------reads in the web graph, change it to CSR format.
	cerr << "reads in the web graph data......" << endl;
	int *src = NULL;
	int *dest = NULL;
	int nnz = 0;
	N = read_webgraph(argv[1], &src, &dest, &nnz) + 1;
	cerr << "There are in total " << N << " nodes in the web graph." << endl;
	col_ind = new int[nnz];
	row_ptr = new int[N+1];
	coord2csr(src, dest, nnz, N, col_ind, row_ptr);
	cerr << "read in the webgraph is done." << endl;

	// get parameters from command line arguments.
	num_threads = atoi(argv[4]);
	num_iter = atoi(argv[3]);

	// prepare the shared data for all threads.
	out_walk = new int[N];
	std::fill(out_walk, out_walk + N, 1);	// initially each node has one mark.
	in_walk = new int*[N];
	for (int i = 0; i < N; i++) {
		in_walk[i] = new int[stride];
		in_walk[i][0] = 0;
	}

	// initialize the mutex in_walk_lock.
	in_walk_lock = new pthread_mutex_t[N];
	for (int i = 0; i < N; i++) {
		pthread_mutex_init(&in_walk_lock[i], NULL);
	}

    if(pthread_barrier_init(&in_walk_barr, NULL, num_threads))
    {
        cerr << "could not create in_walk_barrier in "
        		"parallel_random_walk()" << endl;
        exit(-1);
    }
    if(pthread_barrier_init(&out_walk_barr, NULL, num_threads))
    {
        cerr << "could not create out_walk_barrier in "
        		"parallel_random_walk()" << endl;
        exit(-1);
    }
    if(pthread_barrier_init(&timing_barr, NULL, num_threads))
    {
        cerr << "could not create timing barrier in "
        		"parallel_random_walk()" << endl;
        exit(-1);
    }

	pthread_t *threads_handle = new pthread_t[num_threads];
	// shared data to pass to each thread.
	thread_share_t *share_data = new thread_share_t[num_threads];

	local_top_list_num_marks = new int[TOP_NUM * num_threads];
	local_top_list_nn = new int[TOP_NUM * num_threads];
	local_top_list_pointer = new int[num_threads];

	cerr << "multi-threading random walks are up......" << endl;

	// create threads.
	for(int t = 0; t < num_threads; t++){
		// pack the shared data for the thread.
		share_data[t].tid_ = t;
		share_data[t].col_ind_ = col_ind;
		share_data[t].row_ptr_ = row_ptr;
		share_data[t].out_walk_ = out_walk;
		share_data[t].in_walk_ = in_walk;
		share_data[t].in_walk_lock_ = in_walk_lock;
		int rc = pthread_create(&threads_handle[t], NULL, random_walk_thread,
								(void*)&share_data[t]);
		if (rc) {
			cerr << "fail to create thread " << t
				 << " with the error code " << rc << endl;
			exit(-1);
		}
	}

	// wait for all the threads to complete.
	for(int t = 0; t < num_threads; t++){
		pthread_join(threads_handle[t], NULL);
	}

	cerr << "multi-threading random walks are done." << endl;

	ofstream out(argv[2], ios::trunc);
	out << "Time: " << tElapsed << " seconds when using "
		<< num_threads << " threads for " << num_iter
		<< " steps of random walks." << endl;
	// print out the top TOP_NUM ranked nodes.
	out << "the top " << TOP_NUM << " nodes are:" << endl;
	for (int i = 0; i < TOP_NUM; i++) {
		out << global_top_list_nn[i] << " " << endl;
	}
	cerr << "the solution vector has been saved in file "
		 << argv[4] << endl;
	// print out time.
	cerr << "Time: " << tElapsed << " seconds when using "
		 << num_threads << " threads for " << num_iter
		 << " steps of random walks." << endl;

	// destroy the barriers.
	pthread_barrier_destroy(&in_walk_barr);
	pthread_barrier_destroy(&out_walk_barr);
	pthread_barrier_destroy(&timing_barr);


	// deallocate the storage.
	if (src != NULL) {
		delete [] src;
		src = NULL;
	}
	if (dest != NULL) {
		delete [] dest;
		dest = NULL;
	}
	if (col_ind != NULL) {
		delete [] col_ind;
		col_ind = NULL;
	}
	if (row_ptr != NULL) {
		delete [] row_ptr;
		row_ptr = NULL;
	}
	if (out_walk != NULL) {
		delete [] out_walk;
		out_walk = NULL;
	}
	if (in_walk != NULL) {
		deallocate_pointer_array(in_walk, N);
		in_walk = NULL;
	}
	if (in_walk_lock != NULL) {
		// clean up the mutex in_walk_lock.
		for (int i = 0; i < N; i++) {
			pthread_mutex_destroy(&in_walk_lock[i]);
		}
		delete [] in_walk_lock;
		in_walk_lock = NULL;
	}
	if (threads_handle != NULL) {
		delete [] threads_handle;
		threads_handle = NULL;
	}
	if (share_data != NULL) {
		delete [] share_data;
		share_data = NULL;
	}
	if (local_top_list_num_marks != NULL) {
		delete [] local_top_list_num_marks;
		local_top_list_num_marks = NULL;
	}
	if (local_top_list_nn != NULL) {
		delete [] local_top_list_nn;
		local_top_list_nn = NULL;
	}
	if (local_top_list_pointer != NULL) {
		delete [] local_top_list_pointer;
		local_top_list_pointer = NULL;
	}

	return 0;
}

int main(int argc, char *argv[]) {
	// init the lock for printing purpose.
	pthread_mutex_init(&print_lock, NULL);

	parallel_random_walk(argc, argv);

	// destroy the lock for printing purpose.
	pthread_mutex_destroy(&print_lock);
}
