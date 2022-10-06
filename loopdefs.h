#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

namespace arch{

enum reduceOp { max, min, sum, prod };

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[1], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], thread_data); }

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[2], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], thread_data); }

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[3], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], thread_data); }

template <typename Lambda, typename T>
CUDA_HOSTDEV inline static void lambda_eval(const uint (&idx)[4], T * __restrict__ thread_data, Lambda loop_body) { loop_body(idx[0], idx[1], idx[2], idx[3], thread_data); }

#ifndef USE_CUDA
template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
inline void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) {

  if(nDim == 1){
    uint idx[1];
           
    # pragma omp for
    for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
      lambda_eval(idx, sum, loop_body);
  }
  else if(nDim == 2){
    uint idx[2];
           
    # pragma omp for collapse(2)
    for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
      for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
        lambda_eval(idx, sum, loop_body);
  }
  else if(nDim == 3){
    uint idx[3];
           
    # pragma omp for collapse(3)
    for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
      for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
        for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
          lambda_eval(idx, sum, loop_body);
  }
  else if(nDim == 4){
    uint idx[4];
           
    # pragma omp for collapse(4)
    for (idx[3] = 0; idx[3] < limits[3]; ++idx[3]) 
      for (idx[2] = 0; idx[2] < limits[2]; ++idx[2]) 
        for (idx[1] = 0; idx[1] < limits[1]; ++idx[1]) 
          for (idx[0] = 0; idx[0] < limits[0]; ++idx[0])
            lambda_eval(idx, sum, loop_body);
  }
  else{
    printf("ERROR in %s at line %d: This loop nest dimension is not supported!\n", __FILE__, __LINE__);
    exit(1);
  }
}

#else

#define BLOCKSIZE_R 64

// The macros and functions that can be compiled with a C compiler
#define CUDA_ERR(err) (cuda_error(err, __FILE__, __LINE__))
  inline static void cuda_error(cudaError_t err, const char *file, int line) {
  	if (err != cudaSuccess) {
  		printf("\n\n%s in %s at line %d\n", cudaGetErrorString(err), file, line);
  		// exit(1);
  	}
}

__device__ __forceinline__ static void atomicMax(double *address, double val2) {
  unsigned long long ret = __double_as_longlong(*address);
  while(val2 > __longlong_as_double(ret)) {
    unsigned long long old = ret;
    if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
    break;
  }
}

__device__ __forceinline__ static void atomicMin(double *address, double val2) {
  unsigned long long ret = __double_as_longlong(*address);
  while(val2 < __longlong_as_double(ret)) {
    unsigned long long old = ret;
    if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
    break;
  }
}

template <uint nDim, typename Lambda, typename T>
__device__ __forceinline__ static void loop_eval(const uint idxGlob, const uint * __restrict__ lims, T * __restrict__ thread_data, Lambda loop_body) { 
  uint idx[nDim];
  switch (nDim)
  {
    case 4:
      idx[3] = (idxGlob / (lims[0] * lims[1] * lims[2])) % lims[3];
    case 3:
      idx[2] = (idxGlob / (lims[0] * lims[1])) % lims[2];
    case 2:
      idx[1] = (idxGlob / lims[0]) % lims[1];
    case 1:
      idx[0] = idxGlob % lims[0];  
  }
  lambda_eval(idx, thread_data, loop_body);
}

template <reduceOp op, uint nDim, uint nReduStatic, typename Lambda, typename T>
__global__ static void __launch_bounds__(BLOCKSIZE_R)
reduction_kernel(Lambda loop_body, const T * __restrict__ init_val, T * __restrict__ rslt, const uint * __restrict__ lims, const uint n_total, const uint nReduDynamic, T *thread_data_dynamic)
{
  typedef cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;

  constexpr uint size = nReduStatic ? nReduStatic : 1;

  extern __shared__ typename BlockReduce::TempStorage temp_storage_dynamic[];
  __shared__ typename BlockReduce::TempStorage temp_storage_static[size];
  typename BlockReduce::TempStorage *temp_storage = nReduStatic ? temp_storage_static : temp_storage_dynamic;

  const uint idxGlob = ((blockIdx.x*blockDim.x)+threadIdx.x);

  if (idxGlob < n_total) {

    T thread_data_static[size];
    T *thread_data = nReduStatic ? thread_data_static : &thread_data_dynamic[nReduDynamic * idxGlob];
  
    const uint nReductions = nReduStatic ? nReduStatic : nReduDynamic;  
    for(uint i = 0; i < nReductions; i++)
      thread_data[i] = init_val[i];
  
    // Evaluate the loop body
    loop_eval<nDim>(idxGlob, lims, thread_data, loop_body);
  
    /* Perform reductions */
    for(uint i = 0; i < nReductions; i++){
      /* Compute the block-wide sum for thread0 which stores it */
      if(op == reduceOp::sum){
        T aggregate = BlockReduce(temp_storage[i]).Sum(thread_data[i]);
        if(threadIdx.x == 0) 
          atomicAdd(&rslt[i], aggregate);
      }
      else if(op == reduceOp::max){
        T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Max()); 
        if(threadIdx.x == 0) 
          atomicMax(&rslt[i], aggregate);
      }
      else if(op == reduceOp::min){
        T aggregate = BlockReduce(temp_storage[i]).Reduce(thread_data[i], cub::Min());
        if(threadIdx.x == 0) 
          atomicMin(&rslt[i], aggregate);
      }
      else
        if(threadIdx.x == 0) 
          printf("ERROR at %s:%d: Invalid reduction identifier \"op\".", __FILE__, __LINE__);
    }
  }
}

template <reduceOp op, uint nReduStatic, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T *sum, const uint nReduDynamic) {

  const uint nReductions = nReduStatic ? nReduStatic : nReduDynamic;  

  uint n_total = 1;
  for(uint i = 0; i < nDim; i++)
    n_total *= limits[i];

  const uint blocksize = 64;
  const uint gridsize = (n_total - 1 + blocksize) / blocksize;

  T* d_buf;
  CUDA_ERR(cudaMallocAsync(&d_buf, nReductions*sizeof(T), 0));
  CUDA_ERR(cudaMemcpy(d_buf, sum, nReductions*sizeof(T), cudaMemcpyHostToDevice));
  
  T* d_const_buf;
  CUDA_ERR(cudaMallocAsync(&d_const_buf, nReductions*sizeof(T), 0));
  CUDA_ERR(cudaMemcpy(d_const_buf, d_buf, nReductions*sizeof(T), cudaMemcpyDeviceToDevice));

  uint* d_limits;
  CUDA_ERR(cudaMallocAsync(&d_limits, nDim*sizeof(uint), 0));
  CUDA_ERR(cudaMemcpy(d_limits, limits, nDim*sizeof(uint), cudaMemcpyHostToDevice));

  T* d_thread_data_dynamic;
  if(nReduStatic == 0) {
    constexpr auto cubTempStorageTypeSize = sizeof(typename cub::BlockReduce<T, BLOCKSIZE_R, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1>::TempStorage);
    
    CUDA_ERR(cudaMallocAsync(&d_thread_data_dynamic, nReductions * n_total * sizeof(T), 0));

    reduction_kernel<op, nDim, 0><<<gridsize, blocksize, nReductions * cubTempStorageTypeSize>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, nReductions, d_thread_data_dynamic);

    CUDA_ERR(cudaStreamSynchronize(0));
    CUDA_ERR(cudaFreeAsync(d_thread_data_dynamic, 0));
  }
  else{
    reduction_kernel<op, nDim, nReduStatic><<<gridsize, blocksize>>>(loop_body, d_const_buf, d_buf, d_limits, n_total, nReductions, d_thread_data_dynamic);

    CUDA_ERR(cudaStreamSynchronize(0));
  }
  
  CUDA_ERR(cudaMemcpy(sum, d_buf, nReductions*sizeof(T), cudaMemcpyDeviceToHost));
  CUDA_ERR(cudaFreeAsync(d_buf, 0));
  CUDA_ERR(cudaFreeAsync(d_const_buf, 0));
  CUDA_ERR(cudaFreeAsync(d_limits, 0));
}

#endif

template <reduceOp op, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T &sum) {
  constexpr uint nReductions = 1;
  arch::parallel_reduce<op, nReductions>(limits, loop_body, &sum, nReductions);
}

template <reduceOp op, uint nDim, uint nReductions, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, T (&sum)[nReductions]) { 
  arch::parallel_reduce<op, nReductions>(limits, loop_body, &sum[0], nReductions);
}

template <reduceOp op, uint nDim, typename Lambda, typename T>
__forceinline__ static void parallel_reduce(const uint (&limits)[nDim], Lambda loop_body, std::vector<T> &sum) {
  arch::parallel_reduce<op, 0>(limits, loop_body, sum.data(), sum.size());
}

}