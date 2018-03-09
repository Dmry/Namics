#include "tools.h"  
#include "namics.h"
	
 
#ifdef CUDA
//cublasHandle_t handle;
//cublasStatus_t stat=cublasCreate(&handle);

/*
__global__ void distributeg1(double *G1, double *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY)   {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	if (k < Mz && j < My && i < Mx){
		int pM=jx+jy+1+k+jx*i+jy*j; 
		int pos_r=JX+JY+1;
		int ii;
		int MXm1 = MX-1;
		int MYm1 = MY-1;
		int MZm1 = MZ-1;
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1) ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX;
			if (By[p]+j > MYm1) ii+=(By[p]+j-MY)*JY; else ii+=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1) ii+=(Bz[p]+k-MZ); else ii+=(Bz[p]+k);
			g1[pM]=G1[pos_r+ii];
			pM+=M;
		}
	}
} 

__global__ void collectphi(double* phi, double* GN,double* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY)   {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	int pos_r=(JX+JY+1);//This needs af fix**
	int pM=jx+jy+1+jx*i+jy*j+k;
	int ii,jj,kk;//This needs af fix**
	int MXm1 = MX-1;
	int MYm1 = MY-1;
	int MZm1 = MZ-1;
	if (k < Mz && j < My && i < Mx){
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1)  ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX;
			if (By[p]+j > MYm1)  jj=(By[p]+j-MY)*JY; else jj=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1)  kk=(Bz[p]+k-MZ); else kk=(Bz[p]+k);

			atomicAdd(&phi[pos_r+ii+jj+kk], rho[pM]/GN[p]); //This needs af fix*****************************
			pM+=M;
		}
	}
}
*/

__global__ void dot(double *X, double *Y, double *Z, int M)   {
__shared__ double cache[256];
	int idx=blockDim.x*blockIdx.x+threadIdx.x; 
	int CI=threadIdx.x;
	cache[CI]=0.0;
	while (idx<M) {
		cache[CI] += X[idx]*Y[idx];
		idx += blockDim.x*gridDim.x;
	}
	__syncthreads();
	int i=256/2;
	while (i!=0) {
		if (CI < i) cache[CI] +=cache[CI+i];
		__syncthreads();
		i/=2;
	}
	if (CI==0) Z[blockIdx.x]=cache[0]; 
}


__global__ void sum(double *X, double *Z, int M)   {
__shared__ double cache[256];
	int idx=blockDim.x*blockIdx.x+threadIdx.x; 
	int CI=threadIdx.x;
	cache[CI]=0.0;
	while (idx<M) {
		cache[CI] +=X[idx];
		idx+=blockDim.x*gridDim.x;
	}
	__syncthreads();
	int i=256/2;
	while (i!=0) {
		if (CI < i) cache[CI] +=cache[CI+i];
		__syncthreads();
		i/=2;
	}
	if (CI==0) Z[blockIdx.x]=cache[0]; 
}


__global__ void times(double *P, double *A, double *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(double *P, double *A, double *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(double *P, double C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(double *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void cp (double *P, double *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void yisaminb(double *Y, double *A,double *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yisaplusc(double *Y, double *A, double C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+C;
}
__global__ void yisaplusb(double *Y, double *A,double *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+B[idx];
}
__global__ void yplusisctimesx(double *Y, double *X, double C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void updatealpha(double *Y, double *X, double C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*(X[idx]-1.0);
}
__global__ void picard(double *Y, double *X, double C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = C*Y[idx]+(1-C)*X[idx];
}
__global__ void add(double *P, double *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void dubble(double *P, double *A, double norm, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {if (A[idx]>0) P[idx]*=norm/A[idx]; else {P[idx]=0;}} 
}
__global__ void boltzmann(double *P, double *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void invert(double *SKAM, double *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (MASK[idx]==0.0) SKAM[idx]=1.0; else SKAM[idx]=0.0;
}
__global__ void putalpha(double *g,double *phitot,double *phi_side,double chi,double phibulk,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (phitot[idx]>0) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}
__global__ void div(double *P,double *A,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx] /=A[idx];
}
__global__ void oneminusphitot(double *g, double *phitot, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {if (phitot[idx]>0) g[idx]= 1/phitot[idx]-1; else g[idx]=0;}
}
__global__ void addg(double *g, double *phitot, double *alpha, int M)    {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		if (phitot[idx]>0) g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1; else g[idx]=0;
	}
}

//__global__ void computegn(double *GN, double* G, int M, int n_box)   {
//	int p= blockIdx.x*blockDim.x+threadIdx.x;
//	if (p<n_box){
//		cublasDasum(handle,M,G+p*M,1,GN+p); //in case of float we need the cublasSasum etcetera.
//	}
//}

void TransferDataToHost(double *H, double *D, int M)    {
	cudaMemcpy(H, D, sizeof(double)*M,cudaMemcpyDeviceToHost);
}
void TransferDataToDevice(double *H, double *D, int M)    {
	cudaMemcpy(D, H, sizeof(double)*M,cudaMemcpyHostToDevice);
}
void TransferIntDataToHost(int *H, int *D, int M)   {
	cudaMemcpy(H, D, sizeof(int)*M,cudaMemcpyDeviceToHost);
	if (!CUDA_SUCCESS) cout <<"no success" << endl; else cout <<"success" << endl;
}
void TransferIntDataToDevice(int *H, int *D, int M)     {
	cudaMemcpy(D, H, sizeof(int)*M,cudaMemcpyHostToDevice);
	if (!CUDA_SUCCESS) cout <<"no success" << endl; else cout <<"success" << endl; 
}
#endif

#ifdef CUDA
__global__ void bx(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}


#else
void bx(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
	}
}
void b_x(double *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(double *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(double *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}
#endif

#ifdef CUDA
bool GPU_present()    {
	cudaDeviceReset(); 
	cudaGetLastError(); 
	int deviceCount =0; cuDeviceGetCount(&deviceCount);
   	if (deviceCount ==0) printf("There is no device supporting Cuda.\n"); else cudaDeviceReset();
	//if (deviceCount>0) {
		//stat = cublasCreate(&handle); 
		//if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS handle creation failed \n");}
		//cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	//}
	return deviceCount > 0;
}
#else

bool GPU_present()    {
	return false;
}
#endif


#ifdef CUDA
double *AllOnDev(int N)    {
	double *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(double)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}
int *AllIntOnDev(int N)    {
	int *X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(int)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}
#else
double *AllOnDev(int N)    {
	double *X=NULL;
	printf("Turn on CUDA\n");
	return X;
}
int *AllIntOnDev(int N)    {
	int *X=NULL;
	printf("Turn on CUDA\n");
	return X;
}
#endif

#ifdef CUDA
double Dot(double *x,double *y, int M)   {
	double result=0;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	double *D_XX=AllOnDev(n_blocks);
	dot<<<n_blocks,block_size>>>(x,y,D_XX,M);
	double *H_XX=new double[n_blocks];
	TransferDataToDevice(H_XX, D_XX, n_blocks);
	for (int i=0; i<n_blocks; i++) result+=H_XX[i];
	delete[] H_XX;
	cudaFree(D_XX); 
	return result;
}
#else
double Dot(double *x,double *y, int M)   {
	double result=0.0;
 	for (int i=0; i<M; i++) result +=x[i]*y[i];
	return result;
}
#endif

#ifdef CUDA
double Sum(double *x, int M)   {
	double result=0;
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	double *D_XX=AllOnDev(n_blocks);
	sum<<<n_blocks,block_size>>>(x,D_XX,M);
	double *H_XX=new double[n_blocks];
	TransferDataToDevice(H_XX, D_XX, n_blocks);
	for (int i=0; i<n_blocks; i++) result+=H_XX[i];
	delete[] H_XX;
	cudaFree(D_XX); 
	return result;
}
#else
double Sum(double *x,int M)   {
	double result=0.0;
 	for (int i=0; i<M; i++) result +=x[i];
	return result;
}
#endif

#ifdef CUDA
void AddTimes(double *P, double *A, double *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void AddTimes(double *P, double *A, double *B, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}
#endif

#ifdef CUDA
void Times(double *P, double *A, double *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}
#else
void Times(double *P, double *A, double *B, int M)   {
	for (int i=0; i<M; i++) P[i]=A[i]*B[i];
}
#endif

#ifdef CUDA
void Norm(double *P, double C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
}
#else
void Norm(double *P, double C, int M)   {
	for (int i=0; i<M; i++) P[i] *= C;
}
#endif

#ifdef CUDA
void Zero(double* P, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size>>>(P,M);
}
#else
void Zero(double* P, int M)   {
	for (int i=0; i<M; i++) P[i] =0;
}
#endif

#ifdef CUDA
void Cp(double *P,double *A, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Cp(double *P,double *A, int M)   {
	for (int i=0; i<M; i++) P[i] = A[i];
}
#endif

#ifdef CUDA
void YisAminB(double *Y, double *A, double *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
}
#else
  void YisAminB(double *Y, double *A, double *B, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]-B[i];
}
#endif
#ifdef CUDA
void YisAplusC(double *Y, double *A, double C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusc<<<n_blocks,block_size>>>(Y,A,C,M);
}
#else
void YisAplusC(double *Y, double *A, double C, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+C;
}
#endif
#ifdef CUDA
void YisAplusB(double *Y, double *A, double *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusb<<<n_blocks,block_size>>>(Y,A,B,M);
}
#else
void YisAplusB(double *Y, double *A, double *B, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+B[i];
}
#endif

#ifdef CUDA
void YplusisCtimesX(double *Y, double *X, double C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
}
#else
void YplusisCtimesX(double *Y, double *X, double C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}
#endif

#ifdef CUDA
void UpdateAlpha(double *Y, double *X, double C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	updatealpha<<<n_blocks,block_size>>>(Y,X,C,M);
}
#else
void UpdateAlpha(double *Y, double *X, double C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*(X[i]-1.0);
}
#endif
#ifdef CUDA
  void Picard(double *Y, double *X, double C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	picard<<<n_blocks,block_size>>>(Y,X,C,M);
}
#else
void Picard(double *Y, double *X, double C, int M)    {
	for (int i=0; i<M; i++) Y[i] = C*Y[i]+(1.0-C)*X[i];
}
#endif

#ifdef CUDA
void Add(double *P, double *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Add(double *P, double *A, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i];
}
#endif

#ifdef CUDA
void Dubble(double *P, double *A, double norm,int M)   {
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
}
#else
void Dubble(double *P, double *A, double norm,int M)   {
	for (int i=0; i<M; i++) {if (A[i]>0) P[i]*=norm/A[i]; else P[i]=0;} 
}
#endif

#ifdef CUDA
void Boltzmann(double *P, double *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Boltzmann(double *P, double *A, int M)   {
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}
#endif
#ifdef CUDA
void Invert(double *KSAM, double *MASK, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	invert<<<n_blocks,block_size>>>(KSAM,MASK,M);
}
#else
void Invert(double *KSAM, double *MASK, int M)   {
	for (int i=0; i<M; i++) if (MASK[i]==0.0) KSAM[i]=1.0; else KSAM[i]=0.0;
}
#endif

#ifdef CUDA
void PutAlpha(double *g, double *phitot, double *phi_side, double chi, double phibulk, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
}
#else
void PutAlpha(double *g, double *phitot, double *phi_side, double chi, double phibulk, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}
#endif

#ifdef CUDA
void Div(double *P, double *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	div<<<n_blocks,block_size>>>(P,A,M);
}
#else
void Div(double *P, double *A, int M)   {
	for (int i=0; i<M; i++) if (A[i]>0) P[i]/=A[i]; else P[i]=0;
}
#endif

#ifdef CUDA
void AddG(double *g, double *phitot, double *alpha, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}
#else
void AddG(double *g, double *phitot, double *alpha, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0)  g[i]= g[i] -alpha[i] +1/phitot[i]-1.0; else g[i]=0;
}
#endif

#ifdef CUDA
void OneMinusPhitot(double *g, double *phitot, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	oneminusphitot<<<n_blocks,block_size>>>(g,phitot,M);
}
#else
void OneMinusPhitot(double *g, double *phitot, int M)   {
	for (int i=0; i<M; i++) {if(phitot[i]>0) g[i]= 1/phitot[i]-1; else g[i]=0;}
}
#endif

//#ifdef CUDA
//void ComputeGN(double *GN, double *G, int M, int n_box)   {
//	int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
//	computegn<<<n_blocks,block_size>>>(GN,G,M,n_box);
//}
//#else
//void ComputeGN(double* GN, double* G, int M, int n_box)    {
//	for (int p=0; p<n_box; p++) GN[p]=Sum(G+p*M,M); 
//}
//#endif

void H_Zero(double* H_P, int M)   {//this precedure should act on a host P.
	for (int i=0; i<M; i++) H_P[i] = 0;
}

double H_Sum(double* H, int M){
	double Sum=0;
	for (int i=0; i<M; i++) Sum+=H[i];
	return Sum;
}

#ifdef CUDA
void SetBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#else
void SetBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
	bx(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif


#ifdef CUDA
void RemoveBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#else
void RemoveBoundaries(double *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)    {
	b_x(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}
#endif






