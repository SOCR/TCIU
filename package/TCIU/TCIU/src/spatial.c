
#include <math.h>
#include <R.h>


/* functions for calculating spatial clusters and their properties by thresholding a statistic map/image. A neighbouthood matrix nmat is used to determine connectivity. */


void twovoxtyp(Sint*, Sint*, Sint*, Sint*, Sint*); 


void cluster_mass(float *array, Sint *array_dim, Sint *nmat, Sint *nmat_dim,
float *thresh_value, Sint *ans, float *ans1){

  /* thresholds a 3D array and calculates the number of clusters, the maximum entry in each cluster and the mass above the threshold of each cluster
     
     array - the 3D array
     array_dim - the dimensions of the array
     nmat - the neighbourhood system i.e. a matrix with 3 columns, each row represents an allowable difference between two neighbours
     nmat_dim - the dimensions of nmat
     thresh_value - the value to threshold the array at
     ans - the number of clusters found
     ans1 - an (ans x 6) matrix. First 3 columns give the coordinates of the maximum element in each cluster, column 4 gives the value, column 5 gives the number of elements in the cluster and column 6 gives the cluster mass.*/
 
 
  int i,j,k,num=0,n,row,comp,clust_num;
  Sint *vox_mat,*vox_rel,x,y,z,t;  
  float *vox_mat_vals;
  
  

 x= *(array_dim);
 y= *(array_dim+1);
 z= *(array_dim+2);

 vox_mat=Calloc(3,Sint);
 vox_mat_vals=Calloc(1,float);


 for(i=0;i<x;i++){
    for(j=0;j<y;j++){
      for(k=0;k<z;k++){
	
	if(*(array+i*y*z+j*z+k)>*thresh_value){

	  vox_mat=Realloc(vox_mat,3*(num+1),Sint);
	  vox_mat_vals=Realloc(vox_mat_vals,(num+1),float);
	  
	  *(vox_mat+num*3)=(Sint) i+1;
	  *(vox_mat+num*3+1)=(Sint) j+1;
	  *(vox_mat+num*3+2)=(Sint) k+1;
	  *(vox_mat_vals+num)=*(array+i*y*z+j*z+k);
	  num+=1;
	 
	}
      }
    }
 }
 
    if(num>0){
	vox_rel=Calloc(num*num,Sint);
	
	for(i=0;i<(num-1);i++){
	  for(j=i+1;j<num;j++){
	    
	    twovoxtyp(vox_mat+i*3,vox_mat+j*3,nmat,nmat_dim,vox_rel+i*num+j);
	    *(vox_rel+j*num+i)=*(vox_rel+i*num+j);

	  }
	}
	
	for(i=0;i<num;i++) *(vox_rel+i*num+i)=1;

/* 	for(i=0;i<num;i++){ */
/* 	  for(j=0;j<num;j++){ */
/* 	    *(ans2+i*num+j)=*(vox_rel+i*num+j); */
/* 	  }} */

	n=num;
	row=num;
	
	while(row>1)
	  {

	  comp=row-1;
	  
	  while(comp>0)
	    {
	    t=0;
	    for(i=0;i<num;i++) t+=*(vox_rel+(row-1)*num+i)*(*(vox_rel+(comp-1)*num+i)); 
	    if(t>0)
	      {
		for(i=0;i<num;i++) *(vox_rel+(comp-1)*num+i)+=*(vox_rel+(row-1)*num+i); 
		for(i=0;i<num*(n-row);i++) *(vox_rel+(row-1)*num+i)=*(vox_rel+(row)*num+i);
	  
		n-=1;
		vox_rel=Realloc(vox_rel,num*n,Sint);
		comp=0;
	      }
	    else { comp-=1;}

	    }
	  row-=1;}  
    
	*ans=(Sint) n;
	
	for(i=0;i<n;i++){
	  for(j=0;j<num;j++){
	    if(*(vox_rel+i*num+j)>0) *(vox_rel+i*num+j)=(Sint) i+1;
	  }
	}



	
	for(j=0;j<num;j++){
	  t=0;
	  for(i=0;i<n;i++){
	    t+=*(vox_rel+i*num+j);
	  }
	  *(vox_rel+j)=t;
	}

	for(i=0;i<x;i++){
	  for(j=0;j<y;j++){
	    for(k=0;k<z;k++){
	      *(array+i*y*z+j*z+k)=0.0;
	    }
	  }
	}

	for(j=0;j<num;j++){
	  *(array +(*(vox_mat+j*3)-1)*y*z+(*(vox_mat+j*3+1)-1)*z+
	    (*(vox_mat+j*3+2)-1))=(float) *(vox_rel+j) ;
	}

        for(i=0;i<num;i++){
	  clust_num=(int) *(vox_rel+i); 
	  *(ans1+6*(clust_num-1)+4) += 1;
	  *(ans1+6*(clust_num-1)+5) += *(vox_mat_vals+i)-*thresh_value;
	  
	  if(*(vox_mat_vals+i)>*(ans1+6*(clust_num-1)+3)){
	    *(ans1+6*(clust_num-1))=(float) *(vox_mat+3*i);
	    *(ans1+6*(clust_num-1)+1)=(float) *(vox_mat+3*i+1);
	    *(ans1+6*(clust_num-1)+2)=(float) *(vox_mat+3*i+2);
	    *(ans1+6*(clust_num-1)+3)=*(vox_mat_vals+i); 

	  }
	}

	

       
	Free(vox_rel);
    } else { *ans =0;}

	Free(vox_mat);
	Free(vox_mat_vals);
} 



void twovoxtyp(Sint *voxel1, Sint *voxel2, Sint *nmat, Sint *nmat_dim, Sint *ans){ 
/* tests whether two voxels are neighbours using neighbourhood matrix nmat */


  int i;
  Sint *dif; 
  dif=Calloc(3,Sint); 
 
  *ans=0;
 
  for(i=0;i<3;i++) dif[i]=(*(voxel1+i))-(*(voxel2+i)); 
 
  
  i=0; 
 
  do 
  if(dif[0]==*(nmat+i*3) && dif[1]==*(nmat+i*3+1) && dif[2]==*(nmat+i*3+2)) { 
    *ans=1; 
    i=*nmat_dim;} 
 
  else i+=1; 
  while (i<*nmat_dim); 
  Free(dif); 
} 

void covariance_est(double *array, int *array_dim, int *mask, int *nmat, int *nmat_dim, double *ans){

 
  int i, j, k, i1, j1, k1, x, y, z, n;  
  double mean = 0.0, num_mask = 0.0, cov = 0.0, cov_num = 0.0;
  
  x = *(array_dim);
  y = *(array_dim + 1);
  z = *(array_dim + 2);

  for(i = 0; i < x; i++){
	  for(j = 0; j < y; j++){
		  for(k = 0; k < z; k++){
			  
			  if(*(mask + i * y * z + j * z + k)) {
				  mean += *(array + i * y * z + j * z + k);
				  num_mask += 1.0;
			  }
		  }
	  }
  }
  mean /= num_mask;


  for(i = 0; i < x; i++){
	  for(j = 0; j < y; j++){
		  for(k = 0; k < z; k++){
			  
			  if(*(mask + i * y * z + j * z + k)) {
				  
				  for(n = 0; n < *nmat_dim; n++) {
					  i1 = i + *(nmat + n * 3 + 0);
					  j1 = j + *(nmat + n * 3 + 1);
					  k1 = k + *(nmat + n * 3 + 2);
					  if(i1 >= 0 && i1 < x && j1 >= 0 && j1 < y && k1 >= 0 && k1 < z) {
						  if(*(mask + i1 * y * z + j1 * z + k1)) {
							  cov += (*(array + i * y * z + j * z + k) - mean) * (*(array + i1 * y * z + j1 * z + k1) - mean);
							  cov_num += 1.0;
						  }
					  }
				  }
			  }
		  }
	  }
  }
  cov /= cov_num;

  *ans = cov;

}




