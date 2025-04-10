//construct minimal span generator matrix of hamming code

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
using namespace std;

string func(int n) {
    stringstream result;
    result << n;
    return result.str();
}

int main() {	
	srand(time(0));	
	double t2,t1=time(NULL);	
	int i,i5,i2,i3,i4,j,j2,k,l,*Lx,*Rx,m,n,cnt,flg,**G,**span,**A;
	
	ifstream inf; 
	ofstream outf,outf2; 
	string filename; 
	
	m=4; n=7; inf.open("matrices/7_4_Hamm_G.txt"); filename="7_4_Hamm_G_2.txt"; 
	//m=11; n=15; inf.open("matrices/15_11_Hamm_G.txt"); filename="15_11_Hamm_G_2.txt"; 
	
	G=new int*[m]; 
	for(i=0;i<m;i++) 
		G[i]=new int[n]; //m rows n cols
	
	span=new int*[m]; 
	for(i=0;i<m;i++) 
		span[i]=new int[m]; //m rows n cols
	for(j=0;j<m;j++) 
		for(i=0;i<m;i++) 
			span[j][i]=-1;
			
	A=new int*[n]; 
	for(i=0;i<n;i++) 
		A[i]=new int[m]; //m rows n cols
	for(j=0;j<n;j++) 
		for(i=0;i<m;i++) 
			A[j][i]=-1;
		
	for(j=0;j<m;j++) for(i=0;i<n;i++) inf>>G[j][i]; 
	
	Lx=new int[m];
	Rx=new int[m]; 
	
	
	cnt=0;
	flg=1;
	while(flg && cnt<15) {
		for(j=0;j<m;j++)
			for(i=0;i<n;i++)
				if(G[j][i]) {
					Lx[j]=i+1;
					break;
				}
				
		for(j=0;j<m;j++)
			for(i=n;i>=0;i--)
				if(G[j][i]) {
					Rx[j]=i+1;
					break;
				}
				
		cout<<"Lx: "; for(j=0;j<m;j++) cout<<Lx[j]<<" "; cout<<'\n';
		cout<<"Rx: "; for(j=0;j<m;j++) cout<<Rx[j]<<" "; cout<<'\n';
		
		for(j=0;j<m;j++) {
			for(i=0;i<n;i++) 
				cout<<G[j][i]<<" "; 
			cout<<'\n';
		}
		
		for(i=0;i<m;i++)
			for(j=i+1;j<m;j++)
				if(Rx[i]==Rx[j]) {
					if(Lx[i]>=Lx[j]) {
						for(i2=0;i2<n;i2++)
							G[j][i2]=(G[j][i2]+G[i][i2])%2;
					}
					else {
						for(i2=0;i2<n;i2++)
							G[i][i2]=(G[i][i2]+G[j][i2])%2;
					}
				}
				else if(Lx[i]==Lx[j]) {
					if(Rx[i]<=Rx[j]) {
						for(i2=0;i2<n;i2++)
							G[j][i2]=(G[j][i2]+G[i][i2])%2;
					}
					else {
						for(i2=0;i2<n;i2++)
							G[i][i2]=(G[i][i2]+G[j][i2])%2;
					}
				}
		
		for(i=0;i<m;i++) {
			for(j=i+1;j<m;j++)
				if(Rx[i]==Rx[j] || Lx[i]==Lx[j]) {	
					flg=1;
					break;	
				}
				else
					flg=0;
			if(flg)
				break;
		}
				
		cnt++;
	}	

	//output the minimal span of G
	/*cout<<"G: "<<'\n';
	outf.open(filename.c_str()); 
	for(j=0;j<m;j++) {
		for(i=0;i<n;i++) 
			outf<<G[j][i]<<" "; 
		outf<<'\n';
	}
	outf<<std::endl; 
	outf.close();*/
	
	//finding span of each row of minimal span G
	for(i=0;i<m;i++) {
		cnt=0;
		for(j=Lx[i];j<=Rx[i];j++) {
			span[i][cnt]=j;
			cnt++;
		}
	}
	
	cout<<"span: "<<'\n';
	for(j=0;j<m;j++) {
		for(i=0;i<m;i++) 
			cout<<span[j][i]<<" "; 
		cout<<'\n';
	}
	
	for(i2=1;i2<=n;i2++) {
		cnt=0;
		for(i=0;i<m;i++) 
			for(j=0;j<m;j++) 
				if(span[i][j]==i2) {
					A[i2-1][cnt]=i+1;
					cnt++;
					break;
				}
	
	}
			
	cout<<"A: "<<'\n';
	for(j=0;j<n;j++) {
		for(i=0;i<m;i++) 
			cout<<A[j][i]<<" "; 
		cout<<'\n';
	}
	
	//finding lambda(u) for different cardinality of the A vector and a fixed multiplier
	int a=8, b=3;
	//int a=16, b=4;
	//int a=32, b=5;
	int lambda[a];
	for(i=0;i<a;i++)
    	lambda[i]=0;
    	
	int mult[b]={0,0,1};
    	
	int data[a][b]={
        {0,0,0},
        {0,0,1},
        {0,1,0},
        {0,1,1},
        {1,0,0},
        {1,0,1},
        {1,1,0},
        {1,1,1},
    };
    /*int data[a][b]={
        {0,0,0,0},
        {0,0,0,1},
        {0,0,1,0},
        {0,0,1,1},
        {0,1,0,0},
        {0,1,0,1},
        {0,1,1,0},
        {0,1,1,1},
        {1,0,0,0},
        {1,0,0,1},
        {1,0,1,0},
        {1,0,1,1},
        {1,1,0,0},
        {1,1,0,1},
        {1,1,1,0},
        {1,1,1,1},
    };*/
    /*int data[a][b]={
        {0,0,0,0,0},
        {0,0,0,0,1},
        {0,0,0,1,0},
        {0,0,0,1,1},
        {0,0,1,0,0},
        {0,0,1,0,1},
        {0,0,1,1,0},
        {0,0,1,1,1},
        {0,1,0,0,0},
        {0,1,0,0,1},
        {0,1,0,1,0},
        {0,1,0,1,1},
        {0,1,1,0,0},
        {0,1,1,0,1},
        {0,1,1,1,0},
        {0,1,1,1,1},
        {1,0,0,0,0},
        {1,0,0,0,1},
        {1,0,0,1,0},
        {1,0,0,1,1},
        {1,0,1,0,0},
        {1,0,1,0,1},
        {1,0,1,1,0},
        {1,0,1,1,1},
        {1,1,0,0,0},
        {1,1,0,0,1},
        {1,1,0,1,0},
        {1,1,0,1,1},
        {1,1,1,0,0},
        {1,1,1,0,1},
        {1,1,1,1,0},
        {1,1,1,1,1},
    };
	/*cout<<'\n'<<"data: ";
    for(i=0;i<a;i++) {
    	for(j=0;j<b;j++)
    		cout<<data[i][j]<<" ";
    	cout<<'\n';	
    }*/
    
    for(i=0;i<a;i++) {
    	for(j=0;j<b;j++)
    		lambda[i]+=mult[j]*data[i][j];
    	lambda[i]=lambda[i]%2;
    }
    
    cout<<'\n'<<"mult: ";
    for(i=0;i<b;i++)
    	cout<<mult[i]<<" ";
    	
    cout<<'\n'<<"lambda: ";
    for(i=0;i<a;i++) {
    	cout<<lambda[i]<<" ";
		if(i>0 && !((i+1)%4))
			cout<<", ";
	}
	//t2= time(NULL); cout<<'\n'<<"Executed in "<<t2-t1<<"s"<<'\n';
	cout<<'\n';
}
