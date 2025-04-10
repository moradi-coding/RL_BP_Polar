

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


int main() {

	int **H,m=2049,n=4098,i,j,cnt; //,gama=4,p,flg,cnt;
	//H is the parity-check matrix, m in the number of rows, n is the number of columns

	H= new int*[m]; 
	for(i=0;i<m;i++) 
		H[i]=new int[n];


	ifstream inf;
	inf.open("mat_3_6.txt");

	for(j=0;j<m;j++) 
		for(i=0;i<n;i++) 
			inf >> H[j][i];

	cout<<'\n'<<"VN degrees: ";
	for(i=0;i<n;i++) {
		cnt=0; 
		for(j=0;j<m;j++) if(H[j][i]) cnt++;
		cout<<cnt<<" ";
	}
	
	cout<<'\n'<<'\n'<<"CN degrees: ";
	for(j=0;j<m;j++) {
		cnt=0; 
		for(i=0;i<n;i++) if(H[j][i]) cnt++;
		cout<<cnt<<" "; 
	}


	//checks if the row and column weights are gama and p, respectively
	/*for(j=0;j<n;j++) {
		cnt=0; 
		for(i=0;i<m;i++) if(H[i][j]) cnt++;
		if(cnt!=gama) {flg=1; cout<<'\n'<<"col no.: "<<j; break;}
		else flg=0;
	}
	if(!flg) cout<<'\n'<<"valid col. wt."; 
	else cout<<'\n'<<"invalid col. wt.";

	for(i=0;i<m;i++) {
		cnt=0; 
		for(j=0;j<n;j++) if(H[i][j]) cnt++;
		if(cnt!=p) {flg=1; cout<<'\n'<<"row no.: "<<i; break;}
		else flg=0;
	}
	if(!flg) cout<<'\n'<<"valid row. wt."; 
	else cout<<'\n'<<"invalid row. wt.";*/
	

}
