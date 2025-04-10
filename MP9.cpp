
// sequential Q-learning based GLDPC decoding using Hamming cluster
// one or all clusters decoded per iteration
#include "clus_decode.cpp"
#include "clus_decode_bcjr.cpp"

int MP9(int cnt, int l) { 
		
	long a,i,i2,cw,j,k,cnt2,flg,flg4,s,row;
	double tmp,tmp2;
	
	ifstream inf;
	ofstream outf1, outf2, outf3;

	for(cw=0;cw<CW;cw++) {
		//flg4=1;
		//for(k=0;k<cnt;k++) 
			//if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, cw has been recovered already
		
		//if(flg4) {			
			//generate random scheduling order
			cnt2=1; 
			//if(!l) {
				for(i=0;i<num_cls;i++) j_sv[i]=-1; //refresh
				while(cnt2<=num_cls) {
					flg=0; 
					j=rand()%num_cls;
					for(i=0;i<cnt2;i++) 
						if(j!=j_sv[i]) flg=1; 
						else {flg=0; break;} 
							
					if(flg) {j_sv[cnt2-1]=j; cnt2++;} 
				}
				//for(i=0;i<num_cls;i++) j_sv[i]=i;
			//}
			//cout<<'\n'<<"j_sv: "; for(i=0;i<num_cls;i++) cout<<j_sv[i]<<" "; cout<<endl;
			
			if(meth==4 || meth==6) { //do learned scheduling
				if(!drl) {
					for(a=0;a<num_cls;a++) {
						//cout<<'\n'<<"HERE1 "<<endl;
						//cout<<'\n'<<"a: "<<a<<endl;
						for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a'
							i=vns_cluster[a][j]; //i is a VN of cluster 'a'
							if(i!=-1) {	
								//if(M==2) {
									if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
									else x_hat[cw*n+i]=1;		
								//}
								/*else if(M>2) {
									if(pLR[cw*n+i]<=(rep[0]+rep[1])/2) x_hat[cw*n+i]=0;
									else if(pLR[cw*n+i]>=(rep[M-2]+rep[M-1])/2) x_hat[cw*n+i]=M-1;
									else 
										for(i2=1;i2<M-1;i2++) 
											if(pLR[cw*n+i]>rep[i2] && pLR[cw*n+i]<=rep[i2+1]) x_hat[cw*n+i]=i2;							
								}*/
							}								
						}				
						//cout<<'\n'<<"x_hat: "<<endl; for(i=0;i<n;i++) cout<<x_hat[cw*n+i]<<" "<<endl;
							
						s=0; 
						for(j=0;j<num_vns_cls;j++) {
							i=vns_cluster[a][j]; 
							if(i!=-1) {	
								//s+=x_hat[cw*n+i]*pow(M,num_vns_cls-j-1);
								s+=x_hat[cw*n+i]*pow(2,num_vns_cls-j-1); //binary to decimal conversion
								//inp_vec[j]=x_hat[cw*n+i];
							}
						//else inp_vec[j]=0;
						}
					
					//finding action that gives max. Q value in state s

						row=s; 
						if(meth==4) {
							Q_temp[a]=Q[row*num_cls+a]; //RELDEC SNR mixture
								//if(Q_temp[a])
									//cout<<Q_temp[a]<<endl;
							}
							else if(meth==6) 
								Q_temp[a]=Q2[row*num_cls+a];
							
							/*else if(meth==5) Q_temp[a]=Q2[row*num_cls+a]; //RELDEC SNR specific
							else if(meth==6) Q_temp[a]=Q3[row*num_cls+a]; //M-RELDEC
							else if(meth==7) Q_temp[a]=Q4[row*num_cls+a]; //AM-RELDEC-7
							else if(meth==8) Q_temp[a]=Q5[row*num_cls+a];*/ //AM-RELDEC-75
							//cout<<'\n'<<"Q,row,a: "<<Q[row*num_cls+a]<<" row: "<<row<<" a: "<<a<<endl;

					/*else if(DeepRL && meth!=3) {
						//cout<<"here"<<endl;
						inp.assign(x_hat, x_hat+n);
						const fdeep::tensor inp_tnsr(fdeep::tensor_shape(static_cast<std::size_t>(n)),inp); //converting vector to tensor
						const auto result = model0.predict({inp_tnsr}); const auto y=result[0]; const std::vector<double> vec = y.to_vector(); Q_temp[a]=vec[a];											
					}*/
					//cout<<'\n'<<"s,row: "<<s<<", "<<row;
					}
					//cout<<'\n'<<'\n'<<"l: "<<l<<" ";
					//cout<<'\n'<<"Qtemp: "<<endl; for(i=0;i<num_cls;i++) cout<<Q_temp[i]<<" ";
				
					//creating cluster decoding order
					for(i=0;i<num_cls;i++) 
						indx[i]=i; //refresh
					for(i=0;i<num_cls;i++) {
						for(i2=i+1;i2<num_cls;i2++) {
							if(Q_temp[i]<Q_temp[i2]) {   
								tmp=Q_temp[i];
								Q_temp[i]=Q_temp[i2];
								Q_temp[i2]=tmp;
					
								tmp2=indx[i];
								indx[i]=indx[i2];
								indx[i2]=tmp2;
							}			
						}
					}
				}
				
				else if(drl && l==1) {
					outf3.open("pLR.txt"); 
					for(i=0;i<n;i++) 
						outf3<<pLR[0*n+i]<<" "; 
					outf3.close(); 
		
					system("python predict.py");
					inf.open("out.txt"); 
					for(i2=0;i2<m;i2++) 
						inf >> indx[i2];
						if(indx[i2]>=m)
							indx[i2]=m-1;
				}
							
				//if(!l) {
					//cout<<'\n';
					//cout<<'\n'<<"scheduling order: "; for(i=0;i<num_cls;i++) cout<<indx[i]<<" ";
				//}
			}
			//cout<<'\n'<<"Qtemp_sorted: "<<endl; for(i=0;i<num_cls;i++) cout<<Q_temp[i]<<" ";
			
			for(i2=0;i2<num_cls;i2++) {
				if(meth==4 || meth==6) a=indx[i2];
				else if(meth==3 || meth==5) a=j_sv[i2]; //for random scheduling (5 uses BCJR for GCN, 3 uses BP)
				
				//if(meth==3) cout<<'\n'<<"a: "<<a<<" ";
				
				//cout<<'\n'<<"Q: "; for(i=0;i<m;i++) cout<<Q_temp[i]<<" ";
				//cout<<'\n'<<"Qmax: "<<Qmax<<" "<<endl;
				//cout<<'\n'<<"i2, a: "<<i2<<", "<<a;	
				
				if(meth==5 || meth==6) 
					clus_decode_bcjr(a,l);
				else clus_decode(a);
										
				//check syndrome
				/*for(i=0;i<n;i++) {
					if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
					else x_hat[cw*n+i]=1;
				}
				for(j=0;j<m;j++) syn[j]=0; //refresh
				for(j=0;j<m;j++) {											
					for(i=0;i<n;i++) 
						syn[j]+=H[j*n+i]*x_hat[cw*n+i]; //finding cw syndrome after each iter.
					syn[j]=syn[j]%2;
				}		
				cout<<'\n'<<"i2: "<<i2; 
				cout<<'\n'<<"pLR: "; for(i=0;i<n;i++) cout<<pLR[0*n+i]<<" "; cout<<'\n'<<'\n';
				cout<<'\n'<<"x_hat: "; for(i=0;i<n;i++) cout<<x_hat[0*n+i]<<" "; cout<<'\n';
				cout<<'\n'<<"syn: ";  for(i=0;i<m;i++) cout<<syn[i]<<" "; cout<<'\n';
									
				flg2=0;
				for(j=0;j<m;j++) 
					if(syn[j]) {
						flg2=1; //error 
						break;
					}
				if(!flg2)
					break;*/ //stopping condition reached
				//cout<<'\n'<<"HERE2 "<<endl;
			}	
		//}
	}
	
	if(prep_dataset && cnt_data<dataset_sz) {
		outf1.open("pLR"+func(fn)+".txt", std::ios_base::app); 
		for(i=0;i<n;i++) 
			outf1<<pLR[0*n+i]<<" "; 
		outf1.close(); 
		
		outf2.open("sch_ordr"+func(fn)+".txt", std::ios_base::app); 
		for(i=0;i<m;i++) 
			outf2<<indx[i]<<" "; 
		outf2.close(); 
		
		cnt_data++;
	}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	//cw=0; cout<<'\n'<<"E_c_v_cnt: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v_cnt[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
}

