
#include "bcjr_log_15_11.cpp"
#include "bcjr_log2.cpp"
//#include "bcjr.cpp"

//flooding schedule
void MP0_bcjr(int cnt) { 
	long a,J,i,i2,cw=0,j,k,vn,vidx,cidx,flg=0,flg2,cnt2=0;
	double tmp;	
	
	//for(cw=0;cw<CW;cw++) {
		/*flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;}*/ //if flg4=0, cw has been recovered already
		//if(flg4) {
		
		//horizontal step
		for(a=0;a<num_cls;a++) { //num_cls is the total number of CNs
			for(k=0;k<num_cls-num_gcn;k++) 
				if(a==spc_indx[k]) {
					flg=1;
					break;
				}
				else
					flg=0;
			
			//cout<<'\n'<<"flg: "<<flg<<endl;	
			if(flg) { //if 'a' is SPCN
				cnt2=0; flg2=1;
				for(i2=0;i2<row_wt*cls_sz;i2++)
					vn_sv[i2]=-1;
				for(J=0;J<cls_sz;J++) { 
					j=cns_cluster[a][J]; //j is a SPCN of CN a
					if(j>-1) {
						for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j	
							vidx=vns[j*row_wt+i]; //index of the ith neighboring VN of CN j	
							for(i2=0;i2<cnt2;i2++)
								if(vidx==vn_sv[i2]) {
									flg2=0;
									break;
								}
								else 
									flg2=1;					
							if(vidx>-1) {
								tmp=1; 
								for(k=0;k<row_wt;k++) 
									if(vns[j*row_wt+k]>-1 && k!=i) 
										tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except i
									//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
									for(k=0;k<col_wt;k++) 
										if(cns[vidx*col_wt+k]==j) {
											E_c_v[CW*(vidx*col_wt+k)+cw]=th3(2*atanhf(tmp)); //msg sent by jth CN to ith VN
											if(isnan(E_c_v[CW*(vidx*col_wt+k)+cw]))
												cout<<E_c_v[CW*(vidx*col_wt+k)+cw]<<" ";
											//E_c_v_cnt[CW*(vidx*col_wt+k)+cw]++; 
											if(flg2) {
												E_c_v_cnt2[a]++; 
												vn_sv[cnt2]=vidx;
												cnt2++;
											}
											break;
										} 				
							}
						}
					}
				}
			}
			else { //use BCJR for CN computation
				cnt2=0; flg2=1;
				for(i2=0;i2<row_wt*cls_sz;i2++)
					vn_sv[i2]=-1;
				for(J=0;J<cls_sz;J++) {
					j=cns_cluster[a][J]; //a SPCN of CN 'a'
					//cout<<'\n'<<"j: "<<j<<endl;
					if(j>-1) {
						for(k=0;k<row_wt;k++) {
							vn=vns[j*row_wt+k];
							for(i2=0;i2<cnt2;i2++)
								if(vn==vn_sv[i2]) {
									flg2=0;
									break;
								}
								else 
									flg2=1;
							if(vn>-1 && flg2) {
								L_in[vn]=E_v_c[CW*(row_wt*j+k)+cw]; //input to BCJR
								vn_sv[cnt2]=vn;
								cnt2++;
								//cout<<"E_v_c: "<<E_v_c[CW*(row_wt*j+k)+cw];
							}
						}
					}
				}	
				//cout<<"cnt2: "<<cnt2;
				//cout<<'\n'<<"vn_sv: "; for(i=0;i<num_vns_cls;i++) cout<<vn_sv[i]<< " "; cout<<'\n';
				//cout<<'\n'<<"L_in: "; for(i=0;i<n;i++) cout<<L_in[i]<< " "; //cout<<'\n';
				//if(!l) {cout<<"LLR: "; for(i=0;i<num_vns_cls;i++) cout<<LR[vn_sv[i]]<< " "; cout<<'\n';}
				//cout<<"L_in: "; for(i=0;i<num_vns_cls;i++) cout<<L_in[vn_sv[i]]<< " "; cout<<'\n';
				//bcjr(a); //do BCJR decoding on cluster
				if(fn==-4 || fn==-5 || fn==11) 
					bcjr_log_15_11(a);
				else 
					bcjr_log(a);
				//cout<<"L_out: "; for(i=0;i<num_vns_cls;i++) cout<<L_out[vn_sv[i]]<< " "; cout<<'\n';
				//cout<<'\n'<<"L_out: "; for(i=0;i<n;i++) cout<<L_out[i]<< " "; cout<<'\n';
				
				cnt2=0; flg2=1;	
				for(i2=0;i2<row_wt*cls_sz;i2++)
					vn_sv[i2]=-1;	
				for(J=0;J<cls_sz;J++) {	
					j=cns_cluster[a][J]; //a VN of the CN will receive identical messages from its neighbors
					if(j>-1) {
					for(i=0;i<row_wt;i++) { 
						vidx=vns[j*row_wt+i]; //index of the ith neighboring VN of CN j
						for(i2=0;i2<cnt2;i2++)
							if(vidx==vn_sv[i2]) {
								flg2=0;
								break;
							}
							else 
								flg2=1;
						if(vidx>-1 && flg2) {
							//tmp=1; 		
							//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
							for(k=0;k<col_wt;k++) 
								if(cns[vidx*col_wt+k]==j) {
									E_c_v[CW*(vidx*col_wt+k)+cw]=L_out[vidx]-L_in[vidx]; //E_v_c[CW*(row_wt*j+i)+cw]; //msg sent by jth CN to only 1 VN of its SPCN 
									//E_c_v_cnt[CW*(vidx*col_wt+k)+cw]++; 
									E_c_v_cnt2[a]++; 
									vn_sv[cnt2]=vidx;
									cnt2++;
									//cout<<"vn_sv: "; for(i2=0;i2<num_vns_cls;i2++) cout<<vn_sv[i2]<< " "; cout<<'\n';
									break;
								}	
						}	 			
					}
					}
					
				}
			}
		}	
		//cout<<'\n'<<"cns "<<'\n'; for(i=0;i<n;i++) {for(j=0;j<col_wt;j++) cout<<cns[i*col_wt+j]<< " "; cout<<'\n';}
		//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
		
		//vertical step
		for(i=0;i<n;i++) { //n is no. of VNs
			for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i
				cidx=cns[i*col_wt+j];  //index of jth neighboring CN of VN i
				if(cidx>-1) {
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[i*col_wt+k]>-1 && k!=j) 
							tmp+=E_c_v[CW*(col_wt*i+k)+cw]; //ith VN accumulating msg from all neighboring CNs except j 
					
					for(k=0;k<row_wt;k++) {
						if(vns[cidx*row_wt+k]==i && E_c_v[CW*(i*col_wt+j)+cw]) {
							E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; //msg sent by ith VN to jth CN
							E_v_c_cnt[CW*(row_wt*cidx+k)+cw]++;
							break;
						} 
						else if(vns[cidx*row_wt+k]==i && !E_c_v[CW*(i*col_wt+j)+cw]) {
							E_v_c[CW*(row_wt*cidx+k)+cw]=0;
							break;
						}
					}
				}
			} 

			//updating the aposteriori LLR
			tmp=0; 
			for(k=0;k<col_wt;k++) 
				if(cns[col_wt*i+k]>-1) 
					tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
			pLR[cw*n+i]=LR[cw*n+i]+tmp; 
		}

		//}
	//}
	
	//cout<<"L_out: "; for(i=0;i<n;i++) cout<<L_out[i]<< " "; cout<<'\n'<<'\n';
	//cout<<'\n'<<"vns "<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<vns[i*row_wt+j]<< " "; cout<<'\n';}
	//{cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}}
	//cout<<'\n'<<"l: "<<l<<endl;
}

