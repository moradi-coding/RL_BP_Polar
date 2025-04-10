
int clus_decode_bcjr(long a, int l) { //a is index of scheduled cluster
	long i,I,i2,cw=0,j,k,vn,vidx,cidx,J,flg=0,flg2,cnt2=0;
	double tmp;
		
	for(k=0;k<num_cls-num_gcn;k++) 
		if(a==spc_indx[k]) {
			flg=1;
			break;
		}
		else
			flg=0;
	//cout<<'\n'<<"a, flg: "<<a<<", "<<flg<<endl;

	for(I=0;I<clus_iter;I++) {
	if(flg) { //if 'a' is SPCN
		cnt2=0; flg2=1;
		for(i2=0;i2<row_wt*cls_sz;i2++)
			vn_sv[i2]=-1;
		for(J=0;J<cls_sz;J++) { 
			//cout<<'\n'<<"a: "<<a<<endl;
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
								if(!learn && flg2) {
									E_c_v_cnt2[a]++; 
									vn_sv[cnt2]=vidx;
									cnt2++;
								}
								//E_c_v_cnt[CW*(vidx*col_wt+k)+cw]++; 
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
		//cout<<'\n'<<"L_in: "; for(i=0;i<n;i++) cout<<L_in[i]<< " "; cout<<endl;
		//if(!l) {cout<<"LLR: "; for(i=0;i<num_vns_cls;i++) cout<<LR[vn_sv[i]]<< " "; cout<<'\n';}
		if(fn==-4 || fn==-5 || fn==11) 
			bcjr_log_15_11(a);
		else 
			bcjr_log(a);
		//cout<<"L_out: "; for(i=0;i<num_vns_cls;i++) cout<<L_out[vn_sv[i]]<< " "; cout<<'\n';
		//cout<<'\n';
				
		cnt2=0; flg2=1;	
		for(i2=0;i2<row_wt*cls_sz;i2++)
			vn_sv[i2]=-1;	
		for(J=0;J<cls_sz;J++) {	
			//cout<<"J1: "<<J<<endl;
			j=cns_cluster[a][J]; //a VN of the GCN will receive identical messages from its neighbors
			if(j>-1) {
				for(i=0;i<row_wt;i++) { 
					//cout<<"i: "<<i<<endl;
					//cout<<"j: "<<j<<endl;
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
								E_c_v[CW*(vidx*col_wt+k)+cw]=L_out[vidx]-L_in[vidx]; //E_v_c[CW*(row_wt*j+i)+cw]; //msg sent by jth GCN to only 1 VN of its SPCN
								if(!learn) 
									//E_c_v_cnt[CW*(vidx*col_wt+k)+cw]++; 
									E_c_v_cnt2[a]++; 
								vn_sv[cnt2]=vidx;
								cnt2++;
								break;
							}	
							//cout<<"cnt2: "<<cnt2<<endl;
							//cout<<"vn_sv: "; for(i2=0;i2<row_wt*cls_sz;i2++) cout<<vn_sv[i2]<< " "; cout<<'\n';	
					}	 			
				}
				//cout<<"i loop completed"<<endl;
			}
		}
	}
	//}	
					
	for(J=0;J<num_vns_cls;J++) { 
		//cout<<"J2: "<<J<<endl;
		i=vns_cluster[a][J]; //i is a VN of cluster a					
		if(i!=-1) {
			for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i	
				cidx=cns[i*col_wt+j]; //index of jth neighboring CN of VN i			
				if(cidx>-1) {
					tmp=0; 
					for(k=0;k<col_wt;k++) 
						if(cns[i*col_wt+k]>-1 && k!=j) {
							tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 				
						} //ith VN accumulating msg from all neighboring CNs except j 
						
					for(k=0;k<row_wt;k++) {
							//if(vns[cidx*row_wt+k]==i) {
						if(vns[cidx*row_wt+k]==i && E_c_v[CW*(i*col_wt+j)+cw]) {
							E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; //msg sent by ith VN to jth CN
							E_v_c_cnt[CW*(row_wt*cidx+k)+cw]++;
							break;
						} 
						else if(l>1 && vns[cidx*row_wt+k]==i && !E_c_v[CW*(i*col_wt+j)+cw]) {
							E_v_c[CW*(row_wt*cidx+k)+cw]=0;
							break;
						}
					}
						
					/*for(k=0;k<row_wt;k++) 
						//if(vns[cidx*row_wt+k]==i) {
						if((!flg && vns[cidx*row_wt+k]==i && E_c_v[CW*(i*col_wt+j)+cw]) || (flg && vns[cidx*row_wt+k]==i)) {
							E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; 
							break;
						}*/ //msg sent by ith VN to jth CN			
				}
			} 
				
			//updating the aposteriori LLR
			if(I==clus_iter-1) {
				tmp=0; 
				for(k=0;k<col_wt;k++) 
					if(cns[col_wt*i+k]>-1) 
						tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
				pLR[cw*n+i]=LR[cw*n+i]+tmp; 
			}
		}	
	}
	}
	//cout<<'\n'<<"l, a: "<<l<<", "<<a<<endl;
	//cout<<"L_out: "; for(i=0;i<n;i++) cout<<L_out[i]<< " "; cout<<'\n';
	//cout<<'\n'<<"vns "<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<vns[i*row_wt+j]<< " "; cout<<'\n';}
	//{cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}}
}

