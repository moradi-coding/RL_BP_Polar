
//flooding schedule
void MP0(int cnt) { 
	long a,J,i,i2,cw,j,k,vidx,cidx,flg4,flg2,cnt2=0;
	double tmp;
	
	for(cw=0;cw<CW;cw++) {
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, cw has been recovered already
		
		if(flg4) {
			//horizontal step
			/*for(j=0;j<m;j++) { //m is no. of CNs
				for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN j
	
					vidx=vns[j*row_wt+i]; //index of ith neighboring VN of CN j
	
					if(vidx>-1) {
						tmp=1; 
						for(k=0;k<row_wt;k++) 
							if(vns[j*row_wt+k]>-1 && k!=i) tmp*=tanh(0.5*E_v_c[CW*(row_wt*j+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except i
						
						//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
			
						for(k=0;k<col_wt;k++) 
							if(cns[vidx*col_wt+k]==j) {
								E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); //msg sent by jth CN to ith VN
								E_c_v_cnt[CW*(vidx*col_wt+k)+cw]++; 
								//if(E_c_v_cnt[CW*(vidx*col_wt+k)+cw]>Ifl)
									//cout<<"!!!!"<<'\n';
								break;
							} 				
					}
				} 
			}*/
			
			for(a=0;a<num_cls;a++) {
				cnt2=0; flg2=1;	
				for(i2=0;i2<row_wt*cls_sz;i2++)
					vn_sv[i2]=-1;	
				for(J=0;J<cls_sz;J++) { 
					//cout<<'\n'<<"J: "<<J<<endl;
					j=cns_cluster[a][J]; //j is a CN of cluster a
					//cout<<'\n'<<"j: "<<j<<endl;
								
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
											
								//if(tmp==1) {cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;}
								for(k=0;k<col_wt;k++) 
									if(cns[vidx*col_wt+k]==j) {
										//E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); 
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

			//vertical step
			for(i=0;i<n;i++) { //n is no. of VNs
				for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i
					cidx=cns[i*col_wt+j];  //index of jth neighboring CN of VN i
	
					if(cidx>-1) {
						tmp=0; 
						for(k=0;k<col_wt;k++) 
							if(cns[i*col_wt+k]>-1 && k!=j) {
								tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
							}//ith VN accumulating msg from all neighboring CNs except j 
		
						for(k=0;k<row_wt;k++) 
							if(vns[cidx*row_wt+k]==i) {
								E_v_c[CW*(row_wt*cidx+k)+cw]=tmp+LR[cw*n+i]; 
								E_v_c_cnt[CW*(row_wt*cidx+k)+cw]++;
								break;
							} //msg sent by ith VN to jth CN
							//printf("tmp2=%f\n",tmp2);					
					}
				} 

				//updating the aposteriori LLR
				tmp=0; 
				for(k=0;k<col_wt;k++) 
					if(cns[col_wt*i+k]>-1) tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
				pLR[cw*n+i]=LR[cw*n+i]+tmp; 

			}

		}
	}

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	//{cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}}
	//cout<<'\n'<<"pLR: "; for(cw=0;cw<CW;cw++) {for(i=0;i<n;i++) cout<<pLR[cw*n+i]<<" "; cout<<'\n'<<'\n';}
	
}


