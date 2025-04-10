
int clus_decode(long a) { //a is index of scheduled cluster
		
	long i,i2,I,cw=0,j,k,vidx,cidx,J,flg2,cnt2=0;
	double tmp;
	
	cnt2=0; flg2=1;
	for(i2=0;i2<row_wt*cls_sz;i2++)
		vn_sv[i2]=-1;
					
	for(I=0;I<clus_iter;I++) {
		for(J=0;J<cls_sz;J++) { 
			//cout<<'\n'<<"J: "<<J<<endl;
			//cout<<'\n'<<"cns_cluster[0][0]: "<<cns_cluster[0][0]<<endl;
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
			
						//cout<<'\n'<<"j: "<<j<<" vidx: "<<vidx<<" tmp: "<<tmp<<endl;
						for(k=0;k<col_wt;k++) 
							if(cns[vidx*col_wt+k]==j) {
								E_c_v[CW*(vidx*col_wt+k)+cw]=th3(2*atanhf(tmp)); //msg sent by jth CN to ith VN
								//if(learn) {
									if(isnan(E_c_v[CW*(vidx*col_wt+k)+cw]))
										cout<<E_c_v[CW*(vidx*col_wt+k)+cw]<<" ";
								//}
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
					
		for(J=0;J<num_vns_cls;J++) { 
			i=vns_cluster[a][J]; //i is a VN of cluster a					
			if(i!=-1) {
				for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i	
					cidx=cns[i*col_wt+j]; //index of jth neighboring CN of VN i			
					if(cidx>-1) {
						tmp=0; 
						for(k=0;k<col_wt;k++) 
							if(cns[i*col_wt+k]>-1 && k!=j) 
								tmp+=E_c_v[CW*(col_wt*i+k)+cw]; //ith VN accumulating msg from all neighboring CNs except j 
						
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
}

