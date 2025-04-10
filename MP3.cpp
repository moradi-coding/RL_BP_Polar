
//EDS-LBP alg. of Han-Lui paper
void MP3(int cnt, int l) { 
		
	long a,i,i2,i3,cw,j,j2,k,vidx,cnt2,cidx,cidx2,flg,flg4,J,s,row;
	double tmp,tmp2, prod,val;

	cw=0;
	//for(cw=0;cw<CW;cw++) {
	for(i=0;i<n;i++) 
		for(j=0;j<col_wt;j++) 
			res_v_c_srtd[CW*(col_wt*i+j)+cw]=res_v_c[CW*(col_wt*i+j)+cw]; 
	//cout<<'\n'<<"res_v_c:"<<'\n'; for(i=0;i<n;i++) {for(j=0;j<col_wt;j++) cout<<res_v_c_srtd[CW*(col_wt*i+j)+cw]<<" "; cout<<'\n';}

	for(j=0;j<n;j++) {
		for(i=0;i<col_wt;i++) {
			for(i2=i+1;i2<col_wt;i2++) {
				if(res_v_c_srtd[CW*(col_wt*j+i)+cw]<res_v_c_srtd[CW*(col_wt*j+i2)+cw]) {   
					tmp=res_v_c_srtd[CW*(col_wt*j+i)+cw];
					res_v_c_srtd[CW*(col_wt*j+i)+cw]=res_v_c_srtd[CW*(col_wt*j+i2)+cw];
					res_v_c_srtd[CW*(col_wt*j+i2)+cw]=tmp;
				}			
			}
		}
	}

	for(i=0;i<n;i++) vn_indx[i]=i; //refresh
	for(i=0;i<n;i++) {
		for(i2=i+1;i2<n;i2++) {
			if(res_v_c_srtd[CW*(col_wt*i+0)+cw]<res_v_c_srtd[CW*(col_wt*i2+0)+cw]) {   
				tmp=res_v_c_srtd[CW*(col_wt*i+0)+cw];
				res_v_c_srtd[CW*(col_wt*i+0)+cw]=res_v_c_srtd[CW*(col_wt*i2+0)+cw];
				res_v_c_srtd[CW*(col_wt*i2+0)+cw]=tmp;
					
				tmp2=vn_indx[i];
				vn_indx[i]=vn_indx[i2];
				vn_indx[i2]=tmp2;
			}			
		}
	}
			
	//cout<<'\n'<<"l: "<<l;
	//cout<<'\n'<<"res_max: "; for(i=0;i<n;i++) cout<<res_v_c_srtd[CW*(col_wt*i+0)+cw]<<" "; 
	//cout<<'\n'<<"vn_indx: "; for(i=0;i<n;i++) cout<<vn_indx[i]<<" "; 

	//if(!l) stp=m; else stp=1;
	//stp=1;

	if(!l) a=rand()%n;
	else a=vn_indx[0]; //index of scheduled VN
	//a=l;
	//cout<<'\n'<<"scheduled VN: "<<a;
			
	for(J=0;J<col_wt;J++) { 
		cidx=cns[a*col_wt+J];
		if(cidx>-1) {		
			for(i=0;i<row_wt;i++) { //row_wt is the no. of neighboring VNs of CN cidx	
				vidx=vns[cidx*row_wt+i]; //index of the ith neighboring VN of CN cidx
				//vidx=a;
				
				if(vidx>-1) {
					tmp=1; 
					for(k=0;k<row_wt;k++) 
						if(vns[cidx*row_wt+k]>-1 && k!=i) 
							tmp*=tanh(0.5*E_v_c[CW*(row_wt*cidx+k)+cw]); //jth CN accumulating msgs from all neighboring VNs except vidx
							if(tmp>0.9999)
								tmp=0.9999;
							else if(tmp<-0.9999)
								tmp=-0.9999;
					//cout<<'\n'<<"CN: "<<cidx<<" VN: "<<vidx<<" tmp: "<<tmp<<endl;
					for(k=0;k<col_wt;k++) 
						if(cns[vidx*col_wt+k]==cidx) {
							E_c_v[CW*(vidx*col_wt+k)+cw]=2*atanhf(tmp); //msg sent by CN cidx to ith VN
							E_c_v_cnt[CW*(vidx*col_wt+k)+cw]++; 
							//if(E_c_v_cnt[CW*(vidx*col_wt+k)+cw]>1)
								//cout<<"!!!!"<<'\n';
							break;
						} 				
				}
			} 
		//}
					
		////////////////////////////////////////////////////////
		for(j2=0;j2<row_wt;j2++) { //row_wt is the no. of VNs connected to CN cidx
			i=vns[cidx*row_wt+j2]; 	
			//cout<<'\n'<<"i: "<<i<<endl;					
			if(i>-1) {
				for(j=0;j<col_wt;j++) { //col_wt is the no. of neighboring CNs of VN i	
					cidx2=cns[i*col_wt+j]; //index of jth neighboring CN of VN i			
					if(cidx2>-1) {
						tmp=0; 
						for(k=0;k<col_wt;k++) 
							if(cns[i*col_wt+k]>-1 && k!=j) {
								tmp+=E_c_v[CW*(col_wt*i+k)+cw]; 
							} //ith VN accumulating msg from all neighboring CNs except cidx2 
					
						//cout<<'\n'<<"VN: "<<i<<" CN: "<<cidx2<<" tmp: "<<tmp<<endl;
						//thresholding
						if(tmp>1000)
							tmp=1000;
						else if(tmp<-1000)
							tmp=-1000;
						for(k=0;k<row_wt;k++) 
							if(vns[cidx2*row_wt+k]==i) {
								val=tmp+LR[cw*n+i];	
								//cout<<'\n'<<"val: "<<val<<endl;
								res_v_c[CW*(col_wt*i+j)+cw]=abs((val-E_v_c[CW*(row_wt*cidx2+k)+cw])/(val+E_v_c[CW*(row_wt*cidx2+k)+cw])); //residual update of ith VN
								//if(isnan(res_v_c[CW*(col_wt*i+j)+cw])) 
									//cout<<'\n'<<"val, E_v_c, k: "<<val<<", "<<E_v_c[CW*(row_wt*cidx2+k)+cw]<<", "<<k<<endl;
								E_v_c[CW*(row_wt*cidx2+k)+cw]=val; //msg sent by ith VN to CN cidx2 
								break;
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
		}	
		}
	}
	for(k=0;k<col_wt;k++)
		res_v_c[CW*(col_wt*a+k)+cw]=0;	

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	//cw=0; cout<<'\n'<<"E_c_v_cnt: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v_cnt[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
}

