
//creates GLDPC matrix using a hamming code and a LDPC base matrix

int H_GLDPC() { 
	//lifts the H2 matrix
	int i,j,i2,i3,j2,k,lsft,sum,sum2=0,cnt,flg=1,flg2,ltst_row=0;
	
	//create the GLDPC matrix
	for(i=0;i<m;i++) { //row indices of the base matrix
		cnt=0;
		
		//while(cnt2<=num_spc) {
		flg2=0; 
		for(k=0;k<num_gcn;k++) 
			if(i==spc_idx[k]) {
				flg2=1; //indices of single parity CNs which should be replaced with generalized CNs
				break;
			}
			else 
				flg2=0; 
		
		//if(i<num_spc) {
		if(flg2) {
			for(j=0;j<n;j++) {
				if(Hbase[i][j]) {		
					//replace 1 by a column of a m_sub x n_sub sub-code parity-check matrix
					//for(i2=i*m_sub;i2<(i+1)*m_sub;i2++) {
					for(i2=ltst_row;i2<ltst_row+m_sub;i2++) {
						if(type==1 || (type==2 &&  i<m/2))	
							Hlift[i2][j]=hamming[i2-ltst_row][cnt]; 
						else
							Hlift[i2][j]=hamming2[i2-ltst_row][cnt]; 
						cns_subcode[i][i2-ltst_row]=i2;											
					}	
					cnt++;	
				}
				else {				
					//replace 0 by a zero column vector of length m_sub
					//for(i2=i*m_sub;i2<(i+1)*m_sub;i2++) 
					for(i2=ltst_row;i2<ltst_row+m_sub;i2++) 
						Hlift[i2][j]=0; 	
				}
			}
			ltst_row=i2; //latest row no.
		}
		else { //just use single parity CN
			for(i2=ltst_row;i2<ltst_row+m_sub;i2++) {
				for(j=0;j<n;j++) {
					if(i2==ltst_row) {
						Hlift[i2][j]=Hbase[i][j]; 
						cns_subcode[i][i2-ltst_row]=i2;	
					}
					else
						cns_subcode[i][i2-ltst_row]=-1;
				}										
			}
			ltst_row++;
		}
		//cout<<'\n'<<"cnt: "<<cnt;
			//cout<<'\n'<<"ltst_row: "<<ltst_row;
	}

	
	//finding vns_subcode
	for(i=0;i<tot_r;i++) {	
		//if(i==i/m_sub*m_sub) 
		cnt=0;
		for(j=0;j<n;j++) {
			if(Hlift[i][j]) {
				for(i2=0;i2<row_wt;i2++)
					if(j!=vns_subcode[i/m_sub][i2]) flg=1;
					else {flg=0; break;}
						
				if(flg) {
					//cout<<"Hlift "<<Hlift[i2][j];
					vns_subcode[i/m_sub][cnt]=j;	
					cnt++;
				}
			}
		}
	}

}
