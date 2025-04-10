
#include "init.cpp"
#include "th.cpp"
#include "MP0.cpp"
#include "MP0_bcjr.cpp"
//#include "MP2.cpp"
#include "MP9.cpp"

int BP(int cnt2, int num_dat/*, const auto model0*/) {
		
	long i,i2,k,cw,j,l,strt=0,stp,iter,flg2,flg3,flg4,flg5,flg6,cnt,cnt3,iter_out;

	//E_v_c[j][i][cw]: msg accumulated by CN j from VN i for codeword stream cw, E_v_c[j][i][cw]=E_v_c[CW*(no. of cols *j + i)+cw]
	//E_c_v[i][j][cw]: msg accumulated by VN i from CN j for codeword stream cw
	
	//refresh
	for(i=0;i<CW;i++) excl_cw[i]=-1;

	num=0;
	//while(num<=(L-W)/(mem+1)) {
		init1();
		init2(); 
		//initialization
		//for(cw=0;cw<CW;cw++) for(j=0;j<m;j++) for(k=0;k<row_wt;k++) {
			//E_v_c_mu[CW*(row_wt*j+k)+cw]=0;
		//}

		cnt=0; flg4=1;
		cw=0; //for(cw=0;cw<CW;cw++) 
		/*for(i=0;i<m;i++) 
			for(j=0;j<row_wt;j++) 
				res_c_v[CW*(row_wt*i+j)+cw]=0; //refresh residual
		for(i=0;i<n;i++) 
			for(j=0;j<col_wt;j++) 
				res_v_c[CW*(col_wt*i+j)+cw]=1;*/ //refresh residual

		//if(meth==2) iter=m*Ifl; //for NS 
		//else 
		iter=Ifl; //for others
		if((fn==-1 || fn==-4) && meth==2)
			iter=1;
		
		//cout<<'\n'<<"iter: "<<iter<<endl;
		//cout<<"here0"<<endl;
		//flg6=1; cnt3=0; llr_max=1000;
		//while(flg6) {
		for(l=1;l<=iter;l++) {
			//cout<<"l: "<<l<<'\n';
			if(meth==1) MP0(cnt); //flooding	
			else if(meth==2) MP0_bcjr(cnt); //flooding with BCJR
			else if(meth==3 || meth==4) MP9(cnt,l/*,model0*/); //RL/MRL stuff where meth 3 is random scheduling
			else if(meth==5) MP9(cnt,l); //random seq. with BCJR
			else if(meth==6) MP9(cnt,l); //RL seq. with BCJR
			
			//checking syndromes	
			for(cw=0;cw<CW;cw++) {
				for(j=0;j<m;j++) syn[j]=0; //refresh
				for(k=0;k<cnt;k++) 
					if(cw!=excl_cw[k]) flg4=1;
					else {flg4=0; break;}

				if(flg4) { //if cw was not recovered previously
					for(i=0;i<n;i++) {
						if(pLR[cw*n+i]>=0) x_hat[cw*n+i]=0; 
						else //if(pLR[cw*n+i]<0)
							x_hat[cw*n+i]=1;
						//cout<<x_hat[cw*n+i]<<" ";
					}		
					//check syndrome
					for(j=0;j<m;j++) {											
						for(i=0;i<n;i++) 
							syn[j]+=H[j*n+i]*x_hat[cw*n+i]; //finding cw syndrome after each iter.
						syn[j]=syn[j]%2;
					}	
										
					flg2=0;
					for(j=0;j<m;j++) 
						if(syn[j]) {
							flg2=1; //error 
							break;
						}	
											
					if(!flg2) { //CW has been recovered 
						excl_cw[cnt]=cw;
						cnt++;
					}
				}				
			}		
			iter_out=l;
			if(cnt==CW) break; //all cws were recovered	
			//cout<<'\n'<<"l: "<<l<<endl;
			//cout<<'\n'<<"pLR: "; for(cw=0;cw<CW;cw++) {for(i=0;i<n;i++) cout<<pLR[cw*n+i]<<" "; cout<<'\n'<<'\n';}	
		}
		//cout<<'\n'<<"l: "<<l<<endl;
		//if(!flg2 && !l) l++;	
		
			/*if((meth==2 || meth==5 || meth==6) && flg2 && cnt3<25) { //when doing BCJR, do adaptive thresholding
				llr_max+=1000;
				cnt3++;
				flg6=1;
			}
			else
				flg6=0;
		}*/

		//hard decision
		//if(num<(L-W)/(mem+1) || BC) STP=1;
		//else STP=W/(mem+1); //decision in the last window position

		/*for(j=0;j<STP;j++) {
			if(num<(L-W)/(mem+1))
				strt=num*hshft;
			else if(num==0 && j==0) 
				strt=0;
			else strt+=hshft;

			stp=strt+hshft;*/
			//cout<<'\n'<<"strt: "<<strt<<" stp: "<<stp<<endl;
			if(fn==13) //for polar code
				strt=44;
			else if(fn==14) //for polar code
				strt=72;
			else if(fn==15) //for polar code
				strt=251;
			else if(fn==16) //for polar code
				strt=52;
			else if(fn==17) //for polar code
				strt=123;
			else if(fn==18) //for polar code
				strt=45;
			else if(fn==19) //for polar code
				strt=52;
			else if(fn==20 || fn==21) //for polar code
				strt=127;
			else if(fn==22) //for polar code
				strt=131;
			else if(fn==23) //for polar code
				strt=47;
			else if(fn==24) //for polar code
				strt=46;
			else if(fn==25) //for polar code
				strt=51;
			else if(fn==26) //for polar code
				strt=121;
			else if(fn==27) //for polar code
				strt=275;
			else if(fn==28) //for polar code
				strt=127;
			else if(fn==29) //for polar code
				strt=109;
			else
				strt=0;
				
			cw=0;
			for(i=0;i<n;i++)
				if(isnan(pLR[cw*n+i])) {
					nan_cnt++;
					flg5=1;
					break;
				}
				else
					flg5=0;
					
			for(cw=0;cw<CW;cw++) {
				flg3=0;
				//for(i=strt;i<stp;i++) {
				for(i=strt;i<n;i++) {
					if(!flg5 && pLR[cw*n+i]<0) {
						if(!flg3) {
							if(meth==1) err++; 
							else if(meth==2) err2++;
							else if(meth==3) err3++;
							else if(meth==4) err4++;
							else if(meth==5) err5++;
							else if(meth==6) err6++;
							else if(meth==7) err7++;
							else if(meth==8) err8++;
							flg3=1;
						}
						if(meth==1) biterr++;  //for bit err
						else if(meth==2) biterr2++; 
						else if(meth==3) biterr3++;
						else if(meth==4) biterr4++;
						else if(meth==5) biterr5++;
						else if(meth==6) biterr6++;
						else if(meth==7) biterr7++;
						else if(meth==8) biterr8++; 
						
					}
 					//if(pLR[cw*n+i-strt]>=0) 
						//x_hat[cw*n+i]=0; 
					//else 
						//x_hat[cw*n+i]=1;
				}
				if(!flg5) {
					if(meth==1) dcw++; //no. of decoded cws
					else if(meth==2) dcw2++;
					else if(meth==3) dcw3++;
					else if(meth==4) dcw4++;
					else if(meth==5) dcw5++;
					else if(meth==6) dcw6++;
					else if(meth==7) dcw7++;
					else if(meth==8) dcw8++;
				}

			}
		
			//cout<<'\n'<<"pLR: "<<'\n'; for(cw=0;cw<CW;cw++){ for(i=0;i<n;i++) cout<<pLR[i]<<" "; cout<<'\n'<<'\n';}
			/*if(flg3) {
				cout<<'\n'<<"syn: ";  for(i=0;i<m;i++) cout<<syn[i]<<" "; cout<<'\n';
				cout<<'\n'<<"x_hat: "; for(cw=0;cw<CW;cw++){ for(i=strt;i<stp;i++) cout<<x_hat[cw*n+i-strt]<<" "; cout<<'\n';}
				cout<<'\n'<<"l: "<<l; 
			}*/
			
			//to check what fraction of bits are in error
			if(flg3 && !flg5) { //&& i2==num_dat-1
				//cout<<"abs_sz: "<<abs_sz<<endl;
				//outf.open("xhat_"+func(fn)+"_"+func(i2)+".txt", std::ios_base::app); 
				//ofstream outf,outf2; 
				/*outf.open("xhat_"+func(fn)+".txt"); //outf2.open("pLR"+func(fn)+".txt"); 
				for(i=0;i<n;i++) {
					//outf<<x_hat[i]<<" "; 
					if(x_hat[i]) {
						outf<<i<<" ";
						//outf2<<pLR[i]<<" ";
					}
				}
				outf<<std::endl; //outf2<<std::endl;
				outf.close();*/ //outf2.close();
				
				/*cout<<"CN index: "<<'\n';
				for(i2=0;i2<n;i2++)
					for(i=0;i<num_cls;i++) 
						for(j=0;j<num_vns_cls;j++) 
							if(x_hat[i2] && i2==vns_cluster[i][j]) {
								cout<<i<<" ";
								break;
							}*/
							
				/*ofstream outf; 
				outf.open("y_mat.txt", std::ios_base::app); 
				for(i=0;i<n;i++) 
					outf<<y[i]<<" "; outf<<std::endl; 
				outf.close();*/
				
				/*outf.open("llrs_bcjr.txt", std::ios_base::app); 
				outf<<llr_max<<" "; 
				outf<<std::endl; 
				outf.close();*/
				
				//}
			}

		//}
		num++; //no. of window shifts
	//}
	
	/*if(!flg3 && !flg5) {
		ofstream outf; outf.open("y.txt", std::ios_base::app); for(i=0;i<n;i++) outf<<y[i]<<" "; outf<<std::endl; outf.close();
		not_err++;
		cout<<'\n'<<"not_err: "<<not_err;
	}*/
	
	//if(flg3 && !flg5) {		
		/*cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
		cw=0; cout<<'\n'<<"E_v_c: "<<'\n'; for(j=0;j<m;j++){for(k=0;k<row_wt;k++) cout<<E_v_c[CW*(row_wt*j+k)+cw]<<" "; cout<<'\n';}
		*/
		//cout<<'\n'<<"y: "<<'\n'; for(j=0;j<n;j++)  cout<<y[j]<<" "; cout<<'\n';
	//}
	//cout<<'\n'<<"pLR: "; for(cw=0;cw<CW;cw++) {for(i=0;i<n;i++) cout<<pLR[cw*n+i]<<" "; cout<<'\n';}
	return iter_out;
}
