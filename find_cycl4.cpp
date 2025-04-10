
void find_cycl() {  //cycle length is always even
	long i,i2,i3,j,r,c,l=1,ver=0,hrz=1,rw_wt=0,bck=0,mtch=0,flg=1,flg2=0,flg3=0,flg4=0,flg5=0,flg6=0,rev=0,of=0,strt_L=1,par=0,dnodes=2,dnodes_chk=1,K,idx1,idx2,cnt=0;	

	for(i=col_strt;i<col_stp;i++) if(Hsc[strt_r][i]) rw_wt++; if(rw_wt<2) return;
	snk_r[l-1]=strt_r; snk_c[l-1]=stp_c=strt_c;

	while(l<cyc_len) { 
		//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------hrz
		if(l<cyc_len && hrz) {
			r=snk_r[l-1-of];
			if(l==cyc_len-1 && !Hsc[r][stp_c]) bck=1; 
			else if(par && l==cyc_len-2 && r==stp_r) bck=1; //if this happens, snk can't do ver on last move
			if(!bck) {
				for(i=col_strt;i<col_stp;i++) {//col srch
					//if(col_wt[i]==3 && Hsc[r][i] && (l==cyc_len-1 || (l<cyc_len-1 && i!=strt_c)) ) { //checks col. wt. and prevents sub cycle
					if(Hsc[r][i] && (l==cyc_len-1 || (l<cyc_len-1 && i!=strt_c)) ) { 
						//checking forbidden path
						for(i2=rev-1;i2>=0;i2--) {
							if(fbd_r[i2][l-1-of]==snk_r[l-1-of] && fbd_c[i2][l-1-of]==snk_c[l-1-of] && 
								fbd_r[i2][l-of]==snk_r[l-1-of] && fbd_c[i2][l-of]==i) {
								for(i3=1;i3<cyc_len-1;i3++) {
									idx1=l-1-of-i3; 
									if(idx1>=0) idx2=idx1; else idx2=cyc_len-1+idx1; //causes roll over
									if(idx2!=l-of && fbd_r[i2][idx2]==snk_r[idx2] && fbd_c[i2][idx2]==snk_c[idx2] && 
										fbd_r[i2][idx2+1]==snk_r[idx2+1] && fbd_c[i2][idx2+1]==snk_c[idx2+1]) {
										flg2=1;  
									}
									else if(idx2==l-of && fbd_r[i2][idx2]==snk_r[l-1-of] && fbd_c[i2][idx2]==i && 
										fbd_r[i2][idx2+1]==snk_r[idx2+1] && fbd_c[i2][idx2+1]==snk_c[idx2+1]) {
										flg2=1;  
									}
									else {flg2=0; break;}
								}
							}
							else flg2=0;
							if(flg2) break; //forb. pth exists
						} 
						//cout<<'\n'<<"flg2_h: "<<flg2; 
						//cout<<'\n'<<"forb. pth: "; for(i2=0;i2<rev;i2++) cout<<" "<<frm_r[i2]<<","<<frm_c[i2]<<" to "<<to_r[i2]<<","<<to_c[i2]; cout<<'\n';
						if(!flg2) { //if no fbd path 
							//if(chld_num==2) {cout<<'\n'<<'\n'<<"snake: "; for(i2=0;i2<cyc_len;i2++) cout<<" "<<snk_r[i2]<<","<<snk_c[i2]<<"  "; }
							if(l>=1 || (par && l>=1)) {for(i2=0;i2<cyc_len;i2++) if(r==snk_r[i2] && i==snk_c[i2]) {flg=0; break;} else flg=1; }//snake can't hit itself
							//*if(chld_num==3)*/ cout<<'\n'<<"flg: "<<flg;
							if((flg && l<cyc_len-1) || (flg && l==cyc_len-1 && i==stp_c)) {
								/*if(l==1) cnt++;*/ //wrong
								l++; snk_r[l-1-of]=r; snk_c[l-1-of]=i; ver=1; hrz=0; flg=0; break;
							}
						}
					}
					//cout<<"col. no.: "<<i<<endl;
				}
			}
			if(par && l==strt_L) flg3=1; else flg3=0; //need to change the position of dnodes
			if(!flg3 && hrz) bck=1; else bck=0; //if search is allowed but still in hrz mode, snake has to retreat
			//cout<<'\n'<<"bck: "<<bck;
			//cout<<'\n'<<"snake_h1: "; for(i2=0;i2<cyc_len;i2++) cout<<" "<<snk_r[i2]<<","<<snk_c[i2]<<"  ";
			if(bck) { //cannot retreat from the parent head node
				bck=0;
				if(l==1) {cnt++; hrz=1; ver=0; if(cnt==rw_wt-1) break;} //no cycle can be found from the starting point
				else {
					for(i=0;i<cyc_len;i++) {fbd_r[rev][i]=snk_r[i]; fbd_c[rev][i]=snk_c[i];}  //forbidden path (the entire snake trajectory is a single fb. pth.)
					snk_r[l-1-of]=snk_c[l-1-of]=-1; //erasing retreated point
					l--; rev++; ver=1; hrz=0; 
				} 
			}			
			//*if(chld_num>=3){*/cout<<'\n'<<"L_hrz:"<<l; cout<<" "<<snk_r[l-2-of]<<","<<snk_c[l-2-of]<<" to "<<snk_r[l-1-of]<<","<<snk_c[l-1-of];//}	
			//if(chld_num>=3) {cout<<'\n'<<"fbd_pth: "<<'\n'; for(i2=0;i2<rev;i2++){for(i3=0;i3<cyc_len;i3++) cout<<" "<<fbd_r[i2][i3]<<","<<fbd_c[i2][i3]; cout<<'\n';}}
			//cout<<'\n'<<"snake_h2: "; for(i2=0;i2<cyc_len;i2++) cout<<" "<<snk_r[i2]<<","<<snk_c[i2]<<"  ";
			//cout<<'\n'<<"rev: "<<rev;	
		}

		//--------------------------------------------------------------------------------------------------------------------------------------------------------------------ver
		if(!flg3 && l<cyc_len && ver) {
			c=snk_c[l-1-of];
			if(par && l==cyc_len-1 && !Hsc[stp_r][c]) bck=1; 
			else if(par && l==cyc_len-2 && c==stp_c) bck=1; //if this happens, snk can't do hrz on last move
			//cout<<'\n'<<"rev: "<<rev;
			if(!bck) {
				for(i=row_strt;i<row_stp;i++) {//row srch
					if(Hsc[i][c] && i!=strt_r) {
						for(i2=rev-1;i2>=0;i2--) { //checking forbidden path
							if(fbd_r[i2][l-1-of]==snk_r[l-1-of] && fbd_c[i2][l-1-of]==snk_c[l-1-of] && 
								fbd_r[i2][l-of]==i && fbd_c[i2][l-of]==snk_c[l-1-of]) {
								for(i3=1;i3<cyc_len-1;i3++) {
									idx1=l-1-of-i3; 
									if(idx1>=0) idx2=idx1; else idx2=cyc_len-1+idx1; //causes roll over
									if(idx2!=l-of && fbd_r[i2][idx2]==snk_r[idx2] && fbd_c[i2][idx2]==snk_c[idx2] && 
										fbd_r[i2][idx2+1]==snk_r[idx2+1] && fbd_c[i2][idx2+1]==snk_c[idx2+1]) {
										flg2=1;  
									}
									else if(idx2==l-of && fbd_r[i2][idx2]==i && fbd_c[i2][idx2]==snk_c[l-1-of] && 
										fbd_r[i2][idx2+1]==snk_r[idx2+1] && fbd_c[i2][idx2+1]==snk_c[idx2+1]) {
										flg2=1;  
									}
									else {flg2=0; break;}
								}
							}
							else flg2=0;
							if(flg2) break; //forb. pth exits
						} 
						//cout<<'\n'<<"flg2_v: "<<flg2; 
						//cout<<'\n'<<"forb. pth: "; for(i2=0;i2<rev;i2++) cout<<" "<<frm_r[i2]<<","<<frm_c[i2]<<" to "<<to_r[i2]<<","<<to_c[i2]; cout<<'\n';
						if(!flg2) { //if no fbd path 
							if(l>1 || (par && l>=1)) {for(i2=0;i2<cyc_len;i2++) if(i==snk_r[i2] && c==snk_c[i2]) {flg=0; break;} else flg=1; }//snake can't hit itself
							//*if(chld_num==3)*/ cout<<'\n'<<"flgv: "<<flg;
							if((flg && l<cyc_len-1) || (flg && l==cyc_len-1 && i==stp_r)) {
								l++; snk_r[l-1-of]=i; snk_c[l-1-of]=c; ver=0; hrz=1; flg=0; break;
							}
						}
					}
				}
			}
			if(par && l==strt_L) flg3=1; else flg3=0; 
			if(!flg3 && ver) bck=1; else bck=0; //if search is allowed but still in ver mode, snake has to retreat
			//cout<<'\n'<<"bck: "<<bck;
			if(bck) {
				bck=0;
				for(i=0;i<cyc_len;i++) {fbd_r[rev][i]=snk_r[i]; fbd_c[rev][i]=snk_c[i];}  //forbidden path (the entire snake trajectory is a single fb. pth.)
				snk_r[l-1-of]=snk_c[l-1-of]=-1; //erasing retreated point
				l--; rev++; ver=0; hrz=1; 
			}
			//*if(chld_num>=3){*/cout<<'\n'<<"L_ver:"<<l; cout<<" "<<snk_r[l-2-of]<<","<<snk_c[l-2-of]<<" to "<<snk_r[l-1-of]<<","<<snk_c[l-1-of]; //}
			//if(chld_num>=3) {cout<<'\n'<<"fbd_pth: "<<'\n'; for(i2=0;i2<rev;i2++){for(i3=0;i3<cyc_len;i3++) cout<<" "<<fbd_r[i2][i3]<<","<<fbd_c[i2][i3]; cout<<'\n';}}
			//cout<<'\n'<<"snake_v: "; for(i2=0;i2<cyc_len;i2++) cout<<" "<<snk_r[i2]<<","<<snk_c[i2]<<"  ";
			//if(while_tot==5) break; //testing
		}	
		//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if(l==cyc_len) {
			if(!par) {for(i=0;i<cyc_len;i++) {par_r[i]=snk_r[i]; par_c[i]=snk_c[i];} par=1; par_num++; } //saving cycles
			else {for(i=0;i<cyc_len;i++) {chld_r[chld_num][i]=snk_r[i]; chld_c[chld_num][i]=snk_c[i];} chld_num++; }
			for(i=0;i<cyc_len;i++) {cyc_sv_r[ccnt][i]=snk_r[i]; cyc_sv_c[ccnt][i]=snk_c[i]; if(!(i%2)) cycol[ccnt][i/2]=snk_c[i];}
			ccnt++;
			//cout<<'\n'<<"ccnt: "<<ccnt;
			//cout<<'\n'<<"l:"<<l; cout<<" "<<snk_r[l-2-of]<<","<<snk_c[l-2-of]<<" to "<<snk_r[l-1-of]<<","<<snk_c[l-1-of]<<'\n';
			//cout<<'\n'<<"cycle!!!!!!!!!!!!!!!: "; for(i=0;i<l;i++) cout<<snk_r[i]<<","<<snk_c[i]<<" "; cout<<'\n';
			for(i=0;i<l;i++) {fbd_r[rev][i]=snk_r[i]; fbd_c[rev][i]=snk_c[i];} rev++; //forbidden path
			flg3=1; dnodes_chk=1; //need to find siblings
			//if(chld_num==5) break; //testing						
		}
		//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if(flg3) { 
			flg3=0;
			if(l<cyc_len) { //if no child is found, need to change 'of' or 'dnodes' 
				//dnodes_chk+=1;
				//if(dnodes_chk>dnodes) {
					flg4=1; 
					for(i=0;i<rev;i++) for(i2=0;i2<cyc_len;i2++) fbd_r[i][i2]=fbd_c[i][i2]=-1; rev=0; //refresh forbidden path
					dnodes_chk=1; of++; //'of' determines the starting variable node (eg. in cycle8, of=0 means 1st dnode is the 8th one)
					if(of==cyc_len-dnodes) {of=0; dnodes++;} //need to increase dnodes (no. of nodes that are different from parent)
				//}
				//cout<<'\n'<<"rev: "<<rev;
			}
			//cout<<'\n'<<"of: "<<of<<" dnodes: "<<dnodes<<" dnodes_chk: "<<dnodes_chk<<" flg4: "<<flg4;
			if(dnodes==cyc_len) break; //no more cycles can be found

			if(flg4) { //if change in 'of' or 'dnodes'
				for(i=0;i<cyc_len;i++) {fbd_r[rev][i]=par_r[i]; fbd_c[rev][i]=par_c[i];} //forbidden path
				rev++;
				for(i2=0;i2<chld_num;i2++) { //forb. paths from prev. children
					for(i=0;i<cyc_len;i++) {fbd_r[i2+1][i]=chld_r[i2][i]; fbd_c[i2+1][i]=chld_c[i2][i];} //forbidden path
				}
				rev+=chld_num;
			}

			K=dnodes;
			if(flg4 || !chld_num) { //if change in 'of' or 'dnodes'
				strt_r=par_r[cyc_len-of-K-1]; strt_c=par_c[cyc_len-of-K-1]; //start node is updated
				stp_r=par_r[cyc_len-of-1]; stp_c=par_c[cyc_len-of-1]; //stp node is updated
				for(i2=0;i2<cyc_len;i2++) {snk_r[i2]=par_r[i2]; snk_c[i2]=par_c[i2];} //refreshing the snake
			}
			else {
				strt_r=chld_r[chld_num-1][cyc_len-of-K-1]; strt_c=chld_c[chld_num-1][cyc_len-of-K-1]; //start node is updated
				stp_r=chld_r[chld_num-1][cyc_len-of-1]; stp_c=chld_c[chld_num-1][cyc_len-of-1]; //stop node is updated (stp_c needed in hrz mode, stp_r needed in ver mode)
			}
			for(i=0;i<dnodes;i++) {snk_r[cyc_len-of-1-i]=snk_c[cyc_len-of-1-i]=-1;}
			hrz=!((cyc_len-of-K-1)%2); ver=!hrz; //understanding from strt node index
			l=cyc_len-K; strt_L=l; //snake length at beginning of search
			
			flg4=0;
			//cout<<'\n'<<"l: "<<l<<" hrz:"<<hrz<<" ver: "<<ver<<" strt: "<<strt_r<<","<<strt_c<<" stp_r: "<<stp_r<<" stp_c: "<<stp_c;
			//if(while_tot>=0) {cout<<'\n'<<"fbd_pth: "<<'\n'; for(i2=0;i2<rev;i2++){for(i3=0;i3<cyc_len;i3++) cout<<" "<<fbd_r[i2][i3]<<","<<fbd_c[i2][i3]; cout<<'\n';}}

		}	
		//while_tot++; //cout<<'\n'<<"while_tot: "<<while_tot;
		//rev_tot+=rev;	
	}		
}
	

