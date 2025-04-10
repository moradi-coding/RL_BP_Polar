
//BCJR calclations using the [15,11] Hamming trellis 
void bcjr_log_15_11(int a) { 
	long i,j;
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,
	tmp1b,tmp2b,tmp3b,tmp4b,tmp5b,tmp6b,tmp7b,tmp8b,tmp9b,tmp10b,tmp11b,tmp12b,tmp13b,tmp14b,tmp15b,tmp16b,
	max1,max2,max3,max1b,max2b,max3b,c,y_sum,b,d,num,s000_6_alphalph_old;
	
	double gama1_0, gama1_1, gama2_0, gama2_1, gama3_0, gama3_1, gama4_0, gama4_1, gama5_0, gama5_1, gama6_0, gama6_1, gama7_0, gama7_1, gama8_0, gama8_1, gama9_0, gama9_1, gama10_0, gama10_1, gama11_0, gama11_1, gama12_0, gama12_1, gama13_0, gama13_1, gama14_0, gama14_1, gama15_0, gama15_1; 
	
	double alph, alph1_0, alph1_1, alph2_00, alph2_01, alph2_10, alph2_11, alph3_000, alph3_001, alph3_010, alph3_011, alph3_100, alph3_101, alph3_110, alph3_111,
	alph4_0000, alph4_0001, alph4_0010, alph4_0011, alph4_0100, alph4_0101, alph4_0110, alph4_0111, alph4_1000, alph4_1001, alph4_1010, alph4_1011, alph4_1100, alph4_1101, alph4_1110, alph4_1111,
	alph5_0000, alph5_0001, alph5_0010, alph5_0011, alph5_0100, alph5_0101, alph5_0110, alph5_0111, alph5_1000, alph5_1001, alph5_1010, alph5_1011, alph5_1100, alph5_1101, alph5_1110, alph5_1111,
	alph6_0000, alph6_0001, alph6_0010, alph6_0011, alph6_0100, alph6_0101, alph6_0110, alph6_0111, alph6_1000, alph6_1001, alph6_1010, alph6_1011, alph6_1100, alph6_1101, alph6_1110, alph6_1111,
	alph7_0000, alph7_0001, alph7_0010, alph7_0011, alph7_0100, alph7_0101, alph7_0110, alph7_0111, alph7_1000, alph7_1001, alph7_1010, alph7_1011, alph7_1100, alph7_1101, alph7_1110, alph7_1111,
	alph8_0000, alph8_0001, alph8_0010, alph8_0011, alph8_0100, alph8_0101, alph8_0110, alph8_0111, alph8_1000, alph8_1001, alph8_1010, alph8_1011, alph8_1100, alph8_1101, alph8_1110, alph8_1111,
	alph9_0000, alph9_0001, alph9_0010, alph9_0011, alph9_0100, alph9_0101, alph9_0110, alph9_0111, alph9_1000, alph9_1001, alph9_1010, alph9_1011, alph9_1100, alph9_1101, alph9_1110, alph9_1111,
	alph10_0000, alph10_0001, alph10_0010, alph10_0011, alph10_0100, alph10_0101, alph10_0110, alph10_0111, alph10_1000, alph10_1001, alph10_1010, alph10_1011, alph10_1100, alph10_1101, alph10_1110, alph10_1111,
	alph11_0000, alph11_0001, alph11_0010, alph11_0011, alph11_0100, alph11_0101, alph11_0110, alph11_0111, alph11_1000, alph11_1001, alph11_1010, alph11_1011, alph11_1100, alph11_1101, alph11_1110, alph11_1111,
	alph12_000, alph12_001, alph12_010, alph12_011, alph12_100, alph12_101, alph12_110, alph12_111,
	alph13_00, alph13_01, alph13_10, alph13_11, alph14_0, alph14_1;

	double beta15_0, beta1_0, beta1_1, beta2_00, beta2_01, beta2_10, beta2_11, beta3_000, beta3_001, beta3_010, beta3_011, beta3_100, beta3_101, beta3_110, beta3_111,
	beta4_0000, beta4_0001, beta4_0010, beta4_0011, beta4_0100, beta4_0101, beta4_0110, beta4_0111, beta4_1000, beta4_1001, beta4_1010, beta4_1011, beta4_1100, beta4_1101, beta4_1110, beta4_1111,
	beta5_0000, beta5_0001, beta5_0010, beta5_0011, beta5_0100, beta5_0101, beta5_0110, beta5_0111, beta5_1000, beta5_1001, beta5_1010, beta5_1011, beta5_1100, beta5_1101, beta5_1110, beta5_1111,
	beta6_0000, beta6_0001, beta6_0010, beta6_0011, beta6_0100, beta6_0101, beta6_0110, beta6_0111, beta6_1000, beta6_1001, beta6_1010, beta6_1011, beta6_1100, beta6_1101, beta6_1110, beta6_1111,
	beta7_0000, beta7_0001, beta7_0010, beta7_0011, beta7_0100, beta7_0101, beta7_0110, beta7_0111, beta7_1000, beta7_1001, beta7_1010, beta7_1011, beta7_1100, beta7_1101, beta7_1110, beta7_1111,
	beta8_0000, beta8_0001, beta8_0010, beta8_0011, beta8_0100, beta8_0101, beta8_0110, beta8_0111, beta8_1000, beta8_1001, beta8_1010, beta8_1011, beta8_1100, beta8_1101, beta8_1110, beta8_1111,
	beta9_0000, beta9_0001, beta9_0010, beta9_0011, beta9_0100, beta9_0101, beta9_0110, beta9_0111, beta9_1000, beta9_1001, beta9_1010, beta9_1011, beta9_1100, beta9_1101, beta9_1110, beta9_1111,
	beta10_0000, beta10_0001, beta10_0010, beta10_0011, beta10_0100, beta10_0101, beta10_0110, beta10_0111, beta10_1000, beta10_1001, beta10_1010, beta10_1011, beta10_1100, beta10_1101, beta10_1110, beta10_1111,
	beta11_0000, beta11_0001, beta11_0010, beta11_0011, beta11_0100, beta11_0101, beta11_0110, beta11_0111, beta11_1000, beta11_1001, beta11_1010, beta11_1011, beta11_1100, beta11_1101, beta11_1110, beta11_1111,
	beta12_000, beta12_001, beta12_010, beta12_011, beta12_100, beta12_101, beta12_110, beta12_111,
	beta13_00, beta13_01, beta13_10, beta13_11, beta14_0, beta14_1;
	
	d=0.5; 
	
	c=-1/(2*varn);
	//for(a=0;a<num_cls;a++) {
		//initialization of branch metrics
		for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
			i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'	
			//L_in[i]=L_in[i];
			b=L_in[i]-LR[i]; //gives the extrinsic LLR
			//cout<<"j, b, LLR: "<<j<<", "<<b<<", "<<L_in[i]<<endl;	
			
			if(!j) {
				gama1_0=th(c*pow(abs(y[i]-1),2)+d*1*b); //0 mapped to 1	
				gama1_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); //1 mapped to -1	
			}	
			else if(j==1) {
				gama2_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama2_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==2) {
				gama3_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama3_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==3) {
				gama4_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama4_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==4) {
				gama5_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama5_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==5) {
				gama6_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama6_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==6) {
				gama7_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama7_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==7) {
				gama8_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama8_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==8) {
				gama9_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama9_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==9) {
				gama10_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama10_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==10) {
				gama11_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama11_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==11) {
				gama12_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama12_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==12) {
				gama13_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama13_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==13) {
				gama14_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama14_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==14) {
				gama15_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama15_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
		}
		//cout<<'\n'<<"L_in: "; for(i=0;i<n;i++) cout<<L_in[i]<< " "; cout<<'\n';
		
		//forward metrics (alpha)
		alph=log(1);
		alph1_0=th(alph+gama1_0); //threshold the result
		alph1_1=th(alph+gama1_1); 
		
		alph2_00=th(alph1_0+gama2_0);
		alph2_01=th(alph1_0+gama2_1);
		alph2_10=th(alph1_1+gama2_1);
		alph2_11=th(alph1_1+gama2_0);
		
		alph3_000=th(alph2_00+gama3_0);
		alph3_001=th(alph2_00+gama3_1);
		alph3_010=th(alph2_01+gama3_1);
		alph3_011=th(alph2_01+gama3_0);
		alph3_100=th(alph2_10+gama3_1);
		alph3_101=th(alph2_10+gama3_0);
		alph3_110=th(alph2_11+gama3_0);
		alph3_111=th(alph2_11+gama3_1);
		
		alph4_0000=th(alph3_000+gama4_0);
		alph4_0001=th(alph3_000+gama4_1);
		alph4_0010=th(alph3_001+gama4_0);
		alph4_0011=th(alph3_001+gama4_1);
		alph4_0100=th(alph3_010+gama4_1);
		alph4_0101=th(alph3_010+gama4_0);
		alph4_0110=th(alph3_011+gama4_1);
		alph4_0111=th(alph3_011+gama4_0);
		alph4_1000=th(alph3_100+gama4_1);
		alph4_1001=th(alph3_100+gama4_0);
		alph4_1010=th(alph3_101+gama4_1);
		alph4_1011=th(alph3_101+gama4_0);
		alph4_1100=th(alph3_110+gama4_0);
		alph4_1101=th(alph3_110+gama4_1);
		alph4_1110=th(alph3_111+gama4_0);
		alph4_1111=th(alph3_111+gama4_1);
		
		//cout<<"s001_2_alpha: "<<s001_2_alpha<<endl;
		alph5_0000= th2(alph4_0000+gama5_0,alph4_1000+gama5_1);
		alph5_0001= th2(alph4_0000+gama5_1,alph4_1000+gama5_0);
		alph5_0010= th2(alph4_0001+gama5_1,alph4_1001+gama5_0);
		alph5_0011= th2(alph4_0001+gama5_0,alph4_1001+gama5_1);
		alph5_0100= th2(alph4_0010+gama5_1,alph4_1010+gama5_0);
		alph5_0101= th2(alph4_0010+gama5_0,alph4_1010+gama5_1);
		alph5_0110= th2(alph4_0011+gama5_0,alph4_1011+gama5_1);
		alph5_0111= th2(alph4_0011+gama5_1,alph4_1011+gama5_0);
		alph5_1000= th2(alph4_0100+gama5_1,alph4_1100+gama5_0);
		alph5_1001= th2(alph4_0100+gama5_0,alph4_1100+gama5_1);
		alph5_1010= th2(alph4_0101+gama5_0,alph4_1101+gama5_1);
		alph5_1011= th2(alph4_0101+gama5_1,alph4_1101+gama5_0);
		alph5_1100= th2(alph4_0110+gama5_0,alph4_1110+gama5_1);
		alph5_1101= th2(alph4_0110+gama5_1,alph4_1110+gama5_0);
		alph5_1110= th2(alph4_0111+gama5_1,alph4_1111+gama5_0);
		alph5_1111= th2(alph4_0111+gama5_0,alph4_1111+gama5_1);
		
		alph6_0000= th2(alph5_0000+gama6_0,alph5_0010+gama6_1);
		alph6_0001= th2(alph5_0000+gama6_1,alph5_0010+gama6_0);
		alph6_0010= th2(alph5_0001+gama6_0,alph5_0011+gama6_1);
		alph6_0011= th2(alph5_0001+gama6_1,alph5_0011+gama6_0);
		alph6_0100= th2(alph5_0100+gama6_0,alph5_0110+gama6_1);
		alph6_0101= th2(alph5_0100+gama6_1,alph5_0110+gama6_0);
		alph6_0110= th2(alph5_0101+gama6_0,alph5_0111+gama6_1);
		alph6_0111= th2(alph5_0101+gama6_1,alph5_0111+gama6_0);
		alph6_1000= th2(alph5_1000+gama6_1,alph5_1010+gama6_0);
		alph6_1001= th2(alph5_1000+gama6_0,alph5_1010+gama6_1);
		alph6_1010= th2(alph5_1001+gama6_1,alph5_1011+gama6_0);
		alph6_1011= th2(alph5_1001+gama6_0,alph5_1011+gama6_1);
		alph6_1100= th2(alph5_1100+gama6_1,alph5_1110+gama6_0);
		alph6_1101= th2(alph5_1100+gama6_0,alph5_1110+gama6_1);
		alph6_1110= th2(alph5_1101+gama6_1,alph5_1111+gama6_0);
		alph6_1111= th2(alph5_1101+gama6_0,alph5_1111+gama6_1);
		
		alph7_0000= th2(alph6_0000+gama7_0,alph6_0100+gama7_1);
		alph7_0001= th2(alph6_0000+gama7_1,alph6_0100+gama7_0);
		alph7_0010= th2(alph6_0001+gama7_1,alph6_0101+gama7_0);
		alph7_0011= th2(alph6_0001+gama7_0,alph6_0101+gama7_1);
		alph7_0100= th2(alph6_0010+gama7_0,alph6_0110+gama7_1);
		alph7_0101= th2(alph6_0010+gama7_1,alph6_0110+gama7_0);
		alph7_0110= th2(alph6_0011+gama7_1,alph6_0111+gama7_0);
		alph7_0111= th2(alph6_0011+gama7_0,alph6_0111+gama7_1);
		alph7_1000= th2(alph6_1000+gama7_0,alph6_1100+gama7_1);
		alph7_1001= th2(alph6_1000+gama7_1,alph6_1100+gama7_0);
		alph7_1010= th2(alph6_1001+gama7_1,alph6_1101+gama7_0);
		alph7_1011= th2(alph6_1001+gama7_0,alph6_1101+gama7_1);
		alph7_1100= th2(alph6_1010+gama7_0,alph6_1110+gama7_1);
		alph7_1101= th2(alph6_1010+gama7_1,alph6_1110+gama7_0);
		alph7_1110= th2(alph6_1011+gama7_1,alph6_1111+gama7_0);
		alph7_1111= th2(alph6_1011+gama7_0,alph6_1111+gama7_1);
		
		alph8_0000= th2(alph7_0000+gama8_0,alph7_0010+gama8_1);
		alph8_0001= th2(alph7_0000+gama8_1,alph7_0010+gama8_0);
		alph8_0010= th2(alph7_0001+gama8_1,alph7_0011+gama8_0);
		alph8_0011= th2(alph7_0001+gama8_0,alph7_0011+gama8_1);
		alph8_0100= th2(alph7_0100+gama8_1,alph7_0110+gama8_0);
		alph8_0101= th2(alph7_0100+gama8_0,alph7_0110+gama8_1);
		alph8_0110= th2(alph7_0101+gama8_0,alph7_0111+gama8_1);
		alph8_0111= th2(alph7_0101+gama8_1,alph7_0111+gama8_0);
		alph8_1000= th2(alph7_1000+gama8_1,alph7_1010+gama8_0);
		alph8_1001= th2(alph7_1000+gama8_0,alph7_1010+gama8_1);
		alph8_1010= th2(alph7_1001+gama8_0,alph7_1011+gama8_1);
		alph8_1011= th2(alph7_1001+gama8_1,alph7_1011+gama8_0);
		alph8_1100= th2(alph7_1100+gama8_0,alph7_1110+gama8_1);
		alph8_1101= th2(alph7_1100+gama8_1,alph7_1110+gama8_0);
		alph8_1110= th2(alph7_1101+gama8_1,alph7_1111+gama8_0);
		alph8_1111= th2(alph7_1101+gama8_0,alph7_1111+gama8_1);
		
		alph9_0000= th2(alph8_0000+gama9_0,alph8_0010+gama9_1);
		alph9_0001= th2(alph8_0000+gama9_1,alph8_0010+gama9_0);
		alph9_0010= th2(alph8_0001+gama9_1,alph8_0011+gama9_0);
		alph9_0011= th2(alph8_0001+gama9_0,alph8_0011+gama9_1);
		alph9_0100= th2(alph8_0100+gama9_0,alph8_0110+gama9_1);
		alph9_0101= th2(alph8_0100+gama9_1,alph8_0110+gama9_0);
		alph9_0110= th2(alph8_0101+gama9_1,alph8_0111+gama9_0);
		alph9_0111= th2(alph8_0101+gama9_0,alph8_0111+gama9_1);
		alph9_1000= th2(alph8_1000+gama9_1,alph8_1010+gama9_0);
		alph9_1001= th2(alph8_1000+gama9_0,alph8_1010+gama9_1);
		alph9_1010= th2(alph8_1001+gama9_0,alph8_1011+gama9_1);
		alph9_1011= th2(alph8_1001+gama9_1,alph8_1011+gama9_0);
		alph9_1100= th2(alph8_1100+gama9_1,alph8_1110+gama9_0);
		alph9_1101= th2(alph8_1100+gama9_0,alph8_1110+gama9_1);
		alph9_1110= th2(alph8_1101+gama9_0,alph8_1111+gama9_1);
		alph9_1111= th2(alph8_1101+gama9_1,alph8_1111+gama9_0);
		
		alph10_0000= th2(alph9_0000+gama10_0,alph9_0100+gama10_1);
		alph10_0001= th2(alph9_0000+gama10_1,alph9_0100+gama10_0);
		alph10_0010= th2(alph9_0001+gama10_0,alph9_0101+gama10_1);
		alph10_0011= th2(alph9_0001+gama10_1,alph9_0101+gama10_0);
		alph10_0100= th2(alph9_0010+gama10_0,alph9_0110+gama10_1);
		alph10_0101= th2(alph9_0010+gama10_1,alph9_0110+gama10_0);
		alph10_0110= th2(alph9_0011+gama10_0,alph9_0111+gama10_1);
		alph10_0111= th2(alph9_0011+gama10_1,alph9_0111+gama10_0);
		alph10_1000= th2(alph9_1000+gama10_0,alph9_1100+gama10_1);
		alph10_1001= th2(alph9_1000+gama10_1,alph9_1100+gama10_0);
		alph10_1010= th2(alph9_1001+gama10_0,alph9_1101+gama10_1);
		alph10_1011= th2(alph9_1001+gama10_1,alph9_1101+gama10_0);
		alph10_1100= th2(alph9_1010+gama10_0,alph9_1110+gama10_1);
		alph10_1101= th2(alph9_1010+gama10_1,alph9_1110+gama10_0);
		alph10_1110= th2(alph9_1011+gama10_0,alph9_1111+gama10_1);
		alph10_1111= th2(alph9_1011+gama10_1,alph9_1111+gama10_0);
		
		alph11_0000= th2(alph10_0000+gama11_0,alph10_0100+gama11_1);
		alph11_0001= th2(alph10_0000+gama11_1,alph10_0100+gama11_0);
		alph11_0010= th2(alph10_0001+gama11_0,alph10_0101+gama11_1);
		alph11_0011= th2(alph10_0001+gama11_1,alph10_0101+gama11_0);
		alph11_0100= th2(alph10_0010+gama11_0,alph10_0110+gama11_1);
		alph11_0101= th2(alph10_0010+gama11_1,alph10_0110+gama11_0);
		alph11_0110= th2(alph10_0011+gama11_0,alph10_0111+gama11_1);
		alph11_0111= th2(alph10_0011+gama11_1,alph10_0111+gama11_0);
		alph11_1000= th2(alph10_1000+gama11_0,alph10_1100+gama11_1);
		alph11_1001= th2(alph10_1000+gama11_1,alph10_1100+gama11_0);
		alph11_1010= th2(alph10_1001+gama11_0,alph10_1101+gama11_1);
		alph11_1011= th2(alph10_1001+gama11_1,alph10_1101+gama11_0);
		alph11_1100= th2(alph10_1010+gama11_0,alph10_1110+gama11_1);
		alph11_1101= th2(alph10_1010+gama11_1,alph10_1110+gama11_0);
		alph11_1110= th2(alph10_1011+gama11_0,alph10_1111+gama11_1);
		alph11_1111= th2(alph10_1011+gama11_1,alph10_1111+gama11_0);
		
		alph12_000= th2(alph11_0000+gama12_0,alph11_1000+gama12_1);
		alph12_001= th2(alph11_0001+gama12_1,alph11_1001+gama12_0);
		alph12_010= th2(alph11_0010+gama12_1,alph11_1010+gama12_0);
		alph12_011= th2(alph11_0011+gama12_0,alph11_1011+gama12_1);
		alph12_100= th2(alph11_0100+gama12_1,alph11_1100+gama12_0);
		alph12_101= th2(alph11_0101+gama12_0,alph11_1101+gama12_1);
		alph12_110= th2(alph11_0110+gama12_0,alph11_1110+gama12_1);
		alph12_111= th2(alph11_0111+gama12_1,alph11_1111+gama12_0);
		
		alph13_00= th2(alph12_000+gama13_0,alph12_001+gama13_1);
		alph13_01= th2(alph12_010+gama13_0,alph12_011+gama13_1);
		alph13_10= th2(alph12_100+gama13_0,alph12_101+gama13_1);
		alph13_11= th2(alph12_110+gama13_0,alph12_111+gama13_1);
		
		alph14_0= th2(alph13_00+gama13_0,alph13_10+gama13_1);
		alph14_1= th2(alph13_01+gama13_0,alph13_11+gama13_1);
		
		//alph15_0= th2(alph14_0+gama15_0,alph14_1+gama15_1);
				
		//cout<<"here";
		/*if(abs(tmp1-tmp2)>100) {	
			//cout<<"s000_6_alphalph_old: "<<s000_6_alphalph_old<<endl;	
			cout<<"num: "<<abs(tmp1-tmp2)<<endl;	
			cout<<" num2: "<<exp(-1*abs(tmp1-tmp2))<<endl;	
			cout<<"s000_6_alpha: "<<s000_6_alpha<<endl;	
		}*/	
		
		/*cout<<"alph12_000: "<<alph12_000<<endl;
		cout<<"alph12_001: "<<alph12_001<<endl;
		cout<<"alph12_010: "<<alph12_010<<endl;
		cout<<"alph12_011: "<<alph12_011<<endl;
		cout<<"alph12_100: "<<alph12_100<<endl;
		cout<<"alph12_101: "<<alph12_101<<endl;
		cout<<"alph12_110: "<<alph12_110<<endl;
		cout<<"alph12_111: "<<alph12_111<<endl;
		cout<<endl;
		cout<<"alph13_00: "<<alph13_00<<endl;
		cout<<"alph13_01: "<<alph13_01<<endl;
		cout<<"alph13_10: "<<alph13_10<<endl;
		cout<<"alph13_11: "<<alph13_11<<endl;
		cout<<endl;
		cout<<"alph14_0: "<<alph14_0<<endl;
		cout<<"alph14_1: "<<alph14_1<<endl;
		cout<<endl;*/
	
	
		//backward metrics
		beta15_0=log(1);
		
		beta14_0= th(beta15_0+gama15_0);
		beta14_1= th(beta15_0+gama15_1);
		
		beta13_00= th(beta14_0+gama14_0);
		beta13_01= th(beta14_1+gama14_0);
		beta13_10= th(beta14_0+gama14_1);
		beta13_11= th(beta14_1+gama14_1);
		
		beta12_000= th(beta13_00+gama13_0);
		beta12_001= th(beta13_00+gama13_1);
		beta12_010= th(beta13_01+gama13_0);
		beta12_011= th(beta13_01+gama13_1);
		beta12_100= th(beta13_10+gama13_0);
		beta12_101= th(beta13_10+gama13_1);
		beta12_110= th(beta13_11+gama13_0);
		beta12_111= th(beta13_11+gama13_1);
		
		beta11_0000= th(beta12_000+gama12_0);
		beta11_0001= th(beta12_001+gama12_1);
		beta11_0010= th(beta12_010+gama12_1);
		beta11_0011= th(beta12_011+gama12_0);
		beta11_0100= th(beta12_100+gama12_1);
		beta11_0101= th(beta12_101+gama12_0);
		beta11_0110= th(beta12_110+gama12_0);
		beta11_0111= th(beta12_111+gama12_1);
		beta11_1000= th(beta12_000+gama12_1);
		beta11_1001= th(beta12_001+gama12_0);
		beta11_1010= th(beta12_010+gama12_0);
		beta11_1011= th(beta12_011+gama12_1);
		beta11_1100= th(beta12_100+gama12_0);
		beta11_1101= th(beta12_101+gama12_1);
		beta11_1110= th(beta12_110+gama12_1);
		beta11_1111= th(beta12_111+gama12_0);
		
		beta10_0000= th2(beta11_0000+gama11_0,beta11_0001+gama11_1);
		beta10_0001= th2(beta11_0010+gama11_0,beta11_0011+gama11_1);
		beta10_0010= th2(beta11_0100+gama11_0,beta11_0101+gama11_1);
		beta10_0011= th2(beta11_0110+gama11_0,beta11_0111+gama11_1);
		beta10_0100= th2(beta11_0000+gama11_1,beta11_0001+gama11_0);
		beta10_0101= th2(beta11_0010+gama11_1,beta11_0011+gama11_0);
		beta10_0110= th2(beta11_0100+gama11_1,beta11_0101+gama11_0);
		beta10_0111= th2(beta11_0110+gama11_1,beta11_0111+gama11_0);
		beta10_1000= th2(beta11_1000+gama11_0,beta11_1001+gama11_1);
		beta10_1001= th2(beta11_1010+gama11_0,beta11_1011+gama11_1);
		beta10_1010= th2(beta11_1100+gama11_0,beta11_1101+gama11_1);
		beta10_1011= th2(beta11_1110+gama11_0,beta11_1111+gama11_1);
		beta10_1100= th2(beta11_1000+gama11_1,beta11_1001+gama11_0);
		beta10_1101= th2(beta11_1010+gama11_1,beta11_1011+gama11_0);
		beta10_1110= th2(beta11_1100+gama11_1,beta11_1101+gama11_0);
		beta10_1111= th2(beta11_1110+gama11_1,beta11_1111+gama11_0);
		
		beta9_0000= th2(beta10_0000+gama10_0,beta10_0001+gama10_1);
		beta9_0001= th2(beta10_0010+gama10_0,beta10_0011+gama10_1);
		beta9_0010= th2(beta10_0100+gama10_0,beta10_0101+gama10_1);
		beta9_0011= th2(beta10_0110+gama10_0,beta10_0111+gama10_1);
		beta9_0100= th2(beta10_0000+gama10_1,beta10_0001+gama10_0);
		beta9_0101= th2(beta10_0010+gama10_1,beta10_0011+gama10_0);
		beta9_0110= th2(beta10_0100+gama10_1,beta10_0101+gama10_0);
		beta9_0111= th2(beta10_0110+gama10_1,beta10_0111+gama10_0);
		beta9_1000= th2(beta10_1000+gama10_0,beta10_1001+gama10_1);
		beta9_1001= th2(beta10_1010+gama10_0,beta10_1011+gama10_1);
		beta9_1010= th2(beta10_1100+gama10_0,beta10_1101+gama10_1);
		beta9_1011= th2(beta10_1110+gama10_0,beta10_1111+gama10_1);
		beta9_1100= th2(beta10_1000+gama10_1,beta10_1001+gama10_0);
		beta9_1101= th2(beta10_1010+gama10_1,beta10_1011+gama10_0);
		beta9_1110= th2(beta10_1100+gama10_1,beta10_1101+gama10_0);
		beta9_1111= th2(beta10_1110+gama10_1,beta10_1111+gama10_0);
		
		beta8_0000= th2(beta9_0000+gama9_0,beta9_0001+gama9_1);
		beta8_0001= th2(beta9_0010+gama9_1,beta9_0011+gama9_0);
		beta8_0010= th2(beta9_0000+gama9_1,beta9_0001+gama9_0);
		beta8_0011= th2(beta9_0010+gama9_0,beta9_0011+gama9_1);
		beta8_0100= th2(beta9_0100+gama9_0,beta9_0101+gama9_1);
		beta8_0101= th2(beta9_0110+gama9_1,beta9_0111+gama9_0);
		beta8_0110= th2(beta9_0100+gama9_1,beta9_0101+gama9_0);
		beta8_0111= th2(beta9_0110+gama9_0,beta9_0111+gama9_1);
		beta8_1000= th2(beta9_1000+gama9_1,beta9_1001+gama9_0);
		beta8_1001= th2(beta9_1010+gama9_0,beta9_1011+gama9_1);
		beta8_1010= th2(beta9_1000+gama9_0,beta9_1001+gama9_1);
		beta8_1011= th2(beta9_1010+gama9_1,beta9_1011+gama9_0);
		beta8_1100= th2(beta9_1100+gama9_1,beta9_1101+gama9_0);
		beta8_1101= th2(beta9_1110+gama9_0,beta9_1111+gama9_1);
		beta8_1110= th2(beta9_1100+gama9_0,beta9_1101+gama9_1);
		beta8_1111= th2(beta9_1110+gama9_1,beta9_1111+gama9_0);
		
		beta7_0000= th2(beta8_0000+gama8_0,beta8_0001+gama8_1);
		beta7_0001= th2(beta8_0010+gama8_1,beta8_0011+gama8_0);
		beta7_0010= th2(beta8_0000+gama8_1,beta8_0001+gama8_0);
		beta7_0011= th2(beta8_0010+gama8_0,beta8_0011+gama8_1);
		beta7_0100= th2(beta8_0100+gama8_1,beta8_0101+gama8_0);
		beta7_0101= th2(beta8_0110+gama8_0,beta8_0111+gama8_1);
		beta7_0110= th2(beta8_0100+gama8_0,beta8_0101+gama8_1);
		beta7_0111= th2(beta8_0110+gama8_1,beta8_0111+gama8_0);
		beta7_1000= th2(beta8_1000+gama8_1,beta8_1001+gama8_0);
		beta7_1001= th2(beta8_1010+gama8_0,beta8_1011+gama8_1);
		beta7_1010= th2(beta8_1000+gama8_0,beta8_1001+gama8_1);
		beta7_1011= th2(beta8_1010+gama8_1,beta8_1011+gama8_0);
		beta7_1100= th2(beta8_1100+gama8_0,beta8_1101+gama8_1);
		beta7_1101= th2(beta8_1110+gama8_1,beta8_1111+gama8_0);
		beta7_1110= th2(beta8_1100+gama8_1,beta8_1101+gama8_0);
		beta7_1111= th2(beta8_1110+gama8_0,beta8_1111+gama8_1);
		
		beta6_0000= th2(beta7_0000+gama7_0,beta7_0001+gama7_1);
		beta6_0001= th2(beta7_0010+gama7_1,beta7_0011+gama7_0);
		beta6_0010= th2(beta7_0100+gama7_0,beta7_0101+gama7_1);
		beta6_0011= th2(beta7_0110+gama7_1,beta7_0111+gama7_0);
		beta6_0100= th2(beta7_0000+gama7_1,beta7_0001+gama7_0);
		beta6_0101= th2(beta7_0010+gama7_0,beta7_0011+gama7_1);
		beta6_0110= th2(beta7_0100+gama7_1,beta7_0101+gama7_0);
		beta6_0111= th2(beta7_0110+gama7_0,beta7_0111+gama7_1);
		beta6_1000= th2(beta7_1000+gama7_0,beta7_1001+gama7_1);
		beta6_1001= th2(beta7_1010+gama7_1,beta7_1011+gama7_0);
		beta6_1010= th2(beta7_1100+gama7_0,beta7_1101+gama7_1);
		beta6_1011= th2(beta7_1110+gama7_1,beta7_1111+gama7_0);
		beta6_1100= th2(beta7_1000+gama7_1,beta7_1001+gama7_0);
		beta6_1101= th2(beta7_1010+gama7_0,beta7_1011+gama7_1);
		beta6_1110= th2(beta7_1100+gama7_1,beta7_1101+gama7_0);
		beta6_1111= th2(beta7_1110+gama7_0,beta7_1111+gama7_1);
		
		beta5_0000= th2(beta6_0000+gama6_0,beta6_0001+gama6_1);
		beta5_0001= th2(beta6_0010+gama6_0,beta6_0011+gama6_1);
		beta5_0010= th2(beta6_0000+gama6_1,beta6_0001+gama6_0);
		beta5_0011= th2(beta6_0010+gama6_1,beta6_0011+gama6_0);
		beta5_0100= th2(beta6_0100+gama6_0,beta6_0101+gama6_1);
		beta5_0101= th2(beta6_0110+gama6_0,beta6_0111+gama6_1);
		beta5_0110= th2(beta6_0100+gama6_1,beta6_0101+gama6_0);
		beta5_0111= th2(beta6_0110+gama6_1,beta6_0111+gama6_0);
		beta5_1000= th2(beta6_1000+gama6_1,beta6_1001+gama6_0);
		beta5_1001= th2(beta6_1010+gama6_1,beta6_1011+gama6_0);
		beta5_1010= th2(beta6_1000+gama6_0,beta6_1001+gama6_1);
		beta5_1011= th2(beta6_1010+gama6_0,beta6_1011+gama6_1);
		beta5_1100= th2(beta6_1100+gama6_1,beta6_1101+gama6_0);
		beta5_1101= th2(beta6_1110+gama6_1,beta6_1111+gama6_0);
		beta5_1110= th2(beta6_1100+gama6_0,beta6_1101+gama6_1);
		beta5_1111= th2(beta6_1110+gama6_0,beta6_1111+gama6_1);
		
		beta4_0000= th2(beta5_0000+gama5_0,beta5_0001+gama5_1);
		beta4_0001= th2(beta5_0010+gama5_1,beta5_0011+gama5_0);
		beta4_0010= th2(beta5_0100+gama5_1,beta5_0101+gama5_0);
		beta4_0011= th2(beta5_0110+gama5_0,beta5_0111+gama5_1);
		beta4_0100= th2(beta5_1000+gama5_1,beta5_1001+gama5_0);
		beta4_0101= th2(beta5_1010+gama5_0,beta5_1011+gama5_1);
		beta4_0110= th2(beta5_1100+gama5_0,beta5_1101+gama5_1);
		beta4_0111= th2(beta5_1110+gama5_1,beta5_1111+gama5_0);
		beta4_1000= th2(beta5_0000+gama5_1,beta5_0001+gama5_0);
		beta4_1001= th2(beta5_0010+gama5_0,beta5_0011+gama5_1);
		beta4_1010= th2(beta5_0100+gama5_0,beta5_0101+gama5_1);
		beta4_1011= th2(beta5_0110+gama5_1,beta5_0111+gama5_0);
		beta4_1100= th2(beta5_1000+gama5_0,beta5_1001+gama5_1);
		beta4_1101= th2(beta5_1010+gama5_1,beta5_1011+gama5_0);
		beta4_1110= th2(beta5_1100+gama5_1,beta5_1101+gama5_0);
		beta4_1111= th2(beta5_1110+gama5_0,beta5_1111+gama5_1);
		
		beta3_000= th2(beta4_0000+gama4_0,beta4_0001+gama4_1);
		beta3_001= th2(beta4_0010+gama4_0,beta4_0011+gama4_1);
		beta3_010= th2(beta4_0100+gama4_1,beta4_0101+gama4_0);
		beta3_011= th2(beta4_0110+gama4_1,beta4_0111+gama4_0);
		beta3_100= th2(beta4_1000+gama4_1,beta4_1001+gama4_0);
		beta3_101= th2(beta4_1010+gama4_1,beta4_1011+gama4_0);
		beta3_110= th2(beta4_1100+gama4_0,beta4_1101+gama4_1);
		beta3_111= th2(beta4_1110+gama4_0,beta4_1111+gama4_1);
		
		beta2_00= th2(beta3_000+gama3_0,beta3_001+gama3_1);
		beta2_01= th2(beta3_010+gama3_1,beta3_011+gama3_0);
		beta2_10= th2(beta3_100+gama3_1,beta3_101+gama3_0);
		beta2_11= th2(beta3_110+gama3_0,beta3_111+gama3_1);
		
		beta1_0= th2(beta2_00+gama2_0,beta2_01+gama2_1);
		beta1_1= th2(beta2_10+gama2_1,beta2_11+gama2_0);
		
		/*cout<<"beta12_000: "<<beta12_000<<endl;
		cout<<"beta12_001: "<<beta12_001<<endl;
		cout<<"beta12_010: "<<beta12_010<<endl;
		cout<<"beta12_011: "<<beta12_011<<endl;
		cout<<"beta12_100: "<<beta12_100<<endl;
		cout<<"beta12_101: "<<beta12_101<<endl;
		cout<<"beta12_110: "<<beta12_110<<endl;
		cout<<"beta12_111: "<<beta12_111<<endl;
		cout<<endl;
		cout<<"beta13_00: "<<beta13_00<<endl;
		cout<<"beta13_01: "<<beta13_01<<endl;
		cout<<"beta13_10: "<<beta13_10<<endl;
		cout<<"beta13_11: "<<beta13_11<<endl;
		cout<<endl;
		cout<<"beta14_0: "<<beta14_0<<endl;
		cout<<"beta14_1: "<<beta14_1<<endl;
		cout<<endl;*/
		
		
	//computing LLRs at each stage (messages sent by this CN to its neighboring VNs)
	//cout<<"L_out: "; 
	for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
		i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'
		if(!j) {
			//L_out[i]= log(exp(s000_0_alpha+gama0_00+s000_1_beta)) - log(exp(s000_0_alpha+gama0_01+s010_1_beta));
			tmp1=th(alph+gama1_0+beta1_0);
			tmp1b=th(alph+gama1_1+beta1_1);
			//L_out[i]= log(exp(tmp1)) - log(exp(tmp1b)); 
			//cout<<"j: "<<j<<'\n';
		}
		
		else if(j==1) {
			//L_out[i]= log(exp(s000_1_alpha+gama1_00+s000_2_beta)+exp(s010_1_alpha+gama1_10+s010_2_beta)) - log(exp(s000_1_alpha+gama1_01+s011_2_beta)+exp(s010_1_alpha+gama1_11+s001_2_beta));
			
			tmp1=th(alph1_0+gama2_0+beta2_00);
			tmp2=th(alph1_1+gama2_0+beta2_11);
			
			tmp1b=th(alph1_1+gama2_1+beta2_10);
			tmp2b=th(alph1_0+gama2_1+beta2_01);
					
			/*if(tmp1>tmp2)
				max1=tmp1;
			else
				max1=tmp2;
			if(tmp1b>tmp2b)
				max2=tmp1b;
			else
				max2=tmp2b;*/
					
			//L_out[i]= max1+log(1+exp(-1*abs(tmp1-tmp2))) - (max2+log(1+exp(-1*abs(tmp1b-tmp2b))));
			//L_out[i]= maxx(tmp1,tmp2)+log(1+exp(-1*abs(tmp1-tmp2))) - (maxx(tmp1b,tmp2b)+log(1+exp(-1*abs(tmp1b-tmp2b))));
			//cout<<"j: "<<j<<'\n';
		}
	
		else if(j==2) {
			//L_out[i]= log(exp(s000_2_alpha+gama2_00+s000_3_beta)+exp(s001_2_alpha+gama2_10+s001_3_beta)+exp(s010_2_alpha+gama2_00b+s010_3_beta)+exp(s011_2_alpha+gama2_10b+s011_3_beta)) - log(exp(s000_2_alpha+gama2_01+s001_3_beta)+exp(s001_2_alpha+gama2_11+s000_3_beta)+exp(s010_2_alpha+gama2_01b+s011_3_beta)+exp(s011_2_alpha+gama2_11b+s010_3_beta));
			
			tmp1=th(alph2_00+gama3_0+beta3_000);
			tmp2=th(alph2_01+gama3_0+beta3_011);
			tmp3=th(alph2_10+gama3_0+beta3_101);
			tmp4=th(alph2_11+gama3_0+beta3_110);
			
			tmp1b=th(alph2_00+gama3_1+beta3_001);
			tmp2b=th(alph2_01+gama3_1+beta3_010);
			tmp3b=th(alph2_10+gama3_1+beta3_100);
			tmp4b=th(alph2_11+gama3_1+beta3_111);
			
			//L_out[i]= log(exp(tmp1)+exp(tmp2)+exp(tmp3)+exp(tmp4)) - log(exp(tmp1b)+exp(tmp2b)+exp(tmp3b)+exp(tmp4b));
			//max(a,b,c,d)=max((a,b,c),d)=max(max2(max1(a,b),c)),d)=max(max2,d)=max(max2,tmp4)=max3+log(1+exp(-1*abs(max2-tmp4)))
					
			/*if(tmp1>tmp2)
				max1=tmp1;
			else
				max1=tmp2;
			if(max1>tmp3)
				max2=max1;
			else
				max2=tmp3;	
			if(max2>tmp4)
				max3=max2;
			else
				max3=tmp4;
			
			if(tmp1b>tmp2b)
				max1b=tmp1b;
			else
				max1b=tmp2b;
			if(max1b>tmp3b)
				max2b=max1b;
			else
				max2b=tmp3b;	
			if(max2b>tmp4b)
				max3b=max2b;
			else
				max3b=tmp4b;*/
			//L_out[i]=max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b))));
		}	
		
		else if(j==3) {
			tmp1=th(alph3_000+gama4_0+beta4_0000);
			tmp2=th(alph3_001+gama4_0+beta4_0010);
			tmp3=th(alph3_010+gama4_0+beta4_0101);
			tmp4=th(alph3_011+gama4_0+beta4_0111);
			tmp5=th(alph3_000+gama4_0+beta4_1001);
			tmp6=th(alph3_001+gama4_0+beta4_1011);
			tmp7=th(alph3_010+gama4_0+beta4_1100);
			tmp8=th(alph3_011+gama4_0+beta4_1110);
			
			tmp1b=th(alph3_000+gama4_1+beta4_0001);
			tmp2b=th(alph3_001+gama4_1+beta4_0011);
			tmp3b=th(alph3_010+gama4_1+beta4_0100);
			tmp4b=th(alph3_011+gama4_1+beta4_0110);
			tmp5b=th(alph3_000+gama4_1+beta4_1000);
			tmp6b=th(alph3_001+gama4_1+beta4_1010);
			tmp7b=th(alph3_010+gama4_1+beta4_1101);
			tmp8b=th(alph3_011+gama4_1+beta4_1111);
		}
	
		else if(j==4) {
			tmp1=th(alph4_0000+gama5_0+beta5_0000);
			tmp2=th(alph4_0001+gama5_0+beta5_0011);
			tmp3=th(alph4_0010+gama5_0+beta5_0101);
			tmp4=th(alph4_0011+gama5_0+beta5_0110);
			tmp5=th(alph4_0100+gama5_0+beta5_1001);
			tmp6=th(alph4_0101+gama5_0+beta5_1010);
			tmp7=th(alph4_0110+gama5_0+beta5_1100);
			tmp8=th(alph4_0111+gama5_0+beta5_1111);
			tmp9=th(alph4_1000+gama5_0+beta5_0001);
			tmp10=th(alph4_1001+gama5_0+beta5_0010);
			tmp11=th(alph4_1010+gama5_0+beta5_0100);
			tmp12=th(alph4_1011+gama5_0+beta5_0111);
			tmp13=th(alph4_1100+gama5_0+beta5_1000);
			tmp14=th(alph4_1101+gama5_0+beta5_1011);
			tmp15=th(alph4_1110+gama5_0+beta5_1101);
			tmp16=th(alph4_1111+gama5_0+beta5_1110);
			
			tmp1b=th(alph4_0000+gama5_1+beta5_0001);
			tmp2b=th(alph4_0001+gama5_1+beta5_0010);
			tmp3b=th(alph4_0010+gama5_1+beta5_0100);
			tmp4b=th(alph4_0011+gama5_1+beta5_0111);
			tmp5b=th(alph4_0100+gama5_1+beta5_1000);
			tmp6b=th(alph4_0101+gama5_1+beta5_1011);
			tmp7b=th(alph4_0110+gama5_1+beta5_1101);
			tmp8b=th(alph4_0111+gama5_1+beta5_1110);
			tmp9b=th(alph4_1000+gama5_1+beta5_0000);
			tmp10b=th(alph4_1001+gama5_1+beta5_0011);
			tmp11b=th(alph4_1010+gama5_1+beta5_0101);
			tmp12b=th(alph4_1011+gama5_1+beta5_0110);
			tmp13b=th(alph4_1100+gama5_1+beta5_1001);
			tmp14b=th(alph4_1101+gama5_1+beta5_1010);
			tmp15b=th(alph4_1110+gama5_1+beta5_1100);
			tmp16b=th(alph4_1111+gama5_1+beta5_1111);
		}
	
		else if(j==5) {
			tmp1=th(alph5_0000+gama6_0+beta6_0000);
			tmp2=th(alph5_0001+gama6_0+beta6_0010);
			tmp3=th(alph5_0010+gama6_0+beta6_0001);
			tmp4=th(alph5_0011+gama6_0+beta6_0011);
			tmp5=th(alph5_0100+gama6_0+beta6_0100);
			tmp6=th(alph5_0101+gama6_0+beta6_0110);
			tmp7=th(alph5_0110+gama6_0+beta6_0101);
			tmp8=th(alph5_0111+gama6_0+beta6_0111);
			tmp9=th(alph5_1000+gama6_0+beta6_1001);
			tmp10=th(alph5_1001+gama6_0+beta6_1011);
			tmp11=th(alph5_1010+gama6_0+beta6_1000);
			tmp12=th(alph5_1011+gama6_0+beta6_1010);
			tmp13=th(alph5_1100+gama6_0+beta6_1101);
			tmp14=th(alph5_1101+gama6_0+beta6_1111);
			tmp15=th(alph5_1110+gama6_0+beta6_1100);
			tmp16=th(alph5_1111+gama6_0+beta6_1110);
			
			tmp1b=th(alph5_0000+gama6_1+beta6_0001);
			tmp2b=th(alph5_0001+gama6_1+beta6_0011);
			tmp3b=th(alph5_0010+gama6_1+beta6_0000);
			tmp4b=th(alph5_0011+gama6_1+beta6_0010);
			tmp5b=th(alph5_0100+gama6_1+beta6_0101);
			tmp6b=th(alph5_0101+gama6_1+beta6_0111);
			tmp7b=th(alph5_0110+gama6_1+beta6_0100);
			tmp8b=th(alph5_0111+gama6_1+beta6_0110);
			tmp9b=th(alph5_1000+gama6_1+beta6_1000);
			tmp10b=th(alph5_1001+gama6_1+beta6_1010);
			tmp11b=th(alph5_1010+gama6_1+beta6_1001);
			tmp12b=th(alph5_1011+gama6_1+beta6_1011);
			tmp13b=th(alph5_1100+gama6_1+beta6_1100);
			tmp14b=th(alph5_1101+gama6_1+beta6_1110);
			tmp15b=th(alph5_1110+gama6_1+beta6_1101);
			tmp16b=th(alph5_1111+gama6_1+beta6_1111);
		}
	
		else if(j==6) {
			tmp1=th(alph6_0000+gama7_0+beta7_0000);
			tmp2=th(alph6_0001+gama7_0+beta7_0011);
			tmp3=th(alph6_0010+gama7_0+beta7_0100);
			tmp4=th(alph6_0011+gama7_0+beta7_0111);
			tmp5=th(alph6_0100+gama7_0+beta7_0001);
			tmp6=th(alph6_0101+gama7_0+beta7_0010);
			tmp7=th(alph6_0110+gama7_0+beta7_0101);
			tmp8=th(alph6_0111+gama7_0+beta7_0110);
			tmp9=th(alph6_1000+gama7_0+beta7_1000);
			tmp10=th(alph6_1001+gama7_0+beta7_1011);
			tmp11=th(alph6_1010+gama7_0+beta7_1100);
			tmp12=th(alph6_1011+gama7_0+beta7_1111);
			tmp13=th(alph6_1100+gama7_0+beta7_1001);
			tmp14=th(alph6_1101+gama7_0+beta7_1010);
			tmp15=th(alph6_1110+gama7_0+beta7_1101);
			tmp16=th(alph6_1111+gama7_0+beta7_1111);
			
			tmp1b=th(alph6_0000+gama7_1+beta7_0001);
			tmp2b=th(alph6_0001+gama7_1+beta7_0010);
			tmp3b=th(alph6_0010+gama7_1+beta7_0101);
			tmp4b=th(alph6_0011+gama7_1+beta7_0110);
			tmp5b=th(alph6_0100+gama7_1+beta7_0000);
			tmp6b=th(alph6_0101+gama7_1+beta7_0011);
			tmp7b=th(alph6_0110+gama7_1+beta7_0100);
			tmp8b=th(alph6_0111+gama7_1+beta7_0111);
			tmp9b=th(alph6_1000+gama7_1+beta7_1001);
			tmp10b=th(alph6_1001+gama7_1+beta7_1010);
			tmp11b=th(alph6_1010+gama7_1+beta7_1101);
			tmp12b=th(alph6_1011+gama7_1+beta7_1110);
			tmp13b=th(alph6_1100+gama7_1+beta7_1000);
			tmp14b=th(alph6_1101+gama7_1+beta7_1011);
			tmp15b=th(alph6_1110+gama7_1+beta7_1100);
			tmp16b=th(alph6_1111+gama7_1+beta7_1111);
		}
	
		else if(j==7) {
			tmp1=th(alph7_0000+gama8_0+beta8_0000);
			tmp2=th(alph7_0001+gama8_0+beta8_0011);
			tmp3=th(alph7_0010+gama8_0+beta8_0001);
			tmp4=th(alph7_0011+gama8_0+beta8_0010);
			tmp5=th(alph7_0100+gama8_0+beta8_0101);
			tmp6=th(alph7_0101+gama8_0+beta8_0110);
			tmp7=th(alph7_0110+gama8_0+beta8_0100);
			tmp8=th(alph7_0111+gama8_0+beta8_0111);
			tmp9=th(alph7_1000+gama8_0+beta8_1001);
			tmp10=th(alph7_1001+gama8_0+beta8_1010);
			tmp11=th(alph7_1010+gama8_0+beta8_1000);
			tmp12=th(alph7_1011+gama8_0+beta8_1011);
			tmp13=th(alph7_1100+gama8_0+beta8_1100);
			tmp14=th(alph7_1101+gama8_0+beta8_1111);
			tmp15=th(alph7_1110+gama8_0+beta8_1101);
			tmp16=th(alph7_1111+gama8_0+beta8_1110);
			
			tmp1b=th(alph7_0000+gama8_1+beta8_0001);
			tmp2b=th(alph7_0001+gama8_1+beta8_0010);
			tmp3b=th(alph7_0010+gama8_1+beta8_0000);
			tmp4b=th(alph7_0011+gama8_1+beta8_0011);
			tmp5b=th(alph7_0100+gama8_1+beta8_0100);
			tmp6b=th(alph7_0101+gama8_1+beta8_0111);
			tmp7b=th(alph7_0110+gama8_1+beta8_0101);
			tmp8b=th(alph7_0111+gama8_1+beta8_0110);
			tmp9b=th(alph7_1000+gama8_1+beta8_1000);
			tmp10b=th(alph7_1001+gama8_1+beta8_1011);
			tmp11b=th(alph7_1010+gama8_1+beta8_1001);
			tmp12b=th(alph7_1011+gama8_1+beta8_1010);
			tmp13b=th(alph7_1100+gama8_1+beta8_1101);
			tmp14b=th(alph7_1101+gama8_1+beta8_1110);
			tmp15b=th(alph7_1110+gama8_1+beta8_1100);
			tmp16b=th(alph7_1111+gama8_1+beta8_1111);
		}
	
		else if(j==8) {
			tmp1=th(alph8_0000+gama9_0+beta9_0000);
			tmp2=th(alph8_0001+gama9_0+beta9_0011);
			tmp3=th(alph8_0010+gama9_0+beta9_0001);
			tmp4=th(alph8_0011+gama9_0+beta9_0010);
			tmp5=th(alph8_0100+gama9_0+beta9_0100);
			tmp6=th(alph8_0101+gama9_0+beta9_0111);
			tmp7=th(alph8_0110+gama9_0+beta9_0101);
			tmp8=th(alph8_0111+gama9_0+beta9_0110);
			tmp9=th(alph8_1000+gama9_0+beta9_1001);
			tmp10=th(alph8_1001+gama9_0+beta9_1010);
			tmp11=th(alph8_1010+gama9_0+beta9_1000);
			tmp12=th(alph8_1011+gama9_0+beta9_1011);
			tmp13=th(alph8_1100+gama9_0+beta9_1101);
			tmp14=th(alph8_1101+gama9_0+beta9_1110);
			tmp15=th(alph8_1110+gama9_0+beta9_1100);
			tmp16=th(alph8_1111+gama9_0+beta9_1111);
			
			tmp1b=th(alph8_0000+gama9_1+beta9_0001);
			tmp2b=th(alph8_0001+gama9_1+beta9_0010);
			tmp3b=th(alph8_0010+gama9_1+beta9_0000);
			tmp4b=th(alph8_0011+gama9_1+beta9_0011);
			tmp5b=th(alph8_0100+gama9_1+beta9_0101);
			tmp6b=th(alph8_0101+gama9_1+beta9_0110);
			tmp7b=th(alph8_0110+gama9_1+beta9_0100);
			tmp8b=th(alph8_0111+gama9_1+beta9_0111);
			tmp9b=th(alph8_1000+gama9_1+beta9_1000);
			tmp10b=th(alph8_1001+gama9_1+beta9_1011);
			tmp11b=th(alph8_1010+gama9_1+beta9_1001);
			tmp12b=th(alph8_1011+gama9_1+beta9_1010);
			tmp13b=th(alph8_1100+gama9_1+beta9_1100);
			tmp14b=th(alph8_1101+gama9_1+beta9_1111);
			tmp15b=th(alph8_1110+gama9_1+beta9_1101);
			tmp16b=th(alph8_1111+gama9_1+beta9_1110);
		}
	
		else if(j==9) {
			tmp1=th(alph9_0000+gama10_0+beta10_0000);
			tmp2=th(alph9_0001+gama10_0+beta10_0010);
			tmp3=th(alph9_0010+gama10_0+beta10_0100);
			tmp4=th(alph9_0011+gama10_0+beta10_0110);
			tmp5=th(alph9_0100+gama10_0+beta10_0001);
			tmp6=th(alph9_0101+gama10_0+beta10_0011);
			tmp7=th(alph9_0110+gama10_0+beta10_0101);
			tmp8=th(alph9_0111+gama10_0+beta10_0111);
			tmp9=th(alph9_1000+gama10_0+beta10_1000);
			tmp10=th(alph9_1001+gama10_0+beta10_1010);
			tmp11=th(alph9_1010+gama10_0+beta10_1100);
			tmp12=th(alph9_1011+gama10_0+beta10_1110);
			tmp13=th(alph9_1100+gama10_0+beta10_1001);
			tmp14=th(alph9_1101+gama10_0+beta10_1011);
			tmp15=th(alph9_1110+gama10_0+beta10_1101);
			tmp16=th(alph9_1111+gama10_0+beta10_1111);
			
			tmp1b=th(alph9_0000+gama10_1+beta10_0001);
			tmp2b=th(alph9_0001+gama10_1+beta10_0011);
			tmp3b=th(alph9_0010+gama10_1+beta10_0101);
			tmp4b=th(alph9_0011+gama10_1+beta10_0111);
			tmp5b=th(alph9_0100+gama10_1+beta10_0000);
			tmp6b=th(alph9_0101+gama10_1+beta10_0010);
			tmp7b=th(alph9_0110+gama10_1+beta10_0100);
			tmp8b=th(alph9_0111+gama10_1+beta10_0110);
			tmp9b=th(alph9_1000+gama10_1+beta10_1001);
			tmp10b=th(alph9_1001+gama10_1+beta10_1011);
			tmp11b=th(alph9_1010+gama10_1+beta10_1101);
			tmp12b=th(alph9_1011+gama10_1+beta10_1111);
			tmp13b=th(alph9_1100+gama10_1+beta10_1000);
			tmp14b=th(alph9_1101+gama10_1+beta10_1010);
			tmp15b=th(alph9_1110+gama10_1+beta10_1100);
			tmp16b=th(alph9_1111+gama10_1+beta10_1110);
		}
	
		else if(j==10) {
			tmp1=th(alph10_0000+gama11_0+beta11_0000);
			tmp2=th(alph10_0001+gama11_0+beta11_0010);
			tmp3=th(alph10_0010+gama11_0+beta11_0100);
			tmp4=th(alph10_0011+gama11_0+beta11_0110);
			tmp5=th(alph10_0100+gama11_0+beta11_0001);
			tmp6=th(alph10_0101+gama11_0+beta11_0011);
			tmp7=th(alph10_0110+gama11_0+beta11_0101);
			tmp8=th(alph10_0111+gama11_0+beta11_0111);
			tmp9=th(alph10_1000+gama11_0+beta11_1000);
			tmp10=th(alph10_1001+gama11_0+beta11_1010);
			tmp11=th(alph10_1010+gama11_0+beta11_1100);
			tmp12=th(alph10_1011+gama11_0+beta11_1110);
			tmp13=th(alph10_1100+gama11_0+beta11_1001);
			tmp14=th(alph10_1101+gama11_0+beta11_1011);
			tmp15=th(alph10_1110+gama11_0+beta11_1101);
			tmp16=th(alph10_1111+gama11_0+beta11_1111);
			
			tmp1b=th(alph10_0000+gama11_1+beta11_0001);
			tmp2b=th(alph10_0001+gama11_1+beta11_0011);
			tmp3b=th(alph10_0010+gama11_1+beta11_0101);
			tmp4b=th(alph10_0011+gama11_1+beta11_0111);
			tmp5b=th(alph10_0100+gama11_1+beta11_0000);
			tmp6b=th(alph10_0101+gama11_1+beta11_0010);
			tmp7b=th(alph10_0110+gama11_1+beta11_0100);
			tmp8b=th(alph10_0111+gama11_1+beta11_0110);
			tmp9b=th(alph10_1000+gama11_1+beta11_1001);
			tmp10b=th(alph10_1001+gama11_1+beta11_1011);
			tmp11b=th(alph10_1010+gama11_1+beta11_1101);
			tmp12b=th(alph10_1011+gama11_1+beta11_1111);
			tmp13b=th(alph10_1100+gama11_1+beta11_1000);
			tmp14b=th(alph10_1101+gama11_1+beta11_1010);
			tmp15b=th(alph10_1110+gama11_1+beta11_1100);
			tmp16b=th(alph10_1111+gama11_1+beta11_1110);
		}
		
		else if(j==11) {
			tmp1=th(alph11_0000+gama12_0+beta12_000);
			tmp2=th(alph11_0011+gama12_0+beta12_011);
			tmp3=th(alph11_0101+gama12_0+beta12_101);
			tmp4=th(alph11_0110+gama12_0+beta12_110);
			tmp5=th(alph11_1001+gama12_0+beta12_001);
			tmp6=th(alph11_1010+gama12_0+beta12_010);
			tmp7=th(alph11_1100+gama12_0+beta12_100);
			tmp8=th(alph11_1111+gama12_0+beta12_111);
			
			tmp1b=th(alph11_0001+gama12_1+beta12_001);
			tmp2b=th(alph11_0010+gama12_1+beta12_010);
			tmp3b=th(alph11_0100+gama12_1+beta12_100);
			tmp4b=th(alph11_0111+gama12_1+beta12_111);
			tmp5b=th(alph11_1000+gama12_1+beta12_000);
			tmp6b=th(alph11_1011+gama12_1+beta12_011);
			tmp7b=th(alph11_1101+gama12_1+beta12_101);
			tmp8b=th(alph11_1110+gama12_1+beta12_110);
		}
		
		else if(j==12) {
			tmp1=th(alph12_000+gama13_0+beta13_00);
			tmp2=th(alph12_010+gama13_0+beta13_01);
			tmp3=th(alph12_100+gama13_0+beta13_10);
			tmp4=th(alph12_110+gama13_0+beta13_11);
			
			tmp1b=th(alph12_001+gama13_1+beta13_00);
			tmp2b=th(alph12_011+gama13_1+beta13_01);
			tmp3b=th(alph12_101+gama13_1+beta13_10);
			tmp4b=th(alph12_111+gama13_1+beta13_11);
		}
		
		else if(j==13) {
			tmp1=th(alph13_00+gama14_0+beta14_0);
			tmp2=th(alph13_01+gama14_0+beta14_1);
			
			tmp1b=th(alph13_10+gama14_1+beta14_0);
			tmp2b=th(alph13_11+gama14_1+beta14_1);
		}
		
		else if(j==14) {
			tmp1=th(alph14_0+gama15_0+beta15_0);
			tmp1b=th(alph14_1+gama15_1+beta15_0);
		}
		
		if(j==2 || j==12) {	
			max2=maxx(maxx(tmp1,tmp2),tmp3);	
			max3=maxx(max2,tmp4);	
			max2b=maxx(maxx(tmp1b,tmp2b),tmp3b);	
			max3b=maxx(max2b,tmp4b);
		}
		else if(j==3 || j==11) {
			max2=maxx(maxx(maxx(maxx(maxx(maxx(tmp1,tmp2),tmp3),tmp4),tmp5),tmp6),tmp7);	
			max3=maxx(max2,tmp8);	
			max2b=maxx(maxx(maxx(maxx(maxx(maxx(tmp1b,tmp2b),tmp3b),tmp4b),tmp5b),tmp6b),tmp7b);		
			max3b=maxx(max2b,tmp8b);
			//L_out[i]=max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b))));		
		}
		else { //if(j==4,5,6,7,8,9,10)
			max2=maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(tmp1,tmp2),tmp3),tmp4),tmp5),tmp6),tmp7),tmp8),tmp9),tmp10),tmp11),tmp12),tmp13),
				 tmp14),tmp15);	
			max3=maxx(max2,tmp16);	
			max2b=maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(maxx(tmp1b,tmp2b),tmp3b),tmp4b),tmp5b),tmp6b),tmp7b),tmp8b),tmp9b),tmp10b),tmp11b),
				  tmp12b),tmp13b),tmp14b),tmp15b);			
			max3b=maxx(max2b,tmp16b);
			//L_out[i]=max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b))));	
		}
		
		if(!j || j==14) 
			//L_out[i]= th(log(exp(tmp1)) - log(exp(tmp1b))); 
			L_out[i]= th(tmp1+log(1+exp(-1*abs(tmp1))) - (tmp1b+log(1+exp(-1*abs(tmp1b)))));
		else if(j==1 || j==13) 
			L_out[i]= th(maxx(tmp1,tmp2)+log(1+exp(-1*abs(tmp1-tmp2))) - (maxx(tmp1b,tmp2b)+log(1+exp(-1*abs(tmp1b-tmp2b)))));
		else if(j==2 || j==12) 
			L_out[i]=th(max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b)))));
		else if(j==3 || j==11) 
			L_out[i]=th(max3+log(1+exp(-1*abs(max2-tmp8))) - (max3b+log(1+exp(-1*abs(max2b-tmp8b)))));
		else //for j==4,5,6,7,8,9,10 
			L_out[i]=th(max3+log(1+exp(-1*abs(max2-tmp16))) - (max3b+log(1+exp(-1*abs(max2b-tmp16b)))));
			
	}
	//cout<<"a: "<<a<<endl;	
	//cout<<'\n'<<"L_out: "<<'\n';
	for(i=0;i<n;i++)
		if(isnan(L_out[i])) 
			cout<<"L_out: "<<L_out[i]<<" ";
	
	
}

