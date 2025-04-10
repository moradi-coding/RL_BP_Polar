
//BCJR calclations using the [7,4] Hamming trellis in Vivian paper, calculations based on Ryan-Lin book
//flooding GLDPC decoder for all GCNs
void bcjr_log(int a, int l) { 
		
	long i,j;
	double tmp1,tmp2,tmp3,tmp4,tmp1b,tmp2b,tmp3b,tmp4b,max1,max2,max3,max1b,max2b,max3b,c,y_sum,b,d,num,s000_6_alpha_old;
	d=0.5; //llr_max=25;
	
	if(!l) {
		s000_0_alpha=s000_1_alpha=s000_2_alpha=s000_3_alpha=s000_4_alpha=s000_5_alpha=s000_6_alpha=s000_7_alpha=
s000_0_beta=s000_1_beta=s000_2_beta=s000_3_beta=s000_4_beta=s000_5_beta=s000_6_beta=s000_7_beta=
s001_2_alpha=s001_3_alpha=s001_4_alpha=s001_2_beta=s001_3_beta=s001_4_beta=
s010_1_alpha=s010_2_alpha=s010_3_alpha=s010_4_alpha=s010_5_alpha=
s010_1_beta=s010_2_beta=s010_3_beta=s010_4_beta=s010_5_beta=
s011_2_alpha=s011_3_alpha=s011_4_alpha=s011_2_beta=s011_3_beta=s011_4_beta=
s100_4_alpha=s100_4_beta=
s101_4_alpha=s101_5_alpha=s101_6_alpha=s101_4_beta=s101_5_beta=s101_6_beta=
s110_4_alpha=s110_4_beta=
s111_4_alpha=s111_4_beta=s111_5_alpha=s111_5_beta=0;
//gama0_00=gama0_01=gama1_00=gama1_10=gama1_01=gama1_11=gama2_00=gama2_10=gama2_00b=gama2_10b=gama2_01=gama2_11=gama2_01b=gama2_11b, 
//gama3_00=gama3_00b=gama3_00c=gama3_00d=gama3_01=gama3_01b=gama3_01c=gama3_01d=gama4_00=gama4_00b=gama4_00c=gama4_00d=gama4_01=gama4_11c=gama4_11b,
//gama4_11=gama5_00=gama5_00b=gama5_01b=gama5_01=gama6_00=gama6_01=0;
	}
	
	
	c=-1/(2*varn);
	//for(a=0;a<num_cls;a++) {
		//initialization of branch metrics
		for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
			i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'
			//b=exp(-0.5*L_in[i])/(1+exp(-1*L_in[i]));

			//cout<<"j, b, LLR: "<<j<<", "<<b<<", "<<L_in[i]<<endl;
			//b=L_in[i];
			b=L_in[i]-LR[i]; //gives the extrinsic LLR
			
			if(L_in[i]>llr_max) 
				L_in[i]=llr_max;
			else if(L_in[i]<-1*llr_max) 
				L_in[i]=-1*llr_max;
			
			if(!j) {
				gama0_00=c*pow(abs(y[i]-1),2)+d*1*b; //0 mapped to 1	
				gama0_01=c*pow(abs(y[i]+1),2)+d*-1*b; //1 mapped to -1	
			}	
			else if(j==1) {
				gama1_00=c*pow(abs(y[i]-1),2)+d*1*b; 
				gama1_01=c*pow(abs(y[i]+1),2)+d*-1*b; 
				gama1_10=gama1_00;
				gama1_11=gama1_01;
			}	
			else if(j==2) {
				gama2_00=c*pow(abs(y[i]-1),2)+d*1*b; 
				gama2_01=c*pow(abs(y[i]+1),2)+d*-1*b; 
				gama2_10=gama2_00;
				gama2_00b=gama2_00;
				gama2_10b=gama2_00;
				gama2_11=gama2_01;
				gama2_01b=gama2_01;
				gama2_11b=gama2_01;
			}	
			else if(j==3) {
				gama3_00=c*pow(abs(y[i]-1),2)+d*1*b; 
				gama3_01=c*pow(abs(y[i]+1),2)+d*-1*b; 
				gama3_00b=gama3_00;
				gama3_00c=gama3_00;
				gama3_00d=gama3_00;
				gama3_01b=gama3_01;
				gama3_01c=gama3_01;
				gama3_01d=gama3_01;
			}	
			else if(j==4) {
				gama4_00=c*pow(abs(y[i]-1),2)+d*1*b; 
				gama4_01=c*pow(abs(y[i]+1),2)+d*-1*b; 
				gama4_00b=gama4_00;
				gama4_00c=gama4_00;
				gama4_00d=gama4_00;
				gama4_11b=gama4_01;
				gama4_11c=gama4_01;
				gama4_11=gama4_01;
			}	
			else if(j==5) {
				gama5_00=c*pow(abs(y[i]-1),2)+d*1*b; 
				gama5_01=c*pow(abs(y[i]+1),2)+d*-1*b; 
				gama5_00b=gama5_00;
				gama5_01b=gama5_01;
			}	
			else if(j==6) {
				gama6_00=c*pow(abs(y[i]-1),2)+d*1*b; //0 mapped to 1	
				gama6_01=c*pow(abs(y[i]+1),2)+d*-1*b; //1 mapped to -1	
			}	
		}
		//cout<<"L_in: "; for(i=0;i<n;i++) cout<<L_in[i]<< " "; cout<<'\n'<<'\n';
		
		//forward metrics
		s000_0_alpha=log(1);
		if(abs(s000_1_alpha)<llr_max2) s000_1_alpha=s000_0_alpha+gama0_00;
		if(s000_1_alpha>llr_max) 
			s000_1_alpha=llr_max;
		else if(s000_1_alpha<-1*llr_max) 
			s000_1_alpha=-1*llr_max;
		
		if(abs(s000_2_alpha)<llr_max2) s000_2_alpha=s000_1_alpha+gama1_00;
		if(s000_2_alpha>llr_max) 
			s000_2_alpha=llr_max;
		else if(s000_2_alpha<-1*llr_max) 
			s000_2_alpha=-1*llr_max;
		
		if(abs(s010_1_alpha)<llr_max2) s010_1_alpha=s000_0_alpha+gama0_01;
		if(s010_1_alpha>llr_max) 
			s010_1_alpha=llr_max;
		else if(s010_1_alpha<-1*llr_max) 
			s010_1_alpha=-1*llr_max;
		
		if(abs(s001_2_alpha)<llr_max2) s001_2_alpha=s010_1_alpha+gama1_11;
		if(s001_2_alpha>llr_max) 
			s001_2_alpha=llr_max;
		else if(s001_2_alpha<-1*llr_max) 
			s001_2_alpha=-1*llr_max;
			
		//cout<<"s001_2_alpha: "<<s001_2_alpha<<endl;
		//s000_3_alpha=log(exp(s000_2_alpha+gama2_00) + exp(s001_2_alpha+gama2_11));
		tmp1=s000_2_alpha+gama2_00;
		tmp2=s001_2_alpha+gama2_11; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s000_3_alpha)<llr_max2) s000_3_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_3_alpha>llr_max) 
			s000_3_alpha=llr_max;
		else if(s000_3_alpha<-1*llr_max) 
			s000_3_alpha=-1*llr_max;
		
		//cout<<"s000_3_alpha: "<<s000_3_alpha<<endl;
		
		if(abs(s000_4_alpha)<llr_max2) s000_4_alpha=s000_3_alpha+gama3_00;
		if(s000_4_alpha>llr_max) 
			s000_4_alpha=llr_max;
		else if(s000_4_alpha<-1*llr_max) 
			s000_4_alpha=-1*llr_max;
		
		if(abs(s010_2_alpha)<llr_max2) s010_2_alpha=s010_1_alpha+gama1_10;
		if(s010_2_alpha>llr_max) 
			s010_2_alpha=llr_max;
		else if(s010_2_alpha<-1*llr_max) 
			s010_2_alpha=-1*llr_max;
			
		if(abs(s011_2_alpha)<llr_max2) s011_2_alpha=s000_1_alpha+gama1_01;
		if(s011_2_alpha>llr_max) 
			s011_2_alpha=llr_max;
		else if(s011_2_alpha<-1*llr_max) 
			s011_2_alpha=-1*llr_max;
			
		//s010_3_alpha=log(exp(s010_2_alpha+gama2_00b) + exp(s011_2_alpha+gama2_11b));
		tmp1=s010_2_alpha+gama2_00b;
		tmp2=s011_2_alpha+gama2_11b; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s010_3_alpha)<llr_max2) s010_3_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s010_3_alpha>llr_max) 
			s010_3_alpha=llr_max;
		else if(s010_3_alpha<-1*llr_max) 
			s010_3_alpha=-1*llr_max;
			
		if(abs(s110_4_alpha)<llr_max2) s110_4_alpha=s010_3_alpha+gama3_01c;
		if(s110_4_alpha>llr_max) 
			s110_4_alpha=llr_max;
		else if(s110_4_alpha<-1*llr_max) 
			s110_4_alpha=-1*llr_max;
		
		//s000_5_alpha=log(exp(s000_4_alpha+gama4_00) + exp(s110_4_alpha+gama4_11));	
		tmp1=s000_4_alpha+gama4_00;
		tmp2=s110_4_alpha+gama4_11; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s000_5_alpha)<llr_max2) s000_5_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_5_alpha>llr_max) 
			s000_5_alpha=llr_max;
		else if(s000_5_alpha<-1*llr_max) 
			s000_5_alpha=-1*llr_max;
			
		//s001_3_alpha=log(exp(s001_2_alpha+gama2_10) + exp(s000_2_alpha+gama2_01));
		tmp1=s001_2_alpha+gama2_10;
		tmp2=s000_2_alpha+gama2_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s001_3_alpha)<llr_max2) s001_3_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s001_3_alpha>llr_max) 
			s001_3_alpha=llr_max;
		else if(s001_3_alpha<-1*llr_max) 
			s001_3_alpha=-1*llr_max;
			
		if(abs(s101_4_alpha)<llr_max2) s101_4_alpha=s001_3_alpha+gama3_01b;
		if(s101_4_alpha>llr_max) 
			s101_4_alpha=llr_max;
		else if(s101_4_alpha<-1*llr_max) 
			s101_4_alpha=-1*llr_max;
			
		//s011_3_alpha=log(exp(s011_2_alpha+gama2_10b) + exp(s010_2_alpha+gama2_01b));
		tmp1=s011_2_alpha+gama2_10b;
		tmp2=s010_2_alpha+gama2_01b; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s011_3_alpha)<llr_max2) s011_3_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s011_3_alpha>llr_max) 
			s011_3_alpha=llr_max;
		else if(s011_3_alpha<-1*llr_max) 
			s011_3_alpha=-1*llr_max;
			
		if(abs(s011_4_alpha)<llr_max2) s011_4_alpha=s011_3_alpha+gama3_00d;
		if(s011_4_alpha>llr_max) 
			s011_4_alpha=llr_max;
		else if(s011_4_alpha<-1*llr_max) 
			s011_4_alpha=-1*llr_max;
	
		//s101_5_alpha=log(exp(s101_4_alpha+gama4_00c) + exp(s011_4_alpha+gama4_11c));
		tmp1=s101_4_alpha+gama4_00c;
		tmp2=s011_4_alpha+gama4_11; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s101_5_alpha)<llr_max2) s101_5_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s101_5_alpha>llr_max) 
			s101_5_alpha=llr_max;
		else if(s101_5_alpha<-1*llr_max) 
			s101_5_alpha=-1*llr_max;
			
		if(abs(s010_4_alpha)<llr_max2) s010_4_alpha=s010_3_alpha+gama3_00c;
		if(s010_4_alpha>llr_max) 
			s010_4_alpha=llr_max;
		else if(s010_4_alpha<-1*llr_max) 
			s010_4_alpha=-1*llr_max;
		
		if(abs(s100_4_alpha)<llr_max2) s100_4_alpha=s000_3_alpha+gama3_01;
		if(s100_4_alpha>llr_max) 
			s100_4_alpha=llr_max;
		else if(s100_4_alpha<-1*llr_max) 
			s100_4_alpha=-1*llr_max;
			
		//s010_5_alpha=log(exp(s010_4_alpha+gama4_00b) + exp(s100_4_alpha+gama4_11b));
		tmp1=s010_4_alpha+gama4_00b;
		tmp2=s100_4_alpha+gama4_11b; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s010_5_alpha)<llr_max2) s010_5_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s010_5_alpha>llr_max) 
			s010_5_alpha=llr_max;
		else if(s010_5_alpha<-1*llr_max) 
			s010_5_alpha=-1*llr_max;
			
		//s101_6_alpha=log(exp(s101_5_alpha+gama5_00b) + exp(s010_5_alpha+gama5_01b));
		tmp1=s101_5_alpha+gama5_00b;
		tmp2=s010_5_alpha+gama5_01b; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s101_6_alpha)<llr_max2) s101_6_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s101_6_alpha>llr_max) 
			s101_6_alpha=llr_max;
		else if(s101_6_alpha<-1*llr_max) 
			s101_6_alpha=-1*llr_max;
		
		if(abs(s001_4_alpha)<llr_max2) s001_4_alpha=s001_3_alpha+gama3_00b;
		if(s001_4_alpha>llr_max) 
			s001_4_alpha=llr_max;
		else if(s001_4_alpha<-1*llr_max) 
			s001_4_alpha=-1*llr_max;
			
		if(abs(s111_4_alpha)<llr_max2) s111_4_alpha=s011_3_alpha+gama3_01d;
		if(s111_4_alpha>llr_max) 
			s111_4_alpha=llr_max;
		else if(s111_4_alpha<-1*llr_max) 
			s111_4_alpha=-1*llr_max;
			
		//s111_5_alpha=log(exp(s111_4_alpha+gama4_00d) + exp(s001_4_alpha+gama4_01)); 
		tmp1=s111_4_alpha+gama4_00d;
		tmp2=s001_4_alpha+gama4_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s111_5_alpha)<llr_max2) s111_5_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s111_5_alpha>llr_max) 
			s111_5_alpha=llr_max;
		else if(s111_5_alpha<-1*llr_max) 
			s111_5_alpha=-1*llr_max;
			
		//s000_6_alpha=log(exp(s000_5_alpha+gama5_00) + exp(s111_5_alpha+gama5_01)); 
		tmp1=s000_5_alpha+gama5_00;
		tmp2=s111_5_alpha+gama5_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;	
		
		//s000_6_alpha_old=s000_6_alpha;			
		if(abs(s000_6_alpha)<llr_max2) s000_6_alpha= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_6_alpha>llr_max) 
			s000_6_alpha=llr_max;
		else if(s000_6_alpha<-1*llr_max) 
			s000_6_alpha=-1*llr_max;
		
		
		//cout<<"here";
		/*if(abs(tmp1-tmp2)>100) {	
			//cout<<"s000_6_alpha_old: "<<s000_6_alpha_old<<endl;	
			cout<<"num: "<<abs(tmp1-tmp2)<<endl;	
			cout<<" num2: "<<exp(-1*abs(tmp1-tmp2))<<endl;	
			cout<<"s000_6_alpha: "<<s000_6_alpha<<endl;	
		}*/
			
		/*s000_7_alpha=log(exp(s000_6_alpha+gama6_00) + exp(s101_6_alpha+gama6_01));   
		if(s000_7_alpha>llr_max) 
			s000_7_alpha=llr_max;
		else if(s000_7_alpha<-1*llr_max) 
			s000_7_alpha=-1*llr_max;*/
	
		
		/*cout<<"s000_1_alpha: "<<s000_1_alpha<<endl;
		cout<<"s010_1_alpha: "<<s010_1_alpha<<endl;
		cout<<endl;
		cout<<"s000_2_alpha: "<<s000_2_alpha<<endl;
		cout<<"s001_2_alpha: "<<s001_2_alpha<<endl;
		cout<<"s010_2_alpha: "<<s010_2_alpha<<endl;
		cout<<"s011_2_alpha: "<<s011_2_alpha<<endl;
		cout<<endl;
		cout<<"s000_3_alpha: "<<s000_3_alpha<<endl;
		cout<<"s001_3_alpha: "<<s001_3_alpha<<endl;
		cout<<"s010_3_alpha: "<<s010_3_alpha<<endl;
		cout<<"s011_3_alpha: "<<s011_3_alpha<<endl;
		cout<<endl;
		cout<<"s000_4_alpha: "<<s000_4_alpha<<endl;
		cout<<"s001_4_alpha: "<<s001_4_alpha<<endl;
		cout<<"s010_4_alpha: "<<s010_4_alpha<<endl;
		cout<<"s011_4_alpha: "<<s011_4_alpha<<endl;
		cout<<"s100_4_alpha: "<<s100_4_alpha<<endl;
		cout<<"s101_4_alpha: "<<s101_4_alpha<<endl;
		cout<<"s110_4_alpha: "<<s110_4_alpha<<endl;
		cout<<"s111_4_alpha: "<<s111_4_alpha<<endl;
		cout<<endl;
		cout<<"s000_5_alpha: "<<s000_5_alpha<<endl;
		cout<<"s010_5_alpha: "<<s010_5_alpha<<endl;
		cout<<"s101_5_alpha: "<<s101_5_alpha<<endl;
		cout<<"s111_5_alpha: "<<s111_5_alpha<<endl;
		cout<<endl;
		cout<<"s000_6_alpha: "<<s000_6_alpha<<endl;
		cout<<"s101_6_alpha: "<<s101_6_alpha<<endl;
		cout<<endl;
		cout<<"s000_7_alpha: "<<s000_7_alpha<<endl;
		cout<<endl;*/
	
	
		//backward metrics
		s000_7_beta=log(1);
		
		if(abs(s000_6_beta)<llr_max2) s000_6_beta=s000_7_beta+gama6_00;
		if(s000_6_beta>llr_max) 
			s000_6_beta=llr_max;
		else if(s000_6_beta<-1*llr_max) 
			s000_6_beta=-1*llr_max;
			
		if(abs(s000_5_beta)<llr_max2) s000_5_beta=s000_6_beta+gama5_00;
		if(s000_5_beta>llr_max) 
			s000_5_beta=llr_max;
		else if(s000_5_beta<-1*llr_max) 
			s000_5_beta=-1*llr_max;
			
		if(abs(s000_4_beta)<llr_max2) s000_4_beta=s000_5_beta+gama4_00;
		if(s000_4_beta>llr_max) 
			s000_4_beta=llr_max;
		else if(s000_4_beta<-1*llr_max) 
			s000_4_beta=-1*llr_max;
		
		if(abs(s101_6_beta)<llr_max2) s101_6_beta=s000_7_beta+gama6_01;
		if(s101_6_beta>llr_max) 
			s101_6_beta=llr_max;
		else if(s101_6_beta<-1*llr_max) 
			s101_6_beta=-1*llr_max;
			
		if(abs(s010_5_beta)<llr_max2) s010_5_beta=s101_6_beta+gama5_01b;
		if(s010_5_beta>llr_max) 
			s010_5_beta=llr_max;
		else if(s010_5_beta<-1*llr_max) 
			s010_5_beta=-1*llr_max;
			
		if(abs(s100_4_beta)<llr_max2) s100_4_beta=s010_5_beta+gama4_11b;
		if(s100_4_beta>llr_max) 
			s100_4_beta=llr_max;
		else if(s100_4_beta<-1*llr_max) 
			s100_4_beta=-1*llr_max;
			
		//s000_3_beta=log(exp(s000_4_beta+gama3_00) + exp(s100_4_beta+gama3_01));
		tmp1=s000_4_beta+gama3_00;
		tmp2=s100_4_beta+gama3_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s000_3_beta)<llr_max2) s000_3_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_3_beta>llr_max) 
			s000_3_beta=llr_max;
		else if(s000_3_beta<-1*llr_max) 
			s000_3_beta=-1*llr_max;
			
		if(abs(s101_5_beta)<llr_max2) s101_5_beta=s101_6_beta+gama5_00b;
		if(s101_5_beta>llr_max) 
			s101_5_beta=llr_max;
		else if(s101_5_beta<-1*llr_max) 
			s101_5_beta=-1*llr_max;
		
		if(abs(s111_5_beta)<llr_max2) s111_5_beta=s000_6_beta+gama5_01;
		if(s111_5_beta>llr_max) 
			s111_5_beta=llr_max;
		else if(s111_5_beta<-1*llr_max) 
			s111_5_beta=-1*llr_max;
		
		if(abs(s101_4_beta)<llr_max2) s101_4_beta=s101_5_beta+gama4_00c;
		if(s101_4_beta>llr_max) 
			s101_4_beta=llr_max;
		else if(s101_4_beta<-1*llr_max) 
			s101_4_beta=-1*llr_max;
			
		if(abs(s001_4_beta)<llr_max2) s001_4_beta=s111_5_beta+gama4_01;
		if(s001_4_beta>llr_max) 
			s001_4_beta=llr_max;
		else if(s001_4_beta<-1*llr_max) 
			s001_4_beta=-1*llr_max;
			
		//s001_3_beta=log(exp(s001_4_beta+gama3_00b) + exp(s101_4_beta+gama3_01b));
		tmp1=s001_4_beta+gama3_00b;
		tmp2=s101_4_beta+gama3_01b; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s001_3_beta)<llr_max2) s001_3_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s001_3_beta>llr_max) 
			s001_3_beta=llr_max;
		else if(s001_3_beta<-1*llr_max) 
			s001_3_beta=-1*llr_max;
			
		//s001_2_beta=log(exp(s001_3_beta+gama2_10) + exp(s000_3_beta+gama2_11));
		tmp1=s001_3_beta+gama2_10;
		tmp2=s000_3_beta+gama2_11; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s001_2_beta)<llr_max2) s001_2_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s001_2_beta>llr_max) 
			s001_2_beta=llr_max;
		else if(s001_2_beta<-1*llr_max) 
			s001_2_beta=-1*llr_max;
			
		//s000_2_beta=log(exp(s000_3_beta+gama2_00) + exp(s001_3_beta+gama2_01));
		tmp1=s000_3_beta+gama2_00;
		tmp2=s001_3_beta+gama2_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s000_2_beta)<llr_max2) s000_2_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_2_beta>llr_max) 
			s000_2_beta=llr_max;
		else if(s000_2_beta<-1*llr_max) 
			s000_2_beta=-1*llr_max;
		
		if(abs(s111_4_beta)<llr_max2) s111_4_beta=s111_5_beta+gama4_00d;
		if(s111_4_beta>llr_max) 
			s111_4_beta=llr_max;
		else if(s111_4_beta<-1*llr_max) 
			s111_4_beta=-1*llr_max;
		
		if(abs(s011_4_beta)<llr_max2) s011_4_beta=s101_5_beta+gama4_11c;
		if(s011_4_beta>llr_max) 
			s011_4_beta=llr_max;
		else if(s011_4_beta<-1*llr_max) 
			s011_4_beta=-1*llr_max;
			
		//s011_3_beta=log(exp(s011_4_beta+gama3_00d) + exp(s111_4_beta+gama3_01d));
		tmp1=s011_4_beta+gama3_00d;
		tmp2=s111_4_beta+gama3_01d; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s011_3_beta)<llr_max2) s011_3_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s011_3_beta>llr_max) 
			s011_3_beta=llr_max;
		else if(s011_3_beta<-1*llr_max) 
			s011_3_beta=-1*llr_max;
			
		if(abs(s011_2_beta)<llr_max2) s011_2_beta=s011_3_beta+gama2_10b;
		if(s011_2_beta>llr_max) 
			s011_2_beta=llr_max;
		else if(s011_2_beta<-1*llr_max) 
			s011_2_beta=-1*llr_max;
			
		//s000_1_beta=log(exp(s000_2_beta+gama1_00) + exp(s011_2_beta+gama1_01));
		tmp1=s000_2_beta+gama1_00;
		tmp2=s011_2_beta+gama1_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s000_1_beta)<llr_max2) s000_1_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_1_beta>llr_max) 
			s000_1_beta=llr_max;
		else if(s000_1_beta<-1*llr_max) 
			s000_1_beta=-1*llr_max;
		
		if(abs(s010_4_beta)<llr_max2) s010_4_beta=s010_5_beta+gama4_00b;
		if(s010_4_beta>llr_max) 
			s010_4_beta=llr_max;
		else if(s010_4_beta<-1*llr_max) 
			s010_4_beta=-1*llr_max;
			
		if(abs(s110_4_beta)<llr_max2) s110_4_beta=s000_5_beta+gama4_11;
		if(s110_4_beta>llr_max) 
			s110_4_beta=llr_max;
		else if(s110_4_beta<-1*llr_max) 
			s110_4_beta=-1*llr_max;
			
		//s010_3_beta=log(exp(s010_4_beta+gama3_00c) + exp(s110_4_beta+gama3_01c));
		tmp1=s010_4_beta+gama3_00c;
		tmp2=s110_4_beta+gama3_01c; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s010_3_beta)<llr_max2) s010_3_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s010_3_beta>llr_max) 
			s010_3_beta=llr_max;
		else if(s010_3_beta<-1*llr_max) 
			s010_3_beta=-1*llr_max;
			
		//s010_2_beta=log(exp(s010_3_beta+gama2_00b) + exp(s011_3_beta+gama2_01b));
		tmp1=s010_3_beta+gama2_00b;
		tmp2=s011_3_beta+gama2_01b; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s010_2_beta)<llr_max2) s010_2_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s010_2_beta>llr_max) 
			s010_2_beta=llr_max;
		else if(s010_2_beta<-1*llr_max) 
			s010_2_beta=-1*llr_max;
			
		//s010_1_beta=log(exp(s010_2_beta+gama1_10) + exp(s001_2_beta+gama1_11));
		tmp1=s010_2_beta+gama1_10;
		tmp2=s001_2_beta+gama1_11; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s010_1_beta)<llr_max2) s010_1_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s010_1_beta>llr_max) 
			s010_1_beta=llr_max;
		else if(s010_1_beta<-1*llr_max) 
			s010_1_beta=-1*llr_max;
			
		//s000_0_beta=log(exp(s000_1_beta+gama0_00) + exp(s010_1_beta+gama0_01));
		tmp1=s000_1_beta+gama0_00;
		tmp2=s010_1_beta+gama0_01; 
		if(tmp1>tmp2)
			max1=tmp1;
		else
			max1=tmp2;		
		if(abs(s000_0_beta)<llr_max2) s000_0_beta= max1+log(1+exp(-1*abs(tmp1-tmp2)));
		if(s000_0_beta>llr_max) 
			s000_0_beta=llr_max;
		else if(s000_0_beta<-1*llr_max) 
			s000_0_beta=-1*llr_max;
		
		/*cout<<"s010_1_beta: "<<s010_1_beta<<endl;
		cout<<endl;
		cout<<"s000_2_beta: "<<s000_2_beta<<endl;
		cout<<"s001_2_beta: "<<s001_2_beta<<endl;
		cout<<"s010_2_beta: "<<s010_2_beta<<endl;
		cout<<"s011_2_beta: "<<s011_2_beta<<endl;
		cout<<endl;
		cout<<"s000_3_beta: "<<s000_3_beta<<endl;
		cout<<"s001_3_beta: "<<s001_3_beta<<endl;
		cout<<"s010_3_beta: "<<s010_3_beta<<endl;
		cout<<"s011_3_beta: "<<s011_3_beta<<endl;
		cout<<endl;
		cout<<"s000_4_beta: "<<s000_4_beta<<endl;
		cout<<"s001_4_beta: "<<s001_4_beta<<endl;
		cout<<"s010_4_beta: "<<s010_4_beta<<endl;
		cout<<"s011_4_beta: "<<s011_4_beta<<endl;
		cout<<"s100_4_beta: "<<s100_4_beta<<endl;
		cout<<"s101_4_beta: "<<s101_4_beta<<endl;
		cout<<"s110_4_beta: "<<s110_4_beta<<endl;
		cout<<"s111_4_beta: "<<s111_4_beta<<endl;
		cout<<endl;
		cout<<"s000_5_beta: "<<s000_5_beta<<endl;
		cout<<"s010_5_beta: "<<s010_5_beta<<endl;
		cout<<"s101_5_beta: "<<s101_5_beta<<endl;
		cout<<"s111_5_beta: "<<s111_5_beta<<endl;
		cout<<endl;
		cout<<"s000_6_beta: "<<s000_6_beta<<endl;
		cout<<"s101_6_beta: "<<s101_6_beta<<endl;
		cout<<endl;
		cout<<"s000_7_beta: "<<s000_7_beta<<endl;
		cout<<endl;*/
		
		
	//computing LLRs at each stage (messages sent by this CN to its neighboring VNs)

	//cout<<"L_out: "; 
	for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
		i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'
		if(!j) {
			//L_out[i]= log(exp(s000_0_alpha+gama0_00+s000_1_beta)) - log(exp(s000_0_alpha+gama0_01+s010_1_beta));
			tmp1=s000_0_alpha+gama0_00+s000_1_beta;
			tmp2=s000_0_alpha+gama0_01+s010_1_beta;
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
			L_out[i]= log(exp(tmp1)) - log(exp(tmp2)); 
			
			//cout<<"j: "<<j<<'\n';
			//cout<<exp(tmp1)<<endl;
			//cout<<exp(tmp2)<<endl;
		}
		
		else if(j==1) {
			//L_out[i]= log(exp(s000_1_alpha+gama1_00+s000_2_beta)+exp(s010_1_alpha+gama1_10+s010_2_beta)) - log(exp(s000_1_alpha+gama1_01+s011_2_beta)+exp(s010_1_alpha+gama1_11+s001_2_beta));
			
			tmp1=s000_1_alpha+gama1_00+s000_2_beta;
			tmp2=s010_1_alpha+gama1_10+s010_2_beta;
			tmp1b=s000_1_alpha+gama1_01+s011_2_beta;
			tmp2b=s010_1_alpha+gama1_11+s001_2_beta;
			
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
			if(tmp1b>llr_max) 
				tmp1b=llr_max;
			else if(tmp1b<-1*llr_max) 
				tmp1b=-1*llr_max;
			if(tmp2b>llr_max) 
				tmp2b=llr_max;
			else if(tmp2b<-1*llr_max) 
				tmp2b=-1*llr_max;
					
			if(tmp1>tmp2)
				max1=tmp1;
			else
				max1=tmp2;
			if(tmp1b>tmp2b)
				max2=tmp1b;
			else
				max2=tmp2b;
					
			L_out[i]= max1+log(1+exp(-1*abs(tmp1-tmp2))) - (max2+log(1+exp(-1*abs(tmp1b-tmp2b))));
			//cout<<"j: "<<j<<'\n';
			//cout<<exp(tmp1)<<endl;
			//cout<<exp(tmp2)<<endl;
			//cout<<exp(tmp3)<<endl;
			//cout<<exp(tmp4)<<endl;
		}
	
		else if(j==2) {
			//L_out[i]= log(exp(s000_2_alpha+gama2_00+s000_3_beta)+exp(s001_2_alpha+gama2_10+s001_3_beta)+exp(s010_2_alpha+gama2_00b+s010_3_beta)+exp(s011_2_alpha+gama2_10b+s011_3_beta)) - log(exp(s000_2_alpha+gama2_01+s001_3_beta)+exp(s001_2_alpha+gama2_11+s000_3_beta)+exp(s010_2_alpha+gama2_01b+s011_3_beta)+exp(s011_2_alpha+gama2_11b+s010_3_beta));
			
			tmp1=s000_2_alpha+gama2_00+s000_3_beta;
			tmp2=s001_2_alpha+gama2_10+s001_3_beta;
			tmp3=s010_2_alpha+gama2_00b+s010_3_beta;
			tmp4=s011_2_alpha+gama2_10b+s011_3_beta;
			
			tmp1b=s000_2_alpha+gama2_01+s001_3_beta;
			tmp2b=s001_2_alpha+gama2_11+s000_3_beta;
			tmp3b=s010_2_alpha+gama2_01b+s011_3_beta;
			tmp4b=s011_2_alpha+gama2_11b+s010_3_beta;
			
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
			if(tmp3>llr_max) 
				tmp3=llr_max;
			else if(tmp3<-1*llr_max) 
				tmp3=-1*llr_max;
			if(tmp4>llr_max) 
				tmp4=llr_max;
			else if(tmp4<-1*llr_max) 
				tmp4=-1*llr_max;
			if(tmp1b>llr_max) 
				tmp1b=llr_max;
			else if(tmp1b<-1*llr_max) 
				tmp1b=-1*llr_max;
			if(tmp2b>llr_max) 
				tmp2b=llr_max;
			else if(tmp2b<-1*llr_max) 
				tmp2b=-1*llr_max;
			if(tmp3b>llr_max) 
				tmp3b=llr_max;
			else if(tmp3b<-1*llr_max) 
				tmp3b=-1*llr_max;
			if(tmp4b>llr_max) 
				tmp4b=llr_max;
			else if(tmp4b<-1*llr_max) 
				tmp4b=-1*llr_max;
			//L_out[i]= log(exp(tmp1)+exp(tmp2)+exp(tmp3)+exp(tmp4)) - log(exp(tmp1b)+exp(tmp2b)+exp(tmp3b)+exp(tmp4b));
			//max(a,b,c,d)=max((a,b,c),d)=max(max2(max1(a=b),c)),d)=max(max2,d)=max(max2,tmp4)=max3+log(1+exp(-1*abs(max2-tmp4)))
					
			if(tmp1>tmp2)
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
				max3b=tmp4b;
					
			L_out[i]=max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b))));
		}
	
		else if(j==3) {
			//L_out[i]= log(exp(s000_3_alpha+gama3_00+s000_4_beta)+exp(s001_3_alpha+gama3_00b+s001_4_beta)+exp(s010_3_alpha+gama3_00c+s010_4_beta)+exp(s011_3_alpha+gama3_00d+s011_4_beta)) - log(exp(s000_3_alpha+gama3_01+s100_4_beta)+exp(s001_3_alpha+gama3_01b+s101_4_beta)+exp(s010_3_alpha+gama3_01c+s110_4_beta)+exp(s011_3_alpha+gama3_01d+s111_4_beta));
			
			tmp1=s000_3_alpha+gama3_00+s000_4_beta;
			tmp2=s001_3_alpha+gama3_00b+s001_4_beta; //s010_3_alpha+gama3_00c+s010_4_beta;
			tmp3=s010_3_alpha+gama3_00c+s010_4_beta; //s010_2_alpha+gama2_00b+s010_3_beta;
			tmp4=s011_3_alpha+gama3_00d+s011_4_beta;
			
			tmp1b=s000_3_alpha+gama3_01+s100_4_beta;
			tmp2b=s001_3_alpha+gama3_01b+s101_4_beta;
			tmp3b=s010_3_alpha+gama3_01c+s110_4_beta;
			tmp4b=s011_3_alpha+gama3_01d+s111_4_beta;
			
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
			if(tmp3>llr_max) 
				tmp3=llr_max;
			else if(tmp3<-1*llr_max) 
				tmp3=-1*llr_max;
			if(tmp4>llr_max) 
				tmp4=llr_max;
			else if(tmp4<-1*llr_max) 
				tmp4=-1*llr_max;
			if(tmp1b>llr_max) 
				tmp1b=llr_max;
			else if(tmp1b<-1*llr_max) 
				tmp1b=-1*llr_max;
			if(tmp2b>llr_max) 
				tmp2b=llr_max;
			else if(tmp2b<-1*llr_max) 
				tmp2b=-1*llr_max;
			if(tmp3b>llr_max) 
				tmp3b=llr_max;
			else if(tmp3b<-1*llr_max) 
				tmp3b=-1*llr_max;
			if(tmp4b>llr_max) 
				tmp4b=llr_max;
			else if(tmp4b<-1*llr_max) 
				tmp4b=-1*llr_max;
										
			if(tmp1>tmp2)
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
				max3b=tmp4b;
					
			L_out[i]=max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b))));
		}
	
		else if(j==4) {
			//L_out[i]= log(exp(s000_4_alpha+gama4_00+s000_5_beta)+exp(s010_4_alpha+gama4_00b+s010_5_beta)+exp(s101_4_alpha+gama4_00c+s101_5_beta)+exp(s111_4_alpha+gama4_00d+s111_5_beta)) - log(exp(s001_4_alpha+gama4_01+s111_5_beta)+exp(s011_4_alpha+gama4_11c+s101_5_beta)+exp(s100_4_alpha+gama4_11b+s010_5_beta)+exp(s110_4_alpha+gama4_11+s000_5_beta));
			
			tmp1=s000_4_alpha+gama4_00+s000_5_beta;
			tmp2=s010_4_alpha+gama4_00b+s010_5_beta;
			tmp3=s101_4_alpha+gama4_00c+s101_5_beta;
			tmp4=s111_4_alpha+gama4_00d+s111_5_beta;
			
			tmp1b=s001_4_alpha+gama4_01+s111_5_beta;
			tmp2b=s011_4_alpha+gama4_11c+s101_5_beta;
			tmp3b=s100_4_alpha+gama4_11b+s010_5_beta;
			tmp4b=s110_4_alpha+gama4_11+s000_5_beta;
			
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
			if(tmp3>llr_max) 
				tmp3=llr_max;
			else if(tmp3<-1*llr_max) 
				tmp3=-1*llr_max;
			if(tmp4>llr_max) 
				tmp4=llr_max;
			else if(tmp4<-1*llr_max) 
				tmp4=-1*llr_max;
			if(tmp1b>llr_max) 
				tmp1b=llr_max;
			else if(tmp1b<-1*llr_max) 
				tmp1b=-1*llr_max;
			if(tmp2b>llr_max) 
				tmp2b=llr_max;
			else if(tmp2b<-1*llr_max) 
				tmp2b=-1*llr_max;
			if(tmp3b>llr_max) 
				tmp3b=llr_max;
			else if(tmp3b<-1*llr_max) 
				tmp3b=-1*llr_max;
			if(tmp4b>llr_max) 
				tmp4b=llr_max;
			else if(tmp4b<-1*llr_max) 
				tmp4b=-1*llr_max;
			//L_out[i]= log(exp(tmp1)+exp(tmp2)+exp(tmp3)+exp(tmp4)) - log(exp(tmp1b)+exp(tmp2b)+exp(tmp3b)+exp(tmp4b));
			
			if(tmp1>tmp2)
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
				max3b=tmp4b;
					
			L_out[i]=max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b))));
		}
	
		else if(j==5) {
			//L_out[i]= log(exp(s000_5_alpha+gama5_00+s000_6_beta)+exp(s101_5_alpha+gama5_00b+s101_6_beta)) - log(exp(s010_5_alpha+gama5_01b+s101_6_beta)+exp(s111_5_alpha+gama5_01+s000_6_beta));
			
			tmp1=s000_5_alpha+gama5_00+s000_6_beta;
			tmp2=s101_5_alpha+gama5_00b+s101_6_beta;
			tmp1b=s010_5_alpha+gama5_01b+s101_6_beta;
			tmp2b=s111_5_alpha+gama5_01+s000_6_beta;
			
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
			if(tmp1b>llr_max) 
				tmp1b=llr_max;
			else if(tmp1b<-1*llr_max) 
				tmp1b=-1*llr_max;
			if(tmp2b>llr_max) 
				tmp2b=llr_max;
			else if(tmp2b<-1*llr_max) 
				tmp2b=-1*llr_max;
								
			if(tmp1>tmp2)
				max1=tmp1;
			else
				max1=tmp2;
			if(tmp1b>tmp2b)
				max2=tmp1b;
			else
				max2=tmp2b;
					
			L_out[i]= max1+log(1+exp(-1*abs(tmp1-tmp2))) - (max2+log(1+exp(-1*abs(tmp1b-tmp2b))));
		}

		else if(j==6) {
			//L_out[i]= log(exp(s000_6_alpha+gama6_00+s000_7_beta)) - log(exp(s101_6_alpha+gama6_01+s000_7_beta));
			
			tmp1=s000_6_alpha+gama6_00+s000_7_beta;
			tmp2=s101_6_alpha+gama6_01+s000_7_beta;
			
			if(tmp1>llr_max) 
				tmp1=llr_max;
			else if(tmp1<-1*llr_max) 
				tmp1=-1*llr_max;
			if(tmp2>llr_max) 
				tmp2=llr_max;
			else if(tmp2<-1*llr_max) 
				tmp2=-1*llr_max;
				
			L_out[i]= log(exp(tmp1)) - log(exp(tmp2));
			//cout<<"j: "<<j<<'\n';
			//cout<<exp(tmp1)<<endl;
			//cout<<exp(tmp2)<<endl;
		}
		//cout<<L_out[i]<< " ";
		
		if(L_out[i]>llr_max) 
			L_out[i]=llr_max;
		else if(L_out[i]<-1*llr_max) 
			L_out[i]=-1*llr_max;
			
		//if(isnan(L_out[i]))
			//L_out[i]=0;
	}
	//cout<<"a: "<<a<<endl;
	//cout<<"L_out: "<<endl;
	for(i=0;i<n;i++)
		if(isnan(L_out[i])) 
			cout<<L_out[i]<<" ";
	
	
}

