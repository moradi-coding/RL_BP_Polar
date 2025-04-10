
//BCJR calclations using the [7,4] Hamming trellis in Vivian paper, calculations based on Ryan-Lin book
//flooding GLDPC decoder for all GCNs
void bcjr(int a) { 
		
	long i,j;
	double b,c,d=0.5,y_sum,sum_alpha1,sum_alpha2,sum_alpha3,sum_alpha4,sum_alpha5,sum_alpha6,sum_alpha7,
	sum_beta1,sum_beta2,sum_beta3,sum_beta4,sum_beta5,sum_beta6,sum_beta7,prob1,prob0,prob0a,prob0b,prob0c,prob0d,prob1a,prob1b,prob1c,prob1d,sum;
	
	c=-1/(2*varn);
	//for(a=0;a<num_cls;a*+) {
	//initialization of branch metrics
	/*y_sum=0;
	for(j=0;j<num_vns_cls;j++) {
		//i=vns_cluster[a][j];
		i=j;
		y_sum+=y[i];
	}*/
		
	
		for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
			i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'
			//else L_in[i]=L_out[i]-LR[i];
			
			if(L_in[i]>llr_max) 
				L_in[i]=llr_max;
			else if(L_in[i]<-1*llr_max) 
				L_in[i]=-1*llr_max;
				
			//b=exp(-0.5*L_in[i])/(1+exp(-1*L_in[i]));
			b=L_in[i]-LR[i];
			//cout<<"j, b, LLR: "<<j<<", "<<b<<", "<<L_in[i]<<endl;
			
			if(!j) {
				gama0_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); //0 mapped to 1	
				gama0_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); //1 mapped to -1	
			}	
			else if(j==1) {
				gama1_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama1_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); 
				gama1_10=gama1_00;
				gama1_11=gama1_01;
			}	
			else if(j==2) {
				gama2_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama2_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); 
				gama2_10=gama2_00;
				gama2_00b=gama2_00;
				gama2_10b=gama2_00;
				gama2_11=gama2_01;
				gama2_01b=gama2_01;
				gama2_11b=gama2_01;
			}	
			else if(j==3) {
				gama3_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama3_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); 
				gama3_00b=gama3_00;
				gama3_00c=gama3_00;
				gama3_00d=gama3_00;
				gama3_01b=gama3_01;
				gama3_01c=gama3_01;
				gama3_01d=gama3_01;
			}	
			else if(j==4) {
				gama4_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama4_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); 
				gama4_00b=gama4_00;
				gama4_00c=gama4_00;
				gama4_00d=gama4_00;
				gama4_11b=gama4_01;
				gama4_11c=gama4_01;
				gama4_11=gama4_01;
			}	
			else if(j==5) {
				gama5_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama5_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); 
				gama5_00b=gama5_00;
				gama5_01b=gama5_01;
			}	
			else if(j==6) {
				gama6_00=exp(c*pow(abs(y[i]-1),2)+d*1*b); //0 mapped to 1	
				gama6_01=exp(c*pow(abs(y[i]+1),2)+d*-1*b); //1 mapped to -1	
			}	
		}
		//cout<<"L_in: "; for(i=0;i<n;i++) cout<<L_in[i]<< " "; cout<<'\n'<<'\n';
		//cout<<"gama1_00: "<<gama1_00<<endl;
		
		//forward metrics
		s000_0_alpha=1;
		s000_1_alpha=s000_0_alpha*gama0_00;
		
		s000_2_alpha=s000_1_alpha*gama1_00;
		
		s010_1_alpha=s000_0_alpha*gama0_01;
		
		s001_2_alpha=s010_1_alpha*gama1_11;
		s000_3_alpha=((s000_2_alpha*gama2_00) + (s001_2_alpha*gama2_11));
		s000_4_alpha=s000_3_alpha*gama3_00;
		
		s010_2_alpha=s010_1_alpha*gama1_10;
		s011_2_alpha=s000_1_alpha*gama1_01;
		s010_3_alpha=((s010_2_alpha*gama2_00b) + (s011_2_alpha*gama2_11b));
		s110_4_alpha=s010_3_alpha*gama3_01c;
		
		s000_5_alpha=((s000_4_alpha*gama4_00) + (s110_4_alpha*gama4_11));
		
		s001_3_alpha=((s001_2_alpha*gama2_10) + (s000_2_alpha*gama2_01));
		s101_4_alpha=s001_3_alpha*gama3_01b;
		s011_3_alpha=((s011_2_alpha*gama2_10b) + (s010_2_alpha*gama2_01b));
		s011_4_alpha=s011_3_alpha*gama3_00d;
		s101_5_alpha=((s101_4_alpha*gama4_00c) + (s011_4_alpha*gama4_11c));
		s010_4_alpha=s010_3_alpha*gama3_00c;	
		
		s100_4_alpha=s000_3_alpha*gama3_01;
		s010_5_alpha=((s010_4_alpha*gama4_00b) + (s100_4_alpha*gama4_11b));
		s101_6_alpha=((s101_5_alpha*gama5_00b) + (s010_5_alpha*gama5_01b));
		
		s001_4_alpha=s001_3_alpha*gama3_00b;
		s111_4_alpha=s011_3_alpha*gama3_01d;
		s111_5_alpha=((s111_4_alpha*gama4_00d) + (s001_4_alpha*gama4_01)); 
		s000_6_alpha=((s000_5_alpha*gama5_00) + (s111_5_alpha*gama5_01)); 
		//s000_7_alpha=((s000_6_alpha*gama6_00) + (s101_6_alpha*gama6_01)); 
		
		
		
		//normalizing alphas
		sum_alpha1=s000_1_alpha+s010_1_alpha;
		sum_alpha2=s000_2_alpha+s001_2_alpha+s010_2_alpha+s011_2_alpha;
		sum_alpha3=s000_3_alpha+s001_3_alpha+s010_3_alpha+s011_3_alpha;
		sum_alpha4=s000_4_alpha+s001_4_alpha+s010_4_alpha+s011_4_alpha+s100_4_alpha+s101_4_alpha+s110_4_alpha+s111_4_alpha;
		sum_alpha5=s000_5_alpha+s010_5_alpha+s101_5_alpha+s111_5_alpha;
		sum_alpha6=s000_6_alpha+s101_6_alpha;
		sum_alpha7=s000_7_alpha;
		
		s000_1_alpha/=sum_alpha1;
		s010_1_alpha/=sum_alpha1;
		
		s000_2_alpha/=sum_alpha2;
		s001_2_alpha/=sum_alpha2;
		s010_2_alpha/=sum_alpha2;
		s011_2_alpha/=sum_alpha2;
		
		s000_3_alpha/=sum_alpha3;
		s001_3_alpha/=sum_alpha3;
		s010_3_alpha/=sum_alpha3;
		s011_3_alpha/=sum_alpha3;
		
		s000_4_alpha/=sum_alpha4;
		s001_4_alpha/=sum_alpha4;
		s010_4_alpha/=sum_alpha4;
		s011_4_alpha/=sum_alpha4;
		s100_4_alpha/=sum_alpha4;
		s101_4_alpha/=sum_alpha4;
		s110_4_alpha/=sum_alpha4;
		s111_4_alpha/=sum_alpha4;
		
		s000_5_alpha/=sum_alpha5;
		s010_5_alpha/=sum_alpha5;
		s101_5_alpha/=sum_alpha5;
		s111_5_alpha/=sum_alpha5;
		
		s000_6_alpha/=sum_alpha6;
		s101_6_alpha/=sum_alpha6;
		
		s000_7_alpha/=sum_alpha7;	
		
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
		cout<<endl;*/
	
	
		//backward metrics
		s000_7_beta=1;
		s000_6_beta=s000_7_beta*gama6_00;
		s000_5_beta=s000_6_beta*gama5_00;
		s000_4_beta=s000_5_beta*gama4_00;
		
		s101_6_beta=s000_7_beta*gama6_01;
		
		s010_5_beta=s101_6_beta*gama5_01b;
		s100_4_beta=s010_5_beta*gama4_11b;
		s000_3_beta=((s000_4_beta*gama3_00) + (s100_4_beta*gama3_01));
		
		s101_5_beta=s101_6_beta*gama5_00b;
		
		s111_5_beta=s000_6_beta*gama5_01;
		
		s101_4_beta=s101_5_beta*gama4_00c;
		s001_4_beta=s111_5_beta*gama4_01;
		s001_3_beta=((s001_4_beta*gama3_00b) + (s101_4_beta*gama3_01b));
		
		s001_2_beta=((s001_3_beta*gama2_10) + (s000_3_beta*gama2_11));
		
		s000_2_beta=((s000_3_beta*gama2_00) + (s001_3_beta*gama2_01));
		
		s111_4_beta=s111_5_beta*gama4_00d;
		
		s011_4_beta=s101_5_beta*gama4_11c;
		s011_3_beta=((s011_4_beta*gama3_00d) + (s111_4_beta*gama3_01d));
		s011_2_beta=s011_3_beta*gama2_10b;
		s000_1_beta=((s000_2_beta*gama1_00) + (s011_2_beta*gama1_01));
		
		s010_4_beta=s010_5_beta*gama4_00b;
		s110_4_beta=s000_5_beta*gama4_11;
		s010_3_beta=((s010_4_beta*gama3_00c) + (s110_4_beta*gama3_01c));
		s010_2_beta=((s010_3_beta*gama2_00b) + (s011_3_beta*gama2_01b));
		s010_1_beta=((s010_2_beta*gama1_10) + (s001_2_beta*gama1_11));
		s000_0_beta=((s000_1_beta*gama0_00) + (s010_1_beta*gama0_01));
		
		
		//normalizing betas
		sum_beta1=s000_1_beta+s010_1_beta;
		sum_beta2=s000_2_beta+s001_2_beta+s010_2_beta+s011_2_beta;
		sum_beta3=s000_3_beta+s001_3_beta+s010_3_beta+s011_3_beta;
		sum_beta4=s000_4_beta+s001_4_beta+s010_4_beta+s011_4_beta+s100_4_beta+s101_4_beta+s110_4_beta+s111_4_beta;
		sum_beta5=s000_5_beta+s010_5_beta+s101_5_beta+s111_5_beta;
		sum_beta6=s000_6_beta+s101_6_beta;
		sum_beta7=s000_7_beta;	
		
		s000_1_beta/=sum_beta1;
		s010_1_beta/=sum_beta1;
		
		s000_2_beta/=sum_beta2;
		s001_2_beta/=sum_beta2;
		s010_2_beta/=sum_beta2;
		s011_2_beta/=sum_beta2;
		
		s000_3_beta/=sum_beta3;
		s001_3_beta/=sum_beta3;
		s010_3_beta/=sum_beta3;
		s011_3_beta/=sum_beta3;
		
		s000_4_beta/=sum_beta4;
		s001_4_beta/=sum_beta4;
		s010_4_beta/=sum_beta4;
		s011_4_beta/=sum_beta4;
		s100_4_beta/=sum_beta4;
		s101_4_beta/=sum_beta4;
		s110_4_beta/=sum_beta4;
		s111_4_beta/=sum_beta4;
		
		s000_5_beta/=sum_beta5;
		s010_5_beta/=sum_beta5;
		s101_5_beta/=sum_beta5;
		s111_5_beta/=sum_beta5;
		
		s000_6_beta/=sum_beta6;
		s101_6_beta/=sum_beta6;
		
		s000_7_beta/=sum_beta7;
		
		/*cout<<"s000_1_beta: "<<s000_1_beta<<endl;
		cout<<"s010_1_beta: "<<s010_1_beta<<endl;
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
			prob0=s000_0_alpha*gama0_00*s000_1_beta;
			prob1=s000_0_alpha*gama0_01*s010_1_beta;
			sum=prob1+prob0;
			
			prob1/=sum;
			prob0/=sum;
			L_out[i]= log(prob0/prob1);
			/*cout<<"prob0: "<<prob0<<'\n';
			cout<<"prob1: "<<prob1<<'\n';
			cout<<endl;*/
			//cout<<"sum of prob: "<<prob0+prob1<<'\n';
		}
		else if(j==1) {
			prob0a=s000_1_alpha*gama1_00*s000_2_beta;
			prob0b=s010_1_alpha*gama1_10*s010_2_beta;
			prob1a=s000_1_alpha*gama1_01*s011_2_beta;
			prob1b=s010_1_alpha*gama1_11*s001_2_beta;
			sum=prob1a+prob1b+prob0a+prob0b;
			
			prob1a/=sum;
			prob1b/=sum;
			prob0a/=sum;
			prob0b/=sum;
			/*cout<<"prob0a: "<<prob0a<<'\n';
			cout<<"prob0b: "<<prob0b<<'\n';
			cout<<"prob1a: "<<prob1a<<'\n';
			cout<<"prob1b: "<<prob1b<<'\n';
			cout<<endl;*/
			
			L_out[i]= log((prob0a+prob0b)/(prob1a+prob1b));
		}
	
		else if(j==2) {
			prob0a=s000_2_alpha*gama2_00*s000_3_beta;
			prob0b=s001_2_alpha*gama2_10*s001_3_beta;
			prob0c=s010_2_alpha*gama2_00b*s010_3_beta;
			prob0d=s011_2_alpha*gama2_10b*s011_3_beta;
			
			prob1a=s000_2_alpha*gama2_01*s001_3_beta;
			prob1b=s001_2_alpha*gama2_11*s000_3_beta;
			prob1c=s010_2_alpha*gama2_01b*s011_3_beta;
			prob1d=s011_2_alpha*gama2_11b*s010_3_beta;
			
			sum=prob1a+prob1b+prob1c+prob1d+prob0a+prob0b+prob0c+prob0d;
			
			prob1a/=sum;
			prob1b/=sum;
			prob1c/=sum;
			prob1d/=sum;
			prob0a/=sum;
			prob0b/=sum;
			prob0c/=sum;
			prob0d/=sum;
			
			/*cout<<"prob0a: "<<prob0a<<'\n';
			cout<<"prob0b: "<<prob0b<<'\n';
			cout<<"prob0c: "<<prob0c<<'\n';
			cout<<"prob0d: "<<prob0d<<'\n';
			cout<<"prob1a: "<<prob1a<<'\n';
			cout<<"prob1b: "<<prob1b<<'\n';
			cout<<"prob1c: "<<prob1c<<'\n';
			cout<<"prob1d: "<<prob1d<<'\n';
			cout<<endl;*/
			
			L_out[i]= log(prob0a+prob0b+prob0c+prob0d) - log(prob1a+prob1b+prob1c+prob1d);
		}
	
		else if(j==3) {
			prob0a=s000_3_alpha*gama3_00*s000_4_beta;
			prob0b=s001_3_alpha*gama3_00b*s001_4_beta;
			prob0c=s010_3_alpha*gama3_00c*s010_4_beta;
			prob0d=s011_3_alpha*gama3_00d*s011_4_beta;
			
			prob1a=s000_3_alpha*gama3_01*s100_4_beta;
			prob1b=s001_3_alpha*gama3_01b*s101_4_beta;
			prob1c=s010_3_alpha*gama3_01c*s110_4_beta;
			prob1d=s011_3_alpha*gama3_01d*s111_4_beta;
			
			sum=prob1a+prob1b+prob1c+prob1d+prob0a+prob0b+prob0c+prob0d;
			
			prob1a/=sum;
			prob1b/=sum;
			prob1c/=sum;
			prob1d/=sum;
			prob0a/=sum;
			prob0b/=sum;
			prob0c/=sum;
			prob0d/=sum;
			
			/*cout<<"prob0a: "<<prob0a<<'\n';
			cout<<"prob0b: "<<prob0b<<'\n';
			cout<<"prob0c: "<<prob0c<<'\n';
			cout<<"prob0d: "<<prob0d<<'\n';
			cout<<"prob1a: "<<prob1a<<'\n';
			cout<<"prob1b: "<<prob1b<<'\n';
			cout<<"prob1c: "<<prob1c<<'\n';
			cout<<"prob1d: "<<prob1d<<'\n';
			cout<<endl;*/
			
			L_out[i]= log(prob0a+prob0b+prob0c+prob0d) - log(prob1a+prob1b+prob1c+prob1d);
		}
	
		else if(j==4) {		
			prob0a=s000_4_alpha*gama4_00*s000_5_beta;
			prob0b=s010_4_alpha*gama4_00b*s010_5_beta;
			prob0c=s101_4_alpha*gama4_00c*s101_5_beta;
			prob0d=s111_4_alpha*gama4_00d*s111_5_beta;
			
			prob1a=s001_4_alpha*gama4_01*s111_5_beta;
			prob1b=s011_4_alpha*gama4_11c*s101_5_beta;
			prob1c=s100_4_alpha*gama4_11b*s010_5_beta;
			prob1d=s110_4_alpha*gama4_11*s000_5_beta;
			
			sum=prob1a+prob1b+prob1c+prob1d+prob0a+prob0b+prob0c+prob0d;
			
			prob1a/=sum;
			prob1b/=sum;
			prob1c/=sum;
			prob1d/=sum;
			prob0a/=sum;
			prob0b/=sum;
			prob0c/=sum;
			prob0d/=sum;
			
			/*cout<<"prob0a: "<<prob0a<<'\n';
			cout<<"prob0b: "<<prob0b<<'\n';
			cout<<"prob0c: "<<prob0c<<'\n';
			cout<<"prob0d: "<<prob0d<<'\n';
			cout<<"prob1a: "<<prob1a<<'\n';
			cout<<"prob1b: "<<prob1b<<'\n';
			cout<<"prob1c: "<<prob1c<<'\n';
			cout<<"prob1d: "<<prob1d<<'\n';
			cout<<endl;*/
			
			L_out[i]= log(prob0a+prob0b+prob0c+prob0d) - log(prob1a+prob1b+prob1c+prob1d);
		}
	
		else if(j==5) {
			prob0a=s000_5_alpha*gama5_00*s000_6_beta;
			prob0b=s101_5_alpha*gama5_00b*s101_6_beta;
			prob1a=s010_5_alpha*gama5_01b*s101_6_beta;
			prob1b=s111_5_alpha*gama5_01*s000_6_beta;
			sum=prob1a+prob1b+prob0a+prob0b;
			
			prob1a/=sum;
			prob1b/=sum;
			prob0a/=sum;
			prob0b/=sum;
			
			/*cout<<"prob0a: "<<prob0a<<'\n';
			cout<<"prob0b: "<<prob0b<<'\n';
			cout<<"prob1a: "<<prob1a<<'\n';
			cout<<"prob1b: "<<prob1b<<'\n';
			cout<<endl;*/
			
			L_out[i]= log((prob0a+prob0b)/(prob1a+prob1b));
		}
	
		else if(j==6) {
			prob0=s000_6_alpha*gama6_00*s000_7_beta;
			prob1=s101_6_alpha*gama6_01*s000_7_beta;
			sum=prob1+prob0;
			
			prob1/=sum;
			prob0/=sum;
			
			/*cout<<"prob0: "<<prob0<<'\n';
			cout<<"prob1: "<<prob1<<'\n';
			cout<<endl;*/
			
			L_out[i]= log(prob0/prob1);
			//cout<<"sum of prob: "<<prob0+prob1<<'\n';
		}
		if(L_out[i]>llr_max) 
			L_out[i]=llr_max;
		else if(L_out[i]<-1*llr_max) 
			L_out[i]=-1*llr_max;
		//cout<<L_out[i]<<'\n';
	}
	//cout<<"L_out: "; for(i=0;i<n;i++) cout<<L_out[i]<< " "; cout<<'\n'<<'\n';
	for(i=0;i<n;i++)
		if(isnan(L_out[i])) 
			cout<<L_out[i]<<" ";
	
	//}
	
	
}

