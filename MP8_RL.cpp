
// sequential cluster scheduling LDPC decoding
// only one cluster is decoded per iteration

#include "th.cpp"
#include "bcjr_log_15_11.cpp"
#include "bcjr_log2.cpp"
#include "clus_decode.cpp"
#include "clus_decode_bcjr.cpp"

void MP8_RL(int cnt, int a, int l) { 
	long cw,k,flg4;
	
	for(cw=0;cw<CW;cw++) {
		flg4=1;
		for(k=0;k<cnt;k++) 
			if(cw==excl_cw[k]) {flg4=0; break;} //if flg4=0, cw has been recovered already
		
		if(flg4) {
			//cout<<"a: "<<a<<endl;
			if(!bcjr) clus_decode(a);
			else clus_decode_bcjr(a,l);
		}
	}	

	//cw=0; cout<<'\n'<<"E_c_v: "<<'\n'; for(j=0;j<n;j++){for(k=0;k<col_wt;k++) cout<<E_c_v[CW*(col_wt*j+k)+cw]<<" "; cout<<'\n';}
	
}

