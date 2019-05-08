#include<iostream>
#include<stdio.h>
#include<string>
#include<string.h>
#include<unordered_map>
#include<set>
#include<algorithm>
#include<vector>
#include<stack>
#include<cmath>
#include<sstream>
#include<deque>
#include"19.h"
using namespace std;


int main(){
	//*************** code add by lzx *****************
	DIR *dir=NULL;
	dir = opendir("log");
	if(!dir){
		system("mkdir log");
	}
	
	ofile.open("log/19.log");
	
	prtm = time(NULL); 
	ntm = prtm;
	struct tm* current_time = localtime(&ntm);
	ofile<<"program start...\t"<< ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	//************************************************
	
	if(!freopen("lib.info","r",stdin)){
		perror("lib.info");
	};
	scanf("%d",&FileNum);
	scanf("%d",&lib[0].MinOverLap);
	for(int i=1;i<=FileNum;i++)
	{
		scanf("%s%d%d%d%d%d",lib[i].FilePe,&lib[i].MinSpan,&lib[i].MaxSpan,&lib[i].std,&lib[i].ReadLength, &lib[i].errolp);
		lib[i].MinOverLap = lib[0].MinOverLap;
		lib[i].aveSpan=(lib[i].MinSpan+lib[i].MaxSpan)/2;
		lib[i].MaxSpan=lib[i].aveSpan+3*lib[i].std;
		lib[i].MinSpan=lib[i].aveSpan-3*lib[i].std;
		sprintf(lib[i].FilePe,"pe_%d.txt",i);
		sprintf(lib[i].FileSe,"se_%d.txt",i);
	}
	
	if(!freopen("combine.fasta","r",stdin)){
		perror("combine.fasta");
	};
	string str;
	Pair_Str ps;
	while(cin>>str)
	{
        cin>>str;
		ps.oristr=str;
		ps.revstr=Reverse_Compliment(ps.oristr);
		OriVec.push_back(ps);
		scfs.push_back(deque<int>());
		scfdirs.push_back(deque<int>());
		scfs[OriVec.size()-1].push_back(OriVec.size()-1);
		scfdirs[OriVec.size()-1].push_back(0);
	}
	bused.resize(OriVec.size());
	vmpsptl.resize(OriVec.size());
	vmpsptr.resize(OriVec.size());
	vlen[0].resize(OriVec.size());
	vlen[1].resize(OriVec.size());
	
	for(int i=0;i<bused.size();i++)bused[i] = false;
	
	ifstream ifbs("mapping/bases2.fa");
	while(ifbs>>str)
	{
		bool bsplit = 0;
		for(int i=0;i<str.size();i++){
			if(str[i] == '>'){
				str[i] = ' ';
			}
			if(str[i] == ':'){
				str[i] = '\t';
				bsplit = 1;
				break;
			}
		}
		stringstream ss;
		int ctgno,flk;
		if(bsplit){
			ss<<str;
			ss>>ctgno>>flk;
			ifbs>>str;
			vlen[flk][ctgno] = str.size();
		}else{
			ss<<str;
			ss>>ctgno;
			ifbs>>str;
			vlen[0][ctgno] = str.size();
			vlen[1][ctgno] = str.size();
		}
	}
	
	ifstream ifs("mapping/collect2");
	
	//seq [no] 0:left flank, 1:right flank, 2:full contig
	
	int ctgno1, flk1, pos1, dir1;
	int ctgno2, flk2, pos2, dir2;
	int flkdf1, flkdf2;
	
	int spts = 0;
	while(ifs>>ctgno1>>flk1>>pos1>>dir1){
		ifs>>ctgno2>>flk2>>pos2>>dir2;
		spts++;
		if(ctgno1 == 0 && ctgno2 == 115 || ctgno1 == 115 && ctgno2 == 0){
			// cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
			// cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
			// cout<<"-------------------------"<<endl;
			// continue;
		}
		
		if(ctgno1>OriVec.size()||ctgno2>OriVec.size())cout<<OriVec.size()<<'\t'<<ctgno1<<'\t'<<ctgno2<<endl;
		if(dir1){
			if(dir2){ //dif strand 1 -> -2
				//if(flk1 == 0 || flk2 == 1)continue;
				flkdf1 = vlen[!!flk1][ctgno1]-pos1;
				flkdf2 = vlen[!!flk2][ctgno2]-pos2;
				if(flkdf1<0||flkdf2<0){
					cout<<"1"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno1][ctgno2].vspt1.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 1)); 
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno2][ctgno1].vspt1.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 1));
			}else{ //same strand 1 -> 2
				//if(flk1 == 0 || flk2 == 1)continue;
				flkdf1 = vlen[!!flk1][ctgno1]-pos1;
				flkdf2 = pos2+1; 
				if(flkdf1<0||flkdf2<0){
					cout<<"2"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno1][ctgno2].vspt0.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 0)); 
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno2][ctgno1].vspt0.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 0));
			}
		}else{
			if(dir2){ // same strand 2 -> 1
				//if(flk1 == 1 || flk2 == 0)continue;
				flkdf1 = pos1+1;
				flkdf2 = vlen[!!flk2][ctgno2]-pos2;
				if(flkdf1<0||flkdf2<0){
					cout<<"3"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptl[ctgno2][ctgno1].vspt0.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 0));
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno1][ctgno2].vspt0.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 0));
			}else{ //dif strand -2 -> 1
				//if(flk1 == 1 || flk2 == 0)continue;
				flkdf1 = pos1+1;
				flkdf2 = pos2+1;
				if(flkdf1<0||flkdf2<0){
					cout<<"4"<<endl;
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<endl;
					exit(-1);
				}
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno1][ctgno2].vspt1.push_back(spt(flkdf1, flkdf2, flkdf1 + flkdf2, 1)); 
				if(flkdf1 + flkdf2 <lib[FileNum].MaxSpan)vmpsptr[ctgno2][ctgno1].vspt1.push_back(spt(flkdf2, flkdf1, flkdf1 + flkdf2, 1));
			}
		}
	}
	//exit(0);
	cout<<"spts = "<<spts<<endl;
	//************** debug information ***************
	prtm = ntm; 
	ntm = time(NULL);
	current_time = localtime(&ntm);
	ofile<<"read over, "<<ntm-prtm <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
	ofile<<"------------------------------------------"<<endl;
	prtm = ntm; 
	//************************************************
	
	int spts2 = 0;
	int gt500 = 0;
	for(int i=0;i<OriVec.size();i++){
		for(unordered_map<int, dbvspt >::iterator it = vmpsptl[i].begin(); it!=vmpsptl[i].end(); it++){
			sort(it->second.vspt0.begin(),it->second.vspt0.end(),sptCmp);
			sort(it->second.vspt1.begin(),it->second.vspt1.end(),sptCmp);
			it->second.vspt0_1 = it->second.vspt0;
			it->second.vspt1_1 = it->second.vspt1;
			sort(it->second.vspt0_1.begin(),it->second.vspt0_1.end(),sptCmp2);
			sort(it->second.vspt1_1.begin(),it->second.vspt1_1.end(),sptCmp2);
			
			spts2 += it->second.vspt0.size();
			spts2 += it->second.vspt1.size();
		}
		for(unordered_map<int, dbvspt >::iterator it = vmpsptr[i].begin(); it!=vmpsptr[i].end(); it++){
			sort(it->second.vspt0.begin(),it->second.vspt0.end(),sptCmp);
			sort(it->second.vspt1.begin(),it->second.vspt1.end(),sptCmp);
			it->second.vspt0_1 = it->second.vspt0;
			it->second.vspt1_1 = it->second.vspt1;
			sort(it->second.vspt0_1.begin(),it->second.vspt0_1.end(),sptCmp2);
			sort(it->second.vspt1_1.begin(),it->second.vspt1_1.end(),sptCmp2);
			
			spts2 += it->second.vspt0.size();
			spts2 += it->second.vspt1.size();
			if(it->second.vspt0.size()>300)gt500++;
			if(it->second.vspt1.size()>300)gt500++;
		}
	}
	
	cout<<"spts2 = "<<spts2<<endl;
	cout<<"gt300 = "<<gt500<<endl;
	
	int nwhile = 0;
	time_t ntm2,prtm2,ntm3,prtm3;
	ntm2 = time(NULL);
	ntm3 = time(NULL);
	int lastno, no;
	while(1)
	{
		nwhile++;
		lastno = 0;
		no = 0;
		
		for(int i=0;i<OriVec.size();i++)if(!bused[i])lastno++;
		
		for(int i=0;i<OriVec.size();i++){
			if(!bused[i])
				ofile<<"(i, nContigs) = "<<"("<<i<<", "<<OriVec.size()<<")"<<endl;
			else continue;
			bused[i] = true;
			finalctg.push_back(i);
			bool bRightDone = false;
			while(1){
				deque<target> vtgt;
				
				//right extension  unordered_map<int, vector<spt> >
				int rightedge = scfs[i][scfs[i].size()-1];
				bool brevs = scfdirs[i][scfdirs[i].size()-1];
				
				collectTargets(i, vtgt, rightedge, brevs, notUseAll, true);
				
				ofile<<"final candidates are: \n";
				for(int j=0;j<vtgt.size();j++){
					ofile<<vtgt[j].ctgno<<'\t'<<vtgt[j].partdif<<'\t'<<vtgt[j].gapdis<<'\t'<<vtgt[j].brevs<<'\t'<<OriVec[vtgt[j].ctgno].oristr.size()<<endl;
				}
				if(!scf(vtgt, i, rightedge, 1)){
					if(bRightDone){
						ofile<<"left extension over!"<<endl;
						break;
					}else{
						bRightDone = true;
						ofile<<"right extension over!"<<endl;
						reverse(scfs[i].begin(), scfs[i].end());
						reverse(scfdirs[i].begin(), scfdirs[i].end());
						for(int j = 0;j<scfdirs[i].size();j+=2){
							scfdirs[i][j] = !scfdirs[i][j];
						}
					}
				}
			}
			
			
			
			
			//************** debug information ***************
			prtm2 = ntm2; 
			ntm2 = time(NULL);
			current_time = localtime(&ntm2);
			ofile<<"(i, nContigs) = "<<"("<<i<<", "<<OriVec.size()<<"), "<<ntm2-prtm2 <<" (s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
			ofile<<"------------------------------------------"<<endl;
			//************************************************
		}
		
		for(int i=0;i<OriVec.size();i++)if(!bused[i])no++;
		//************** debug information ***************
		prtm3 = ntm3; 
		ntm3 = time(NULL);
		current_time = localtime(&ntm3);
		ofile<<"scaffold and gap round "<<nwhile<<" over,\tscaffold success "<<lastno-no<<" times, "<<ntm3-prtm3 <<"(s)\t"<<current_time->tm_year+1900<<"/"<<current_time->tm_mon<<"/"<<current_time->tm_mday<<" "<<current_time->tm_hour<<":"<<current_time->tm_min<<":"<<current_time->tm_sec<<endl;
		ofile<<"------------------------------------------"<<endl;
		//************************************************
		//if(0 == lastno-no)break;
		break;
	} 
	
	for(int q = 0;q < finalctg.size(); q++){
		int r = finalctg[q];
		if(scfs[r].size() == 1){
			deque<target> vtgt, vtgt2;
			//right extension  unordered_map<int, vector<spt> >
			int rightedge = r;
			bool brevs = false;
			
			collectTargets(rightedge, vtgt, rightedge, false, useAll, false);
			if(vtgt.size())ofile<<"vtgt, final candidates are: \n";
			for(int j=0;j<vtgt.size();j++){
				ofile<<vtgt[j].ctgno<<'\t'<<vtgt[j].partdif<<'\t'<<vtgt[j].gapdis<<'\t'<<vtgt[j].brevs<<'\t'<<OriVec[vtgt[j].ctgno].oristr.size()<<endl;
			}
			
			if(vtgt.size() == 1 && vtgt[0].ctgno != r){
				cout<<"1 : "<<r<<' '<<vtgt[0].ctgno<<endl;
				insert(vtgt[0], r, 1);
				finalctg.erase(finalctg.begin()+q);
				continue;
			}
			
			collectTargets(rightedge, vtgt2, rightedge, true, useAll, false);
			if(vtgt2.size())ofile<<"vtgt2,final candidates are: \n";
			for(int j=0;j<vtgt2.size();j++){
				ofile<<vtgt2[j].ctgno<<'\t'<<vtgt2[j].partdif<<'\t'<<vtgt2[j].gapdis<<'\t'<<vtgt2[j].brevs<<'\t'<<OriVec[vtgt2[j].ctgno].oristr.size()<<endl;
			}
			if(vtgt2.size() == 1 && vtgt2[0].ctgno != r){
				cout<<"2: "<<r<<' '<<vtgt2[0].ctgno<<endl;
				//insert(vtgt2[0], r, 0);
				//finalctg.erase(finalctg.begin()+q);
			}
		}
	}
	cout<<"process insert over!"<<endl;
	
	ofstream ofscf("scaffold.fasta");
	ofstream oforder("order19.info");
	int cnt = 0;
	for(int q = 0;q < finalctg.size(); q++){
		int r = finalctg[q];
		ofscf<<">"<<cnt++<<endl;
		oforder<<"scafffold "<<r<<" :"<<endl;
		for(int h= 0;h < scfs[r].size();h++){
			if(h%2 == 0){ // contig
				oforder<<scfs[r][h];
				if(scfdirs[r][h]) {
					oforder<<'-';
					ofscf<<OriVec[scfs[r][h]].revstr;
				}else {
					oforder<<'+';
					ofscf<<OriVec[scfs[r][h]].oristr;
				}
				oforder<<'\t';
			}else{ // gap N
				while(scfdirs[r][h]--)ofscf<<'N';
			}
		}
		oforder<<endl;
		ofscf<<endl;
	}
	
	return 0;
}