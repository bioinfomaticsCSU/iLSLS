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
using namespace std;

//************** debug information ***************
#include<fstream>
#include<dirent.h>
ofstream ofile;
time_t ntm;
time_t prtm;
//************************************************

const bool useAll = 1;
const bool notUseAll = 0;

const int LimitSz=9;

struct tm* current_time;
struct spt{
	int posSelf;
	int posOther;
	int partdif;
	bool brvs;
	spt(){};
	spt(int ps,int po, int pd, int br):posSelf(ps), posOther(po), partdif(pd), brvs(br){}
};

struct target{
	int ctgno;
	int partdif;
	int gapdis;
	bool brevs;
	target(){};
	target(int c, int p, bool b):ctgno(c), partdif(p), brevs(b){};
};

typedef struct{
	deque<spt> vspt0;
	deque<spt> vspt0_1; //和vspt0的内容一样，只是元素顺序不一样
	deque<spt> vspt1;  //目标是反向互补的
	deque<spt> vspt1_1;
}dbvspt;

vector<bool> bused;
vector<int> finalctg;
vector<unordered_map<int, dbvspt> > vmpsptl;
vector<unordered_map<int, dbvspt> > vmpsptr;

vector<deque<int> > scfs;
vector<deque<int> > scfdirs;

struct Lib
{
	char FilePe[18];
	char FileSe[18];
	int MinSpan;//min_insertsize
	int MaxSpan;//max_insertsize
	int aveSpan;//ave_insertsize
	int MinOverLap;//K值
	int ReadLength;//读数长度
	int errolp;
	int std;//standard deviation
};

int FileNum;
const int MaxFileNum=10;
Lib lib[MaxFileNum];//存储文库信息
struct Pair_Str 
{
	string oristr;
	string revstr;
};
vector<Pair_Str> OriVec;
vector<int> vlen[2];

bool sptCmp(spt e1,spt e2)
{
	return e1.posSelf > e2.posSelf;
}

bool sptCmp2(spt e1,spt e2)
{
	return e1.posOther > e2.posOther;
}

bool targetCmp(target e1,target e2)
{
	return e1.gapdis < e2.gapdis;
}

string Reverse_Compliment(const string& ori)
{
	string ans=ori;
	int len=ori.size();
	for(int i=len-1;i>=0;i--)
	{
		if(ori[i]=='A') ans[len-i-1]='T';
		else if(ori[i]=='C') ans[len-i-1]='G';
		else if(ori[i]=='G') ans[len-i-1]='C';
		else if(ori[i]=='T') ans[len-i-1]='A';
	}
	return ans;
}

int main(){
	//*************** code add by lzx *****************
	DIR *dir=NULL;
	dir = opendir("log");
	if(!dir){
		system("mkdir log");
	}
	
	ofile.open("log/pairs.log");
	
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
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<" ~ "<<flkdf1<<' '<<vlen[!!flk1][ctgno1]<<' '<<pos1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<" ~ "<<flkdf2<<' '<<vlen[!!flk2][ctgno2]<<' '<<pos2<<endl;
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
					cout<<ctgno1<<' '<<flk1<<' '<<pos1<<' '<<dir1<<" ~ "<<flkdf1<<' '<<vlen[!!flk1][ctgno1]<<' '<<pos1<<endl;
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<" ~ "<<endl;
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
					cout<<ctgno2<<' '<<flk2<<' '<<pos2<<' '<<dir2<<" ~ "<<flkdf2<<' '<<vlen[!!flk2][ctgno2]<<' '<<pos2<<endl;
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
		ofile<<"==================================\n==================================\nvmpsptl:"<<endl;
		ofile<<"contig "<<i<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptl[i].begin(); it!=vmpsptl[i].end();){
			sort(it->second.vspt0.begin(),it->second.vspt0.end(),sptCmp);
			sort(it->second.vspt1.begin(),it->second.vspt1.end(),sptCmp);
			it->second.vspt0_1 = it->second.vspt0;
			it->second.vspt1_1 = it->second.vspt1;
			sort(it->second.vspt0_1.begin(),it->second.vspt0_1.end(),sptCmp2);
			sort(it->second.vspt1_1.begin(),it->second.vspt1_1.end(),sptCmp2);
			
			spts2 += it->second.vspt0.size();
			spts2 += it->second.vspt1.size();
			if(it->second.vspt0.size()>=LimitSz){
				ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
				for(deque<spt>::iterator itv = it->second.vspt0.begin();itv != it->second.vspt0.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			
			if(it->second.vspt1.size()>=LimitSz){
				ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
				for(deque<spt>::iterator itv = it->second.vspt1.begin();itv != it->second.vspt1.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			it++;
		}
		ofile<<"\n==================================\n==================================\n"<<endl<<endl;
		
		ofile<<"==================================\n==================================\nvmpsptr:"<<endl;
		ofile<<"contig "<<i<<endl;
		for(unordered_map<int, dbvspt >::iterator it = vmpsptr[i].begin(); it!=vmpsptr[i].end();){
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
			
			if(it->second.vspt0.size()>=LimitSz){
				ofile<<"candidate "<<it->first<<", brevs = false, surpports "<<it->second.vspt0.size()<<endl;
				for(deque<spt>::iterator itv = it->second.vspt0.begin();itv != it->second.vspt0.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			
			if(it->second.vspt1.size()>=LimitSz){
				ofile<<"candidate "<<it->first<<", brevs = true, surpports "<<it->second.vspt1.size()<<endl;
				for(deque<spt>::iterator itv = it->second.vspt1.begin();itv != it->second.vspt1.end(); itv++){
					ofile<<'('<<itv->posSelf<<','<<itv->posOther<<','<<itv->partdif<<"), ";
				}
				ofile<<endl;
			}
			it++;
		}
		ofile<<"==================================\n==================================\n"<<endl<<endl;
	}
	
	cout<<"spts2 = "<<spts2<<endl;
	cout<<"gt300 = "<<gt500<<endl;

	return 0;
}