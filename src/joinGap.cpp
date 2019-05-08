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

struct tm* current_time;

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
const int LimitSz=9;
const int LimitSz2=9;
Lib lib[MaxFileNum];//存储文库信息
struct Pair_Str 
{
	string oristr;
	string revstr;
};
vector<Pair_Str> OriVec;

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
	
	ofile.open("log/joinGap.log");
	
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
	}
	
	ifstream ifctg("ctgpos.info");
	
	int pre = 0;
	int ctgno, len, gap;
	char c3('1'), c4(':');
	string ctgname;
	gap = 0;
	//int cnt = 0;
	unordered_map<string, unordered_map<int, int>> ctginfo; // <ctg name in ref, <ctgno, pos in ref> >
	unordered_map<string, unordered_map<int, int>> ctginfo2;
	unordered_map<int, string> ctgno2seq;
	while(ifctg>>ctgno>>len>>gap>>ctgname){
		replace(ctgname.begin(), ctgname.end(), ':', ' ');
		string seqname;
		int gb, startpos, endpos;
		stringstream ss(ctgname);
		ss>>seqname>>gb>>startpos>>endpos;
		ctginfo[seqname][ctgno] = startpos;
		ctginfo2[seqname][ctgno] = endpos;
		ctgno2seq[ctgno] = seqname;
	}
	
	ofstream ofscf("scaffold.fasta");
	ifstream ifs("order.info");
	int cnt = 0;
	
	int prectg = -1;
	int notSameSeq = 0;
	int wrongDis = 0;
	int corrt = 0;
	while(getline(ifs, str)){
		getline(ifs, str);
		stringstream ss(str);
		vector<string> vs;
		string strtmp;
		while(ss>>strtmp){
			vs.push_back(strtmp);
		}
		int ctg, nOfn, gap(0);
		bool brevs;
		ofscf<<">"<<cnt++<<endl;
		for(int i = 0; i < vs.size(); i+=3){
			nOfn = 0;
			if(i+2<vs.size()){
				ctg = atoi(vs[i].c_str());
				if(vs[i+1][0] == '-')brevs = 1;
				nOfn = atoi(vs[i+2].c_str());
			}else{
				ctg = atoi(vs[i].c_str());
				if(vs[i+1][0] == '-')brevs = 1;
			}
			if(brevs){
				if(gap>=0)ofscf<<OriVec[ctg].revstr;
				else{
					for(int h = abs(gap); h < OriVec[ctg].revstr.size();h++)ofscf<<OriVec[ctg].revstr[h];
				}
			}else{
				if(gap>=0)ofscf<<OriVec[ctg].oristr;
				else{
					for(int h = abs(gap); h < OriVec[ctg].oristr.size();h++)ofscf<<OriVec[ctg].oristr[h];
				}
			} 
			if(prectg<0){
				prectg = ctg;
			};
			string seq1, seq2;
			seq1 = ctgno2seq[prectg];
			seq2 = ctgno2seq[ctg];
			int end, start;
			if(seq1 == seq2){
				if(vs[i+1][0] == '-'){
					end = ctginfo2[seq1][ctg];
					start = ctginfo[seq1][prectg];
				}else{
					end = ctginfo2[seq1][prectg];
					start = ctginfo[seq1][ctg];
				}
				if(vs[i+1][0] == '-')gap = start - end;
				else gap = end - start;
				if(abs(nOfn - gap) > 1000)wrongDis++;
				else corrt++;
				cout<<nOfn<<'\t'<<gap<<'\t'<<nOfn - gap<<"\t# "<<ctg<<' '<<ctginfo[seq1][ctg]<<' '<<ctginfo2[seq1][ctg]<<endl;
			}else{
				notSameSeq++;
				gap = 0;
			}
			
			
			if(gap<0)gap = 0;
			while(gap--)ofscf<<'n';
			prectg = ctg;
		}
		ofscf<<endl;
	}
	
	cout<<"\n==============================="<<endl;
	cout<<notSameSeq<<'\t'<<wrongDis<<'\t'<<corrt<<endl;
	return 0;
}