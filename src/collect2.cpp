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

typedef struct mapinfo{
	int chrm;
	int pos;
	int mapq;
	bool brevs;
	mapinfo(){pos = -1;};
	mapinfo(int c, int p, bool b, int m):chrm(c), pos(p), brevs(b), mapq(m){};
}mapinfo;

vector<vector<mapinfo> > vvmapinfo;
vector<vector<mapinfo> > vvmapinfo2;

int FileNum;
const int MaxFileNum=10;
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
Lib lib[MaxFileNum];//存储文库信息
int MaxLength;

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
	
	ofile.open("log/collect.log");
	
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
		MaxLength = lib[i].aveSpan + 5*lib[i].std;
	}
	
	if(!freopen("combine.fasta","r",stdin)){
		perror("combine.fasta");
	};
	string str;
	int nctgs = 0;
	vector<int> vsz;
	while(cin>>str)
	{
        cin>>str;
		nctgs++;
		vsz.push_back(str.size());
	}
	
	int covl[nctgs][MaxLength];
	for(int i = 0; i < nctgs; i++){
		for(int j = 0; j < MaxLength; j++)covl[i][j] = 0;
	}
	
	int covr[nctgs][MaxLength];
	for(int i = 0; i < nctgs; i++){
		for(int j = 0; j < MaxLength; j++)covr[i][j] = 0;
	}
	//cout<<"init covl[][] covr[][] over\n"<<endl;
	
	int readsnum;
	ifstream ifnum("mapping/reads_num.info");
	ifnum>>readsnum;
	ifnum.close();
	vvmapinfo.resize(readsnum);
	vvmapinfo2.resize(readsnum);
	
	int readsmapped1(0), totalreads1(0);
	ifstream ifsam("mapping/single.2.1.1.sam");
	
	int chrm, seqno, flag, pos, mapq;
	int mul1(0);
	bool brevs;
	const int szbuff = 1024;
	char buff[szbuff];
	
	while(ifsam.getline(buff,szbuff)){
		stringstream ss(buff);
		ss>>seqno>>flag>>chrm>>pos>>mapq;
		//cout<<seqno<<'\t'<<flag<<'\t'<<chrm<<'\t'<<pos<<endl;
		totalreads1++;
		if(! (flag & 0x04) ){
			brevs = flag & 0x10;
			readsmapped1++;
			vvmapinfo[seqno].push_back(mapinfo(chrm, pos, brevs, mapq));
		}
	}
	ifsam.close();
	cout<<"#1 totalreads1 = "<<totalreads1<<endl;
	cout<<"#1 readsmapped1 = "<<readsmapped1<<endl;
	
	int readsmapped2(0), totalreads2(0);
	int mul2(0);
	ifsam.open("mapping/single.2.1.2.sam");
	
	while(ifsam.getline(buff,szbuff)){
		stringstream ss(buff);
		ss>>seqno>>flag>>chrm>>pos>>mapq;
		totalreads2++;
		if(! (flag & 0x04) ){
			brevs = flag & 0x10;
			readsmapped2++;
			vvmapinfo2[seqno].push_back(mapinfo(chrm, pos, brevs, mapq));
		}
	}
	ifsam.close();
	cout<<"#2 totalreads2 = "<<totalreads2<<endl;
	cout<<"#2 readsmapped2 = "<<readsmapped2<<endl;
	
	
	ofstream ofclct("mapping/collect2");
	for(int i = 0; i < vvmapinfo.size(); i++){
		if(vvmapinfo[i].size()>1){
			mul1++;
		}
		if(vvmapinfo2[i].size()>1){
			mul2++;
		}
		if(vvmapinfo[i].empty() || vvmapinfo2[i].empty())continue;
		for(int p = 0; p < vvmapinfo[i].size(); p++){
			for(int q = 0; q < vvmapinfo2[i].size(); q++){
				if(vvmapinfo[i][p].chrm == vvmapinfo2[i][q].chrm)continue;
				ofclct<<vvmapinfo[i][p].chrm<<'\t'<<"2\t"<<vvmapinfo[i][p].pos<<'\t'<<vvmapinfo[i][p].brevs<<endl;
				ofclct<<vvmapinfo2[i][q].chrm<<'\t'<<"2\t"<<vvmapinfo2[i][q].pos<<'\t'<<vvmapinfo2[i][q].brevs<<endl;
			}
		}
	}
	
	cout<<endl;
	cout<<"#2 mul1 = "<<mul1<<endl;
	cout<<"#2 mul2 = "<<mul2<<endl;
	
	for(int i = 0; i < vvmapinfo.size(); i++){
		for(vector<mapinfo>::iterator it = vvmapinfo[i].begin(); it != vvmapinfo[i].end(); it++){
			if(vsz[it->chrm] - it->pos <= MaxLength){ // right side
				int start = vsz[it->chrm] - it->pos;
				while((start - vsz[it->chrm] - it->pos <= lib[FileNum].ReadLength) && start <= MaxLength){
					covr[it->chrm][start]++;
					start++;
				}
			}
			
			if(it->pos <= MaxLength){ // left side
				int start = it->pos;
				while((start - it->pos <= lib[FileNum].ReadLength) && start <= MaxLength){
					covl[it->chrm][start]++;
					start++;
				}
			}
		}
	}
	
	ofstream ofcov("mapping/cov");
	for(int i = 0; i < vsz.size(); i++){
		for(int j = 0; j < MaxLength; j++)ofcov<<covl[i][j]<<' ';
		ofcov<<endl;
		for(int j = 0; j < MaxLength; j++)ofcov<<covr[i][j]<<' ';
		ofcov<<endl;
	}
	return 0;
}